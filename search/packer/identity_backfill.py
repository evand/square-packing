#!/usr/bin/env python3
"""One-off (10-11): fill the store's `minimum` / `copy` tables (pk/identity.py) from existing verdicts, and relabel
basin-cache rows whose stored status shows a weak certificate (PSD with zero modes / inconclusive) as cls 'weak'.

Sources of verdicts (exactsolve status strings, keyed by exact S):
  basin table rows (cls certified, S stored; status may be 'import' = written by an unsaved 10-09 command, status lost)
  runs/**/exact.json (explorer certify), runs/**/census_exact.json (census3), runs/census/n*/*/class.json,
  runs/census/n*/descend.jsonl (descended not-min points -> copy links)
A certified basin row whose status is lost and not recoverable from the files is kept as certified with status
'import (second order not recorded)' and counted separately in the report (re-solve to settle it).

  identity_backfill.py [--dry]      (writes to the store unless --dry; snapshot the store first)
"""
import argparse, collections, glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pk.store import Store
from pk import identity as ID

LOST = 'import (second order not recorded)'


def status_index():
    """min_key -> list of status strings seen for that exact S anywhere in the saved verdict files."""
    idx = collections.defaultdict(list)
    def put(S, st):
        if S and st:
            idx[ID.min_key(S)].append(st)
    for f in glob.glob(f'{HERE}/runs/**/exact.json', recursive=True):
        if '/ext/' in f: continue
        try: d = json.load(open(f))
        except Exception: continue
        for v in (d.values() if isinstance(d, dict) else d):
            if isinstance(v, dict): put(v.get('S'), v.get('status'))
    for f in glob.glob(f'{HERE}/runs/**/census_exact.json', recursive=True):
        try: d = json.load(open(f))
        except Exception: continue
        for r in d:
            e = r.get('ex') if isinstance(r, dict) else None
            if isinstance(e, dict): put(e.get('S'), e.get('status'))
    for f in glob.glob(f'{HERE}/runs/census/n*/*/class.json'):
        try: r = json.load(open(f))
        except Exception: continue
        put(r.get('S'), r.get('status'))
    for f in glob.glob(f'{HERE}/runs/census/n*/descend.jsonl'):
        for l in open(f):
            r = json.loads(l); put(r.get('S'), r.get('status'))
    return idx


def verdict(status, idx_list):
    """Best verdict: the row's own status if informative, else the best one among the indexed statuses."""
    cands = [status] if status and not status.startswith(('import', 'known110')) else []
    cands += idx_list
    best = None
    for st in cands:
        c = ID.classify_status(st)
        if c in (ID.CERTIFIED, ID.WEAK) and (best is None or ID.RANK[c] > ID.RANK[best[0]]):
            best = (c, st)
    return best


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dry', action='store_true'); a = ap.parse_args()
    st = Store()
    idx = status_index()
    print(f'status index: {len(idx)} exact sides with a saved status')
    C = collections.Counter()
    rows = st.db.execute("SELECT n, key, cls, S, status FROM basin WHERE cls='certified'").fetchall()
    for r in rows:
        if not r['S']:
            C['certified basin row without S (skipped)'] += 1; continue
        v = verdict(r['status'], idx.get(ID.min_key(r['S']), []))
        if v is None:
            cls, stt = ID.CERTIFIED, LOST; C['status lost -> certified (unverified class)'] += 1
        else:
            cls, stt = v; C[f'basin -> {cls}'] += 1
        if not a.dry:
            if cls == ID.WEAK:
                st.db.execute("UPDATE basin SET cls='weak' WHERE n=? AND key=?", (r['n'], r['key']))
            mk, new = st.add_minimum(r['n'], r['S'], cls, stt, src='backfill:basin', commit=False)
            st.add_copy(r['n'], r['key'], mk, via='backfill:basin', commit=False)
    for f in sorted(glob.glob(f'{HERE}/runs/census/n*/descend.jsonl')):
        n = int(f.split('/n')[-1].split('/')[0])
        for l in open(f):
            r = json.loads(l)
            if not r.get('S') or r.get('cls') != 'certified':
                continue
            c = ID.classify_status(r.get('status'))
            if c not in (ID.CERTIFIED, ID.WEAK):
                c = ID.CERTIFIED; C['descend: status not classifiable -> certified'] += 1
            C[f'descend -> {c}'] += 1
            if not a.dry:
                mk, _ = st.add_minimum(n, r['S'], c, r.get('status'), src='backfill:descend', pid=r.get('pid_out'), commit=False)
                st.add_copy(n, r['key'], mk, via='descend', commit=False)
    if not a.dry:
        st.db.commit()
    for k, v in sorted(C.items()):
        print(f'  {k}: {v}')
    for n in sorted({r['n'] for r in st.db.execute('SELECT DISTINCT n FROM minimum')}):
        print(f'  n={n}: minimum rows', dict(collections.Counter(m['cls'] for m in st.minima(n))))


if __name__ == '__main__':
    main()
