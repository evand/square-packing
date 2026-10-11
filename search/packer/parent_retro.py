#!/usr/bin/env python3
"""Retrospective (10-11, Evan's question): are partially descended (not-min) states useful parents?

Observational, from logged explorer runs at n = 110 (runs/lineage/L1, L2, runs/insert110/x*): every sub-11 archive entry
gets the store's verdict for its side key (basin: certified / weak / not-min / unresolved) and, via the copy table, the
exact minimum it is (certified) or descends to (not-min points descended by descend_census).  Each proposal is credited to
its parent's class.  Child outcomes (child = 'new?' result below 11, looked up by side key):
  cert   the child is a certified minimum
  moved  ... and a different minimum from the parent's own (parents with a known minimum only)
  disc   ... and the first time this run reaches that minimum (start entries' minima excluded): a within-run discovery
Rates per 1000 proposals and per worker-hour (logged 'sec'; these runs predate measured CPU), by parent class and side gap.
Matched: not-min parents vs the certified parent that is their own descended minimum, in the same run.
Confounds: parents were chosen by elite / value rules, not at random; per-run archives differ.  Use it to decide whether
census v2 needs a controlled arm, not as a verdict.

  parent_retro.py [--json runs/parent_retro.json]
"""
import argparse, collections, glob, json, os, sqlite3, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REC = 10.9967833966316
BINS = [(0, 1e-4), (1e-4, 1e-3), (1e-3, 3e-3), (3e-3, 1.0)]


def gbin(s):
    if s >= 11: return 'above 11'
    g = s - REC
    for lo, hi in BINS:
        if lo <= g < hi: return f'{lo:g}-{hi:g}'
    return 'other'


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--json', default=f'{HERE}/runs/parent_retro.json'); a = ap.parse_args()
    db = sqlite3.connect(f'file:{HERE}/runs/store.sqlite?mode=ro', uri=True)
    B = dict(db.execute('SELECT key, cls FROM basin WHERE n=110').fetchall())
    CP = collections.defaultdict(set)
    for sk, mk in db.execute('SELECT side_key, min_key FROM copy WHERE n=110'):
        CP[sk].add(mk)
    def lab(s):
        if s >= 11: return 'above-11', None
        k = f'{s:.9f}'
        m = CP.get(k)
        return B.get(k, 'unlabelled'), (next(iter(m)) if m and len(m) == 1 else None)

    T = collections.defaultdict(lambda: collections.Counter())      # (class, bin) -> counters
    matched = collections.Counter()
    runs = sorted(d for d in glob.glob(f'{HERE}/runs/lineage/L*/*') + glob.glob(f'{HERE}/runs/insert110/x*')
                  if os.path.exists(f'{d}/proposals.jsonl'))
    for d in runs:
        E = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
        L = {e['i']: lab(e['s']) for e in E}
        seen = {L[e['i']][1] for e in E if e['parent'] < 0 and L[e['i']][1]}          # start entries' minima
        P = [json.loads(l) for l in open(f'{d}/proposals.jsonl')]
        per_parent = collections.defaultdict(lambda: collections.Counter())
        for p in P:
            par = E[p['par']] if p['par'] < len(E) else None
            if par is None: continue
            pc, pm = L[par['i']]
            c = T[(pc, gbin(par['s']))]
            c['props'] += 1; c['sec'] += p.get('sec', 0.0)
            pp = per_parent[par['i']]; pp['props'] += 1
            if p.get('st') != 'new?' or p.get('s') is None or p['s'] >= 11:
                continue
            cc, cm = lab(p['s'])
            if cc != 'certified':
                continue
            c['cert'] += 1; pp['cert'] += 1
            if pm is not None and cm is not None:
                c['known_par_min'] += 1
                if cm != pm: c['moved'] += 1; pp['moved'] += 1
            if cm is not None and cm not in seen:
                seen.add(cm); c['disc'] += 1; pp['disc'] += 1
        # matched pairs: a not-min parent and the certified parent that is its own minimum, both expanded here
        cert_by_min = {}
        for e in E:
            c_, m_ = L[e['i']]
            if c_ == 'certified' and m_: cert_by_min.setdefault(m_, e['i'])
        for e in E:
            c_, m_ = L[e['i']]
            if c_ == 'not-min' and m_ in cert_by_min:
                a_, b_ = per_parent.get(e['i']), per_parent.get(cert_by_min[m_])
                if a_ and b_ and a_['props'] and b_['props']:
                    matched['pairs'] += 1
                    for k in ('props', 'cert', 'moved', 'disc'):
                        matched[f'nm_{k}'] += a_[k]; matched[f'c_{k}'] += b_[k]
    out = {}
    print(f'{len(runs)} runs; rates per 1000 proposals (per worker-hour in brackets)')
    print(f'{"parent class":12} {"gap bin":12} {"props":>7} {"cert":>13} {"moved":>13} {"disc":>13}  moved/known')
    for (pc, gb), c in sorted(T.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        if c['props'] < 200: continue
        h = c['sec'] / 3600 or 1
        f = lambda k: f'{1000 * c[k] / c["props"]:5.1f} ({c[k] / h:5.1f})'
        mk = f'{c["moved"]}/{c["known_par_min"]}' if c['known_par_min'] else '-'
        print(f'{pc:12} {gb:12} {c["props"]:7d} {f("cert"):>13} {f("moved"):>13} {f("disc"):>13}  {mk}')
        out[f'{pc}|{gb}'] = dict(c)
    for pc in ('certified', 'weak', 'not-min', 'above-11'):
        c = sum((T[k] for k in T if k[0] == pc), collections.Counter())
        if c['props']:
            h = c['sec'] / 3600 or 1
            print(f'  all {pc:10}: {c["props"]:6d} props, cert {1000*c["cert"]/c["props"]:.1f}, moved {1000*c["moved"]/c["props"]:.1f}, '
                  f'disc {1000*c["disc"]/c["props"]:.2f} per 1000 ({c["disc"]/h:.2f} disc per worker-h)')
    if matched['pairs']:
        m = matched
        print(f'matched: {m["pairs"]} (not-min parent, its own certified minimum) pairs expanded in the same run')
        for k in ('cert', 'moved', 'disc'):
            print(f'  {k}: not-min {1000*m["nm_"+k]/m["nm_props"]:.1f} vs certified {1000*m["c_"+k]/m["c_props"]:.1f} per 1000 '
                  f'({m["nm_"+k]} / {m["nm_props"]} vs {m["c_"+k]} / {m["c_props"]})')
    out['matched'] = dict(matched)
    json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
