#!/usr/bin/env python3
"""Record the nlpdesc pass over the s(110) unresolved side keys and weak minima (10-15) into the store.

Inputs: runs/resolve110_nlp/{unresolved,weak}/<key>/ (in.txt = the store representative, nlp.out.txt = nlpdesc output,
nlp.out.txt.ex/ = exactsolve on it) and final.json (one row per target: kind, dir, cls, S, s_in = output side, s0, gain).
Writes, as descend_census does: the landed packing (source 'nlpdesc:resolve110'), the minimum row (certified / weak,
merged within 1e-17), the copy link input side key -> minimum.  Basin rows of the unresolved side keys: 'not-min' when
the descent moved the side by > 1e-9 (status names the landing S) or exactsolve found a descent there; otherwise the
landed minimum's class (the side key was an unpolished copy).  Weak minimum rows are not changed (the representative
copies are linked by a 1e-9 side key only, so a descent from one is evidence, not proof, about the weak point).

  python3 resolve110_nlp.py [--dry]
"""
import argparse, collections, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pk.store import Store, now
from pk.packing import read
from pk import identity as ID

ROOT = os.path.join(HERE, 'runs', 'resolve110_nlp')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dry', action='store_true'); a = ap.parse_args()
    st = Store(); R = json.load(open(f'{ROOT}/final.json')); C = collections.Counter()
    for r in R:
        d = f"{ROOT}/{r['kind']}/{r['dir']}"
        landed = r['cls'] in (ID.CERTIFIED, ID.WEAK) and r.get('S')
        mk = None
        if landed:
            known = st.find_minimum(110, r['S']) is not None
            C[f"{r['kind']}:{r['cls']}-{'known' if known else 'new'}"] += 1
            if not a.dry:
                pid2, _ = st.add(read(f'{d}/nlp.out.txt'), source='nlpdesc:resolve110', finder='ours',
                                 ref=f"from {r['kind']} {r['dir']}", status=r['cls'], commit=False)
                mk, _ = st.add_minimum(110, r['S'], r['cls'], r.get('status'), src='nlpdesc:resolve110', pid=pid2, commit=False)
                st.add_copy(110, ID.side_key(float(r['s_in'])), mk, via='nlpdesc', commit=False)
                if r['kind'] == 'unresolved':
                    st.add_copy(110, r['dir'], mk, via='nlpdesc', commit=False)
        else:
            C[f"{r['kind']}:{r['cls']}"] += 1
        if r['kind'] != 'unresolved' or a.dry:
            continue
        if r['cls'] == 'not-min':
            row = ('not-min', None, f"NOT a local minimum: nlpdesc descent to {r['s_in']} (still not-min there)")
        elif landed and r['gain'] > 1e-9:
            row = ('not-min', None, f"NOT a local minimum: nlpdesc descent (-{r['gain']:.2e}) to {r['cls']} S = {r['S'][:22]}")
        elif landed:
            row = (r['cls'], r['S'], r.get('status'))
        else:
            continue                                  # still unresolved: leave the basin row
        st.db.execute('UPDATE basin SET cls=?, S=?, status=?, added=? WHERE n=110 AND key=?', (*row, now(), r['dir']))
    if not a.dry:
        st.db.commit()
    print(dict(sorted(C.items())))
    print('basin:', [tuple(x) for x in st.db.execute("SELECT cls, count(*) FROM basin WHERE n=110 GROUP BY cls")])
    print('minimum < 11:', [tuple(x) for x in st.db.execute(
        "SELECT cls, count(*) FROM minimum WHERE n=110 AND CAST(S AS REAL) < 11 GROUP BY cls")])


if __name__ == '__main__':
    main()
