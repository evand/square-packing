#!/usr/bin/env python3
"""DUP2 (10-18): duplicate-half of CURRENT asymmetric records (the de Winter / Schadt direction, from the record).
Non-integer store-best records with tilts that are not Md / Ma invariant, ranked by Hungarian mirror mismatch under
their better diagonal mirror; top K.  Per record: dup-lo, dup-hi; symmetric soft-to-rigid relax (r0 0.05 / 0.1 / 0.2 x 2
seeds) + free polish.  Released side vs the record (frontier).
  symdup2.py --out runs/symsub/dup2 [--top 60]"""
import sys, os, json, argparse, tempfile, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import symtt
from symdup import job


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--top', type=int, default=60)
    a = ap.parse_args()
    from multiprocessing import Pool
    from pk.store import Store
    sys.path.insert(0, 'rational'); import ratgen
    fr = ratgen.frontier()
    os.makedirs(a.out, exist_ok=True)
    sc = json.load(open('runs/symsub/symscan.json')); st = Store()
    cand = []
    for n, r in sc.items():
        n = int(n)
        if n < 30 or not r['ntilt'] or abs(r['s'] - round(r['s'])) < 1e-9 or {'Md', 'Ma'} & set(r['inv']):
            continue
        p = st.get(r['pid'])
        el, d = symtt.best_diag_mirror(p.sq, p.s)
        cand.append((d, n, r['pid'], p, el))
    cand.sort(key=lambda c: c[0])
    print(len(cand), 'asymmetric records; top', a.top, 'mismatch', f'{cand[0][0]:.3f} .. {cand[a.top - 1][0]:.3f}', flush=True)
    T = []
    for d, n, pid, p, el in cand[:a.top]:
        rec = min(fr.get(n, p.s), p.s)
        for how in ('dup-lo', 'dup-hi'):
            for r0 in (0.05, 0.1, 0.2):
                for seed in (1, 2):
                    T.append((n, rec, pid, p.s, p.sq.tolist(), how, el, r0, seed))
    with Pool(a.procs) as pool:
        R = pool.map(job, T)
    json.dump(R, open(os.path.join(a.out, 'runs.json'), 'w'))
    g = collections.defaultdict(list)
    for r in R:
        if r['s'] is not None:
            g[(r['n'], r['how'])].append(r['s'] - r['rec'])
    mm = {n: d for d, n, *_ in cand}
    rows = []
    for n in sorted({k[0] for k in g}):
        b = min(min(g[(n, 'dup-lo')], default=9), min(g[(n, 'dup-hi')], default=9))
        rows.append((b, n))
    print('best released side - record, by n (mismatch):')
    for b, n in sorted(rows):
        print(n, f'mismatch {mm[n]:.3f}', f'lo {min(g[(n, "dup-lo")], default=float("nan")):+.5f} hi {min(g[(n, "dup-hi")], default=float("nan")):+.5f}',
              '<<< BELOW' if b < -1e-9 else ('= tie' if abs(b) <= 1e-7 else ''))
