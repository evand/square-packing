#!/usr/bin/env python3
"""Spot check (10-11): do tour legs between differently oriented minima overstate the expansion?  morph.match2 matches
squares without box symmetries.  For the K mixed-orientation tour legs with the largest expansion: re-orient B by the box
symmetry minimising match2's squared-corner-movement cost, rerun build._barrier (same bisection), compare.

  python3 landscape/orient_check.py --k 5 --procs 5 --tmp DIR        (run from packer/; ~1 min per pair per proc)"""
import argparse, json, os, sys
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
from scipy.optimize import linear_sum_assignment
import morph, build
from mapdata import sym, SYMS
RUN = os.path.join(os.path.dirname(HERE), 'runs/landscape110')


def loaddeg(p):
    L = [l.split() for l in open(p) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1]); return s, np.array([[float(v) for v in l[:3]] for l in L[1:n + 1]])


def cost(A, B):
    CA, CB = morph.corners(A[1]), morph.corners(B[1])
    d = np.min([(np.linalg.norm(CA[:, None] - np.roll(CB, k, axis=1)[None], axis=-1) ** 2).sum(-1) for k in range(4)], 0)
    r, c = linear_sum_assignment(d); return d[r, c].sum()


def best_orient(pa, pb, tmp):
    A = morph.load(pa); s, sq = loaddeg(pb); res = []
    for g in SYMS:
        q = sym(sq, s, *g); B = (s, [(x, y, np.radians(t)) for x, y, t in q]); res.append((cost(A, B), g, q))
    res.sort(key=lambda r: r[0]); c0 = [r for r in res if r[1] == (0, 0, 0)][0][0]
    c, g, q = res[0]
    out = os.path.join(tmp, os.path.basename(pb).replace('.txt', f'_g{g[0]}{g[1]}{g[2]}.txt'))
    with open(out, 'w') as f:
        f.write(f'110 {float(s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in q))
    return g, c0, c, out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--k', type=int, default=5); ap.add_argument('--procs', type=int, default=5)
    ap.add_argument('--tmp', required=True); a = ap.parse_args()
    R = json.load(open(f'{RUN}/reps.json')); cyc = json.load(open(f'{RUN}/cycle.json'))
    path = lambda i: R[i]['path'] if os.path.isabs(R[i]['path']) else os.path.join(os.path.dirname(HERE), R[i]['path'])
    legs = []
    for l in cyc['legs']:
        g, c0, c, out = best_orient(path(l['i']), path(l['j']), a.tmp)
        if g != (0, 0, 0) and l.get('s'):
            legs.append((l['s'] - max(float(R[l['i']]['S']), float(R[l['j']]['S'])), l, g, c0, c, out))
    legs.sort(key=lambda r: -r[0])
    print(f'{len(legs)} tour legs where re-orienting B lowers the matching cost; checking the {a.k} largest expansions')
    jobs = [(l['i'], l['j'], path(l['i']), out) for _, l, *_, out in legs[:a.k]]
    with Pool(a.procs) as P:
        new = P.map(build._barrier, jobs)
    for (ex, l, g, c0, c, out), r in zip(legs[:a.k], new):
        base = r.get('base'); nex = None if r['s'] is None else r['s'] - base
        print(f"leg {l['i']:2d}->{l['j']:2d}: stored expansion {ex:.4f}; matching cost {c0:.2f} -> {c:.2f} with sym {g}; "
              f"re-oriented expansion {'fail' if nex is None else f'{nex:.4f}'}")


if __name__ == '__main__':
    main()
