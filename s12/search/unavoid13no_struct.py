#!/usr/bin/env python3
"""unavoid13-no step 3: h on structured subfamilies (axis-parallel, +45deg, +near-axis, +Pythagorean
angles), with the shape of the optimal sets.

    python3 search/unavoid13no_struct.py --stages A0,A1,A2,A3 --grid 10 --out runs/unavoid13no_struct

Stages (each D4-closed):
  A0  axis-parallel squares with centres on the 1/grid lattice (contains the tiles and every shifted
      4x3 / 3x4 / 3x3 tiling at that shift resolution);
  A1  + 45deg squares (u = tan 22.5deg is irrational: u = 5/12 -> theta = 45.24deg and its mirror
      44.76deg are used, both rational), centres on the lattice clamped into the admissible range;
  A2  + near-axis angles u = 1/57, 1/28, 1/19 (theta = 2.0, 4.1, 6.0deg) and mirrors;
  A3  + Pythagorean angles u = 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3 (rational cos, sin).
"""
import sys, os, time, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13no_lib as N


def angle_grid(m, u, K):
    c, s = L.trig(u); w = abs(c) + abs(s)
    lo, hi = w / 2, m - w / 2
    g = [lo] + [F(k, K) for k in range(int(lo * K) + 1, int(hi * K) + 1)] + [hi]
    g = sorted(set(x for x in g if lo <= x <= hi))
    return [(u, x, y) for x in g for y in g]


STAGES = {
    'A0': [F(0)],
    'A1': [F(5, 12)],
    'A2': [F(1, 57), F(1, 28), F(1, 19)],
    'A3': [F(1, 8), F(1, 7), F(1, 6), F(1, 5), F(1, 4), F(1, 3), F(2, 5), F(1, 2), F(3, 5), F(2, 3)],
}


def describe(P):
    P = np.asarray(P, float)
    lat = sum(1 for x, y in P if abs(x - round(x)) < 1e-6 and abs(y - round(y)) < 1e-6)
    edge = sum(1 for x, y in P if (abs(x - round(x)) < 1e-6) != (abs(y - round(y)) < 1e-6))
    return f"{len(P)} points: {lat} at lattice points, {edge} on tile edges, {len(P)-lat-edge} interior"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', default='4'); ap.add_argument('--stages', default='A0,A1,A2,A3')
    ap.add_argument('--grid', type=int, default=10); ap.add_argument('--out', default='runs/unavoid13no_struct')
    ap.add_argument('--threads', type=int, default=4); ap.add_argument('--ip-time', type=float, default=3600)
    ap.add_argument('--k', type=int, default=13, help='also run the feasibility IP for k points')
    a = ap.parse_args()
    m = F(a.m); os.makedirs(a.out, exist_ok=True)
    Fam = []
    for st in a.stages.split(','):
        for u in STAGES[st]:
            Fam += angle_grid(m, u, a.grid)
        Fam = L.d4_closure(m, Fam)
        tag = f"{st}_g{a.grid}"
        L.write_family(f"{a.out}/{tag}_family.txt", m, Fam, header=f"stage {st}, grid 1/{a.grid}")
        t0 = time.time()
        sq, V, B = N.maximal_cells(m, Fam, threads=a.threads)
        lp, _ = N.lp_value(B, threads=a.threads)
        print(f"[{tag}] {len(Fam)} squares, {B.shape[0]} cells, LP {lp:.4f}", flush=True)
        t1 = time.time()
        r = L.solve_hitting_set(B, threads=a.threads, time_limit=a.ip_time)
        print(f"[{tag}] h = {r['obj']} ({r['status']}, bound {r['bound']}, {time.time()-t1:.0f}s)", flush=True)
        if r['x'] is not None:
            P = V[np.nonzero(r['x'] > 0.5)[0]]
            print(f"[{tag}] optimal set: {describe(P)}: " + " ".join(f"({x:.3f},{y:.3f})" for x, y in sorted(map(tuple, P))), flush=True)
            np.save(f"{a.out}/{tag}_P.npy", P)
        json.dump(dict(stage=st, n=len(Fam), K=int(B.shape[0]), lp=lp, h=r['obj'], status=r['status'],
                       bound=r['bound'], t=time.time() - t0), open(f"{a.out}/{tag}.json", 'w'))


if __name__ == '__main__':
    main()
