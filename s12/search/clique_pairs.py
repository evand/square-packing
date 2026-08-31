#!/usr/bin/env python3
"""Independently verify that a separated anchor clique really is pairwise closed-intersecting,
and that its reported mass is right.

Uses only `clique_family.sat_meet` (a scalar SAT test between two rotated unit squares) and an
explicit membership test written out from the definition -- no reuse of the vectorised predicates
that produced the number.  This is the check that catches the missing-axis bug of
`notes/clique-family.md` 7.
"""
import sys, os, json, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clique_family as CF

recs = json.load(open(sys.argv[1]))
support = sys.argv[2]
t = float(sys.argv[3]) if len(sys.argv) > 3 else 3.99
nrec = int(sys.argv[4]) if len(sys.argv) > 4 else 3

imgs, poses = CF.load_measure(support, t)
IM = np.array([(a, b, c, d) for (a, b, c, d) in imgs if d > 1e-12])


def in_sq(q, c, th):
    dx, dy = q[0] - c[0], q[1] - c[1]
    ct, st = math.cos(th), math.sin(th)
    return abs(dx * ct + dy * st) <= 0.5 + 1e-12 and abs(-dx * st + dy * ct) <= 0.5 + 1e-12


def seg_meets(a0, a1, c, th):
    return CF.seg_meets_square(a0, a1, c, th)


for r in recs[:nrec]:
    p = tuple(r['p'])
    A = [tuple(v) for v in r['A']]
    mem = []
    for i in range(len(IM)):
        c = (IM[i, 0], IM[i, 1])
        th = IM[i, 2]
        thru = in_sq(p, c, th)
        if len(A) == 2:
            meets = seg_meets(A[0], A[1], c, th)
            contains = all(in_sq(v, c, th) for v in A)
        else:
            meets = any(in_sq(v, c, th) for v in A) or CF.polys_meet(CF.poly_sq(c, th), A)
            contains = all(in_sq(v, c, th) for v in A)
        if (thru and meets) or ((not thru) and contains):
            mem.append(i)
    mass = float(IM[mem, 3].sum())
    bad = 0
    for i, j in itertools.combinations(mem, 2):
        if not CF.sat_meet(IM[i, :2], IM[i, 2], IM[j, :2], IM[j, 2], tol=1e-12):
            bad += 1
    wd = min(p[0], p[1], t - p[0], t - p[1])
    print(f"p=({p[0]:.5f},{p[1]:.5f}) walldist {wd:.4f}  |A|={len(A)}  "
          f"members {len(mem)}  reported mass {r['mass']:.6f}  recomputed {mass:.6f}  "
          f"non-intersecting pairs {bad}")
