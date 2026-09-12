#!/usr/bin/env python3
"""rung2_bound.py -- the m^2 obstruction of search/RUNG2.md, and its check on a certificate.

A *monotone witness certificate* for a pose box B is a set S of cover points with total weight
>= 1 such that EVERY p in S lies in the closed unit square Q(c, theta) at EVERY admissible pose
(c, theta) of B.  CORE, P1 and the ADM primitive of search/zeromargin.py all produce exactly this
kind of certificate, and so does any union of them; the triangle lemma TRI does NOT (it is
disjunctive: it says *some* vertex is captured, not that a fixed one is).

THEOREM (search/RUNG2.md sec 2).  Fix a tile (i, j) of [0,m]^2 and signs sx, sy with sx = +1 if
i < m-1 (the centre may move right) and sx = -1 if i = m-1 (it may only move left), likewise sy.
The pose (i+1/2, j+1/2, 0) is admissible and, the leaves being closed, finitely many and space-
filling, one leaf box B contains it together with poses (i+1/2 + sx e, j+1/2 + sy e, theta) for
arbitrarily small e, theta > 0.  B's monotone witness set S therefore lies in

    PIN(i, j, sx, sy) = intersection of the closed unit squares over that whole family,

so w(PIN(i, j, sx, sy)) >= 1 for every tile and every admissible sign pair.  Minimising the total
weight subject to those constraints is an LP; on a 1/20 or 1/40 lattice of candidate points its
value is exactly m^2 = 16, and the LP value on a lattice is an upper bound for the continuum one,
so the computed 16 is the bound for lattice covers (every cover this project builds is one: the
columns of search/closed4.py live on a 1/1000 lattice).

So no cover of [0,4]^2 with total weight below 16 -- in particular none with W < 13 -- can be
certified by CORE / P1 / ADM and their unions alone, at any subdivision depth.  Rung 2 needs a
DISJUNCTIVE primitive (the weighted analogue of TRI / the two-region argument).

Modes
  bound  [--m 4] [--pitch 20] [--octants]   the LP behind the theorem, on a lattice: with the
         one-sided (octant) constraints the value is exactly m^2; without them (i.e. if only the
         *two-sided* information at each tile pose is used) it drops to m^2/2, which is why the
         one-sided reading matters.
  check  CERT                               the direct test on a certificate-format cover: the
         exact (Fraction) weight of PIN(i,j,sx,sy) for every tile and sign pair, and the list of
         those below 1 -- each of which is a pose the box checker can never certify.

usage: python3 search/rung2_bound.py bound [--m 4] [--pitch 20] [--octants]
       python3 search/rung2_bound.py check runs/closed4_best_x103.txt
"""
import argparse, math, sys
from fractions import Fraction as F
import numpy as np


def pin_mask(X, Y, m, i, j, thetas, sx=0, sy=0, tol=1e-12):
    """Points captured at every pose of the admissible curve through the tile-(i,j) pose
    (i+1/2, j+1/2, 0).  sx/sy != 0 additionally move the centre in that direction, i.e. restrict
    to the poses of a leaf box that touches the tile pose from one side."""
    ok = np.ones(len(X), dtype=bool)
    for th in thetas:
        w = math.cos(th) + math.sin(th)
        for dcx in ((0.0,) if sx == 0 else (0.0, sx * 1e-3)):
            for dcy in ((0.0,) if sy == 0 else (0.0, sy * 1e-3)):
                cx = min(max(i + 0.5 + dcx, w / 2), m - w / 2)
                cy = min(max(j + 0.5 + dcy, w / 2), m - w / 2)
                c, s = math.cos(th), math.sin(th)
                dx, dy = X - cx, Y - cy
                ok &= (np.abs(dx * c + dy * s) <= 0.5 + tol) & (np.abs(-dx * s + dy * c) <= 0.5 + tol)
    return ok


def bound(a):
    from scipy.optimize import linprog
    m, K = a.m, a.pitch
    g = np.arange(0, m * K + 1) / K
    GX, GY = np.meshgrid(g, g, indexing='ij')
    X, Y = GX.ravel(), GY.ravel()
    thetas = [0.0, 1e-4, 3e-4, 1e-3, 3e-3]
    rows, labels = [], []
    for i in range(m):
        for j in range(m):
            sxs = [0] if not a.octants else ([1] if i == 0 else [-1] if i == m - 1 else [1, -1])
            sys_ = [0] if not a.octants else ([1] if j == 0 else [-1] if j == m - 1 else [1, -1])
            for sx in sxs:
                for sy in sys_:
                    msk = pin_mask(X, Y, m, i, j, thetas, sx, sy)
                    rows.append(msk.astype(float)); labels.append((i, j, sx, sy, int(msk.sum())))
    A = np.array(rows)
    res = linprog(c=np.ones(len(X)), A_ub=-A, b_ub=-np.ones(len(rows)), bounds=(0, None), method='highs')
    print(f"container [0,{m}]^2, lattice pitch 1/{K} ({len(X)} candidate points), "
          f"{len(rows)} PIN constraints, octants={a.octants}")
    print(f"LP value (lower bound on any monotone-witness-certifiable cover on this lattice) = {res.fun:.6f}")
    w = res.x; sup = np.nonzero(w > 1e-9)[0]
    print(f"support: {len(sup)} points, total {w.sum():.6f}")
    for k in sup[np.argsort(-w[sup])][:20]:
        print(f"   ({X[k]:.3f}, {Y[k]:.3f}) w={w[k]:.6f}")


def check(path):
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    ms = F(sn, sd)
    assert ms.denominator == 1, "container side must be an integer for the tile test"
    m = int(ms)
    Xi = np.empty(n, dtype=np.int64); Yi = np.empty(n, dtype=np.int64); Wi = [0] * n
    for k in range(n):
        x, y, w = map(int, tok[5 + 3 * k: 8 + 3 * k])
        Xi[k] = x; Yi[k] = y; Wi[k] = w
    tot = F(sum(Wi), W)
    Xf = Xi / D; Yf = Yi / D
    print(f"container [0,{m}]^2, {n} points, total weight {tot} = {float(tot):.9f}")
    print("w(PIN(i,j,sx,sy)) for every tile and every admissible sign pair -- each must be >= 1")
    print("for a monotone witness certificate to exist at the pose (i+1/2, j+1/2, 0):")
    thetas = [0.0, 1e-5, 1e-4, 1e-3]
    bad = []
    for j in range(m - 1, -1, -1):
        cells = []
        for i in range(m):
            sxs = [1] if i == 0 else [-1] if i == m - 1 else [1, -1]
            sys_ = [1] if j == 0 else [-1] if j == m - 1 else [1, -1]
            worst = None
            for sx in sxs:
                for sy in sys_:
                    msk = pin_mask(Xf, Yf, m, i, j, thetas, sx, sy)
                    v = F(sum(Wi[k] for k in np.nonzero(msk)[0]), W)
                    if worst is None or v < worst[0]: worst = (v, sx, sy)
                    if v < 1: bad.append((i, j, sx, sy, v))
            cells.append(worst[0])
        print("  j=%d: " % j + "  ".join(f"{float(v):8.5f}" for v in cells))
    print(f"\nPIN sets below 1: {len(bad)}")
    for i, j, sx, sy, v in bad:
        print(f"   tile ({i},{j}) sx={sx:+d} sy={sy:+d}: weight {v} = {float(v):.6f}  -> the pose "
              f"({i + 0.5}, {j + 0.5}, 0) can never be certified by a monotone witness set")
    print("VERDICT:", "no monotone-witness obstruction found at the tile poses" if not bad else
          "OBSTRUCTED -- this cover cannot be certified by CORE/P1/ADM at any depth")


def dual(a):
    """The explicit dual certificate of Theorem 1, verified on a complete set of cell
    representatives.  PIN-membership depends only on the position of a point relative to the
    half-integer grid, so it is constant on every cell of that arrangement, and the 1/4-lattice
    contains a representative of every cell (2-cells at (k+1/2)/2, 1-cells at their midpoints,
    0-cells themselves).  Checking `multiplicity <= 1` on the 1/4-lattice therefore checks it on
    ALL of [0,m]^2 -- the bound that follows holds for every cover, not only lattice ones."""
    m = a.m
    sx = [1, 1] + [-1] * (m - 2)       # the sign pattern of Theorem 1: the +/- switch must sit
                                       # between two INTERIOR columns (Lemma 2, case (iv))
    sy = list(sx)
    if a.switch is not None:
        sx = [1] * (a.switch + 1) + [-1] * (m - a.switch - 1); sy = list(sx)
    g = np.arange(0, m * 4 + 1) / 4
    GX, GY = np.meshgrid(g, g, indexing='ij')
    X, Y = GX.ravel(), GY.ravel()
    thetas = [0.0, 1e-5, 1e-4, 1e-3]
    mult = np.zeros(len(X), dtype=int)
    for i in range(m):
        for j in range(m):
            mult += pin_mask(X, Y, m, i, j, thetas, sx[i], sy[j]).astype(int)
    print(f"container [0,{m}]^2; dual certificate lambda = 1 on the {m*m} constraints")
    print(f"  sigma_x per column: {sx}")
    print(f"  sigma_y per row:    {sy}")
    print(f"sum of lambda = {m*m}")
    print(f"max multiplicity over the 1/4-lattice ({len(X)} cell representatives) = {mult.max()}")
    bad = np.nonzero(mult > 1)[0]
    for k in bad[:10]: print(f"   OVERLAP at ({X[k]}, {Y[k]}): multiplicity {mult[k]}")
    print("VERDICT:", f"the {m*m} PIN sets are pairwise disjoint  =>  every monotone-witness-"
          f"certifiable cover of [0,{m}]^2 has total weight >= {m*m}" if mult.max() <= 1 else
          "NOT disjoint for this sign pattern")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['bound', 'check', 'dual'])
    ap.add_argument('cert', nargs='?')
    ap.add_argument('--m', type=int, default=4)
    ap.add_argument('--pitch', type=int, default=20)
    ap.add_argument('--octants', action='store_true')
    ap.add_argument('--switch', type=int, default=None,
                    help='dual: the last column with sigma_x = +1 (default 1: the switch sits '
                         'between columns 1 and 2, both interior, which Lemma 2 requires)')
    a = ap.parse_args()
    if a.mode == 'bound': bound(a)
    elif a.mode == 'dual': dual(a)
    else: check(a.cert)


if __name__ == '__main__':
    main()
