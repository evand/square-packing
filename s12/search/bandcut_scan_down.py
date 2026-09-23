#!/usr/bin/env python3
"""bandcut_scan_down.py -- follow the tilted-band family of `s(T^2-T) < T` packings
DOWN from `T = 12` to `T = 4`   (2026-09-20, tasks/bandcut-scan).

Input: a real member of the family, read off a published packing with
`search/bandcut_scan_svg.py` (the `T = 12` and `T = 11` files parse to exactly `132` and
`110` squares and verify `delta > 0`).

One step `T -> T-1` is a **seam contraction**: pick a vertical grid line `a` and a horizontal
one `b`, delete every square whose centre lies in the unit strip `[a, a+1]` or `[b, b+1]`,
translate what is left across the seams, and re-close the container to `[0, T-1]^2`.  That is
the operation the family itself is built on: the big axis-parallel block loses one row and one
column, and the band has to be re-threaded.  The count is then brought back to exactly
`n' = (T-1)^2 - (T-1)` (drop the worst square, or add one at the roomiest free spot) and the
configuration is re-optimised by the disjunctive LP of `S6_SKELETON.md` §3.1 in the centres
with the angles held fixed (`bandcut_scan.polish`), with a relocation loop around it.

Everything reported is a feasible point -> a LOWER bound on `delta*_{T}`.

    python3 search/bandcut_scan_down.py --start p132.json --T 12 --steps 8 --tries 24
"""
import argparse
import json
import math
import os
import sys

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bandcut_scan import float_value, polish                          # noqa: E402


def load(path):
    d = json.load(open(path))
    sq = d['sq']
    X = np.array([s[0] for s in sq]); Y = np.array([s[1] for s in sq])
    TH = np.array([s[2] for s in sq])
    return X, Y, TH, d['s']


def clearances(X, Y, TH, T):
    """per-square: min over (its four wall slacks, its gaps to the other squares)."""
    n = len(X)
    c, s = np.cos(TH), np.sin(TH)
    p = 0.5 * (np.abs(c) + np.abs(s))
    out = np.minimum.reduce([X - p, T - X - p, Y - p, T - Y - p])
    for i in range(n):
        dx, dy = X - X[i], Y - Y[i]
        g = np.maximum.reduce([np.abs(c[i] * dx + s[i] * dy), np.abs(-s[i] * dx + c[i] * dy),
                               np.abs(c * dx + s * dy), np.abs(-s * dx + c * dy)])
        d = TH - TH[i]
        m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
        gg = g - m
        gg[i] = 1e9
        out[i] = min(out[i], gg.min())
        out = np.minimum(out, np.where(np.arange(n) == i, 1e9, gg))
    return out


def best_free_spot(X, Y, TH, T, angles, step=0.05):
    """the roomiest place for one more square, over a grid of centres x angles."""
    gs = np.arange(0.5, T - 0.5 + 1e-9, step)
    best = (-1e9, None)
    c, s = np.cos(TH), np.sin(TH)
    for a in angles:
        ca, sa = math.cos(a), math.sin(a)
        pa = 0.5 * (abs(ca) + abs(sa))
        d = TH - a
        m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
        for x in gs:
            if x - pa < -1e-12 or T - x - pa < -1e-12:
                continue
            dx = X - x
            for y in gs:
                if y - pa < -1e-12 or T - y - pa < -1e-12:
                    continue
                dy = Y - y
                g = np.maximum.reduce([np.abs(ca * dx + sa * dy), np.abs(-sa * dx + ca * dy),
                                       np.abs(c * dx + s * dy), np.abs(-s * dx + c * dy)]) - m
                v = min(g.min(), x - pa, T - x - pa, y - pa, T - y - pa)
                if v > best[0]:
                    best = (v, (x, y, a))
    return best


def contract(X, Y, TH, T, a, b):
    """delete the squares whose centres sit in the seams [a,a+1] x R and R x [b,b+1],
    then close up: container [0,T]^2 -> [0,T-1]^2."""
    keep = ~(((X > a) & (X < a + 1)) | ((Y > b) & (Y < b + 1)))
    X2, Y2, TH2 = X[keep].copy(), Y[keep].copy(), TH[keep].copy()
    X2 = np.where(X2 >= a + 1, X2 - 1.0, X2)
    Y2 = np.where(Y2 >= b + 1, Y2 - 1.0, Y2)
    return X2, Y2, TH2


def fix_count(X, Y, TH, T, n, rng, angles):
    while len(X) > n:
        cl = clearances(X, Y, TH, T)
        k = int(np.argmin(cl))
        X = np.delete(X, k); Y = np.delete(Y, k); TH = np.delete(TH, k)
    while len(X) < n:
        v, spot = best_free_spot(X, Y, TH, T, angles)
        if spot is None:
            return None
        X = np.append(X, spot[0]); Y = np.append(Y, spot[1]); TH = np.append(TH, spot[2])
    return X, Y, TH


def improve(X, Y, TH, T, rng, rounds=6, cut=3.2, angles=None, verbose=False):
    """alternate LP polish with relocation of the worst square."""
    best_v, bX, bY, bTH = float_value(X, Y, TH, T)[0], X.copy(), Y.copy(), TH.copy()
    v, X, Y = polish(X, Y, TH, T, cut=cut, iters=30, rounds=3)
    if v > best_v:
        best_v, bX, bY, bTH = v, X.copy(), Y.copy(), TH.copy()
    for r in range(rounds):
        cl = clearances(X, Y, TH, T)
        k = int(np.argmin(cl))
        X2 = np.delete(X, k); Y2 = np.delete(Y, k); TH2 = np.delete(TH, k)
        vv, spot = best_free_spot(X2, Y2, TH2, T, angles)
        if spot is None:
            break
        X2 = np.append(X2, spot[0]); Y2 = np.append(Y2, spot[1]); TH2 = np.append(TH2, spot[2])
        v2, X2, Y2 = polish(X2, Y2, TH2, T, cut=cut, iters=30, rounds=3)
        if v2 > best_v + 1e-12:
            best_v, bX, bY, bTH = v2, X2.copy(), Y2.copy(), TH2.copy()
            X, Y, TH = X2, Y2, TH2
        else:
            break
    return best_v, bX, bY, bTH


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', required=True)
    ap.add_argument('--T', type=int, required=True)
    ap.add_argument('--steps', type=int, default=8)
    ap.add_argument('--tries', type=int, default=16)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=None)
    ap.add_argument('--cut', type=float, default=3.2)
    a = ap.parse_args()

    X, Y, TH, s = load(a.start)
    T = float(a.T)
    lam = T / s
    X, Y = X * lam, Y * lam
    v0 = float_value(X, Y, TH, T)[0]
    print('start  T=%2d  n=%3d  delta = %+.9e' % (a.T, len(X), v0), flush=True)
    rng = np.random.default_rng(a.seed)
    chain = [(a.T, len(X), v0, X.copy(), Y.copy(), TH.copy())]

    for step in range(a.steps):
        Tc = chain[-1][0]
        if Tc <= 2:
            break
        Xc, Yc, THc = chain[-1][3], chain[-1][4], chain[-1][5]
        Tn = Tc - 1
        nn = Tn * Tn - Tn
        angles = sorted(set(np.round(THc % (math.pi / 2), 9))) + [0.0]
        angles = sorted(set([float(x) for x in angles]))
        best = None
        cuts = [(aa, bb) for aa in range(Tc - 1) for bb in range(Tc - 1)]
        rng.shuffle(cuts)
        for (aa, bb) in cuts[:a.tries]:
            X2, Y2, TH2 = contract(Xc, Yc, THc, Tc, float(aa), float(bb))
            r = fix_count(X2, Y2, TH2, float(Tn), nn, rng, angles)
            if r is None:
                continue
            X2, Y2, TH2 = r
            v, X3, Y3, TH3 = improve(X2, Y2, TH2, float(Tn), rng, angles=angles, cut=a.cut)
            if best is None or v > best[0]:
                best = (v, X3, Y3, TH3, aa, bb)
            print('    T=%2d cut(%d,%d)  delta = %+.6e   (best %+.6e)'
                  % (Tn, aa, bb, v, best[0]), flush=True)
        if best is None:
            print('  T=%2d : no candidate' % Tn); break
        v, X3, Y3, TH3, aa, bb = best
        ntilt = int(np.sum(np.abs(((TH3 + 1e-12) % (math.pi / 2))) > 1e-4))
        print('T=%2d  n=%3d  best delta = %+.9e   (seam %d,%d; %d tilted, %d near-axis)'
              % (Tn, len(X3), v, aa, bb, ntilt, len(X3) - ntilt), flush=True)
        chain.append((Tn, len(X3), v, X3, Y3, TH3))

    print()
    print('%-4s %-5s %-18s' % ('T', 'n', 'delta (feasible point)'))
    for (Tc, nc, vc, _, _, _) in chain:
        print('%-4d %-5d %+.9e' % (Tc, nc, vc))
    if a.out:
        json.dump({'chain': [{'T': c[0], 'n': c[1], 'delta': c[2],
                              'sq': [[float(x), float(y), float(t)]
                                     for x, y, t in zip(c[3], c[4], c[5])]} for c in chain]},
                  open(a.out, 'w'))
        print('wrote', a.out)


if __name__ == '__main__':
    main()
