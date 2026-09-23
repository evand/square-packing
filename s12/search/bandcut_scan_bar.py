#!/usr/bin/env python3
"""bandcut_scan_bar.py -- the band family with the band forced to be *wall-to-wall*,
scanned `T = 12 ... 4`   (2026-09-20, tasks/bandcut-scan).

Why this and not `bandcut_scan_family.py`: with the band merely present, the LP in the
centres slides it into a harmless corner and reports the axis-parallel plateau value
(`delta = 0` at `T = 4`).  What makes the band family a *counterexample* family is that the
band is threaded across the container so that no wall-to-wall chain of `T` near-axis squares
survives (`ARCH_TLEDGER.md` §5.4: `A4` is the lemma that fails at `T = 17`).  So here the bar
is anchored: its first square is on the `lo` wall and its last on the `hi` wall, both at slack
exactly `delta` -- the two wall rows of the chain bound of `ARCH_TLEDGER.md` §2.1, imposed as
equalities in the LP.

    band(p, t)      p unit squares at a common tilt t in a straight chain, anchored in x
    cross(p,q,t)    that, plus a q-chain at 90 - t anchored in y  (cuts rows AND columns)

Everything else is axis-parallel and free.  Output is a feasible point, i.e. a lower bound
on `delta*_T` restricted to the family.

    python3 search/bandcut_scan_bar.py --T 4 --deg 24:66:1 --p 2,3,4,5 --cross
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
from bandcut_scan import float_value, polish                            # noqa: E402
from bandcut_scan_family import fill_grid, fill, worst_index, chain_report   # noqa: E402


def bar_x(p, t, T, y0):
    """p squares at tilt t in a straight chain spanning the container in x, centred on y0."""
    c, s = math.cos(t), math.sin(t)
    half = 0.5 * (c + s)
    gap = (T - 2 * half - (p - 1) * c) / max(p - 1, 1) if p > 1 else 0.0
    step = c + max(gap, 0.0) * c / max(c, 1e-9) if False else None
    # spread the p centres evenly between the two wall-touching extremes along e
    x0, x1 = half, T - half
    if p == 1:
        return [(0.5 * T, y0, t)]
    d = (x1 - x0) / (p - 1)
    out = []
    for k in range(p):
        out.append((x0 + k * d, y0 + (k - (p - 1) / 2.0) * d * (s / c), t))
    return out


def bar_y(q, t, T, x0):
    tt = math.pi / 2 - t
    return [(y, x, math.pi / 2 - a) for (x, y, a) in bar_x(q, t, T, x0)]


def refine_anchored(X, Y, TH, T, anchors, lock, rounds=3, cut=3.2, step=0.1):
    v, X, Y = polish(X, Y, TH, T, cut=cut, iters=30, rounds=3, anchors=anchors)
    best = (v, X.copy(), Y.copy(), TH.copy())
    for _ in range(rounds):
        k, cl = worst_index(X, Y, TH, T)
        if k < lock:
            cl2 = cl.copy(); cl2[:lock] = 1e9
            k = int(np.argmin(cl2))
        X2 = np.delete(X, k); Y2 = np.delete(Y, k); TH2 = np.delete(TH, k)
        X2, Y2, TH2 = fill(X2, Y2, TH2, T, 1, step=step)
        v2, X2, Y2 = polish(X2, Y2, TH2, T, cut=cut, iters=30, rounds=3, anchors=anchors)
        if v2 > best[0] + 1e-13:
            best = (v2, X2.copy(), Y2.copy(), TH2.copy())
            X, Y, TH = X2, Y2, TH2
        else:
            break
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--T', type=int, required=True)
    ap.add_argument('--deg', default='24:66:2')
    ap.add_argument('--p', default='2,3,4,5')
    ap.add_argument('--cross', action='store_true')
    ap.add_argument('--offstep', type=float, default=0.5)
    ap.add_argument('--step', type=float, default=0.1)
    ap.add_argument('--cut', type=float, default=3.2)
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--topk', type=int, default=10)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    T, n = a.T, a.T * a.T - a.T
    lo, hi, st = [float(x) for x in a.deg.split(':')]
    degs = [lo + k * st for k in range(int(round((hi - lo) / st)) + 1)]
    ps = [int(x) for x in a.p.split(',')]
    offs = list(np.arange(0.8, T - 0.7, a.offstep))
    res = []
    print('# T = %d  n = %d   bar lengths %s   %d angles   cross=%s'
          % (T, n, ps, len(degs), a.cross), flush=True)
    for dg in degs:
        t = math.radians(dg)
        bestd = None
        for p in ps:
            c, s = math.cos(t), math.sin(t)
            if p * c + s > T + 1e-12:            # the bar cannot span: too long for the wall
                continue
            for y0 in offs:
                b = bar_x(p, t, T, y0)
                sets = [('bar p=%d y0=%.1f' % (p, y0), b, [(0, p - 1, 'x')])]
                if a.cross:
                    for q in ps:
                        if q * c + s > T + 1e-12:
                            continue
                        for x0 in offs:
                            b2 = bar_y(q, t, T, x0)
                            sets.append(('cross p=%d q=%d y0=%.1f x0=%.1f' % (p, q, y0, x0),
                                         b + b2, [(0, p - 1, 'x'), (p, p + q - 1, 'y')]))
                for (tag, tilted, anch) in sets:
                    nb = len(tilted)
                    if nb > n:
                        continue
                    X = np.array([z[0] for z in tilted]); Y = np.array([z[1] for z in tilted])
                    TH = np.array([z[2] for z in tilted])
                    if nb > 1 and float_value(X, Y, TH, T)[0] < -1e-9:
                        continue
                    X, Y, TH = fill_grid(X, Y, TH, T, n - nb)
                    if len(X) < n:
                        X, Y, TH = fill(X, Y, TH, T, n - len(X), step=a.step)
                    if len(X) < n:
                        continue
                    v, X2, Y2, TH2 = refine_anchored(X, Y, TH, T, anch, nb,
                                                     rounds=a.rounds, cut=a.cut, step=a.step)
                    res.append((v, dg, tag, nb, X2, Y2, TH2))
                    if bestd is None or v > bestd[0]:
                        bestd = (v, tag)
        if bestd:
            print('  t = %7.3f   best delta = %+.9e   [%s]' % (dg, bestd[0], bestd[1]),
                  flush=True)
    res.sort(key=lambda z: -z[0])
    print()
    print('# top %d at T = %d' % (a.topk, T))
    for (v, dg, tag, nb, X2, Y2, TH2) in res[:a.topk]:
        ch, nw, npr = chain_report(X2, Y2, TH2, T, v)
        print('  delta = %+.9e  t = %7.3f  band = %2d  k = %3d  tight chain x/y = %d/%d  [%s]'
              % (v, dg, nb, n - nb, ch['x'], ch['y'], tag), flush=True)
    if a.out and res:
        v, dg, tag, nb, X2, Y2, TH2 = res[0]
        json.dump({'T': T, 'n': n, 'delta': v, 'tdeg': dg, 'tag': tag, 'band': nb, 's': T,
                   'sq': [[float(x), float(y), float(th)]
                          for x, y, th in zip(X2, Y2, TH2)]}, open(a.out, 'w'))
        print('wrote', a.out)


if __name__ == '__main__':
    main()
