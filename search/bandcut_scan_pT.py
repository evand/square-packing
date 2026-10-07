#!/usr/bin/env python3
"""bandcut_scan_pT.py -- the sharp end of the band family: one **full-width `(p,1)` stack**,
scanned in the tilt, at every `T`   (2026-09-20, tasks/bandcut-scan).

The whole mechanism of `s(T^2-T) < T` is one inequality (`ARCH_TLEDGER.md` §4, and
Arslanov-Mustafin-Shangitbayev §2 in the same words): a stack of `p` unit squares at a common
tilt `t` spans `p cos t + sin t`, which is `< p` as soon as `t > t_p = 2 arctan(1/p)`.  Lay
such a stack wall-to-wall across the container and every row it crosses gains slack; the price
is its height `cos t + p sin t`, which the `T` free cells of `n = T^2 - T` have to pay for.

This script places exactly that -- `p` squares at tilt `t`, in a straight chain, with the
first anchored on the left wall and the last on the right wall (equality rows in the LP, so
the band cannot slide away) -- fills the rest with axis-parallel squares, and maximises the
margin over the centres with the disjunctive LP of `S6_SKELETON.md` §3.1.  A feasible point,
so a LOWER bound on `delta*_T`.

    python3 search/bandcut_scan_pT.py --T 4 --p 4 --fine
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
from bandcut_scan import float_value                                    # noqa: E402
from bandcut_scan_family import fill_grid, fill                         # noqa: E402
from bandcut_scan_bar import bar_x, refine_anchored                     # noqa: E402


def run(T, p, tdegs, y0s, step, cut, rounds):
    n = T * T - T
    best = (-9.0, None, None, None, None)
    for tdeg in tdegs:
        t = math.radians(tdeg)
        c, s = math.cos(t), math.sin(t)
        if p * c + s > T + 1e-12:
            continue
        for y0 in y0s:
            b = bar_x(p, t, T, y0)
            X = np.array([z[0] for z in b]); Y = np.array([z[1] for z in b])
            TH = np.array([z[2] for z in b])
            if p > 1 and float_value(X, Y, TH, T)[0] < -1e-9:
                continue
            X, Y, TH = fill_grid(X, Y, TH, T, n - p)
            if len(X) < n:
                X, Y, TH = fill(X, Y, TH, T, n - len(X), step=step)
            if len(X) < n:
                continue
            v, X2, Y2, TH2 = refine_anchored(X, Y, TH, T, [(0, p - 1, 'x')], p,
                                             rounds=rounds, cut=cut, step=step)
            if v > best[0]:
                best = (v, tdeg, y0, (X2, Y2, TH2), None)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--T', type=int, required=True)
    ap.add_argument('--p', type=int, default=None, help='stack length (default T, T-1, T-2)')
    ap.add_argument('--step', type=float, default=0.1)
    ap.add_argument('--cut', type=float, default=3.2)
    ap.add_argument('--rounds', type=int, default=3)
    ap.add_argument('--coarse', type=float, default=0.5, help='tilt step, degrees')
    ap.add_argument('--ystep', type=float, default=0.2)
    ap.add_argument('--fine', action='store_true', help='second pass, 0.01 deg / 0.02 offset')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    T = a.T
    ps = [a.p] if a.p else [T, T - 1, T - 2, T - 3]
    print('# T = %d  n = %d   window angles t_p = 2 arctan(1/p):' % (T, T * T - T),
          {p: round(math.degrees(2 * math.atan(1.0 / p)), 4) for p in ps if p > 1}, flush=True)
    overall = None
    for p in ps:
        if p < 2:
            continue
        t0 = math.degrees(2 * math.atan(1.0 / p))
        tds = list(np.arange(t0 + 0.02, min(t0 + 24.0, 70.0), a.coarse))
        ys = list(np.arange(0.6, T - 0.5, a.ystep))
        v, tdeg, y0, cfg, _ = run(T, p, tds, ys, a.step, a.cut, a.rounds)
        if cfg is None:
            print('  p = %2d : no wall-to-wall stack fits' % p, flush=True)
            continue
        if a.fine:
            tds = list(np.arange(tdeg - a.coarse, tdeg + a.coarse, 0.02))
            ys = list(np.arange(max(0.6, y0 - a.ystep), min(T - 0.5, y0 + a.ystep), 0.02))
            v2, td2, y2, cfg2, _ = run(T, p, tds, ys, a.step, a.cut, a.rounds)
            if v2 > v:
                v, tdeg, y0, cfg = v2, td2, y2, cfg2
        X2, Y2, TH2 = cfg
        val, kp, wall = float_value(X2, Y2, TH2, T)
        ntilt = int(np.sum(np.abs(TH2 % (math.pi / 2)) > 1e-9))
        print('  p = %2d : best delta = %+.9e   t = %8.4f deg (t_p = %8.4f, excess %+.4f)  '
              'y0 = %.2f  tilted = %d  k = %d  binding pair %s  wall %.3e'
              % (p, v, tdeg, t0, tdeg - t0, y0, ntilt, len(X2) - ntilt, kp, wall), flush=True)
        if overall is None or v > overall[0]:
            overall = (v, p, tdeg, X2, Y2, TH2)
    if overall:
        print('T = %2d   best over the full-width-stack family : delta = %+.9e  (p = %d, t = %.4f)'
              % (T, overall[0], overall[1], overall[2]), flush=True)
        if a.out:
            json.dump({'T': T, 'n': T * T - T, 'delta': overall[0], 'p': overall[1],
                       'tdeg': overall[2], 's': T,
                       'sq': [[float(x), float(y), float(t)] for x, y, t
                              in zip(overall[3], overall[4], overall[5])]}, open(a.out, 'w'))
            print('wrote', a.out)


if __name__ == '__main__':
    main()
