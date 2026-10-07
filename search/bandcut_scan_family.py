#!/usr/bin/env python3
"""bandcut_scan_family.py -- the band family, defined for every `T`, and scanned
`T = 12 ... 4`   (2026-09-20, tasks/bandcut-scan).

The family, as Arslanov-Mustafin-Shangitbayev build it (EJC 28(4) 2021 #P4.22, sec 2):

  > "In Figure 2 we see one of the main ideas for packing unit squares: using of stacks (4,1)
  >  tilted by an angle alpha = arcsin(8/17).  The main idea for squeezing a packing follows
  >  from it: tilting stacks (4,1) by an angle alpha + eps so that the stack (4,1) is located
  >  in a vertical strip of width 4 - delta, where eps and delta are sufficiently small."

A **bar** = a `(p,1)` stack = `p` unit squares in a straight chain at a common tilt `t`,
consecutive centres differing by the unit edge vector.  Bounding box
`(p cos t + sin t) x (cos t + p sin t)`; the first factor is `<= p` exactly when
`tan(t/2) >= 1/p`, with equality at the window angle `t_p = 2 arctan(1/p)`
(`ARCH_TLEDGER.md` §4; `arcsin(8/17) = 2 arctan(1/4)` is the case `p = 4`).  Tilting past
`t_p` buys slack `p - (p cos t + sin t) > 0` in the `p` columns the bar crosses.

A **band** = bars at a common tilt, consecutive bars at perpendicular distance exactly 1
(edge to edge), with a free drift along the bar direction.  Bar lengths follow a profile
(the published ones are staircases, e.g. `2,3,4,4,3,2` in the `(8,4)` rectangle of Figure 6).
Bands come in two chiralities, `+t` (bars near-horizontal, slack in `x`) and `90 - t`
(bars near-vertical, slack in `y`); the published packings use one of each, meeting at an
elbow, so that every row *and* every column is cut.

Everything else is an axis-parallel square on (or near) the integer grid.  The whole
configuration is then handed to the disjunctive LP in the centres at fixed angles
(`S6_SKELETON.md` §3.1, via `bandcut_scan.polish`), which is exact for its assignment; what
comes back is a **feasible point**, i.e. a LOWER bound on `delta*_T`.

    python3 search/bandcut_scan_family.py --T 4 --deg 28:60:2 --profiles 4 22 33 232
    python3 search/bandcut_scan_family.py --T 12 --deg 28:34:1 --profiles 234432 --both
"""
import argparse
import itertools
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


# --------------------------------------------------------------------------- geometry
def make_band(t, profile, x0, y0, dx, mirror=False):
    """bars of lengths `profile` at tilt `t`, perpendicular spacing 1, drift `dx` per bar.

    `mirror=True` reflects the whole band in the line `y = x`, which sends the tilt `t` to
    `90 - t` and turns a slack-in-`x` band into a slack-in-`y` band."""
    c, s = math.cos(t), math.sin(t)
    dy = (1.0 + dx * s) / c                     # perpendicular spacing exactly 1
    out = []
    for r, p in enumerate(profile):
        bx, by = x0 + r * dx, y0 + r * dy
        for k in range(p):
            out.append((bx + k * c, by + k * s, t))
    if mirror:
        out = [(y, x, math.pi / 2 - a) for (x, y, a) in out]
    return out


def clear_map(X, Y, TH, T, ang, step):
    """clearance of an axis-parallel-or-`ang` unit square centred at each lattice point."""
    gs = np.arange(0.5, T - 0.5 + 1e-12, step)
    ca, sa = math.cos(ang), math.sin(ang)
    pa = 0.5 * (abs(ca) + abs(sa))
    GX, GY = np.meshgrid(gs, gs, indexing='ij')
    out = np.minimum.reduce([GX - pa, T - GX - pa, GY - pa, T - GY - pa])
    d = TH - ang
    m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
    for i in range(len(X)):
        dx, dy = GX - X[i], GY - Y[i]
        ci, si = math.cos(TH[i]), math.sin(TH[i])
        g = np.maximum.reduce([np.abs(ca * dx + sa * dy), np.abs(-sa * dx + ca * dy),
                               np.abs(ci * dx + si * dy), np.abs(-si * dx + ci * dy)]) - m[i]
        np.minimum(out, g, out=out)
    return gs, out


def fill_grid(X, Y, TH, T, n_need):
    """fast bulk fill: axis-parallel squares on the integer cells, best clearance first."""
    X, Y, TH = list(X), list(Y), list(TH)
    cells = [(ix + 0.5, iy + 0.5) for ix in range(T) for iy in range(T)]
    Xa, Ya, THa = np.array(X), np.array(Y), np.array(TH)
    scored = []
    for (cx, cy) in cells:
        if len(Xa):
            dx, dy = Xa - cx, Ya - cy
            ci, si = np.cos(THa), np.sin(THa)
            m = 0.5 + 0.5 * (np.abs(np.cos(THa)) + np.abs(np.sin(THa)))
            g = np.maximum.reduce([np.abs(dx), np.abs(dy),
                                   np.abs(ci * dx + si * dy),
                                   np.abs(-si * dx + ci * dy)]) - m
            v = float(g.min())
        else:
            v = 1e9
        scored.append((v, cx, cy))
    scored.sort(key=lambda z: -z[0])
    used = []
    for (v, cx, cy) in scored:
        if len(used) >= n_need:
            break
        if any(max(abs(cx - u), abs(cy - w)) < 1.0 - 1e-9 for (u, w) in used):
            continue
        used.append((cx, cy))
        X.append(cx); Y.append(cy); TH.append(0.0)
    return np.array(X), np.array(Y), np.array(TH)


def fill(X, Y, TH, T, n_need, step=0.05, angs=(0.0,)):
    """greedily add `n_need` further squares at the roomiest lattice points."""
    X, Y, TH = list(X), list(Y), list(TH)
    for _ in range(n_need):
        best = (-1e18, None)
        for a in angs:
            gs, cm = clear_map(np.array(X), np.array(Y), np.array(TH), T, a, step)
            k = int(np.argmax(cm))
            i, j = divmod(k, len(gs))
            if cm[i, j] > best[0]:
                best = (float(cm[i, j]), (float(gs[i]), float(gs[j]), a))
        if best[1] is None:
            break
        X.append(best[1][0]); Y.append(best[1][1]); TH.append(best[1][2])
    return np.array(X), np.array(Y), np.array(TH)


def worst_index(X, Y, TH, T):
    n = len(X)
    c, s = np.cos(TH), np.sin(TH)
    p = 0.5 * (np.abs(c) + np.abs(s))
    cl = np.minimum.reduce([X - p, T - X - p, Y - p, T - Y - p])
    for i in range(n):
        dx, dy = X - X[i], Y - Y[i]
        g = np.maximum.reduce([np.abs(c[i] * dx + s[i] * dy), np.abs(-s[i] * dx + c[i] * dy),
                               np.abs(c * dx + s * dy), np.abs(-s * dx + c * dy)])
        d = TH - TH[i]
        m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
        gg = g - m
        gg[i] = 1e9
        cl = np.minimum(cl, gg)
        cl[i] = min(cl[i], gg.min())
    return int(np.argmin(cl)), cl


def refine(X, Y, TH, T, rounds=4, cut=3.2, step=0.05, angs=(0.0,), lock=0):
    """LP polish, then relocate the tightest *unlocked* square, repeat.

    `lock` protects the first `lock` squares (the band) from relocation: without it the
    relocation loop simply dissolves the band back into the axis-parallel plateau, which is
    the right answer for `delta*` but tells us nothing about the band family."""
    v, X, Y = polish(X, Y, TH, T, cut=cut, iters=30, rounds=3)
    best = (v, X.copy(), Y.copy(), TH.copy())
    for _ in range(rounds):
        k, _cl = worst_index(X, Y, TH, T)
        if k < lock:
            cl = _cl.copy(); cl[:lock] = 1e9
            k = int(np.argmin(cl))
        X2 = np.delete(X, k); Y2 = np.delete(Y, k); TH2 = np.delete(TH, k)
        X2, Y2, TH2 = fill(X2, Y2, TH2, T, 1, step=step, angs=angs)
        v2, X2, Y2 = polish(X2, Y2, TH2, T, cut=cut, iters=30, rounds=3)
        if v2 > best[0] + 1e-13:
            best = (v2, X2.copy(), Y2.copy(), TH2.copy())
            X, Y, TH = X2, Y2, TH2
        else:
            break
    return best


# --------------------------------------------------------------------------- tight rows
def tight_report(X, Y, TH, T, delta, tol=1e-7):
    """which rows attain `delta` -- the binding constraints."""
    n = len(X)
    c, s = np.cos(TH), np.sin(TH)
    p = 0.5 * (np.abs(c) + np.abs(s))
    walls = []
    for i in range(n):
        for val, tag in ((X[i] - p[i], 'L'), (T - X[i] - p[i], 'R'),
                         (Y[i] - p[i], 'B'), (T - Y[i] - p[i], 'T')):
            if val <= delta + tol:
                walls.append((i, tag))
    pairs = []
    for i in range(n):
        dx, dy = X - X[i], Y - Y[i]
        g = np.maximum.reduce([np.abs(c[i] * dx + s[i] * dy), np.abs(-s[i] * dx + c[i] * dy),
                               np.abs(c * dx + s * dy), np.abs(-s * dx + c * dy)])
        d = TH - TH[i]
        m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
        gg = g - m
        for j in range(i + 1, n):
            if gg[j] <= delta + tol:
                pairs.append((i, j))
    return walls, pairs


def chain_report(X, Y, TH, T, delta, tol=1e-7):
    """look for a wall-to-wall chain among the tight rows: squares linked left-to-right
    (or bottom-to-top) by tight pair rows, anchored on both walls."""
    n = len(X)
    walls, pairs = tight_report(X, Y, TH, T, delta, tol)
    out = {}
    for axis, (lo, hi) in (('x', ('L', 'R')), ('y', ('B', 'T'))):
        C = X if axis == 'x' else Y
        adj = {i: [] for i in range(n)}
        for (i, j) in pairs:
            a, b = (i, j) if C[i] <= C[j] else (j, i)
            adj[a].append(b)
        src = set(i for (i, tg) in walls if tg == lo)
        snk = set(i for (i, tg) in walls if tg == hi)
        depth = {i: (1 if i in src else -10**9) for i in range(n)}
        for i in sorted(range(n), key=lambda k: C[k]):
            if depth[i] < 0:
                continue
            for j in adj[i]:
                depth[j] = max(depth[j], depth[i] + 1)
        best = max([depth[i] for i in snk], default=-1)
        out[axis] = best
    return out, len(walls), len(pairs)


# --------------------------------------------------------------------------- the scan
def specs(T, t, profiles, both, drifts, anchors):
    for prof_s in profiles:
        prof = [int(ch) for ch in prof_s]
        for dx in drifts:
            for (x0, y0) in anchors:
                b1 = make_band(t, prof, x0, y0, dx, mirror=False)
                yield ('one %s dx=%.2f @(%.2f,%.2f)' % (prof_s, dx, x0, y0), b1)
                if both:
                    for (x1, y1) in anchors:
                        b2 = make_band(t, prof, x1, y1, dx, mirror=True)
                        yield ('two %s dx=%.2f @(%.2f,%.2f)+(%.2f,%.2f)'
                               % (prof_s, dx, x0, y0, x1, y1), b1 + b2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--T', type=int, required=True)
    ap.add_argument('--deg', default='28:62:2', help='lo:hi:step in degrees')
    ap.add_argument('--profiles', nargs='+', default=['4'],
                    help='bar-length profiles, e.g. 4  22  234432')
    ap.add_argument('--both', action='store_true', help='also try a mirrored second band')
    ap.add_argument('--drifts', default='0,0.5,1.0')
    ap.add_argument('--step', type=float, default=0.05)
    ap.add_argument('--cut', type=float, default=3.2)
    ap.add_argument('--rounds', type=int, default=4)
    ap.add_argument('--anchor-step', type=float, default=0.5)
    ap.add_argument('--topk', type=int, default=8)
    ap.add_argument('--lock', action='store_true', help='protect the band from relocation')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    T = a.T
    n = T * T - T
    lo, hi, st = [float(x) for x in a.deg.split(':')]
    degs = [lo + k * st for k in range(int(round((hi - lo) / st)) + 1)]
    drifts = [float(x) for x in a.drifts.split(',')]
    anchors = [(x, y) for x in np.arange(0.3, T - 0.8, a.anchor_step)
               for y in np.arange(0.3, T - 0.8, a.anchor_step)]
    results = []
    print('# T = %d, n = %d, %d angles, %d profiles' % (T, n, len(degs), len(a.profiles)),
          flush=True)
    for dg in degs:
        t = math.radians(dg)
        bestd = None
        for (tag, tilted) in specs(T, t, a.profiles, a.both, drifts, anchors):
            if len(tilted) > n:
                continue
            X = np.array([q[0] for q in tilted]); Y = np.array([q[1] for q in tilted])
            TH = np.array([q[2] for q in tilted])
            if len(X) and (min(X.min(), Y.min()) < 0.2 or
                           max(X.max(), Y.max()) > T - 0.2):
                continue
            v0 = float_value(X, Y, TH, T)[0] if len(X) > 1 else 1.0
            if v0 < -1e-9:
                continue                                  # bars overlap: not a band
            nb = len(X)
            if n - nb > 24:
                X, Y, TH = fill_grid(X, Y, TH, T, n - nb)
            if len(X) < n:
                X, Y, TH = fill(X, Y, TH, T, n - len(X), step=a.step)
            if len(X) < n:
                continue
            v, X2, Y2, TH2 = refine(X, Y, TH, T, rounds=a.rounds, cut=a.cut, step=a.step,
                                    lock=(nb if a.lock else 0))
            results.append((v, dg, tag, X2, Y2, TH2))
            if bestd is None or v > bestd[0]:
                bestd = (v, tag)
        if bestd:
            print('  t = %7.3f deg   best delta = %+.9e   [%s]' % (dg, bestd[0], bestd[1]),
                  flush=True)
    results.sort(key=lambda z: -z[0])
    print()
    print('# top %d over the whole sweep at T = %d' % (a.topk, T))
    for (v, dg, tag, X2, Y2, TH2) in results[:a.topk]:
        ch, nw, npr = chain_report(X2, Y2, TH2, T, v)
        ntilt = int(np.sum(np.abs(TH2 % (math.pi / 2)) > 1e-6))
        print('  delta = %+.9e  t = %7.3f  k(near-axis) = %3d  tilted = %3d  '
              'longest tight chain x/y = %d/%d  tight walls %d pairs %d   [%s]'
              % (v, dg, n - ntilt, ntilt, ch['x'], ch['y'], nw, npr, tag), flush=True)
    if a.out and results:
        v, dg, tag, X2, Y2, TH2 = results[0]
        json.dump({'T': T, 'n': n, 'delta': v, 'tdeg': dg, 'tag': tag,
                   'sq': [[float(x), float(y), float(th)]
                          for x, y, th in zip(X2, Y2, TH2)]}, open(a.out, 'w'))
        print('wrote', a.out)


if __name__ == '__main__':
    main()
