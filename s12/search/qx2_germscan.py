#!/usr/bin/env python3
"""qx2_germscan.py -- exact scan of the tile-germ blow-up poses of a box cover (task quadrant-exact-w2).

At a centre c0 on the grid where lines lie on the edges of the axis-parallel square, the mass at (c0 + u t, u) for
u -> 0 converges to a piecewise-linear function of t (the germ blow-up, LINE_COVER.md sec 2).  This evaluates the
EXACT mass (zm_mixed.exact_mass, Fractions) at u = 10^-12 and both tilt directions (theta and 90deg - theta, the latter
as the diagonal image), for c0 on the 1/10-grid of [1/2, 7/2]^2 and t on a grid of [-T, T]^2.  A diagnostic, not a
certificate: it reports the minimum and every pose below 1.
usage: python3 qx2_germscan.py COVER [--pitch 1/10] [--tmax 6/5] [--tstep 1/10] [--hi 7/2] [--nproc 4]
(--hi: upper end of the centre grid; 7/2 for the R = 2 box, R + 3/2 in general)
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import multiprocessing as mp
import zm_mixed as ZM, mixed_cover as MC, zeromargin as zm

_COV = None


def work(args):
    c0x, c0y, ts, u = args
    c, s = zm.trig(u); hw = (c + s) / 2
    out = []
    best = None
    for tx in ts:
        for ty in ts:
            for (X, Y) in ((c0x + u * tx, c0y + u * ty),):
                if X < hw or Y < hw or X > _COV.m - hw or Y > _COV.m - hw: continue
                v = ZM.exact_mass(_COV, X, Y, u)
                if best is None or v < best[0]: best = (v, X, Y)
                if v < 1: out.append((v, X, Y))
    return (c0x, c0y), best, out


def main():
    global _COV
    ap = argparse.ArgumentParser(); ap.add_argument('cover')
    ap.add_argument('--pitch', default='1/10'); ap.add_argument('--tmax', default='6/5'); ap.add_argument('--tstep', default='1/10')
    ap.add_argument('--u', default='1/1000000000000'); ap.add_argument('--hi', default='7/2'); ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    _COV = ZM.Cover(MC.load(a.cover))
    p = F(a.pitch); T = F(a.tmax); dt = F(a.tstep); u = F(a.u)
    n = int(T / dt)
    ts = [k * dt for k in range(-n, n + 1)]
    g = [F(1, 2) + k * p for k in range(int((F(a.hi) - F(1, 2)) / p) + 1)]
    jobs = []
    for sgn in (1, -1):
        for x in g:
            for y in g:
                # sgn = -1: the tilt 90deg - theta, i.e. the diagonal image (y, x) at +theta (the cover is D4-invariant)
                jobs.append((x, y, ts, u) if sgn == 1 else (y, x, ts, u))
    with mp.get_context('fork').Pool(a.nproc) as pool:
        res = pool.map(work, jobs, chunksize=4)
    best = min((r[1] for r in res if r[1] is not None), key=lambda t: t[0])
    bad = [b for r in res for b in r[2]]
    print(f"germ scan {a.cover}: {len(jobs)} centres x {len(ts)**2} offsets at u = {u}: min exact mass {float(best[0]):.12f} "
          f"at ({float(best[1]):.12f}, {float(best[2]):.12f}); poses below 1: {len(bad)}")
    for v, X, Y in sorted(bad)[:20]:
        print(f"   {float(v):.9f} at cx = {float(X):.13f}, cy = {float(Y):.13f}")


if __name__ == '__main__':
    main()
