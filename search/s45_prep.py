#!/usr/bin/env python3
"""s45_prep.py -- inputs for the m = 7 cover LP (search/S45_COVER.md; also usable at any integer m).

  transplant IN.txt M OUT.txt
      stretch a D4-symmetric certificate-format cover of [0,m0]^2 (m0 odd, m0 >= 5) to [0,M]^2 by
      replicating its centre column/row of unit cells: a coordinate X in [0, m0] (units of 1/D) maps to
      X if X <= h0*D, to X + (M - m0)*D if X >= (h0+1)*D (h0 = (m0-1)/2), and a point strictly inside the
      centre strip (h0, h0+1) goes to every copy X + k*D, k = 0..M-m0.  Points ON the lines x = h0,
      h0+1 map to themselves / their shifted copy; the lines in between (h0+1 .. h0+M-m0) get a copy
      of the line x = h0+1.  Corner and wall-band cells are then transplanted unchanged, the centre
      cell is tiled over the new interior.  Weights are unchanged (HEURISTIC warm start: not a cover).
  neartile M OUT.txt [--rad 0.012 --step 0.003 --th0 0.02 --th1 1.5 --dth 0.02]
      the near-tile pose pool of S21_COVER.md sec 1 for [0,M]^2: the M^2 tile centres +- rad
      (step `step`) x theta = th0..th1 deg (step dth), clipped into the admissible box.  With D4-
      symmetric covers one cell per C4 orbit is written (C4 fixes theta; --all for every tile).
"""
import sys, math, argparse
import numpy as np


def transplant(inp, M, out):
    t = open(inp).read().split(); it = iter(t)
    sn, sd, D, WD, n = [int(next(it)) for _ in range(5)]
    assert sd == 1 and sn % 2 == 1 and sn >= 5, "needs an odd integer container"
    m0 = sn; h0 = (m0 - 1) // 2; ext = M - m0; assert ext >= 0 and ext % 2 == 0
    def mapc(X):
        if X <= h0 * D: return [X]
        if X >= (h0 + 1) * D:
            return [X + k * D for k in range(ext + 1)] if X == (h0 + 1) * D else [X + ext * D]
        return [X + k * D for k in range(ext + 1)]
    pts = {}
    for _ in range(n):
        X, Y, W = int(next(it)), int(next(it)), int(next(it))
        for a in mapc(X):
            for b in mapc(Y):
                pts[(a, b)] = pts.get((a, b), 0) + W
    with open(out, 'w') as f:
        f.write(f"{M} 1\n{D}\n{WD}\n{len(pts)}\n")
        for (a, b), W in sorted(pts.items()): f.write(f"{a} {b} {W}\n")
    tot = sum(pts.values()) / WD
    print(f"transplanted {inp} (s={m0}, {n} points) -> {out} (s={M}, {len(pts)} points, total {tot:.6f})")


def neartile(M, out, rad=0.012, step=0.003, th0=0.02, th1=1.5, dth=0.02, allcells=False):
    offs = np.arange(-rad, rad + 1e-12, step)
    ths = np.arange(th0, th1 + 1e-9, dth)
    n = 0
    with open(out, 'w') as f:
        f.write(f"# cx cy theta_rad -- near-tile poses: tile centres +-{rad} (step {step}), theta {th0}..{th1} deg "
                f"(step {dth}); s45_prep.py, S45_COVER.md\n")
        for i in range(M):
            for j in range(M):
                rot = [(i, j), (M - 1 - j, i), (M - 1 - i, M - 1 - j), (j, M - 1 - i)]
                if not allcells and (i, j) != min(rot): continue      # one cell per C4 orbit (C4 fixes theta)
                for d in ths:
                    th = math.radians(d); w = (math.cos(th) + math.sin(th)) / 2
                    lo, hi = w + 1e-7, M - w - 1e-7
                    for dx in offs:
                        for dy in offs:
                            cx = min(max(i + 0.5 + dx, lo), hi); cy = min(max(j + 0.5 + dy, lo), hi)
                            f.write("%.17g %.17g %.17g\n" % (cx, cy, th)); n += 1
    print(f"{n} near-tile poses written to {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['transplant', 'neartile'])
    ap.add_argument('args', nargs='+')
    ap.add_argument('--rad', type=float, default=0.012); ap.add_argument('--step', type=float, default=0.003)
    ap.add_argument('--th0', type=float, default=0.02); ap.add_argument('--th1', type=float, default=1.5)
    ap.add_argument('--dth', type=float, default=0.02); ap.add_argument('--all', action='store_true')
    a = ap.parse_args()
    if a.mode == 'transplant': transplant(a.args[0], int(a.args[1]), a.args[2])
    else: neartile(int(a.args[0]), a.args[1], a.rad, a.step, a.th0, a.th1, a.dth, a.all)


if __name__ == '__main__':
    main()
