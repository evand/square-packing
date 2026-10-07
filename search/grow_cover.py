#!/usr/bin/env python3
"""grow_cover.py -- warm starts for the m -> m+1 rung of the cover LP (search/S32_COVER.md).  HEURISTIC inputs only.

  transplant IN OUT
      grow a D4-symmetric certificate-format cover of [0,m]^2 (m odd or even) into one of [0,m+1]^2 by
      duplicating the middle cell column and row: a coordinate x of the old cover goes to
      {x : x < c+1} u {x + 1 : x >= c}, c = floor(m/2) (so the old centre cell [c, c+1] appears twice and
      the corner / wall bands keep their positions relative to the nearest wall).  The result is D4-symmetric,
      but NOT a valid cover (seams at the duplicated cell); it is meant as the column set of
      `closed4.py run --from-cert` (its weights only choose the initial rows).
  neartile OUT --s S [--rad 0.012] [--step 0.003] [--dmax 1.5] [--dstep 0.02]
      the near-tile pose pool of S21_COVER.md sec 1: every unit-tile centre (i+1/2, j+1/2) +- rad on a
      step grid, theta = dstep..dmax degrees, clipped into the open admissible box; `cx cy theta_rad` rows
      for `closed4.py --seed-rows [--lazy-seed]`.
"""
import sys, math, argparse
import numpy as np


def transplant(inp, out):
    tok = open(inp).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    assert sd == 1, "integer container side expected"
    m = sn; c = m // 2; K2 = (m + 1) * D
    acc = {}
    for i in range(n):
        X, Y, w = map(int, tok[5 + 3 * i: 8 + 3 * i])
        xs = ([X] if X < (c + 1) * D else []) + ([X + D] if X >= c * D else [])
        ys = ([Y] if Y < (c + 1) * D else []) + ([Y + D] if Y >= c * D else [])
        for x in xs:
            for y in ys: acc[(x, y)] = acc.get((x, y), 0) + w
    # D4 check
    for (x, y), w in acc.items():
        for (a, b) in ((K2 - x, y), (y, x)):
            assert acc.get((a, b)) == w, "transplant broke the D4 symmetry"
    with open(out, 'w') as f:
        f.write(f"{m + 1} 1\n{D}\n{W}\n{len(acc)}\n")
        for (x, y), w in sorted(acc.items()): f.write(f"{x} {y} {w}\n")
    print(f"{inp} (s={m}, {n} points, total {sum(map(int, tok[7::3][:n])) / W:.6f}) -> {out} "
          f"(s={m + 1}, {len(acc)} points, total {sum(acc.values()) / W:.6f})")


def neartile(out, s, rad=0.012, step=0.003, dmax=1.5, dstep=0.02):
    off = np.arange(-rad, rad + 1e-12, step)
    degs = np.arange(dstep, dmax + 1e-9, dstep)
    rows = []
    for d in degs:
        th = math.radians(d); w = (math.cos(th) + math.sin(th)) / 2; lo, hi = w + 1e-7, s - w - 1e-7
        for i in range(int(s)):
            for j in range(int(s)):
                for ox in off:
                    for oy in off:
                        rows.append((min(max(i + .5 + ox, lo), hi), min(max(j + .5 + oy, lo), hi), th))
    with open(out, 'w') as f:
        f.write(f"# cx cy theta_rad -- near-tile poses: tile centres +-{rad} (step {step}), theta {dstep}..{dmax} deg "
                f"(step {dstep}); grow_cover.py\n")
        for cx, cy, th in rows: f.write("%.17g %.17g %.17g\n" % (cx, cy, th))
    print(f"{len(rows)} near-tile poses written to {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['transplant', 'neartile'])
    ap.add_argument('args', nargs='+')
    ap.add_argument('--s', type=float, default=None)
    ap.add_argument('--rad', type=float, default=0.012); ap.add_argument('--step', type=float, default=0.003)
    ap.add_argument('--dmax', type=float, default=1.5); ap.add_argument('--dstep', type=float, default=0.02)
    a = ap.parse_args()
    if a.mode == 'transplant': transplant(a.args[0], a.args[1])
    else: neartile(a.args[0], a.s, a.rad, a.step, a.dmax, a.dstep)


if __name__ == '__main__':
    main()
