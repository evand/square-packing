#!/usr/bin/env python3
"""family_rows.py -- emit the pose families that `search/closed4.py`'s row lattice steps over.

`search/RUNG2.md` §4.3: the 163 boxes the exact checker could not certify on
`runs/closed4_best_x103.txt` turned out to be a genuine violation of the cover, at

    (c_x, c_y, theta) = (3/2, 1461/2000, 2 arctan(1/40000))   capture 9420217/10000000 = 0.9420217

and its 8 dihedral images.  This is `ZEROMARGIN.md` §4 item 3: an axis-parallel square
`[1,2] x [c_y-1/2, c_y+1/2]` with BOTH vertical edges on the heavy grid lines `x = 1` and `x = 2`,
sliding freely in `c_y`.  `closed4.py` places its rows on a `--pitch` lattice (0.02 by default) in
the rotated frame at each angle of `angle_list()`, so it samples `c_x = 3/2` but steps over the
dip in `c_y`, which is a few thousandths wide; the LP therefore never sees the constraint.

This script writes those poses explicitly, in the `cx cy theta_rad` format that
`closed4.py --seed-rows` reads (and that `zeromargin.py --oracle` writes):

  * **item 3**: `c_x` on each half-integer line (both square edges then sit on grid lines), `c_y`
    on a fine pitch across the admissible range, at a geometric ladder of small angles and their
    reflections; and the transpose (`c_y` fixed, `c_x` free).
  * **tilted**: a local cluster around a named pose and its dihedral images -- used for the
    `(3.3965, 1.4235, 76.43deg)` family of `FAMILY.md` §2, whose `D4` images are where the
    right-hand half of the residue sits.

usage: python3 search/family_rows.py OUT.txt [--s 4] [--pitch 0.0025] [--tilted]
"""
import argparse, math
import numpy as np


def d4(x, y, s):
    return [(x, y), (s - x, y), (x, s - y), (s - x, s - y),
            (y, x), (s - y, x), (y, s - x), (s - y, s - x)]


def item3(s, pitch, angles):
    """squares with both edges of one axis on grid lines, sliding along the other axis"""
    out = []
    lines = [k + 0.5 for k in range(int(s))]          # 0.5, 1.5, ..., s-0.5
    for th in angles:
        w = (math.cos(th) + math.sin(th)) / 2
        lo, hi = w + 1e-7, s - w - 1e-7
        if hi <= lo: continue
        free = np.arange(lo, hi + 1e-12, pitch)
        for c in lines:
            cc = min(max(c, lo), hi)
            for v in free:
                out.append((cc, v, th))              # c_x fixed, c_y free
                out.append((v, cc, th))              # and the transpose
    return out


def tilted(s, poses, nring, rad, arad):
    """a local cluster around each named pose and each of its dihedral images"""
    out = []
    rng = np.random.default_rng(7)
    for (cx, cy, th) in poses:
        for (ax, ay) in d4(cx, cy, s):
            for k in range(nring):
                r = rad * (k + 1) / nring
                for a in np.linspace(0, 2 * math.pi, 12, endpoint=False):
                    for dt in (-arad, 0.0, arad):
                        t = th + dt
                        if not (0 <= t <= math.pi / 2): continue
                        w = (math.cos(t) + math.sin(t)) / 2
                        lo, hi = w + 1e-7, s - w - 1e-7
                        if hi <= lo: continue
                        out.append((min(max(ax + r * math.cos(a), lo), hi),
                                    min(max(ay + r * math.sin(a), lo), hi), t))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--s', type=float, default=4.0)
    ap.add_argument('--pitch', type=float, default=0.0025)
    ap.add_argument('--tilted', action='store_true')
    a = ap.parse_args()
    angles = [0.0] + [10 ** e for e in (-6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5)]
    rows = item3(a.s, a.pitch, angles)
    if a.tilted:
        rows += tilted(a.s, [(3.3965, 1.4235, math.radians(76.4282))], 3, 0.02, math.radians(0.4))
        rows += tilted(a.s, [(1.5, 0.7305, 2.86e-5)], 3, 0.01, 1e-5)
    with open(a.out, 'w') as f:
        f.write("# cx cy theta_rad -- search/family_rows.py\n")
        for cx, cy, th in rows:
            f.write("%.17g %.17g %.17g\n" % (float(cx), float(cy), float(th)))
    print(f"{len(rows)} poses written to {a.out} "
          f"(item-3 pitch {a.pitch}, {len(angles)} angles, tilted={a.tilted})")


if __name__ == '__main__':
    main()
