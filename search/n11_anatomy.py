#!/usr/bin/env python3
"""Anatomy of a D4-symmetrised fractional packing measure of closed unit squares in [0,t]^2.

Reads either
  * a `packing_dual.py` support file   `pose cx cy theta_deg mu`   (floats), or
  * a `dual_exact.py` support file     `pose p q cx cy mass`       (exact; theta = 2 atan(p/q)),
both with the `# t=...` header and the convention that a pose carries `mu/8` on each of its 8
dihedral images (`leaf_ceiling.d4_images` / `packing_dual.pose_images`).

Reports, for the measure as a whole (all 8*len(poses) images):
  mass, concentration, angle histogram, mass by region of the *centre*, and the region
  decomposition of the coverage integral  int_R cov = sum_i m_i * area(Q_i ^ R)  over
  the four corner unit boxes, the four wall strips and the central box, plus a coverage map.

Everything here is descriptive: it decides nothing and certifies nothing (the certified numbers
are `dual_exact.py check`'s `M` and `L`).  Areas are computed by exact convex polygon clipping
in floating point (Sutherland-Hodgman against the four half-planes of the box).

    python3 search/n11_anatomy.py runs/dual_exact_3.82_support.txt [--t 191/50] [--map 40]
"""
import argparse
import math
import os
import sys
from fractions import Fraction as Fr


def read_any(path, t_override=None):
    """-> (t, [(cx, cy, theta_rad, mu)])  one entry per POSE (mass mu, spread over 8 images)."""
    t = t_override
    out = []
    for line in open(path):
        s = line.strip()
        if s.startswith('#'):
            if t is None and 't=' in s:
                t = float(Fr(s.split('t=')[1].split()[0]))
            elif t is None and 't =' in s:
                t = float(Fr(s.split('t =')[1].split()[0]))
            continue
        q = s.split('#')[0].split()
        if not q:
            continue
        if q[0] == 't' and t is None:
            t = float(Fr(q[1]))
        if q[0] != 'pose':
            continue
        if len(q) == 5:                                   # packing_dual: cx cy theta_deg mu
            out.append((float(q[1]), float(q[2]), math.radians(float(q[3])), float(q[4])))
        elif len(q) == 6:                                 # dual_exact: p q cx cy mass
            p, qq = int(q[1]), int(q[2])
            out.append((float(Fr(q[3])), float(Fr(q[4])), 2.0 * math.atan2(p, qq), float(Fr(q[5]))))
        else:
            raise SystemExit(f'{path}: unrecognised pose line with {len(q)-1} fields')
    if t is None:
        raise SystemExit(f'{path}: no container side; pass --t')
    return t, out


def images(cx, cy, th, t):
    """the 8 dihedral images, the convention of packing_dual.pose_images / leaf_ceiling.d4_images"""
    return [(cx, cy, th), (t - cx, cy, -th), (cx, t - cy, -th), (t - cx, t - cy, th),
            (cy, cx, -th), (t - cy, cx, th), (cy, t - cx, th), (t - cy, t - cx, -th)]


def corners_of(cx, cy, th):
    c, s = math.cos(th), math.sin(th)
    return [(cx + c * a - s * b, cy + s * a + c * b) for a, b in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]


def clip_box(poly, x0, y0, x1, y1):
    """Sutherland-Hodgman clip of a convex polygon against the axis box."""
    for k, (sign, lim) in enumerate(((1, x0), (-1, -x1), (1, y0), (-1, -y1))):
        ax = k < 2
        out = []
        n = len(poly)
        if n == 0:
            return out
        for i in range(n):
            p, q = poly[i], poly[(i + 1) % n]
            dp = sign * (p[0] if ax else p[1]) - lim
            dq = sign * (q[0] if ax else q[1]) - lim
            if dp >= 0:
                out.append(p)
            if (dp > 0) != (dq > 0) and dp != dq:
                r = dp / (dp - dq)
                out.append((p[0] + r * (q[0] - p[0]), p[1] + r * (q[1] - p[1])))
        poly = out
    return poly


def area(poly):
    a = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        a += x0 * y1 - x1 * y0
    return abs(a) / 2.0


def inside(px, py, cx, cy, th, tol=1e-9):
    c, s = math.cos(th), math.sin(th)
    dx, dy = px - cx, py - cy
    return abs(dx * c + dy * s) <= .5 + tol and abs(-dx * s + dy * c) <= .5 + tol


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--t', default=None)
    ap.add_argument('--map', type=int, default=0, help='ASCII coverage map at this resolution')
    ap.add_argument('--top', type=int, default=16)
    a = ap.parse_args()
    tov = float(Fr(a.t)) if a.t else None
    for path in a.files:
        t, poses = read_any(path, tov)
        mass = sum(p[3] for p in poses)
        print(f"\n=== {path}   t = {t:.6f}   poses = {len(poses)}   images = {8*len(poses)}   mass = {mass:.7f}")
        # ---------------------------------------------------------------- concentration
        mu = sorted((p[3] for p in poses), reverse=True)
        cum = 0.0
        marks = [.25, .5, .75, .9, .99]
        j = 0
        for i, v in enumerate(mu):
            cum += v
            while j < len(marks) and cum >= marks[j] * mass:
                print(f"  {marks[j]*100:5.1f}% of the mass on the {i+1} heaviest poses")
                j += 1
        # ---------------------------------------------------------------- angles
        print("  angle histogram (pose mass by 5 deg bin of |theta| in [0,45]):")
        for lo in range(0, 45, 5):
            sel = [p for p in poses if lo <= abs(math.degrees(p[2])) < lo + 5 or (lo == 40 and abs(math.degrees(p[2])) >= 45 - 1e-9)]
            m = sum(p[3] for p in sel)
            if m > 1e-9:
                print(f"    [{lo:2d},{lo+5:2d}) deg : mass {m:8.4f} ({100*m/mass:5.1f}%)  poses {len(sel)}")
        for thr in (0.01, 0.1, 0.5, 1.0, 2.0, 5.0):
            m = sum(p[3] for p in poses if abs(math.degrees(p[2])) < thr)
            print(f"    |theta| < {thr:4.2f} deg : mass {m:8.4f} ({100*m/mass:5.1f}%)")
        m45 = sum(p[3] for p in poses if abs(abs(math.degrees(p[2])) - 45) < 1.0)
        print(f"    |theta-45| < 1 deg   : mass {m45:8.4f} ({100*m45/mass:5.1f}%)")
        # ---------------------------------------------------------------- centres
        print("  pose centre, distance to the nearest wall (canonical representative):")
        for lo, hi in ((0, .52), (.52, .8), (.8, 1.2), (1.2, 1.6), (1.6, 2.4)):
            sel = [p for p in poses if lo <= min(p[0], t - p[0], p[1], t - p[1]) < hi]
            m = sum(p[3] for p in sel)
            if m > 1e-9:
                print(f"    [{lo:4.2f},{hi:4.2f}) : mass {m:8.4f} ({100*m/mass:5.1f}%)  poses {len(sel)}")
        # ---------------------------------------------------------------- regions
        # int_R cov = sum over IMAGES i of m_i * area(Q_i ^ R);  m_i = mu/8.
        w = t - 2.0                                        # wall-strip length
        regions = [("corner box  [0,1]^2", [(0, 0, 1, 1), (t - 1, 0, t, 1), (0, t - 1, 1, t), (t - 1, t - 1, t, t)]),
                   (f"wall strip  [1,{t-1:.3f}]x[0,1]", [(1, 0, t - 1, 1), (1, t - 1, t - 1, t),
                                                         (0, 1, 1, t - 1), (t - 1, 1, t, t - 1)]),
                   (f"centre box  [1,{t-1:.3f}]^2", [(1, 1, t - 1, t - 1)])]
        print(f"  coverage integral by region (sum over regions = mass; area(C) = {t*t:.4f}):")
        tot_check = 0.0
        for name, boxes in regions:
            acc = 0.0
            ar = sum((b[2] - b[0]) * (b[3] - b[1]) for b in boxes)
            for cx, cy, th, m in poses:
                for icx, icy, ith in images(cx, cy, th, t):
                    poly = corners_of(icx, icy, ith)
                    for b in boxes:
                        acc += (m / 8.0) * area(clip_box(poly, *b))
            tot_check += acc
            print(f"    {name:28s} area {ar:7.4f} (x{len(boxes)})  mass {acc:8.4f} ({100*acc/mass:5.1f}%)  mean cov {acc/ar:.4f}")
        print(f"    {'sum':28s} {'':16s}  mass {tot_check:8.4f}  (should equal {mass:.4f})")
        # ---------------------------------------------------------------- heaviest poses
        print(f"  heaviest {a.top} poses (cx cy theta_deg  mu   mu/8 per image):")
        for cx, cy, th, m in sorted(poses, key=lambda p: -p[3])[:a.top]:
            d = min(cx, t - cx, cy, t - cy)
            print(f"    {cx:9.5f} {cy:9.5f} {math.degrees(th):10.5f}  {m:8.5f}  {m/8:7.5f}   wall-dist {d:.4f}")
        # ---------------------------------------------------------------- map
        if a.map:
            n = a.map
            print(f"  coverage map ({n}x{n} cell centres; digit = floor(10*cov), # = cov in [0.95,1], * = >1):")
            imgs = [(icx, icy, ith, m / 8.0) for cx, cy, th, m in poses for (icx, icy, ith) in images(cx, cy, th, t)]
            for jj in range(n - 1, -1, -1):
                row = ''
                py = (jj + .5) * t / n
                for ii in range(n):
                    px = (ii + .5) * t / n
                    c = sum(mm for (icx, icy, ith, mm) in imgs if inside(px, py, icx, icy, ith))
                    row += '*' if c > 1.0000001 else ('#' if c >= 0.95 else str(min(9, int(c * 10))))
                print('   ' + row)


if __name__ == '__main__':
    main()
