#!/usr/bin/env python3
"""Level-2 region bookkeeping: split a fractional packing measure (a branch dual or a
`packing_dual.py` / `clique_ceiling.py` support file) into corner boxes, wall slots and the
interior, by the CENTRE of the pose.

Regions for a container [0,t]^2 and a corner size r (default 1):
  C_j   the four corner boxes [0,r]^2 and images                                      (4)
  W_j   the eight wall slots [r, t/2] x [0, r], [t/2, t-r] x [0, r] and images         (8)
  I     the interior [r, t-r]^2                                                        (1)
These partition the container (up to shared boundaries), so the masses add up to the total.

Usage:  python3 search/level2_regions.py FILE [--t T] [--r R]
FILE is either a branch dual (`cx cy theta h flag y`, mass y on the ORBIT) or a support file
(`pose cx cy theta_deg mu`) or an exact support (`pose p q cx cy mass`).
"""
import argparse, math, sys
from fractions import Fraction

import numpy as np


def read_measure(path):
    """-> (t, [(cx, cy, theta_rad, mass_of_orbit)]).  mass is the mass of the whole D4 orbit."""
    t = None
    out = []
    kind = None
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            for tok in s.replace(',', ' ').split():
                if tok.startswith('s=') and '/' in tok:
                    a, b = tok[2:].split('/')
                    t = float(a) / float(b)
                elif tok.startswith('t='):
                    try:
                        t = float(tok[2:])
                    except ValueError:
                        pass
                elif tok == 't' :
                    kind = kind
            if 't =' in s and t is None:
                q = s.split('t =')[1].split()[0]
                if '/' in q:
                    a, b = q.split('/')
                    t = float(a) / float(b)
                else:
                    t = float(q.rstrip(';'))
            continue
        q = s.split()
        if q[0] == 'pose':
            if len(q) == 5:                      # cx cy theta_deg mu
                out.append((float(q[1]), float(q[2]), math.radians(float(q[3])), float(q[4])))
            elif len(q) == 6:                    # p q cx cy mass  (theta = 2 arctan(p/q))
                p, qq = Fraction(q[1]), Fraction(q[2])
                th = 2 * math.atan(float(p) / float(qq))
                out.append((float(Fraction(q[3])), float(Fraction(q[4])), th, float(Fraction(q[5]))))
        elif len(q) == 6:                        # branch dual: cx cy theta h flag y
            out.append((float(q[0]), float(q[1]), float(q[2]), float(q[5])))
    if t is None:
        raise SystemExit(f'{path}: could not read t from the header')
    return t, out


def images(cx, cy, th, t):
    return [(cx, cy, th), (t - cx, cy, -th), (cx, t - cy, -th), (t - cx, t - cy, th),
            (cy, cx, -th), (t - cy, cx, th), (cy, t - cx, th), (t - cy, t - cx, -th)]


def classify(cx, cy, t, r):
    """region label of a centre: ('C', j) / ('W', j) / ('I', 0); slots numbered 0..7 anticlockwise
    from the one on the bottom wall next to the origin corner."""
    lo, hi, h = r, t - r, t / 2
    inx = lo <= cx <= hi
    iny = lo <= cy <= hi
    if inx and iny:
        return ('I', 0)
    if not inx and not iny:                                     # corner box
        j = (1 if cx > h else 0) + (2 if cy > h else 0)
        return ('C', j)
    if iny:                                                     # left or right wall
        if cx < h:                                              # left wall
            return ('W', 7 if cy < h else 6)
        return ('W', 2 if cy < h else 3)
    if cx < h:                                                  # bottom-left / top-left
        return ('W', 0 if cy < h else 5)
    return ('W', 1 if cy < h else 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('FILE')
    ap.add_argument('--t', type=float, default=None)
    ap.add_argument('--r', type=float, default=1.0)
    ap.add_argument('--top', type=int, default=12, help='heaviest poses per region class to print')
    a = ap.parse_args()

    t, meas = read_measure(a.FILE)
    if a.t:
        t = a.t
    r = a.r
    total = sum(m for *_, m in meas)
    mc = np.zeros(4); mw = np.zeros(8); mi = 0.0
    detail = {'C': [], 'W': [], 'I': []}
    for (cx, cy, th, m) in meas:
        for (x, y, u) in images(cx, cy, th, t):
            k, j = classify(x, y, t, r)
            if k == 'C':
                mc[j] += m / 8
            elif k == 'W':
                mw[j] += m / 8
            else:
                mi += m / 8
        k, j = classify(cx, cy, t, r)
        detail[k].append((m, cx, cy, math.degrees(th)))

    print(f'{a.FILE}\n  t = {t:.6f}  r = {r}  total mass = {total:.6f}')
    print(f'  corners  {mc.sum():8.4f}   per box {np.round(mc, 4).tolist()}')
    print(f'  wall     {mw.sum():8.4f}   per slot {np.round(mw, 4).tolist()}')
    print(f'  interior {mi:8.4f}')
    print(f'  check    {mc.sum() + mw.sum() + mi:8.4f}')
    for k, name in (('C', 'corner'), ('W', 'slot'), ('I', 'interior')):
        d = sorted(detail[k], reverse=True)[:a.top]
        print(f'  --- heaviest orbits with representative in a {name} region ---')
        for (m, cx, cy, th) in d:
            print(f'      mu={m:9.5f}  ({cx:.4f}, {cy:.4f}, {th:6.2f} deg)')


if __name__ == '__main__':
    main()
