#!/usr/bin/env python3
"""m7_zuncert.py LOG [LOG ...] > POSES -- turn zmx2 UNCERT boxes into LP row poses `cx cy theta_rad`
(search/S45_COVER.md sec 10).  For every UNCERT box: its float-min pose, plus the box centre and the 4 (x,y) corners
at the box's middle angle (the boxes are tiny, so this is a local cut family).  HEURISTIC helper: the poses are
only used as extra LP rows; nothing is certified by this script."""
import sys, re, math

pat = re.compile(r"x\[([-\d.]+),([-\d.]+)\] y\[([-\d.]+),([-\d.]+)\] u\[([-\d.]+),([-\d.]+)\].*float-min ([-\d.]+) at "
                 r"\(([-\d.]+),([-\d.]+),u ([-\d.]+)\)")
seen = set()
for path in sys.argv[1:]:
    for line in open(path):
        if not line.startswith('UNCERT'): continue
        m = pat.search(line)
        if not m: continue
        x0, x1, y0, y1, u0, u1, fm, fx, fy, fu = map(float, m.groups())
        um = (u0 + u1) / 2
        out = [(fx, fy, 2 * math.atan(fu)), ((x0 + x1) / 2, (y0 + y1) / 2, 2 * math.atan(um))]
        out += [(x, y, 2 * math.atan(um)) for x in (x0, x1) for y in (y0, y1)]
        for p in out:
            k = tuple(round(v, 7) for v in p)
            if k in seen: continue
            seen.add(k); print(f"{p[0]:.9f} {p[1]:.9f} {p[2]:.12f}   # float-min {fm}")
