#!/usr/bin/env python3
"""Render parsed squares to PNG (matplotlib) for visual checks against the source SVG."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_svg import parse
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
src, out = sys.argv[1], sys.argv[2]
s, sq = parse(src); s = float(s)
fig, ax = plt.subplots(figsize=(7, 7)); ax.set_aspect('equal'); ax.set_xlim(-0.05*s, 1.05*s); ax.set_ylim(-0.05*s, 1.05*s)
ax.add_patch(Rectangle((0, 0), s, s, fill=False, lw=1.5))
for i, (cx, cy, th) in enumerate(sq):
    cx, cy, th = float(cx), float(cy), float(th); c, sn = math.cos(math.radians(th)), math.sin(math.radians(th))
    pts = [(cx + c*dx - sn*dy, cy + sn*dx + c*dy) for dx, dy in ((-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5))]
    ax.add_patch(Polygon(pts, closed=True, fc='#B2B2B2' if th < 1e-9 else '#7fb3d5', ec='k', lw=0.4, alpha=0.6))
    if len(sq) <= 60: ax.text(cx, cy, str(i), ha='center', va='center', fontsize=6)
ax.set_title('%s: %d squares, s=%.6f' % (os.path.basename(src), len(sq), s)); ax.axis('off')
fig.savefig(out, dpi=110, bbox_inches='tight'); print(out, len(sq))
