#!/usr/bin/env python3
"""Why do count-families split into two islands? (10-11)  For chosen families: fingerprint-distance clusters (threshold
0.02), the medoid of each cluster, drawn side by side after mapping it to the box symmetry that best matches the record.

  ./landscape/mapenv.sh python3 landscape/splits.py 21/19/67/3 19/18/69/4 20/19/67/4     (run from packer/)"""
import collections, json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pk import fingerprint as fp
O = 'runs/landscape110/map'
M = json.load(open(f'{O}/minima.json')); D = np.load(f'{O}/fpD.npy')
fam = ['/'.join(map(str, m['fam'][:4])) + (f"/{m['fam'][4]}" if m['fam'][4] else '') for m in M]
S = [float(m['S']) for m in M]


def dec(m): q = m['sq']; return np.array([(q[i] / 1e4, q[i + 1] / 1e4, q[i + 2] / 100) for i in range(0, len(q), 3)])


def sym(sq, s, t, fy, fx):
    x, y, th = sq[:, 0].copy(), sq[:, 1].copy(), sq[:, 2].copy()
    if t: x, y, th = y, x, 90 - th
    if fy: y, th = s - y, -th
    if fx: x, th = s - x, -th
    return np.c_[x, y, th % 90]


rec = int(np.argmin(S)); R = fp.image(S[rec], dec(M[rec]))


def orient(i):
    sq = dec(M[i]); s = S[i]
    best = min(((t, fy, fx) for t in (0, 1) for fy in (0, 1) for fx in (0, 1)),
               key=lambda c: np.abs(fp.image(s, sym(sq, s, *c)) - R).sum())
    return sym(sq, s, *best)


def draw(ax, sq, s, title):
    ax.add_patch(plt.Rectangle((0, 0), s, s, fill=False, lw=1.2))
    for x, y, t in sq:
        r = np.radians(t); c, sn = np.cos(r) / 2, np.sin(r) / 2
        a = t % 90; off = min(a, 90 - a)
        col = '#E6E9E5' if off < 1.5 else ('#0B7A9B' if a < 40 else '#B98A3C' if a > 50 else '#C0392B')
        ax.add_patch(Polygon([(x + c - sn, y + sn + c), (x - c - sn, y - sn + c), (x - c + sn, y - sn - c), (x + c + sn, y + sn - c)],
                             fc=col, ec='#3A4448', lw=.5))
    ax.set_xlim(-.2, 11.2); ax.set_ylim(-.2, 11.2); ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, fontsize=9)


fams = sys.argv[1:]
fig, axs = plt.subplots(len(fams), 3, figsize=(15, 5 * len(fams)))
for row, f in zip(np.atleast_2d(axs), fams):
    m = np.array([i for i, g in enumerate(fam) if g == f])
    k = connected_components(D[np.ix_(m, m)] < .02, directed=False)[1]
    cl = [c for c, _ in collections.Counter(k).most_common(2)]
    for ax, c in zip(row, cl):
        mm = m[k == c]; med = mm[np.argmin(D[np.ix_(mm, mm)].sum(1))]
        draw(ax, orient(med), S[med], f'{f}, island of {len(mm)}: medoid side {S[med]:.7f}\n(best in island {min(S[i] for i in mm):.7f})')
    draw(row[2], orient(rec), S[rec], f'record (22/21/67/0) for reference, side {S[rec]:.7f}')
fig.suptitle('Count-family splits: medoid of each fingerprint cluster, mapped to the box symmetry closest to the record.  '
             'Blue ~28 deg band, gold ~62 deg band, red 40-50 deg junction squares', fontsize=10)
plt.tight_layout(); plt.savefig(f'{O}/splits.png', dpi=70); print(f'{O}/splits.png')
