#!/usr/bin/env python3
"""High-res renders of chosen sweep layouts (10-11): colour = side (record to 11), family labels per island.

  ./landscape/mapenv.sh python3 landscape/render.py [--align] 'umap_p=0.5_k=250_md=2_spread=4' ...   (run from packer/)
--align rotates the principal axis horizontal (a layout's orientation is arbitrary).
"""
import collections, json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maps import RUN, REC, TOP, sidebar
OUT = f'{RUN}/map'
M = json.load(open(f'{OUT}/minima.json')); fam = [tuple(m['fam']) for m in M]; S = np.array([float(m['S']) for m in M])
Z = np.load(f'{OUT}/sweep.npz'); res = {r[0].replace(' ', '_'): r[1:] for r in json.load(open(f'{OUT}/sweep.json'))}


def short(f): return f'{f[0]}/{f[1]}/{f[2]}/{f[3]}' + (f'/{f[4]}' if f[4] else '')


ALIGN = '--align' in sys.argv


def align(E):
    """Rotate so the principal axis is horizontal (distances unchanged)."""
    E = E - E.mean(0); w, v = np.linalg.eigh(np.cov(E.T)); return E @ v[:, ::-1]


for name in [a for a in sys.argv[1:] if not a.startswith('--')]:
    E = align(Z[name]) if ALIGN else Z[name]; w = (E.max(0) - E.min(0)).max()
    fig, ax = plt.subplots(figsize=(15, 13))
    o = np.argsort(-S)                       # draw low sides last (on top)
    sc = ax.scatter(E[o, 0], E[o, 1], c=S[o], cmap='viridis_r', vmin=REC, vmax=TOP, s=14, lw=.2, edgecolor='0.3')
    b = int(np.argmin(S)); ax.scatter(*E[b], marker='*', s=260, c='red', edgecolor='k', zorder=5)
    for f, n in collections.Counter(fam).items():
        m = np.array([i for i, g in enumerate(fam) if g == f])
        if n < 4: continue
        P = E[m]; A = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)) < .04 * w
        k, lab = connected_components(A, directed=False)
        for c in range(k):
            q = P[lab == c]
            if len(q) < 4: continue
            ax.annotate(short(f), (np.median(q[:, 0]), q[:, 1].max()), xytext=(0, 4), textcoords='offset points',
                        fontsize=7.5, ha='center', color='0.15',
                        bbox=dict(boxstyle='round,pad=.12', fc='white', ec='none', alpha=.7))
    c, tr, gl = res[name]
    ax.set_title(f'{name}: crowd {c:.1%}, trust {tr:.3f}, global {gl:+.2f}   ({len(S)} certified minima; star = record; '
                 'labels L/B/axis/J per island of >= 4)', fontsize=10)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal', adjustable='datalim'); sidebar(sc, ax)
    plt.tight_layout(); fn = f'{OUT}/r_{name.replace("=", "").replace(".", "")}{"_al" if ALIGN else ""}.png'; plt.savefig(fn, dpi=100); plt.close()
    print(fn)
