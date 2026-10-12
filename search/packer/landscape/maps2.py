#!/usr/bin/env python3
"""Hover-oriented map variants (10-11): inline family labels; two-level layout (families as non-overlapping discs,
members evenly spaced inside); crowding metric; within-family structure check.

  ./landscape/mapenv.sh python3 landscape/maps2.py        (run from packer/, after mapdata.py and maps.py)
"""
import collections, json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maps import RUN, REC, TOP, mds, sidebar

OUT = f'{RUN}/map'


def short(f):
    return f'{f[0]}/{f[1]}/{f[2]}/{f[3]}' + (f'/{f[4]}' if f[4] else '')


def crowding(E, frac=0.01):
    """Fraction of points whose nearest neighbour is closer than frac x the plot's longer side; median NN distance."""
    E = np.asarray(E, float); w = (E.max(0) - E.min(0)).max()
    d = np.sqrt(((E[:, None] - E[None]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
    nn = d.min(1) / w
    return float((nn < frac).mean()), float(np.median(nn))


def sunflower(k):
    i = np.arange(k) + .5
    r = np.sqrt(i / k); a = i * np.pi * (3 - np.sqrt(5))
    return np.c_[r * np.cos(a), r * np.sin(a)]


def twolevel(D, fam, iters=3000, fill=.45, seed=0):
    F = sorted(set(fam), key=lambda f: -fam.count(f))
    idx = {f: np.array([i for i, g in enumerate(fam) if g == f]) for f in F}
    med = {f: idx[f][np.argmin(D[np.ix_(idx[f], idx[f])].sum(1))] for f in F}
    Dm = np.array([[D[med[f], med[g]] for g in F] for f in F])
    P0 = mds(Dm); P0 /= np.abs(P0).max()
    cnt = np.array([len(idx[f]) for f in F], float)
    r = np.sqrt(cnt); r *= np.sqrt(fill * 4 / (np.pi * (r ** 2).sum()))      # discs cover ~fill of the [-1,1]^2 box
    r = np.maximum(r, .02)
    P = P0.copy()
    for it in range(iters):                     # push overlapping discs apart (+ small gap), weak spring to MDS
        dP = P[:, None] - P[None]; dist = np.sqrt((dP ** 2).sum(-1)) + 1e-12
        need = r[:, None] + r[None] + .015; ov = np.maximum(need - dist, 0); np.fill_diagonal(ov, 0)
        P += .25 * (ov[:, :, None] * dP / dist[:, :, None]).sum(1)
        P += .002 * (P0 - P)
    E = np.zeros((len(fam), 2))
    for k, f in enumerate(F):
        m = idx[f]; g = sunflower(len(m)) * r[k] * .92
        if len(m) >= 3:
            L = mds(D[np.ix_(m, m)]); L /= (np.abs(L).max() or 1)
            L = L * np.sqrt((g ** 2).sum(1).mean() / max((L ** 2).sum(1).mean(), 1e-12))
            row, col = linear_sum_assignment(((L[:, None] - g[None]) ** 2).sum(-1))
            E[m[row]] = P[k] + g[col]
        else:
            E[m] = P[k] + g
    return E, F, P, r


def within(D, fam, S):
    """Is there structure inside families?  Per family >= 15: share of within-family distance variance on the first 2 MDS
    axes, Spearman of axis 1 with side; within vs between nearest-neighbour distances."""
    rows = []
    for f, n in collections.Counter(fam).most_common():
        if n < 15: break
        m = np.array([i for i, g in enumerate(fam) if g == f]); Dw = D[np.ix_(m, m)]
        J = np.eye(n) - 1 / n; B = -.5 * J @ (Dw ** 2) @ J; w = np.sort(np.linalg.eigvalsh(B))[::-1]
        L = mds(Dw); rho = spearmanr(L[:, 0], S[m])[0]
        rows.append((short(f), n, np.median(Dw[np.triu_indices(n, 1)]), w[:2].sum() / w[w > 0].sum(), rho))
    return rows


def labelled(ax, E, fam, S, title, kmin=5, by='S'):
    T = minimum_spanning_tree(D).tocoo()
    for i, j in zip(T.row, T.col):
        ax.plot(E[[i, j], 0], E[[i, j], 1], color='0.85', lw=.5, zorder=0)
    sc = ax.scatter(E[:, 0], E[:, 1], c=S, cmap='viridis_r', vmin=REC, vmax=TOP, s=10, zorder=2)
    b = int(np.argmin(S)); ax.scatter(*E[b], marker='*', s=200, c='red', edgecolor='k', zorder=5)
    cnt = collections.Counter(fam)
    for f, n in cnt.items():
        if n < kmin: continue
        m = [i for i, g in enumerate(fam) if g == f]; c = np.median(E[m], 0)
        top = E[m][:, 1].max()
        ax.text(c[0], top, short(f), fontsize=7, ha='center', va='bottom', zorder=6,
                bbox=dict(boxstyle='round,pad=.15', fc='white', ec='none', alpha=.75))
    ax.set_title(title, fontsize=10); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal', adjustable='datalim')
    return sc


M = json.load(open(f'{OUT}/minima.json')); D = np.load(f'{OUT}/fpD.npy')
fam = [tuple(m['fam']) for m in M]; S = np.array([float(m['S']) for m in M])
Z = dict(np.load(f'{OUT}/layouts.npz'))
E2, F, P, r = twolevel(D, fam)
np.save(f'{OUT}/twolevel.npy', E2)
print('crowding: share with nearest neighbour < 1% of width, median NN / width')
for k, E in list(Z.items()) + [('two-level', E2)]:
    c, med = crowding(E); print(f'  {k:28s} {c:6.1%}  {med:.4f}')
print('within-family structure (families >= 15): median dist, var on 2 MDS axes, Spearman(axis1, side)')
for row in within(D, fam, S):
    print('  %-14s n=%3d  med %.4f  var2 %.2f  rho %+.2f' % row)
nn = np.sort(D, 1)[:, 1]
same = np.array([[fam[i] == fam[j] for j in range(len(fam))] for i in range(len(fam))])
print('median within-family distance %.4f, between %.4f, NN %.4f' % (np.median(D[same & ~np.eye(len(D), dtype=bool)]),
                                                                     np.median(D[~same]), np.median(nn)))
fig, ax = plt.subplots(1, 2, figsize=(22, 11))
labelled(ax[0], Z['t-SNE_perplexity_150'], fam, S, 't-SNE perplexity 150, families >= 5 labelled (L/B/axis/J)')
sc = labelled(ax[1], E2, fam, S, 'Two-level: discs per family (area ~ count, placed by family distance, overlap removed); '
              'members evenly spaced inside', kmin=3)
for k, f in enumerate(F):
    ax[1].add_patch(plt.Circle(P[k], r[k], fill=False, ec='0.75', lw=.6))
sidebar(sc, ax[1])
fig.suptitle('s(110): 740 certified minima, colour = side (record to 11), grey = MST of fingerprint distance')
plt.tight_layout(); plt.savefig(f'{OUT}/labelled.png', dpi=72)
print('->', f'{OUT}/labelled.png')
