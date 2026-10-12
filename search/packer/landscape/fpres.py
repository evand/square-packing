#!/usr/bin/env python3
"""Fingerprint resolution / blur vs exact matching distance (10-11), on runs/count_sep/k150.json (150 minima, all pairs).
Scores: Spearman on all pairs; Spearman on near pairs (exact <= its 25th percentile); same nearest neighbour; mean overlap
of each minimum's 10 nearest neighbours (what UMAP uses).  Blur in box units = sig * s / g.

  python3 landscape/fpres.py        (run from packer/; ~1 min)
"""
import json, sys, os, time
sys.path.insert(0, '.')
import numpy as np
from scipy.stats import spearmanr
from pk.store import Store
from pk import refset, fingerprint as fp

d = json.load(open('runs/count_sep/k150.json')); st = Store()
R = {m['key']: m for m in refset.load(d['ref'])['minima']}
P = [st.get(R[k]['rep']) for k in d['keys']]
P = [(float(p.s), [tuple(map(float, r)) for r in p.sq]) for p in P]
K = len(P); E = np.full((K, K), np.inf)
for i, j, x in d['pairs']:
    E[i, j] = E[j, i] = x
iu = np.triu_indices(K, 1); ex = E[iu]; near = ex <= np.percentile(ex, 25)
kn = lambda M: np.argsort(M, 1)[:, :10]
EK = kn(E)
print(f'{"grid":>5} {"sig":>5} {"blur":>6} | {"rho all":>7} {"rho near":>8} {"sameNN":>6} {"10NN ovl":>8} | s/pair')
for g, sig in ((12, .3), (24, .6), (24, 1.2), (48, 1.2), (48, .6), (96, 1.2), (96, .6), (192, 1.2), (192, .6)):
    t = time.time()
    I = np.stack([fp.image(s, sq, g, sig) for s, sq in P]).reshape(K, -1).astype(np.float32)
    SY = np.stack([fp.syms(fp.image(s, sq, g, sig)).reshape(8, -1) for s, sq in P]).astype(np.float32)
    F = np.empty((K, K))
    for i in range(K):
        F[i] = np.abs(SY[None, i] - I[:, None]).sum(-1).min(-1) / 110
    F = np.minimum(F, F.T); np.fill_diagonal(F, np.inf)
    f = F[iu]; dt = (time.time() - t) / (K * (K - 1) / 2)
    ovl = np.mean([len(set(a) & set(b)) / 10 for a, b in zip(EK, kn(F))])
    print(f'{g:5d} {sig:5.2f} {sig * 11 / g:6.3f} | {spearmanr(ex, f)[0]:7.3f} {spearmanr(ex[near], f[near])[0]:8.3f} '
          f'{np.mean(E.argmin(1) == F.argmin(1)):6.0%} {ovl:8.0%} | {dt * 1e6:.0f} us', flush=True)
