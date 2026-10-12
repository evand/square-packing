#!/usr/bin/env python3
"""Angle-aware fingerprint variants vs exact matching distance (10-11), same harness / scores as fpres.py.

  cur        current: centres, channels axis / tilted (> 3 deg)
  ang a      centres; channels density, a*cos(4 theta), a*sin(4 theta) (weighted by each square's blob); the sin
             channel changes sign under the 4 reflections among the 8 box symmetries
  cur+ang a  current 2 channels + a*cos4, a*sin4
  corners    the 4 corners of each square (weight 1/4), one channel
  corners+c  corners + centres channels

  python3 landscape/fpangle.py        (run from packer/; ~1 min)
"""
import json, sys, time
sys.path.insert(0, '.')
import numpy as np
from scipy.stats import spearmanr
from pk.store import Store
from pk import refset

G, SIG = 24, .6
d = json.load(open('runs/count_sep/k150.json')); st = Store()
R = {m['key']: m for m in refset.load(d['ref'])['minima']}
P = [st.get(R[k]['rep']) for k in d['keys']]
P = [(float(p.s), np.array([tuple(map(float, r)) for r in p.sq])) for p in P]
K = len(P); E = np.full((K, K), np.inf)
for i, j, x in d['pairs']:
    E[i, j] = E[j, i] = x
iu = np.triu_indices(K, 1); ex = E[iu]; near = ex <= np.percentile(ex, 25)
EK = np.argsort(E, 1)[:, :10]


def blobs(s, xy, w):
    """Channels: sum over points of w[:, c] * Gaussian blob at xy (box coords / s).  Returns C x G x G."""
    u = np.clip(xy[:, 0] / s, 0, 1) * G - .5; v = np.clip(xy[:, 1] / s, 0, 1) * G - .5; c = np.arange(G)
    wx = np.exp(-.5 * ((c[None] - u[:, None]) / SIG) ** 2); wy = np.exp(-.5 * ((c[None] - v[:, None]) / SIG) ** 2)
    wx /= wx.sum(1, keepdims=True); wy /= wy.sum(1, keepdims=True)
    return np.stack([(wy * w[:, k:k + 1]).T @ wx for k in range(w.shape[1])])


def corners(sq):
    t = np.radians(sq[:, 2]); c, s_ = np.cos(t) * .5, np.sin(t) * .5
    out = []
    for a, b in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        out.append(np.c_[sq[:, 0] + a * c - b * s_, sq[:, 1] + a * s_ + b * c])
    return np.concatenate(out)


def variant(name, a=1.0):
    def f(s, sq):
        th = np.radians(4 * (sq[:, 2] % 90)); off = np.minimum(sq[:, 2] % 90, 90 - sq[:, 2] % 90)
        tilt = (off > 3).astype(float); one = np.ones(len(sq))
        if name == 'cur':
            return blobs(s, sq[:, :2], np.c_[1 - tilt, tilt]), []
        if name == 'ang':
            return blobs(s, sq[:, :2], np.c_[one, a * np.cos(th), a * np.sin(th)]), [2]
        if name == 'cur+ang':
            return blobs(s, sq[:, :2], np.c_[1 - tilt, tilt, a * np.cos(th), a * np.sin(th)]), [3]
        if name == 'corners':
            return blobs(s, corners(sq), np.full((4 * len(sq), 1), .25)), []
        if name == 'corners+c':
            return np.concatenate([blobs(s, corners(sq), np.full((4 * len(sq), 1), .25)), blobs(s, sq[:, :2], one[:, None])]), []
    return f


def syms(im, odd):
    out = []
    for t in (0, 1):
        z0 = np.transpose(im, (0, 2, 1)) if t else im
        for fy in (0, 1):
            for fx in (0, 1):
                z = z0[:, ::-1, :] if fy else z0
                z = z[:, :, ::-1] if fx else z
                if (t + fy + fx) % 2 and odd:
                    z = z.copy(); z[odd] *= -1
                out.append(z)
    return np.stack(out)


print(f'{"variant":>16} | {"rho all":>7} {"rho near":>8} {"sameNN":>6} {"10NN ovl":>8}')
for name, a in (('cur', 1), ('ang', .5), ('ang', 1), ('ang', 2), ('cur+ang', .5), ('cur+ang', 1), ('cur+ang', 2),
                ('corners', 1), ('corners+c', 1)):
    f = variant(name, a)
    ims = [f(s, sq) for s, sq in P]; odd = ims[0][1]
    I = np.stack([im for im, _ in ims]).reshape(K, -1)
    SY = np.stack([syms(im, odd).reshape(8, -1) for im, _ in ims])
    F = np.empty((K, K))
    for i in range(K):
        F[i] = np.abs(SY[None, i] - I[:, None]).sum(-1).min(-1)
    F = np.minimum(F, F.T); np.fill_diagonal(F, np.inf); fv = F[iu]
    ovl = np.mean([len(set(x) & set(y)) / 10 for x, y in zip(EK, np.argsort(F, 1)[:, :10])])
    lab = name + (f' a={a}' if 'ang' in name else '')
    print(f'{lab:>16} | {spearmanr(ex, fv)[0]:7.3f} {spearmanr(ex[near], fv[near])[0]:8.3f} '
          f'{np.mean(E.argmin(1) == F.argmin(1)):6.0%} {ovl:8.0%}', flush=True)
