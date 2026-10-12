"""Cheap configuration fingerprint + distance (10-11; for an in-explorer "far from everything known" parent rule).

Exact matching distance (pk.crossing.match_distance) costs ~0.2 s per pair: too slow to compare each new archive entry
with ~1000 others.  Fingerprint = Gaussian-smoothed G x G images of the square centres (coordinates / s), one channel for
axis squares and one for tilted squares (> 3 deg off axis).  Channel images are invariant to which band is "L" or "B", so
the 8 box symmetries are plain image flips / transposes.  Distance = min over the 8 symmetries of the L1 image difference / n.
Calibrated against exact matching distances in count_sep.py output (`fingerprint.py --calibrate runs/count_sep/k150.json`).
"""
from __future__ import annotations
import numpy as np

G = 24
SIG = 0.6          # in grid cells


def image(s, sq, g=G, sig=SIG):
    P = np.asarray([(x, y, t) for x, y, t in sq], float)
    a = np.mod(P[:, 2], 90.0); off = np.minimum(a, 90.0 - a)
    u = np.clip(P[:, 0] / s, 0, 1) * g - 0.5; v = np.clip(P[:, 1] / s, 0, 1) * g - 0.5
    c = np.arange(g)
    wx = np.exp(-0.5 * ((c[None, :] - u[:, None]) / sig) ** 2)      # n x g
    wy = np.exp(-0.5 * ((c[None, :] - v[:, None]) / sig) ** 2)
    wx /= wx.sum(1, keepdims=True); wy /= wy.sum(1, keepdims=True)
    tilt = off > 3.0
    return np.stack([(wy[~tilt].T @ wx[~tilt]), (wy[tilt].T @ wx[tilt])])     # 2 x g x g (rows = y)


def syms(im):
    out = []
    for t in (im, np.transpose(im, (0, 2, 1))):
        for fy in (False, True):
            for fx in (False, True):
                z = t[:, ::-1, :] if fy else t
                out.append(z[:, :, ::-1] if fx else z)
    return np.stack(out)          # 8 x 2 x g x g


def dist(A, B_syms, n):
    """A: one image; B_syms: syms() of the other.  Min over symmetries of L1 / n."""
    return float(np.abs(B_syms - A[None]).sum(axis=(1, 2, 3)).min() / n)


def nearest(A, bank, n):
    """Distance from image A to every image in bank (k x 2 x g x g), min over A's 8 symmetries: vectorised."""
    S = syms(A)                                              # 8 x ...
    if len(bank) == 0:
        return np.array([])
    d = np.abs(bank[:, None] - S[None]).sum(axis=(2, 3, 4)) / n     # k x 8
    return d.min(1)


def calibrate(path):
    import json, sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    from pk.store import Store
    from pk import refset
    d = json.load(open(path)); st = Store()
    R = {m['key']: m for m in refset.load(d['ref'])['minima']}
    ims = []
    for k in d['keys']:
        p = st.get(R[k]['rep']); ims.append(image(float(p.s), [tuple(map(float, r)) for r in p.sq]))
    n = len(st.get(R[d['keys'][0]]['rep']).sq)
    SY = [syms(im) for im in ims]
    ex = np.array([x for _, _, x in d['pairs']])
    fp = np.array([dist(ims[i], SY[j], n) for i, j, _ in d['pairs']])
    from scipy.stats import spearmanr, pearsonr
    print(f'{len(ex)} pairs: Spearman {spearmanr(ex, fp)[0]:.3f}, Pearson {pearsonr(ex, fp)[0]:.3f}')
    # nearest-neighbour agreement: does the fingerprint's NN equal the exact NN?
    K = len(ims); E = np.full((K, K), np.inf); F = np.full((K, K), np.inf)
    for (i, j, x), f in zip(d['pairs'], fp):
        E[i, j] = E[j, i] = x; F[i, j] = F[j, i] = f
    ennd = E.min(1); fnnd = F.min(1)
    print(f'NN distance per minimum: Spearman exact vs fingerprint {spearmanr(ennd, fnnd)[0]:.3f}; '
          f'same NN {np.mean(E.argmin(1) == F.argmin(1)):.0%}; top-quartile overlap '
          f'{len(set(np.argsort(-ennd)[:K // 4]) & set(np.argsort(-fnnd)[:K // 4])) / (K // 4):.0%}')
    for q in (10, 50, 90):
        print(f'  exact p{q} {np.percentile(ex, q):.4f}  fingerprint p{q} {np.percentile(fp, q):.4f}')


def bank(holdout, out):
    """Images of a holdout's visible minima (representatives), for explore.py --fp-bank."""
    import json, sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    from pk.store import Store
    from pk import refset
    H = json.load(open(f'{refset.ROOT}/{holdout}/holdout.json')); R = refset.load(H['ref']); st = Store()
    vis = set(H['visible']); ims = []; keys = []
    for m in R['minima']:
        if m['key'] in vis and m['rep']:
            p = st.get(m['rep']); ims.append(image(float(p.s), [tuple(map(float, r)) for r in p.sq])); keys.append(m['key'])
    np.savez_compressed(out, images=np.stack(ims), keys=np.array(keys))
    print(f'{out}: {len(ims)} visible minima of {holdout}')


if __name__ == '__main__':
    import sys
    if '--calibrate' in sys.argv:
        calibrate(sys.argv[sys.argv.index('--calibrate') + 1])
    elif '--bank' in sys.argv:
        i = sys.argv.index('--bank'); bank(sys.argv[i + 1], sys.argv[i + 2])
    else:
        raise SystemExit('fingerprint.py --calibrate COUNT_SEP.json | --bank HOLDOUT OUT.npz')
