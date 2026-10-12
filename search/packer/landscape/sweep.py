#!/usr/bin/env python3
"""Parameter sweep for the hover map (10-11): t-SNE / UMAP on fingerprint distance ** p.  Scores per layout:
crowd = share of points with nearest neighbour < 0.3 % of the plot width; fill = share of 12 x 12 bounding-box cells
occupied (principal axis aligned; reported by --fill); trust = sklearn trustworthiness (k = 10) vs the
untransformed distance; glob = Spearman(layout distance, D) on all pairs.  Saves every layout to sweep.npz.

  ./landscape/mapenv.sh python3 landscape/sweep.py [p,k,md,spread ...]   (run from packer/; full grid a few minutes)
"""
import json, os, sys, time, warnings
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
from sklearn.manifold import TSNE, trustworthiness
import umap
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maps import RUN
OUT = f'{RUN}/map'
D = np.load(f'{OUT}/fpD.npy'); iu = np.triu_indices(len(D), 1)


def fill(E, g=12):
    """Share of g x g cells of the bounding box (after aligning the principal axis) holding at least one point."""
    E = E - E.mean(0); E = E @ np.linalg.eigh(np.cov(E.T))[1]
    u = (E - E.min(0)) / (E.max(0) - E.min(0) + 1e-12) * (g - 1e-9)
    return len(set(map(tuple, u.astype(int)))) / g ** 2


def score(E):
    w = (E.max(0) - E.min(0)).max(); d = np.sqrt(((E[:, None] - E[None]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
    crowd = float((d.min(1) / w < .003).mean())
    tr = trustworthiness(D, E, n_neighbors=10, metric='precomputed')
    gl = spearmanr(pdist(E), D[iu])[0]
    return crowd, tr, gl


def grid():
    runs = {}
    for p in (1, .5, .3):
        Dp = D ** p
        for perp in (30, 100, 300):
            for ee in (12, 4):
                runs[f'tsne p={p} perp={perp} ee={ee}'] = lambda Dp=Dp, perp=perp, ee=ee: TSNE(
                    perplexity=perp, early_exaggeration=ee, metric='precomputed', init='random', random_state=0).fit_transform(Dp)
        for k in (30, 100, 250):
            for md, sp in ((.5, 1), (1, 2), (2, 4)):
                runs[f'umap p={p} k={k} md={md} spread={sp}'] = umap_run(p, k, md, sp)
    return runs


def umap_run(p, k, md, sp, rep=1.0, neg=5):
    return lambda: umap.UMAP(n_neighbors=k, min_dist=md, spread=sp, repulsion_strength=rep, negative_sample_rate=neg,
                             metric='precomputed', random_state=0).fit_transform(D ** p)


def report_fill():
    Es = np.load(f'{OUT}/sweep.npz'); res = {r[0].replace(' ', '_'): r for r in json.load(open(f'{OUT}/sweep.json'))}
    for k in Es.files:
        if k.startswith('umap'):
            r = res[k]; print(f'{r[0]:44s} fill {fill(Es[k]):5.1%} crowd {r[1]:5.1%} trust {r[2]:.3f} glob {r[3]:+.2f}')


def main():
    """--fill: fill score for every saved UMAP layout.No arguments: the full grid (overwrites).  Arguments 'p,k,md,spread[,repulsion[,negative_sample_rate]]' ...: those UMAP runs, appended."""
    if sys.argv[1:] == ['--fill']:
        return report_fill()
    if len(sys.argv) > 1:
        runs = {}
        for a in sys.argv[1:]:
            f = a.split(','); p, k, md, sp = f[:4]; rep = float(f[4]) if len(f) > 4 else 1.0; neg = int(f[5]) if len(f) > 5 else 5
            name = f'umap p={p} k={k} md={md} spread={sp}' + (f' rep={f[4]}' if rep != 1 else '') + (f' neg={neg}' if neg != 5 else '')
            runs[name] = umap_run(float(p), int(k), float(md), float(sp), rep, neg)
        Es = dict(np.load(f'{OUT}/sweep.npz')); res = json.load(open(f'{OUT}/sweep.json'))
    else:
        runs, Es, res = grid(), {}, []
    for name, f in runs.items():
        t = time.time(); E = f(); Es[name.replace(' ', '_')] = E
        c, tr, gl = score(E); res = [r for r in res if r[0] != name] + [(name, c, tr, gl)]
        print(f'{name:44s} fill {fill(E):5.1%} crowd {c:5.1%} trust {tr:.3f} glob {gl:+.2f}  ({time.time() - t:.0f} s)', flush=True)
    np.savez(f'{OUT}/sweep.npz', **Es)
    json.dump(res, open(f'{OUT}/sweep.json', 'w'))


if __name__ == '__main__':
    main()
