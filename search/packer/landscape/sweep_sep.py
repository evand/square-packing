#!/usr/bin/env python3
"""Fill vs family separation (10-11): UMAP k = 250 over p x (min_dist, spread) x 5 seeds.
fill = share of 12 x 12 cells of a square frame occupied; sep = median distance to the nearest other-family point /
median distance to the nearest same-family point; crowd as sweep.py; glob = Spearman(layout distance, D).
Saves every layout to sep.npz (key 'p_md_sp_seed').

  ./landscape/mapenv.sh python3 landscape/sweep_sep.py        (run from packer/; ~2 min)"""
import json, os, sys, warnings
import numpy as np, umap
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maps import RUN
OUT = f'{RUN}/map'
D = np.load(f'{OUT}/fpD.npy'); iu = np.triu_indices(len(D), 1)
fam = [tuple(m['fam']) for m in json.load(open(f'{OUT}/minima.json'))]
same = np.array([[a == b for b in fam] for a in fam])


def fillsq(E, g=12):
    E = E - E.mean(0); w = (E.max(0) - E.min(0)).max(); c = (E.max(0) + E.min(0)) / 2
    return len(set(map(tuple, (((E - c) / w + .5) * (g - 1e-9)).astype(int)))) / g ** 2


def metrics(E):
    d = np.sqrt(((E[:, None] - E[None]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
    w = (E.max(0) - E.min(0)).max()
    own = np.where(same, d, np.inf).min(1); oth = np.where(~same, d, np.inf).min(1)
    ok = np.isfinite(own)
    return fillsq(E), np.median(oth) / np.median(own[ok]), (d.min(1) / w < .003).mean(), spearmanr(pdist(E), D[iu])[0]


def main():
    Es = {}
    print(f'{"p":>4} {"md":>4} {"sp":>3} |  fill   sep  crowd  glob   (mean over 5 seeds; fill range)')
    for p in (.5, .6, .7, .85):
        for md, sp in ((2, 4), (1, 4), (.5, 4), (1, 2)):
            R = []
            for seed in range(5):
                E = umap.UMAP(n_neighbors=250, min_dist=md, spread=sp, metric='precomputed', random_state=seed).fit_transform(D ** p)
                Es[f'{p}_{md}_{sp}_{seed}'] = E; R.append(metrics(E))
            R = np.array(R); m = R.mean(0)
            print(f'{p:4} {md:4} {sp:3} | {m[0]:5.0%} {m[1]:5.2f} {m[2]:5.1%} {m[3]:5.2f}   ({R[:, 0].min():.0%}-{R[:, 0].max():.0%})', flush=True)
    np.savez(f'{OUT}/sep.npz', **Es)


if __name__ == '__main__':
    main()
