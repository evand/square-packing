#!/usr/bin/env python3
"""The page's map layout (10-11, Evan's baseline): UMAP on fingerprint distance ** 0.5, 250 neighbours, min_dist 1,
spread 2; of seeds 0-4 the one filling the most cells of a 12 x 12 square frame (the rule that picked seed 3 on the
740-minimum snapshot).  Writes map/layout.json: seed, fill, xy normalised to [0,1]^2 (equal aspect), in minima.json order.

  ./landscape/mapenv.sh python3 landscape/layout.py        (run from packer/, after mapdata.py)"""
import json, os, sys, warnings
import numpy as np, umap
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maps import RUN
from sweep_sep import fillsq, metrics
P, K, MD, SP, SEEDS = .5, 250, 1.0, 2.0, range(5)
OUT = f'{RUN}/map'


def main():
    D = np.load(f'{OUT}/fpD.npy'); best = None
    for seed in SEEDS:
        E = umap.UMAP(n_neighbors=K, min_dist=MD, spread=SP, metric='precomputed', random_state=seed).fit_transform(D ** P)
        f = fillsq(E); print(f'seed {seed}: fill {f:.0%}', flush=True)
        if best is None or f > best[0]: best = (f, seed, E)
    f, seed, E = best
    E = E - (E.max(0) + E.min(0)) / 2; E = E / (E.max(0) - E.min(0)).max() + .5
    json.dump(dict(p=P, k=K, min_dist=MD, spread=SP, seed=seed, fill=f, xy=np.round(E, 5).tolist()), open(f'{OUT}/layout.json', 'w'))
    np.savez(f'{OUT}/sweep_page.npz', **{f'umap_p={P}_k={K}_md={MD:g}_spread={SP:g}_seed={seed}': E})
    print(f'-> layout.json: seed {seed}, fill {f:.0%}, n {len(E)}')


if __name__ == '__main__':
    main()
