#!/usr/bin/env python3
"""Standalone interactive mock-up of the s(110) connectivity map (10-11): hover a dot to see the packing, click to hold.
Data embedded (opens from file://).  Layouts from sweep.npz; island labels as render.py.

  ./landscape/mapenv.sh python3 landscape/mock.py [--out runs/landscape110/map/mock.html] NAME ...   (run from packer/)
"""
import argparse, collections, json, os, sys
import numpy as np
from scipy.sparse.csgraph import connected_components
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maps import RUN, REC, TOP
OUT = f'{RUN}/map'


def short(f): return f'{f[0]}/{f[1]}/{f[2]}/{f[3]}' + (f'/{f[4]}' if f[4] else '')


def norm(E):
    E = E - (E.max(0) + E.min(0)) / 2; w = (E.max(0) - E.min(0)).max(); return E / w + .5     # equal aspect in [0,1]^2


def labels(E, fam, kmin=4, link=.04):
    out = []
    for f, n in collections.Counter(fam).items():
        if n < kmin: continue
        m = np.array([i for i, g in enumerate(fam) if g == f]); P = E[m]
        k, lab = connected_components(np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)) < link, directed=False)
        for c in range(k):
            q = P[lab == c]
            if len(q) >= kmin:
                out.append([round(float(np.median(q[:, 0])), 4), round(float(q[:, 1].max()), 4), short(f)])
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=f'{OUT}/mock.html'); ap.add_argument('names', nargs='+')
    a = ap.parse_args()
    M = json.load(open(f'{OUT}/minima.json')); Z = np.load(f'{OUT}/sweep.npz')
    fam = [tuple(m['fam']) for m in M]
    lay = []
    for nm in a.names:
        E = norm(Z[nm]); lay.append(dict(name=nm.replace('_', ' '), xy=[round(v, 4) for v in E.ravel().tolist()], labels=labels(E, fam)))
    data = dict(rec=min(float(m['S']) for m in M), top=TOP, minima=[dict(S=m['S'], fam=short(tuple(m['fam'])), sq=m['sq']) for m in M], layouts=lay)
    html = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mock_template.html')).read()
    open(a.out, 'w').write(html.replace('/*DATA*/null', json.dumps(data, separators=(',', ':'))))
    print(a.out, os.path.getsize(a.out) // 1024, 'KB')


if __name__ == '__main__':
    main()
