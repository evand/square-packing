#!/usr/bin/env python3
"""Which store-best packings are invariant under a D4 element (10-18, stage 0 of the symmetric-subspace test).
Elements: Md (x,y)->(y,x), Ma (x,y)->(s-y,s-x), C2, Mx, My.  Angle maps mod 90: mirrors th -> -th, C2 keeps th.
  symscan.py [--lo 30 --hi 324] [--tol 1e-6]   -> runs/symsub/symscan.json"""
import sys, os, json, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.spatial import cKDTree

EL = {'Md': (lambda x, y, s: (y, x), -1), 'Ma': (lambda x, y, s: (s - y, s - x), -1),
      'C2': (lambda x, y, s: (s - x, s - y), 1), 'Mx': (lambda x, y, s: (s - x, y), -1),
      'My': (lambda x, y, s: (x, s - y), -1)}


def angdiff(a, b):
    d = (a - b) % 90
    return np.minimum(d, 90 - d)


def invariant(sq, s, el, tol):
    f, sg = EL[el]
    x2, y2 = f(sq[:, 0], sq[:, 1], s)
    t2 = (sg * sq[:, 2]) % 90
    tr = cKDTree(sq[:, :2])
    d, k = tr.query(np.c_[x2, y2])
    return bool((d < tol).all() and (angdiff(t2, sq[k, 2]) < tol * 100).all()), float(np.mean(d < 1e-3))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--lo', type=int, default=30); ap.add_argument('--hi', type=int, default=324)
    ap.add_argument('--tol', type=float, default=1e-6)
    a = ap.parse_args()
    from pk.store import Store
    st = Store(); pb = st.public_best()
    out = {}
    for n in range(a.lo, a.hi + 1):
        if n not in pb:
            continue
        p = st.get(pb[n][1]); sq = p.sq
        ntilt = int((np.abs(((sq[:, 2] + 45) % 90) - 45) > 1e-6).sum())
        r = {el: invariant(sq, p.s, el, a.tol) for el in EL}
        out[n] = dict(s=p.s, pid=pb[n][1], ntilt=ntilt, inv=[el for el in EL if r[el][0]],
                      match={el: r[el][1] for el in EL})
    os.makedirs('runs/symsub', exist_ok=True)
    json.dump(out, open('runs/symsub/symscan.json', 'w'), indent=0)
    md = [n for n, r in out.items() if 'Md' in r['inv'] and r['ntilt'] and abs(r['s'] - round(r['s'])) > 1e-9]
    print('Md-invariant, tilted, non-integer:', len(md), md)
    near = [n for n, r in out.items() if not r['inv'] and r['ntilt'] and max(r['match'].values()) >= 0.95]
    print('asymmetric but >=95% mirror-matched:', len(near))


if __name__ == '__main__':
    main()
