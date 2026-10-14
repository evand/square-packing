#!/usr/bin/env python3
"""Family map of the frontier on the triangle (10-14; Evan: "blur jlevy's triangle view and squint: which packings are
out of place?").  Triangle position: k = ceil(sqrt n), c = k^2 - n (column; c = 0 is the full grid).

Per store-best packing: tilted count T, angle families (pk.bands.families: sizes + mean angle), bands, and the
thinness of the tilted set (smallest / largest principal spread of the tilted centres).  Code (heuristic, for
squinting, not a theorem):
  G  grid (s = k, or no tilted squares)
  D  45-degree block (>= 80 % of tilted within 1 degree of 45, thick)
  S  45-degree strip(s) (same, thin: thickness / length < 0.25)
  V  one low-angle band (tilted angles away from 45, one band)
  2  two or more low-angle bands
  J  jumble (>= 3 angle families, or tilted angles spread with no dominant family)

  python3 familymap.py [--nmax 324] [--json out.json]
"""
from __future__ import annotations
import argparse, json, math, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def classify(n, s, sq):
    from pk import bands as B
    sq = np.asarray(sq, float)
    k = math.isqrt(n - 1) + 1
    t = B.tilted(sq)
    T = int(t.sum())
    a = sq[t, 2] % 90.0
    dev45 = np.abs(a - 45.0) if T else np.array([])
    bs, _ = B.bands(sq)
    fam = B.families(sq)
    fam_sizes = fam[0] if isinstance(fam, tuple) else fam
    info = dict(n=n, k=k, c=k * k - n, s=float(s), T=T, bands=list(map(int, bs)),
                fams=[int(x) for x in fam_sizes] if fam_sizes is not None else [])
    if s >= k - 1e-9 or T < 5:
        code = 'G'
    else:
        P = sq[t, :2]
        ev = np.sort(np.linalg.eigvalsh(np.cov(P.T))) if T >= 3 else np.array([0, 1])
        thin = math.sqrt(max(ev[0], 0) / max(ev[1], 1e-12))
        info['thin'] = round(thin, 3)
        info['a_med'] = round(float(np.median(np.minimum(a, 90 - a))), 2)
        if np.mean(dev45 < 1.0) >= 0.8:
            code = 'S' if thin < 0.25 else 'D'
        elif len(info['fams']) >= 3:
            code = 'J'
        elif len(bs) >= 2:
            code = '2'
        else:
            code = 'V'
    info['code'] = code
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmin', type=int, default=2)
    ap.add_argument('--nmax', type=int, default=324)
    ap.add_argument('--json')
    a = ap.parse_args()
    from pk.store import Store
    st = Store()
    rows = []
    for n in range(a.nmin, a.nmax + 1):
        r = st.db.execute("SELECT id, s FROM packing WHERE n=? AND status != 'infeasible' AND (pen IS NULL OR pen <= 1e-9)"
                          " ORDER BY s LIMIT 1", (n,)).fetchone()
        if r is None:
            continue
        info = classify(n, r['s'], st.get(r['id']).sq)
        info['id'] = r['id']
        rows.append(info)
    if a.json:
        json.dump(rows, open(a.json, 'w'), indent=0)
    by = {(r['k'], r['c']): r for r in rows}
    K = max(r['k'] for r in rows)
    print('rows k, columns c = k^2 - n (c = 0 at the right edge)')
    for k in range(2, K + 1):
        cells = []
        for c in range(2 * k - 2, -1, -1):
            r = by.get((k, c))
            cells.append(r['code'] if r else ' ')
        print(f'k={k:2d} {(k - 1) ** 2 + 1:4d}..{k * k:4d}  ' + ' ' * (2 * (K - k)) + ' '.join(cells))
    print()
    for r in rows:
        if r['code'] not in 'G':
            print(f"{r['n']:4d} k{r['k']:2d} c{r['c']:2d} {r['code']} s {r['s']:.6f} T {r['T']:3d} fams {r['fams']} bands {r['bands']}"
                  f" thin {r.get('thin')} a {r.get('a_med')}")


if __name__ == '__main__':
    main()
