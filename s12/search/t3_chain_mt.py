#!/usr/bin/env python3
"""t3_chain_mt: the architecture's literal claim, tested correctly.

`notes/proof-architecture.md` sec 0a item 14 / `notes/review-2026-09-20c.md`:
   "every configuration has a tight wall-to-wall chain of T squares, tilted or not,
    on which MT gives delta <= 0".

MT (`notes/bandcut-cost.md` sec 1) is valid ONLY for a chain whose links are separated along a
COMMON unit normal e = (cos phi, sin phi), cos phi > 0:

    delta <= MT = [ cos(phi) (T - ubar) + sin(phi) R - M ] / (k - 1 + 2 cos(phi))
    M = sum_i m(theta_{i+1} - theta_i),  ubar = (u_1 + u_k)/2,  R = (c_k - c_1)_y .

`search/bandcut_k.py chain_certificate` does NOT enforce the common normal: it builds the chain
from the x/y-DAG (which at tilt 45 deg is degenerate, |cos| = |sin|) and then evaluates B_T at the
MEAN tilt.  At 45 deg that produces values BELOW the true margin, i.e. invalid bounds (sec 5.2 of
`notes/t3-chain.md`).  This script enumerates chains link-normal by link-normal instead, so every
number it prints is a valid upper bound on delta.
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                    # noqa: E402

TOL = 1e-7


def sep_rows(z, n, lvl, tol=TOL):
    """-> dict (i, j) -> list of (nx, ny, m): every unit normal along which the ordered pair
    i -> j is separated at level lvl (j on the + side)."""
    X, Y, TH = np.array(z[1::3]), np.array(z[2::3]), np.array(z[3::3])
    out = defaultdict(list)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = min(i, j), max(i, j)
            D = TH[b] - TH[a]
            m = 0.5 + 0.5 * (abs(math.cos(D)) + abs(math.sin(D)))
            dd = np.array([X[j] - X[i], Y[j] - Y[i]])
            for o in (i, j):
                C, S = math.cos(TH[o]), math.sin(TH[o])
                for (d0, d1) in ((C, S), (-S, C)):
                    for sg in (1, -1):
                        nx, ny = sg * d0, sg * d1
                        if nx * dd[0] + ny * dd[1] - m >= lvl - tol:
                            out[(i, j)].append((nx, ny, m))
    return out


def best_mt(z, n, T, lvl, k=None, ntol=1e-6):
    """min over (common normal e, chain of k along e, axis) of MT.  A valid upper bound on delta.
    The `axis` is the wall pair used; `e` must have positive component along it."""
    k = k or int(round(T))
    X, Y, TH = np.array(z[1::3]), np.array(z[2::3]), np.array(z[3::3])
    u = np.abs(np.cos(TH)) + np.abs(np.sin(TH))
    S = sep_rows(z, n, lvl)
    normals = []
    for lst in S.values():
        for (nx, ny, m) in lst:
            if not any(abs(nx - a) < ntol and abs(ny - b) < ntol for (a, b) in normals):
                normals.append((nx, ny))
    best = (math.inf, None)
    for (nx, ny) in normals:
        for ax in ('x', 'y'):
            cphi = nx if ax == 'x' else ny          # component along the wall pair's axis
            sphi = ny if ax == 'x' else nx          # the transverse component
            if cphi <= 1e-9:
                continue
            adj = {i: [] for i in range(n)}
            for (i, j), lst in S.items():
                for (a, b, m) in lst:
                    if abs(a - nx) < ntol and abs(b - ny) < ntol:
                        adj[i].append((j, m))
                        break
            # every path on k vertices
            def walk(path, M):
                if len(path) == k:
                    i0, ik = path[0], path[-1]
                    perp = (Y if ax == 'x' else X)
                    R = float(perp[ik] - perp[i0]) * (1.0 if ax == 'x' else 1.0)
                    ub = 0.5 * (u[i0] + u[ik])
                    val = (cphi * (T - ub) + sphi * R - M) / (k - 1 + 2 * cphi)
                    yield val, (ax, list(path), round(math.degrees(math.atan2(ny, nx)), 3),
                                round(R, 5))
                    return
                for (j, m) in adj[path[-1]]:
                    if j not in path:
                        yield from walk(path + [j], M + m)
            for i in range(n):
                for val, info in walk([i], 0.0):
                    if val < best[0]:
                        best = (val, info)
    return best


def main():
    recs = []
    for f in sys.argv[1:]:
        recs += [json.loads(line) for line in open(f)]
    fams = defaultdict(list)
    print('# MT on chains with a COMMON separating normal.  Every value is a valid upper bound on')
    print('# delta for the configuration it is computed at (so it must be >= delta*).')
    bad = []
    invalid = 0
    for r in recs:
        z = np.array(r['z'])
        v, info = best_mt(z, 6, 3.0, r['delta'])
        r['MT'] = v
        r['MTinfo'] = info
        if v < r['delta'] - 1e-9:
            invalid += 1
        fams[r['fam']].append(r)
        if v > 1e-9:
            bad.append(r)
    print(f'# {"family":13s} {"N":>3s} {"median MT":>12s} {"max MT":>12s} {"MT<=0":>8s} '
          f'{"H<=0":>8s} {"no common-normal chain":>22s}')
    for fam in sorted(fams):
        rs = fams[fam]
        M = sorted(r['MT'] for r in rs)
        print(f'  {fam:13s} {len(rs):3d} {M[len(M)//2]:+12.4e} {M[-1]:+12.4e} '
              f'{sum(1 for r in rs if r["MT"] <= 1e-9):4d}/{len(rs):<3d} '
              f'{sum(1 for r in rs if r["H"] <= 1e-9):4d}/{len(rs):<3d} '
              f'{sum(1 for r in rs if r["MT"] == math.inf):22d}')
    print(f'\n# MT <= 0 at {len(recs)-len(bad)}/{len(recs)} sampled optima; '
          f'{invalid} values below delta* (must be 0 -- MT is an upper bound)')
    for r in sorted(bad, key=lambda r: -r['MT'])[:20]:
        tl = ' '.join(f'{math.degrees(s6skel.norm_tilt(t)):+6.1f}' for t in r['theta'])
        print(f"  {r['fam']:12s} d*={r['delta']:+.4e} MT={r['MT']:+.4e} H={r['H']:+.4e} "
              f"{r['MTinfo']}  tilts[{tl}]")


if __name__ == '__main__':
    main()
