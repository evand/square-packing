#!/usr/bin/env python3
"""break_h1_cells: for a jsonl of recorded optima (lab, k, eps, delta, pattern, z), take the
best configuration per (label-family, k, eps), re-solve it at every eps by continuation (same
angle pattern, tilts rescaled), and report delta*, the CW / H / HP / H1 / H2 ladder, the
certificate depth, the weight-ladder depth and the fitted exponent p in delta* = -gamma eps^p.

Imports only; nothing in search/ is modified."""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import break_h1 as B                                             # noqa: E402
import s6skel                                                    # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('--T', type=float, required=True)
    ap.add_argument('--eps', default='0.25,0.5,1,2,5,10')
    ap.add_argument('--top', type=int, default=3)
    ap.add_argument('--maxchains', type=int, default=200)
    ap.add_argument('--tol', type=float, default=1e-10)
    ap.add_argument('--njit', type=int, default=40)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    T = a.T
    epss = sorted([float(s) for s in a.eps.split(',')], reverse=True)
    rs = [json.loads(line) for line in open(a.file)]
    rs = [r for r in rs if r.get('delta') is not None and r.get('cfree', True)]
    byk = {}
    for r in rs:
        key = (r['k'], r['eps'])
        byk.setdefault(key, []).append(r)
    # seeds: per k, the record with the largest delta*/eps^2 (the least degenerate branch)
    seeds = {}
    for (k, e), lst in byk.items():
        b = max(lst, key=lambda r: r['delta'] / math.radians(r['eps']) ** 2)
        if k not in seeds or (b['delta'] / math.radians(b['eps']) ** 2
                              > seeds[k]['delta'] / math.radians(seeds[k]['eps']) ** 2):
            seeds[k] = b
    rng = np.random.default_rng(23)
    out = open(a.out, 'w')
    for k in sorted(seeds):
        base = seeds[k]
        z0 = np.array(base['z'])
        n = (len(z0) - 1) // 3
        pat = [tuple(c) for c in base['pattern']] if base.get('pattern') else None
        X0, Y0, TH0 = z0[1::3], z0[2::3], z0[3::3]
        print(f"\n## T={T:g} n={n} k={k}  seed {base['lab']} at eps={base['eps']:g} "
              f"(delta* {base['delta']:+.4e})", flush=True)
        prev = None
        rows = []
        for epsdeg in epss:
            eps = math.radians(epsdeg)
            th = [(0.0 if abs(s6skel.norm_tilt(t)) < 1e-9 else (eps if t > 0 else -eps))
                  for t in TH0]
            seedsXY = [(X0.copy(), Y0.copy())]
            bx, by = prev if prev is not None else (X0, Y0)
            seedsXY.append((bx.copy(), by.copy()))
            for _ in range(a.njit):
                seedsXY.append((bx + rng.normal(0, 0.03, n), by + rng.normal(0, 0.03, n)))
            v, z, hit = B.delta_star_at(th, T, rng, seeds=seedsXY, limit=0, jitters=0, extra=0,
                                        pattern=pat, k=(k if pat else 0))
            if z is None:
                print(f'  eps={epsdeg:g}: infeasible', flush=True)
                continue
            cf = B.chainfree_ok(z, n, T, v, range(k)) if pat else None
            prev = (z[1::3].copy(), z[2::3].copy())
            o = B.bounds(z, n, T, v, tol=a.tol, maxchains=a.maxchains, eps_rad=eps,
                         kmax=7, do_cyc=False, nch2=1, cap2=1500)
            rows.append((epsdeg, v))
            print(f"  eps={epsdeg:6g}  delta*={v:+.8e}  /e^2 {-v/eps**2:.5f} "
                  f"/e^3 {-v/eps**3:.4f} /e^4 {-v/eps**4:.3f}  cfree={cf} hit={hit}", flush=True)
            print(f"      CW {o['CW']:+.4e}  H {o['H']:+.4e}  HP {o['HP']:+.4e}  "
                  f"H1 {o['H1']:+.4e}  H2 {o['H2']:+.4e}  F {o['F']:+.6e}  "
                  f"cdepth={o['cdepth']} wdepth={o['wdepth']} mu={o['mu']} "
                  f"nchains={o['nchains']} ({o['secs']:.0f}s)", flush=True)
            print('      wlevels ' + ' '.join(f'{l}:{c}x{w:.3e}' for l, c, w in o['wlevels']),
                  flush=True)
            rec = dict(lab=base['lab'], T=T, n=n, k=k, eps=epsdeg, delta=float(v),
                       cfree=(None if cf is None else bool(cf)), hit=int(hit),
                       pattern=base.get('pattern'), z=[float(q) for q in z])
            for key in ('CW', 'H', 'HP', 'H1', 'H2', 'F', 'cdepth', 'wdepth', 'wlevels',
                        'mu', 'E', 'V', 'cls', 'nchains', 'supp'):
                rec[key] = o[key]
            out.write(json.dumps(rec) + '\n'); out.flush()
        if len(rows) >= 2:
            p, g = B.fit_exponent([e for e, _v in rows], [v for _e, v in rows])
            ps = [(rows[i][0], math.log(rows[i][1] / rows[i + 1][1])
                   / math.log(math.radians(rows[i][0]) / math.radians(rows[i + 1][0])))
                  for i in range(len(rows) - 1) if rows[i][1] < 0 and rows[i + 1][1] < 0]
            print(f"  FIT  delta* = -{g:.4f} eps^{p:.3f}   local slopes "
                  + ' '.join(f'{e:g}:{q:.3f}' for e, q in ps), flush=True)
    out.close()


if __name__ == '__main__':
    main()
