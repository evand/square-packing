#!/usr/bin/env python3
"""M(tb, tc, H): the chain relaxation minimised over heights and bands at fixed tilts (2026-10-21)."""
import math, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import differential_evolution, minimize
from relax import Lmin
from frontier4 import Lt, F

def M(args):
    tb, tc, H = args
    f = lambda w: Lmin([tb, tc, *np.clip(w, 0, 1)], H)
    best = math.inf
    for s in range(3):
        r = differential_evolution(f, [(0, 1)] * 4, seed=s, tol=1e-12, maxiter=1500, polish=False,
                                   init=np.vstack([np.random.default_rng(s).random((60, 4)),
                                                   [[1, 0, 0, 1], [0, 1, 1, 0]]]))
        r2 = minimize(f, r.x, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-14})
        best = min(best, r.fun, r2.fun)
    return tb, tc, best

if __name__ == '__main__':
    H = float(sys.argv[1]) if len(sys.argv) > 1 else 1.93
    ts = np.radians([0, 10, 20, 30, 35, 38, 40.5, 43, 45, 47, 49.5, 52, 55, 60, 70, 80, 89.99])
    jobs = [(a, b, H) for a in ts for b in ts]
    with ProcessPoolExecutor(8) as ex:
        res = list(ex.map(M, jobs, chunksize=4))
    f, tstar = F(H)
    print(f'H={H} F={f:.6f} t*={tstar:.2f}')
    print('tb\\tc ' + ' '.join(f'{math.degrees(t):6.1f}' for t in ts))
    k = 0
    for a in ts:
        row = []
        for b in ts:
            row.append(res[k][2]); k += 1
        print(f'{math.degrees(a):6.1f} ' + ' '.join(f'{v - f:6.3f}' for v in row))
