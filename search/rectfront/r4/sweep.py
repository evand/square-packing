#!/usr/bin/env python3
"""Sweep of the chain relaxation vs R_4's frontier (2026-10-21).  Populations seeded with the frontier
configuration (both mirror forms) plus random points; reports relaxed min - min(4, F(H))."""
import math, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import differential_evolution, minimize
from relax import Lmin
from frontier4 import F

UB = np.array([math.pi / 2, math.pi / 2, 1, 1, 1, 1])

def one(H):
    f, tdeg = F(H)
    t = math.radians(tdeg)
    seeds = [[t, t, 1, 0, 0, 1], [math.pi / 2 - t, math.pi / 2 - t, 0, 1, 1, 0]]
    best = (math.inf, None)
    for s in range(8):
        rng = np.random.default_rng(s)
        init = np.vstack([rng.random((200, 6)) * UB] + ([np.array(seeds)] if s % 2 == 0 else []))
        r = differential_evolution(Lmin, list(zip([0] * 6, UB)), args=(H,), init=init, seed=s, tol=1e-13,
                                   maxiter=4000, polish=False)
        r2 = minimize(Lmin, r.x, args=(H,), method='Nelder-Mead',
                      options={'xatol': 1e-13, 'fatol': 1e-15, 'maxiter': 40000})
        for v in (r.fun, r2.fun):
            if v < best[0]: best = (v, np.clip(r2.x if v == r2.fun else r.x, 0, UB))
    return H, best[0], min(4, f), best[1]

if __name__ == '__main__':
    Hs = [1.2, 1.5, 1.7, 1.8, 1.85, 1.87, 1.88] + list(np.round(np.linspace(1.8817, 1.9995, 25), 5))
    with ProcessPoolExecutor(8) as ex:
        for H, v, tgt, x in ex.map(one, Hs):
            print(f'H={H:.5f} relaxed={v:.12f} target={tgt:.12f} gap={v - tgt:+.2e} '
                  f'tb={math.degrees(x[0]):.2f} tc={math.degrees(x[1]):.2f} y=({x[2]:.2f},{x[3]:.2f}) lo=({x[4]:.2f},{x[5]:.2f})',
                  flush=True)
