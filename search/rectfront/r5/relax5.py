#!/usr/bin/env python3
"""R_5, strip H < 2: is the chain relaxation exact? (2026-10-21)
Squares sorted a..e by centre x.  Kept: left shelf on b (band lo), right shelf on d (band lo2), b-c and c-d
disjoint (x order), optionally b-d disjoint.  Target min(5, 4 + 6 sqrt2 - 4H) (Friedman's S_4 pieces 4.1, 4.2)."""
import math, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import differential_evolution, minimize
sys.path.insert(0, '../r4')
from relax import minx_window, kk, sigma

UB = np.array([math.pi / 2] * 3 + [1] * 5)

def Lmin(v, H, bd=False):
    v = np.clip(v, 0, UB)
    tb, tc, td, ub, uc, ud, a, a2 = v
    ks = [kk(t) for t in (tb, tc, td)]
    if max(ks) * 2 > H: return 10.0
    yb, yc, yd = (k + u * (H - 2 * k) for k, u in zip(ks, (ub, uc, ud)))
    lo, lo2 = a * (H - 1), a2 * (H - 1)
    m = minx_window(tb, lo - yb, lo + 1 - yb)
    xb = ks[0] if m is None else max(ks[0], 1 - m)
    xc = xb + max(0.0, sigma(tb, tc, yc - yb))
    xd = xc + max(0.0, sigma(tc, td, yd - yc))
    if bd:
        xd = max(xd, xb + max(0.0, sigma(tb, td, yd - yb)))
    m2 = minx_window(-td, lo2 - yd, lo2 + 1 - yd)
    rd = ks[2] if m2 is None else max(ks[2], 1 - m2)
    return xd + rd

def target(H):
    return min(5.0, 4 + 6 * math.sqrt(2) - 4 * H)

def one(args):
    H, bd = args
    q = math.pi / 4
    seeds = np.array([[q, q, q, 1, 0, 1, 0, 0]])
    best = math.inf; bx = None
    for s in range(6):
        rng = np.random.default_rng(s)
        init = np.vstack([rng.random((240, 8)) * UB] + ([seeds] if s % 2 == 0 else []))
        r = differential_evolution(Lmin, list(zip([0] * 8, UB)), args=(H, bd), init=init, seed=s, tol=1e-13,
                                   maxiter=4000, polish=False)
        r2 = minimize(Lmin, r.x, args=(H, bd), method='Nelder-Mead',
                      options={'xatol': 1e-13, 'fatol': 1e-15, 'maxiter': 40000})
        for val, x in ((r.fun, r.x), (r2.fun, r2.x)):
            if val < best: best, bx = val, np.clip(x, 0, UB)
    return H, bd, best, bx

if __name__ == '__main__':
    Hs = [1.3, 1.6, 1.8, 1.85, 1.87, 1.88, 1.9, 1.93, 1.96, 1.99]
    with ProcessPoolExecutor(10) as ex:
        for H, bd, v, x in ex.map(one, [(H, bd) for bd in (False, True) for H in Hs]):
            t = target(H)
            print(f"{'b-d ' if bd else 'plain'} H={H:.3f} relaxed={v:.9f} target={t:.9f} gap={v - t:+.2e} "
                  f"tilts=({', '.join(f'{math.degrees(a):.1f}' for a in x[:3])}) y=({', '.join(f'{a:.2f}' for a in x[3:6])}) "
                  f"lo=({x[6]:.2f},{x[7]:.2f})", flush=True)
