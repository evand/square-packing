#!/usr/bin/env python3
"""One-square reduction of R_4 (2026-10-21).
phi(m, H) = min over a square b right of a left shelf (tilt, height, band) of its support in direction (1, m).
G(m) = 2 phi(m) - m H is a lower bound for L whenever the b-c separating line has slope m."""
import math, sys
import numpy as np
from scipy.optimize import differential_evolution, minimize, minimize_scalar
from relax import minx_window, kk
from frontier4 import Lt, F

def supp_dir(t, m):
    C, S = math.cos(t), math.sin(t)
    return (abs(C + m * S) + abs(-S + m * C)) / 2

def obj(v, m, H):
    t, uy, a = np.clip(v, 0, [math.pi / 2, 1, 1])
    k = kk(t)
    if 2 * k > H: return 10.0
    y = k + uy * (H - 2 * k); lo = a * (H - 1)
    mm = minx_window(t, lo - y, lo + 1 - y)
    x = k if mm is None else max(k, 1 - mm)
    return x + m * y + supp_dir(t, m)

def phi(m, H):
    best = math.inf
    for s in range(4):
        r = differential_evolution(obj, [(0, math.pi / 2), (0, 1), (0, 1)], args=(m, H), seed=s, tol=1e-13,
                                   maxiter=2000, polish=False)
        r2 = minimize(obj, r.x, args=(m, H), method='Nelder-Mead', options={'xatol': 1e-13, 'fatol': 1e-15})
        best = min(best, r.fun, r2.fun)
    return best

if __name__ == '__main__':
    H = float(sys.argv[1]) if len(sys.argv) > 1 else 1.93
    f, ts = F(H)
    print(f'H={H}  min(4,F)={min(4, f):.9f}  t*={ts:.3f}  m*=-cot t*={-1 / math.tan(math.radians(ts)):.4f}')
    for m in [0, -0.2, -0.5, -0.8, -1.0, -1.1, -1.15, -1 / math.tan(math.radians(ts)), -1.25, -1.4, -1.7, -2.5, -5]:
        g = 2 * phi(m, H) - m * H
        tm = math.atan2(1, -m) if m < 0 else math.pi / 2
        lt = Lt(tm, H) if 0 < tm < math.pi / 2 else float('inf')
        print(f'm={m:8.4f}  G={g:.9f}  G-min(4,F)={g - min(4, f):+.2e}  t_m={math.degrees(tm):6.2f}  '
              f'G-min(4,Lt(t_m))={g - min(4, lt):+.2e}', flush=True)
