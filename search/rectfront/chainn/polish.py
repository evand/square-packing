#!/usr/bin/env python3
"""Continuous polish of the chain relaxation from the DP chain (all tilts, heights, both shelf bands free)."""
import sys, math
import numpy as np
from scipy.optimize import minimize
from dp import dp, minx_window
from targets import target
import paths as _p  # noqa  (reuse path extraction)

def sig(tb, tc, dy):
    best = math.inf
    for t in (tb, tc):
        C, S = math.cos(t), math.sin(t)
        for nx, ny in ((C, S), (-C, -S), (-S, C), (S, -C)):
            if nx > 1e-12:
                sb = abs(nx * math.cos(tb) + ny * math.sin(tb)) + abs(-nx * math.sin(tb) + ny * math.cos(tb))
                sc = abs(nx * math.cos(tc) + ny * math.sin(tc)) + abs(-nx * math.sin(tc) + ny * math.cos(tc))
                best = min(best, ((sb + sc) / 2 - ny * dy) / nx)
    return max(best, 0.0)

def Lrel(v, H, m):
    t = np.clip(v[:m], 0, math.pi / 2); u = np.clip(v[m:2 * m], 0, 1); lo, lo2 = np.clip(v[2 * m:], 0, 1) * (H - 1)
    k = (np.cos(t) + np.sin(t)) / 2
    if np.any(2 * k > H): return 1e3
    y = k + u * (H - 2 * k)
    mm = minx_window(t[0], lo - y[0], lo + 1 - y[0]); x = k[0] if mm is None else max(k[0], 1 - mm)
    for i in range(m - 1): x += sig(t[i], t[i + 1], y[i + 1] - y[i])
    mm = minx_window(-t[-1], lo2 - y[-1], lo2 + 1 - y[-1]); x += k[-1] if mm is None else max(k[-1], 1 - mm)
    return x

def polish(H, n, seed=0):
    out, (T, Y, Sg, St, En, back) = dp(H, n, 91, 41)
    V = back[n - 3]; j = int(np.argmin(V + En)); seq = [j]
    for r in range(n - 4, -1, -1):
        j = int(np.argmin(back[r] + Sg[:, j])); seq.append(j)
    seq.reverse(); m = n - 2
    t0 = T[seq]; k = (np.cos(t0) + np.sin(t0)) / 2; u0 = (Y[seq] - k) / (H - 2 * k)
    rng = np.random.default_rng(seed)
    best = None
    def bestlo(t, y):
        g = np.linspace(0, 1, 401)
        f = lambda a: (lambda mm: (k0 if mm is None else max(k0, 1 - mm)))(minx_window(t, a * (H - 1) - y, a * (H - 1) + 1 - y))
        k0 = (abs(math.cos(t)) + abs(math.sin(t))) / 2
        return g[int(np.argmin([f(a) for a in g]))]
    lo0 = np.array([bestlo(t0[0], Y[seq[0]]), bestlo(-t0[-1], Y[seq[-1]])])
    for rep in range(8):
        v0 = np.concatenate([t0, u0, lo0]) + (rng.normal(0, 0.01, 2 * m + 2) if rep else 0)
        r = minimize(Lrel, v0, args=(H, m), method='Powell', options={'xtol': 1e-12, 'ftol': 1e-14, 'maxfev': 200000})
        r = minimize(Lrel, r.x, args=(H, m), method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-15, 'maxfev': 200000})
        if best is None or r.fun < best.fun: best = r
    return out[n], best.fun, target(n, H)

if __name__ == '__main__':
    from concurrent.futures import ProcessPoolExecutor
    jobs = [(H, n) for H in (1.84, 1.87, 1.9, 1.95, 1.99) for n in (6, 7, 8, 9, 11, 12)]
    with ProcessPoolExecutor(8) as ex:
        for (H, n), (d, p, t) in zip(jobs, ex.map(polish, *zip(*jobs))):
            print(f"H={H} n={n}: DP {d:.6f}  polished {p:.9f}  target {t:.9f}  gap {p - t:+.2e}", flush=True)
