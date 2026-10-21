#!/usr/bin/env python3
"""Does the chain relaxation reach R_4's frontier?  (2026-10-21)
Relaxed problem for strip height H < 2: squares b, c (tilts tb, tc, heights yb, yc), a left shelf band
[lo, lo+1] on b (points of b there have x >= 1), a right shelf band [lo2, lo2+1] on c (x <= L - 1), b and c
disjoint with x_b <= x_c.  Least L = xb_min + sigma(b, c) + (L - xc)_min; minimised over the 6 parameters."""
import math, sys
import numpy as np
from scipy.optimize import differential_evolution, minimize
from frontier4 import F

def verts(t):
    C, S = math.cos(t), math.sin(t)
    e, f = np.array([C, S]) / 2, np.array([-S, C]) / 2
    return [e + f, -e + f, -e - f, e - f]

def kk(t):
    return (abs(math.cos(t)) + abs(math.sin(t))) / 2

def minx_window(t, y0, y1):
    """min x of the centred square restricted to heights [y0, y1]; None if empty."""
    V = verts(t); best = None
    for i in range(4):
        p, q = V[i], V[(i + 1) % 4]
        cands = []
        if y0 <= p[1] <= y1: cands.append(p[0])
        for yy in (y0, y1):
            if (p[1] - yy) * (q[1] - yy) < 0:
                s = (yy - p[1]) / (q[1] - p[1]); cands.append(p[0] + s * (q[0] - p[0]))
        for x in cands:
            best = x if best is None else min(best, x)
    return best

def supp(t, n):
    C, S = math.cos(t), math.sin(t)
    return abs(n[0] * C + n[1] * S) + abs(-n[0] * S + n[1] * C)

def sigma(tb, tc, dy):
    best = math.inf
    for t in (tb, tc):
        C, S = math.cos(t), math.sin(t)
        for n in ((C, S), (-C, -S), (-S, C), (S, -C)):
            if n[0] > 1e-12:
                best = min(best, ((supp(tb, n) + supp(tc, n)) / 2 - n[1] * dy) / n[0])
    return best

def Lmin(v, H):
    v = np.clip(v, 0, [math.pi / 2, math.pi / 2, 1, 1, 1, 1])
    tb, tc, ub, uc, a, a2 = v
    kb, kc = kk(tb), kk(tc)
    if 2 * kb > H or 2 * kc > H: return 10.0
    yb = kb + ub * (H - 2 * kb); yc = kc + uc * (H - 2 * kc)
    lo = a * (H - 1); lo2 = a2 * (H - 1)
    m = minx_window(tb, lo - yb, lo + 1 - yb)
    xb = kb if m is None else max(kb, 1 - m)
    m2 = minx_window(-tc, lo2 - yc, lo2 + 1 - yc)   # mirror of c
    rc = kc if m2 is None else max(kc, 1 - m2)
    return xb + max(0.0, sigma(tb, tc, yc - yb)) + rc

if __name__ == '__main__':
    bounds = [(0, math.pi / 2), (0, math.pi / 2), (0, 1), (0, 1), (0, 1), (0, 1)]
    for H in [float(x) for x in sys.argv[1:]] or [1.85, 1.88, 1.9, 1.93, 1.96, 1.99]:
        best = None
        for seed in range(6):
            r = differential_evolution(Lmin, bounds, args=(H,), seed=seed, tol=1e-12, maxiter=3000,
                                       popsize=40, polish=False)
            r2 = minimize(Lmin, r.x, args=(H,), method='Nelder-Mead',
                          options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 20000})
            if best is None or r2.fun < best.fun: best = r2
        f, tdeg = F(H)
        tb, tc, ub, uc, a, a2 = np.clip(best.x, 0, [math.pi / 2, math.pi / 2, 1, 1, 1, 1])
        print(f'H={H}: relaxed {best.fun:.9f}  frontier min(4, {f:.9f})  gap {best.fun - min(4, f):+.2e}  '
              f'tb={math.degrees(tb):.2f} tc={math.degrees(tc):.2f} ub={ub:.3f} uc={uc:.3f} lo={a:.3f} lo2={a2:.3f}',
              flush=True)
