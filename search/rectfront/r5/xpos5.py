#!/usr/bin/env python3
"""R_5 strip: X(c) = min over separating lines (t, sign) of B_b(t, sign) - m y + h_c(m)/2,
a lower bound on c's centre x given a, b to its left (2026-10-21).  Compare with lines along c's edges.
Line at angle d: direction (cos d, sin d), f = x + m y with m = -cot d  (d in (0,90): sign -1, t = d;
d in (90,180): sign +1, t = 180 - d)."""
import math
import numpy as np
from scipy.optimize import minimize_scalar
from split5 import Bb

def hsup(th, m):
    P, Q = math.cos(th), math.sin(th)
    return abs(P + m * Q) + abs(-Q + m * P)

def expr(t, s, th, y, H):
    m = s / math.tan(t)
    return Bb(t, s, H) - m * y + hsup(th, m) / 2

def X(th, y, H):
    best = (math.inf, None)
    for s in (-1, 1):
        ts = np.radians(np.linspace(0.05, 89.95, 3600))
        vals = [expr(t, s, th, y, H) for t in ts]
        i = int(np.argmin(vals))
        r = minimize_scalar(lambda t: expr(t, s, th, y, H), bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                            method='bounded', options={'xatol': 1e-12})
        if r.fun < best[0]: best = (r.fun, (math.degrees(r.x), s))
    return best

def aligned(th, y, H):
    out = []
    for d in (th, th + math.pi / 2):
        d = d % math.pi
        if 1e-9 < d < math.pi / 2 - 1e-9: out.append(expr(d, -1, th, y, H))
        if math.pi / 2 + 1e-9 < d < math.pi - 1e-9: out.append(expr(math.pi - d, 1, th, y, H))
    return min(out) if out else math.inf

if __name__ == '__main__':
    for H in [1.9, 1.95]:
        for thd in [10, 30, 45, 60, 80]:
            th = math.radians(thd)
            k = (math.cos(th) + math.sin(th)) / 2
            for y in [k, H / 2, H - k]:
                v, arg = X(th, y, H)
                print(f'H={H} th={thd} y={y:.3f}: X={v:.6f} at t={arg[0]:.3f} s={arg[1]:+d}   aligned={aligned(th, y, H):.6f}', flush=True)
