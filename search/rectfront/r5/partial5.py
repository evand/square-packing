#!/usr/bin/env python3
"""R_5 strip: partial decoupling of the middle square (2026-10-21).
V1: c's tilt shared, heights separate per half.  V2: heights shared, tilts separate.  Grid in (t1, t2)."""
import math
import numpy as np
from split5 import Bb, Bd, target

THS = np.linspace(0, math.pi / 2, 361)

def hsup(th, m):
    P, Q = math.cos(th), math.sin(th)
    return abs(P + m * Q) + abs(-Q + m * P)

def V1(m1, m2, H):
    best = math.inf
    for th in THS:
        P, Q = math.cos(th), math.sin(th); k = (P + Q) / 2
        if 2 * k > H: continue
        a = hsup(th, m1) / 2 + min(-m1 * k, -m1 * (H - k))
        b = hsup(th, m2) / 2 + min(m2 * k, m2 * (H - k))
        best = min(best, a + b)
    return best

def V2(m1, m2, H):
    best = math.inf
    for kk in np.linspace(0.5, math.sqrt(2) / 2, 60):       # shared height via shared k-range: y in {k, H-k}
        pass
    # heights shared: y free in [k_min, H - k_min]; each half picks its tilt with k <= min(y, H - y)
    for y in np.linspace(0.5, H - 0.5, 121):
        a = b = math.inf
        for th in THS:
            P, Q = math.cos(th), math.sin(th); k = (P + Q) / 2
            if k > y + 1e-12 or k > H - y + 1e-12: continue
            a = min(a, hsup(th, m1) / 2 - m1 * y)
            b = min(b, hsup(th, m2) / 2 + m2 * y)
        best = min(best, a + b)
    return best

ts = np.radians(np.linspace(1, 89, 45))
for H in [1.6, 1.8, 1.87, 1.9, 1.95, 1.99]:
    for nm, V in (('V1 shared tilt', V1), ('V2 shared y', V2)):
        best = (math.inf, None)
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                for t1 in ts:
                    m1 = s1 / math.tan(t1); b1 = Bb(t1, s1, H)
                    for t2 in ts:
                        m2 = s2 / math.tan(t2)
                        g = b1 + V(m1, m2, H) + Bd(t2, s2, H)
                        if g < best[0]: best = (g, (s1, s2, round(math.degrees(t1)), round(math.degrees(t2))))
        print(f'H={H} {nm}: min {best[0]:.6f} target {target(H):.6f} gap {best[0] - target(H):+.2e} at {best[1]}', flush=True)
