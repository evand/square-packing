#!/usr/bin/env python3
"""Final5, middle square aligned with line 1 (2026-10-21).  Relaxations: h2(a1) >= 1, then the right half
depends on c only through y; its min over line 2, rho(y), is concave, so rho >= chord of rho(1/2), rho(H - 1/2).
Check  min over t1 [left(t1, y*) + chord(y*)] >= target, y* = k1 (mu >= 0) or H - k1 (mu < 0), both tried."""
import math
import numpy as np
from cases5 import bound

R2 = math.sqrt(2)
TS = np.radians(np.concatenate([np.linspace(0.02, 89.98, 1500), [45.0]]))

def right(t, e, br, H, y):
    C, S = math.cos(t), math.sin(t); v = C if e else -C
    return (bound(C, S, H, br) + (0 if e else C * H)) / S + 1 / (2 * S) + v / S * y

def left(t, e, br, H, y):
    C, S = math.cos(t), math.sin(t); v = C if e else -C
    return (bound(C, S, H, br) + (C * H if e else 0)) / S + 1 / (2 * S) - v / S * y

worst = math.inf
for H in np.concatenate([np.linspace(1.0, 1.9999, 61), [(6 * R2 - 1) / 4]]):
    tgt = min(5.0, 4 + 6 * R2 - 4 * H)
    r_lo = min(right(t, e, br, H, 0.5) for t in TS for e in (False, True) for br in 'UP')
    r_hi = min(right(t, e, br, H, H - 0.5) for t in TS for e in (False, True) for br in 'UP')
    chord = lambda y: r_lo + (r_hi - r_lo) * (y - 0.5) / (H - 1) if H > 1 else r_lo
    for e in (False, True):
        for br in 'UP':
            for t in TS:
                C, S = math.cos(t), math.sin(t); k = (C + S) / 2
                if 2 * k > H: continue
                for y in (k, H - k):
                    val = left(t, e, br, H, y) + chord(y) - tgt
                    if val < worst: worst = val; arg = (round(H, 4), e, br, round(math.degrees(t), 2), round(y, 3))
print('aligned candidate, relaxed: worst slack', worst, 'at', arg)
