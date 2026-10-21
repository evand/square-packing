#!/usr/bin/env python3
"""Diagnose aligned5: exact rho(y) (min over line 2 at that y) with h2 >= 1, no chord."""
import math
import numpy as np
from aligned5 import left, right, R2
TS = np.radians(np.concatenate([np.linspace(0.05, 89.95, 600), [45.0]]))
worst = math.inf
for H in np.concatenate([np.linspace(1.0, 1.9999, 21), [(6 * R2 - 1) / 4]]):
    tgt = min(5.0, 4 + 6 * R2 - 4 * H)
    rho = lambda y: min(right(t, e, br, H, y) for t in TS for e in (False, True) for br in 'UP')
    for e in (False, True):
        for br in 'UP':
            for t in TS[::3]:
                C, S = math.cos(t), math.sin(t); k = (C + S) / 2
                if 2 * k > H: continue
                for y in (k, H - k):
                    val = left(t, e, br, H, y) + rho(y) - tgt
                    if val < worst: worst = val; arg = (round(H, 4), e, br, round(math.degrees(t), 2), round(y, 3))
print('aligned, h2>=1, exact rho: worst slack', worst, 'at', arg)
