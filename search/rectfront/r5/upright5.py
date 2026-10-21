#!/usr/bin/env python3
"""Final5, middle square upright (2026-10-21): with y fixed (1/2 or H - 1/2) the bound splits into
left(t1, eps1, br1) + right(t2, eps2, br2).  Check  min left + min right >= target  for each y, ignoring mu's sign."""
import math, itertools
import numpy as np
from cases5 import bound

R2 = math.sqrt(2)
def left(t, e, br, H, y):
    C, S = math.cos(t), math.sin(t); v = C if e else -C
    return (bound(C, S, H, br) + (C * H if e else 0)) / S + (C + S) / (2 * S) - v / S * y
def right(t, e, br, H, y):
    C, S = math.cos(t), math.sin(t); v = C if e else -C
    return (bound(C, S, H, br) + (0 if e else C * H)) / S + (C + S) / (2 * S) + v / S * y

ts = np.radians(np.concatenate([np.linspace(0.02, 89.98, 2000), [45.0]]))
worst = math.inf
for H in np.concatenate([np.linspace(1.0, 1.9999, 81), [(6 * R2 - 1) / 4]]):
    tgt = min(5.0, 4 + 6 * R2 - 4 * H)
    for y in (0.5, H - 0.5):
        lm = min(left(t, e, br, H, y) for t in ts for e in (False, True) for br in 'UP')
        rm = min(right(t, e, br, H, y) for t in ts for e in (False, True) for br in 'UP')
        if lm + rm - tgt < worst: worst = lm + rm - tgt; arg = (round(H, 4), y, lm, rm)
print('upright candidate, decoupled: worst slack', worst, 'at', arg)
