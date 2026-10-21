#!/usr/bin/env python3
"""R_5 strip: fully decoupled halves (2026-10-21).
left(t1)  = B_b(t1) + min_{c tilt, y} [h(m1)/2 - m1 y]
right(t2) = B_d(t2) + min_{c tilt, y} [h(m2)/2 + m2 y]
Need min left + min right >= target; reports each half's min and the target split."""
import math
import numpy as np
from split5 import Bb, Bd, target

def halfw(m, sgn_y, H):
    best = math.inf
    for th in np.linspace(0, math.pi / 2, 2001):
        P, Q = math.cos(th), math.sin(th); k = (P + Q) / 2
        if 2 * k > H: continue
        h = abs(P + m * Q) + abs(-Q + m * P)
        best = min(best, h / 2 + min(sgn_y * m * y for y in (k, H - k)))
    return best

ts = np.radians(np.concatenate([np.linspace(0.2, 89.8, 241), [45.0]]))
for H in [1.3, 1.6, 1.8, 1.85, 1.87, 1.8713, 1.88, 1.9, 1.95, 1.99]:
    L = R = (math.inf, None)
    for t in ts:
        for s in (-1, 1):
            m = s * math.cos(t) / math.sin(t)
            l = Bb(t, s, H) + halfw(m, -1, H)
            r = Bd(t, s, H) + halfw(m, +1, H)
            if l < L[0]: L = (l, (round(math.degrees(t), 2), s))
            if r < R[0]: R = (r, (round(math.degrees(t), 2), s))
    print(f'H={H}: left {L[0]:.6f} at {L[1]}, right {R[0]:.6f} at {R[1]}, sum {L[0] + R[0]:.6f}, '
          f'target {target(H):.6f}, gap {L[0] + R[0] - target(H):+.2e}', flush=True)
