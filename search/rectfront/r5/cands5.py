#!/usr/bin/env python3
"""R_5 strip (2026-10-21): the middle square's width W(theta) is concave on each piece between kinks, so its
min over theta is at c aligned with line 1, aligned with line 2, or upright.  Check against a fine theta grid,
and report the minimum of B1 + W + B2 per candidate."""
import math
import numpy as np
from split5 import Bb, Bd, target

def W(th, m1, m2, H):
    P, Q = math.cos(th), math.sin(th); k = (P + Q) / 2
    h = lambda m: abs(P + m * Q) + abs(-Q + m * P)
    y = k if m2 >= m1 else H - k
    return (m2 - m1) * y + (h(m1) + h(m2)) / 2

def line_tilt(t, s):
    """the c tilt in [0, 90) aligned with a line of slope parameter m = s cot t (direction angle d)."""
    d = t if s < 0 else math.pi - t           # direction angle of the line
    return d % (math.pi / 2)

ts = np.radians(np.linspace(1, 89, 45))
worst_gap = 0
for H in [1.6, 1.85, 1.9, 1.95, 1.99]:
    res = {'grid': math.inf, 'al1': math.inf, 'al2': math.inf, 'up': math.inf}
    arg = {}
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            for t1 in ts:
                m1 = s1 / math.tan(t1); b1 = Bb(t1, s1, H)
                for t2 in ts:
                    m2 = s2 / math.tan(t2); b2 = Bd(t2, s2, H)
                    grid = min(W(th, m1, m2, H) for th in np.linspace(0, math.pi / 2, 721))
                    cands = {'al1': W(line_tilt(t1, s1), m1, m2, H), 'al2': W(line_tilt(t2, s2), m1, m2, H),
                             'up': W(0.0, m1, m2, H)}
                    worst_gap = min(worst_gap, min(cands.values()) - grid)
                    for kk, v in [('grid', grid)] + list(cands.items()):
                        tot = b1 + v + b2
                        if tot < res[kk]: res[kk] = tot; arg[kk] = (s1, s2, round(math.degrees(t1)), round(math.degrees(t2)))
    print(f'H={H} target={target(H):.6f} ' + '  '.join(f'{kk}={res[kk] - target(H):+.2e}@{arg[kk]}' for kk in res), flush=True)
print('min over candidates minus grid min (should be <= ~0):', worst_gap)
