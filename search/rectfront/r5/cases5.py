#!/usr/bin/env python3
"""Final5 by case (2026-10-21): eps1, eps2, middle candidate x, one-square branch for b and for d.
Worst slack of the conclusion over (t1, t2, H) per case, and where."""
import math, itertools
import numpy as np
from final5_check import cr, dt, midW

R2 = math.sqrt(2)
def bound(C, S, H, br):
    if br == 'U': return 2 * S - C * (H - 1)
    return (C * (S - C) + 1 + C + S - H) / C if C > 1e-15 else math.inf

def least_L(e1, e2, t1, t2, H, xi, br1, br2):
    C1, S1, C2, S2 = math.cos(t1), math.sin(t1), math.cos(t2), math.sin(t2)
    a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
    v1 = C1 if e1 else -C1; v2 = C2 if e2 else -C2
    s1 = bound(C1, S1, H, br1) + (C1 * H if e1 else 0)
    x = [(1.0, 0.0), a1, a2, (0.0, 1.0)][xi]
    mu = v2 / S2 - v1 / S1
    y = (x[0] + x[1]) / 2 if mu >= 0 else H - (x[0] + x[1]) / 2
    m = midW(mu, 1 / (2 * S1), 1 / (2 * S2), a1, a2, x, y)
    return s1 / S1 + m + (bound(C2, S2, H, br2) + (0 if e2 else C2 * H)) / S2

ts = np.radians(np.concatenate([np.linspace(0.05, 89.95, 181), [45.0]]))
Hs = np.concatenate([np.linspace(1.0, 1.999, 41), [1.8713]])
rows = []
if __name__ == "__main__":
  for e1, e2, xi, b1, b2 in itertools.product((False, True), (False, True), (0, 1, 2), 'UP', 'UP'):
      worst = (math.inf, None)
      for H in Hs:
          tgt = min(5.0, 4 + 6 * R2 - 4 * H)
          for t1 in ts:
              for t2 in ts:
                  L = least_L(e1, e2, t1, t2, H, xi, b1, b2)
                  if L - tgt < worst[0]: worst = (L - tgt, (round(H, 3), round(math.degrees(t1), 2), round(math.degrees(t2), 2)))
      rows.append((worst[0], e1, e2, ['up', 'a1', 'a2'][xi], b1, b2, worst[1]))
  for r in sorted(rows, key=lambda r: r[0]):
      print(f'slack {r[0]:+.4f}  eps=({int(r[1])},{int(r[2])}) x={r[3]} branches={r[4]}{r[5]}  at (H,t1,t2)={r[6]}')
