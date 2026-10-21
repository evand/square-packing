#!/usr/bin/env python3
"""R_4 tilted frontier (2026-10-21).  Point-symmetric family: uprights in opposite corners, b on the ceiling
and c on the floor, both at tilt t; contacts a-corner/b-edge, b-c edge (normal -e_perp), c-edge/d-corner.
   L(t, H) = 2 + (1 - 2(H - k - 1) S)/C + (1 - (H - 2k) C)/S,   C = cos t, S = sin t, k = (C + S)/2.
Frontier: F(H) = min_t L(t, H)."""
import math
from scipy.optimize import minimize_scalar

def Lt(t, H):
    C, S = math.cos(t), math.sin(t); k = (C + S) / 2
    return 2 + (1 - 2 * (H - k - 1) * S) / C + (1 - (H - 2 * k) * C) / S

def F(H):
    r = minimize_scalar(lambda t: Lt(t, H), bounds=(0.05, math.pi / 2 - 0.05), method='bounded',
                        options={'xatol': 1e-13})
    return r.fun, math.degrees(r.x)

if __name__ == '__main__':
    print('check measured 3.912287 x 1.912287 @39.95:', Lt(math.radians(39.95), 1.912287))
    for H in [1.86, 1.88, 1.9, 1.912287, 1.93, 1.95, 1.97, 1.98, 1.99, 2.0]:
        print(H, *F(H))
