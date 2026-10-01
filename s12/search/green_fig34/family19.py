"""Friedman's n = 19 set of Figure 34 (DS7), reconstructed (see GREEN_FIG34.md sec. 4).

Side t = 6*sqrt(2) - 4 = 2 + 3d with d = 2*sqrt(2) - 2.  Rows y = 1 + k d (k = 0..3).
  rows 0 and 3:  x = 1, 1+d, 1+2d, 1+3d                       (4 points each)
  row 1:         x = 1, 1+u, 1+u+d, 1+u+2d, 1+3d              (5 points)
  row 2:         180-degree image of row 1: x = 1, 1+d-u, 1+2d-u, 1+3d-u, 1+3d
18 points in all.  Friedman's figure has u ~ 0.554 (21 px of 37.9 px/unit).
"""
import math

def friedman19(t, u, unit=1.0):
    """t: side; d = (t - 2*unit)/3; u: offset (absolute units).  unit = wall distance."""
    d = (t - 2 * unit) / 3
    o = unit
    R0 = [o, o + d, o + 2 * d, o + 3 * d]
    R1 = [o, o + u, o + u + d, o + u + 2 * d, o + 3 * d]
    R2 = [o, o + d - u, o + 2 * d - u, o + 3 * d - u, o + 3 * d]
    P = []
    for k, row in enumerate((R0, R1, R2, R0)):
        P += [(x, o + k * d) for x in row]
    return P
