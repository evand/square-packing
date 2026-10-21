#!/usr/bin/env python3
"""Validity of Friedman's R_4 configuration at tilt t, height H, width L = Lt(t, H) (2026-10-21).
a = [0,1]^2, d = [L-1,L] x [H-1,H], b at tilt t on the ceiling touching a's corner (1,1) with its lower-left
edge, c = b rotated by 180 deg about the centre.  Checks containment and pairwise interior-disjointness (SAT)."""
import math
import numpy as np
from frontier4 import Lt

def square(c, th):
    C, S = math.cos(th), math.sin(th); e, f = np.array([C, S]) / 2, np.array([-S, C]) / 2
    return np.array([c + e + f, c - e + f, c - e - f, c + e - f]), [np.array([C, S]), np.array([-S, C])]

def sep_margin(A, B):
    """max over the axes of both squares of the separation gap (>= 0 means interiors disjoint)."""
    best = -math.inf
    for n in A[1] + B[1]:
        pa, pb = A[0] @ n, B[0] @ n
        best = max(best, pb.min() - pa.max(), pa.min() - pb.max())
    return best

def config(t, H):
    C, S = math.cos(t), math.sin(t); k = (C + S) / 2
    L = Lt(t, H)
    yb = H - k; xb = 1 + (0.5 - (H - k - 1) * S) / C
    a = square(np.array([0.5, 0.5]), 0.0); d = square(np.array([L - 0.5, H - 0.5]), 0.0)
    b = square(np.array([xb, yb]), t); c = square(np.array([L - xb, H - yb]), t)
    return L, [a, b, c, d]

def check(t, H, tol=1e-12):
    L, sqs = config(t, H)
    inside = all(v[0] >= -tol and v[0] <= L + tol and v[1] >= -tol and v[1] <= H + tol for s in sqs for v in s[0])
    m = min(sep_margin(sqs[i], sqs[j]) for i in range(4) for j in range(i + 1, 4))
    return L, inside, m

if __name__ == '__main__':
    bad = []
    for H in np.linspace(1.0, 1.9999, 200):
        for td in np.linspace(0.5, 89.5, 357):
            L, ins, m = check(math.radians(td), H)
            if L < 4 and (not ins or m < -1e-9):
                bad.append((H, td, L, ins, m))
    print('configs with Lt < 4 that are invalid:', len(bad))
    for b in bad[:10]: print(b)
    if bad:
        print('t range of bad:', min(b[1] for b in bad), max(b[1] for b in bad), ' H range:', min(b[0] for b in bad), max(b[0] for b in bad))
    print('symmetry Lt(t)=Lt(90-t)?', [round(Lt(math.radians(x), 1.93) - Lt(math.radians(90 - x), 1.93), 9) for x in (30, 40, 44)])
