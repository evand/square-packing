#!/usr/bin/env python3
"""Exact limits of the odd-n strip method (2026-10-24).
 H_C: least H with KapC (Lean, StripOdd) true.  Aligned-lines reduction (StripC docstring):
      f(t, H) = p + C^2 - S^2 - C^2 H - (4 sqrt2 - 1 - 2H) S C >= 0, linear in H; H_C = max_t root.  Grid check of the
      full KapC margin alongside.
 H_b: least H for which n-independent parity potentials work: d_P(m) = inf over prefixes ending at a line of index
      parity P of the excess sum (min-plus fixed point on 2 x M states); condition: every cheap step (c < a) between
      lines m, m' has d_P(m) + c + d_P(m') >= 2g, P = 0, 1 (prefix and suffix have the same parity for odd n).
      Then L >= min(n a, 2g) + n lam = min(n, F_n) for every odd n.  Grid values over-estimate d (fewer routes), so
      H_b from the grid is a slight under-estimate; refine the slope grid to see convergence.
 Thresholds H*(n) (F_n = n): odd n, H*(n) = (2(n-2) sqrt2 + 4 - n)/(n - 1)."""
import math, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from kwalk import Kmat

R2 = math.sqrt(2)

def HC_aligned():
    def root(t):
        C, S = math.cos(t), math.sin(t); p = C + S
        A = p + C * C - S * S - (4 * R2 - 1) * S * C; B = 2 * S * C - C * C   # f = A + H B
        return -A / B if B > 0 else -np.inf
    r = minimize_scalar(lambda t: -root(t), bounds=(math.atan(0.6), math.atan(10)), method='bounded',
                        options={'xatol': 1e-13})
    return -r.fun, math.tan(r.x)

def grids(H, phis):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, _, _ = Kmat(H, ms)
    lam = 2 * R2 - H; a = 1 - lam; beta = 1 - H / 2; g = a + beta
    return ms, K - lam, lam, a, beta, g

def kapc_margin(H, phis):
    ms, c, lam, a, beta, g = grids(H, phis); am = np.abs(ms)
    return (c + lam - (2 * lam - 1 + beta * (am[None, :] - am[:, None]))).min()

def parity_margin(H, phis):
    ms, c, lam, a, beta, g = grids(H, phis)
    z = int(np.argmin(np.abs(ms)))
    d = [np.full(len(ms), np.inf), c[z].copy()]          # d[P]: parity P = index mod 2; line 1 from the wall
    for _ in range(400):
        n0 = np.minimum(d[0], (d[1][:, None] + c).min(0))
        n1 = np.minimum(d[1], (d[0][:, None] + c).min(0))
        if np.array_equal(n0, d[0]) and np.array_equal(n1, d[1]): break
        d = [n0, n1]
    cheap = c < a - 1e-12
    return min(np.where(cheap, d[P][:, None] + c + d[P][None, :], np.inf).min() for P in (0, 1)) - 2 * g

def Hstar(n): return (2 * (n - 2) * R2 + 4 - n) / (n - 1)

if __name__ == '__main__':
    hc, m = HC_aligned()
    print(f"H_C (aligned reduction) = {hc:.10f} at tan t = {m:.6f}", flush=True)
    phis = np.arange(-85, 85 + 1e-9, 0.5)
    h1 = brentq(lambda H: kapc_margin(H, phis), 1.86, 1.875, xtol=1e-7)
    print(f"slope grid 0.5 deg: KapC margin zero at H = {h1:.7f}", flush=True)
    # parity margin is 0 where it holds (zigzag tight), < 0 below the limit: tabulate
    for step in (0.5, 0.25):
        phis = np.arange(-85, 85 + 1e-9, step)
        row = [(H, parity_margin(H, phis)) for H in (1.84, 1.845, 1.8475, 1.85, 1.8525, 1.855, 1.86, 1.87)]
        print(f"slope grid {step} deg, parity margin:", " ".join(f"{H}:{v:+.1e}" for H, v in row), flush=True)
    print("H*(n), odd n:", " ".join(f"{n}:{Hstar(n):.5f}" for n in range(3, 22, 2)), f" limit {2*R2-1:.5f}")
