#!/usr/bin/env python3
"""Shortest-path form of the odd-n strip theorem (2026-10-22).
c(m, m') = K(m, m') - lam >= 0.  Odd n: L >= min(n, F_n)  <=  sum c >= min(n a, 2g).
Sufficient: every odd walk with a cheap step (c < a) has sum c >= 2g.  With nonnegative c this is
   min over cheap steps (m at parity j -> m' at 1-j) of  dS_j(m) + c(m, m') + dE_{1-j}(m')  >= 2g,
dS_j(m) = least excess of a prefix from the left wall (index 0) to a line m at index parity j,
dE_j(m) = least excess of a suffix from a line m at index parity j to the right wall (index n odd)."""
import math, sys
import numpy as np
from kwalk import Kmat

def dists(H, ms, K):
    lam = 2 * math.sqrt(2) - H
    c = K - lam
    z = int(np.argmin(np.abs(ms)))
    INF = 1e9
    # prefix: index 1 line has parity 1, first step wall -> m: c(0, m)
    dS = np.full((2, len(ms)), INF); dS[1] = c[z]
    for _ in range(400):
        new = dS.copy()
        new[0] = np.minimum(new[0], (dS[1][:, None] + c).min(0))
        new[1] = np.minimum(new[1], (dS[0][:, None] + c).min(0))
        if np.allclose(new, dS, atol=0, rtol=0): break
        dS = new
    # suffix to the right wall at odd index: last line at index n-1 even -> parity 0, step c(m, 0)
    dE = np.full((2, len(ms)), INF); dE[0] = c[:, z]
    for _ in range(400):
        new = dE.copy()
        new[1] = np.minimum(new[1], (c + dE[0][None, :]).min(1))
        new[0] = np.minimum(new[0], (c + dE[1][None, :]).min(1))
        if np.allclose(new, dE, atol=0, rtol=0): break
        dE = new
    return c, dS, dE, z

def check(H, phis):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, _, _ = Kmat(H, ms)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    c, dS, dE, z = dists(H, ms, K)
    worst = np.inf; arg = None
    for j in (0, 1):
        T = dS[j][:, None] + c + dE[1 - j][None, :]
        T = np.where(c < a - 1e-12, T, np.inf)
        i, k = np.unravel_index(np.argmin(T), T.shape)
        if T[i, k] < worst: worst, arg = T[i, k], (j, phis[i], phis[k])
    # also the whole-walk minimum excess (any odd walk, cheap or not), for reference
    allmin = min((dS[j][:, None] + c + dE[1 - j][None, :]).min() for j in (0, 1))
    return worst - 2 * g, arg, g, a, dS, dE, ms

if __name__ == '__main__':
    phis = np.linspace(-84, 84, 337)
    for H in [float(x) for x in sys.argv[1:]] or [1.83, 1.835, 1.84, 1.85, 1.87, 1.9, 1.93, 1.96, 1.99]:
        gap, arg, g, a, dS, dE, ms = check(H, phis)
        print(f"H={H:.3f}: min over cheap steps of (dS + c + dE) - 2g = {gap:+.3e} at parity {arg[0]} "
              f"{arg[1]:+.1f} -> {arg[2]:+.1f}   (2g = {2*g:.4f}, a = {a:.4f})", flush=True)
