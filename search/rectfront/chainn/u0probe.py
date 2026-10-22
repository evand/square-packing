#!/usr/bin/env python3
"""E5: structure of the even-parity prefix potential u0(m') = min_{m,t} [u1(m) + kap(t; m, m') - lam],
u1(m) = 1 + |m|(1 - H/2) - lam (one-step wall cost, upright square).  Also check u1 = full dS1 and u0 = full dS0."""
import math, sys, numpy as np
from excess import dists
from kwalk import Kmat
phis = np.linspace(-84, 84, 337)
ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
ts = np.radians(np.linspace(0, 90, 1801))
for H in (1.85, 1.9, 1.96):
    lam = 2 * math.sqrt(2) - H
    u1 = 1 + np.abs(ms) * (1 - H / 2) - lam
    K, _, _ = Kmat(H, ms); c, dS, dE, z = dists(H, ms, K)
    print(f"H={H}: max |u1 - dS1| = {np.abs(u1 - dS[1]).max():.2e}")
    C, S = np.cos(ts), np.sin(ts); k = (C + S) / 2; ok = 2 * k <= H
    g = lambda m: np.abs(C + m * S) + np.abs(m * C - S)
    for ph2 in (0, 20, 35, 45, 52, 60, 75):
        mp = math.tan(math.radians(ph2))
        best = (np.inf,)
        for i, m in enumerate(ms):
            v = u1[i] + (g(m) + g(mp)) / 2 - abs(m - mp) * (H / 2 - k) - lam
            v = np.where(ok, v, np.inf); j = int(np.argmin(v))
            if v[j] < best[0]: best = (v[j], phis[i], math.degrees(ts[j]))
        jj = int(np.argmin(np.abs(phis - ph2)))
        print(f"   m'={ph2:3d}deg: u0 = {best[0]:.5f} (dS0 {dS[0][jj]:.5f}) via m={best[1]:+.1f}deg, t={best[2]:.1f}deg")
