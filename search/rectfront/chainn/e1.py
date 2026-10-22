#!/usr/bin/env python3
"""E1: per-square excess kappa(t; m, m') - lam >= 0, for 2 sqrt2 - 1 <= H <= 2, 2 k_t <= H.
kappa = (g_t(m) + g_t(m'))/2 - |m - m'| (H/2 - k_t),  lam = 2 sqrt2 - H."""
import math, numpy as np
rng = np.random.default_rng(1)
def kap(t, m, mp, H):
    C, S = np.cos(t), np.sin(t); k = (C + S) / 2
    g = lambda m: np.abs(C + m * S) + np.abs(m * C - S)
    return (g(m) + g(mp)) / 2 - np.abs(m - mp) * (H / 2 - k) - (2 * math.sqrt(2) - H), 2 * k <= H
best = (np.inf, None)
for it in range(200):
    N = 2_000_000
    t = rng.uniform(0, math.pi / 2, N)
    m = np.tan(rng.uniform(-1.55, 1.55, N)); mp = np.tan(rng.uniform(-1.55, 1.55, N))
    H = rng.uniform(2 * math.sqrt(2) - 1, 2, N)
    v, ok = kap(t, m, mp, H)
    v = np.where(ok, v, np.inf)
    i = int(np.argmin(v))
    if v[i] < best[0]: best = (v[i], (math.degrees(t[i]), math.degrees(math.atan(m[i])), math.degrees(math.atan(mp[i])), H[i]))
print("min excess", best)
