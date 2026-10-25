#!/usr/bin/env python3
"""Do the conditions follow from the per-square lemmas already in Lean?  (2026-10-24)
Replace the kernel by LB(m, m') = max(kap_B, kap_S, kap_same [same sign], kap_wall [m or m' = 0]) and recheck
W, S10, S01, C0, C1 for the tangent potentials (lincheck.exact_pots).  Also report where the true kernel is needed."""
import sys, math, numpy as np
import lincheck as lc
from kwalk import Kmat

def LB(H, ms):
    lam = 2 * math.sqrt(2) - H; b = 1 - H / 2; Z = 2 + 2 * math.sqrt(2) - 2 * H
    x = np.abs(ms)[:, None]; y = np.abs(ms)[None, :]; m = ms[:, None]; mp = ms[None, :]
    B = 2 * (math.sqrt(2) - 1) + b * (x + y)
    S = np.minimum(1 + b * np.abs(x - y), Z - b * (x + y))
    same = np.where(m * mp >= 0, 1 + b * np.abs(x - y), -np.inf)
    wall = np.where((m == 0) | (mp == 0), 1 + b * np.maximum(x, y), -np.inf)
    return np.maximum.reduce([B, S, same, wall])

def conds(H, c, u0, u1, z, a, g):
    cheap = c < a - 1e-12
    r = {'W': (c[z] - u1).min(), 'S10': (u1[:, None] + c - u0[None, :]),
         'S01': (u0[:, None] + c - u1[None, :]),
         'C0': np.where(cheap, u0[:, None] + c + u0[None, :] - 2 * g, np.inf),
         'C1': np.where(cheap, u1[:, None] + c + u1[None, :] - 2 * g, np.inf)}
    return r

if __name__ == '__main__':
    ms = lc.ms; phis = lc.phis; z = int(np.argmin(np.abs(ms)))
    for H in [1.849, 1.86, 1.8699]:
        lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2; g = a + b
        u0, u1 = lc.exact_pots(H)
        cl = LB(H, ms) - lam
        K, _, _ = Kmat(H, ms, 3601); ck = K - lam
        print(f"H={H}: LB <= K violations: {(cl - ck).max():+.1e}")
        for nm, M in conds(H, cl, u0, u1, z, a, g).items():
            if nm == 'W': print(f"   W {M:+.2e}"); continue
            i, j = np.unravel_index(np.argmin(M), M.shape)
            print(f"   {nm}: min {M[i, j]:+.2e} at ({phis[i]:+.2f}, {phis[j]:+.2f})", flush=True)
