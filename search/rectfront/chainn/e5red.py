#!/usr/bin/env python3
"""E5 reduced to per-square inequalities (exact kap, random samples), beta = 1 - H/2, a = 1 - lam:
 (A) c < a  =>  c >= beta (2 - |m| - |m'|)
 (B) c >= beta (|m| + |m'| - 2)
 (C') c >= beta (|m'| - |m|) - a
 (S1) c >= min(a + beta (|m'| - |m|), beta (2 - |m| - |m'|))"""
import math, sys, numpy as np
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
Hlo = float(sys.argv[1]) if len(sys.argv) > 1 else 1.87
worst = {k: (np.inf, None) for k in "A B C S1".split()}
for it in range(100):
    N = 2_000_000
    t = rng.uniform(0, math.pi / 2, N)
    # mix: uniform angles, plus samples near the edge slopes and +-1
    phi = rng.uniform(-1.5, 1.5, N); phi2 = rng.uniform(-1.5, 1.5, N)
    sel = rng.random(N) < 0.3; phi[sel] = (t[sel] - (rng.random(sel.sum()) < 0.5) * math.pi / 2) + rng.normal(0, 0.02, sel.sum())
    sel = rng.random(N) < 0.3; phi2[sel] = (t[sel] - (rng.random(sel.sum()) < 0.5) * math.pi / 2) + rng.normal(0, 0.02, sel.sum())
    m, mp = np.tan(phi), np.tan(phi2)
    H = rng.uniform(Hlo, 2, N)
    C, S = np.cos(t), np.sin(t); k = (C + S) / 2; ok = 2 * k <= H
    g = lambda m: np.abs(C + m * S) + np.abs(m * C - S)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2
    c = (g(m) + g(mp)) / 2 - np.abs(m - mp) * (H / 2 - k) - lam
    am, amp = np.abs(m), np.abs(mp)
    tests = {"A": np.where((c < a) & ok, c - b * (2 - am - amp), np.inf),
             "B": np.where(ok, c - b * (am + amp - 2), np.inf),
             "C": np.where(ok, c - (b * (amp - am) - a), np.inf),
             "S1": np.where(ok, c - np.minimum(a + b * (amp - am), b * (2 - am - amp)), np.inf)}
    for kk, v in tests.items():
        i = int(np.argmin(v))
        if v[i] < worst[kk][0]:
            worst[kk] = (v[i], (round(math.degrees(t[i]), 2), round(math.degrees(math.atan(m[i])), 2),
                                round(math.degrees(math.atan(mp[i])), 2), round(H[i], 4)))
for kk, v in worst.items(): print(kk, f"{v[0]:+.2e}", "(t, phi, phi', H) =", v[1])
