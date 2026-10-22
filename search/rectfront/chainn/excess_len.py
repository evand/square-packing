#!/usr/bin/env python3
"""Length-aware check near the threshold: for odd n, min excess over odd walks of length n that contain a cheap step,
vs 2g, only for n with n a >= 2g (else the grid bound follows by monotonicity in H).  Also: how many steps the
shortest prefixes dS use (structure of the distance functions)."""
import math, sys
import numpy as np
from kwalk import Kmat
phis = np.linspace(-84, 84, 337)
ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
z = int(np.argmin(np.abs(ms)))
for H in [1.83, 1.835, 1.84, 1.845, 1.85, 1.9]:
    K, _, _ = Kmat(H, ms)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    c = K - lam; cheap = c < a - 1e-12
    N = math.ceil(2 * g / a)
    # forward DP over length with flag "seen a cheap step": A = no cheap yet, B = cheap seen
    A = c[z].copy(); B = np.where(cheap[z], c[z], np.inf); A = np.where(cheap[z], np.inf, A)
    worst = []
    for n in range(2, max(N + 40, 60)):
        endB = min((B + c[:, z]).min(), (A + np.where(cheap[:, z], c[:, z], np.inf)).min())
        if n % 2 == 1 and n >= N:
            worst.append((n, endB - 2 * g))
        nA = (A[:, None] + np.where(cheap, np.inf, c)).min(0)
        nB = np.minimum((B[:, None] + c).min(0), (A[:, None] + np.where(cheap, c, np.inf)).min(0))
        A, B = nA, nB
    w = min(worst, key=lambda t: t[1])
    print(f"H={H}: N={N}; over odd n in [N, {worst[-1][0]}]: min (excess - 2g) = {w[1]:+.3e} at n={w[0]};"
          f" at n=N: {worst[0][1]:+.3e}", flush=True)
