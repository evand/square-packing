#!/usr/bin/env python3
"""E5 probe: cheap-step condition with parity-free prefix distance d = min(dS0, dS1), and slack landscape."""
import math, sys, numpy as np
from excess import dists
from kwalk import Kmat
phis = np.linspace(-84, 84, 337)
ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
for H in [float(x) for x in sys.argv[1:]] or [1.85, 1.87, 1.9, 1.93, 1.96, 1.99]:
    K, _, _ = Kmat(H, ms); c, dS, dE, z = dists(H, ms, K)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    d = np.minimum(dS[0], dS[1]); d1 = dS[1]
    cheap = c < a - 1e-12
    for name, u in (("parity-free min(dS0,dS1)", d), ("one-step c(0,m)", d1), ("const g... (u=min(g, c(0,m)))", np.minimum(g, d1))):
        T = np.where(cheap, u[:, None] + c + u[None, :], np.inf)
        i, k = np.unravel_index(np.argmin(T), T.shape)
        # validity as a prefix lower bound: u(m') <= u(m) + c(m,m') and u <= c(0,.)
        viol = (u[None, :] - u[:, None] - c).max(); wall = (u - c[z]).max()
        print(f"H={H:.2f} {name:32s}: cheap min - 2g = {T[i,k]-2*g:+.2e} at {phis[i]:+.1f}->{phis[k]:+.1f};"
              f" step-ineq viol {viol:+.1e}, wall viol {wall:+.1e}")
