#!/usr/bin/env python3
"""E5 with piecewise-linear u0: replace the aligned route by tangent lines at |m| = 1 (and optionally more),
u0 = min(R_up, max_k tangent_k).  Checks the four conditions as in e5cert.py."""
import math, sys, numpy as np
from kwalk import Kmat
from e5cert import pots

def Ral(H, m):
    lam = 2 * math.sqrt(2) - H; phi = np.arctan(np.abs(m)); C, S = np.cos(phi), np.sin(phi)
    return 1 + (C / S) * (1 - H / 2) + (2 * (C + S) - H) / (2 * S * C) - 2 * lam

def check(H, phis, tang_pts):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, _, _ = Kmat(H, ms)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    c = K - lam; z = int(np.argmin(np.abs(ms)))
    u1, u0ex = pots(H, ms)
    Rup = 2 * a + np.abs(ms) * (1 - H / 2)
    am = np.abs(ms); h = 1e-6
    lines = []
    for x0 in tang_pts:
        f0 = Ral(H, x0); d = (Ral(H, x0 + h) - Ral(H, x0 - h)) / (2 * h)
        lines.append(f0 + d * (am - x0))
    u0 = np.minimum(Rup, np.max(lines, axis=0))
    S1 = (u1[:, None] + c - u0[None, :]).min()
    S2 = (u0[:, None] + c - u1[None, :]).min()
    cheap = c < a - 1e-12
    Cm = min(np.where(cheap, u[:, None] + c + u[None, :], np.inf).min() for u in (u0, u1)) - 2 * g
    return S1, S2, Cm, (u0 - u0ex).min()

if __name__ == '__main__':
    phis = np.linspace(-85, 85, 681)
    for pts in ([1.0], [1.0, 1.3], [0.8, 1.0, 1.3]):
        print("tangents at |m| =", pts)
        for H in [1.87, 1.88, 1.9, 1.93, 1.96, 1.99]:
            S1, S2, Cm, lo = check(H, phis, pts)
            print(f"  H={H:.2f}: S1 {S1:+.1e} S2 {S2:+.1e} C {Cm:+.1e}  (u0 - exact) min {lo:+.1e}", flush=True)
