#!/usr/bin/env python3
"""Closed-form parity potentials for odd n, H_b <= H < 2 (2026-10-24).
 E(m)  = a + beta |cot phi| + (2(C + S) - H)/(2 S C) - lam   (m = tan phi, C = |cos phi|, S = |sin phi|):
         one upright to the perpendicular line, then a square aligned with both lines.
 u0(m) = min(2a + beta|m|, E(m))      (even lines)
 u1(m) = min(a + beta|m|, a + E(m))   (odd lines)
Conditions (c = K - lam, K the kernel min over tilts):
 (W) c(0, m) >= u1(m);  (S10) u0(m') <= u1(m) + c;  (S01) u1(m') <= u0(m) + c;
 (C) c < a  =>  u_P(m) + c + u_P(m') >= 2g, P = 0, 1."""
import math, sys
import numpy as np
from kwalk import Kmat

R2 = math.sqrt(2)

def pots(H, ms):
    lam = 2 * R2 - H; a = 1 - lam; beta = 1 - H / 2
    phi = np.arctan(np.abs(ms)); C, S = np.cos(phi), np.sin(phi); am = np.abs(ms)
    with np.errstate(divide='ignore', invalid='ignore'):
        E = a + beta * (C / S) + (2 * (C + S) - H) / (2 * S * C) - lam
    E = np.where((S * C > 1e-12) & (C + S <= H), E, np.inf)
    return np.minimum(2 * a + beta * am, E), np.minimum(a + beta * am, a + E)

def check(H, phis, nt=3601):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, _, _ = Kmat(H, ms, nt)
    lam = 2 * R2 - H; a = 1 - lam; beta = 1 - H / 2; g = a + beta
    c = K - lam; z = int(np.argmin(np.abs(ms)))
    u0, u1 = pots(H, ms)
    r = {'W': (c[z] - u1).min(),
         'S10': (u1[:, None] + c - u0[None, :]).min(),
         'S01': (u0[:, None] + c - u1[None, :]).min()}
    cheap = c < a - 1e-12
    for P, u in ((0, u0), (1, u1)):
        r[f'C{P}'] = np.where(cheap, u[:, None] + c + u[None, :], np.inf).min() - 2 * g
    return r

if __name__ == '__main__':
    phis = np.arange(-85, 85 + 1e-9, float(sys.argv[1]) if len(sys.argv) > 1 else 0.25)
    for H in [1.8457, 1.846, 1.85, 1.86, 1.8685, 1.87, 1.9, 1.95, 1.99]:
        r = check(H, phis)
        print(f"H={H:.4f} " + "  ".join(f"{k} {v:+.1e}" for k, v in r.items())
              + ("" if min(r.values()) > -1e-12 else "   FAIL"), flush=True)
