#!/usr/bin/env python3
"""E5 closed-form certificate check (grid; K from kwalk, fine t grid).
u1(m) = 1 + |m|(1 - H/2) - lam;  u0(m) = min(2 - 2 lam + |m|(1 - H/2), aligned(phi)), m = tan phi,
aligned(phi) = 1 + |cot phi|(1 - H/2) + (2(C+S) - H)/(2 S C) - 2 lam  (C = |cos phi|, S = |sin phi|).
Conditions (c = K - lam):
 (W)  c(0, m) >= u1(m)
 (S1) u0(m') <= u1(m) + c(m, m')       (S2) u1(m') <= u0(m) + c(m, m')
 (C)  c(m, m') < a  =>  u_j(m) + c(m, m') + u_j(m') >= 2g  (j = 0, 1);  wall: c(0, m) < a => c(0, m) + u0(m) >= 2g."""
import math, sys, numpy as np
from kwalk import Kmat

def pots(H, ms):
    lam = 2 * math.sqrt(2) - H
    phi = np.arctan(np.abs(ms)); C, S = np.cos(phi), np.sin(phi)
    u1 = 1 + np.abs(ms) * (1 - H / 2) - lam
    with np.errstate(divide='ignore', invalid='ignore'):
        al = 1 + (C / S) * (1 - H / 2) + (2 * (C + S) - H) / (2 * S * C) - 2 * lam
    al = np.where(S * C > 1e-12, al, np.inf)
    u0 = np.minimum(2 - 2 * lam + np.abs(ms) * (1 - H / 2), al)
    return u1, u0

def check(H, phis):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, _, _ = Kmat(H, ms)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    c = K - lam; z = int(np.argmin(np.abs(ms)))
    u1, u0 = pots(H, ms)
    W = (c[z] - u1).min()
    S1 = (u1[:, None] + c - u0[None, :]).min()
    S2 = (u0[:, None] + c - u1[None, :]).min()
    cheap = c < a - 1e-12
    Cm = min(np.where(cheap, u[:, None] + c + u[None, :], np.inf).min() for u in (u0, u1)) - 2 * g
    Cw = np.where(c[z] < a - 1e-12, c[z] + u0, np.inf).min() - 2 * g
    return W, S1, S2, Cm, Cw

if __name__ == '__main__':
    phis = np.linspace(-85, 85, 681)
    for H in [float(x) for x in sys.argv[1:]] or [1.83, 1.85, 1.86, 1.87, 1.88, 1.9, 1.93, 1.96, 1.99]:
        W, S1, S2, Cm, Cw = check(H, phis)
        print(f"H={H:.3f}: W {W:+.1e}  S1 {S1:+.1e}  S2 {S2:+.1e}  C {Cm:+.1e}  Cwall {Cw:+.1e}", flush=True)
