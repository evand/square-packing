#!/usr/bin/env python3
"""D-case reduction (2026-10-21): tau_D = S - C + sigma+ + gamma (1 + P + Q - H)/P,  P = cos th, Q = sin th,
C = cos t, S = sin t, sigma = SP - CQ = sin(t - th), gamma = SQ + CP = cos(t - th).
 D1 (th <= t): tau_D - lam Psi - (1 - lam) Up >= 0, lam = CQ/(SP)  (H drops out)
 D2 (th >= t): tau_D - Psi >= 0 at H = 2 (coefficient of H is S sigma/(C P) <= 0)."""
import math
import numpy as np

def tauD(P, Q, C, S, H):
    sg = S * P - C * Q; gm = S * Q + C * P
    return S - C + max(sg, 0) + gm * (1 + P + Q - H) / P

def Psi(C, S, H): return S - C + (1 + C + S - H) / C
def Up(C, S, H): return 2 * S - C * (H - 1)

w1 = w2 = math.inf; a1 = a2 = None
for t in np.radians(np.linspace(0.2, 89.8, 449)):
    C, S = math.cos(t), math.sin(t)
    for th in np.radians(np.linspace(0, 89.8, 450)):
        P, Q = math.cos(th), math.sin(th)
        if th <= t:
            lam = C * Q / (S * P)
            v = [tauD(P, Q, C, S, H) - lam * Psi(C, S, H) - (1 - lam) * Up(C, S, H) for H in (1.0, 1.5, 2.0)]
            assert max(v) - min(v) < 1e-9, v
            if v[0] < w1: w1, a1 = v[0], (math.degrees(th), math.degrees(t))
        else:
            v = tauD(P, Q, C, S, 2.0) - Psi(C, S, 2.0)
            if v < w2: w2, a2 = v, (math.degrees(th), math.degrees(t))
print('D1 min', w1, 'at (th, t)', a1)
print('D2 min', w2, 'at (th, t)', a2)
