#!/usr/bin/env python3
"""Per-square conditions split by slope sign (2026-10-24), tangent potentials (lincheck.exact_pots), 1.849 <= H < 1.87.
 same sign: does kap_same / kap_B / kap_S / kap_wall (LB) suffice for S01, S10 and the cheap condition?
 opposite sign: true kernel; cheap condition WITHOUT the cheap restriction (c >= 2g - u_P(m) - u_P(m') everywhere)?"""
import sys, math, numpy as np
import lincheck as lc
from lbcheck import LB
from kwalk import Kmat

ms = lc.ms; phis = lc.phis; P1, P2 = np.meshgrid(ms, ms, indexing='ij')
same = P1 * P2 >= 0; opp = P1 * P2 < 0
for H in [1.849, 1.855, 1.86, 1.865, 1.8699]:
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2; g = a + b
    u0, u1 = lc.exact_pots(H)
    cl = LB(H, ms) - lam
    K, _, _ = Kmat(H, ms, 7201); ck = K - lam
    out = []
    for nm, c, mask in (('same/LB', cl, same), ('opp/K', ck, opp)):
        S01 = np.where(mask, u0[:, None] + c - u1[None, :], np.inf).min()
        S10 = np.where(mask, u1[:, None] + c - u0[None, :], np.inf).min()
        cheap = c < a - 1e-12
        Cr = min(np.where(mask & cheap, u[:, None] + c + u[None, :] - 2 * g, np.inf).min() for u in (u0, u1))
        Cu = min(np.where(mask, u[:, None] + c + u[None, :] - 2 * g, np.inf).min() for u in (u0, u1))
        out.append(f"{nm}: S01 {S01:+.1e} S10 {S10:+.1e} cheap {Cr:+.1e} cheap-unrestricted {Cu:+.1e}")
    print(f"H={H}  " + " | ".join(out), flush=True)
