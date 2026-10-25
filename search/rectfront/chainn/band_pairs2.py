#!/usr/bin/env python3
"""Breakpoint pairs for the final band design (2026-10-25): H in [H*(9), 1.87], H*(9) = (14 sqrt2 - 5)/8;
T through (8/7, Lz(8/7)) with slope H/16 - 1/10; exact region pieces on both sides:
  u0: A0 on [0, x2], Lz on [x2, 8/7], T on [8/7, inf),  x2 = 1 - a/(2b)
  u1: O on [0, y1], a + T on [y1, inf),                  y1 = T0/(b - T1) (>= 8/7 iff H >= H*(9))."""
import math, sys
import numpy as np
from band_dec import consts
from band_pairs import cands

R2 = math.sqrt(2); Q = 8 / 7; H9 = (14 * R2 - 5) / 8

def pieces(H):
    lam, a, b, g, _, T1 = consts(H)
    Lz = (a + 2 * b, -b); T0 = Lz[0] + Lz[1] * Q - T1 * Q; T = (T0, T1)
    x2 = 1 - a / (2 * b); y1 = T0 / (b - T1)
    A0 = (2 * a, b); O = (a, b)
    u0 = [('A0', 0.0, x2, A0), ('Lz', x2, Q, Lz), ('T', Q, math.inf, T)]
    u1 = [('O', 0.0, y1, O), ('aT', y1, math.inf, (a + T0, T1))]
    return dict(lam=lam, a=a, b=b, g=g, u0=u0, u1=u1, y1=y1, x2=x2)

def conds(H, C, S, P):
    lam, g = P['lam'], P['g']; p = C + S
    fx = (p / 2, S - H / 2); kx = (C, S); gy = (p / 2, C - H / 2); ky = (S, C)
    res = {}
    def pair(name, X, Y, const):
        for lx, vx in X:
            for ly, vy in Y:
                res[(name, lx, ly)] = vx + vy + const
    z = [('-', 0.0, math.inf, (0.0, 0.0))]
    pair('Sa', cands(fx, kx, +1, P['u0']), cands(gy, ky, -1, P['u1']), -lam)
    pair('Sb', cands(fx, kx, +1, P['u1']), cands(gy, ky, -1, P['u0']), -lam)
    pair('N', cands(fx, kx, +1, z, hi=0.5), cands(gy, ky, +1, z), -1.0)
    pair('C0', cands(fx, kx, +1, P['u0'], lo=0.5), cands(gy, ky, +1, P['u0'], lo=0.5), -lam - 2 * g)
    pair('C1', cands(fx, kx, +1, P['u1'], lo=0.5), cands(gy, ky, +1, P['u1'], lo=0.5), -lam - 2 * g)
    return res

if __name__ == '__main__':
    worst = {}
    for H in np.linspace(H9, 1.87, 9):
        P = pieces(H)
        for t in np.linspace(0, math.pi / 2, 3601):
            C, S = math.cos(t), math.sin(t)
            for k, v in conds(H, C, S, P).items():
                if k not in worst or v < worst[k][0]: worst[k] = (v, H, math.degrees(t))
    by = {}
    for (nm, lx, ly), (v, H, td) in worst.items(): by.setdefault(nm, []).append((v, lx, ly, H, td))
    tot = 0
    for nm, rows in by.items():
        rows.sort(); tot += len(rows)
        print(f"{nm}: {len(rows)} pairs, min {rows[0][0]:+.2e}; below 1e-3: "
              + "; ".join(f"{v:+.1e} {lx}/{ly} (H {H:.4f}, t {td:.2f})" for v, lx, ly, H, td in rows if v < 1e-3))
    print("total", tot)
