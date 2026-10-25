#!/usr/bin/env python3
"""Band theorem with fixed rational breakpoints (2026-10-25).

Change of potentials: T passes through (q, Lz(q)), q = 8/7, with the tangent slope T1 = H/16 - 1/10, so L = max(Lz, T)
switches exactly at q for every H.  Branch splits on the - side are free choices (each branch is an upper bound of the
potential everywhere); fixed rationals y1 (u1: O below, a + L above) and y2 (u0: A0 below, L above).  With all
breakpoints fixed, every leaf inequality is affine in H.

Evaluation by exact candidates (band_pairs.cands): min over a piece is at its ends or at the frame kink.
Reports, per H, the feasible windows for y1 and y2 (all sub-inequalities of Sa resp. Sb >= 0), and checks C0, C1, N."""
import math, sys
import numpy as np
from band_dec import consts
from band_pairs import cands

Q = 8 / 7

def lines(H):
    lam, a, b, g, T0_, T1 = consts(H)
    Lz = (a + 2 * b, -b)
    Lzq = Lz[0] + Lz[1] * Q
    T = (Lzq - T1 * Q, T1)
    return lam, a, b, g, Lz, T

def frames(H, n=1801):
    for t in np.linspace(0, math.pi / 2, n):
        C, S = math.cos(t), math.sin(t)
        if C + S <= H: yield C, S, t

def side(C, S, H, xside):
    p = C + S
    return ((p / 2, S - H / 2), (C, S)) if xside else ((p / 2, C - H / 2), (S, C))

def worst_sum(H, X, Y, const, n=1801):
    """X, Y: lists of (pieces, sign, lo, hi) groups; condition holds iff for every frame
    min over all X groups' candidates + min over all Y groups' candidates + const >= 0 -- but '+ side min' groups
    are separate requirements, so X is a list of alternatives that must EACH hold: we return min over x-groups,
    y-groups of (min_x + min_y)."""
    w = math.inf
    for C, S, t in frames(H, n):
        fx, kx = side(C, S, H, True); gy, ky = side(C, S, H, False)
        for pcsx, sx, lox, hix in X:
            vx = min(v for _, v in cands(fx, kx, sx, pcsx, lox, hix))
            for pcsy, sy, loy, hiy in Y:
                vy = min(v for _, v in cands(gy, ky, sy, pcsy, loy, hiy))
                w = min(w, vx + vy + const)
    return w

def conds(H, y1, y2, n=1801):
    lam, a, b, g, Lz, T = lines(H)
    A0 = (2 * a, b); O = (a, b); aLz = (a + Lz[0], Lz[1]); aT = (a + T[0], T[1])
    inf = math.inf
    # + side: u0 >= A0 everywhere, and u0 = min(A0, L) >= ... both branches: A0 globally; L as Lz on [0,Q], T on [Q,inf)
    u0plus = [([('A0', 0, inf, A0)], +1), ([('Lz', 0, Q, Lz), ('T', Q, inf, T)], +1)]
    u1plus = [([('O', 0, inf, O)], +1), ([('aLz', 0, Q, aLz), ('aT', Q, inf, aT)], +1)]
    # - side, branch per region
    def minus_u1(y1):
        pcs = [('O', 0, y1, O)]
        if y1 < Q: pcs.append(('aLz', y1, Q, aLz))
        pcs.append(('aT', max(y1, Q), inf, aT))
        return [(pcs, -1)]
    def minus_u0(y2):
        pcs = [('A0', 0, y2, A0)]
        if y2 < Q: pcs.append(('Lz', y2, Q, Lz))
        pcs.append(('T', max(y2, Q), inf, T))
        return [(pcs, -1)]
    G = lambda groups, lo=0.0, hi=inf: [(p, s, lo, hi) for p, s in groups]
    r = {}
    r['Sa'] = worst_sum(H, G(u0plus), G(minus_u1(y1)), -lam, n)
    r['Sb'] = worst_sum(H, G(u1plus), G(minus_u0(y2)), -lam, n)
    r['C0'] = worst_sum(H, G(u0plus, 0.5), G(u0plus, 0.5), -lam - 2 * g, n)
    r['C1'] = worst_sum(H, G(u1plus, 0.5), G(u1plus, 0.5), -lam - 2 * g, n)
    zero = [([('-', 0, inf, (0.0, 0.0))], +1)]
    r['N'] = worst_sum(H, G(zero, 0.0, 0.5), G(zero), -1.0, n)
    return r

def window(H, which, grid, n=361):
    """Feasible split values (Sa for y1, Sb for y2) on the grid."""
    ok = []
    for y in grid:
        r = conds_one(H, which, y, n)
        if r >= -1e-12: ok.append(y)
    return ok

def conds_one(H, which, y, n):
    lam, a, b, g, Lz, T = lines(H)
    A0 = (2 * a, b); O = (a, b); aLz = (a + Lz[0], Lz[1]); aT = (a + T[0], T[1]); inf = math.inf
    if which == 'Sa':
        X = [([('A0', 0, inf, A0)], +1, 0.0, inf), ([('Lz', 0, Q, Lz), ('T', Q, inf, T)], +1, 0.0, inf)]
        pcs = [('O', 0, y, O)] + ([('aLz', y, Q, aLz)] if y < Q else []) + [('aT', max(y, Q), inf, aT)]
    else:
        X = [([('O', 0, inf, O)], +1, 0.0, inf), ([('aLz', 0, Q, aLz), ('aT', Q, inf, aT)], +1, 0.0, inf)]
        pcs = [('A0', 0, y, A0)] + ([('Lz', y, Q, Lz)] if y < Q else []) + [('T', max(y, Q), inf, T)]
    return worst_sum(H, X, [(pcs, -1, 0.0, inf)], -lam, n)

if __name__ == '__main__':
    Hs = [float(h) for h in sys.argv[1:]] or [1.849, 1.852, 1.855, 1.858, 1.861, 1.864, 1.867, 1.87]
    grid = np.round(np.arange(0.5, 2.0001, 1 / 32), 6)
    for H in Hs:
        r = conds(H, 1.25, 0.8, n=721)
        wa = window(H, 'Sa', grid); wb = window(H, 'Sb', grid)
        fmt = lambda w: f"[{min(w):.4f}, {max(w):.4f}] ({len(w)} pts)" if w else "none"
        print(f"H={H}: C0 {r['C0']:+.1e} C1 {r['C1']:+.1e} N {r['N']:+.1e} | y1 ok {fmt(wa)} | y2 ok {fmt(wb)}", flush=True)
