#!/usr/bin/env python3
"""Band theorem: breakpoint pairs and their margins (2026-10-25).

Each condition is min_x Phi(x) + min_y Psi(y) >= 0 with Phi, Psi piecewise linear: on each potential piece the
function is affine + (Cx - S)^+ (resp. (Sy - C)^+), so its minimum over the piece is at a piece end or at the frame kink
S/C (C/S) when inside.  Potentials by active piece (exact, region-based):
  u0: A0 on [0, x2], Lz on [x2, xLT], T on [xLT, inf)            x2 = 1 - a/(2b), xLT = Lz/T crossing
  u1: O on [0, y1], a + L beyond (a + Lz up to xLT if y1 < xLT, then a + T)   y1 = max(1 + a/2b, T0/(b - T1))
Conditions (opposite slopes, c = kap - lam):  Sa: c >= u1(y) - u0(x);  Sb: c >= u0(y) - u1(x);
  N: x <= 1/2 => kap >= 1;  C_P: x, y >= 1/2 => c + u_P(x) + u_P(y) >= 2g.
For every (x-candidate, y-candidate) pair: least value over frames and H (only frames where both candidates exist).
Pairs with small least value need exact algebra in Lean; the rest have slack."""
import math, sys
import numpy as np
from band_dec import consts

def pieces(H):
    lam, a, b, g, T0, T1 = consts(H)
    x2 = 1 - a / (2 * b); xLT = (a + 2 * b - T0) / (b + T1); y1 = max(1 + a / (2 * b), T0 / (b - T1))
    A0 = (2 * a, b); Lz = (a + 2 * b, -b); T = (T0, T1); O = (a, b)
    add = lambda p, s: (p[0] + s, p[1])
    u0 = [('A0', 0.0, x2, A0), ('Lz', x2, xLT, Lz), ('T', xLT, math.inf, T)]
    if y1 < xLT:
        u1 = [('O', 0.0, y1, O), ('aLz', y1, xLT, add(Lz, a)), ('aT', xLT, math.inf, add(T, a))]
    else:
        u1 = [('O', 0.0, y1, O), ('aT', y1, math.inf, add(T, a))]
    return dict(lam=lam, a=a, b=b, g=g, u0=u0, u1=u1)

def cands(lin, kink, sign, pcs, lo=0.0, hi=math.inf):
    """Candidates for min over x in [lo, hi] of lin(x) + kink-part + sign * piece(x).
    lin = (const, slope) of the frame part without the kink, kink = (Ck, Sk) for (Ck x - Sk)^+.
    Returns list of (label, value) and whether the function decreases to -inf on an unbounded piece."""
    out = []
    Ck, Sk = kink; k = Sk / Ck if Ck > 1e-15 else math.inf
    for nm, l, h, (p0, p1) in pcs:
        l = max(l, lo); h = min(h, hi)
        if l > h: continue
        f = lambda x: lin[0] + lin[1] * x + max(Ck * x - Sk, 0) + sign * (p0 + p1 * x)
        out.append((f'{nm}@lo', f(l)))
        if h < math.inf: out.append((f'{nm}@hi', f(h)))
        else:
            slope = lin[1] + Ck + sign * p1
            if slope < -1e-12: out.append((f'{nm}@inf', -math.inf))
        if l < k < h: out.append((f'{nm}@kink', f(k)))
    return out

def conds(H, C, S, P):
    lam, a, b, g = P['lam'], P['a'], P['b'], P['g']
    p = C + S
    fx = (p / 2, S - H / 2); kx = (C, S)       # half of the constant p on each side
    gy = (p / 2, C - H / 2); ky = (S, C)
    res = {}
    def pair(name, X, Y, const=0.0):
        for lx, vx in X:
            for ly, vy in Y:
                res[(name, lx, ly)] = vx + vy + const
    zero = [('none', 0.0, 0.0, math.inf, (0.0, 0.0))]
    pair('Sa', cands(fx, kx, +1, P['u0']), cands(gy, ky, -1, P['u1']), -lam)
    pair('Sb', cands(fx, kx, +1, P['u1']), cands(gy, ky, -1, P['u0']), -lam)
    pair('N', cands(fx, kx, +1, [('-', 0.0, math.inf, (0.0, 0.0))], hi=0.5),
         cands(gy, ky, +1, [('-', 0.0, math.inf, (0.0, 0.0))]), -1.0)
    pair('C0', cands(fx, kx, +1, P['u0'], lo=0.5), cands(gy, ky, +1, P['u0'], lo=0.5), -lam - 2 * g)
    pair('C1', cands(fx, kx, +1, P['u1'], lo=0.5), cands(gy, ky, +1, P['u1'], lo=0.5), -lam - 2 * g)
    return res

if __name__ == '__main__':
    worst = {}
    for H in np.linspace(1.849, 1.8699, 8):
        P = pieces(H)
        for t in np.linspace(0, math.pi / 2, 1801):
            C, S = math.cos(t), math.sin(t)
            if C + S > H: continue
            for k, v in conds(H, C, S, P).items():
                if k not in worst or v < worst[k][0]:
                    worst[k] = (v, H, math.degrees(t))
    by = {}
    for (nm, lx, ly), (v, H, td) in worst.items():
        by.setdefault(nm, []).append((v, lx, ly, H, td))
    tot = 0
    for nm, rows in by.items():
        rows.sort()
        tight = [r for r in rows if r[0] < 0.01]
        tot += len(rows)
        print(f"{nm}: {len(rows)} pairs, {len(tight)} with least value < 0.01, min {rows[0][0]:+.2e}")
        for v, lx, ly, H, td in tight:
            print(f"     {v:+.2e}  x {lx:10s} y {ly:10s}  (H {H:.4f}, tilt {td:.2f})")
    print("total pairs", tot)
