#!/usr/bin/env python3
"""Vertex enumeration for the band theorem (odd n, 1.849 <= H < 1.87), 2026-10-24.

Opposite slopes m = x >= 0, m' = -y <= 0, frame (C, S) >= 0, C^2 + S^2 = 1, C + S <= H, k = (C + S)/2:
  kap = Tm(x) + Ty(y) + beta (x - y)                                   (identity, StripC.kap_C)
  Tm(x) = (C + x S + |x C - S|)/2 - x (1 - k)                          kink x = S/C
  Ty(y) = (|C - y S| + y C + S)/2 - y (H - 1 - k)                      kink y = C/S
Potentials in x = |m| (lincheck.exact_pots):
  A0 = 2a + beta x, Lz = a + beta (2 - x), T = 7H/12 + 14/3 - 4 sqrt2 + (H/16 - 1/10)(x - 4/3), L = max(Lz, T)
  u0 = min(A0, L), u1 = min(a + beta x, a + L).
Conditions (c = kap - lam), each separable Phi(x) + Psi(y) >= 0:
  (Sa) c >= u1(y) - u0(x):   Phi = Tm + beta x + u0(x),  Psi = Ty - beta y - lam - u1(y)
  (Sb) c >= u0(y) - u1(x):   Phi = Tm + beta x + u1(x),  Psi = Ty - beta y - lam - u0(y)
  (N)  x <= 1/2 => kap >= 1: Phi = Tm + beta x (x in [0, 1/2]), Psi = Ty - beta y - 1
  (C_P) x, y >= 1/2 => c + u_P(x) + u_P(y) >= 2g:  Phi = Tm + beta x + u_P(x) (x >= 1/2), Psi = Ty - beta y - lam + u_P(y) - 2g (y >= 1/2)
The x-swapped versions follow by kap_mirror (frame (S, C)).  Each Phi, Psi is piecewise linear; its minimum over the
range is at a breakpoint (frame kink, potential kinks, range end) unless it decreases to infinity (reported).
Output: per condition, the least value of min Phi + min Psi over a fine frame grid and H grid, and which breakpoint
pair attains it."""
import math, sys
import numpy as np

R2 = math.sqrt(2)

def pieces(H):
    lam = 2 * R2 - H; a = 1 - lam; b = 1 - H / 2; g = a + b
    T1 = H / 16 - 1 / 10; T0 = 7 * H / 12 + 14 / 3 - 4 * R2 - T1 * 4 / 3
    A0 = lambda x: 2 * a + b * x
    Lz = lambda x: a + b * (2 - x)
    T = lambda x: T0 + T1 * x
    L = lambda x: max(Lz(x), T(x))
    u0 = lambda x: min(A0(x), L(x))
    u1 = lambda x: min(a + b * x, a + L(x))
    # potential kinks: pairwise crossings of the affine pieces a+bx, 2a+bx, Lz, a+Lz, T, a+T
    lines = [(a, b), (2 * a, b), (a + 2 * b, -b), (2 * a + 2 * b, -b), (T0, T1), (a + T0, T1)]
    ks = set()
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            (p, q), (r, s) = lines[i], lines[j]
            if abs(q - s) > 1e-15:
                x = (r - p) / (q - s)
                if x > 0: ks.add(x)
    return dict(lam=lam, a=a, b=b, g=g, u0=u0, u1=u1, kinks=sorted(ks))

def Tm(C, S, H, x):
    k = (C + S) / 2
    return (C + x * S + abs(x * C - S)) / 2 - x * (1 - k)

def Ty(C, S, H, y):
    k = (C + S) / 2
    return (abs(C - y * S) + y * C + S) / 2 - y * (H - 1 - k)

def pl_min(f, cands, lo, hi_slope_pt):
    """Min of a piecewise-linear f over [lo, inf) given candidate breakpoints; also checks the slope at infinity
    (f(big) - f(big/2)); returns (value, argmin label)."""
    best = (math.inf, None)
    for lab, x in cands:
        if x >= lo - 1e-15:
            v = f(max(x, lo))
            if v < best[0]: best = (v, lab)
    X = hi_slope_pt
    if f(2 * X) < f(X) - 1e-12: best = (-math.inf, 'inf')
    return best

def conds(H, C, S):
    P = pieces(H); lam, a, b, g, u0, u1 = P['lam'], P['a'], P['b'], P['g'], P['u0'], P['u1']
    kx = [('S/C', S / C if C > 1e-12 else 1e9)] + [(f'k{i}', x) for i, x in enumerate(P['kinks'])] + [('0', 0.0), ('1/2', 0.5)]
    ky = [('C/S', C / S if S > 1e-12 else 1e9)] + [(f'k{i}', x) for i, x in enumerate(P['kinks'])] + [('0', 0.0), ('1/2', 0.5)]
    out = {}
    def sep(name, Phi, Psi, xlo=0.0, ylo=0.0, xhi=None):
        if xhi is not None:   # x in [0, xhi]: candidates clipped, no infinity
            cx = [(l, min(x, xhi)) for l, x in kx] + [('xhi', xhi)]
            vx = min((Phi(x), l) for l, x in cx if x >= xlo)
        else:
            vx = pl_min(Phi, kx, xlo, 1e4)
        vy = pl_min(Psi, ky, ylo, 1e4)
        out[name] = (vx[0] + vy[0], vx[1], vy[1])
    sep('Sa', lambda x: Tm(C, S, H, x) + b * x + u0(x), lambda y: Ty(C, S, H, y) - b * y - lam - u1(y))
    sep('Sb', lambda x: Tm(C, S, H, x) + b * x + u1(x), lambda y: Ty(C, S, H, y) - b * y - lam - u0(y))
    sep('N', lambda x: Tm(C, S, H, x) + b * x, lambda y: Ty(C, S, H, y) - b * y - 1, xhi=0.5)
    for Pn, u in (('0', u0), ('1', u1)):
        sep('C' + Pn, lambda x, u=u: Tm(C, S, H, x) + b * x + u(x),
            lambda y, u=u: Ty(C, S, H, y) - b * y - lam + u(y) - 2 * g, xlo=0.5, ylo=0.5)
    return out

if __name__ == '__main__':
    Hs = [1.849, 1.852, 1.855, 1.86, 1.865, 1.8699]
    ts = np.linspace(0, math.pi / 2, 9001)
    worst = {}
    for H in Hs:
        for t in ts:
            C, S = math.cos(t), math.sin(t)
            if C + S > H: continue
            for nm, (v, lx, ly) in conds(H, C, S).items():
                if nm not in worst or v < worst[nm][0]:
                    worst[nm] = (v, H, math.degrees(t), lx, ly)
    for nm, (v, H, td, lx, ly) in worst.items():
        print(f"{nm:3s}: min {v:+.3e} at H {H}, tilt {td:.3f} deg, x at {lx}, y at {ly}")
