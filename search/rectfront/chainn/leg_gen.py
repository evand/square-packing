#!/usr/bin/env python3
"""Pair table for the band theorem, legibility version (2026-10-25).

Every band condition is separable: kap(x, -y) = C + S + kh(C, S, x) + kh(S, C, y), kh(c, s, z) = (s - H/2) z + (c z - s)^+.
For one square and one line, the potential plus kh is a min of convex one-kink pieces (or minus: a max of concave ones), so
its extreme over the line's slope is at the start (0 or 1/2), at the potential's breakpoints, or at the ALIGNED slope
(t = S/C on the left, t' = C/S on the right: the line runs along a side of the square).  So each condition reduces to
finitely many PAIRS (left candidate, right candidate); at each pair the potentials are single lines once we know which
piece is active (region), and the max atoms are split by sign.  Each sub-case: 0 <= P(C, S, H, r) under region
hypotheses; P is the pair inequality times positive denominators.

  u0: A0 = 2a + b z on [0, x2],  Lz = a + b(2 - z) on [x2, 8/7],  T = T0 + T1 z on [8/7, inf)
  u1: A1 = a + b z on [0, y1],   aT = a + T on [y1, inf)          (x2 <= 8/7 <= y1 on the band)
"""
import math, itertools
import numpy as np
import sympy as sp
C, S, H, r = sp.symbols('C S H r', real=True)
H9 = (14 * r - 5) / 8; H87 = sp.Rational(187, 100)
lam = 2 * r - H; a = 1 - lam; b = 1 - H / 2; g = a + b
Q = sp.Rational(8, 7); T1 = H / 16 - sp.Rational(1, 10); T0 = a + b * (2 - Q) - T1 * Q
LINES = {'A0': (2 * a, b), 'Lz': (a + 2 * b, -b), 'T': (T0, T1), 'A1': (a, b), 'aT': (a + T0, T1)}
PT = {'0': (sp.Integer(0), sp.Integer(1)), '1/2': (sp.Integer(1), sp.Integer(2)), '8/7': (sp.Integer(8), sp.Integer(7)),
      'x2': (2 * b - a, 2 * b), 'y1': (T0, b - T1), 't': (S, C), "t'": (C, S)}
U = {'u0': [('A0', None, 'x2'), ('Lz', 'x2', '8/7'), ('T', '8/7', None)], 'u1': [('A1', None, 'y1'), ('aT', 'y1', None)]}
# the piece active at a fixed point (boundary points: either piece; pick the lower one)
FIXED = {('u0', '0'): 'A0', ('u0', '1/2'): 'A0', ('u0', 'x2'): 'A0', ('u0', '8/7'): 'Lz', ('u0', 'y1'): 'T',
         ('u1', '0'): 'A1', ('u1', '1/2'): 'A1', ('u1', '8/7'): 'A1', ('u1', 'x2'): 'A1', ('u1', 'y1'): 'A1'}
CONDS = {  # (left potential, left points, right potential, right points, kind)
    'Sa': ('u0', ['0', '8/7', 't'], 'u1', ['0', "t'", 'y1'], 'step'),
    'Sb': ('u1', ['0', '8/7', 't'], 'u0', ['0', 'x2', '8/7', "t'"], 'step'),
    'C0': ('u0', ['1/2', '8/7', 't'], 'u0', ['1/2', '8/7', "t'"], 'cheap'),
    'C1': ('u1', ['1/2', '8/7', 't'], 'u1', ['1/2', '8/7', "t'"], 'cheap')}

def le(z1, z2):
    """z1 <= z2 for fractions, as 'expr >= 0' (denominators positive)."""
    (n1, d1), (n2, d2) = PT[z1], PT[z2]
    return sp.expand(n2 * d1 - n1 * d2)

def pieces_at(u, z):
    if z in ('t', "t'"):
        out = []
        for nm, lo, hi in U[u]:
            hyps = []
            if lo: hyps.append(('reg', le(lo, z)))
            if hi: hyps.append(('reg', le(z, hi)))
            out.append((nm, hyps))
        return out
    return [(FIXED[(u, z)], [])]

def side(z, c, s, nm, sign):
    """d * (sign * line(z) + kh(c, s, z)) as (poly without max, max-arg or None, d)."""
    n, d = PT[z]; pc, ps = LINES[nm]
    lin = sign * (pc * d + ps * n) + (s - H / 2) * n
    marg = sp.expand(c * n - s * d)
    if (z == 't' and (c, s) == (C, S)) or (z == "t'" and (c, s) == (S, C)): marg = None   # aligned: kink, atom = 0
    if z == '0': marg = None   # c*0 - s*1 = -s <= 0
    return sp.expand(lin), marg, d

def subcases():
    out = []
    for cond, (uf, xs, ug, ys, kind) in CONDS.items():
        for zx, zy in itertools.product(xs, ys):
            for (px, hx), (py, hy) in itertools.product(pieces_at(uf, zx), pieces_at(ug, zy)):
                Lx, mx, dx = side(zx, C, S, px, +1)
                Ly, my, dy = side(zy, S, C, py, -1 if kind == 'step' else +1)
                const = C + S - lam - (2 * g if kind == 'cheap' else 0)
                hyps = hx + hy
                if kind == 'cheap':
                    if zx == 't': hyps.append(('reg', le('1/2', 't')))
                    if zy == "t'": hyps.append(('reg', le('1/2', "t'")))
                for sx in ([None] if mx is None else [+1, -1]):
                    for sy in ([None] if my is None else [+1, -1]):
                        P = Lx * dy + Ly * dx + const * dx * dy
                        h2 = list(hyps)
                        if sx == +1: P += mx * dy; h2.append(('max', mx))
                        if sx == -1: h2.append(('max', -mx))
                        if sy == +1: P += my * dx; h2.append(('max', my))
                        if sy == -1: h2.append(('max', -my))
                        tag = f"{cond}_{zx}_{zy}_{px}_{py}" + (f"_{'p' if sx > 0 else 'n'}" if sx else '') + (f"{'p' if sy > 0 else 'n'}" if sy else '')
                        out.append(dict(cond=cond, zx=zx, zy=zy, px=px, py=py, P=sp.expand(P), hyps=h2, tag=tag))
    return out

def numcheck(sc, nth=721, nH=7):
    """min of P / (dx dy) over the region (sampled); None if region empty."""
    f = sp.lambdify((C, S, H, r), sc['P']); hs = [sp.lambdify((C, S, H, r), e) for _, e in sc['hyps']]
    R = math.sqrt(2); best = None
    for Hv in np.linspace((14 * R - 5) / 8, 1.87, nH):
        for th in np.linspace(1e-7, math.pi / 2 - 1e-7, nth):
            c, s = math.cos(th), math.sin(th)
            if all(h(c, s, Hv, R) >= -1e-12 for h in hs):
                v = f(c, s, Hv, R)
                if best is None or v < best[0]: best = (v, Hv, math.degrees(th))
    return best

if __name__ == '__main__':
    scs = subcases(); print(len(scs), 'sub-cases')
    n_empty = 0
    for sc in scs:
        res = numcheck(sc)
        if res is None: n_empty += 1; continue
        flag = '  <-- NEG' if res[0] < -1e-9 else ('  tight' if res[0] < 1e-9 else '')
        if res[0] < 2e-3: print(f"  {sc['tag']:28s} min {res[0]:+.2e} H={res[1]:.4f} {res[2]:.2f}deg{flag}")
    print('empty regions:', n_empty)
