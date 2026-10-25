#!/usr/bin/env python3
"""Band theorem generator (2026-10-25): odd n, H*(9) <= H <= 187/100, H*(9) = (14 sqrt2 - 5)/8.

Opposite slopes x, y >= 0 (positive / negative line), frame (C, S):
  kap(x, -y) = Fx(x) + Gy(y),  Fx = p/2 + (S - H/2) x + max(C x - S, 0),  Gy = p/2 + (C - H/2) y + max(S y - C, 0).
Potentials, exact region pieces (x = |m|):
  u0: A0 = 2a + b x on [0, x2],  Lz = a + b(2 - x) on [x2, Q],  T on [Q, inf)      x2 = (2b - a)/(2b), Q = 8/7
  u1: O = a + b x on [0, y1],    aT = a + T on [y1, inf)                         y1 = T0/(b - T1)
  T(x) = Lz(Q) + T1 (x - Q), T1 = H/16 - 1/10.
Conditions (each: 0 <= Phi(x) + Psi(y) on every cell = x-piece x y-piece):
  Sa: Fx + u0(x) - lam + Gy - u1(y);  Sb: Fx + u1(x) - lam + Gy - u0(y);
  C0 (x, y >= 1/2): Fx + u0(x) + Gy + u0(y) - lam - 2g;  C1 likewise with u1;  N (x <= 1/2): Fx + Gy - 1.
Cell proof (Lean `plx`, `plx_inf`): a convex piece f = A + B x + max(Cx - S, 0) on [ln/ld, hn/hd] is >= 0 (times m,
plus E) if it is at both ends and at the kink when inside; unbounded: lower end, kink, and slope B + C >= 0.
Nested x then y, the leaves are polynomial inequalities with at most two max atoms.

This module builds the leaves (sympy), splits the max atoms, and checks every sub-leaf numerically."""
import math, sys, itertools
import numpy as np
import sympy as sp

C, S, H, r = sp.symbols('C S H r', real=True)
R2 = math.sqrt(2)
H9 = (14 * r - 5) / 8
H87 = sp.Rational(187, 100)
lam = 2 * r - H; a = 1 - lam; b = 1 - H / 2; g = a + b
Q = sp.Rational(8, 7)
T1 = H / 16 - sp.Rational(1, 10)
LzQ = a + b * (2 - Q)
T0 = LzQ - T1 * Q
p = C + S

LINES = {'A0': (2 * a, b), 'Lz': (a + 2 * b, -b), 'T': (T0, T1), 'O': (a, b), 'aT': (a + T0, T1), '0': (0, 0)}
X2 = (2 * b - a, 2 * b)          # fractions (numerator, denominator), denominator > 0 on the band
Y1 = (T0, b - T1)
ZERO = (sp.Integer(0), sp.Integer(1)); HALF = (sp.Integer(1), sp.Integer(2)); QF = (sp.Integer(8), sp.Integer(7))
U0 = [('A0', ZERO, X2), ('Lz', X2, QF), ('T', QF, None)]
U1 = [('O', ZERO, Y1), ('aT', Y1, None)]

def clip(pieces, lo=None, hi=None):
    """restrict pieces to [lo, hi] (lo, hi fixed fractions known to lie in the first / last piece)"""
    out = []
    for nm, l, h in pieces:
        out.append([nm, l, h])
    if lo is not None: out[0][1] = lo
    if hi is not None:
        out = [out[0]]; out[0][2] = hi
    return [tuple(x) for x in out]

CONDS = {
    'Sa': (U0, +1, U1, -1, -lam, None, None),
    'Sb': (U1, +1, U0, -1, -lam, None, None),
    'C0': (clip(U0, lo=HALF), +1, clip(U0, lo=HALF), +1, -lam - 2 * g, None, None),
    'C1': (clip(U1, lo=HALF), +1, clip(U1, lo=HALF), +1, -lam - 2 * g, None, None),
    'N': ([('0', ZERO, HALF)], +1, [('0', ZERO, None)], +1, sp.Integer(-1), None, None),
}

def xform(piece, sign, const):
    """x side: f = A + B x + max(Cx - S, 0)"""
    nm, l, h = piece; pc, ps = LINES[nm]
    return p / 2 + sign * pc + const, S - H / 2 + sign * ps

def yform(piece, sign):
    nm, l, h = piece; pc, ps = LINES[nm]
    return p / 2 + sign * pc, C - H / 2 + sign * ps

class Leaf:
    """0 <= expr, with max atoms [(arg, ...)], under extra hypotheses (list of (expr >= 0 / > 0 strings))."""
    def __init__(self, cond, cell, tag, expr, maxes, hyps):
        self.cond, self.cell, self.tag, self.expr, self.maxes, self.hyps = cond, cell, tag, expr, maxes, hyps

M1, M2 = sp.symbols('M1 M2', real=True)   # max atoms

def cell_leaves(cond, i, j, xp, yp, xs, ys, const):
    """Nested plx: x first (m = 1, E = Psi), then y."""
    A, B = xform(xp, xs, const); A2, B2 = yform(yp, ys)
    _, (ln, ld), xh = xp; _, (lyn, lyd), yh = yp
    out = []
    # x-ends: E_x = A ld + B ln + max(C ln - S ld, 0) (scaled by ld);  x-kink: E_x = A C + B S (scaled by C)
    xends = [('xlo', ld, A * ld + B * ln + M1, [(M1, C * ln - S * ld)], [])]
    if xh is not None:
        hn, hd = xh
        xends.append(('xhi', hd, A * hd + B * hn + M1, [(M1, C * hn - S * hd)], []))
        xends.append(('xk', C, A * C + B * S, [], [('pos', C), ('ge', S * ld - ln * C), ('ge', hn * C - S * hd)]))
    else:
        xends.append(('xk', C, A * C + B * S, [], [('pos', C), ('ge', S * ld - ln * C)]))
    for xt, mx, Ex, mxs, hx in xends:
        # y-lemma with m = mx, E = Ex:  0 <= mx (A2 + B2 y + max(S y - C, 0)) + Ex
        yends = [('ylo', lyd, mx * (A2 * lyd + B2 * lyn + M2) + Ex * lyd, [(M2, S * lyn - C * lyd)], [])]
        if yh is not None:
            yn, yd = yh
            yends.append(('yhi', yd, mx * (A2 * yd + B2 * yn + M2) + Ex * yd, [(M2, S * yn - C * yd)], []))
            yends.append(('yk', S, mx * (A2 * S + B2 * C) + Ex * S, [], [('pos', S), ('ge', C * lyd - lyn * S), ('ge', yn * S - C * yd)]))
        else:
            yends.append(('yk', S, mx * (A2 * S + B2 * C) + Ex * S, [], [('pos', S), ('ge', C * lyd - lyn * S)]))
        for yt, _, ex, mys, hy in yends:
            out.append(Leaf(cond, (i, j), f'{xt}_{yt}', ex, mxs + mys, hx + hy))
        if yh is None:   # slope condition for the unbounded y piece: 0 <= B2 + S  (times mx > 0)
            out.append(Leaf(cond, (i, j), f'{xt}_yslope', B2 + S, [], hx))
    if xh is None:
        out.append(Leaf(cond, (i, j), 'xslope', B + C, [], []))
    return out

def all_leaves():
    L = []
    for cond, (xps, xs, yps, ys, const, _, _) in CONDS.items():
        for i, xp in enumerate(xps):
            for j, yp in enumerate(yps):
                L += cell_leaves(cond, i, j, xp, yp, xs, ys, const)
    return L

def subcases(leaf):
    """split each max atom M = max(arg, 0): case 'pos' (arg >= 0, M = arg) or 'neg' (arg <= 0, M = 0)."""
    out = []
    for choice in itertools.product((0, 1), repeat=len(leaf.maxes)):
        e = leaf.expr; hyps = list(leaf.hyps)
        for (Mv, arg), ch in zip(leaf.maxes, choice):
            if ch: e = e.subs(Mv, arg); hyps.append(('ge', arg))
            else: e = e.subs(Mv, 0); hyps.append(('ge', -arg))
        out.append((choice, sp.expand(e), hyps))
    return out

def numeric_check(expr, hyps, nt=721, nH=5):
    f = sp.lambdify((C, S, H, r), expr, 'numpy')
    hf = [(k, sp.lambdify((C, S, H, r), h, 'numpy')) for k, h in hyps]
    t = np.linspace(0, math.pi / 2, nt); Cv, Sv = np.cos(t), np.sin(t)
    best = (math.inf, None)
    H9v = (14 * R2 - 5) / 8
    for Hv in np.linspace(H9v, 1.87, nH):
        val = np.broadcast_to(f(Cv, Sv, Hv, R2), Cv.shape).astype(float)
        ok = np.ones_like(Cv, bool)
        for k, hfun in hf:
            hv = np.broadcast_to(hfun(Cv, Sv, Hv, R2), Cv.shape)
            ok &= (hv > 1e-12) if k == 'pos' else (hv >= -1e-12)
        if ok.any():
            k = np.argmin(np.where(ok, val, np.inf))
            if val[k] < best[0]: best = (val[k], (Hv, math.degrees(t[k])))
    return best

if __name__ == '__main__':
    L = all_leaves()
    print(len(L), 'leaves')
    n_sub = 0; worst = []
    for lf in L:
        for ch, e, hy in subcases(lf):
            n_sub += 1
            v, at = numeric_check(e, hy)
            worst.append((v, lf.cond, lf.cell, lf.tag, ch, at))
    worst.sort(key=lambda w: w[0])
    print(n_sub, 'sub-leaves;', sum(1 for w in worst if w[0] == math.inf), 'empty regions')
    for w in worst[:25]:
        print(f"  {w[0]:+.3e} {w[1]} cell {w[2]} {w[3]} max-case {w[4]} at {w[5]}")
    print("negative:", sum(1 for w in worst if w[0] < -1e-9))
