#!/usr/bin/env python3
"""Symbolic cell analysis for (S) (and similar): in each cell the target f is linear in the slopes; enumerate the cell
polygon's vertices symbolically, evaluate f there, and report each vertex value's minimum over (t, H) numerically.
Also reports recession directions (unbounded cells) and the coefficient signs."""
import itertools, math, sys
import numpy as np
import sympy as sp

C, S, H, r = sp.symbols('C S H r', positive=True)   # r = sqrt2
m, y = sp.symbols('m y', real=True)
b = 1 - H / 2; Z = 2 + 2 * r - 2 * H; k = (C + S) / 2; c0 = H / 2 - k

def kap_opp(s1, s2):
    """m >= 0, m' = -y <= 0.  s1 = sign(mC - S), s2 = sign(C - yS)."""
    return (C + m * S + s1 * (m * C - S) + s2 * (C - y * S) + (y * C + S)) / 2 - (m + y) * c0

def kap_same(s1, s2):
    """0 <= m' = y <= m.  s1 = sign(mC - S), s2 = sign(yC - S)."""
    return (C + m * S + s1 * (m * C - S) + C + y * S + s2 * (y * C - S)) / 2 - (m - y) * c0

def numeric_min(expr, feas=(), Hlo=2 * math.sqrt(2) - 1):
    f = sp.lambdify((C, S, H, r), expr, 'numpy')
    fs = [sp.lambdify((C, S, H, r), g, 'numpy') for g in feas]
    ts = np.linspace(1e-4, math.pi / 2 - 1e-4, 1500); Hs = np.linspace(Hlo, 2, 60)
    T, HH = np.meshgrid(ts, Hs)
    Cv, Sv = np.cos(T), np.sin(T)
    with np.errstate(all='ignore'):
        v = f(Cv, Sv, HH, math.sqrt(2)) * np.ones_like(T)
    ok = Cv + Sv <= HH
    with np.errstate(all='ignore'):
        for g in fs:
            ok &= (g(Cv, Sv, HH, math.sqrt(2)) * np.ones_like(T)) >= -1e-9
    if not ok.any(): return (math.inf, None, None)
    v = np.where(ok, v, np.nan)
    i = np.nanargmin(v)
    return np.nanmin(v), math.degrees(T.flat[i]), HH.flat[i]

def analyse(name, f, cons):
    """cons: list of linear forms >= 0 in (m, y)."""
    print(f"== {name}")
    fm, fy = sp.simplify(sp.diff(f, m)), sp.simplify(sp.diff(f, y))
    print(f"   df/dm = {sp.factor(fm)},  df/dy = {sp.factor(fy)}")
    verts = []
    for g1, g2 in itertools.combinations(cons, 2):
        sol = sp.solve([g1, g2], [m, y], dict=True)
        if not sol: continue
        sol = sol[0]
        if m not in sol or y not in sol: continue
        verts.append(((g1, g2), sol))
    for (g1, g2), sol in verts:
        val = sp.simplify(f.subs(sol))
        # feasibility of the vertex: all constraints >= 0 (check numerically at a few points; report)
        feas = [sp.lambdify((C, S, H, r), sp.simplify(g.subs(sol)), 'numpy') for g in cons]
        mn = numeric_min(val, [sp.simplify(g.subs(sol)) for g in cons])
        if mn[1] is None: continue
        print(f"   vertex m={sp.simplify(sol[m])}, y={sp.simplify(sol[y])}  [{g1} = 0 = {g2}]:  min {mn[0]:+.4e} at t={mn[1]:.2f}, H={mn[2]:.4f}"
              f"\n        f = {sp.factor(sp.simplify(val))}")

if __name__ == '__main__':
    for s1, s2 in itertools.product((1, -1), repeat=2):
        f = kap_opp(s1, s2)
        cons = [y, m - y, s1 * (m * C - S), s2 * (C - y * S)]
        Gb = 1 + b * (m - y); Zb = Z - b * (m + y)
        analyse(f"opp s1={s1} s2={s2} grid (G <= Zb)", sp.expand(f - Gb), cons + [Zb - Gb])
        analyse(f"opp s1={s1} s2={s2} zig  (Zb <= G)", sp.expand(f - Zb), cons + [Gb - Zb])
