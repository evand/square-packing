#!/usr/bin/env python3
"""Find a Positivstellensatz-style certificate by LP:  target = sum_i lam_i * prod_i + (C^2+S^2-1) q1 + (r^2-2) q2,
lam_i >= 0, prod_i = products (degree <= D) of nonnegative atoms; q1, q2 free polynomials of bounded degree.
Exact rational solution attempted afterwards (linear solve on the support)."""
import itertools, sys
import sympy as sp
import numpy as np
from scipy.optimize import linprog

def find(target, atoms, gens, D=2, qdeg=3, extra_eqs=()):
    prods = []
    for d in range(0, D + 1):
        for combo in itertools.combinations_with_replacement(range(len(atoms)), d):
            e = sp.Integer(1)
            for i in combo: e *= atoms[i][1]
            prods.append((" * ".join(atoms[i][0] for i in combo) or "1", sp.expand(e)))
    eqs = [sp.expand(e) for e in extra_eqs]
    # free multipliers for equalities: all monomials up to qdeg
    mons = sorted(sp.itermonomials(gens, qdeg), key=sp.default_sort_key)
    cols = [p for _, p in prods] + [sp.expand(m * e) for e in eqs for m in mons] + [sp.expand(-m * e) for e in eqs for m in mons]
    nlam = len(prods)
    allp = [sp.Poly(c, *gens) for c in cols] + [sp.Poly(target, *gens)]
    monset = sorted(set(m for p in allp for m in p.monoms()))
    idx = {m: i for i, m in enumerate(monset)}
    A = np.zeros((len(monset), len(cols)))
    for j, p in enumerate(allp[:-1]):
        for m, c in p.terms(): A[idx[m], j] = float(c)
    bvec = np.zeros(len(monset))
    for m, c in allp[-1].terms(): bvec[idx[m]] = float(c)
    res = linprog(np.r_[np.ones(nlam), np.zeros(len(cols) - nlam)] * 0 + np.r_[np.ones(nlam), 1e-3 * np.ones(len(cols) - nlam)],
                  A_eq=A, b_eq=bvec, bounds=[(0, None)] * len(cols), method='highs')
    if res.status != 0:
        return None
    lam = res.x[:nlam]
    used = [(prods[i][0], lam[i]) for i in range(nlam) if lam[i] > 1e-9]
    return used

if __name__ == '__main__':
    C, S, H, r = sp.symbols('C S H r')
    d = 2 - H; e = 1 + 2 * r - 2 * H
    target = sp.expand(S * d * (C + S - 1) + e * S * (S - 1) + C * d * (C + 1 - H))
    atoms = [("C", C), ("S", S), ("(1-C)", 1 - C), ("(1-S)", 1 - S), ("(C+S-1)", C + S - 1),
             ("(H-2r+1)", H - 2 * r + 1), ("(2-H)", 2 - H), ("e", e), ("u", S * d - e * C), ("(r-1)", r - 1)]
    for D in (2, 3):
        out = find(target, atoms, (C, S, H, r), D=D, qdeg=2, extra_eqs=(C**2 + S**2 - 1, r**2 - 2))
        print("D =", D, ":", out)
        if out: break
