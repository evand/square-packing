#!/usr/bin/env python3
"""Final5 aligned leaves (x = a1): numerator polynomials and Handelman LP certificates (2026-10-21).
For each leaf: E = L_bound - target = N/D on its region; N' = sign(D) N must be >= 0.  Region: 0 <= u, v <= 1,
mu's sign, the kink's sign.  Certificate: N' = sum lam_a prod g_i^{a_i}, lam >= 0, generators g = {u, 1-u, v, 1-v,
region factors}.  sqrt 2 = r: certify at r_lo and r_hi (linear in r)."""
import itertools, math, sys
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from leaves5 import leaf, u, v, H, r

R_LO, R_HI = sp.Rational(141421, 100000), sp.Rational(141422, 100000)
H0 = (6 * r - 1) / 4

def numer(expr):
    n, d = sp.fraction(sp.together(expr))
    return sp.expand(n), sp.factor(d)

def sample_sign(poly, region, pts=200):
    """sign of poly on the region by sampling (for the denominator, which has no zero inside)."""
    f = sp.lambdify((u, v), poly)
    rng = np.random.default_rng(0); sgn = set()
    for _ in range(pts):
        a, b = rng.uniform(0.01, 0.99, 2)
        if all(sp.lambdify((u, v), g)(a, b) >= 0 for g in region):
            val = f(a, b)
            if abs(val) > 1e-12: sgn.add(val > 0)
    return sgn

def _polyarr(expr, dmax):
    a = np.zeros((dmax + 1, dmax + 1))
    for (i, j), c in sp.Poly(expr, u, v).terms(): a[i, j] = float(c)
    return a

def _mul(a, b, dmax):
    out = np.zeros((dmax + 1, dmax + 1))
    for i, j in zip(*np.nonzero(a)):
        out[i:, j:] += a[i, j] * b[:dmax + 1 - i, :dmax + 1 - j]
    return out

def _exps(degs, deg):
    """exponent vectors with sum a_i deg_i <= deg."""
    if not degs:
        yield ()
        return
    for a in range(deg // degs[0] + 1):
        for rest in _exps(degs[1:], deg - a * degs[0]):
            yield (a,) + rest

def handelman(N, gens, deg):
    """LP: N = sum lam * prod gens^alpha (degree <= deg), lam >= 0.  Returns (ok, certificate)."""
    dN = sp.Poly(N, u, v).total_degree()
    dmax = max(deg, dN)
    garr = [_polyarr(g, dmax) for g in gens]
    gdeg = [max(1, sp.Poly(g, u, v).total_degree()) for g in gens]
    cols, alphas = [], []
    for alpha in _exps(gdeg, max(deg, dN)):
        p = np.zeros((dmax + 1, dmax + 1)); p[0, 0] = 1
        for g, a in zip(garr, alpha):
            for _ in range(a): p = _mul(p, g, dmax)
        cols.append(p.ravel()); alphas.append(alpha)
    A = np.array(cols).T
    b = _polyarr(N, dmax).ravel()
    keep = np.any(A != 0, axis=1) | (b != 0)
    res = linprog(np.ones(A.shape[1]), A_eq=A[keep], b_eq=b[keep], bounds=[(0, None)] * A.shape[1], method='highs')
    if res.status != 0: return False, None
    return True, [(alphas[i], res.x[i]) for i in range(len(alphas)) if res.x[i] > 1e-12]

def regions(e1, e2, mu_pos, kink_pos):
    """region polynomials (>= 0) in u, v for mu's sign and the kink's sign (x = a1)."""
    C1, S1 = (1 - u**2), 2 * u      # times (1 + u^2)
    C2, S2 = (1 - v**2), 2 * v
    a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
    v1 = C1 if e1 else -C1; v2 = C2 if e2 else -C2
    mu_num = sp.factor(v2 * S1 - v1 * S2)              # sign of mu (S's positive)
    kink = sp.factor(a2[0] * a1[1] - a2[1] * a1[0])     # cr a2 a1 (up to positive factors)
    out = [mu_num if mu_pos else -mu_num, kink if kink_pos else -kink]
    return out

if __name__ == '__main__':
    deg = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    base = [u, 1 - u, v, 1 - v]
    stats = []
    for e1, e2, b1, b2, mp, kp in itertools.product((False, True), (False, True), 'UP', 'UP', (True, False), (True, False)):
        reg = regions(e1, e2, mp, kp)
        Lb = leaf(e1, e2, 1, b1, b2, mp, kp)
        for Hname, Hv, tgt in [('1', 1, 5), ('h0', H0, 5), ('2', 2, 6 * r - 4)]:
            N, D = numer(Lb.subs(H, Hv) - tgt)
            sg = sample_sign(sp.expand(D).subs(r, math.sqrt(2)), reg)
            if not sg: stats.append((e1, e2, b1, b2, mp, kp, Hname, 'empty')); continue
            if len(sg) > 1: stats.append((e1, e2, b1, b2, mp, kp, Hname, 'den-sign?')); continue
            Np = N if sg == {True} else -N
            # factor out region generators: try to certify N' at r_lo and r_hi
            gens = base + [g for g in reg]
            # use linear/irreducible factors of the region polys as generators
            fac = []
            for g in reg:
                for f, m in sp.factor_list(g)[1]:
                    fac.append(f)
            gens = base + [f for f in fac if sp.Poly(f, u, v).total_degree() >= 1]
            ok = True
            for rv in (R_LO, R_HI):
                Nr = sp.expand(Np.subs(r, rv))
                if not handelman(Nr, gens, deg)[0]: ok = False; break
            stats.append((e1, e2, b1, b2, mp, kp, Hname, 'OK' if ok else 'FAIL'))
            print(stats[-1], flush=True)
    from collections import Counter
    print(Counter(s[-1] for s in stats))
