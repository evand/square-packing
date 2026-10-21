#!/usr/bin/env python3
"""Certificates for Final5 leaves (2026-10-21).
A leaf: (eps1, eps2, br1, br2, cand, mu_pos, kink_pos, Hint).  In half angles u (line 1), v (line 2), with
H and r = sqrt 2 as variables, E = B1'/S1 + M + B2'/S2 and the target (5 on [1, h0], 4 + 6r - 4H on [h0, 2]).
P = numerator of E - target, made >= 0 on the region by the sign of the denominator.  Certificate: products of
generators (u, 1-u, v, 1-v, H-interval, r bounds, region factors) found by a float LP; Lean's linarith checks."""
import itertools, json, math, sys
import numpy as np
import sympy as sp
from scipy.optimize import linprog

u, v, H, r = sp.symbols('u v H r', real=True)
R_LO, R_HI = sp.Rational(1414213, 1000000), sp.Rational(1414214, 1000000)
H0 = (6 * r - 1) / 4

def CS(w):
    return (1 - w**2) / (1 + w**2), 2 * w / (1 + w**2)

def bnd(C, S, br):
    return 2 * S - C * (H - 1) if br == 'U' else (C * (S - C) + 1 + C + S - H) / C

def parts(e1, e2, br1, br2, cand, mu_pos, kink_pos):
    C1, S1 = CS(u); C2, S2 = CS(v)
    a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
    v1 = C1 if e1 else -C1; v2 = C2 if e2 else -C2
    x = {'up': (1, 0), 'a1': a1, 'a2': a2}[cand]
    cr = lambda p, q: p[0] * q[1] - p[1] * q[0]
    dt = lambda p, q: p[0] * q[0] + p[1] * q[1]
    mu = v2 / S2 - v1 / S1
    y = (x[0] + x[1]) / 2 if mu_pos else H - (x[0] + x[1]) / 2
    def absv(z, kind):           # kind: 'z' zero, '+' / '-' known sign, 'k' the kink with sign kink_pos
        return {'z': 0, '+': z, '-': -z, 'k': z if kink_pos else -z}[kind]
    k1 = cr(a1, x); k2 = cr(a2, x)
    kind1 = {'up': '-', 'a1': 'z', 'a2': 'k'}[cand]
    kind2 = {'up': '-', 'a1': 'k', 'a2': 'z'}[cand]
    h1 = absv(k1, kind1) + dt(a1, x); h2 = absv(k2, kind2) + dt(a2, x)
    M = mu * y + h1 / (2 * S1) + h2 / (2 * S2)
    E = (bnd(C1, S1, br1) + (C1 * H if e1 else 0)) / S1 + M + (bnd(C2, S2, br2) + (0 if e2 else C2 * H)) / S2
    # region conditions (>= 0), as polynomials up to positive factors
    mu_num = sp.factor(sp.together(mu * S1 * S2 * (1 + u**2) * (1 + v**2)))
    region = [mu_num if mu_pos else -mu_num]
    kink = {'up': None, 'a1': k2, 'a2': k1}[cand]
    if kink is not None:
        kn = sp.factor(sp.together(kink * (1 + u**2) * (1 + v**2)))
        region.append(kn if kink_pos else -kn)
    return E, region

def poly_factors(region):
    out = []
    for g in region:
        c, fl = sp.factor_list(sp.numer(sp.together(g)))
        sign = sp.sign(c)
        facs = [f for f, m in fl if m % 2 == 1]
        out.append((sign, facs))
    return out

if __name__ == '__main__':
    case = (False, False, 'U', 'U', 'a1', True, True)
    E, region = parts(*case)
    for Hint in ('lo', 'hi'):
        tgt = 5 if Hint == 'lo' else 4 + 6 * r - 4 * H
        N, D = sp.fraction(sp.together(E - tgt))
        print(Hint, 'deg N', sp.Poly(sp.expand(N), u, v, H, r).total_degree(), ' D =', sp.factor(D))
        print('   region:', [sp.factor(g) for g in region])

# ---------------------------------------------------------------- LP over products (dict polynomials)
VARS = (u, v, H, r)

def todict(expr):
    return {m: float(c) for m, c in sp.Poly(sp.expand(expr), *VARS).terms()}

def dmul(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = tuple(x + y for x, y in zip(ma, mb))
            out[m] = out.get(m, 0.0) + ca * cb
    return out

def gens_for(region, Hint):
    """named generators: (name, sympy expr, kind) with kind in {'uv', 'H', 'r'}."""
    G = [('u', u, 'uv'), ('1-u', 1 - u, 'uv'), ('v', v, 'uv'), ('1-v', 1 - v, 'uv')]
    for i, g in enumerate(region):
        c, fl = sp.factor_list(sp.numer(sp.together(g)))
        sgn = 1 if c > 0 else -1
        lin = [f for f, m in fl if m % 2 == 1]
        # a region condition g >= 0 with g = c * prod f^m: keep it as one generator when it is linear or
        # a product of a linear factor and a positive factor (u v + 1 > 0)
        pos = [f for f in lin if f.free_symbols and all(sp.Poly(f, u, v).coeffs()[k] > 0 for k in range(len(sp.Poly(f, u, v).coeffs()))) and sp.Poly(f, u, v).total_degree() > 1]
        rest = [f for f in lin if f not in pos]
        if len(rest) == 1:
            G.append((f'reg{i}', sp.expand(sgn * rest[0]), 'uv'))
        elif len(rest) > 1:
            G.append((f'reg{i}', sp.expand(sgn * sp.Mul(*rest)), 'uv'))
    if Hint == 'lo':
        G += [('H-1', H - 1, 'H'), ('h0-H', H0 - H, 'H')]
    else:
        G += [('H-h0', H - H0, 'H'), ('2-H', 2 - H, 'H')]
    G += [('r-rlo', r - R_LO, 'r'), ('rhi-r', R_HI - r, 'r')]
    return G

def certify(P, G, deg, margin=0.0, const_min=0.0):
    uvG = [g for g in G if g[2] == 'uv']; HG = [g for g in G if g[2] == 'H']; rG = [g for g in G if g[2] == 'r']
    uvd = [todict(g[1]) for g in uvG]; uvdeg = [max(1, sp.Poly(g[1], u, v).total_degree()) for g in uvG]
    def exps(i, budget):
        if i == len(uvG):
            yield ()
            return
        for a in range(budget // uvdeg[i] + 1):
            for rest in exps(i + 1, budget - a * uvdeg[i]):
                yield (a,) + rest
    cols, labels = [], []
    pow_cache = {}
    for alpha in exps(0, deg):
        p = {(0, 0, 0, 0): 1.0}
        for gi, a in enumerate(alpha):
            for _ in range(a): p = dmul(p, uvd[gi])
        for hg in [None] + HG:
            ph = p if hg is None else dmul(p, todict(hg[1]))
            for rg in [None] + rG:
                pr = ph if rg is None else dmul(ph, todict(rg[1]))
                cols.append(pr); labels.append((alpha, hg[0] if hg else None, rg[0] if rg else None))
    target = todict(P)
    monos = sorted({m for c in cols for m in c} | set(target))
    idx = {m: i for i, m in enumerate(monos)}
    A = np.zeros((len(monos), len(cols)))
    for j, c in enumerate(cols):
        for m, val in c.items(): A[idx[m], j] = val
    b = np.zeros(len(monos))
    for m, val in target.items(): b[idx[m]] = val
    # margin: require the constant-term slack  P - margin = combination (scale-free check of robustness)
    if margin:
        b[idx[(0, 0, 0, 0)]] -= margin if (0, 0, 0, 0) in idx else 0
    bnds = [(0, None)] * len(cols)
    if const_min:
        for j, lab in enumerate(labels):
            if not any(lab[0]) and lab[1] is None and lab[2] is None: bnds[j] = (const_min, None)
    res = linprog(np.ones(len(cols)), A_eq=A, b_eq=b, bounds=bnds, method='highs',
                  options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
    if res.status != 0:
        return None
    resid = np.abs(A @ res.x - b).max()
    used = [(labels[j], res.x[j]) for j in range(len(cols)) if res.x[j] > 1e-10]
    return used, resid

def leaf_P(case, Hint):
    E, region = parts(*case)
    tgt = 5 if Hint == 'lo' else 4 + 6 * r - 4 * H
    N, D = sp.fraction(sp.together(E - tgt))
    return sp.expand(N), sp.factor(D), region

# ---------------------------------------------------------------- Lean emission
def lean_expr(e):
    s = sp.sstr(sp.expand(e), order='lex').replace('**', '^')
    return s

def hyp_name(gname):
    return {'u': 'g_u', '1-u': 'g_1u', 'v': 'g_v', '1-v': 'g_1v', 'H-1': 'g_H1', 'h0-H': 'g_h0H',
            'H-h0': 'g_Hh0', '2-H': 'g_2H', 'r-rlo': 'g_rlo', 'rhi-r': 'g_rhi',
            'u-1/2': 'g_uh', '1/2-u': 'g_hu', 'v-1/2': 'g_vh', '1/2-v': 'g_hv'}.get(gname, 'g_' + gname.replace('-', '_').replace('/', 'o'))

def lean_hint(label, G):
    alpha, hn, rn = label
    uvG = [g for g in G if g[2] == 'uv']
    factors = []
    for g, a in zip(uvG, alpha):
        if a == 1: factors.append(hyp_name(g[0]))
        elif a > 1: factors.append(f'(pow_nonneg {hyp_name(g[0])} {a})')
    if hn: factors.append(hyp_name(hn))
    if rn: factors.append(hyp_name(rn))
    if not factors: return None
    t = factors[0]
    for f in factors[1:]: t = f'(mul_nonneg {t} {f})'
    return t

def lean_leaf(name, P, G, used):
    hyps = ' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in G)
    hints = [h for h in (lean_hint(l, G) for l, _ in used) if h]
    return (f'theorem {name} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(P)} := by\n'
            f'  linarith [{", ".join(hints)}]\n')

# ---------------------------------------------------------------- fallbacks: subdivision
def gens_box(G, box):
    """add box generators for a sub-box {u in [ua, ub], v in [va, vb]} (only the non-trivial sides)."""
    (ua, ub), (va, vb) = box
    extra = []
    if ua > 0: extra.append((f'u-{ua}', u - ua, 'uv'))
    if ub < 1: extra.append((f'{ub}-u', ub - u, 'uv'))
    if va > 0: extra.append((f'v-{va}', v - va, 'uv'))
    if vb < 1: extra.append((f'{vb}-v', vb - v, 'uv'))
    return G + extra

def certify_split(P, G, degs=(6, 8, 10)):
    """try the whole box, then halves in u, in v, then quarters.  Returns list of (box, used) or None."""
    half = sp.Rational(1, 2)
    plans = [[((0, 1), (0, 1))],
             [((0, half), (0, 1)), ((half, 1), (0, 1))],
             [((0, 1), (0, half)), ((0, 1), (half, 1))],
             [((a, b), (c, d)) for (a, b) in ((0, half), (half, 1)) for (c, d) in ((0, half), (half, 1))]]
    for plan in plans:
        out = []
        for box in plan:
            Gb = gens_box(G, box)
            got = None
            for deg in degs:
                res = certify(P, Gb, deg)
                if res: got = (box, res[0], Gb); break
            if not got: break
            out.append(got)
        if len(out) == len(plan): return out
    return None

# ---------------------------------------------------------------- tight leaves: exact r^2 = 2
def certify_exact(P, G, deg, eqdeg=None):
    """like certify, plus free columns (r^2 - 2) * monomial (both signs) for monomials of degree <= eqdeg."""
    eqdeg = deg if eqdeg is None else eqdeg
    base = certify_cols(P, G, deg)
    cols, labels, target = base
    eqcols, eqlabels = [], []
    for i in range(eqdeg + 1):
        for j in range(eqdeg + 1 - i):
            for k in range(2):
                for l in range(2):
                    m = {(i, j, k, l + 2): 1.0, (i, j, k, l): -2.0}
                    eqcols.append(m); eqlabels.append(('eq', (i, j, k, l)))
    allc = cols + eqcols
    monos = sorted({mm for c in allc for mm in c} | set(target))
    idx = {mm: n for n, mm in enumerate(monos)}
    A = np.zeros((len(monos), len(allc)))
    for jj, c in enumerate(allc):
        for mm, val in c.items(): A[idx[mm], jj] = val
    b = np.zeros(len(monos))
    for mm, val in target.items(): b[idx[mm]] = val
    nb = len(cols)
    bounds = [(0, None)] * nb + [(None, None)] * len(eqcols)
    cost = np.concatenate([np.ones(nb), np.zeros(len(eqcols))])
    res = linprog(cost, A_eq=A, b_eq=b, bounds=bounds, method='highs')
    if res.status != 0: return None
    used = [(labels[jj], res.x[jj]) for jj in range(nb) if res.x[jj] > 1e-10]
    eqused = [(eqlabels[jj - nb][1], res.x[jj]) for jj in range(nb, len(allc)) if abs(res.x[jj]) > 1e-10]
    return used, eqused, np.abs(A @ res.x - b).max()

def certify_cols(P, G, deg):
    uvG = [g for g in G if g[2] == 'uv']; HG = [g for g in G if g[2] == 'H']; rG = [g for g in G if g[2] == 'r']
    uvd = [todict(g[1]) for g in uvG]; uvdeg = [max(1, sp.Poly(g[1], u, v).total_degree()) for g in uvG]
    def exps(i, budget):
        if i == len(uvG):
            yield ()
            return
        for a in range(budget // uvdeg[i] + 1):
            for rest in exps(i + 1, budget - a * uvdeg[i]):
                yield (a,) + rest
    cols, labels = [], []
    for alpha in exps(0, deg):
        p = {(0, 0, 0, 0): 1.0}
        for gi, a in enumerate(alpha):
            for _ in range(a): p = dmul(p, uvd[gi])
        for hg in [None] + HG:
            ph = p if hg is None else dmul(p, todict(hg[1]))
            for rg in [None] + rG:
                pr = ph if rg is None else dmul(ph, todict(rg[1]))
                cols.append(pr); labels.append((alpha, hg[0] if hg else None, rg[0] if rg else None))
    return cols, labels, todict(P)

# ---------------------------------------------------------------- strip edge factors
EDGE = None
def strip_edges(Pexpr):
    """P = c * prod(f^m) * core with f among u, v, 1-u, 1-v, 1+u, 1+v, 1+u^2, 1+v^2 (all >= 0 on the box).
    Returns (factor_expr (>= 0 on the box), core) with P = factor_expr * core."""
    c, fl = sp.factor_list(sp.expand(Pexpr))
    nice = {u: u, v: v, u - 1: 1 - u, v - 1: 1 - v, u + 1: 1 + u, v + 1: 1 + v, u**2 + 1: 1 + u**2, v**2 + 1: 1 + v**2}
    fac, core = sp.Integer(1), sp.Integer(1) * c
    for f, m in fl:
        f = sp.expand(f)
        hit = None
        for key, pos in nice.items():
            if sp.expand(f - key) == 0: hit = pos; break
        if hit is not None:
            fac *= hit ** m
            if sp.expand(hit - f) != 0: core *= (-1) ** m
        else:
            core *= f ** m
    return fac, sp.expand(core)
