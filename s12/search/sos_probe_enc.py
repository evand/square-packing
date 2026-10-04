#!/usr/bin/env python3
"""sos_probe_enc.py -- polynomial encoding of a fixed separation type at n = 6, T = 3 (2026-10-03).

Task `tasks/sos-probe/README.md`, write-up `search/SOS_PROBE.md`.

Variables (in this order): T, then per square i: x_i, y_i, c_i, s_i  with
    c_i = cos(theta_i), s_i = sin(theta_i),  c_i^2 + s_i^2 = 1,  c_i >= |s_i|   (theta_i in [-45, 45 deg]:
    every orientation mod 90 once, and normal 0 / normal 1 of every square stay within 45 deg of x / y,
    so a typed separating direction is honestly 'near x' or 'near y').
(The brief's u_i = tan(theta_i/2) gives the same set; (c, s) keeps every constraint QUADRATIC, while
the u-form has degree 3 walls and degree 5 pair rows.  sos_probe_enc.u_form() checks the equivalence.)

Constraints of one separation type (all closed, so the feasible set is the closure of the packings of
that type; by scaling, s(6) = 3 is equivalent to  T >= 3  on every type's feasible set):
  walls   x_i - (c_i+-s_i)/2 >= 0,  T - x_i - (c_i+-s_i)/2 >= 0,  same in y   (bounding box, half-width (c+|s|)/2)
  pair    (i, j) separated by edge normal k of square o, sign sg  (n = sg * n_{o,k}):
            n.(c_q - c_o) - 1/2 - cosD/2 -+ sinD/2 >= 0        (two rows: |sin D| = max(+-sin D))
          q the other square, cosD = n_{o,0}.e_q >= 0 automatically (|D| <= 90 deg), sinD = n_{o,0}.f_q   (separating
          axis theorem: square q's half-extent along n is (|cos D| + |sin D|)/2).
A configuration is a list of 3x3 tiling cells (row, col); a type is a rule assigning (o, k, sg) to
each pair from the cells.  'far' pairs (cell distance 2 along an axis) are optional.
"""
import itertools
import math
from fractions import Fraction as Fr

import numpy as np

# ---------------------------------------------------------------- sparse polynomials (dict exp->Fr)


class Poly:
    __slots__ = ('nv', 't')

    def __init__(self, nv, terms=None):
        self.nv = nv
        self.t = {k: v for k, v in (terms or {}).items() if v != 0}

    @staticmethod
    def var(nv, i, c=1):
        e = [0] * nv; e[i] = 1
        return Poly(nv, {tuple(e): Fr(c)})

    @staticmethod
    def const(nv, c):
        return Poly(nv, {(0,) * nv: Fr(c)})

    def __add__(self, o):
        if not isinstance(o, Poly):
            o = Poly.const(self.nv, o)
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = r.get(k, 0) + v
        return Poly(self.nv, r)

    __radd__ = __add__

    def __neg__(self):
        return Poly(self.nv, {k: -v for k, v in self.t.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, Poly) else -Fr(o))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Poly):
            return Poly(self.nv, {k: v * Fr(o) for k, v in self.t.items()})
        r = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                k = tuple(a + b for a, b in zip(k1, k2))
                r[k] = r.get(k, 0) + v1 * v2
        return Poly(self.nv, r)

    __rmul__ = __mul__

    def deg(self):
        return max((sum(k) for k in self.t), default=0)

    def vars(self):
        return sorted({i for k in self.t for i, e in enumerate(k) if e})

    def subs(self, m):
        """substitute variables: m = {var index: Poly}."""
        out = Poly(self.nv)
        for k, v in self.t.items():
            term = Poly(self.nv, {tuple(0 if i in m else e for i, e in enumerate(k)): v})
            for i, e in enumerate(k):
                if e and i in m:
                    for _ in range(e):
                        term = term * m[i]
            out = out + term
        return out

    def ev(self, z):
        return sum(float(v) * np.prod([z[i] ** e for i, e in enumerate(k) if e]) for k, v in self.t.items())


# ---------------------------------------------------------------- configurations and types

ROW_A = [(0, 0), (0, 1), (0, 2), (1, 1), (2, 0), (2, 2)]     # row of three + the S6_LOCAL j=3 shape
ROW_B = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 1), (2, 2)]     # row of three + a staircase
ROW_C = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 1)]
PINWHEEL = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]   # 3x3 minus the main diagonal


def pair_type(ci, cj, i, j, diag='x', owner='low', far=True):
    """(o, k, sg) for cells ci=(r,c) of square i and cj of square j, or None (unconstrained).
    diag: axis used for diagonal neighbours ('x' or 'y'); owner: whose normal ('low' = the square
    with smaller coordinate along the axis, 'high' = the other, 'i' = square i)."""
    dr, dc = cj[0] - ci[0], cj[1] - ci[1]
    if max(abs(dr), abs(dc)) >= 2 and not far:
        return None
    if dr == 0:
        ax = 'x'
    elif dc == 0:
        ax = 'y'
    elif abs(dr) != abs(dc):
        ax = 'x' if abs(dc) > abs(dr) else 'y'
    else:
        ax = diag
    d = dc if ax == 'x' else dr                       # j minus i along the axis
    lo, hi = (i, j) if d > 0 else (j, i)
    o = {'low': lo, 'high': hi, 'i': i}[owner]
    k = 0 if ax == 'x' else 1
    # direction pointing from o towards the other square
    sg = 1 if o == lo else -1
    return (o, k, sg)


def make_type(cells, diag='x', owner='low', far=True, overrides=None):
    typ = {}
    for i, j in itertools.combinations(range(len(cells)), 2):
        t = pair_type(cells[i], cells[j], i, j, diag, owner, far)
        if overrides and (i, j) in overrides:
            t = overrides[(i, j)]
        if t is not None:
            typ[(i, j)] = t
    return typ


# ---------------------------------------------------------------- the polynomial system


def system(n, typ, Tfix=None, extra=True, orth=None):
    """returns (nv, ineqs, eqs, names).  Variable 0 is T (or the constant Tfix).
    orth: optional string of '+', '-', '0' per square fixing the sign of s_i; in a fixed orthant the
    half-width is (c_i + |s_i|)/2 with a known sign, so one wall row per side suffices, and the
    'extra' rows (c <= 1, |s| <= 1) are dropped (2(1-c) = (1-c)^2 + s^2 mod the ideal)."""
    nv = 1 + 4 * n
    V = lambda i: Poly.var(nv, i)
    T = V(0) if Tfix is None else Poly.const(nv, Tfix)
    X = [V(1 + 4 * i) for i in range(n)]
    Y = [V(2 + 4 * i) for i in range(n)]
    C = [V(3 + 4 * i) for i in range(n)]
    S = [V(4 + 4 * i) for i in range(n)]
    half = Fr(1, 2)
    G, H, names = [], [], []
    for i in range(n):
        o = orth[i] if orth else '0'
        sgns = {'+': (1,), '-': (-1,), '0': (1, -1)}[o]
        for sgn in sgns:                          # half-width p_i = (c_i + |s_i|)/2
            P = (C[i] + S[i] * sgn) * half
            for nm, g in (('xlo', X[i] - P), ('xhi', T - X[i] - P), ('ylo', Y[i] - P), ('yhi', T - Y[i] - P)):
                G.append(g); names.append((nm, i, sgn))
        if o == '0':
            G.append(C[i] - S[i]); names.append(('c>=s', i))
            G.append(C[i] + S[i]); names.append(('c>=-s', i))
            if extra:   # redundant but valid: c <= 1, |s| <= 1
                G.append(1 - C[i]); names.append(('c<=1', i))
                G.append(1 - S[i]); names.append(('s<=1', i))
                G.append(1 + S[i]); names.append(('s>=-1', i))
        else:
            sg = 1 if o == '+' else -1
            G.append(C[i] - S[i] * sg); names.append(('c>=|s|', i))
            G.append(S[i] * sg); names.append(('orth', i, o))
        H.append(C[i] * C[i] + S[i] * S[i] - 1)
    for (i, j), (o, k, sg) in sorted(typ.items()):
        q = j if o == i else i
        nx, ny = (C[o], S[o]) if k == 0 else (-S[o], C[o])
        nx, ny = nx * sg, ny * sg
        proj = nx * (X[q] - X[o]) + ny * (Y[q] - Y[o])
        cosD = C[o] * C[q] + S[o] * S[q]
        sinD = S[q] * C[o] - C[q] * S[o]
        for pm in (1, -1):
            G.append(proj - half - cosD * half - sinD * (half * pm))
            names.append(('pair', i, j, o, k, sg, pm))
    return nv, G, H, names


# ---------------------------------------------------------------- numerical sanity (local min of T)


def numeric_minT(cells, typ, rng, nstart=200, tilt=0.6, jit=0.15, fix_theta=None, tmin=None, orthant_plus=False):
    """multistart SLSQP: minimise T over the type's feasible set (angles theta, centres).
    fix_theta: dict i -> fixed angle.  Returns best T and the point."""
    from scipy.optimize import minimize
    n = len(cells)
    _nv, G, H, _nm = system(n, typ)
    # compile constraints to functions of z = (T, x, y, theta)
    def full(z):
        T = z[0]; w = [T]
        for i in range(n):
            th = z[1 + 2 * n + i] if fix_theta is None or i not in fix_theta else fix_theta[i]
            w += [z[1 + i], z[1 + n + i], math.cos(th), math.sin(th)]
        return w
    comp = [[(float(v), [(i, e) for i, e in enumerate(k) if e]) for k, v in g.t.items()] for g in G]

    def allg(z):
        w = full(z)
        return np.array([sum(v * math.prod(w[i] ** e for i, e in ie) for v, ie in cg) for cg in comp])
    cons = [{'type': 'ineq', 'fun': allg}]
    bounds = [(1.0, 4.0)] + [(0, 4)] * (2 * n) + [(0 if orthant_plus else -math.pi / 4, math.pi / 4)] * n
    best = (9, None)
    for _ in range(nstart):
        z0 = [3.0]
        z0 += [c[1] + 0.5 + rng.uniform(-jit, jit) for c in cells]
        z0 += [c[0] + 0.5 + rng.uniform(-jit, jit) for c in cells]
        z0 += list(rng.uniform(0 if orthant_plus else -tilt, tilt, n))
        r = minimize(lambda z: z[0], z0, method='SLSQP', constraints=cons, bounds=bounds,
                     options={'maxiter': 500, 'ftol': 1e-12})
        if r.success:
            viol = allg(r.x).min()
            if viol > -1e-8 and r.x[0] < best[0]:
                best = (r.x[0], r.x)
    return best


def u_form_check(rng, trials=200):
    """the (c,s) walls/pair rows agree with direct geometry (vertex projection) at random angles."""
    worst = 0.0
    for _ in range(trials):
        to, tq = rng.uniform(-math.pi / 4, math.pi / 4, 2)
        n = np.array([math.cos(to), math.sin(to)])
        verts = [np.array([math.cos(tq), math.sin(tq)]) * a / 2 + np.array([-math.sin(tq), math.cos(tq)]) * b / 2
                 for a in (1, -1) for b in (1, -1)]
        ext = max(abs(n @ v) for v in verts)
        cosD = math.cos(to) * math.cos(tq) + math.sin(to) * math.sin(tq)
        sinD = math.sin(tq) * math.cos(to) - math.cos(tq) * math.sin(to)
        worst = max(worst, abs(ext - (abs(cosD) + abs(sinD)) / 2), abs((math.cos(to) + abs(math.sin(to))) / 2
                    - max(abs(np.array([1, 0]) @ v) for v in [np.array([math.cos(to), math.sin(to)]) * a / 2
                    + np.array([-math.sin(to), math.cos(to)]) * b / 2 for a in (1, -1) for b in (1, -1)])))
    return worst
