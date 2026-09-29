#!/usr/bin/env python3
"""qx2_zm.py -- zero-margin exact checker for the qx2 box covers (task quadrant-exact-w2, 2026-09-28).

Statement checked: for every closed unit square Q in [0, m]^2 (any centre, any angle) mu(Q) >= 1, where mu is a
FORMAT.md v1 mixed cover made of axis-parallel uniform segments and ONE polygon: the square U = [a, m-a]^2 with
density exactly 1 (Lebesgue).  No points.  Proofs: search/QUADRANT_EXACT.md sec 4 (Lemmas U, K, E, Z).

This file IMPORTS search/zm_mixed.py (pinned, not edited; its sha256 and zeromargin.py's are printed) and runs its
MixedChecker box recursion with three new primitives tried first on every box:
  LEB   (Lemma U)  every admissible Q of the box lies in U: mu(Q) >= area(Q) = 1.
  CAP   (Lemma K)  Q pokes out of U through the lines y = a / x = a by a depth <= the minimum line density there.
  EXACT (Lemma E)  exact minimisation at fixed u by concavity on the cells of an arrangement of "bad" breaklines,
                   with every vertex a rational curve in u, and the certification of mass >= 1 along each curve by
                   exact univariate Bernstein positivity (zero margin allowed, e.g. at u -> 0).
and falls back to zm_mixed's own primitives (Lemmas S, T, L, R; polygon Lemma S(b) for U) otherwise.
The theta = 0 face is certified separately (qx2_exact.py axis, Lemma Z); EXACT certifies u in (u0, u1].
"""
import sys, os, time, math, argparse, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
from fractions import Fraction as F
from math import comb
import numpy as np
import zm_mixed as ZM
import zeromargin as zm
import mixed_cover as MC

HALF = F(1, 2)
# ============================================================================================== polynomials in u
# coefficient lists, low degree first, Fraction entries


def ptrim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p


def padd(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def psub(p, q):
    n = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0) for i in range(n)])


def pscale(p, k):
    if k == 0: return []
    return [c * k for c in p]


def pmul(p, q):
    if not p or not q: return []
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0: continue
        for j, b in enumerate(q):
            if b: r[i + j] += a * b
    return ptrim(r)


def ppow(p, k):
    r = [F(1)]
    for _ in range(k): r = pmul(r, p)
    return r


def peval(p, x):
    v = F(0)
    for c in reversed(p): v = v * x + c
    return v


def pevalf(p, x):
    v = 0.0
    for c in reversed(p): v = v * x + float(c)
    return v


def lowdeg(p):
    """index of the lowest nonzero coefficient (the order of the root at u = 0), None for the zero poly"""
    for i, c in enumerate(p):
        if c != 0: return i
    return None


def shift_scale(p, u0, h):
    """coefficients of p(u0 + h t) in t"""
    n = len(p)
    q = [F(0)] * n
    # Horner with (u0 + h t)
    for c in reversed(p):
        # q = q * (u0 + h t) + c
        nq = [F(0)] * n
        for i in range(n):
            if q[i]:
                nq[i] += q[i] * u0
                if i + 1 < n: nq[i + 1] += q[i] * h
        nq[0] += c
        q = nq
    return q


def bernstein(p, u0, u1):
    """Bernstein coefficients of p on [u0, u1] (degree len(p)-1)."""
    p = ptrim(p)
    if not p: return [F(0)]
    n = len(p) - 1
    c = shift_scale(p, u0, u1 - u0)
    return [sum((F(comb(i, k), comb(n, k)) * c[k] for k in range(i + 1)), F(0)) for i in range(n + 1)]


def nonneg(p, u0, u1, depth=6):
    """True only if p(u) >= 0 for every u in [u0, u1] (Bernstein convex-hull property, bisection)."""
    p = ptrim(p)
    if not p: return True
    b = bernstein(p, u0, u1)
    if all(x >= 0 for x in b): return True
    if b[0] < 0 or b[-1] < 0: return False          # b0 = p(u0), bn = p(u1)
    if depth == 0: return False
    m = (u0 + u1) / 2
    return nonneg(p, u0, m, depth - 1) and nonneg(p, m, u1, depth - 1)


def positive(p, u0, u1, depth=6):
    """True only if p(u) > 0 for every u in [u0, u1]."""
    p = ptrim(p)
    if not p: return False
    b = bernstein(p, u0, u1)
    if all(x > 0 for x in b): return True
    if b[0] <= 0 or b[-1] <= 0: return False
    if depth == 0: return False
    m = (u0 + u1) / 2
    return positive(p, u0, m, depth - 1) and positive(p, m, u1, depth - 1)


def sign_on(p, u0, u1, depth=6):
    """+1 / -1 if p is certified strictly positive / negative on [u0, u1], else 0."""
    if positive(p, u0, u1, depth): return 1
    if positive(pscale(p, -1), u0, u1, depth): return -1
    return 0


PN = [F(1), F(0), F(1)]         # N = 1 + u^2
PC = [F(1), F(0), F(-1)]        # C = 1 - u^2
PS = [F(0), F(2)]               # S = 2u
PN2 = [HALF, F(0), HALF]        # N/2
ONE = [F(1)]


# ============================================================================================== rational functions
# num / (C^eC S^eS N^eN D^eD) with D = the current candidate determinant (a polynomial set per evaluation context);
# C, N > 0 on [0,1); S > 0 on (0,1); the sign of D on the bin is certified before use.
class RF:
    __slots__ = ('num', 'e')

    def __init__(self, num, e=(0, 0, 0, 0)):
        self.num = ptrim(num); self.e = tuple(e)

    @staticmethod
    def const(v): return RF([F(v)])


def rf_lift(x, e, Dp):
    """rewrite x with the (larger) exponent vector e"""
    num = x.num
    for k, (a, b) in enumerate(zip(x.e, e)):
        if b < a: raise ValueError("rf_lift: exponent decrease")
        if b > a: num = pmul(num, ppow((PC, PS, PN, Dp)[k], b - a))
    return num


def rf_add(x, y, Dp, sy=1):
    e = tuple(max(a, b) for a, b in zip(x.e, y.e))
    nx = rf_lift(x, e, Dp); ny = rf_lift(y, e, Dp)
    return RF(padd(nx, pscale(ny, sy)), e)


def rf_scale(x, k):
    return RF(pscale(x.num, F(k)), x.e)


def rf_sum(xs, Dp):
    if not xs: return RF([])
    e = tuple(max(x.e[k] for x in xs) for k in range(4))
    num = []
    for x in xs: num = padd(num, rf_lift(x, e, Dp))
    return RF(num, e)


def rf_nonneg(x, u0, u1, sD, depth=6):
    """x >= 0 for all u in (u0, u1] (u0 >= 0), given sign(D) = sD on the bin.  Only the numerator's sign matters
    (the denominator is > 0 up to sD^eD); at u0 = 0 the denominator may vanish (S), which does not matter for
    u in (0, u1]; the numerator is tested on the closed [u0, u1] (sufficient by continuity)."""
    s = 1 if (x.e[3] % 2 == 0 or sD > 0) else -1
    return nonneg(pscale(x.num, s), u0, u1, depth)


def rf_positive(x, u0, u1, sD, depth=6):
    """x > 0 for all u in (u0, u1] -- via numerator > 0 on [u0,u1], or, when u0 = 0 and the numerator vanishes at
    0, numerator/u^k > 0 on [0, u1] (u^k > 0 on (0, u1])."""
    s = 1 if (x.e[3] % 2 == 0 or sD > 0) else -1
    p = pscale(x.num, s)
    if u0 == 0:
        k = lowdeg(p)
        if k is None: return False
        p = p[k:]
    return positive(p, u0, u1, depth)


def rf_float(x, u, Dp):
    v = pevalf(x.num, u)
    den = (1 - u * u) ** x.e[0] * (2 * u) ** x.e[1] * (1 + u * u) ** x.e[2] * pevalf(Dp, u) ** x.e[3]
    return v / den


# ============================================================================================== chord-end options
# Each option: t (the line parameter of the chord end) satisfies  F t = a0 + ax cx + ay cy,  F in {C, S}.
#   H line y = eta (t = x): k0 up C, k1 lo C, k2 lo S, k3 up S;  V line x = xi (t = y): k0 up S, k1 lo S, k2 up C, k3 lo C.
# (QUADRANT_EXACT.md sec 4.3; derived from |X|, |Y| <= 1/2 with X = ((p-c).(C,S))/N, Y = ((p-c).(-S,C))/N.)
def options(orient, p):
    if orient == 'H':
        eta = p
        return [('up', 'C', padd(PN2, pscale(PS, -eta)), PC, PS),
                ('lo', 'C', psub(pscale(PN2, -1), pscale(PS, eta)), PC, PS),
                ('lo', 'S', psub(pscale(PC, eta), PN2), PS, pscale(PC, -1)),
                ('up', 'S', padd(pscale(PC, eta), PN2), PS, pscale(PC, -1))]
    xi = p
    return [('up', 'S', psub(PN2, pscale(PC, xi)), PC, PS),
            ('lo', 'S', psub(pscale(PN2, -1), pscale(PC, xi)), PC, PS),
            ('up', 'C', padd(PN2, pscale(PS, xi)), pscale(PS, -1), PC),
            ('lo', 'C', padd(pscale(PN2, -1), pscale(PS, xi)), pscale(PS, -1), PC)]


def opt_float(opt, cx, cy, u):
    typ, Fk, a0, ax, ay = opt
    Fv = (1 - u * u) if Fk == 'C' else 2 * u
    return (pevalf(a0, u) + pevalf(ax, u) * cx + pevalf(ay, u) * cy) / Fv


# ============================================================================================== line density profiles
class Profile:
    """piecewise-constant density on a line: breakpoints b[0] < ... < b[n], density r[i] on (b[i], b[i+1]);
    0 outside.  G(t) = mass on (-inf, t].  Piece index i = 0..n: piece 0 = (-inf, b0] (G = 0), piece n+1 =
    [b_n, inf) (G = total), piece i+1 = [b_i, b_{i+1}]."""

    def __init__(self, segs):
        bs = sorted(set([s[0] for s in segs] + [s[1] for s in segs]))
        r = []
        for i in range(len(bs) - 1):
            lo, hi = bs[i], bs[i + 1]
            r.append(sum((dn for t0, t1, dn in segs if t0 <= lo and t1 >= hi), F(0)))
        self.b = bs; self.r = r
        g = [F(0)]
        for i in range(len(r)): g.append(g[-1] + r[i] * (bs[i + 1] - bs[i]))
        self.g = g                       # G(b[i])
        self.tot = g[-1]
        self.bf = [float(x) for x in bs]

    def pieces(self):
        """(lo, hi, rho, G(lo)) for all pieces including the two constant tails (lo/hi None = infinite)"""
        out = [(None, self.b[0], F(0), F(0))]
        for i in range(len(self.r)): out.append((self.b[i], self.b[i + 1], self.r[i], self.g[i]))
        out.append((self.b[-1], None, F(0), self.tot))
        return out

    def density_left_right(self, i):
        """densities left / right of breakpoint b[i]"""
        left = self.r[i - 1] if i >= 1 else F(0)
        right = self.r[i] if i < len(self.r) else F(0)
        return left, right

    def G(self, t):
        if t <= self.b[0]: return F(0)
        if t >= self.b[-1]: return self.tot
        import bisect
        i = bisect.bisect_right(self.b, t) - 1
        return self.g[i] + self.r[i] * (t - self.b[i])

    def Gf(self, t):
        import bisect
        if t <= self.bf[0]: return 0.0
        if t >= self.bf[-1]: return float(self.tot)
        i = bisect.bisect_right(self.bf, t) - 1
        return float(self.g[i]) + float(self.r[i]) * (t - self.bf[i])

    def piece_index_f(self, t):
        import bisect
        if t < self.bf[0]: return 0
        if t > self.bf[-1]: return len(self.r) + 1
        return min(bisect.bisect_right(self.bf, t), len(self.r)) if t >= self.bf[0] else 0


# ============================================================================================== c-lines, corners
# A c-line is (A, B, K): A(u) cx + B(u) cy + K(u) = 0, polynomial coefficients.
def cline_opt_eq_b(opt, b):
    typ, Fk, a0, ax, ay = opt
    Fp = PC if Fk == 'C' else PS
    return (ax, ay, psub(a0, pscale(Fp, b)))


def cline_opt_eq_opt(o1, o2):
    """t_1 = t_2 :  (a0_1 + ax_1 cx + ay_1 cy) F_2 - (a0_2 + ...) F_1 = 0"""
    F1 = PC if o1[1] == 'C' else PS
    F2 = PC if o2[1] == 'C' else PS
    return (psub(pmul(o1[3], F2), pmul(o2[3], F1)), psub(pmul(o1[4], F2), pmul(o2[4], F1)),
            psub(pmul(o1[2], F2), pmul(o2[2], F1)))


def coord_rep(kind, val):
    """a rectangle-edge coordinate as (num, den) polynomials: 'R' constant val, 'W' w/2 = (C+S)/(2N)"""
    if kind == 'R': return ([F(val)], ONE)
    return (padd(PC, PS), pscale(PN, 2))


def cline_edge(axis, kind, val):
    """the rectangle edge cx = coord (axis 0) or cy = coord (axis 1) as a c-line"""
    n, d = coord_rep(kind, val)
    if axis == 0: return (d, [], pscale(n, -1))
    return ([], d, pscale(n, -1))


def cline_at(cl, xr, yr):
    """numerator of the c-line's value at the corner (x, y) = (xn/xd, yn/yd), times xd*yd (> 0)"""
    A, B, K = cl
    (xn, xd), (yn, yd) = xr, yr
    return padd(padd(pmul(A, pmul(xn, yd)), pmul(B, pmul(yn, xd))), pmul(K, pmul(xd, yd)))


def trivial(cl):
    return not ptrim(cl[0]) and not ptrim(cl[1])


# ============================================================================================== EXACT (Lemma E)
class Exact:
    """Lemma E: exact certification of  mu(Q) >= tau  on a pose box with u in (u0, u1]  (see QUADRANT_EXACT.md 4.4)."""

    def __init__(self, cov, Ua, Ub, maxdepth_u=32, budget=96, verbose=False):
        if cov.points: raise ValueError("EXACT: covers with points are not supported")
        self.cov = cov; self.m = cov.m; self.Ua = Ua; self.Ub = Ub
        self.L = []
        for key, L in cov.lines.items():
            if key[0] not in ('H', 'V'): raise ValueError("EXACT: non-axis segment")
            prof = Profile(L['segs'])
            self.L.append(dict(key=key, o=key[0], p=key[1], prof=prof, opts=options(key[0], key[1]),
                               f=L['f']))
        self.maxdepth_u = maxdepth_u; self.budget = budget; self.verbose = verbose
        self.anchor = 'hi'
        self.stat = dict(cand=0, combos=0)

    # ------------------------------------------------------------------ geometry of the box
    def rect(self, box):
        """lower edges of the admissible centre rectangle on the bin: ('R', cx0) or ('W', w/2) per axis, or None
        when the wall crosses the box side within the bin (the driver then splits the bin)."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        if u1 >= F(2, 5): return None
        m = self.m
        c1, s1 = zm.trig(u1); c0, s0 = zm.trig(u0)
        w0, w1 = c0 + s0, c1 + s1                 # w increasing on [0, 45 deg]
        if cx1 + w1 / 2 > m or cy1 + w1 / 2 > m: return None
        out = []
        for lo in (cx0, cy0):
            if w0 / 2 >= lo: out.append(('W', None))
            elif w1 / 2 <= lo: out.append(('R', lo))
            else: return None
        return out

    def lines_in_reach(self, box):
        cx0, cx1, cy0, cy1 = (float(v) for v in box[:4])
        cxm, cym = (cx0 + cx1) / 2, (cy0 + cy1) / 2
        R = 0.7072 + 0.5 * math.hypot(cx1 - cx0, cy1 - cy0) + 1e-9
        out = []
        for L in self.L:
            px, py, dx, dy, tmin, tmax = L['f']
            if L['o'] == 'H':
                if abs(py - cym) > R: continue
                if float(tmax) < cxm - R or float(tmin) > cxm + R: continue
            else:
                if abs(px - cxm) > R: continue
                if float(tmax) < cym - R or float(tmin) > cym + R: continue
            out.append(L)
        return out

    def u_regime(self, box):
        """the Lebesgue square U = [a,b]^2: 'none' if Q ∩ U is empty at every pose of the box; otherwise a list of
        (axis, mode) for the lines x = a (axis 0) / y = a (axis 1) that Q may cross, with
            area(Q ∩ U) >= Phi_0 + Phi_1 - 1,   Phi_axis <= area(Q ∩ {axis-coordinate >= a})   (Lemma E', 4.4)
        mode 'std': Phi = 1 - ghat(d), ghat the cap area continued linearly beyond the third vertex (a concave
        minorant, exact while d <= cos theta); mode 'tan': d >= cos theta on the whole box, Phi = the tangent at a
        fixed depth d* of the convex branch (s + c - d)^2 / (2 s c) (an affine minorant).  Requires Q ⊂ {x, y <= b}."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        a, b = self.Ua, self.Ub
        c1, s1 = zm.trig(u1); c0, s0 = zm.trig(u0)
        w1 = c1 + s1; w0 = c0 + s0                        # w increasing on the bin (u1 < 2/5)
        if cx1 + w1 / 2 <= a or cy1 + w1 / 2 <= a: return 'none'
        if cx1 + w1 / 2 > b or cy1 + w1 / 2 > b: return None
        caps = []
        for axis, lo, hi in ((0, cx0, cx1), (1, cy0, cy1)):
            dmax = a - lo + w1 / 2
            if dmax <= 0: continue
            dmin = a - hi + w0 / 2
            if dmin >= c0:
                dstar = min(F(1), (dmin + dmax) / 2)
                caps.append((axis, 'tan', dstar))
            else:
                caps.append((axis, 'std', None))
        if len(caps) == 2 and all(m_ == 'std' for _, m_, _ in caps) and self.corner_ok(box):
            caps.append((2, 'corner', self.anchor))
        return caps

    def corner_ok(self, box):
        """the lower-left configuration of Lemma E'' on the whole box: Q ∩ {x < a, y < a} is the quadrilateral
        (BL, bottom edge ∩ {x = a}, (a, a), left edge ∩ {y = a}):  x_BR >= a, y_TL >= a, (a, a) in Q, alpha, beta >= 0,
        beta C >= alpha S  (alpha = (a - x_BL) N, beta = (a - y_BL) N)."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        a = self.Ua
        c1, s1 = zm.trig(u1); c0, s0 = zm.trig(u0)
        if cx0 + (c0 + s0) / 2 < a: return False                 # x_BR = cx + (c + s)/2 >= a
        if cy0 + (c1 - s1) / 2 < a: return False                 # y_TL = cy + (c - s)/2 >= a
        CmS = psub(PC, PS); CpS = padd(PC, PS)
        al_hi = padd(pscale(PN, a - cx0), pscale(CmS, HALF)); al_lo = padd(pscale(PN, a - cx1), pscale(CmS, HALF))
        be_hi = padd(pscale(PN, a - cy0), pscale(CpS, HALF)); be_lo = padd(pscale(PN, a - cy1), pscale(CpS, HALF))
        checks = [al_lo, be_lo,
                  psub(pmul(be_lo, PC), pmul(al_hi, PS)),                               # beta C - alpha S >= 0
                  psub(PN2, padd(pscale(PC, a - cx0), pscale(PS, a - cy0))),              # X(a,a) <= 1/2
                  psub(PN2, psub(pscale(PC, a - cy0), pscale(PS, a - cx1)))]              # Y(a,a) <= 1/2
        return all(nonneg(p, u0, u1, 4) for p in checks)

    def corner_term(self, box, Dp, X, Y, anchor='hi'):
        """Lemma E'': a concave minorant of area(Q ∩ {x < a, y < a}) = (2 C alpha beta + S (beta^2 - alpha^2)) / (2 C N^2),
        with alpha beta >= alpha_lo beta + beta_lo alpha - alpha_lo beta_lo (McCormick) and S beta^2 >= S (2 b* beta - b*^2)."""
        cx0, cx1, cy0, cy1 = box[:4]
        a = self.Ua
        CmS = psub(PC, PS); CpS = padd(PC, PS)
        # alpha = (a - cx) N + (C - S)/2 = an / D,  beta = (a - cy) N + (C + S)/2 = bn / D
        an = padd(pmul(psub(pscale(Dp, a), X), PN), pscale(pmul(CmS, Dp), HALF))
        bn = padd(pmul(psub(pscale(Dp, a), Y), PN), pscale(pmul(CpS, Dp), HALF))
        # McCormick anchor: 'lo' uses (alpha - alpha_lo)(beta - beta_lo) >= 0 (exact at the box corner (cx1, cy1)),
        # 'hi' uses (alpha_hi - alpha)(beta_hi - beta) >= 0 (exact at (cx0, cy0)); both are the same formula with the
        # anchor values substituted: alpha beta >= A beta + B alpha - A B.
        if anchor == 'lo':
            al_lo = padd(pscale(PN, a - cx1), pscale(CmS, HALF)); be_lo = padd(pscale(PN, a - cy1), pscale(CpS, HALF))
        else:
            al_lo = padd(pscale(PN, a - cx0), pscale(CmS, HALF)); be_lo = padd(pscale(PN, a - cy0), pscale(CpS, HALF))
        bstar = padd(pscale(PN, a - (cy0 + cy1) / 2), pscale(CpS, HALF))
        D2 = pmul(Dp, Dp)
        t1 = pscale(pmul(PC, padd(padd(pmul(pmul(al_lo, bn), Dp), pmul(pmul(be_lo, an), Dp)), pscale(pmul(pmul(al_lo, be_lo), D2), -1))), 2)
        t2 = pmul(PS, psub(pscale(pmul(pmul(bstar, bn), Dp), 2), pmul(pmul(bstar, bstar), D2)))
        t3 = pscale(pmul(PS, pmul(an, an)), -1)
        num = pscale(padd(padd(t1, t2), t3), HALF)
        return RF(num, (1, 0, 2, 2))

    def cap_lines(self, axis):
        """c-lines d = 0 (Q touches the line) and d = sin(theta) (the second vertex crosses it)"""
        a = self.Ua
        n0 = padd(pscale(PN, 2 * a), padd(PC, PS)); n1 = padd(pscale(PN, 2 * a), psub(PC, PS))
        out = []
        for n in (n0, n1):
            if axis == 0: out.append((pscale(PN, 2), [], pscale(n, -1)))
            else: out.append(([], pscale(PN, 2), pscale(n, -1)))
        return out

    def phi_term(self, cap, Dp, X, Y, u0, u1, sD):
        """Phi_axis along the candidate curve: an RF lower bound of area(Q ∩ {coordinate >= a})."""
        axis, mode, dstar = cap
        if mode == 'std':
            return rf_add(RF.const(1), self.cap_term(axis, Dp, X, Y, u0, u1, sD), Dp, -1)
        a = self.Ua
        Z = X if axis == 0 else Y
        dn = pscale(padd(pmul(pscale(PN, 2), psub(pscale(Dp, a), Z)), pmul(padd(PC, PS), Dp)), HALF)   # d = dn/(N D)
        k = psub(padd(PS, PC), pscale(PN, dstar))                          # (s + c - d*) N
        # T = [k^2 D / 2 - k (dn - d* N D)] / (S C D)
        num = psub(pscale(pmul(pmul(k, k), Dp), HALF), pmul(k, psub(dn, pscale(pmul(PN, Dp), dstar))))
        return RF(num, (1, 1, 0, 1))

    def cap_term(self, axis, Dp, X, Y, u0, u1, sD):
        """an RF upper bound of ghat(d) along the candidate curve (ghat: 0 for d <= 0, d^2/(2sc) on [0, s],
        (d - s/2)/c for d >= s; ghat >= the true cap area for every d)."""
        a = self.Ua
        Z = X if axis == 0 else Y
        # d = a - z + (C+S)/(2N) = dn / (N D),  dn = (2N (a D - Z) + (C + S) D) / 2
        dn = pscale(padd(pmul(pscale(PN, 2), psub(pscale(Dp, a), Z)), pmul(padd(PC, PS), Dp)), HALF)
        d = RF(dn, (0, 0, 1, 1))
        if rf_nonneg(rf_scale(d, -1), u0, u1, sD, 4):                    # d <= 0: no cap
            return RF([])
        s_ = RF(PS, (0, 0, 1, 0))                                        # sin theta = S / N
        if rf_nonneg(rf_add(d, s_, Dp, -1), u0, u1, sD, 4):              # d >= sin: linear regime
            return RF(psub(dn, pscale(pmul(PS, Dp), HALF)), (1, 0, 0, 1))   # (d N - S/2)/C
        return RF(pscale(pmul(dn, dn), HALF), (1, 1, 0, 2))               # d^2 N^2/(2 S C) >= g for all d <= cos

    # ------------------------------------------------------------------ main entry
    def certify(self, box, tau=F(1), lam=F(0)):
        """True only if mu(Q) >= tau at every admissible pose of the box with u in (u0, u1] (u in [u0,u1] if u0>0).
        In the corner regime (Lemma E'') both McCormick anchors are tried."""
        self.anchor = 'hi'
        ok, why = self._certify(box, tau, lam)
        if not ok and self._had_corner:
            self.anchor = 'lo'
            ok, why = self._certify(box, tau, lam)
            self.anchor = 'hi'
        return ok, why

    def _certify(self, box, tau, lam):
        self._had_corner = False
        rect = self.rect(box)
        if rect is None: return False, 'rect'
        reg = self.u_regime(box)
        if reg is None: return False, 'U-regime'
        self._had_corner = reg != 'none' and any(cp[1] == 'corner' for cp in reg)
        cx0, cx1, cy0, cy1, u0, u1 = box
        xlo = coord_rep(rect[0][0], cx0); xhi = coord_rep('R', cx1)
        ylo = coord_rep(rect[1][0], cy0); yhi = coord_rep('R', cy1)
        corners = [(xr, yr) for xr in (xlo, xhi) for yr in (ylo, yhi)]
        edges = [cline_edge(0, rect[0][0], cx0), cline_edge(0, 'R', cx1),
                 cline_edge(1, rect[1][0], cy0), cline_edge(1, 'R', cy1)]
        lines = self.lines_in_reach(box)

        def crosses(cl):
            """False only if the c-line is certified to miss R(u) for every u of the bin"""
            if trivial(cl): return False
            sg = [sign_on(cline_at(cl, xr, yr), u0, u1, 4) for xr, yr in corners]
            if all(s == 1 for s in sg) or all(s == -1 for s in sg): return False
            return True
        # ---- breaklines: bad kinks of each option, chord-empty pairs
        brk = []
        for L in lines:
            prof = L['prof']; opts = L['opts']
            for i, (typ, Fk, a0, ax, ay) in enumerate(opts):
                for bi, b in enumerate(prof.b):
                    left, right = prof.density_left_right(bi)
                    bad = (right > left) if typ == 'up' else (right < left)
                    if not bad: continue
                    cl = cline_opt_eq_b(opts[i], b)
                    if crosses(cl): brk.append(cl)
            ups = [o for o in opts if o[0] == 'up']; los = [o for o in opts if o[0] == 'lo']
            for ou in ups:
                for ol in los:
                    cl = cline_opt_eq_opt(ou, ol)
                    if crosses(cl): brk.append(cl)
        if reg != 'none':
            for axis, mode, _ in reg:
                if mode != 'std' or axis == 2: continue
                for cl in self.cap_lines(axis):
                    if crosses(cl): brk.append(cl)
        allc = edges + brk
        self.nbrk = len(brk)
        # ---- candidates
        for i in range(len(allc)):
            for j in range(i + 1, len(allc)):
                if i < 4 and j < 4 and (i // 2 == j // 2): continue      # parallel rectangle edges
                ok, why = self._candidate(allc[i], allc[j], box, rect, lines, reg, (tau, lam), corners)
                if not ok:
                    return False, why
        return True, f"EXACT[{len(brk)} brk]"

    # ------------------------------------------------------------------ one candidate vertex curve
    def _candidate(self, l1, l2, box, rect, lines, reg, tau, corners):
        A1, B1, K1 = l1; A2, B2, K2 = l2
        Dp = psub(pmul(A1, B2), pmul(A2, B1))
        if not Dp: return True, 'parallel'
        X = psub(pmul(B1, K2), pmul(B2, K1)); Y = psub(pmul(A2, K1), pmul(A1, K2))
        u0, u1 = box[4], box[5]
        self._left = self.budget                   # sub-bin evaluations allowed for this candidate
        return self._cand_bin(Dp, X, Y, box, rect, lines, reg, tau, u0, u1, 0)

    def _cand_bin(self, Dp, X, Y, box, rect, lines, reg, tau, u0, u1, depth):
        cx0, cx1, cy0, cy1 = box[:4]
        # sign of D on (u0, u1]
        k = lowdeg(Dp) if u0 == 0 else 0
        sD = sign_on(Dp[k:], u0, u1, 5)
        if sD == 0:
            return self._split(Dp, X, Y, box, rect, lines, reg, tau, u0, u1, depth, 'D sign')
        cxr = RF(X, (0, 0, 0, 1)); cyr = RF(Y, (0, 0, 0, 1))
        # outside the admissible rectangle for every u of the sub-bin?  (strict, certified)
        cons = []
        for axis, (kind, val), hi in ((0, rect[0], cx1), (1, rect[1], cy1)):
            lo_val = box[0] if axis == 0 else box[2]
            cc = cxr if axis == 0 else cyr
            if kind == 'R': lo = RF(psub(cc.num, pscale(Dp, lo_val)), (0, 0, 0, 1))
            else: lo = RF(psub(pmul(pscale(PN, 2), cc.num), pmul(padd(PC, PS), Dp)), (0, 0, 1, 1))
            hi_ = RF(psub(pscale(Dp, hi), cc.num), (0, 0, 0, 1))
            cons += [lo, hi_]
        for c in cons:
            if rf_positive(rf_scale(c, -1), u0, u1, sD, 4):
                return True, 'outside'
        # float guidance at a few u's
        us = [float(u0) + (float(u1) - float(u0)) * t for t in (0.02, 0.5, 1.0)]
        if float(u0) == 0: us[0] = float(u1) * 1e-3
        pos = []
        for u in us:
            dv = pevalf(Dp, u)
            pos.append((pevalf(X, u) / dv, pevalf(Y, u) / dv, u))
        terms = []           # list of lists of RF alternatives (a combo picks one per list; the sum must be >= 0)
        const = []           # RF terms added to every combo
        tau0, lam = tau                      # the claim checked: mass >= tau0 + lam * u
        const.append(RF([-F(tau0), -F(lam)]))
        if reg != 'none':
            # area(Q ∩ U) >= sum_axis Phi_axis - (#caps - 1) [+ area(Q ∩ {x<a, y<a}) when the corner regime holds]
            lcaps = [cp for cp in reg if cp[1] != 'corner']
            const.append(RF.const(1 - len(lcaps)))
            for cap in lcaps:
                const.append(self.phi_term(cap, Dp, X, Y, u0, u1, sD))
            if len(lcaps) < len(reg):
                const.append(self.corner_term(box, Dp, X, Y, [cp for cp in reg if cp[1] == 'corner'][0][2]))
        for L in lines:
            r = self._line_terms(L, Dp, X, Y, pos, u0, u1, sD)
            if r is None: continue
            ualts, lalts = r
            ualts = self._prune_dom(ualts, u0, u1, sD, Dp, pos)
            nl = self._prune_dom([rf_scale(x, -1) for x in lalts], u0, u1, sD, Dp, pos)
            if len(ualts) == 1: const.append(ualts[0])
            else: terms.append(ualts)
            if len(nl) == 1: const.append(nl[0])
            else: terms.append(nl)
        base = rf_sum(const, Dp)
        # combos: the bound is the sum over lines of a minimum over alternatives = the minimum over combos
        ncomb = 1
        for alts in terms: ncomb *= len(alts)
        if ncomb > 64:
            return self._split(Dp, X, Y, box, rect, lines, reg, tau, u0, u1, depth, 'combos')
        combos = [base]
        for alts in terms:
            combos = [rf_add(c, a, Dp) for c in combos for a in alts]
        self.stat['cand'] += 1; self.stat['combos'] += len(combos)
        # constraints of R(u) that the candidate may violate on this sub-bin (for the S-procedure below)
        open_cons = [g for g in cons if not rf_nonneg(g, u0, u1, sD, 3)]
        for c in combos:
            if not rf_nonneg(c, u0, u1, sD, 5) and not self._sproc(c, open_cons, u0, u1, sD, Dp):
                if depth >= self.maxdepth_u and self.verbose:
                    um = float(u0 + u1) / 2
                    self.lastfail = dict(u=(float(u0), float(u1)), pos=[(pevalf(X, uu) / pevalf(Dp, uu), pevalf(Y, uu) / pevalf(Dp, uu), uu)
                                                                       for uu in (float(u0) + 1e-12, um, float(u1))],
                                         combo=[rf_float(cc, um, Dp) for cc in combos], nterms=len(terms))
                return self._split(Dp, X, Y, box, rect, lines, reg, tau, u0, u1, depth, 'value')
        return True, 'ok'

    @staticmethod
    def _sproc(c, gs, u0, u1, sD, Dp):
        """S-procedure: c - lam * g >= 0 on the sub-bin with lam >= 0 constant implies c >= 0 wherever g >= 0,
        i.e. wherever the candidate is inside that side of R(u) (outside it the candidate is no vertex of R(u))."""
        us = [float(u0) + (float(u1) - float(u0)) * t for t in (0.0, 0.25, 0.5, 0.75, 1.0)]
        if float(u0) == 0: us[0] = float(u1) * 1e-6
        for g in gs:
            # a lam that makes c - lam g >= 0 at the float sample points, if any
            lam = 0.0; okf = True
            for u in us:
                cv = rf_float(c, u, Dp); gv = rf_float(g, u, Dp)
                if cv >= 0: continue
                if gv >= 0: okf = False; break
                lam = max(lam, cv / gv)
            if not okf: continue
            for f_ in (1.01, 1.5, 3, 10):
                L = F(lam * f_ + 1e-12).limit_denominator(10 ** 9)
                if rf_nonneg(rf_add(c, rf_scale(g, L), Dp, -1), u0, u1, sD, 4): return True
        return False

    def _split(self, Dp, X, Y, box, rect, lines, reg, tau, u0, u1, depth, why):
        self._left -= 1
        if depth >= self.maxdepth_u or self._left <= 0:
            return False, why
        m = (u0 + u1) / 2
        ok, w = self._cand_bin(Dp, X, Y, box, rect, lines, reg, tau, u0, m, depth + 1)
        if not ok: return ok, w
        return self._cand_bin(Dp, X, Y, box, rect, lines, reg, tau, m, u1, depth + 1)

    # ------------------------------------------------------------------ one line at the candidate
    def _line_terms(self, L, Dp, X, Y, pos, u0, u1, sD):
        """lower bound of the line's mass along the candidate curve on the sub-bin, as (up alternatives,
        lo alternatives): mass >= min(up alts) - max(lo alts); or None (line dropped: bound 0).
        The up alternatives cover every (option, piece) that the up end can use on the sub-bin (options certified
        dominated are dropped), and likewise for lo."""
        prof = L['prof']; opts = L['opts']
        # float chord at the sample positions
        fl = []
        nonempty = False
        for (cx, cy, u) in pos:
            vals = [opt_float(o, cx, cy, u) for o in opts]
            up = min(v for o, v in zip(opts, vals) if o[0] == 'up')
            lo = max(v for o, v in zip(opts, vals) if o[0] == 'lo')
            if prof.Gf(up) - prof.Gf(lo) > 1e-12: nonempty = True
            fl.append(vals)
        if not nonempty: return None
        # the options' exact values along the curve
        tv = []
        for (typ, Fk, a0, ax, ay) in opts:
            num = padd(padd(pmul(a0, Dp), pmul(ax, X)), pmul(ay, Y))
            tv.append(RF(num, (1, 0, 0, 1) if Fk == 'C' else (0, 1, 0, 1)))
        out = []
        for side in ('up', 'lo'):
            idx = [i for i, o in enumerate(opts) if o[0] == side]
            # float-best option per sample; keep the ones that are best somewhere, certify the others dominated
            best = set()
            for vals in fl:
                vs = [(vals[i], i) for i in idx]
                best.add(min(vs)[1] if side == 'up' else max(vs)[1])
            keep = set(best)
            for i in idx:
                if i in keep: continue
                dom = False
                for j in keep:
                    d = rf_add(tv[i], tv[j], Dp, -1)          # t_i - t_j
                    if side == 'lo': d = rf_scale(d, -1)      # lo: need t_i <= t_j
                    if rf_nonneg(d, u0, u1, sD, 4): dom = True; break
                if not dom: keep.add(i)
            alts = []
            for i in sorted(keep):
                for pc in self._pieces(prof, tv[i], [vals[i] for vals in fl], Dp, u0, u1, sD):
                    lo_, hi_, rho, g = pc
                    # A_p(t) = g + rho (t - lo_)   (tails: constant)
                    if rho == 0: alts.append(RF.const(g))
                    else:
                        base_ = g - rho * lo_
                        alts.append(rf_add(RF.const(base_), rf_scale(tv[i], rho), Dp))
            out.append(self._prune_const(alts, side))
        return out[0], out[1]

    @staticmethod
    def _prune_dom(alts, u0, u1, sD, Dp, pos):
        """drop alternatives certified >= another kept alternative on the sub-bin (the minimum is unchanged)."""
        if len(alts) <= 1: return alts
        us = [p[2] for p in pos]
        order = sorted(range(len(alts)), key=lambda i: min(rf_float(alts[i], u, Dp) for u in us))
        keep = []
        for i in order:
            if any(rf_nonneg(rf_add(alts[i], alts[j], Dp, -1), u0, u1, sD, 3) for j in keep): continue
            keep.append(i)
        return [alts[i] for i in keep]

    @staticmethod
    def _prune_const(alts, side):
        """identical alternatives collapse; among constants only the extreme one matters"""
        seen = []; consts = []
        for a in alts:
            if len(a.num) <= 1 and a.e == (0, 0, 0, 0): consts.append(a.num[0] if a.num else F(0)); continue
            if any(a.num == b.num and a.e == b.e for b in seen): continue
            seen.append(a)
        if consts: seen.append(RF.const(min(consts) if side == 'up' else max(consts)))
        return seen

    def _pieces(self, prof, t, tfloats, Dp, u0, u1, sD):
        """pieces of the profile that the value t(u) may occupy on the sub-bin: the float range widened until the
        bounding breakpoints are certified (t >= b_lo and t <= b_hi on the sub-bin)."""
        P = prof.pieces()
        tmin, tmax = min(tfloats), max(tfloats)
        lo_i = prof.piece_index_f(tmin); hi_i = prof.piece_index_f(tmax)
        nb = len(prof.b)
        # certify t >= b[lo_i - 1] (lower boundary of piece lo_i)  -- piece i (1..n) is [b[i-1], b[i]]
        while lo_i > 0:
            b = prof.b[lo_i - 1]
            if rf_nonneg(rf_add(t, RF.const(b), Dp, -1), u0, u1, sD, 4): break
            lo_i -= 1
        while hi_i < nb:
            b = prof.b[hi_i]
            if rf_nonneg(rf_add(RF.const(b), t, Dp, -1), u0, u1, sD, 4): break
            hi_i += 1
        return P[lo_i:hi_i + 1]


# ============================================================================================== LEB (Lemma U) and CAP (Lemma K)
def u_square(cov):
    """the cover's single polygon must be an axis-parallel square [a, b]^2 with density exactly 1."""
    if len(cov.polys) != 1: raise ValueError("expected exactly one polygon (the Lebesgue square)")
    P = cov.polys[0]
    xs = sorted(set(v[0] for v in P['V'])); ys = sorted(set(v[1] for v in P['V']))
    if len(P['V']) != 4 or xs != ys or len(xs) != 2: raise ValueError("the polygon is not a square [a,b]^2")
    a, b = xs
    if P['w'] != (b - a) ** 2 or P['area'] != (b - a) ** 2: raise ValueError("the polygon's density is not 1")
    if not (0 < a < b <= cov.m): raise ValueError("bad square")
    return a, b


def line_min_density(cov, key, t0, t1):
    """min over [t0, t1] of the (summed) density of the cover's segments on the line `key`; 0 if not covered."""
    L = cov.lines.get(key)
    if L is None: return F(0)
    prof = Profile(L['segs'])
    if t0 < prof.b[0] or t1 > prof.b[-1]: return F(0)
    best = None
    for i in range(len(prof.r)):
        lo, hi = prof.b[i], prof.b[i + 1]
        if hi <= t0 or lo >= t1: continue
        best = prof.r[i] if best is None else min(best, prof.r[i])
    return best if best is not None else F(0)


def cert_leb(box, B, a, b):
    """Lemma U: every admissible Q of the box lies in [a, b]^2 (where mu >= Lebesgue): mu(Q) >= 1."""
    cx0, cx1, cy0, cy1 = box[:4]
    h = B['whi'] / 2
    return cx0 - h >= a and cy0 - h >= a and cx1 + h <= b and cy1 + h <= b


def cert_cap(box, B, cov, a, b):
    """Lemma K (QUADRANT_EXACT.md 4.2): Q may leave U = [a,b]^2 only through the lines y = a and x = a, by depths
    d_y, d_x <= the minimum density of the cover's segments on those lines over the reachable chord range, and the
    centre is inside the half-planes (cy >= a when Q crosses y = a; cx >= a when it crosses x = a).  Then
    mu(Q) >= area(Q ∩ U) + mass(y = a) + mass(x = a) >= 1 - c_y d_y - c_x d_x + rho_y c_y + rho_x c_x >= 1."""
    cx0, cx1, cy0, cy1 = box[:4]
    h = B['whi'] / 2
    if cx1 + h > b or cy1 + h > b: return False
    ok_any = False
    for axis in (0, 1):
        lo0 = cy0 if axis == 0 else cx0          # the coordinate across the line (y for the line y = a)
        d = a - lo0 + h                          # max depth below the line over the box
        if d <= 0: continue                      # never crosses this line
        if lo0 < a: return False                 # the width condition needs the centre on U's side
        o0, o1 = (cx0, cx1) if axis == 0 else (cy0, cy1)
        key = ('H', a) if axis == 0 else ('V', a)
        rho = line_min_density(cov, key, o0 - h, o1 + h)
        if not (d <= rho): return False
        ok_any = True
    return True


# ============================================================================================== the checker
class QXChecker(ZM.MixedChecker):
    """zm_mixed.MixedChecker with LEB / CAP tried first on every box and EXACT tried before a box is split
    (when u1 <= exact_umax).  Everything else is zm_mixed's (imported, unedited)."""

    def __init__(self, cover, exact_umax=F(1, 20), exact_from=0, use_exact=True, **kw):
        super().__init__(cover, **kw)
        self.Ua, self.Ub = u_square(cover)
        self.exact = Exact(cover, self.Ua, self.Ub) if use_exact else None
        self.exact_umax = exact_umax; self.exact_from = exact_from

    def run_box(self, root):
        zc = self.zc
        stats = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'PIECE': 0, 'EMPTY': 0,
                 'UNCERT': 0, 'boxes': 0, 'maxdepth': 0, 'THR': 0, 'LIN': 0, 'TPTS': 0, 'SPLIT': 0, 'cpu': 0.0,
                 'LEB': 0, 'CAP': 0, 'EXACT': 0, 'EXACT0': 0, 'AXIS': 0}
        unc = []
        _t0 = time.process_time(); leaves = []
        n = len(zc.P)
        base = np.zeros(n, dtype=bool); base[self.ph] = True
        stack = [(root, 0, None, F(0))]
        bincache = {}
        while stack:
            box, depth, inh, Lpar = stack.pop()
            stats['boxes'] += 1
            stats['maxdepth'] = max(stats['maxdepth'], depth)
            cx0, cx1, cy0, cy1, u0, u1 = box
            if self.clip and u1 > u0:
                cu1 = zm.clip_bin(box, self.m)
                if cu1 < u1:
                    u1 = cu1; box = (cx0, cx1, cy0, cy1, u0, u1)
            key = (u0, u1)
            if key not in bincache:
                B = zm.bin_data(u0, u1)
                bincache[key] = (B, None)
            B = bincache[key][0]
            lo = B['wlo'] / 2
            if cx1 < lo or cx0 > self.m - lo or cy1 < lo or cy0 > self.m - lo:
                stats['EMPTY'] += 1
                if self.dump: leaves.append((box, 'EMPTY', None))
                continue
            kind = wit = None
            # ---- new primitives
            if u1 == 0:
                kind, wit = 'AXIS', None          # only theta = 0 poses left (after clip_bin): Lemma Z
            elif cert_leb(box, B, self.Ua, self.Ub):
                kind, wit = 'LEB', None
            elif cert_cap(box, B, self.cov, self.Ua, self.Ub):
                kind, wit = 'CAP', None
            thr = False
            L = F(0)
            if kind is None:
                # ---- zm_mixed's primitives (verbatim logic of MixedChecker.run_box; no points in these covers,
                #      so the point primitives only see the phantom)
                lparts = None
                if self.use_pieces:
                    L, thr = self.piece_bound(box, B)
                    lparts = self._lparts
                    if Lpar > L: L = Lpar
                if L >= 1:
                    kind, wit = 'PIECE', str(L)
                else:
                    Lr = self._set_phantom(L)
                    inh2 = base.copy() if inh is None else (inh | base)
                    zc._last_inT = None
                    kind, wit = zc.cert_adm(box, B, inh2)
                    if kind is None: kind, wit = zc.cert_p1(box, B)
                    if kind is None: kind, wit = zc.cert_mix(box, B, inh2)
                    if kind is None and self.use_chain and depth >= self.chain_from:
                        kind, wit = zc.cert_chain(box, B, inh=inh2)
                    if kind is None and self.use_split and lparts is not None:
                        kind, wit = self.cert_split(box, B, inh, L, lparts)
            if kind is None and self.exact is not None and u1 <= self.exact_umax and depth >= self.exact_from \
                    and u1 > u0:
                ok, why = self.exact.certify(box)
                if ok:
                    kind = 'EXACT0' if u0 == 0 else 'EXACT'; wit = why
            if kind is not None:
                stats[kind] += 1
                if thr == 'T': stats['THR'] += 1
                if thr == 'L': stats['LIN'] += 1
                if self.dump: leaves.append((box, kind, wit))
                continue
            if depth >= self.max_depth:
                stats['UNCERT'] += 1; unc.append(box); continue
            dx, dy, du = cx1 - cx0, cy1 - cy0, 2 * (u1 - u0)
            if self.theta_bias > 1 and u0 * 2 < B['wlo'] and (
                    cx0 < B['whi'] / 2 or cx1 > self.m - B['whi'] / 2 or
                    cy0 < B['whi'] / 2 or cy1 > self.m - B['whi'] / 2):
                du = du * self.theta_bias
            kid = zc._last_inT if zc._last_inT is not None else inh
            if kid is not None: kid = kid.copy(); kid[self.ph] = False
            if dx >= dy and dx >= du:
                mid = (cx0 + cx1) / 2
                stack.append(((cx0, mid, cy0, cy1, u0, u1), depth + 1, kid, L))
                stack.append(((mid, cx1, cy0, cy1, u0, u1), depth + 1, kid, L))
            elif dy >= du:
                mid = (cy0 + cy1) / 2
                stack.append(((cx0, cx1, cy0, mid, u0, u1), depth + 1, kid, L))
                stack.append(((cx0, cx1, mid, cy1, u0, u1), depth + 1, kid, L))
            else:
                mid = (u0 + u1) / 2
                stack.append(((cx0, cx1, cy0, cy1, u0, mid), depth + 1, kid, L))
                stack.append(((cx0, cx1, cy0, cy1, mid, u1), depth + 1, kid, L))
        stats['cpu'] = time.process_time() - _t0
        return stats, unc, leaves


# ============================================================================================== Lemma Z on the box cover
def axis_face(cov, lo=HALF, hi=None, verbose=True):
    """Lemma Z (QUADRANT_EXACT.md 4.1) on the box cover itself: at theta = 0 the mass of the closed square
    [cx-1/2, cx+1/2] x [cy-1/2, cy+1/2] is, on each open cell of the grid of breakpoints (every segment end / density
    change / line position, and the Lebesgue square's sides, each +-1/2), a multilinear function of (cx, cy); its
    values at points of the grid lines are >= the adjacent one-sided limits (upper semicontinuity: closed squares).
    So inf over the theta = 0 poses with centre in [lo, hi]^2 = min over the one-sided limits at the grid corners,
    computed here exactly.  Returns (min, argmin, count, tight count)."""
    m = cov.m
    if hi is None: hi = m / 2
    a, b = u_square(cov)
    xs = set([a, b]); prof = {}
    for key, L in cov.lines.items():
        P = Profile(L['segs']); prof[key] = P
        xs.update(P.b); xs.add(key[1])
    grid = sorted(set(v + s for v in xs for s in (-HALF, HALF) if lo <= v + s <= hi) | {lo, hi})

    def inside(p, c, s):
        l, h = c - HALF, c + HALF
        if s > 0: return l < p <= h
        if s < 0: return l <= p < h
        return l <= p <= h

    def ov(a0, a1, b0, b1): return max(F(0), min(a1, b1) - max(a0, b0))
    best = None; n = 0; tight = 0
    for cx in grid:
        for cy in grid:
            for sx in (-1, 1):
                for sy in (-1, 1):
                    if (cx == lo and sx < 0) or (cx == hi and sx > 0) or (cy == lo and sy < 0) or (cy == hi and sy > 0):
                        continue
                    v = ov(a, b, cx - HALF, cx + HALF) * ov(a, b, cy - HALF, cy + HALF)
                    for key, P in prof.items():
                        if key[0] == 'H':
                            if inside(key[1], cy, sy): v += P.G(cx + HALF) - P.G(cx - HALF)
                        else:
                            if inside(key[1], cx, sx): v += P.G(cy + HALF) - P.G(cy - HALF)
                    n += 1
                    if v == 1: tight += 1
                    if best is None or v < best[0]: best = (v, (cx, cy, sx, sy))
    if verbose:
        print(f"Lemma Z (theta = 0 face) on the box cover: {n} one-sided limit corners on a {len(grid)}^2 grid over "
              f"[{lo},{hi}]^2; min = {best[0]} = {float(best[0]):.15f} at {tuple(str(t) for t in best[1])}; "
              f"exactly tight: {tight}; {'OK' if best[0] >= 1 else '*** VIOLATED ***'}")
    return best, n, tight


# ============================================================================================== driver
KEYS = ('ADM', 'CORE', 'P1', 'MIX', 'CHAIN', 'TRI', 'PIECE', 'EMPTY', 'UNCERT', 'boxes', 'maxdepth', 'THR', 'LIN',
        'TPTS', 'SPLIT', 'cpu', 'LEB', 'CAP', 'EXACT', 'EXACT0', 'AXIS')
_CHK = None


def _work(root):
    return root, _CHK.run_box(root)


def file_sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def main():
    global _CHK
    ap = argparse.ArgumentParser()
    ap.add_argument('path')
    ap.add_argument('--depth', type=int, default=18); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--pitch', default='1/10'); ap.add_argument('--ubins', type=int, default=8)
    ap.add_argument('--exact-umax', default='1/20'); ap.add_argument('--exact-from', type=int, default=0)
    ap.add_argument('--no-exact', action='store_true')
    ap.add_argument('--cx-lo'); ap.add_argument('--cx-hi'); ap.add_argument('--cy-lo'); ap.add_argument('--cy-hi')
    ap.add_argument('--u-lo'); ap.add_argument('--u-hi')
    ap.add_argument('--resume', default=None); ap.add_argument('--progress', type=int, default=50)
    ap.add_argument('--unc-out', default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    shas = {n: file_sha(os.path.join(here, n)) for n in ('qx2_zm.py', 'zm_mixed.py', 'zeromargin.py', 'mixed_cover.py')}
    shas['input'] = file_sha(a.path)
    zs = ZM.zm_sha()
    print(f"zeromargin.py sha256 {zs} ({'= pinned' if zs == ZM.ZM_SHA_PINNED else '*** DIFFERS FROM PINNED ***'})")
    for k, v in shas.items(): print(f"sha256 {v}  {k}")
    print('argv:', ' '.join(sys.argv), flush=True)
    cov = ZM.Cover(MC.load(a.path))
    if not cov.symmetric_d4(): print("ERROR: measure not D4-invariant"); sys.exit(2)
    Ua, Ub = u_square(cov)
    print(f"container [0,{cov.m}]^2, {len(cov.points)} points, {sum(len(L['segs']) for L in cov.lines.values())} segments "
          f"on {len(cov.lines)} lines, Lebesgue square [{Ua},{Ub}]^2; total {cov.total} = {float(cov.total):.12f}")
    R0 = zm.d4_roots(cov.m, F(a.pitch), a.ubins)

    def keep(r):
        x0, x1, y0, y1, u0, u1 = r
        if a.cx_lo is not None and x1 <= F(a.cx_lo): return False
        if a.cx_hi is not None and x0 >= F(a.cx_hi): return False
        if a.cy_lo is not None and y1 <= F(a.cy_lo): return False
        if a.cy_hi is not None and y0 >= F(a.cy_hi): return False
        if a.u_lo is not None and u1 <= F(a.u_lo): return False
        if a.u_hi is not None and u0 >= F(a.u_hi): return False
        return True
    partial = any(v is not None for v in (a.cx_lo, a.cx_hi, a.cy_lo, a.cy_hi, a.u_lo, a.u_hi))
    R0 = [r for r in R0 if keep(r)]
    print(f"D4 roots [0,m/2]^2 x u in [0,1/2]: {len(R0)} roots{' (PARTIAL)' if partial else ''}; depth {a.depth}; "
          f"EXACT {'off' if a.no_exact else 'on for u1 <= ' + a.exact_umax}", flush=True)
    chk = QXChecker(cov, exact_umax=F(a.exact_umax), exact_from=a.exact_from, use_exact=not a.no_exact,
                    max_depth=a.depth, use_chain=False, cert_mode=False)
    _CHK = chk
    done = set()
    if a.resume and os.path.exists(a.resume):
        for ln in open(a.resume):
            d = json.loads(ln)
            if d.get('kind') == 'header':
                if d['sha256'] != shas: print("ERROR: resume header shas differ"); sys.exit(2)
                continue
            done.add(tuple(d['root']))
    fres = open(a.resume, 'a') if a.resume else None
    if fres and not done: fres.write(json.dumps(dict(kind='header', sha256=shas, argv=sys.argv)) + "\n"); fres.flush()
    todo = [r for r in R0 if tuple(str(v) for v in r) not in done]
    tot = {k: 0 for k in KEYS}; unc_all = []
    if a.resume and done:
        for ln in open(a.resume):
            d = json.loads(ln)
            if d.get('kind') == 'header': continue
            for k in KEYS: tot[k] = max(tot[k], d['st'].get(k, 0)) if k == 'maxdepth' else tot[k] + d['st'].get(k, 0)
            unc_all += [tuple(F(v) for v in b) for b in d['unc']]
    t0 = time.time()
    import multiprocessing as mp
    pool = mp.get_context('fork').Pool(a.nproc) if a.nproc > 1 else None
    it = pool.imap_unordered(_work, todo, chunksize=1) if pool else map(_work, todo)
    for i, (root, (st, unc, leaves)) in enumerate(it):
        for k in KEYS: tot[k] = max(tot[k], st.get(k, 0)) if k == 'maxdepth' else tot[k] + st.get(k, 0)
        unc_all += unc
        if fres:
            fres.write(json.dumps(dict(root=[str(v) for v in root], st=st, unc=[[str(v) for v in b] for b in unc])) + "\n")
            fres.flush()
        if (i + 1) % a.progress == 0:
            print(f"  {i+1}/{len(todo)} roots, {tot['boxes']} boxes, uncert {tot['UNCERT']}, {time.time()-t0:.0f}s", flush=True)
    if pool: pool.close(); pool.join()
    print(f"done in {time.time()-t0:.0f}s: boxes {tot['boxes']}, max depth {tot['maxdepth']}, CPU {tot['cpu']:.0f} s")
    print("  leaves: " + "  ".join(f"{k} {tot[k]}" for k in ('LEB', 'CAP', 'EXACT', 'EXACT0', 'AXIS', 'PIECE', 'ADM', 'P1', 'MIX',
                                                               'SPLIT', 'EMPTY', 'UNCERT')))
    if unc_all:
        print("uncertified boxes (cx0 cx1 cy0 cy1 theta0 theta1 deg):")
        for b in sorted(unc_all)[:60]:
            print("  ", *[f"{float(v):.6f}" for v in b[:4]],
                  f"{math.degrees(2*math.atan(float(b[4]))):.5f} {math.degrees(2*math.atan(float(b[5]))):.5f}")
        if len(unc_all) > 60: print(f"  ... {len(unc_all)} in total")
        if a.unc_out:
            with open(a.unc_out, 'w') as f:
                for b in unc_all: f.write(" ".join(str(v) for v in b) + "\n")
    tag = "VERIFIED-D4 (u > 0; u = 0 face: qx2_exact.py axis)" if tot['UNCERT'] == 0 else "NOT VERIFIED"
    print(tag + (" (PARTIAL)" if partial else ""))


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == 'axis':
        cov = ZM.Cover(MC.load(sys.argv[2]))
        best, n, t = axis_face(cov)
        sys.exit(0 if best[0] >= 1 else 1)
    main()
