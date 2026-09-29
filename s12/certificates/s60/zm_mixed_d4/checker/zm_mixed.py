#!/usr/bin/env python3
"""Exact checker for MIXED covers (points + uniform segments + uniform convex polygons), format v1.

Statement checked (tasks/line-cover/FORMAT.md):  for every closed unit square Q contained in [0,s]^2, at every
centre and angle,  mu(Q) >= 1,  mu = sum of the point masses, the uniform-by-length segment measures and the
uniform-by-area polygon measures; Q closed (a segment on an edge of Q counts in full).

This file IMPORTS search/zeromargin.py (never edits it; its sha256 is printed and compared with the pinned value)
and reuses its pose boxes, its admissibility machinery (Lemma A/B/C of RUNG2.md: `_adm_specs`, `_xn`,
`Checker._cond_poly`, `clip_bin`, `bin_data`) and its point primitives ADM / P1 / MIX / CHAIN unchanged.
What is new (soundness argument with proofs: search/ZM_MIXED.md):

  * PIECE bound.  For a pose box, a certified lower bound L on the segment + polygon mass captured at EVERY
    admissible pose of the box:
      - Lemma S (core of a line): the set of points p of a line that satisfy containment inequality k at every
        admissible pose of the box contains the interval of parameters t where all Bernstein coefficients of
        the Lemma-B polynomial G_k(u; p(t)) are <= 0; those coefficients are AFFINE in t, so the interval is
        exact rational.  Over the Lemma-B bound choices ('R' / wall) the convex hull of the intervals is still
        certified (the true certified set is convex).  Intersecting the four conditions gives J4, and
        mass(seg ∩ J4) is certified.  Polygons: the same with half-planes (Bernstein coefficients affine in p).
      - Lemma T (threshold / germ lemma): for two axis-parallel lines at distance exactly 1 (x = xi and
        x = xi + 1, or y = eta and y = eta + 1) every pose with theta in [0, 90deg) has ONE threshold T in the
        extended reals such that the up-line's singular inequality holds at t >= T and the down-line's at
        t <= T + u0.  Hence the pair's mass is >= min over T of a continuous piecewise-linear function, which
        is computed exactly.  This is what certifies tile germs (theta -> 0 with edges on grid lines), where
        the mass is discontinuous in pose and no fixed-witness (core) argument can terminate.
      - Corollary T' (second attempt): points lying on the germ lines join the groups as steps and are withheld
        from the point primitives for that box.
      - Lemma L / L' (axis lines): at fixed u every chord end is affine in the centre, so the lines' bound
        core + rho_up min(t_up - a_up, Delta) + rho_lo min(b_lo - t_lo, Delta) is concave in the centre: its minimum
        is at the corners of (a rectangle containing) the admissible centre rectangle, and in u it is bounded by
        Bernstein ratios.  First-order exact off the germs (the core bound loses O(box) at every chord end).
  * The points are certified by zeromargin's own primitives with a PHANTOM point of weight L that is declared
    in Q at every admissible pose (zeromargin's `inh` mechanism: a set of points already proved to lie in Q
    at every admissible pose of the box).  Every one of its primitives then proves  mu_pts(Q) + L >= 1.

Usage:
  python3 search/zm_mixed.py cert FILE [--d4 | --full] [--disj --chain-from 0] [--depth 18] [--nproc 5] ...
  python3 search/zm_mixed.py pose FILE --cx C --cy C --u U         exact mu(Q) at one rational pose
  python3 search/zm_mixed.py stress FILE LEAFDUMP [--per 20]        float/exact re-check of certified leaves
  python3 search/zm_mixed.py selftest                               randomised tests of Lemmas S and T
"""
import sys, os, time, math, argparse, hashlib, bisect, random
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import zeromargin as zm
import mixed_cover as MC

ZM_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeromargin.py')
ZM_SHA_PINNED = '640fe453c1a32f4aa580ca2b1261c6406923a4d7c131f65604a432c7fc2086ab'
HALF = F(1, 2)
REACH = F(7072, 10000)            # > sqrt(2)/2: every point of Q(c, theta) is within this of c


def zm_sha():
    return hashlib.sha256(open(ZM_PATH, 'rb').read()).hexdigest()


# =========================================================================================== the cover
class Cover:
    """Exact (Fraction) form of a mixed cover.  Segments are grouped by supporting line; each line has a base
    point P0 and a direction d (P(t) = P0 + t d), and its segments are parameter intervals [t0, t1] with a
    density dens = mass / (t1 - t0) per unit of t.  Axis lines use t = y (vertical, d = (0,1)) or t = x
    (horizontal, d = (1,0)), so the two lines of a germ pair share their parameter."""

    def __init__(self, cv):
        MC.validate(cv)
        self.raw = cv
        D, W = cv['D'], cv['W']
        self.D, self.Wd = D, W
        self.m = F(cv['s_num'], cv['s_den'])
        self.S = cv['s_num'] * D // cv['s_den']            # container side in coordinate units
        self.points = [(F(X, D), F(Y, D), F(w, W)) for X, Y, w in cv['points']]
        self.total = MC.total(cv)
        # ---- lines
        self.lines = {}
        for (X0, Y0, X1, Y1, w) in cv['segments']:
            if w == 0: continue
            if X0 == X1:
                key = ('V', F(X0, D)); P0 = (F(X0, D), F(0)); d = (F(0), F(1))
                t0, t1 = sorted((F(Y0, D), F(Y1, D)))
            elif Y0 == Y1:
                key = ('H', F(Y0, D)); P0 = (F(0), F(Y0, D)); d = (F(1), F(0))
                t0, t1 = sorted((F(X0, D), F(X1, D)))
            else:
                a, b = X1 - X0, Y1 - Y0
                g = math.gcd(a, b); a //= g; b //= g
                if a < 0 or (a == 0 and b < 0): a, b = -a, -b
                c = b * X0 - a * Y0                          # line: b X - a Y = c (integer coordinates)
                key = ('G', a, b, c)
                if key in self.lines:
                    P0, d = self.lines[key]['P0'], self.lines[key]['d']
                else:
                    P0 = (F(X0, D), F(Y0, D)); d = (F(a, D), F(b, D))
                dd = d[0] * d[0] + d[1] * d[1]
                ta = ((F(X0, D) - P0[0]) * d[0] + (F(Y0, D) - P0[1]) * d[1]) / dd
                tb = ((F(X1, D) - P0[0]) * d[0] + (F(Y1, D) - P0[1]) * d[1]) / dd
                t0, t1 = min(ta, tb), max(ta, tb)
            L = self.lines.setdefault(key, dict(key=key, P0=P0, d=d, segs=[]))
            L['segs'].append((t0, t1, F(w, W) / (t1 - t0)))
        for L in self.lines.values():
            L['segs'].sort()
            self._prefix(L)
            P0, d = L['P0'], L['d']
            L['f'] = (float(P0[0]), float(P0[1]), float(d[0]), float(d[1]),
                      min(s[0] for s in L['segs']), max(s[1] for s in L['segs']))
        # points lying exactly on an axis line that carries segments (for Lemma T with points), index = position
        # in self.points (= the point's index in the zeromargin checker)
        self.line_points = {}
        for i, (px, py, w) in enumerate(self.points):
            for key in (('V', px), ('H', py)):
                if key in self.lines and w > 0:
                    t = py if key[0] == 'V' else px
                    self.line_points.setdefault(key, []).append((t, w, i))
        self.V = {k[1]: L for k, L in self.lines.items() if k[0] == 'V'}
        self.H = {k[1]: L for k, L in self.lines.items() if k[0] == 'H'}
        # ---- polygons
        self.polys = []
        for (w, vs) in cv['polygons']:
            if w == 0: continue
            check_simple_convex(vs)
            V = [(F(X, D), F(Y, D)) for X, Y in vs]
            A = poly_area(V)
            xs = [float(v[0]) for v in V]; ys = [float(v[1]) for v in V]
            self.polys.append(dict(V=V, w=F(w, W), area=A, bb=(min(xs), max(xs), min(ys), max(ys))))
        self.has_pieces = bool(self.lines) or bool(self.polys)

    @staticmethod
    def _prefix(L):
        L['tot'] = sum((s[1] - s[0]) * s[2] for s in L['segs'])

    # ---- symmetry of the MEASURE (canonical multisets: sufficient for invariance)
    def _canon(self, g):
        from collections import Counter
        P = Counter()
        for X, Y, w in self.raw['points']:
            if w: P[g(X, Y)] += w
        Sg = Counter()
        for X0, Y0, X1, Y1, w in self.raw['segments']:
            if w: Sg[(frozenset((g(X0, Y0), g(X1, Y1))), w)] += 1
        Pg = Counter()
        for w, vs in self.raw['polygons']:
            if w: Pg[(frozenset(g(X, Y) for X, Y in vs), w)] += 1
        return P, Sg, Pg

    def invariant(self, g):
        return self._canon(lambda X, Y: (X, Y)) == self._canon(g)

    def symmetric_d2(self):
        S = self.S
        return self.invariant(lambda X, Y: (S - X, Y)) and self.invariant(lambda X, Y: (X, S - Y))

    def symmetric_d4(self):
        S = self.S
        return self.invariant(lambda X, Y: (S - X, Y)) and self.invariant(lambda X, Y: (Y, X))


def swapped(cv):
    """the image of a cover under the reflection (x, y) -> (y, x)."""
    return dict(cv, points=[(Y, X, w) for X, Y, w in cv['points']],
                segments=[(Y0, X0, Y1, X1, w) for X0, Y0, X1, Y1, w in cv['segments']],
                polygons=[(w, [(Y, X) for X, Y in reversed(vs)]) for w, vs in cv['polygons']])


# =========================================================================================== small exact geometry
def poly_area(V):
    a = 0
    for i in range(len(V)):
        x0, y0 = V[i]; x1, y1 = V[(i + 1) % len(V)]
        a += x0 * y1 - x1 * y0
    return a / 2


def check_simple_convex(vs):
    """mixed_cover.validate checks that every turn is a left turn (or straight) and that the signed area is
    positive; that still admits a polygon winding twice round.  Require in addition: no repeated vertex and the
    vertex set, without its collinear vertices, is exactly the convex hull in CCW order."""
    if len(set(vs)) != len(vs): raise ValueError('polygon with a repeated vertex')
    k = len(vs)
    strict = []
    for i in range(k):
        (x0, y0), (x1, y1), (x2, y2) = vs[i - 1], vs[i], vs[(i + 1) % k]
        if (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1) != 0: strict.append(vs[i])
    H = convex_hull(vs)
    if len(H) != len(strict) or set(H) != set(strict): raise ValueError('polygon is not a simple convex polygon')
    j = strict.index(H[0])
    if strict[j:] + strict[:j] != H: raise ValueError('polygon vertices are not in CCW hull order')


def convex_hull(pts):
    """Andrew's monotone chain, CCW, collinear points dropped; exact for Fractions / ints."""
    P = sorted(set(pts))
    if len(P) <= 2: return P
    def cross(o, a, b): return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, hi = [], []
    for p in P:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(P):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], p) <= 0: hi.pop()
        hi.append(p)
    return lo[:-1] + hi[:-1]


def clip_halfplane(V, a, b, c):
    """convex polygon V ∩ {a x + b y + c <= 0}  (Sutherland-Hodgman, exact)."""
    out = []
    n = len(V)
    for i in range(n):
        P = V[i]; Q = V[(i + 1) % n]
        fp = a * P[0] + b * P[1] + c; fq = a * Q[0] + b * Q[1] + c
        if fp <= 0: out.append(P)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    return out


def clip_convex(V, K):
    """V ∩ K for convex polygons (K CCW): clip by every edge half-plane of K."""
    n = len(K)
    for i in range(n):
        if len(V) < 3: return []
        (x0, y0), (x1, y1) = K[i], K[(i + 1) % n]
        # inside of a CCW edge: cross(e, p - p0) >= 0  <=>  -(x1-x0)(y-y0) + (y1-y0)(x-x0) <= 0
        a = y1 - y0; b = -(x1 - x0); c = -(a * x0 + b * y0)
        V = clip_halfplane(V, a, b, c)
    return V if len(V) >= 3 else []


# =========================================================================================== intervals
# An interval is (lo, hi) with lo/hi a Fraction or None (= -inf / +inf); the empty set is None.
def iv_and(a, b):
    if a is None or b is None: return None
    lo = a[0] if b[0] is None else (b[0] if a[0] is None else max(a[0], b[0]))
    hi = a[1] if b[1] is None else (b[1] if a[1] is None else min(a[1], b[1]))
    if lo is not None and hi is not None and lo > hi: return None
    return (lo, hi)


def iv_hull(a, b):
    if a is None: return b
    if b is None: return a
    lo = None if (a[0] is None or b[0] is None) else min(a[0], b[0])
    hi = None if (a[1] is None or b[1] is None) else max(a[1], b[1])
    return (lo, hi)


def seg_len_in(t0, t1, iv):
    """length of [t0, t1] ∩ iv."""
    if iv is None: return 0
    lo = t0 if iv[0] is None else max(t0, iv[0])
    hi = t1 if iv[1] is None else min(t1, iv[1])
    return hi - lo if hi > lo else 0


# =========================================================================================== Bernstein (Lemma C)
_C = [[1, 0, 0, 0, 0], [1, 1, 0, 0, 0], [1, 2, 1, 0, 0], [1, 3, 3, 1, 0], [1, 4, 6, 4, 1]]
_BW = [[F(_C[i][j], _C[4][j]) for j in range(5)] for i in range(5)]


def bern(a, u0, h):
    """The 5 Bernstein coefficients (degree 4) on [u0, u0 + h] of sum_k a[k] u^k (deg <= 4): the polynomial is
    a convex combination of them at every u of the interval (RUNG2.md Lemma C).  LINEAR in a."""
    b = []
    hp = F(1)
    for j in range(5):
        sj = F(0)
        up = F(1)
        for k in range(j, 5):
            if a[k]: sj += a[k] * _C[k][j] * up
            up *= u0
        b.append(sj * hp); hp *= h
    return [sum(_BW[i][j] * b[j] for j in range(i + 1)) for i in range(5)]


def cond_poly(px, py, xk, yk, cond, m):
    """zeromargin's Lemma-B polynomial G_cond(u) for the point (px, py) and the bound choice (xk, yk), exactly
    as zeromargin.Checker._adm_cond_ok_ref builds it: G <= 0 on the bin  =>  the inequality holds at every
    admissible pose of the box (Lemma A).  G is AFFINE in (px, py)."""
    Xn = zm._xn(xk[0], xk[1], m); Yn = zm._xn(yk[0], yk[1], m)
    U = (2 * px - Xn[0], -Xn[1], 2 * px - Xn[2])
    V = (2 * py - Yn[0], -Yn[1], 2 * py - Yn[2])
    return zm.Checker._cond_poly(U, V, cond, xk[0] == 'R' and yk[0] == 'R')


def _slots(specs, cond):
    Axs, Bxs, Ays, Bys = specs
    return ((Axs, Ays), (Bxs, Bys), (Bxs, Ays), (Axs, Bys))[cond]


def line_cond_iv(P0, d, specs, cond, u0, h, m):
    """Lemma S: an interval of parameters t such that P0 + t d satisfies inequality `cond` at every admissible
    pose of the box.  For each bound choice the Bernstein coefficients are beta_i(t) = b0_i + t bd_i, so
    {t : all beta_i(t) <= 0} is an exact interval; the hull over the choices is certified because the true
    certified set is convex and contains each of them."""
    xs, ys = _slots(specs, cond)
    res = None
    for xk in xs:
        for yk in ys:
            g0 = cond_poly(P0[0], P0[1], xk, yk, cond, m)
            g1 = cond_poly(P0[0] + d[0], P0[1] + d[1], xk, yk, cond, m)
            b0 = bern(g0, u0, h)
            bd = bern([g1[i] - g0[i] for i in range(5)], u0, h)
            iv = (None, None)
            for i in range(5):
                if bd[i] > 0: iv = iv_and(iv, (None, -b0[i] / bd[i]))
                elif bd[i] < 0: iv = iv_and(iv, (-b0[i] / bd[i], None))
                elif b0[i] > 0: iv = None
                if iv is None: break
            if iv is not None: res = iv_hull(res, iv)
    return res


def cond_region(specs, cond, u0, h, m, base, frame):
    """Lemma S for polygons: a convex polygon of points satisfying inequality `cond` at every admissible pose
    of the box: the hull over the bound choices of {p in frame : all Bernstein coefficients <= 0}; each
    coefficient is affine in p, so each choice gives an intersection of 5 half-planes."""
    xs, ys = _slots(specs, cond)
    pts = []
    bx, by = base
    for xk in xs:
        for yk in ys:
            g0 = cond_poly(bx, by, xk, yk, cond, m)
            gx = cond_poly(bx + 1, by, xk, yk, cond, m)
            gy = cond_poly(bx, by + 1, xk, yk, cond, m)
            b0 = bern(g0, u0, h)
            bX = bern([gx[i] - g0[i] for i in range(5)], u0, h)
            bY = bern([gy[i] - g0[i] for i in range(5)], u0, h)
            Vp = list(frame)
            for i in range(5):
                # b0 + (x - bx) bX + (y - by) bY <= 0
                if bX[i] == 0 and bY[i] == 0:
                    if b0[i] > 0: Vp = []
                else:
                    Vp = clip_halfplane(Vp, bX[i], bY[i], b0[i] - bX[i] * bx - bY[i] * by)
                if len(Vp) < 3: Vp = []; break
            pts += Vp
    if len(pts) < 3: return []
    H = convex_hull(pts)
    return H if len(H) >= 3 else []


# =========================================================================================== line mass functions
class LineMass:
    """mass of (segments of one line) ∩ J ∩ (-inf, y]  as a continuous piecewise-linear function of y, with
    exact prefix sums: G(y) = sum_i dens_i |[e0_i, e1_i] ∩ (-inf, y]|, the e's already clipped to J."""

    def __init__(self, pieces):
        self.p = pieces
        a = sorted((e0, dn) for e0, e1, dn in pieces)
        b = sorted((e1, dn) for e0, e1, dn in pieces)
        self.e0 = [x for x, _ in a]; self.e1 = [x for x, _ in b]
        self.D0 = [F(0)]; self.M0 = [F(0)]
        for x, dn in a: self.D0.append(self.D0[-1] + dn); self.M0.append(self.M0[-1] + dn * x)
        self.D1 = [F(0)]; self.M1 = [F(0)]
        for x, dn in b: self.D1.append(self.D1[-1] + dn); self.M1.append(self.M1[-1] + dn * x)
        self.tot = sum((e1 - e0) * dn for e0, e1, dn in pieces)
        self.bps = self.e0 + self.e1

    def below(self, y):
        """G(y); y may be None meaning +inf... use below_inf / 0 for the infinite cases."""
        i = bisect.bisect_left(self.e0, y)          # e0 < y
        j = bisect.bisect_left(self.e1, y)          # e1 < y
        return (y * self.D0[i] - self.M0[i]) - (y * self.D1[j] - self.M1[j])

    def above(self, y):
        return self.tot - self.below(y)


def line_pieces(L, J, lo_t, hi_t):
    """the segments of line L restricted to the interval J and to the reach window [lo_t, hi_t] (floats)."""
    out = []
    if J is None: return out
    for t0, t1, dn in L['segs']:
        if float(t1) < lo_t or float(t0) > hi_t: continue
        a = t0 if J[0] is None else max(t0, J[0])
        b = t1 if J[1] is None else min(t1, J[1])
        if b > a: out.append((a, b, dn))
    return out


# =========================================================================================== Lemma L (linear chord ends)
# Polynomials in u are coefficient lists (low degree first).  S = 2u, C = 1 - u^2 (both > 0 for 0 < u < 1).
PS = [F(0), F(2)]
PC = [F(1), F(0), F(-1)]


def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def pscale(p, k):
    return [k * a for a in p]


def pmul(p, q):
    r = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q): r[i + j] += a * b
    return r


def peval(p, x):
    v = 0.0
    for a in reversed(p): v = v * x + float(a)
    return v


# threshold numerators of the four inequalities on an axis line at a corner (see ZM_MIXED.md sec 2, Lemma L):
#   V line x = xi, t = y, a = xi - cx:  t_k = cy + num_k / den_k
#   H line y = eta, t = x, b = eta - cy: t_k = cx + num_k / den_k
def thr_num(fam, k, a):
    h = HALF
    if fam == 'V':
        return {0: ([h - a, F(0), h + a], 'S'), 1: ([-h - a, F(0), a - h], 'S'),
                2: ([h, 2 * a, h], 'C'), 3: ([-h, 2 * a, -h], 'C')}[k]
    return {0: ([h, -2 * a, h], 'C'), 1: ([-h, -2 * a, -h], 'C'),
            2: ([a - h, F(0), -a - h], 'S'), 3: ([a + h, F(0), h - a], 'S')}[k]


UPLO = {'V': ({0: 'S', 2: 'C'}, {1: 'S', 3: 'C'}), 'H': ({0: 'C', 3: 'S'}, {1: 'C', 2: 'S'})}


def end_minorant(pcs, LM, a, D, up, xmid):
    """For the gain g(x) = mass in [a, a + x] (up) or [a - x, a] (down), 0 <= x <= D, an affine l(x) = s x + i with
    s >= 0 and l <= g on [0, D]: the edge of the lower convex hull of g's vertices (g is piecewise linear, so its
    greatest convex minorant is that hull) that contains xmid (a float guess of where the chord end sits).
    Returns (s, i, g(D)).  Used in Lemma L: gain >= min(l(x), g(D))."""
    if D <= 0: return F(0), F(0), F(0)
    if up:
        g = lambda x: LM.below(a + x) - LM.below(a)
        xs = [e - a for e in LM.bps if a < e < a + D]
    else:
        g = lambda x: LM.below(a) - LM.below(a - x)
        xs = [a - e for e in LM.bps if a - D < e < a]
    P = sorted(set([F(0), D] + xs))
    V = [(x, g(x)) for x in P]
    H = []
    for pnt in V:                      # lower hull, left to right
        while len(H) >= 2 and (H[-1][0] - H[-2][0]) * (pnt[1] - H[-2][1]) - (H[-1][1] - H[-2][1]) * (pnt[0] - H[-2][0]) <= 0:
            H.pop()
        H.append(pnt)
    xm = min(max(xmid, 0.0), float(D))
    for i in range(len(H) - 1):
        if float(H[i + 1][0]) >= xm or i == len(H) - 2:
            (x0, y0), (x1, y1) = H[i], H[i + 1]
            sl = (y1 - y0) / (x1 - x0)
            return sl, y0 - sl * x0, V[-1][1]
    return F(0), F(0), V[-1][1]


def _tk_float(fam, k, line, cx, cy, u):
    num, dt = thr_num(fam, k, F(0)); n1, _ = thr_num(fam, k, F(1))
    a = float(line) - (cx if fam == 'V' else cy)
    nv = [float(num[i]) + a * float(n1[i] - num[i]) for i in range(3)]
    val = nv[0] + nv[1] * u + nv[2] * u * u
    den = 2 * u if dt == 'S' else 1 - u * u
    return (cy if fam == 'V' else cx) + val / den


def move_bound(fam, line, ks, ref, up, box, m):
    """an exact upper bound, over the box, of  min_{k in ks} t_k - ref  (up) or  ref - max_{k in ks} t_k  (lo):
    min over k of the max over sub-bins and corners of a containing rectangle (t_k is affine in the centre) of
    the Bernstein-ratio upper bound.  None if unavailable."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    ch = [corner_choices(cx0, True, m, u0, u1), corner_choices(cx1, False, m, u0, u1),
          corner_choices(cy0, True, m, u0, u1), corner_choices(cy1, False, m, u0, u1)]
    pts = sorted(set([u0, u1] + [q for c in ch for (a, b, _) in c for q in (a, b)]))
    res = None
    for k in ks:
        mk = None
        for j in range(len(pts) - 1):
            v0, v1 = pts[j], pts[j + 1]
            X0, X1, Y0, Y1 = [next(f for (a, b, f) in c if a <= v0 and v1 <= b) for c in ch]
            for X in (X0, X1):
                for Y in (Y0, Y1):
                    t = thr_rf(fam, k, line, X, Y)
                    e = rf_add(t, {K0: [-ref]}) if up else rf_add({K0: [ref]}, t, -1)
                    ub = rf_bound(e, v0, v1, upper=True)
                    if ub is None: mk = None; break
                    mk = ub if mk is None else max(mk, ub)
                if mk is None: break
            if mk is None: break
        if mk is not None and (res is None or mk < res): res = mk
    return res


def lemma_l_data(L, key, I, w, box, m=None):
    """eligibility and data of Lemma L for one axis line in one box, or None."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    fam = key[0]
    ups, los = UPLO[fam]
    typed_up = [k for k, dt in ups.items() if dt == 'C' or u0 > 0]
    typed_lo = [k for k, dt in los.items() if dt == 'C' or u0 > 0]
    a_up = None; b_lo = None
    for k in typed_up:
        iv = I[k]
        if iv is None or iv[0] is not None or iv[1] is None: return None
        a_up = iv[1] if a_up is None else min(a_up, iv[1])
    for k in typed_lo:
        iv = I[k]
        if iv is None or iv[1] is not None or iv[0] is None: return None
        b_lo = iv[0] if b_lo is None else max(b_lo, iv[0])
    Delta = 2 * (cx1 - cx0 + cy1 - cy0) + 8 * (u1 - u0) + F(1, 10 ** 6)
    vertex = not (b_lo < a_up)
    Dup = Dlo = Delta
    if not vertex and m is not None:        # tight caps: how far the chord ends can actually move (any Delta > 0 is sound)
        mu = move_bound(fam, key[1], typed_up, a_up, True, box, m)
        ml = move_bound(fam, key[1], typed_lo, b_lo, False, box, m)
        if mu is not None: Dup = max(min(Delta, mu), F(1, 10 ** 12))
        if ml is not None: Dlo = max(min(Delta, ml), F(1, 10 ** 12))
    if vertex:                  # Lemma L' (short chord, e.g. a vertex of Q crossing the line): see ZM_MIXED.md
        lo_r, hi_r = min(a_up, b_lo) - Delta, max(a_up, b_lo) + Delta
    else:
        lo_r, hi_r = b_lo - Dlo, a_up + Dup
    for k in range(4):
        if k in typed_up or k in typed_lo: continue
        iv = I[k]                                   # an untyped inequality must hold on the whole range
        if iv is None or (iv[0] is not None and iv[0] > lo_r) or (iv[1] is not None and iv[1] < hi_r): return None
    pcs = line_pieces(L, (None, None), w[0], w[1])
    if not pcs: return None
    if vertex:                  # Lemma L' with hull minorant / majorant (ZM_MIXED.md): mass >= l1(min(r_up,B)-p) - l2(max(r_lo,A)-p)
        LMv = LineMass(pcs)
        cxm, cym, um = float(cx0 + cx1) / 2, float(cy0 + cy1) / 2, max(float(u0 + u1) / 2, 1e-9)
        rf_up = min(_tk_float(fam, k, key[1], cxm, cym, um) for k in typed_up)
        rf_lo = max(_tk_float(fam, k, key[1], cxm, cym, um) for k in typed_lo)
        pst = F((min(max(rf_up, float(lo_r)), float(hi_r)) + min(max(rf_lo, float(lo_r)), float(hi_r))) / 2).limit_denominator(10 ** 9)
        pst = min(max(pst, lo_r), hi_r)
        g = lambda x: LMv.below(pst + x) - LMv.below(pst)
        xs = sorted(set([lo_r - pst, hi_r - pst] + [e - pst for e in LMv.bps if lo_r < e < hi_r]))
        V = [(x, g(x)) for x in xs]
        def hull(V, lower):
            H = []
            for q in V:
                while len(H) >= 2:
                    cr = (H[-1][0] - H[-2][0]) * (q[1] - H[-2][1]) - (H[-1][1] - H[-2][1]) * (q[0] - H[-2][0])
                    if (cr <= 0) if lower else (cr >= 0): H.pop()
                    else: break
                H.append(q)
            return H
        def edge(H, xf):
            for i in range(len(H) - 1):
                if float(H[i + 1][0]) >= xf or i == len(H) - 2:
                    (x0, y0), (x1, y1) = H[i], H[i + 1]
                    sl = (y1 - y0) / (x1 - x0); return sl, y0 - sl * x0
        s1, i1 = edge(hull(V, True), rf_up - float(pst))
        s2, i2 = edge(hull(V, False), rf_lo - float(pst))
        if s1 <= 0 and s2 <= 0: return None
        return dict(fam=fam, line=key[1], up=typed_up, lo=typed_lo, vertex=True, A=lo_r, Bv=hi_r, p=pst,
                    s1=s1, i1=i1, s2=s2, i2=i2, rho=s1, core=F(0))
    LM = LineMass(pcs)
    core = LM.below(a_up) - LM.below(b_lo)
    cxm, cym, um = float(cx0 + cx1) / 2, float(cy0 + cy1) / 2, max(float(u0 + u1) / 2, 1e-9)
    xu = min(_tk_float(fam, k, key[1], cxm, cym, um) for k in typed_up) - float(a_up)
    xl = float(b_lo) - max(_tk_float(fam, k, key[1], cxm, cym, um) for k in typed_lo)
    su, iu, cu = end_minorant(pcs, LM, a_up, Dup, True, xu)
    sl, il, cl = end_minorant(pcs, LM, b_lo, Dlo, False, xl)
    return dict(fam=fam, line=key[1], up=typed_up, lo=typed_lo, a_up=a_up, b_lo=b_lo, Dup=Dup, Dlo=Dlo,
                core=core, s_up=su, i_up=iu, cap_up=cu, s_lo=sl, i_lo=il, cap_lo=cl)


# Rational functions of u as dicts {(eS, eC, eN): poly}: value = sum_key poly(u) / (S^eS C^eC N^eN),
# S = 2u, C = 1 - u^2, N = 1 + u^2 (all > 0 on a bin with 0 < u < 1; S only used when u0 > 0).
PN = [F(1), F(0), F(1)]
K0 = (0, 0, 0)


def rf_add(x, y, sgn=1):
    r = {k: list(v) for k, v in x.items()}
    for k, p in y.items(): r[k] = padd(r.get(k, [F(0)]), pscale(p, sgn))
    return r


def rf_scale(x, c):
    return {k: pscale(p, c) for k, p in x.items()}


def rf_mulpoly_den(x, q, den):
    """x * q / den, den a key (eS, eC, eN) of exponents to add."""
    return {(k[0] + den[0], k[1] + den[1], k[2] + den[2]): pmul(p, q) for k, p in x.items()}


def rf_eval(x, u):
    v = 0.0
    S, C, N = 2 * u, 1 - u * u, 1 + u * u
    for (a, b, c), p in x.items():
        if not any(p): continue
        d = (S ** a) * (C ** b) * (N ** c)
        v += peval(p, u) / d if d > 0 else float('inf')
    return v


def bern_n(a, u0, h, n):
    """degree-n Bernstein coefficients on [u0, u0 + h] of sum_k a[k] u^k (deg <= n)."""
    a = list(a) + [F(0)] * (n + 1 - len(a))
    b = []
    hp = F(1)
    for j in range(n + 1):
        sj = F(0); up = F(1)
        for k in range(j, n + 1):
            if a[k]: sj += a[k] * math.comb(k, j) * up
            up *= u0
        b.append(sj * hp); hp *= h
    return [sum(F(math.comb(i, j), math.comb(n, j)) * b[j] for j in range(i + 1)) for i in range(n + 1)]


def rf_bound(x, u0, u1, upper=False):
    """exact lower (upper) bound of the rational function over [u0, u1]: P/Den with Den = S^A C^B N^E the common
    denominator; if every degree-n Bernstein coefficient of Den is > 0, then P - lam Den = sum_i B_i (beta_i(P) -
    lam beta_i(Den)) >= 0 for lam = min_i beta_i(P)/beta_i(Den) (<= 0 for the max).  None if unavailable."""
    keys = [k for k, p in x.items() if any(p)]
    if not keys: return F(0)
    A = max(k[0] for k in keys); Bc = max(k[1] for k in keys); E = max(k[2] for k in keys)
    if A > 0 and u0 <= 0: return None
    def mono(a, b, c):
        r = [F(1)]
        for _ in range(a): r = pmul(r, PS)
        for _ in range(b): r = pmul(r, PC)
        for _ in range(c): r = pmul(r, PN)
        return r
    Den = mono(A, Bc, E)
    P = [F(0)]
    for k in keys: P = padd(P, pmul(x[k], mono(A - k[0], Bc - k[1], E - k[2])))
    n = max(len(P), len(Den)) - 1
    h = u1 - u0
    bP = bern_n(P, u0, h, n); bD = bern_n(Den, u0, h, n)
    if any(v <= 0 for v in bD): return None
    rs = [bP[i] / bD[i] for i in range(n + 1)]
    return max(rs) if upper else min(rs)


def corner_choices(c, lo_side, m, u0, u1):
    """For the centre bound max(c0, w/2) (lower side) or min(c1, m - w/2) (upper side): BOTH the constant c and the
    wall term are sound replacements (each gives a rectangle CONTAINING the admissible one at every u), so the
    choice only affects tightness.  Returns the sub-bins of [u0, u1] (split near where w(u)/2 crosses the
    constant, rational split points) with the tighter choice on each, as rational functions of u."""
    const = {K0: [c]}
    wall = {(0, 0, 1): [HALF, F(1), -HALF]} if lo_side else {K0: [m], (0, 0, 1): [-HALF, F(-1), HALF]}
    target = 2 * float(c) if lo_side else 2 * float(m - c)        # w(u) = target at the crossing
    splits = []
    if 0 < target / math.sqrt(2) < 1:
        for th in (math.asin(target / math.sqrt(2)) - math.pi / 4, math.pi / 4 - math.asin(target / math.sqrt(2)) + math.pi / 2):
            if 0 < th < math.pi / 2:
                uc = F(math.tan(th / 2)).limit_denominator(10 ** 9)
                if u0 < uc < u1: splits.append(uc)
    pts = [u0] + sorted(set(splits)) + [u1]
    out = []
    for i in range(len(pts) - 1):
        um = (float(pts[i]) + float(pts[i + 1])) / 2
        wv = ((1 - um * um) + 2 * um) / (1 + um * um)
        if lo_side: ch = const if float(c) >= wv / 2 else wall
        else: ch = const if float(c) <= float(m) - wv / 2 else wall
        out.append((pts[i], pts[i + 1], ch))
    return out


def thr_rf(fam, k, line, X, Y):
    """t_k as a rational function of u at the corner (X, Y) (rational functions): base + num_k(a)/den_k with
    num_k = q0 + a q1 (see thr_num) and a = line - (X for V, Y for H), base = Y for V, X for H."""
    n0, dt = thr_num(fam, k, F(0))
    n1, _ = thr_num(fam, k, F(1))
    q1 = [n1[i] - n0[i] for i in range(3)]
    den = (1, 0, 0) if dt == 'S' else (0, 1, 0)
    coord, base = (X, Y) if fam == 'V' else (Y, X)
    r = dict(base)
    r = rf_add(r, {den: padd(n0, pscale(q1, line))})          # (q0 + line q1) / den
    r = rf_add(r, rf_mulpoly_den(coord, q1, den), -1)          # - coord q1 / den
    return r


def _line_opts(d, X, Y):
    if d.get('vertex'):         # l1(min(min_up t_k, B) - p) - l2(max(max_lo t_j, A) - p)
        p_, s1, i1, s2, i2 = d['p'], d['s1'], d['i1'], d['s2'], d['i2']
        up = [{K0: [s1 * (d['Bv'] - p_) + i1]}] + \
             [rf_add(rf_scale(thr_rf(d['fam'], k, d['line'], X, Y), s1), {K0: [i1 - s1 * p_]}) for k in d['up']]
        lo = [{K0: [-(s2 * (d['A'] - p_) + i2)]}] + \
             [rf_add(rf_scale(thr_rf(d['fam'], k, d['line'], X, Y), -s2), {K0: [s2 * p_ - i2]}) for k in d['lo']]
        return [up, lo]
    out = []
    for side, ks in (('up', d['up']), ('lo', d['lo'])):
        s_, i_, cap = (d['s_up'], d['i_up'], d['cap_up']) if side == 'up' else (d['s_lo'], d['i_lo'], d['cap_lo'])
        if s_ == 0 and i_ == 0 and cap == 0: continue
        opts = [{K0: [cap]}]
        for k in ks:
            t = thr_rf(d['fam'], k, d['line'], X, Y)
            if side == 'up': x = rf_add(t, {K0: [-d['a_up']]})          # x = t_k - a_up
            else: x = rf_add({K0: [d['b_lo']]}, t, -1)                   # x = b_lo - t_k
            opts.append(rf_add(rf_scale(x, s_), {K0: [i_]}))             # s x + i
        out.append(opts)
    return out


def lemma_l_joint(datas, box, B, m):
    """Lemma L, joint bound (Corollary L): on each sub-bin, min over the 4 corners of a rectangle containing the
    admissible centre rectangle at fixed u (box sides or wall terms, rational in u) of a lower bound over u of
    sum_l [core_l + min(up options) + min(lo options)].  None if some Bernstein bound is unavailable."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    ch = [corner_choices(cx0, True, m, u0, u1), corner_choices(cx1, False, m, u0, u1),
          corner_choices(cy0, True, m, u0, u1), corner_choices(cy1, False, m, u0, u1)]
    pts = sorted(set([u0, u1] + [q for c in ch for (a, b, _) in c for q in (a, b)]))
    core = sum((d['core'] for d in datas), F(0))
    best = None
    for j in range(len(pts) - 1):
        v0, v1 = pts[j], pts[j + 1]
        pick = [next(f for (a, b, f) in c if a <= v0 and v1 <= b) for c in ch]
        um = (float(v0) + float(v1)) / 2
        for X in (pick[0], pick[1]):
            for Y in (pick[2], pick[3]):
                tot = {K0: [core]}
                slack = F(0)
                for d in datas:
                    for opts in _line_opts(d, X, Y):
                        vals = [rf_eval(o, um) for o in opts]
                        i0 = min(range(len(opts)), key=lambda i: vals[i])
                        tot = rf_add(tot, opts[i0])
                        sl = F(0)
                        for i, o in enumerate(opts):
                            if i == i0: continue
                            mu = rf_bound(rf_add(opts[i0], o, -1), v0, v1, upper=True)
                            if mu is None: return None
                            if mu > sl: sl = mu
                        slack += sl
                lb = rf_bound(tot, v0, v1)
                if lb is None: return None
                v = lb - slack
                if best is None or v < best: best = v
    return best

# =========================================================================================== Lemma R (region-wise)
def g_coeffs(kind, px, py):
    """zeromargin's violation polynomial of a point, G = alpha(u) + beta(u) cx + gamma(u) cy (G <= 0 <=> the
    inequality `kind` holds; zeromargin.Checker._gcoef with a = px - cx, b = py - cy)."""
    Cp, Sp, Np = PC, PS, PN
    if kind == 0:      # 2(aC + bS) - N
        al = padd(padd(pscale(Cp, 2 * px), pscale(Sp, 2 * py)), pscale(Np, -1)); be = pscale(Cp, -2); ga = pscale(Sp, -2)
    elif kind == 1:    # -2(aC + bS) - N
        al = padd(padd(pscale(Cp, -2 * px), pscale(Sp, -2 * py)), pscale(Np, -1)); be = pscale(Cp, 2); ga = pscale(Sp, 2)
    elif kind == 2:    # -2aS + 2bC - N
        al = padd(padd(pscale(Sp, -2 * px), pscale(Cp, 2 * py)), pscale(Np, -1)); be = pscale(Sp, 2); ga = pscale(Cp, -2)
    else:              # 2aS - 2bC - N
        al = padd(padd(pscale(Sp, 2 * px), pscale(Cp, -2 * py)), pscale(Np, -1)); be = pscale(Sp, -2); ga = pscale(Cp, 2)
    return al, be, ga


def g_at(gc, X, Y):
    al, be, ga = gc
    return rf_add(rf_add({K0: al}, rf_mulpoly_den(X, be, K0)), rf_mulpoly_den(Y, ga, K0))


def _div_by(x, poly):
    """x / poly for poly = k*C or k*S (the beta/gamma of g_coeffs); None otherwise."""
    if len(poly) == 3 and poly[1] == 0 and poly[0] == -poly[2] and poly[0] != 0:      # k (1 - u^2)
        return rf_mulpoly_den(x, [1 / poly[0]], (0, 1, 0))
    if len(poly) == 2 and poly[0] == 0 and poly[1] != 0:                               # 2k u = k S
        return rf_mulpoly_den(x, [2 / poly[1]], (1, 0, 0))
    return None


def phi_at(datas, X, Y, v0, v1):
    """Corollary L's lower bound over [v0, v1] of the lines' Lemma L right-hand sides at the centre (X(u), Y(u))."""
    um = (float(v0) + float(v1)) / 2
    tot = {K0: [sum((d['core'] for d in datas), F(0))]}
    slack = F(0)
    for d in datas:
        for opts in _line_opts(d, X, Y):
            vals = [rf_eval(o, um) for o in opts]
            i0 = min(range(len(opts)), key=lambda i: vals[i])
            tot = rf_add(tot, opts[i0])
            sl = F(0)
            for i, o in enumerate(opts):
                if i == i0: continue
                mu = rf_bound(rf_add(opts[i0], o, -1), v0, v1, upper=True)
                if mu is None: return None
                if mu > sl: sl = mu
            slack += sl
    lb = rf_bound(tot, v0, v1)
    return None if lb is None else lb - slack


def rf_mul(x, y):
    r = {}
    for k1, p1 in x.items():
        for k2, p2 in y.items():
            k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
            r[k] = padd(r.get(k, [F(0)]), pmul(p1, p2))
    return r


def _mono(a, b, c):
    r = [F(1)]
    for _ in range(a): r = pmul(r, PS)
    for _ in range(b): r = pmul(r, PC)
    for _ in range(c): r = pmul(r, PN)
    return r


def _trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p


def rf_combine(x):
    """the same rational function over one common denominator S^A C^B N^E (a single key)."""
    ks = [k for k, p in x.items() if any(p)]
    if not ks: return {K0: [F(0)]}
    A = max(k[0] for k in ks); Bc = max(k[1] for k in ks); E = max(k[2] for k in ks)
    P = [F(0)]
    for k in ks: P = padd(P, pmul(x[k], _mono(A - k[0], Bc - k[1], E - k[2])))
    return {(A, Bc, E): P}


def rf_div(x, d):
    """x / d for a rational function d that is a single term kappa S^i C^j N^e / (S^a C^b N^c); None otherwise."""
    d = rf_combine(d)
    ks = [k for k, p in d.items() if any(p)]
    if len(ks) != 1: return None
    (a, b, c) = ks[0]; p = _trim(d[ks[0]])
    for i in range(4):
        for j in range(4):
            for e in range(4):
                mo = _mono(i, j, e)
                if len(mo) != len(p): continue
                kap = p[-1] / mo[-1]
                if kap != 0 and all(p[t] == kap * mo[t] for t in range(len(p))):
                    num = pscale(_mono(a, b, c), 1 / kap)
                    return rf_mulpoly_den(x, num, (i, j, e))
    return None


def rf_is_zero(x):
    return all(not any(p) for p in rf_combine(x).values())


def gc_to_rf(gc):
    al, be, ga = gc
    return ({K0: al}, {K0: be}, {K0: ga})


def aff_at(D, X, Y):
    D0, DX, DY = D
    return rf_add(rf_add(D0, rf_mul(DX, X)), rf_mul(DY, Y))


def region_phi(datas, box, m, cons):
    """Lemma R: a lower bound on the lines' mass over the admissible poses of the box that also satisfy the
    constraints cons = [(D, side)], D = (D0, DX, DY) rational functions of u with D(c, u) = D0 + DX c_x + DY c_y
    (affine in the centre at fixed u), side 'le': D <= 0, 'ge': D >= 0; the constraint lines must be pairwise
    parallel (no line-line vertices).  At fixed u the centres form (a subset of) rect(u) ∩ half-planes, a convex
    polygon whose vertices are among the rectangle corners and the constraint lines ∩ edge lines; the concave Lemma L
    bound is >= its minimum over any finite set whose hull contains the polygon, i.e. over every candidate not PROVEN
    (rf_bound over the sub-bin) to violate a constraint or to leave its edge.  'EMPTY' if every candidate is excluded
    on every sub-bin; None if a bound is unavailable."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    ch = [corner_choices(cx0, True, m, u0, u1), corner_choices(cx1, False, m, u0, u1),
          corner_choices(cy0, True, m, u0, u1), corner_choices(cy1, False, m, u0, u1)]
    pts = sorted(set([u0, u1] + [q for c in ch for (a, b, _) in c for q in (a, b)]))
    best = 'EMPTY'
    for j in range(len(pts) - 1):
        v0, v1 = pts[j], pts[j + 1]
        X0, X1, Y0, Y1 = [next(f for (a, b, f) in c if a <= v0 and v1 <= b) for c in ch]

        def excluded(V, edge=None):
            X, Y = V
            for D, side in cons:
                g = aff_at(D, X, Y)
                if side == 'le':
                    lb = rf_bound(g, v0, v1)
                    if lb is not None and lb > 0: return True
                else:
                    ub = rf_bound(g, v0, v1, upper=True)
                    if ub is not None and ub < 0: return True
            if edge is not None:
                coord, lo, hi = edge
                a_ = rf_bound(rf_add(coord, hi, -1), v0, v1)
                if a_ is not None and a_ > 0: return True
                b_ = rf_bound(rf_add(lo, coord, -1), v0, v1)
                if b_ is not None and b_ > 0: return True
            return False
        cand = []
        for X in (X0, X1):
            for Y in (Y0, Y1):
                if not excluded((X, Y)): cand.append((X, Y))
        for D, side in cons:
            D0, DX, DY = D
            if not rf_is_zero(DX):
                for Y in (Y0, Y1):          # edge cy = Y:  cx = -(D0 + DY Y)/DX
                    X = rf_div(rf_scale(rf_add(D0, rf_mul(DY, Y)), -1), DX)
                    if X is None or any(k[0] and u0 <= 0 for k in X): return None
                    if not excluded((X, Y), (X, X0, X1)): cand.append((X, Y))
            if not rf_is_zero(DY):
                for X in (X0, X1):          # edge cx = X:  cy = -(D0 + DX X)/DY
                    Y = rf_div(rf_scale(rf_add(D0, rf_mul(DX, X)), -1), DY)
                    if Y is None or any(k[0] and u0 <= 0 for k in Y): return None
                    if not excluded((X, Y), (Y, Y0, Y1)): cand.append((X, Y))
        for V in cand:
            v = phi_at(datas, V[0], V[1], v0, v1)
            if v is None: return None
            if best == 'EMPTY' or v < best: best = v
    return best


def vertex_split(datas, box, m):
    """Lemma V: for the Lemma L' (short-chord) line with the largest density, split the box's poses by the sign of
    t_k - t_j (k, j its up / lo inequalities active at the box centre), an affine function of the centre at fixed u.
    Where t_k <= t_j the line is counted as 0; where t_k >= t_j with its Lemma L' term.  Returns the min of the two
    region bounds of the joint Lemma L expression (other short-chord lines dropped, i.e. counted 0), or None."""
    vs = [d for d in datas if d.get('vertex')]
    if not vs: return None
    d = max(vs, key=lambda d: d['rho'])
    others = [e for e in datas if not e.get('vertex')]
    cx0, cx1, cy0, cy1, u0, u1 = box
    if u0 <= 0: return None
    cxm, cym, um = float(cx0 + cx1) / 2, float(cy0 + cy1) / 2, float(u0 + u1) / 2
    k = min(d['up'], key=lambda k: _tk_float(d['fam'], k, d['line'], cxm, cym, um))
    j = max(d['lo'], key=lambda k: _tk_float(d['fam'], k, d['line'], cxm, cym, um))
    Z0 = {K0: [F(0)]}; O1 = {K0: [F(1)]}
    def Dat(X, Y): return rf_add(thr_rf(d['fam'], k, d['line'], X, Y), thr_rf(d['fam'], j, d['line'], X, Y), -1)
    D0 = Dat(Z0, Z0)
    D = (D0, rf_add(Dat(O1, Z0), D0, -1), rf_add(Dat(Z0, O1), D0, -1))
    ra = region_phi(others, box, m, [(D, 'le')]) if others else 'EMPTY0'
    rb = region_phi(others + [d], box, m, [(D, 'ge')])
    if ra is None or rb is None: return None
    if ra == 'EMPTY0':
        ra = region_phi([], box, m, [(D, 'le')])
        if ra is None: return None
    vals = [v for v in (ra, rb) if v != 'EMPTY']
    return min(vals) if vals else None


# =========================================================================================== the checker
class MixedChecker:
    def __init__(self, cover, max_depth=18, use_chain=False, chain_from=0, theta_bias=4, clip=True,
                 use_thr=True, use_pieces=True, dump=False, use_lin=True, use_split=True, cert_mode=False):
        # cert_mode: the reduced trusted surface for certificates (audit S5): Corollary T' (points on germ lines)
        # and the polygon code are unreachable, and a cover containing polygons is refused.
        if cert_mode and cover.polys:
            raise ValueError("--cert-mode: the cover contains polygons (the polygon code is disabled in cert mode)")
        self.cert_mode = cert_mode
        self.cov = cover
        self.m = cover.m
        self.max_depth = max_depth
        self.use_chain = use_chain
        self.chain_from = chain_from
        self.theta_bias = theta_bias
        self.clip = clip
        self.use_thr = use_thr
        self.use_lin = use_lin
        self.use_split = use_split
        self.use_pieces = use_pieces and cover.has_pieces
        self.dump = dump
        pts = [(x, y) for x, y, _ in cover.points]
        ws = [w for _, _, w in cover.points]
        # the PHANTOM (last index): a point far outside, so no geometric test ever accepts it; its weight is set
        # per box to the certified piece bound L and it enters every primitive through `inh` only.
        pts.append((F(-1000), F(-1000))); ws.append(F(1, cover.Wd))
        self.zc = zm.Checker(self.m, pts, ws, use_tri=False, max_depth=max_depth, use_adm=True,
                             theta_bias=theta_bias, use_chain=use_chain, chain_from=chain_from, clip=clip)
        self.ph = len(pts) - 1
        zc = self.zc
        self.Wden = zc.Wden
        order = [int(k) for k in zc.order if int(k) != self.ph]
        zc.order = np.array([self.ph] + order, dtype=np.int64)
        self._set_phantom(F(0))
        self.stat_thr = 0

    def _set_phantom(self, L):
        zc = self.zc; k = self.ph
        if self.Wden:
            num = (L.numerator * self.Wden) // L.denominator        # floor: L rounded DOWN to a multiple of 1/Wden
            Lr = F(num, self.Wden)
            zc.Wnum[k] = num
        else:
            Lr = L
        zc.W[k] = Lr
        zc.Wf[k] = float(Lr)
        return Lr

    # ---------------------------------------------------------------- PIECE bound
    def piece_bound(self, box, B, with_pts=False, moved=None):
        """certified lower bound on (segment + polygon) mass at every admissible pose of the box, and whether the
        threshold lemma raised it above the core bound."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        assert 0 <= u0 <= u1 < 1, 'Lemma T needs theta in [0, 90deg)'
        m = self.m
        specs = self.zc._adm_specs(box, B)
        h = u1 - u0
        cxm, cym = float(cx0 + cx1) / 2, float(cy0 + cy1) / 2
        R = 0.7072 + 0.5 * math.hypot(float(cx1 - cx0), float(cy1 - cy0)) + 1e-9
        cov = self.cov
        self._lparts = None
        ivcache = {}

        def ivs(L):
            k = L['key']
            if k not in ivcache:
                ivcache[k] = [line_cond_iv(L['P0'], L['d'], specs, c, u0, h, m) for c in range(4)]
            return ivcache[k]

        def window(L):
            """reach window of the line's parameter, or None if the line is out of reach (floats; dropping a piece
            only lowers the bound, so any error here is on the safe side for soundness)."""
            px, py, dx, dy, tmin, tmax = L['f']
            dd = dx * dx + dy * dy
            tc = ((cxm - px) * dx + (cym - py) * dy) / dd
            qx, qy = px + tc * dx - cxm, py + tc * dy - cym
            dist2 = qx * qx + qy * qy
            if dist2 > R * R: return None
            half = math.sqrt(max(R * R - dist2, 0.0) / dd) + 1e-9
            lo, hi = tc - half, tc + half
            if hi < float(tmin) or lo > float(tmax): return None
            return lo, hi

        used = set()
        groups = []                 # (keys, Lemma T bound, core bound)
        # ---- Lemma T groups: axis-parallel pairs at distance exactly 1
        if self.use_thr:
            for fam, lines, cm, up_cond, down_cond in (('V', cov.V, cxm, 1, 0), ('H', cov.H, cym, 2, 3)):
                # V: up = x = xi (inequality 1: X >= -1/2 <=> y >= T), down = x = xi + 1 (inequality 0)
                # H: up = y = eta + 1 (inequality 2), down = y = eta (inequality 3)
                cands = set()
                for v in lines:
                    fv = float(v)
                    if fam == 'V':
                        if abs(fv + 0.5 - cm) <= 0.3: cands.add(v)          # v is the up line, xi = v
                        if abs(fv - 0.5 - cm) <= 0.3: cands.add(v - 1)      # v is the down line
                    else:
                        if abs(fv - 0.5 - cm) <= 0.3: cands.add(v - 1)      # v is the up (top) line, eta = v - 1
                        if abs(fv + 0.5 - cm) <= 0.3: cands.add(v)          # v is the down (bottom) line
                for base in sorted(cands):
                    upk, dnk = (base, base + 1) if fam == 'V' else (base + 1, base)
                    Lu = lines.get(upk); Ld = lines.get(dnk)
                    if Lu is not None and (fam, upk) in used: Lu = None
                    if Ld is not None and (fam, dnk) in used: Ld = None
                    wu = window(Lu) if Lu is not None else None
                    wd = window(Ld) if Ld is not None else None
                    if wu is None and wd is None: continue
                    bound, core = self._group(Lu, wu, Ld, wd, ivs, up_cond, down_cond, u0, with_pts, moved)
                    keys = []
                    if Lu is not None: used.add((fam, upk)); keys.append((fam, upk))
                    if Ld is not None: used.add((fam, dnk)); keys.append((fam, dnk))
                    groups.append((keys, bound, core))
        # ---- per-line data: core bound (Lemma S) and, for axis lines, the Lemma L data
        lcore = {}; ldata = {}
        for key, L in cov.lines.items():
            w = window(L)
            if w is None: continue
            I = ivs(L)
            J4 = I[0]
            for c in (1, 2, 3): J4 = iv_and(J4, I[c])
            lcore[key] = sum(((b - a) * dn for a, b, dn in line_pieces(L, J4, w[0], w[1])), F(0))
            if self.use_lin and key[0] in ('V', 'H'):
                dL = lemma_l_data(L, key, I, w, box, self.m)
                if dL is not None: ldata[key] = dL
        gkeys = set(k for keys, _, _ in groups for k in keys)
        others = [k for k in lcore if not (k[0] in ('V', 'H') and (k[0], k[1]) in gkeys)]

        def lin_or_core(keys):
            """(value, Lemma L gained?, parts); parts = (Lemma L datas, core value of those lines, L value) so that
            value = rest + max(core, L) with rest = value - max(core, L)."""
            keys = list(keys)
            el = [k for k in keys if k in ldata]
            base_ = sum((lcore[k] for k in keys if k not in ldata), F(0))
            cs = sum((lcore[k] for k in el), F(0))
            if not el: return base_, False, None
            best = None                                         # (lj, datas, core of those datas)
            lj = lemma_l_joint([ldata[k] for k in el], box, B, m)
            if lj is not None: best = (lj, [ldata[k] for k in el], cs)
            nv = [k for k in el if not ldata[k].get('vertex')]
            if len(nv) < len(el) and nv:  # the short-chord terms can be negative: also try without them
                lj2 = lemma_l_joint([ldata[k] for k in nv], box, B, m)
                if lj2 is not None:
                    csv = sum((lcore[k] for k in el if ldata[k].get('vertex')), F(0))
                    if best is None or lj2 + csv > best[0]:
                        best = (lj2 + csv, [ldata[k] for k in nv], cs - csv)
            if any(ldata[k].get('vertex') for k in el):          # Lemma V: split on the short chord's existence
                vv = vertex_split([ldata[k] for k in el], box, m)
                if vv is not None and (best is None or vv > best[0]):
                    best = (vv, [ldata[k] for k in el], cs)
            if best is None: return base_ + cs, False, None
            ljv, datas, csd = best
            # value of the part carried by `datas`: ljv - (vertex cores added) ; keep it simple:
            part_L = ljv - (cs - csd)          # = Lemma L joint value of `datas`
            part = max(part_L, csd)
            val = base_ + (cs - csd) + part
            return val, part_L > csd, (datas, csd, part_L)

        A = sum((max(bd, cr) for _, bd, cr in groups), F(0))
        tgain = any(bd > cr for _, bd, cr in groups)
        oth, lgain, parts = lin_or_core(others)
        A += oth
        total = A
        if groups and ldata:
            Bv, lg2, parts2 = lin_or_core([k for keys, _, _ in groups for k in keys if k in lcore] + others)
            if Bv > total: total = Bv; tgain = False; lgain = lg2; parts = parts2
        # for the region-wise bound (Lemma R): total = rest + max(core_part, L_part), rest excludes the datas' lines
        self._lparts = None if parts is None else (parts[0], parts[1], total - max(parts[1], parts[2]))
        thr_gain = 'T' if tgain else ('L' if lgain else False)
        # ---- polygons (Lemma S, 2-D)
        if cov.polys:
            assert not self.cert_mode, "polygon code reached in cert mode"
            fx0, fx1, fy0, fy1 = cxm - R, cxm + R, cym - R, cym + R
            near = [P for P in cov.polys if not (P['bb'][1] < fx0 or P['bb'][0] > fx1 or
                                                  P['bb'][3] < fy0 or P['bb'][2] > fy1)]
            if near:
                r = REACH
                frame = [(cx0 - r, cy0 - r), (cx1 + r, cy0 - r), (cx1 + r, cy1 + r), (cx0 - r, cy1 + r)]
                K = frame
                for c in range(4):
                    Kc = cond_region(specs, c, u0, h, m, (cx0, cy0), frame)
                    K = clip_convex(K, Kc) if Kc else []
                    if not K: break
                if K:
                    for P in near:
                        Vi = clip_convex(P['V'], K)
                        if Vi:
                            pa = P['w'] * poly_area(Vi) / P['area']
                            total += pa
                            if self._lparts is not None:
                                self._lparts = (self._lparts[0], self._lparts[1], self._lparts[2] + pa)
        return total, thr_gain

    def _group(self, Lu, wu, Ld, wd, ivs, up_cond, down_cond, u0, with_pts=False, moved=None):
        """Lemma T: inf over T of  F_up(min(T, t_up)) + G_down(max(T + u0, t_dn))  [+ the point steps when
        with_pts: a point of the up line at parameter t_p in J3 counts iff min(T, t_up) <= t_p, one of the down
        line iff t_p <= max(T + u0, t_dn); such points are appended to `moved` and must then be withheld from the
        point primitives].  Also returns the plain core bound of the segments."""
        terms = []          # (kind, LineMass, tstar, points [(t_p, w)])
        core = F(0)
        for L, wnd, cond, kind in ((Lu, wu, up_cond, 'up'), (Ld, wd, down_cond, 'down')):
            if L is None or wnd is None: continue
            I = ivs(L)
            J3 = (None, None)
            for c in range(4):
                if c != cond: J3 = iv_and(J3, I[c])
            sing = I[cond]
            if kind == 'up':    # certified singular set must be an up-set [t*, inf); else ignored (t* = +inf)
                tstar = None if (sing is None or sing[1] is not None) else sing[0]
            else:               # a down-set (-inf, t*]; else ignored (t* = -inf)
                tstar = None if (sing is None or sing[0] is not None) else sing[1]
            if sing is not None and sing[0] is None and sing[1] is None: tstar = 'all'
            LM = LineMass(line_pieces(L, J3, wnd[0], wnd[1]))
            pts = []
            if with_pts and J3 is not None:
                for tp, w, idx in self.cov.line_points.get(L['key'], ()):
                    if idx in moved: continue          # a grid-crossing point is counted by ONE group only
                    if (J3[0] is None or tp >= J3[0]) and (J3[1] is None or tp <= J3[1]):
                        pts.append((tp, w)); moved.append(idx)
            terms.append((kind, LM, tstar, pts))
            J4 = iv_and(J3, sing)
            core += sum((b - a) * dn for a, b, dn in line_pieces(L, J4, wnd[0], wnd[1]))
        # f = PL(T) + STEP(T): PL continuous piecewise linear, STEP a sum of indicator steps; both change only at
        # the breakpoints, so inf f = min over breakpoints b of f(b) and over the open intervals between
        # consecutive breakpoints (and the two tails) of STEP(midpoint) + min(PL at the two ends).
        cand = set()
        for kind, LM, ts, pts in terms:
            if ts == 'all': continue
            for e in LM.bps + [tp for tp, _ in pts]:
                if kind == 'up' and (ts is None or e <= ts): cand.add(e)
                if kind == 'down' and (ts is None or e >= ts): cand.add(e - u0)
            if ts is not None: cand.add(ts if kind == 'up' else ts - u0)
        cand = sorted(cand) or [F(0)]

        def PL(T):
            v = F(0)
            for kind, LM, ts, pts in terms:
                if ts == 'all': v += LM.tot; continue
                if kind == 'up': v += LM.above(T if ts is None else min(T, ts))
                else: v += LM.below(T + u0 if ts is None else max(T + u0, ts))
            return v

        def STEP(T):
            v = F(0)
            for kind, LM, ts, pts in terms:
                if not pts: continue
                if ts == 'all': v += sum(w for _, w in pts); continue
                if kind == 'up':
                    x = T if ts is None else min(T, ts)
                    v += sum((w for tp, w in pts if x <= tp), F(0))
                else:
                    y = T + u0 if ts is None else max(T + u0, ts)
                    v += sum((w for tp, w in pts if tp <= y), F(0))
            return v
        pl = [PL(T) for T in cand]
        best = min(pl[i] + STEP(T) for i, T in enumerate(cand))
        mids = [(cand[0] - 1, pl[0], pl[0]), (cand[-1] + 1, pl[-1], pl[-1])]
        mids += [((cand[i] + cand[i + 1]) / 2, pl[i], pl[i + 1]) for i in range(len(cand) - 1)]
        for T, pa, pb in mids:
            v = STEP(T) + min(pa, pb)
            if v < best: best = v
        return best, core

    # ---------------------------------------------------------------- SPLIT: region-wise pieces + points (Lemma R)
    def cert_split(self, box, B, inh, Lbox, lparts, maxchain=400, diag=None):
        """Disjunctive certification coupling pieces and points region by region (ZM_MIXED.md Lemma R).
        T = points in Q at every admissible pose (ADM / P1 / inherited).  For each swing kind: a chain of swing
        points q_1..q_k with G_{q_1} <= ... <= G_{q_k} on the box (exact); every admissible pose lies in a region
        r = 0..k:  G_{q_r} <= 0 (r >= 1) and G_{q_{r+1}} > 0 (r < k).  In region r the points of T, of the
        down-set (G_p <= G_{q_r} on the box) and of the up-set (G_p + lam G_{q_{r+1}} <= 0 on the box) are in Q,
        and the pieces carry >= rest + max(core, region_phi(...)) (Lemma R).  Certified if every region reaches 1."""
        if lparts is None: lparts = ([], F(0), Lbox)        # no Lemma L lines: the pieces are the constant Lbox
        datas, core_part, rest = lparts
        zc = self.zc
        cx0, cx1, cy0, cy1, u0, u1 = box
        if u0 <= 0: return (None, None)
        specs = zc._adm_specs(box, B)
        ma, cmasks = zc._adm_mask(specs, u0, u1, per_cond=True)
        mp = zc._p1_mask(box, B)
        t = 1 - B['whi'] / 2
        n = len(zc.P) - 1                                       # the phantom (last) is not a point here
        inT = np.zeros(n, dtype=bool); wT = F(0)
        sel = ma[:n] | mp[:n]
        if inh is not None: sel = sel | inh[:n]
        for k in np.nonzero(sel)[0]:
            px, py = zc.P[k]
            if (inh is not None and inh[k]) or (mp[k] and zc._p1_exact(px, py, box, t)) or \
               (ma[k] and zc._adm_exact(px, py, specs, u0, u1)):
                inT[k] = True; wT += zc.W[k]
        base = wT + rest
        if base + core_part >= 1 and diag is None: return ('SPLIT', 'T')
        cxm, cym = float((cx0 + cx1) / 2), float((cy0 + cy1) / 2)
        rad = 0.7072 + 0.5 * math.hypot(float(cx1 - cx0), float(cy1 - cy0))
        reach = (((zc.Pxf[:n] - cxm) ** 2 + (zc.Pyf[:n] - cym) ** 2) <= rad * rad) & (~inT) & (zc.Wf[:n] > 0)
        nfail = np.zeros(n, dtype=np.int8)
        for cm in cmasks: nfail += (~cm[:n])
        reach &= (nfail <= 1)
        cand = []
        for k in np.nonzero(reach)[0]:
            px, py = zc.P[k]
            bad = None
            for c in range(4):
                if not cmasks[c][k]:
                    if bad is not None: bad = -1; break
                    bad = c; continue
                if not zc._adm_cond_ok(px, py, specs, c, u0, u1):
                    if bad is not None: bad = -1; break
                    bad = c
            if bad is None or bad < 0: continue
            cand.append((int(k), bad))
        if not cand: return (None, None)
        ctx = zc._box_ctx(box)
        # zeromargin's _gle0 docstring speaks of positive weights; here it is called with weights +1 and -1.  Its
        # corner argument (the G's are affine in the centre at fixed u, so the max over the centre box is at a
        # corner) holds for any real weights (audit nit N4).
        def le(terms): return zc._gle0(terms, ctx)
        m = self.m
        kinds = sorted({kd for _, kd in cand}, key=lambda kd: -sum(zc.Wf[k] for k, d in cand if d == kd))
        for kd in kinds:
            grp = [ck for ck in cand if ck[1] == kd]
            def proxy(ck):
                g = zc._gcoef(kd, zc.P[ck[0]][0] - F(cxm), zc.P[ck[0]][1] - F(cym))
                return float(g[0] + g[1] * (u0 + u1) / 2 + g[2] * ((u0 + u1) / 2) ** 2)
            grp.sort(key=proxy)
            ch = []
            for ck in grp:
                if not ch or le([(1, ch[-1][0], kd), (-1, ck[0], kd)]): ch.append(ck)
                if len(ch) >= maxchain: break
            kk = len(ch)
            # down[r] / up[r] membership (exact integer tests, binary searches as in zeromargin's CHAIN)
            W = {k: zc.W[k] for k, _ in cand}
            dlo = {}; ulo = {}
            for (k, kind) in cand:
                lo, hi = 1, kk + 1
                while lo < hi:
                    mid = (lo + hi) // 2
                    q = ch[mid - 1]
                    if q[0] == k or le([(1, k, kind), (-1, q[0], q[1])]): hi = mid
                    else: lo = mid + 1
                dlo[k] = lo
                lo, hi = 0, kk
                while lo < hi:
                    mid = (lo + hi + 1) // 2
                    q = ch[mid - 1]
                    good = q[0] != k and any(le([(a, k, kind), (b, q[0], q[1])]) for a, b in ((1, 1), (2, 1), (1, 2)))
                    if good: lo = mid
                    else: hi = mid - 1
                ulo[k] = lo
            gcs = [gc_to_rf(g_coeffs(kd, *zc.P[q[0]])) for q in ch]
            ok = True
            regs = []
            for r in range(kk + 1):
                Wr = sum((W[k] for k, _ in cand if dlo[k] <= r or ulo[k] >= r + 1), F(0))
                cons = []
                if r >= 1: cons.append((gcs[r - 1], 'le'))
                if r < kk: cons.append((gcs[r], 'ge'))
                if diag is None and base + Wr + (Lbox - rest) >= 1: continue
                ph = region_phi(datas, box, m, cons)
                pb = rest + (core_part if ph in (None, 'EMPTY') else max(core_part, ph))
                if diag is not None:
                    regs.append(dict(r=r, pts=[int(k) for k in np.nonzero(inT)[0]] +
                                     [k for k, _ in cand if dlo[k] <= r or ulo[k] >= r + 1], pieces=pb, empty=(ph == 'EMPTY')))
                if ph == 'EMPTY': continue
                if wT + Wr + pb < 1:
                    ok = False
                    if diag is None: break
            if diag is not None:
                diag.append(dict(kind=kd, chain=[q[0] for q in ch], regions=regs, ok=ok))
            if ok and diag is None: return ('SPLIT', ('kind', kd, 'k', kk))
        return (None, None)

    # ---------------------------------------------------------------- recursion (mirrors zeromargin.run_box)
    def run_box(self, root):
        zc = self.zc
        stats = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'PIECE': 0, 'EMPTY': 0,
                 'UNCERT': 0, 'boxes': 0, 'maxdepth': 0, 'THR': 0, 'LIN': 0, 'TPTS': 0, 'SPLIT': 0, 'cpu': 0.0}
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
                Bf = tuple(float(B[k]) for k in ('c0', 's0', 'c1', 's1', 'cD', 'sD'))
                bincache[key] = (B, Bf)
            B, Bf = bincache[key]
            lo = B['wlo'] / 2
            if cx1 < lo or cx0 > self.m - lo or cy1 < lo or cy0 > self.m - lo:
                stats['EMPTY'] += 1
                if self.dump: leaves.append((box, 'EMPTY', None))
                continue
            L = F(0); thr = False
            lparts = None
            if self.use_pieces:
                L, thr = self.piece_bound(box, B)
                lparts = self._lparts
                if Lpar > L: L = Lpar                    # the parent's bound holds on this sub-box too
            kind = wit = None
            if L >= 1:
                kind, wit = 'PIECE', str(L)
            else:
                Lr = self._set_phantom(L)
                inh2 = base.copy() if inh is None else (inh | base)
                zc._last_inT = None
                kind, wit = zc.cert_adm(box, B, inh2)
                if kind is None:
                    kind, wit = zc.cert_p1(box, B)
                if kind is None:
                    kind, wit = zc.cert_mix(box, B, inh2)
                if kind is None and self.use_chain and depth >= self.chain_from:
                    kind, wit = zc.cert_chain(box, B, inh=inh2)
                if kind is None and self.use_thr and self.cov.line_points and not self.cert_mode:
                    # second attempt (Lemma T with points): the points on the germ-pair lines are accounted for
                    # inside the groups and withheld (weight 0) from the point primitives
                    moved = []
                    L2, _ = self.piece_bound(box, B, with_pts=True, moved=moved)
                    if moved:
                        if Lpar > L2: L2 = Lpar
                        if L2 >= 1:
                            kind, wit = 'PIECE', str(L2)
                        else:
                            save = [(k, zc.W[k], zc.Wf[k], zc.Wnum[k]) for k in moved]
                            try:
                                for k in moved:
                                    zc.W[k] = F(0); zc.Wf[k] = 0.0; zc.Wnum[k] = 0
                                Lr = self._set_phantom(L2)
                                zc._last_inT = None
                                kind, wit = zc.cert_adm(box, B, inh2)
                                if kind is None:
                                    kind, wit = zc.cert_p1(box, B)
                                if kind is None:
                                    kind, wit = zc.cert_mix(box, B, inh2)
                                if kind is None and self.use_chain and depth >= self.chain_from:
                                    kind, wit = zc.cert_chain(box, B, inh=inh2)
                            finally:
                                for k, a1, a2, a3 in save:
                                    zc.W[k] = a1; zc.Wf[k] = a2; zc.Wnum[k] = a3
                        if kind is not None:
                            stats['TPTS'] += 1
                            if self.dump: wit = ('moved', sorted(moved), wit)
                if kind is None and self.use_split and lparts is not None:
                    kind, wit = self.cert_split(box, B, inh, L, lparts)
                if kind is not None and self.dump:
                    wit = (str(Lr), wit)
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


# =========================================================================================== exact / float pose mass
def exact_mass(cov, cx, cy, u):
    """exact mu(Q(c, theta)) at a rational pose, theta = 2 atan u (Fractions; closed square)."""
    c, s = zm.trig(u)
    tot = F(0)
    for px, py, w in cov.points:
        if zm.in_rot_square(px - cx, py - cy, c, s): tot += w
    for L in cov.lines.values():
        (P0x, P0y), (dx, dy) = L['P0'], L['d']
        # X(t) = (P0x + t dx - cx) c + (P0y + t dy - cy) s,  Y(t) = -(...) s + (...) c ; |X|, |Y| <= 1/2
        iv = (None, None)
        for a0, a1 in (((P0x - cx) * c + (P0y - cy) * s, dx * c + dy * s),
                       (-(P0x - cx) * s + (P0y - cy) * c, -dx * s + dy * c)):
            for sg in (1, -1):                       # sg (a0 + t a1) <= 1/2
                b0, b1 = sg * a0 - HALF, sg * a1
                if b1 > 0: iv = iv_and(iv, (None, -b0 / b1))
                elif b1 < 0: iv = iv_and(iv, (-b0 / b1, None))
                elif b0 > 0: iv = None
                if iv is None: break
            if iv is None: break
        for t0, t1, dn in L['segs']:
            tot += seg_len_in(t0, t1, iv) * dn
    if cov.polys:
        Q = [(cx + (sx * c - sy * s) / 2, cy + (sx * s + sy * c) / 2)
             for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
        for P in cov.polys:
            Vi = clip_convex(P['V'], Q)
            if Vi: tot += P['w'] * poly_area(Vi) / P['area']
    return tot


def admissible(m, cx, cy, u):
    c, s = zm.trig(u)
    w = abs(c) + abs(s)
    return w / 2 <= cx <= m - w / 2 and w / 2 <= cy <= m - w / 2


# =========================================================================================== driver
_MC = None
def _worker1(root):
    return root, _MC.run_box(root)


def run_sweep(chk, R, nproc, chunksize, progress, label='', resume=None, done=None):
    """resume: a jsonl file getting one line per finished root (root, census, uncertified boxes); roots already
    in `done` (read from it) are skipped and their census added."""
    import json
    global _MC
    _MC = chk
    prev = []
    if done:
        prev = [done[k] for k in [(label, tuple(str(v) for v in r)) for r in R] if k in done]
        R = [r for r in R if (label, tuple(str(v) for v in r)) not in done]
    fres = open(resume, 'a') if resume else None
    tot = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'PIECE': 0, 'EMPTY': 0, 'UNCERT': 0,
           'boxes': 0, 'maxdepth': 0, 'THR': 0, 'LIN': 0, 'TPTS': 0, 'SPLIT': 0, 'cpu': 0.0}
    unc_all = []; leaves_all = []
    t0 = time.time()
    import multiprocessing as _mp
    if nproc <= 1:
        it = map(_worker1, R)
        pool = None
    else:
        pool = _mp.get_context('fork').Pool(nproc)
        it = pool.imap_unordered(_worker1, R, chunksize=chunksize)
    for st, unc in prev:
        for k in tot:
            tot[k] = max(tot[k], st.get(k, 0)) if k == 'maxdepth' else tot[k] + st.get(k, 0)
        unc_all += unc
    t1 = time.time()
    for i, (root, (st, unc, leaves)) in enumerate(it):
        for k in tot:
            tot[k] = max(tot[k], st[k]) if k == 'maxdepth' else tot[k] + st[k]
        unc_all += unc; leaves_all += leaves
        if fres:
            fres.write(json.dumps(dict(label=label, root=[str(v) for v in root], st=st,
                                       unc=[[str(v) for v in b] for b in unc], cpu=None)) + "\n")
            fres.flush()
        if (i + 1) % progress == 0:
            print(f"  {label}{i+1}/{len(R)} roots, {tot['boxes']} boxes, uncert {tot['UNCERT']}, "
                  f"{time.time()-t0:.0f}s", flush=True)
    if pool is not None:
        pool.close(); pool.join()
    if fres: fres.close()
    return tot, unc_all, leaves_all, time.time() - t0


def file_sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def run_header(a, cov, input_path):
    """Provenance of a cert run (audit S2/S3): the sha256 of every file the verdict depends on and every setting
    that changes what is checked.  Written as the first line of the --resume jsonl and into the --manifest;
    --resume refuses a file whose header differs.  (--nproc, --chunksize, --progress, --dump, --resume,
    --manifest do not change the verdict and are not part of it.)"""
    mode = 'D4' if a.d4 else ('full' if a.full else 'D2')
    return dict(
        kind='zm_mixed cert header', version=2,
        sha256={'zm_mixed.py': file_sha(os.path.abspath(__file__)), 'mixed_cover.py': file_sha(MC.__file__),
                'zeromargin.py': zm_sha(), 'input': file_sha(input_path)},
        input=input_path,
        total=str(cov.total),
        settings=dict(mode=mode, depth=a.depth, pitch=str(F(a.pitch)), ubins=a.ubins,
                      chain=bool(a.disj), chain_from=a.chain_from, theta_bias=a.theta_bias, clip=not a.no_clip,
                      lemma_T=not a.no_thr, lemma_L=not a.no_lin, split=not a.no_split,
                      cert_mode=bool(a.cert_mode), tprime=(not a.no_thr) and not a.cert_mode,
                      split_maxchain=400, reach=str(REACH),
                      region=dict(cx_lo=a.cx_lo, cx_hi=a.cx_hi, cy_lo=a.cy_lo, cy_hi=a.cy_hi,
                                  u_lo=a.u_lo, u_hi=a.u_hi)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what', help='cert | pose | stress | selftest')
    ap.add_argument('path', nargs='?')
    ap.add_argument('leafdump', nargs='?')
    ap.add_argument('--disj', action='store_true', help='enable zeromargin CHAIN for the points')
    ap.add_argument('--chain-from', type=int, default=0)
    ap.add_argument('--theta-bias', type=int, default=4)
    ap.add_argument('--depth', type=int, default=18)
    ap.add_argument('--nproc', type=int, default=5)
    ap.add_argument('--pitch', type=str, default='1/10')
    ap.add_argument('--ubins', type=int, default=8)
    ap.add_argument('--chunksize', type=int, default=1)
    ap.add_argument('--progress', type=int, default=100)
    ap.add_argument('--dump', type=str, default=None)
    ap.add_argument('--resume', type=str, default=None,
                    help='jsonl file: one line per finished root; on restart the roots in it are skipped (their '
                         'census and uncertified boxes are counted)')
    ap.add_argument('--full', action='store_true',
                    help='no symmetry assumed: the cover and its (x,y)->(y,x) image, all centres, u in [0,1/2]')
    ap.add_argument('--d4', action='store_true')
    ap.add_argument('--cert-mode', action='store_true',
                    help='certificate mode: Corollary T\' and the polygon code disabled; covers with polygons refused')
    ap.add_argument('--manifest', type=str, default=None,
                    help='write a JSON manifest (header: shas + settings; result: census + verdict) at the end')
    ap.add_argument('--no-thr', action='store_true', help='disable Lemma T (for comparison)')
    ap.add_argument('--no-lin', action='store_true', help='disable Lemma L (for comparison)')
    ap.add_argument('--no-split', action='store_true', help='disable the region-wise SPLIT primitive (Lemma R)')
    ap.add_argument('--no-clip', action='store_true')
    ap.add_argument('--cx-lo', type=str, default=None); ap.add_argument('--cx-hi', type=str, default=None)
    ap.add_argument('--cy-lo', type=str, default=None); ap.add_argument('--cy-hi', type=str, default=None)
    ap.add_argument('--u-lo', type=str, default=None, help='keep only root u-bins with u1 > this (PARTIAL)')
    ap.add_argument('--u-hi', type=str, default=None, help='keep only root u-bins with u0 < this (PARTIAL)')
    ap.add_argument('--cx', type=str); ap.add_argument('--cy', type=str); ap.add_argument('--u', type=str)
    ap.add_argument('--per', type=int, default=20)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()

    sha = zm_sha()
    if a.what == 'selftest':
        import zm_mixed_test
        sys.exit(zm_mixed_test.selftest(a.per * 15, a.seed))
    cv = MC.load(a.path)
    cov = Cover(cv)
    if a.what == 'pose':
        cx, cy, u = F(a.cx), F(a.cy), F(a.u)
        v = exact_mass(cov, cx, cy, u)
        print(f"pose cx={cx} cy={cy} u={u} (theta {math.degrees(2*math.atan(float(u))):.9f} deg), admissible "
              f"{admissible(cov.m, cx, cy, u)}: EXACT mu(Q) = {v} = {float(v):.12f} "
              f"({'OK' if v >= 1 else '*** < 1 ***'})")
        return
    if a.what == 'stress':
        import zm_mixed_test
        sys.exit(zm_mixed_test.stress(cov, a.leafdump, a.per, a.seed))
    assert a.what == 'cert'
    import json
    print(f"zeromargin.py sha256 {sha} ({'= pinned' if sha == ZM_SHA_PINNED else '*** DIFFERS FROM PINNED ***'})")
    hdr = run_header(a, cov, a.path)
    for k, v in hdr['sha256'].items():
        print(f"sha256 {v}  {k}" + (f" ({a.path})" if k == 'input' else ''))
    print('argv:', ' '.join(sys.argv))
    print('settings:', json.dumps(hdr['settings'], sort_keys=True), flush=True)
    if a.cert_mode:
        if cov.polys:
            print("ERROR: --cert-mode: the cover contains polygons; refused"); sys.exit(2)
        print(f"CERT MODE: Corollary T' and polygon code disabled ({sum(len(v) for v in cov.line_points.values())} "
              f"points on segment lines are treated as ordinary points)")
    print(f"container [0,{cov.m}]^2, {len(cov.points)} points, {sum(len(L['segs']) for L in cov.lines.values())} "
          f"segments on {len(cov.lines)} lines, {len(cov.polys)} polygons; total {cov.total} = {float(cov.total):.9f}")
    partial = any(v is not None for v in (a.cx_lo, a.cx_hi, a.cy_lo, a.cy_hi, a.u_lo, a.u_hi))
    pitch = F(a.pitch)
    covers = [(cov, '')]
    if a.d4:
        if not cov.symmetric_d4():
            print("ERROR: --d4: the measure is not D4-invariant (canonical multisets differ)"); sys.exit(2)
        print("D4: measure invariant under x -> m-x and x <-> y; roots [0,m/2]^2 x u in [0,1/2]")
        R0 = zm.d4_roots(cov.m, pitch, a.ubins)
    elif a.full:
        print("FULL: no symmetry; the cover on u in [0,1/2] and its (x,y)->(y,x) image on u in [0,1/2]")
        R0 = zm.roots(cov.m, pitch, a.ubins, cy_max=cov.m)
        covers.append((Cover(swapped(cv)), 'swap:'))
    else:
        if not cov.symmetric_d2():
            print("measure not invariant under x -> m-x and y -> m-y: use --full"); sys.exit(2)
        print("D2: measure invariant under x -> m-x and y -> m-y; roots cy <= m/2, u in [0,1/2]")
        R0 = zm.roots(cov.m, pitch, a.ubins)

    def keep(r):
        x0, x1, y0, y1, u0, u1 = r
        if a.cx_lo is not None and x1 <= F(a.cx_lo): return False
        if a.cx_hi is not None and x0 >= F(a.cx_hi): return False
        if a.cy_lo is not None and y1 <= F(a.cy_lo): return False
        if a.cy_hi is not None and y0 >= F(a.cy_hi): return False
        if a.u_lo is not None and u1 <= F(a.u_lo): return False
        if a.u_hi is not None and u0 >= F(a.u_hi): return False
        return True
    R0 = [r for r in R0 if keep(r)]
    if partial:
        print(f"PARTIAL SWEEP: cx [{a.cx_lo}, {a.cx_hi}], cy [{a.cy_lo}, {a.cy_hi}], u [{a.u_lo}, {a.u_hi}] -- "
              f"not a verification of the whole container")
    print(f"{len(R0)} root boxes x {len(covers)} cover(s), depth limit {a.depth}, pitch {pitch}, u-bins {a.ubins}, "
          f"CHAIN {'on' if a.disj else 'off'}, Lemma T {'off' if a.no_thr else 'on'}, Lemma L {'off' if a.no_lin else 'on'}", flush=True)
    grand = None; unc_all = []; leaves_all = []; wall = 0
    for cc, label in covers:
        chk = MixedChecker(cc, max_depth=a.depth, use_chain=a.disj, chain_from=a.chain_from,
                           theta_bias=a.theta_bias, clip=not a.no_clip, use_thr=not a.no_thr, use_lin=not a.no_lin, use_split=not a.no_split,
                           dump=a.dump is not None, cert_mode=a.cert_mode)
        done = {}
        if a.resume and os.path.exists(a.resume) and os.path.getsize(a.resume) > 0:
            with open(a.resume) as fh:
                lines = [ln for ln in fh if ln.strip()]
            h0 = json.loads(lines[0])
            if h0.get('kind') != hdr['kind']:
                print(f"ERROR: --resume {a.resume}: no header line (written by an older zm_mixed.py?); refused")
                sys.exit(2)
            if {k: h0.get(k) for k in ('sha256', 'settings', 'total')} != \
                    {k: hdr[k] for k in ('sha256', 'settings', 'total')}:
                print(f"ERROR: --resume {a.resume}: its header (shas / settings) differs from this run; refused")
                for k in ('sha256', 'settings', 'total'):
                    if h0.get(k) != hdr[k]: print(f"  {k}: file {json.dumps(h0.get(k))}\n  {k}: now  {json.dumps(hdr[k])}")
                sys.exit(2)
            for ln in lines[1:]:
                d = json.loads(ln)
                done[(d['label'], tuple(d['root']))] = (d['st'], [tuple(F(v) for v in b) for b in d['unc']])
            print(f"resume: {len(done)} finished roots read from {a.resume} (header matches)", flush=True)
        elif a.resume and label == '':
            with open(a.resume, 'w') as fh:
                fh.write(json.dumps(hdr) + "\n")
        tot, unc, leaves, dt = run_sweep(chk, R0, a.nproc, a.chunksize, a.progress, label, a.resume, done)
        wall += dt
        unc_all += [(label, b) for b in unc]; leaves_all += [(label, l) for l in leaves]
        if grand is None: grand = tot
        else:
            for k in grand: grand[k] = max(grand[k], tot[k]) if k == 'maxdepth' else grand[k] + tot[k]
    tot = grand
    print(f"done in {wall:.0f}s: boxes {tot['boxes']}, max depth {tot['maxdepth']}, CPU {tot['cpu']:.0f} s")
    print(f"  leaves: PIECE {tot['PIECE']}  ADM {tot['ADM']}  P1 {tot['P1']}  MIX {tot['MIX']}  CHAIN {tot['CHAIN']}  SPLIT {tot['SPLIT']}  "
          f"EMPTY {tot['EMPTY']}  UNCERTIFIED {tot['UNCERT']}   (leaves where Lemma T / Lemma L raised the piece bound: {tot['THR']} / {tot['LIN']}; closed by Lemma T with line points: {tot['TPTS']})")
    if unc_all:
        print("uncertified boxes (cover cx0 cx1 cy0 cy1 theta0 theta1 deg):")
        for lab, b in sorted(unc_all, key=lambda t: t[1])[:40]:
            print("  ", lab or '-', *[f"{float(v):.6f}" for v in b[:4]],
                  f"{math.degrees(2*math.atan(float(b[4]))):.4f} {math.degrees(2*math.atan(float(b[5]))):.4f}")
        if len(unc_all) > 40: print(f"  ... {len(unc_all)} in total")
    if a.dump:
        with open(a.dump, 'w') as f:
            f.write(f"# zm_mixed leaves; container {cov.m}; cover ; cx0 cx1 cy0 cy1 u0 u1 ; kind ; witness\n")
            for lab, (box, kind, wit) in leaves_all:
                f.write((lab or '-') + " ; " + " ".join(str(v) for v in box) + f" ; {kind} ; {wit}\n")
            for lab, box in unc_all:
                f.write((lab or '-') + " ; " + " ".join(str(v) for v in box) + " ; UNCERT ; -\n")
        print(f"leaves written to {a.dump}")
    tag = "VERIFIED" + ("-D4" if a.d4 else "") if tot['UNCERT'] == 0 else "NOT VERIFIED"
    print(tag + (" (PARTIAL)" if partial else ""))
    if a.manifest:
        man = dict(header=hdr, argv=sys.argv, created=time.strftime('%Y-%m-%d %H:%M:%S'),
                   result=dict(roots=len(R0) * len(covers), census=tot, uncertified=len(unc_all), wall_s=round(wall, 1),
                               verdict=tag + (" (PARTIAL)" if partial else "")))
        if a.resume:
            man['records'] = dict(path=a.resume, sha256=file_sha(a.resume))
        with open(a.manifest, 'w') as fh:
            json.dump(man, fh, indent=1, sort_keys=False); fh.write("\n")
        print(f"manifest written to {a.manifest}")


if __name__ == '__main__':
    main()
