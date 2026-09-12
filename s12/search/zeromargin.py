#!/usr/bin/env python3
"""Exact zero-margin checker for closed covers at the container [0,m]^2 (task E, rung 1).

Statement checked (unit weights):  every closed unit square Q, at every angle, contained in the
closed square [0,m]^2, contains at least one point of the set P.

Method: adaptive subdivision of pose space (cx, cy, u = tan(theta/2)) into boxes with rational
endpoints; a box is certified by one of

  CORE  the point p lies in the exact core of the box: for every centre c of the box's centre
        rectangle and every angle in the box's angle bin, p is in the closed unit square (c, theta).
        The core over an angle bin [t0, t1] at a fixed centre is the intersection of the rotated
        squares R_t Q, which equals  R_t0 Q  ∩  R_t1 Q  ∩  {x : dir(x) mod 90deg in [t0,t1] => |x| <= 1/2}
        (support-function argument; see search/ZEROMARGIN.md).  It is convex, so testing the four
        corners of the rectangle p - rect suffices.  All arithmetic is in Fractions.
  P1    the "2x2 box lemma": if |c - p|_inf <= 1 - w(theta)/2 with w = |cos| + |sin| (i.e. the
        square lies in p + [-1,1]^2), then p in S.  Proof: in the rotated frame the coordinate of
        p - c is at most (1 - w/2) w = w - w^2/2 <= 1/2 since (w-1)^2 >= 0.  Applied to the
        admissible poses of the box only (a wall makes the constraint on that side automatic).
        Exact at the wall/corner tight poses, where CORE never terminates (quadratic margin).
  TRI   (optional) the triangle lemma: a triangle with vertices in P and sides <= 1 containing the
        box's centre rectangle certifies it for every angle (Friedman/Stromquist).
  EMPTY the box contains no admissible pose.

Only admissible poses matter (cx, cy in [w/2, m - w/2]); the centre rectangle is clipped to the
admissible range at the bin's smallest w before CORE, which over-tests inadmissible poses in the
sliver (sound), and P1 handles the slivers at the tight wall poses.

Symmetry: if P is invariant under x -> m-x and y -> m-y (checked exactly), the fundamental domain
is theta in [0, 45deg] (x-reflection sends theta to -theta) and cy <= m/2 (the 180deg rotation).
The root domain uses u in [0, 1/2] (theta up to 53deg) to keep endpoints rational.

Usage:  python3 search/zeromargin.py friedman14 [--tri] [--depth 16] [--nproc 4] [--dump FILE]
        python3 search/zeromargin.py cert certificates/xxx.txt ...   (weighted: sum of certified weights >= 1)
"""
import sys, os, time, argparse, itertools, math
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
from fractions import Fraction as F
from multiprocessing import Pool
import numpy as np

HALF = F(1, 2)

def trig(u):
    """cos, sin of theta = 2 atan u, exact."""
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d

def bin_data(u0, u1):
    c0, s0 = trig(u0); c1, s1 = trig(u1)
    cD = c1 * c0 + s1 * s0          # cos(t1 - t0)
    sD = s1 * c0 - c1 * s0          # sin(t1 - t0)  (>= 0)
    w0, w1 = c0 + s0, c1 + s1
    # max of w = cos+sin over the bin: at 45deg if the bin contains it (u = sqrt2-1: u^2+2u-1 = 0)
    if u0 * u0 + 2 * u0 - 1 < 0 < u1 * u1 + 2 * u1 - 1:
        whi = F(14143, 10000)       # > sqrt 2
    else:
        whi = max(w0, w1)
    wlo = min(w0, w1)
    return dict(c0=c0, s0=s0, c1=c1, s1=s1, cD=cD, sD=sD, whi=whi, wlo=wlo)

def clip_bin(box, m, steps=28):
    """Shrink a box's angle bin to a sub-bin that still contains every ADMISSIBLE pose of the box.

    A pose (c, theta) of the box is admissible iff  w(theta)/2 <= c_x <= m - w(theta)/2  and the
    same in y, i.e. iff  w(theta) <= K  with

        K = 2 min(cx1, m - cx0, cy1, m - cy0).

    On [0, 45deg] w is strictly increasing, so {u : w(u) <= K} = [0, u-] and every pose of the box
    with u > u- is inadmissible and constrains nothing.  Returns a rational u* >= u- found by
    bisection (so no admissible pose is ever dropped -- the clip is sound in the direction that
    matters), or u1 unchanged when the bin is not inside [0, 45deg], where w is not monotone.

    Without this, a box pressed against a wall -- e.g. cx in [7/2, 7/2 + 1/160] at m = 4, whose
    only admissible poses are (7/2, cy, 0) -- is tested over its whole bin, where every pose is
    inadmissible, and no primitive can certify it however deep the subdivision goes.  See
    search/RUNG2.md sec 4.6.
    """
    cx0, cx1, cy0, cy1, u0, u1 = box
    if u1 * u1 + 2 * u1 - 1 > 0:          # the bin reaches past 45deg: w is not monotone on it
        return u1
    K = 2 * min(cx1, m - cx0, cy1, m - cy0)
    if K <= 0: return u1
    def w(u):
        c, s = trig(u)
        return c + s
    if w(u1) <= K: return u1              # the whole bin can hold admissible poses
    lo, hi = u0, u1                       # w(lo) <= K < w(hi) is the invariant (w(u0) <= K, else EMPTY)
    if w(lo) > K: return u0               # nothing admissible at all; EMPTY will catch it
    if w(lo) == K: return lo              # u- = u0 exactly: the admissible bin is the single angle
    for _ in range(steps):                # u0, and only a DEGENERATE bin lets the primitives see
        mid = (lo + hi) / 2               # it -- at a wall this is the whole content of the box
        if w(mid) <= K: lo = mid
        else: hi = mid
    return lo if w(lo) == K else hi       # hi >= u-, so no admissible pose is ever dropped


# ------------------------------------------------------------------ ADM primitive (Lemma A/B/C)
# A point p is certified for a pose box (rect x u-bin) when four scalar inequalities hold for every
# u in the bin; each inequality, after clearing the positive denominator (1+u^2)^2, is a polynomial
# of degree <= 4 in u with rational coefficients.  See search/RUNG2.md Lemma A (monotone corners),
# Lemma B (the wall bounds) and Lemma C (the Bernstein enclosure).
COMB = [[1, 0, 0, 0, 0], [1, 1, 0, 0, 0], [1, 2, 1, 0, 0], [1, 3, 3, 1, 0], [1, 4, 6, 4, 1]]

def _max_quad(a0, a1, a2, u0, u1):
    """exact max of a0 + a1 u + a2 u^2 over [u0, u1] (Fractions)."""
    v = a0 + a1 * u0 + a2 * u0 * u0
    v1 = a0 + a1 * u1 + a2 * u1 * u1
    if v1 > v: v = v1
    if a2 < 0:
        uv = -a1 / (2 * a2)
        if u0 < uv < u1:
            vv = a0 + a1 * uv + a2 * uv * uv
            if vv > v: v = vv
    return v

def _max_bern(a, u0, u1):
    """Bernstein (convex-hull) upper bound for max of sum_k a[k] u^k, deg <= 4, over [u0, u1].
    Exact at the endpoints; the overestimate is O((u1-u0)^2).  Sound: the polynomial is a convex
    combination of its Bernstein coefficients at every point of the interval."""
    h = u1 - u0
    b = []
    hp = F(1)
    for j in range(5):
        sj = F(0)
        for k in range(j, 5):
            if a[k]:
                sj += a[k] * COMB[k][j] * (u0 ** (k - j))
        b.append(sj * hp)
        hp *= h
    best = b[0]
    for i in range(1, 5):
        bi = F(0)
        for j in range(i + 1):
            bi += F(COMB[i][j], COMB[4][j]) * b[j]
        if bi > best: best = bi
    return best

def _poly_ok(a, u0, u1):
    """True if max_{[u0,u1]} sum a[k] u^k <= 0, by the exact quadratic max when deg <= 2 and by the
    Bernstein bound otherwise (sound in both cases)."""
    if a[3] == 0 and a[4] == 0:
        return _max_quad(a[0], a[1], a[2], u0, u1) <= 0
    return _max_bern(a, u0, u1) <= 0

# The centre bounds, as c = Xn(u) / (2 (1+u^2)) with Xn a quadratic in u:
#   'R' a: a constant rational bound a          -> Xn = 2a (1+u^2)
#   'W':   the near-wall bound  c >= w(theta)/2 -> Xn = (1-u^2) + 2u
#   'M':   the far-wall bound   c <= m - w/2    -> Xn = 2m(1+u^2) - (1-u^2) - 2u
def _xn(kind, val, m):
    if kind == 'R': return (2 * val, F(0), 2 * val)
    if kind == 'W': return (F(1), F(2), F(-1))
    return (2 * m - 1, F(-2), 2 * m + 1)

def in_rot_square(a, b, c, s):
    """(a,b) in R_theta Q  (closed unit square centred at 0, angle theta with cos c, sin s)."""
    x = a * c + b * s
    y = -a * s + b * c
    return -HALF <= x <= HALF and -HALF <= y <= HALF

def in_core(a, b, B):
    """(a,b) in the intersection of R_t Q over t in [t0, t1]."""
    if not in_rot_square(a, b, B['c0'], B['s0']): return False
    if not in_rot_square(a, b, B['c1'], B['s1']): return False
    # direction of (a,b) relative to t0, folded mod 90deg into the first quadrant
    x = a * B['c0'] + b * B['s0']
    y = -a * B['s0'] + b * B['c0']
    if x >= 0 and y >= 0: p, q = x, y
    elif x < 0 and y >= 0: p, q = y, -x
    elif x < 0 and y < 0: p, q = -x, -y
    else: p, q = -y, x
    in_sector = (p == 0) or (q * B['cD'] <= p * B['sD'])
    if in_sector:
        return a * a + b * b <= F(1, 4)
    return True

def in_core_f(a, b, Bf, tol=1e-9):
    """float version, lenient by tol (a superset of the exact test)."""
    c0, s0, c1, s1, cD, sD = Bf
    h = 0.5 + tol
    x = a * c0 + b * s0; y = -a * s0 + b * c0
    if not (-h <= x <= h and -h <= y <= h): return False
    x1 = a * c1 + b * s1; y1 = -a * s1 + b * c1
    if not (-h <= x1 <= h and -h <= y1 <= h): return False
    if x >= 0 and y >= 0: p, q = x, y
    elif x < 0 and y >= 0: p, q = y, -x
    elif x < 0 and y < 0: p, q = -x, -y
    else: p, q = -y, x
    if p <= tol or q * cD <= p * sD + tol:
        return a * a + b * b <= 0.25 + tol
    return True

def in_core_f_vec(a, b, Bf, tol=1e-9):
    """numpy-vectorised version of in_core_f: a, b are arrays; returns a bool array, a superset
    (lenient by tol) of the exact core test, used only as a pre-filter -- every candidate it
    passes is re-checked exactly in Fractions before being counted."""
    c0, s0, c1, s1, cD, sD = Bf
    h = 0.5 + tol
    x = a * c0 + b * s0; y = -a * s0 + b * c0
    ok = (np.abs(x) <= h) & (np.abs(y) <= h)
    x1 = a * c1 + b * s1; y1 = -a * s1 + b * c1
    ok &= (np.abs(x1) <= h) & (np.abs(y1) <= h)
    # fold (x,y) into the first quadrant mod 90deg -> (p,q)
    p = np.where(x >= 0, np.where(y >= 0, x, -y), np.where(y >= 0, y, -x))
    q = np.where(x >= 0, np.where(y >= 0, y, x), np.where(y >= 0, -x, -y))
    in_sector = (p <= tol) | (q * cD <= p * sD + tol)
    r2ok = (a * a + b * b) <= 0.25 + tol
    ok &= np.where(in_sector, r2ok, True)
    return ok

class Checker:
    def __init__(self, m, points, weights=None, use_tri=False, max_depth=16, dump=None,
                 use_adm=True, theta_bias=1, use_chain=False, chain_from=0, clip=True):
        self.m = F(m)
        self.P = [(F(x), F(y)) for x, y in points]
        self.W = [F(w) for w in weights] if weights else [F(1)] * len(self.P)
        self.Pf = [(float(x), float(y)) for x, y in self.P]
        self.use_tri = use_tri
        self.use_adm = use_adm
        self.clip = clip
        self._last_inT = None
        self.use_chain = use_chain
        self.chain_from = chain_from
        self.theta_bias = theta_bias
        self.max_depth = max_depth
        self.tris = self._triangles() if use_tri else []
        self.dump = dump
        # vectorised (numpy) point/weight arrays for the float pre-filter; a fixed order by
        # decreasing weight so the exact confirmation loop reaches total >= 1 in as few
        # Fraction operations as possible.
        n = len(self.P)
        self.Pxf = np.array([p[0] for p in self.Pf], dtype=float) if n else np.zeros(0)
        self.Pyf = np.array([p[1] for p in self.Pf], dtype=float) if n else np.zeros(0)
        self.Wf = np.array([float(w) for w in self.W], dtype=float) if n else np.zeros(0)
        # exact integer numerators over a common denominator: lets CHAIN total the weight of a
        # UNION of index sets with numpy (one int64 sum) instead of adding Fractions one by one,
        # which is what makes the per-region test affordable -- and correct, since the regions'
        # witness sets overlap and their weights must NOT simply be added.
        den = 1
        for w in self.W:
            d = w.denominator
            den = den * d // math.gcd(den, d)
            if den > 10 ** 15: den = 0; break
        self.Wden = den
        self.Wnum = (np.array([int(w * den) for w in self.W], dtype=np.int64) if (den and n)
                     else np.zeros(n, dtype=np.int64))
        self.order = np.argsort(-self.Wf) if n else np.zeros(0, dtype=np.int64)

    def _triangles(self):
        n = len(self.P); T = []
        for i, j, k in itertools.combinations(range(n), 3):
            ok = True
            for a, b in ((i, j), (j, k), (i, k)):
                dx = self.P[a][0] - self.P[b][0]; dy = self.P[a][1] - self.P[b][1]
                if dx * dx + dy * dy > 1: ok = False; break
            if ok: T.append((self.P[i], self.P[j], self.P[k]))
        return T

    # ---- symmetry -------------------------------------------------------------------------
    def symmetric(self):
        S = set(self.P)
        m = self.m
        return all((m - x, y) in S for x, y in S) and all((x, m - y) in S for x, y in S)

    # ---- certification tests ----------------------------------------------------------------
    def cert_core(self, box, B, Bf):
        cx0, cx1, cy0, cy1 = box[:4]
        # clip to the widest admissible range over the bin (w minimal)
        lo = B['wlo'] / 2; hi = self.m - lo
        cx0 = max(cx0, lo); cx1 = min(cx1, hi); cy0 = max(cy0, lo); cy1 = min(cy1, hi)
        if cx0 > cx1 or cy0 > cy1: return ('EMPTY', None)
        cxf0, cxf1 = float(cx0), float(cx1); cyf0, cyf1 = float(cy0), float(cy1)
        # numpy float pre-filter over ALL points at once (superset of the exact test: every
        # point that could pass the exact test passes this one, by construction of in_core_f_vec).
        mask = in_core_f_vec(self.Pxf - cxf0, self.Pyf - cyf0, Bf)
        mask &= in_core_f_vec(self.Pxf - cxf0, self.Pyf - cyf1, Bf)
        mask &= in_core_f_vec(self.Pxf - cxf1, self.Pyf - cyf0, Bf)
        mask &= in_core_f_vec(self.Pxf - cxf1, self.Pyf - cyf1, Bf)
        if not mask.any() or self.Wf[mask].sum() < 1.0:
            return (None, None)          # even the lenient float superset can't reach weight 1
        # exact confirmation only on the (few) candidates, heaviest first
        total = F(0); used = []
        for k in self.order:
            if not mask[k]: continue
            px, py = self.P[k]
            ok = True
            for a in (px - cx0, px - cx1):
                for b in (py - cy0, py - cy1):
                    if not in_core(a, b, B): ok = False; break
                if not ok: break
            if ok:
                total += self.W[k]; used.append(int(k))
                if total >= 1: return ('CORE', used)
        return (None, None)

    # ---- ADM: the admissible-box primitive (RUNG2.md Lemmas A-C) ---------------------------
    def _adm_specs(self, box, B):
        """For each of the four centre bounds return the list of sound polynomial bounds available
        on this bin: 'R' (the box side), 'W'/'M' (the near/far container wall).  max(cx0, w/2) is
        not a polynomial in u, so when the two cross inside the bin both are offered and each
        inequality may use whichever of them certifies it (each is separately a sound bound)."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        m = self.m; whi2 = B['whi'] / 2; wlo2 = B['wlo'] / 2
        def lo(c0):
            if c0 >= whi2: return (('R', c0),)          # the wall bound is dominated: skip it
            return (('R', c0), ('W', None))
        def hi(c1):
            if c1 <= m - whi2: return (('R', c1),)
            return (('R', c1), ('M', None))
        return lo(cx0), hi(cx1), lo(cy0), hi(cy1)

    @staticmethod
    def _cond_poly(U, V, cond, both_rect):
        """The polynomial G with  G(u) <= 0 on the bin  <=>  the condition holds for every
        theta = 2 atan u in the bin.  U = 2(1+u^2) p_x - Xn, V = 2(1+u^2) p_y - Yn."""
        U0, U1, U2 = U; V0, V1, V2 = V
        if both_rect:                      # U = U0 (1+u^2), V = V0 (1+u^2): divide out (1+u^2) > 0
            if cond == 0: g = (U0 - 1, 2 * V0, -U0 - 1)
            elif cond == 1: g = (-U0 - 1, -2 * V0, U0 - 1)
            elif cond == 2: g = (V0 - 1, -2 * U0, -V0 - 1)
            else: g = (-V0 - 1, 2 * U0, V0 - 1)
            return (g[0], g[1], g[2], 0, 0)
        if cond == 0:   # X <= 1/2
            return (U0 - 1, U1 + 2 * V0, U2 - U0 + 2 * V1 - 2, -U1 + 2 * V2, -U2 - 1)
        if cond == 1:   # X >= -1/2
            return (-U0 - 1, -U1 - 2 * V0, U0 - U2 - 2 * V1 - 2, U1 - 2 * V2, U2 - 1)
        if cond == 2:   # Y <= 1/2
            return (V0 - 1, V1 - 2 * U0, V2 - V0 - 2 * U1 - 2, -V1 - 2 * U2, -V2 - 1)
        return (-V0 - 1, 2 * U0 - V1, V0 - V2 + 2 * U1 - 2, V1 + 2 * U2, V2 - 1)

    def _adm_cond_ok(self, px, py, specs, cond, u0, u1):
        """does inequality `cond` hold for p at every admissible pose of the box?  (Lemma A)"""
        m = self.m
        Axs, Bxs, Ays, Bys = specs
        # condition -> (x-slot, y-slot): X<=1/2 at (Ax,Ay); X>=-1/2 at (Bx,By);
        #              Y<=1/2 at (Bx,Ay);  Y>=-1/2 at (Ax,By)
        xs, ys = ((Axs, Ays), (Bxs, Bys), (Bxs, Ays), (Axs, Bys))[cond]
        for xk in xs:
            Xn = _xn(xk[0], xk[1], m)
            U = (2 * px - Xn[0], -Xn[1], 2 * px - Xn[2])
            for yk in ys:
                Yn = _xn(yk[0], yk[1], m)
                V = (2 * py - Yn[0], -Yn[1], 2 * py - Yn[2])
                g = self._cond_poly(U, V, cond, xk[0] == 'R' and yk[0] == 'R')
                if _poly_ok(g, u0, u1): return True
        return False

    def _adm_exact(self, px, py, specs, u0, u1):
        for cond in range(4):
            if not self._adm_cond_ok(px, py, specs, cond, u0, u1): return False
        return True

    def _adm_mask(self, specs, u0, u1, per_cond=False):
        """numpy pre-filter: a lenient superset of _adm_exact over the whole point set.  With
        per_cond, also return the four per-condition masks; a point failing one of those FAILS the
        exact test for that condition too (the float bound is a lenient over-estimate), which is
        what makes the CHAIN candidate screen cheap."""
        f0, f1 = float(u0), float(u1); h = f1 - f0
        Px, Py = self.Pxf, self.Pyf; mf = float(self.m)
        tol = 1e-9
        def xn(k, v):
            if k[0] == 'R': fv = float(k[1]); return (2 * fv, 0.0, 2 * fv)
            if k[0] == 'W': return (1.0, 2.0, -1.0)
            return (2 * mf - 1, -2.0, 2 * mf + 1)
        def bound(a, both_rect):
            if both_rect:
                v = np.maximum(a[0] + a[1] * f0 + a[2] * f0 * f0, a[0] + a[1] * f1 + a[2] * f1 * f1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    uv = np.where(a[2] < 0, -a[1] / (2 * a[2]), f0 - 1.0)
                inr = (uv > f0) & (uv < f1)
                return np.where(inr, a[0] + a[1] * uv + a[2] * uv * uv, v)
            b = []; hp = 1.0
            for j in range(5):
                sj = 0.0
                for k in range(j, 5):
                    if k == j: sj = sj + a[k]
                    else: sj = sj + a[k] * COMB[k][j] * (f0 ** (k - j))
                b.append(sj * hp); hp *= h
            best = b[0]
            for i in range(1, 5):
                bi = 0.0
                for j in range(i + 1):
                    bi = bi + (COMB[i][j] / COMB[4][j]) * b[j]
                best = np.maximum(best, bi)
            return best
        Axs, Bxs, Ays, Bys = specs
        mask = np.ones(len(Px), dtype=bool); conds = []
        for cond, (xs, ys) in enumerate(((Axs, Ays), (Bxs, Bys), (Bxs, Ays), (Axs, Bys))):
            best = None
            for xk in xs:
                Xn = xn(xk, None)
                U = (2 * Px - Xn[0], -Xn[1], 2 * Px - Xn[2])
                for yk in ys:
                    Yn = xn(yk, None)
                    V = (2 * Py - Yn[0], -Yn[1], 2 * Py - Yn[2])
                    br = (xk[0] == 'R' and yk[0] == 'R')
                    g = self._cond_poly(U, V, cond, br)
                    bv = bound(g, br)
                    best = bv if best is None else np.minimum(best, bv)
            cm = (best <= tol)
            if per_cond: conds.append(cm)
            mask &= cm
            if not per_cond and not mask.any(): break
        return (mask, conds) if per_cond else mask

    def cert_adm(self, box, B, inh=None):
        """`inh` is a set of points already proved (at an ancestor box) to lie in Q at every
        admissible pose there.  A sub-box has a SUBSET of those poses, so the proof is inherited
        verbatim -- this is the per-subtree cache, and it is exact, not a heuristic."""
        u0, u1 = box[4], box[5]
        specs = self._adm_specs(box, B)
        mask = self._adm_mask(specs, u0, u1)
        cand = mask if inh is None else (mask | inh)
        if not cand.any() or self.Wf[cand].sum() < 1.0: return (None, None)
        total = F(0); used = []
        for k in self.order:
            if not cand[k]: continue
            if (inh is not None and inh[k]) or self._adm_exact(*self.P[k], specs, u0, u1):
                total += self.W[k]; used.append(int(k))
                if total >= 1: return ('ADM', used)
        return (None, None)

    def _p1_mask(self, box, B):
        cx0, cx1, cy0, cy1 = box[:4]
        t = 1 - B['whi'] / 2
        tf = float(t); mf = float(self.m); eps = 1e-9
        cxf0, cxf1, cyf0, cyf1 = float(cx0), float(cx1), float(cy0), float(cy1)
        cond1 = (self.Pxf - 1 <= eps) | (cxf0 >= self.Pxf - tf - eps)
        cond2 = (self.Pxf + 1 >= mf - eps) | (cxf1 <= self.Pxf + tf + eps)
        cond3 = (self.Pyf - 1 <= eps) | (cyf0 >= self.Pyf - tf - eps)
        cond4 = (self.Pyf + 1 >= mf - eps) | (cyf1 <= self.Pyf + tf + eps)
        return cond1 & cond2 & cond3 & cond4

    def _p1_exact(self, px, py, box, t):
        cx0, cx1, cy0, cy1 = box[:4]; m = self.m
        return ((px - 1 <= 0 or cx0 >= px - t) and (px + 1 >= m or cx1 <= px + t)
                and (py - 1 <= 0 or cy0 >= py - t) and (py + 1 >= m or cy1 <= py + t))

    def cert_p1(self, box, B):
        mask = self._p1_mask(box, B)
        if not mask.any() or self.Wf[mask].sum() < 1.0:
            return (None, None)
        t = 1 - B['whi'] / 2
        total = F(0); used = []
        for k in self.order:
            if not mask[k]: continue
            px, py = self.P[k]
            if not self._p1_exact(px, py, box, t): continue
            total += self.W[k]; used.append(int(k))
            if total >= 1: return ('P1', used)
        return (None, None)

    def cert_mix(self, box, B, inh=None):
        """ADM and P1 certify DIFFERENT points of the same box; their union is a legitimate
        witness set (each of its points lies in Q for every admissible pose of the box), so a box
        neither primitive can carry alone is still certified if the union reaches weight 1."""
        u0, u1 = box[4], box[5]
        specs = self._adm_specs(box, B)
        ma = self._adm_mask(specs, u0, u1)
        mp = self._p1_mask(box, B)
        mask = ma | mp if inh is None else (ma | mp | inh)
        if not mask.any() or self.Wf[mask].sum() < 1.0: return (None, None)
        t = 1 - B['whi'] / 2
        total = F(0); used = []
        for k in self.order:
            if not mask[k]: continue
            px, py = self.P[k]
            if (inh is not None and inh[k]) or \
               (mp[k] and self._p1_exact(px, py, box, t)) or \
               (ma[k] and self._adm_exact(px, py, specs, u0, u1)):
                total += self.W[k]; used.append(int(k))
                if total >= 1: return ('MIX', used)
        return (None, None)

    # ---- CHAIN: the disjunctive primitive (RUNG2.md Lemmas E-G) -----------------------------
    # Violation polynomials.  With a = p_x - cx, b = p_y - cy, C = 1-u^2, S = 2u, N = 1+u^2, the
    # four containment inequalities of p are  G_k <= 0,  k = 0..3, where 2 N^2 X = ... :
    #   k=0  X <=  1/2 violated:   G = 2a C + 2b S - N = (2a-1) + 4b u + (-2a-1) u^2
    #   k=1  X >= -1/2 violated:   G = -2a C - 2b S - N = (-2a-1) - 4b u + (2a-1) u^2
    #   k=2  Y <=  1/2 violated:   G = -2a S + 2b C - N = (2b-1) - 4a u + (-2b-1) u^2
    #   k=3  Y >= -1/2 violated:   G = 2a S - 2b C - N = (-2b-1) + 4a u + (2b-1) u^2
    # Each is AFFINE in (cx, cy) and QUADRATIC in u, so any nonnegative combination of them is
    # too, and its maximum over a pose box is attained at one of the four corners of the centre
    # rectangle -- there exactly, by the quadratic vertex test.  No Bernstein slack, no depth.
    @staticmethod
    def _gcoef(kind, a, b):
        if kind == 0: return (2 * a - 1, 4 * b, -2 * a - 1)
        if kind == 1: return (-2 * a - 1, -4 * b, 2 * a - 1)
        if kind == 2: return (2 * b - 1, -4 * a, -2 * b - 1)
        return (-2 * b - 1, 4 * a, 2 * b - 1)

    def _gmax(self, terms, box):
        """exact max over the whole pose box of  sum_i lam_i * G_{p_i, kind_i}  (terms =
        [(lam, pointindex, kind), ...]).  Sound for the admissible poses too, since they are a
        subset of the box."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        best = None
        for cx in (cx0, cx1):
            for cy in (cy0, cy1):
                c0 = c1 = c2 = F(0)
                for lam, k, kind in terms:
                    px, py = self.P[k]
                    g = self._gcoef(kind, px - cx, py - cy)
                    c0 += lam * g[0]; c1 += lam * g[1]; c2 += lam * g[2]
                v = _max_quad(c0, c1, c2, u0, u1)
                if best is None or v > best: best = v
        return best

    def cert_chain(self, box, B, lams=(F(1), F(1, 2), F(2)), inh=None):
        """Disjunctive certification by monotone chains of pivots (RUNG2.md secs 6-7).

        T = the points ADM/P1 certify for the whole box.  A *swing* point is one whose four
        containment inequalities all hold on the box except exactly one, `G_p`.  Among the swing
        points of one kind pick a chain q_1, ..., q_k with G_{q_1} <= ... <= G_{q_k} everywhere on
        the box (verified exactly on consecutive pairs; transitivity does the rest).  The k+1 sets
        {G_{q_1} > 0}, {G_{q_r} <= 0 < G_{q_{r+1}}}, {G_{q_k} <= 0} partition the box, and on the
        r-th of them every q_j with j <= r is captured, as is every point a with
        G_a + lam G_{q_{r+1}} <= 0 on the box for some lam > 0.

        One chain suffices at a wall pose (one cut slides).  At an interior tile pose two cuts
        slide independently, so two chains of different kinds are used and the regions are the
        PRODUCT of the two partitions; a product region can be empty, which is certified by
        exhibiting lam > 0 with max_B (G_A + lam G_B) <= 0 -- the two pivots cannot both be
        violated at once -- and an empty region needs no witness."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        specs = self._adm_specs(box, B)
        # --- T: everything ADM or P1 certifies outright
        ma, cmasks = self._adm_mask(specs, u0, u1, per_cond=True)
        mp = self._p1_mask(box, B)
        t = 1 - B['whi'] / 2
        inT = np.zeros(len(self.P), dtype=bool)
        wT = F(0)
        for k in self.order:
            if not (ma[k] or mp[k] or (inh is not None and inh[k])): continue
            px, py = self.P[k]
            if (inh is not None and inh[k]) or \
               (mp[k] and self._p1_exact(px, py, box, t)) or \
               (ma[k] and self._adm_exact(px, py, specs, u0, u1)):
                inT[k] = True; wT += self.W[k]
        self._last_inT = inT              # the per-subtree cache handed to this box's children
        if wT >= 1: return ('ADM', [int(k) for k in np.nonzero(inT)[0]])
        # --- reachable: |p - c| <= sqrt2/2 + (box diagonal)/2 for some centre c of the box
        cxm, cym = float((cx0 + cx1) / 2), float((cy0 + cy1) / 2)
        rad = 0.7072 + 0.5 * math.hypot(float(cx1 - cx0), float(cy1 - cy0))
        reach = (((self.Pxf - cxm) ** 2 + (self.Pyf - cym) ** 2) <= rad * rad) & (~inT) & (self.Wf > 0)
        if float(wT) + self.Wf[reach].sum() < 1.0: return (None, None)
        # --- candidate swing points: exactly one of the four inequalities can fail on the box
        # a condition whose FLOAT bound already fails is exactly failed, so a point with two or
        # more float-failing conditions is never a single-swing candidate: screen those out with
        # no Fraction arithmetic at all.
        nfail = np.zeros(len(self.P), dtype=np.int8)
        for cm in cmasks: nfail += (~cm)
        reach &= (nfail <= 1)
        cand = []
        for k in np.nonzero(reach)[0]:
            px, py = self.P[k]
            bad = None
            for c in range(4):
                if not cmasks[c][k]:
                    if bad is not None: bad = -1; break
                    bad = c; continue
                if not self._adm_cond_ok(px, py, specs, c, u0, u1):
                    if bad is not None: bad = -1; break
                    bad = c
            if bad is None or bad < 0: continue
            cand.append((int(k), bad))
        if len(cand) < 2: return (None, None)
        if float(wT) + sum(self.Wf[k] for k, _ in cand) < 1.0: return (None, None)

        def gmax1(k1, d1, k2, d2, lam):
            return self._gmax([(F(1), k1, d1), (lam, k2, d2)], box)

        # --- one chain per swing kind, heaviest kind first
        kinds = sorted({kd for _, kd in cand}, key=lambda kd: -sum(self.Wf[k] for k, d in cand if d == kd))
        chains = []
        for kd in kinds:
            grp = [ck for ck in cand if ck[1] == kd]
            if len(grp) < 2: continue
            def proxy(ck):
                g = self._gcoef(ck[1], self.P[ck[0]][0] - F(cxm), self.P[ck[0]][1] - F(cym))
                return float(g[0] + g[1] * (u0 + u1) / 2)
            grp.sort(key=proxy)
            ch = []
            for ck in grp:
                if not ch: ch.append(ck); continue
                if self._gmax([(F(1), ch[-1][0], ch[-1][1]), (F(-1), ck[0], ck[1])], box) <= 0:
                    ch.append(ck)
            if len(ch) >= 1: chains.append(ch)
        if not chains: return (None, None)

        def analyse(ch):
            """prefix weights of the chain's down-sets, and the up-set weight/members per pivot"""
            kk = len(ch)
            n = len(self.P)
            down = [np.zeros(n, dtype=bool) for _ in range(kk + 1)]
            for r in range(1, kk + 1):
                down[r] = down[r - 1].copy(); down[r][ch[r - 1][0]] = True
            up = [np.zeros(n, dtype=bool) for _ in range(kk + 2)]
            for (k, kind) in cand:
                lo, hi = 0, kk
                while lo < hi:
                    mid = (lo + hi + 1) // 2
                    kq, kdq = ch[mid - 1]
                    good = kq != k and any(gmax1(k, kind, kq, kdq, lam) <= 0 for lam in lams)
                    if good: lo = mid
                    else: hi = mid - 1
                for r in range(1, lo + 1): up[r][k] = True
            return down, up

        if not self.Wden: return (None, None)         # no common denominator: CHAIN unavailable
        def enough(*masks):
            u = masks[0].copy()
            for mk in masks[1:]: u |= mk
            return int(self.Wnum[u].sum()) >= self.Wden
        info = [analyse(ch) for ch in chains]
        # --- single chain
        for ci, ch in enumerate(chains):
            down, up = info[ci]; kk = len(ch)
            if all(enough(inT, down[r], up[r + 1]) for r in range(kk + 1)):
                used = inT | down[kk]
                for r in range(1, kk + 1): used = used | up[r]
                return ('CHAIN', sorted(int(v) for v in np.nonzero(used)[0]))
        # --- product of two chains of different kinds
        for i in range(len(chains)):
            for j in range(i + 1, len(chains)):
                A, Bc = chains[i], chains[j]
                (downA, upA), (downB, upB) = info[i], info[j]
                ka, kb = len(A), len(Bc)
                # empty[r] = the largest s such that region (r, s) is provably empty (a staircase:
                # G_A and G_B are both non-decreasing along their chains, so is G_A + lam G_B)
                empt = [0] * (ka + 1)
                for r in range(ka):
                    lo, hi = 0, kb
                    while lo < hi:
                        mid = (lo + hi + 1) // 2
                        good = any(gmax1(A[r][0], A[r][1], Bc[mid - 1][0], Bc[mid - 1][1], lam) <= 0
                                   for lam in lams)
                        if good: lo = mid
                        else: hi = mid - 1
                    empt[r] = lo
                ok = True
                for r in range(ka + 1):
                    for sdx in range(kb + 1):
                        if r < ka and sdx < kb and empt[r] >= sdx + 1: continue   # region is empty
                        if not enough(inT, downA[r], downB[sdx], upA[r + 1], upB[sdx + 1]):
                            ok = False; break
                    if not ok: break
                if ok:
                    used = inT | downA[ka] | downB[kb]
                    for r in range(1, ka + 1): used = used | upA[r]
                    for r in range(1, kb + 1): used = used | upB[r]
                    return ('CHAIN', sorted(int(v) for v in np.nonzero(used)[0]))
        return (None, None)

    def cert_tri(self, box):
        cx0, cx1, cy0, cy1 = box[:4]
        corners = [(cx0, cy0), (cx1, cy0), (cx0, cy1), (cx1, cy1)]
        for (A, Bp, C) in self.tris:
            if all(self._in_tri(c, A, Bp, C) for c in corners):
                wmin = min(self.W[self.P.index(v)] for v in (A, Bp, C))
                if wmin >= 1: return ('TRI', (A, Bp, C))
        return (None, None)

    @staticmethod
    def _in_tri(c, A, B, C):
        def cross(o, p, q): return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])
        d1 = cross(A, B, c); d2 = cross(B, C, c); d3 = cross(C, A, c)
        return (d1 >= 0 and d2 >= 0 and d3 >= 0) or (d1 <= 0 and d2 <= 0 and d3 <= 0)

    # ---- recursion --------------------------------------------------------------------------
    def run_box(self, root):
        """Certify one root box; returns stats dict and the list of uncertified boxes."""
        stats = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'EMPTY': 0,
                 'UNCERT': 0, 'boxes': 0, 'maxdepth': 0}
        unc = []; leaves = []
        stack = [(root, 0, None)]
        bincache = {}
        while stack:
            box, depth, inh = stack.pop()
            stats['boxes'] += 1
            stats['maxdepth'] = max(stats['maxdepth'], depth)
            cx0, cx1, cy0, cy1, u0, u1 = box
            if self.clip and u1 > u0:
                cu1 = clip_bin(box, self.m)
                if cu1 < u1:
                    u1 = cu1; box = (cx0, cx1, cy0, cy1, u0, u1)
            key = (u0, u1)
            if key not in bincache:
                B = bin_data(u0, u1)
                Bf = tuple(float(B[k]) for k in ('c0', 's0', 'c1', 's1', 'cD', 'sD'))
                bincache[key] = (B, Bf)
            B, Bf = bincache[key]
            # no admissible pose at all?
            lo = B['wlo'] / 2
            if cx1 < lo or cx0 > self.m - lo or cy1 < lo or cy0 > self.m - lo:
                stats['EMPTY'] += 1; leaves.append((box, 'EMPTY', None)); continue
            self._last_inT = None
            if self.use_adm:
                kind, wit = self.cert_adm(box, B, inh)
            else:
                kind, wit = self.cert_core(box, B, Bf)
            if kind is None:
                kind, wit = self.cert_p1(box, B)
            if kind is None and self.use_adm:
                kind, wit = self.cert_mix(box, B, inh)
            if kind is None and self.use_chain and depth >= self.chain_from:
                kind, wit = self.cert_chain(box, B, inh=inh)
            if kind is None and self.use_tri:
                kind, wit = self.cert_tri(box)
            if kind is not None:
                stats[kind] += 1
                if self.dump: leaves.append((box, kind, wit))
                continue
            if depth >= self.max_depth:
                stats['UNCERT'] += 1; unc.append(box); continue
            # split the longest dimension (theta measured in radians ~ 2 du); near a wall with
            # theta small the admissible width margin w(theta)/2 - 1/2 grows linearly in theta, so
            # precision in the centre needs precision in the angle: bias the split towards u there
            # (ZEROMARGIN.md sec 7 / FAMILY.md sec 2b diagnosis).
            dx, dy, du = cx1 - cx0, cy1 - cy0, 2 * (u1 - u0)
            if self.theta_bias > 1 and u0 * 2 < B['wlo'] and (
                    cx0 < B['whi'] / 2 or cx1 > self.m - B['whi'] / 2 or
                    cy0 < B['whi'] / 2 or cy1 > self.m - B['whi'] / 2):
                du = du * self.theta_bias
            kid = self._last_inT if self._last_inT is not None else inh
            if dx >= dy and dx >= du:
                mid = (cx0 + cx1) / 2
                stack.append(((cx0, mid, cy0, cy1, u0, u1), depth + 1, kid))
                stack.append(((mid, cx1, cy0, cy1, u0, u1), depth + 1, kid))
            elif dy >= du:
                mid = (cy0 + cy1) / 2
                stack.append(((cx0, cx1, cy0, mid, u0, u1), depth + 1, kid))
                stack.append(((cx0, cx1, mid, cy1, u0, u1), depth + 1, kid))
            else:
                mid = (u0 + u1) / 2
                stack.append(((cx0, cx1, cy0, cy1, u0, mid), depth + 1, kid))
                stack.append(((cx0, cx1, cy0, cy1, mid, u1), depth + 1, kid))
        return stats, unc, leaves

def _worker(args):
    chk, root = args
    return chk.run_box(root)

def roots(m, pitch=F(1, 10), ubins=8, cy_max=None, cx_lo=None, cx_hi=None):
    """Root boxes aligned to the pitch grid (so tight poses sit on box boundaries).
    cx_lo/cx_hi restrict the sweep to a band of centre-x columns -- a *partial* run, useful for
    timing or for isolating a region; it proves nothing about the poses it skips, and the summary
    says so."""
    m = F(m); cy_max = m / 2 if cy_max is None else F(cy_max)
    nx = int(m / pitch); ny = int(cy_max / pitch)
    R = []
    for i in range(nx):
        if cx_lo is not None and (i + 1) * pitch <= F(cx_lo): continue
        if cx_hi is not None and i * pitch >= F(cx_hi): continue
        for j in range(ny):
            for k in range(ubins):
                R.append((i * pitch, (i + 1) * pitch, j * pitch, (j + 1) * pitch,
                          F(k, 2 * ubins), F(k + 1, 2 * ubins)))
    return R

FRIEDMAN14 = [(1, 1), (F(8, 5), 1), (F(12, 5), 1), (3, 1),
              (1, F(9, 5)), (2, F(9, 5)), (3, F(9, 5)),
              (1, F(11, 5)), (2, F(11, 5)), (3, F(11, 5)),
              (1, 3), (F(8, 5), 3), (F(12, 5), 3), (3, 3)]

def read_cert(path):
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    pts, ws = [], []
    for i in range(n):
        X, Y, w = map(int, tok[5 + 3 * i: 8 + 3 * i])
        pts.append((F(X, D), F(Y, D))); ws.append(F(w, W))
    return F(sn, sd), pts, ws

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what', help='friedman14 | cert | pose')
    ap.add_argument('path', nargs='?')
    ap.add_argument('--tri', action='store_true')
    ap.add_argument('--no-adm', action='store_true', help='use the old CORE primitive instead of ADM')
    ap.add_argument('--no-clip', action='store_true',
                    help='do not shrink a box angle bin to its admissible sub-bin (clip_bin)')
    ap.add_argument('--disj', action='store_true',
                    help='enable the CHAIN disjunctive primitive (RUNG2.md sec 6): a monotone chain '
                         'of pivot inequalities partitions the box and each region gets its own '
                         'witness set.  Required for any cover of total weight below m^2.')
    ap.add_argument('--chain-from', type=int, default=0,
                    help='only attempt CHAIN at subdivision depth >= this (it is the most expensive '
                         'primitive; the cheap ones handle almost every box)')
    ap.add_argument('--theta-bias', type=int, default=4,
                    help='near a container wall with theta small, weight the u-dimension by this '
                         'factor when choosing which dimension to halve (1 = the old rule)')
    ap.add_argument('--depth', type=int, default=16)
    ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--pitch', type=str, default='1/10')
    ap.add_argument('--ubins', type=int, default=8)
    ap.add_argument('--dump', type=str, default=None)
    ap.add_argument('--oracle', type=str, default=None,
                     help='write the uncertified boxes worst-case poses (centre + 4 corners x 2 angle '
                          'endpoints, floats derived from the exact box) as (cx, cy, theta_rad) rows, '
                          'one per line, for use as separating cutting planes in a cover LP')
    ap.add_argument('--full', action='store_true', help='no symmetry reduction (cy up to m, u up to 1)')
    ap.add_argument('--cx-lo', type=str, default=None,
                    help='restrict the sweep to root columns with cx >= this (a PARTIAL run: it '
                         'proves nothing about the columns it skips)')
    ap.add_argument('--cx-hi', type=str, default=None, help='... and cx <= this')
    ap.add_argument('--u', type=str, default=None, help='pose mode: u = tan(theta/2), as a Fraction-parseable string ("1/3", "0.3333")')
    ap.add_argument('--cx', type=str, default=None, help='pose mode: centre x, Fraction-parseable')
    ap.add_argument('--cy', type=str, default=None, help='pose mode: centre y, Fraction-parseable')
    a = ap.parse_args()
    if a.what == 'pose':
        # exact captured weight at one rational pose (cx, cy, theta = 2*atan(u)): a rigorous,
        # non-adaptive confirmation that a specific pose found by a float search (stress test,
        # local polish, ...) really is (or is not) a violation of the cover -- no subdivision,
        # no depth limit, pure Fraction arithmetic; the pose's u must be rational for cos/sin to
        # be rational (trig(u)), so snap any float u to a nearby fraction first.
        m, pts, ws = read_cert(a.path)
        u = F(a.u); cx = F(a.cx); cy = F(a.cy)
        c, s = trig(u)
        total = F(0); used = []
        for k, (px, py) in enumerate(pts):
            dx, dy = px - cx, py - cy
            x = dx * c + dy * s; y = -dx * s + dy * c
            if -HALF <= x <= HALF and -HALF <= y <= HALF:
                total += ws[k]; used.append((k, px, py, ws[k]))
        th_deg = math.degrees(2 * math.atan(float(u)))
        w_adm = abs(c) + abs(s)
        adm = (w_adm / 2 <= cx <= m - w_adm / 2) and (w_adm / 2 <= cy <= m - w_adm / 2)
        print(f"container [0,{m}]^2; pose cx={cx} ({float(cx):.9f}) cy={cy} ({float(cy):.9f}) "
              f"u={u} theta={th_deg:.9f} deg; admissible: {adm}")
        print(f"EXACT captured weight = {total} = {float(total):.12f}  "
              f"({'OK: >= 1' if total >= 1 else '*** VIOLATION: < 1 ***'})")
        print(f"captured points ({len(used)}): {[(k, str(px), str(py), str(w)) for k, px, py, w in used]}")
        sys.exit(0)
    if a.what == 'friedman14':
        m, pts, ws = 4, FRIEDMAN14, None
    else:
        m, pts, ws = read_cert(a.path)
    chk = Checker(m, pts, ws, use_tri=a.tri, max_depth=a.depth, dump=a.dump,
                  use_adm=not a.no_adm, theta_bias=a.theta_bias,
                  use_chain=a.disj, chain_from=a.chain_from, clip=not a.no_clip)
    sym = chk.symmetric()
    print(f"container [0,{m}]^2, {len(pts)} points, total weight {float(sum(chk.W)):.6f}, "
          f"symmetric under x->m-x and y->m-y: {sym}, triangles: {len(chk.tris) if a.tri else 'off'}")
    if not sym and not a.full:
        print("point set not symmetric: use --full"); sys.exit(2)
    pitch = F(a.pitch)
    if a.full:
        R = roots(m, pitch, a.ubins * 2, cy_max=m, cx_lo=a.cx_lo, cx_hi=a.cx_hi)
        R = [(x0, x1, y0, y1, u0 * 2, u1 * 2) for (x0, x1, y0, y1, u0, u1) in R]  # u in [0,1]: theta to 90deg
    else:
        R = roots(m, pitch, a.ubins, cx_lo=a.cx_lo, cx_hi=a.cx_hi)
    if a.cx_lo is not None or a.cx_hi is not None:
        print(f"PARTIAL SWEEP: cx restricted to [{a.cx_lo}, {a.cx_hi}] -- this is not a verification "
              f"of the whole container")
    print(f"{len(R)} root boxes, depth limit {a.depth}, pitch {pitch}, u-bins {a.ubins}")
    t0 = time.time()
    tot = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'EMPTY': 0, 'UNCERT': 0,
           'boxes': 0, 'maxdepth': 0}
    unc_all = []; leaves_all = []
    with Pool(a.nproc) as pool:
        for i, (st, unc, leaves) in enumerate(pool.imap_unordered(_worker, [(chk, r) for r in R], chunksize=4)):
            for k in tot:
                tot[k] = max(tot[k], st[k]) if k == 'maxdepth' else tot[k] + st[k]
            unc_all += unc; leaves_all += leaves
            if (i + 1) % 500 == 0:
                print(f"  {i+1}/{len(R)} roots, {tot['boxes']} boxes, uncert {tot['UNCERT']}, {time.time()-t0:.0f}s", flush=True)
    print(f"done in {time.time()-t0:.0f}s: boxes {tot['boxes']}, max depth {tot['maxdepth']}")
    print(f"  leaves: ADM {tot['ADM']}  CORE {tot['CORE']}  P1 {tot['P1']}  MIX {tot['MIX']}  "
          f"CHAIN {tot['CHAIN']}  TRI {tot['TRI']}  EMPTY {tot['EMPTY']}  "
          f"UNCERTIFIED {tot['UNCERT']}")
    if unc_all:
        print("uncertified boxes (cx0 cx1 cy0 cy1 u0 u1 -> theta0 theta1 deg):")
        for b in sorted(unc_all)[:40]:
            print("  ", *[f"{float(v):.6f}" for v in b[:4]],
                  f"{math.degrees(2*math.atan(float(b[4]))):.3f} {math.degrees(2*math.atan(float(b[5]))):.3f}")
        if len(unc_all) > 40: print(f"  ... {len(unc_all)} in total")
    if a.dump:
        with open(a.dump, 'w') as f:
            f.write(f"# container {m}; box: cx0 cx1 cy0 cy1 u0 u1 ; kind ; witness\n")
            for box, kind, wit in leaves_all:
                f.write(" ".join(str(v) for v in box) + f" ; {kind} ; {wit}\n")
        print(f"leaves written to {a.dump}")
    if a.oracle:
        with open(a.oracle, 'w') as f:
            f.write(f"# container {m}; uncertified box worst-case poses: cx cy theta_rad\n")
            for box in unc_all:
                bx0, bx1, by0, by1, bu0, bu1 = box
                th0, th1 = 2 * math.atan(float(bu0)), 2 * math.atan(float(bu1))
                cxs = sorted({float(bx0), float(bx1), float((bx0 + bx1) / 2)})
                cys = sorted({float(by0), float(by1), float((by0 + by1) / 2)})
                for cxv in cxs:
                    for cyv in cys:
                        for thv in (th0, th1, (th0 + th1) / 2):
                            f.write(f"{cxv!r} {cyv!r} {thv!r}\n")
        print(f"oracle rows (uncertified-box corner/centre poses) written to {a.oracle}: "
              f"{len(unc_all)} boxes -> up to {len(unc_all) * 27} poses")
    partial = a.cx_lo is not None or a.cx_hi is not None
    print(("VERIFIED" if tot['UNCERT'] == 0 else "NOT VERIFIED")
          + (f" (PARTIAL: cx in [{a.cx_lo}, {a.cx_hi}] only)" if partial else ""))

if __name__ == '__main__':
    main()
