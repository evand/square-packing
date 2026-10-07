"""Leaf bounds for the wall-strip demonstration: lower bounds on 'node values' over a box of poses.

Exact mode (`exact=True`, the default) uses fractions.Fraction throughout: every returned number is
a rational that is a PROVEN lower bound.  Float mode (`exact=False`) runs the same formulas in
floating point; it is used only to guide the search (search.py) and proves nothing.

Poses.  A pose is (theta, c): a closed unit square at angle theta and centre height c.  Angles are
taken mod 90 deg, so theta in [0, 90 deg].  Parametrise by t = tan(theta/2) in [0, 1]:

    C = cos theta = (1 - t^2)/(1 + t^2),   S = sin theta = 2t/(1 + t^2)     (rational in t)

C decreases and S increases with t.  A *box* is [t0, t1] x [c0, c1] with rational endpoints.

Half-chords.  For the line y = z put e = z - c.  For 0 < theta < 90 deg the chord of the square on
that line is [x - lL, x + lR] (x = centre abscissa) with   (notes/chord-lemma.md, proof of Lemma 1)

    lR(e) = min( (1/2 - e S)/C , (e C + 1/2)/S )        lL(e) = min( (1/2 + e S)/C , (1/2 - e C)/S )

and the bounding-box half-width is p = (C + S)/2.  The NODE VALUE of a square whose left separation
line is y = z1 and right separation line is y = z2 is

    V(z1, z2) = lL(z1) + lR(z2),   V(END, z2) = p + lR(z2),   V(z1, END) = lL(z1) + p.

Two independent lower bounds are computed and the larger is returned:

 (J) joint:  expanding the two minima, V(z1, z2) = min of four closed forms (d = z1 - z2)
        A  = (1 + d S)/C                      B  = (1 - d C)/S
        Cc = (p - c + z2 + d S^2)/(C S)       Dd = (p + c - z1 + d S^2)/(C S)
     (Cc, Dd use C^2 + S^2 = 1).  Each form is linear in c, so its minimum over c is at c = cl or
     c = c1; over theta it is bounded by interval arithmetic (numerator low, denominator high).
 (S) separate: lower-bound lL(z1) and lR(z2) separately.  Each of the four terms is
        T1(k) = (1/2 - k S)/C   or   T2(k) = (k C + 1/2)/S     (k = +-e)
     whose exact minimum over a theta-interval is at an endpoint or at the interior critical point
     S = 2k (resp. C = -2k), where the value is sqrt(1 - 4k^2)/2 (rounded DOWN to a rational).
 (G) good-line floor: if every pose in the box has |c - z| <= D(theta) (the chord band of
     notes/chord-lemma.md Cor. 2) then the chord on y = z is >= 1, so V(z, z) >= 1 and
     V(END, z), V(z, END) >= 1 (because lL, lR <= p).

The two endpoint angles theta = 0 and theta = 90 deg (t = 0, t = 1) make S or C vanish; there the
formulas are read as limits, and check.py separately verifies that each bound is <= the true value
at those two poses (both are axis-parallel squares with V = 1 for every line meeting the interior).
"""
from fractions import Fraction as Fr
import math

NEG = None   # "no bound" (minus infinity)


def CS(t):
    d = 1 + t * t
    return (1 - t * t) / d, 2 * t / d


def Dof(u):
    """Chord-band half-height D as a function of u = C + S (notes/chord-lemma.md, Cor. 2)."""
    return (u - u * u + 1) / 2


SQRT2_UB = Fr(14142136, 10 ** 7)   # > sqrt 2 (asserted in self_test)


def sqrt_lb_exact(x):
    """A rational r with 0 <= r <= sqrt(x), x a nonnegative Fraction, error < 1e-15."""
    if x <= 0:
        return Fr(0)
    Q = 10 ** 15
    return Fr(math.isqrt(x.numerator * Q * Q // x.denominator), Q)


def contains45(t0, t1):
    """t0 <= tan 22.5deg = sqrt2 - 1 <= t1, decided exactly: sqrt2 - 1 <= t  iff  (t+1)^2 >= 2."""
    return (t0 + 1) ** 2 <= 2 <= (t1 + 1) ** 2


def _min(vals):
    if any(v is NEG for v in vals):
        return NEG
    return min(vals)


def _div(num, dmin, dmax):
    """Lower bound of n/d over n >= num, d in [dmin, dmax] (0 <= dmin <= dmax, dmax > 0)."""
    if num >= 0:
        return num / dmax
    return num / dmin if dmin > 0 else NEG


class Box:
    def __init__(self, t0, t1, c0, c1, exact=True):
        conv = Fr if exact else float
        self.exact = exact
        self.t0, self.t1, self.c0, self.c1 = conv(t0), conv(t1), conv(c0), conv(c1)
        assert 0 <= self.t0 < self.t1 <= 1 and self.c0 < self.c1
        self.C1, self.S0 = CS(self.t0)        # cos is largest, sin smallest at t0
        self.C0, self.S1 = CS(self.t1)
        u_end = [self.C1 + self.S0, self.C0 + self.S1]
        self.ins = contains45(self.t0, self.t1)
        # u = C + S = sqrt2 cos(theta - 45deg): min at an endpoint, max sqrt2 at 45deg
        self.plb = min(u_end) / 2
        u_ub = (SQRT2_UB if exact else math.sqrt(2)) if self.ins else max(u_end)
        self.Dlb = Dof(u_ub)                  # D is decreasing in u on [1, sqrt2]
        cs_end = [self.C1 * self.S0, self.C0 * self.S1]   # C S = sin(2 theta)/2: unimodal
        self.CSmin = min(cs_end)
        self.CSmax = (Fr(1, 2) if exact else 0.5) if self.ins else max(cs_end)
        self.cl = max(self.c0, self.plb)      # c >= p >= plb for every pose inside the container
        self.sqrt_lb = sqrt_lb_exact if exact else (lambda x: math.sqrt(max(x, 0)))
        self.half = Fr(1, 2) if exact else 0.5

    def key(self):
        return (self.t0, self.t1, self.c0, self.c1)

    def nonempty(self):
        return self.cl < self.c1

    def may_be_H(self, a):
        """True if the box may contain a pose with c > a + D(theta)."""
        return self.c1 > a + self.Dlb

    def good(self, z):
        """Every pose in the box has |c - z| <= D(theta)  (so its chord on y = z has length >= 1)."""
        return self.cl >= z - self.Dlb and self.c1 <= z + self.Dlb

    # ---- separate half-chord bounds (S) ------------------------------------------------------
    def _T1min(self, k):
        """min over the theta-interval of (1/2 - k S)/C."""
        C0, C1, S0, S1, h = self.C0, self.C1, self.S0, self.S1, self.half
        vals = [(h - k * S0) / C1]
        vals.append((h - k * S1) / C0 if C0 > 0 else (None if h - k * S1 <= 0 else float('inf')))
        if vals[-1] is None:
            return NEG
        if S0 < 2 * k < S1:
            vals.append(self.sqrt_lb(1 - 4 * k * k) / 2)
        return min(v for v in vals if v != float('inf'))

    def _T2min(self, k):
        """min over the theta-interval of (k C + 1/2)/S."""
        C0, C1, S0, S1, h = self.C0, self.C1, self.S0, self.S1, self.half
        vals = [(k * C0 + h) / S1]
        vals.append((k * C1 + h) / S0 if S0 > 0 else (None if k * C1 + h <= 0 else float('inf')))
        if vals[-1] is None:
            return NEG
        if C0 < -2 * k < C1:
            vals.append(self.sqrt_lb(1 - 4 * k * k) / 2)
        return min(v for v in vals if v != float('inf'))

    def half_lb(self, z, side):
        """Lower bound of lR(z) (side 'R') or lL(z) (side 'L') over the box."""
        vals = []
        for c in (self.cl, self.c1):
            e = z - c
            k = e if side == 'R' else -e      # lR: (1/2 - eS)/C, (eC + 1/2)/S ; lL: same with k = -e
            vals += [self._T1min(k), self._T2min(k)]
        return _min(vals)

    # ---- joint bound (J) ---------------------------------------------------------------------
    def joint_lb(self, z1, z2):
        C0, C1, S0, S1, h = self.C0, self.C1, self.S0, self.S1, self.half
        vals = []
        for c in (self.cl, self.c1):
            if z1 is None or z2 is None:
                z = z2 if z1 is None else z1
                e = z - c
                k = e if z1 is None else -e
                v1 = _div(min(h - k * S0, h - k * S1), C0, C1)
                v2 = _div(min(k * C0 + h, k * C1 + h), S0, S1)
                vals += [NEG if v is NEG else self.plb + v for v in (v1, v2)]
                continue
            d = z1 - z2
            dS2 = min(d * S0 * S0, d * S1 * S1)
            vals.append(_div(min(1 + d * S0, 1 + d * S1), C0, C1))
            vals.append(_div(min(1 - d * C0, 1 - d * C1), S0, S1))
            vals.append(_div(self.plb - c + z2 + dS2, self.CSmin, self.CSmax))
            vals.append(_div(self.plb + c - z1 + dS2, self.CSmin, self.CSmax))
        return _min(vals)

    def V(self, z1, z2):
        """Proven lower bound of the node value V(z1, z2) over the box (None = END)."""
        assert not (z1 is None and z2 is None)
        cands = [self.joint_lb(z1, z2)]
        l = self.plb if z1 is None else self.half_lb(z1, 'L')
        r = self.plb if z2 is None else self.half_lb(z2, 'R')
        cands.append(NEG if (l is NEG or r is NEG) else l + r)
        one = Fr(1) if self.exact else 1.0
        if z1 is not None and z2 is not None and z1 == z2 and self.good(z1):
            cands.append(one)
        if (z1 is None) != (z2 is None) and self.good(z2 if z1 is None else z1):
            cands.append(one)
        cands = [v for v in cands if v is not NEG]
        return max(cands) if cands else NEG


def self_test():
    assert SQRT2_UB * SQRT2_UB > 2
    assert contains45(Fr(0), Fr(1)) and not contains45(Fr(0), Fr(2, 5)) and not contains45(Fr(1, 2), Fr(1))
    r = sqrt_lb_exact(Fr(2))
    assert r * r <= 2 and (r + Fr(1, 10 ** 14)) ** 2 > 2
    return True
