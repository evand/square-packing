"""Vectorised interval forward-mode AD (numpy), the batch version of iad.py (2026-10-21).

VI: arrays lo, hi (shape (N,)).  Directed rounding by error-free transforms: TwoSum for sums, Dekker's TwoProduct
(Veltkamp splitting) for products, so exact results stay exact and inexact ones move outward by one ulp.
cos / sin: endpoint values padded by 2 ulps (glibc < 1 ulp), exact at input 0; interior extrema by integer tests.
VD: value v, one derivative g (the ray parameter rho), base-point value v0; enc() = v meet (v0 + RHO g).
Semantics identical to iad.py (r3local.py's scalar checker), which stays the reference.
"""
import math
import numpy as np

_INF = np.inf
_SPLIT = 134217729.0      # 2^27 + 1


def _dn(x):
    return np.nextafter(x, -_INF)


def _up(x):
    return np.nextafter(x, _INF)


def _two_sum_err(a, b, s):
    bb = s - a
    return (a - (s - bb)) + (b - bb)


def _split(a):
    c = _SPLIT * a
    ah = c - (c - a)
    return ah, a - ah


def _two_prod_err(a, b, p):
    ah, al = _split(a); bh, bl = _split(b)
    return ((ah * bh - p) + ah * bl + al * bh) + al * bl


class VI:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = np.asarray(lo, dtype=float)
        self.lo = lo
        self.hi = lo if hi is None else np.asarray(hi, dtype=float)

    @staticmethod
    def c(x, n):
        a = np.full(n, float(x)); return VI(a, a.copy())

    def __add__(s, o):
        if not isinstance(o, VI): o = VI(np.broadcast_to(np.asarray(o, float), s.lo.shape))
        lo = s.lo + o.lo; hi = s.hi + o.hi
        lo = np.where(_two_sum_err(s.lo, o.lo, lo) < 0, _dn(lo), lo)
        hi = np.where(_two_sum_err(s.hi, o.hi, hi) > 0, _up(hi), hi)
        return VI(lo, hi)
    __radd__ = __add__

    def __neg__(s):
        return VI(-s.hi, -s.lo)

    def __sub__(s, o):
        if not isinstance(o, VI): o = VI(np.broadcast_to(np.asarray(o, float), s.lo.shape))
        return s + (-o)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        if not isinstance(o, VI): o = VI(np.broadcast_to(np.asarray(o, float), s.lo.shape))
        A = np.stack([s.lo, s.lo, s.hi, s.hi]); Bm = np.stack([o.lo, o.hi, o.lo, o.hi])
        P = A * Bm
        E = _two_prod_err(A, Bm, P)              # true product = P + E
        lo = P.min(axis=0); hi = P.max(axis=0)
        lo_bad = ((P == lo) & (E < 0)).any(axis=0)
        hi_bad = ((P == hi) & (E > 0)).any(axis=0)
        return VI(np.where(lo_bad, _dn(lo), lo), np.where(hi_bad, _up(hi), hi))
    __rmul__ = __mul__

    def meet(s, o):
        lo = np.maximum(s.lo, o.lo); hi = np.minimum(s.hi, o.hi)
        ok = lo <= hi
        return VI(np.where(ok, lo, s.lo), np.where(ok, hi, s.hi))


def _cosb(t):
    c = np.cos(t)
    exact = t == 0.0
    return np.where(exact, 1.0, _dn(_dn(c))), np.where(exact, 1.0, _up(_up(c)))


def _sinb(t):
    c = np.sin(t)
    exact = t == 0.0
    return np.where(exact, 0.0, _dn(_dn(c))), np.where(exact, 0.0, _up(_up(c)))


def vcos(x):
    a, b = x.lo, x.hi
    la, ha = _cosb(a); lb, hb = _cosb(b)
    lo = np.minimum(la, lb); hi = np.maximum(ha, hb)
    k0 = np.ceil(_dn(a) / math.pi); k1 = np.floor(_up(b) / math.pi)
    has = k1 >= k0
    # an even multiple of pi inside -> max 1; odd -> min -1 (intervals here are far shorter than 2 pi)
    even_in = has & ((np.mod(k0, 2) == 0) | (k1 > k0))
    odd_in = has & ((np.mod(k0, 2) == 1) | (k1 > k0))
    hi = np.where(even_in, 1.0, hi); lo = np.where(odd_in, -1.0, lo)
    return VI(np.maximum(lo, -1.0), np.minimum(hi, 1.0))


def vsin(x):
    a, b = x.lo, x.hi
    la, ha = _sinb(a); lb, hb = _sinb(b)
    lo = np.minimum(la, lb); hi = np.maximum(ha, hb)
    # extrema at pi/2 + k pi (k even: +1, odd: -1); conservative inclusion with outward-rounded ends
    k0 = np.ceil((_dn(a) - math.pi / 2) / math.pi - 1e-12); k1 = np.floor((_up(b) - math.pi / 2) / math.pi + 1e-12)
    has = k1 >= k0
    even_in = has & ((np.mod(k0, 2) == 0) | (k1 > k0))
    odd_in = has & ((np.mod(k0, 2) == 1) | (k1 > k0))
    hi = np.where(even_in, 1.0, hi); lo = np.where(odd_in, -1.0, lo)
    return VI(np.maximum(lo, -1.0), np.minimum(hi, 1.0))


def vdiv(a, b):
    """a / b, b excluding 0 (checked by the caller: raises if any b contains 0)."""
    if np.any((b.lo <= 0) & (b.hi >= 0)):
        raise ZeroDivisionError('vdiv by an interval containing 0')
    los, his = [], []
    for x in (a.lo, a.hi):
        for y in (b.lo, b.hi):
            q = x / y
            p = q * y
            e = _two_prod_err(q, y, p)          # q y = p + e exactly
            r = (x - p) - e                      # x - q y (x - p exact by Sterbenz: p is within a factor 2 of x)
            d = r / y                            # sign of the true x/y - q
            los.append(np.where(d < 0, _dn(q), q)); his.append(np.where(d > 0, _up(q), q))
    return VI(np.min(los, axis=0), np.max(his, axis=0))


class VD:
    __slots__ = ('v', 'g', 'v0', 'R')

    def __init__(self, v, g, v0, R):
        self.v, self.g, self.v0, self.R = v, g, v0, R      # R: VI of rho ranges (None: no base enclosure)

    def enc(self):
        if self.R is None: return self.v
        return self.v.meet(self.v0 + self.R * self.g)

    @staticmethod
    def c(x, n, R=None):
        v = x if isinstance(x, VI) else VI.c(x, n)
        return VD(v, VI.c(0.0, n), v, R)

    def __add__(s, o):
        if not isinstance(o, VD): return VD(s.v + o, s.g, s.v0 + o, s.R)
        return VD(s.v + o.v, s.g + o.g, s.v0 + o.v0, s.R if s.R is not None else o.R)
    __radd__ = __add__

    def __neg__(s):
        return VD(-s.v, -s.g, -s.v0, s.R)

    def __sub__(s, o):
        return s + (-o)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        if not isinstance(o, VD):
            return VD(s.v * o, s.g * o, s.v0 * o, s.R)
        sv, ov = s.enc(), o.enc()
        return VD(sv * ov, s.g * ov + sv * o.g, s.v0 * o.v0, s.R if s.R is not None else o.R)
    __rmul__ = __mul__


def dcos(x):
    v = x.enc()
    return VD(vcos(v), -vsin(v) * x.g, vcos(x.v0), x.R)


def dsin(x):
    v = x.enc()
    return VD(vsin(v), vcos(v) * x.g, vsin(x.v0), x.R)


def dabs(x):
    e = x.enc()
    v0 = x.v0
    a0lo = np.where(v0.lo >= 0, v0.lo, np.where(v0.hi <= 0, -v0.hi, 0.0))
    a0hi = np.where(v0.lo >= 0, v0.hi, np.where(v0.hi <= 0, -v0.lo, np.maximum(-v0.lo, v0.hi)))
    pos = e.lo >= 0; neg = e.hi <= 0; mix = ~(pos | neg)
    vlo = np.where(pos, e.lo, np.where(neg, -e.hi, 0.0))
    vhi = np.where(pos, e.hi, np.where(neg, -e.lo, np.maximum(-e.lo, e.hi)))
    # derivative: +g, -g, or the Clarke hull g * [-1, 1]
    gm = np.maximum(np.abs(x.g.lo), np.abs(x.g.hi))
    glo = np.where(pos, x.g.lo, np.where(neg, -x.g.hi, -gm))
    ghi = np.where(pos, x.g.hi, np.where(neg, -x.g.lo, gm))
    # base value of |x|: exact branch when the base value's sign is known
    b0lo = np.where(pos & (v0.lo >= 0), v0.lo, np.where(neg & (v0.hi <= 0), -v0.hi, a0lo))
    b0hi = np.where(pos & (v0.lo >= 0), v0.hi, np.where(neg & (v0.hi <= 0), -v0.lo, a0hi))
    return VD(VI(vlo, vhi), VI(glo, ghi), VI(b0lo, b0hi), x.R)


def dfsq(x):
    c, s = dcos(x), dsin(x)
    return (dabs(c) + dabs(s)) * 0.5


def ddiv(x, y):
    xv, yv = x.enc(), y.enc()
    q = vdiv(xv, yv)
    y2 = yv * yv
    return VD(q, vdiv(x.g * yv - xv * y.g, y2), vdiv(x.v0, y.v0), x.R if x.R is not None else y.R)
