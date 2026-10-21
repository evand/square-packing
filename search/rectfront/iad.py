"""Interval forward-mode AD (value + gradient, all intervals), float64 with outward padding (2026-10-20).

I(lo, hi): interval.  D(v, g): v an I, g a tuple of I (partial derivatives).  Functions cos, sin, abs_s (|x| with
the sign of x known on the domain, else Clarke hull), fsq (|cos x| + |sin x|)/2 (generalised gradient hull across
its kinks at multiples of pi/2).  Angles in radians.  Used by r3local.py; margins are required to exceed 1e-9, far
above the padding (1e-15 per op relative).
"""
import math

PAD = 1e-15
_INF = math.inf


def _dn(x):
    return math.nextafter(x, -_INF)


def _up(x):
    return math.nextafter(x, _INF)


def _add_lo(a, b):
    s = a + b
    bb = s - a; err = (a - (s - bb)) + (b - bb)          # TwoSum: a + b = s + err exactly
    return s if err >= 0 else _dn(s)


def _add_hi(a, b):
    s = a + b
    bb = s - a; err = (a - (s - bb)) + (b - bb)
    return s if err <= 0 else _up(s)


def _mul_lo(a, b):
    p = a * b
    err = math.fma(a, b, -p)                              # a b = p + err exactly
    return p if err >= 0 else _dn(p)


def _mul_hi(a, b):
    p = a * b
    err = math.fma(a, b, -p)
    return p if err <= 0 else _up(p)


def _o(lo, hi):
    """Outward padding for results of library functions (cos, sin: < 1 ulp in glibc; we pad 2 ulps)."""
    return _dn(_dn(lo)), _up(_up(hi))


class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        if hi is None: hi = lo
        self.lo, self.hi = float(lo), float(hi)

    def __add__(s, o):
        o = o if isinstance(o, I) else I(o)
        return I(_add_lo(s.lo, o.lo), _add_hi(s.hi, o.hi))
    __radd__ = __add__

    def __neg__(s):
        return I(-s.hi, -s.lo)

    def __sub__(s, o):
        o = o if isinstance(o, I) else I(o)
        return s + (-o)

    def __rsub__(s, o):
        return I(o) - s

    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o)
        pairs = ((s.lo, o.lo), (s.lo, o.hi), (s.hi, o.lo), (s.hi, o.hi))
        return I(min(_mul_lo(a, b) for a, b in pairs), max(_mul_hi(a, b) for a, b in pairs))
    __rmul__ = __mul__

    def __repr__(s):
        return f'[{s.lo:.6g}, {s.hi:.6g}]'

    def hull(s, o):
        return I(min(s.lo, o.lo), max(s.hi, o.hi))


def _cos_bounds(t):
    if t == 0.0: return 1.0, 1.0
    c = math.cos(t)
    return _dn(_dn(c)), _up(_up(c))


def _sin_bounds(t):
    if t == 0.0: return 0.0, 0.0
    c = math.sin(t)
    return _dn(_dn(c)), _up(_up(c))


def icos(x):
    """cos over [lo, hi] (radians); interior extrema at k pi.  |x| < 2 pi assumed for the extremum test."""
    a, b = x.lo, x.hi
    if b - a >= 2 * math.pi: return I(-1, 1)
    la, ha = _cos_bounds(a); lb, hb = _cos_bounds(b)
    lo, hi = min(la, lb), max(ha, hb)
    k0, k1 = math.ceil(_dn(a) / math.pi), math.floor(_up(b) / math.pi)
    for k in range(k0, k1 + 1):
        v = 1.0 if k % 2 == 0 else -1.0
        lo, hi = min(lo, v), max(hi, v)
    return I(max(lo, -1.0), min(hi, 1.0))


def isin(x):
    """sin over [lo, hi]; interior extrema at pi/2 + k pi."""
    a, b = x.lo, x.hi
    if b - a >= 2 * math.pi: return I(-1, 1)
    la, ha = _sin_bounds(a); lb, hb = _sin_bounds(b)
    lo, hi = min(la, lb), max(ha, hb)
    k0, k1 = math.ceil((_dn(a) - math.pi / 2) / math.pi) - 1, math.floor((_up(b) - math.pi / 2) / math.pi) + 1
    for k in range(k0, k1 + 1):
        c = math.pi / 2 + k * math.pi              # candidate extremum (float approx; include if possibly inside)
        if _dn(a) <= _up(c) and _dn(c) <= _up(b):
            v = 1.0 if k % 2 == 0 else -1.0
            lo, hi = min(lo, v), max(hi, v)
    return I(max(lo, -1.0), min(hi, 1.0))


RHO = None   # I([0, rho0]): the ray parameter range, set by the caller (enables base-point enclosures)


def _meet(a, b):
    lo, hi = max(a.lo, b.lo), min(a.hi, b.hi)
    return I(lo, hi) if lo <= hi else a


class D:
    """value v over the box, gradient g (tuple of I), and v0: the value at the base point rho = 0 (same box in the
    non-radial variables).  With one AD variable rho in RHO, the value is enclosed by v and by v0 + RHO * g."""
    __slots__ = ('v', 'g', 'v0')

    def __init__(self, v, g, v0=None):
        self.v, self.g = v, tuple(g)
        self.v0 = v if v0 is None else v0

    def enc(self):
        if RHO is None or len(self.g) != 1: return self.v
        return _meet(self.v, self.v0 + RHO * self.g[0])

    @staticmethod
    def const(c, n):
        c = c if isinstance(c, I) else I(c)
        return D(c, [I(0)] * n, c)

    def __add__(s, o):
        if not isinstance(o, D): return D(s.v + o, s.g, s.v0 + o)
        return D(s.v + o.v, [a + b for a, b in zip(s.g, o.g)], s.v0 + o.v0)
    __radd__ = __add__

    def __neg__(s):
        return D(-s.v, [-a for a in s.g], -s.v0)

    def __sub__(s, o):
        return s + (-o)

    def __rsub__(s, o):
        return (-s) + o

    def __mul__(s, o):
        if not isinstance(o, D):
            return D(s.v * o, [a * o for a in s.g], s.v0 * o)
        sv, ov = s.enc(), o.enc()
        return D(sv * ov, [a * ov + sv * b for a, b in zip(s.g, o.g)], s.v0 * o.v0)
    __rmul__ = __mul__


def dcos(x):
    v = x.enc()
    return D(icos(v), [-isin(v) * a for a in x.g], icos(x.v0))


def dsin(x):
    v = x.enc()
    return D(isin(v), [icos(v) * a for a in x.g], isin(x.v0))


def dabs(x):
    """|x|: exact branch if the sign is known on the box (via the base-point enclosure too), else Clarke hull."""
    e = x.enc()
    a0 = I(abs(x.v0.lo), abs(x.v0.hi)) if (x.v0.lo >= 0 or x.v0.hi <= 0) else I(0, max(-x.v0.lo, x.v0.hi))
    if e.lo >= 0: return D(e, x.g, x.v0 if x.v0.lo >= 0 else a0)
    if e.hi <= 0: return D(-e, [-a for a in x.g], -x.v0 if x.v0.hi <= 0 else a0)
    v = I(0, max(-e.lo, e.hi))
    return D(v, [a * I(-1, 1) for a in x.g], a0)


def dfsq(x):
    """(|cos x| + |sin x|)/2."""
    c, s = dcos(x), dsin(x)
    return (dabs(c) + dabs(s)) * 0.5


def idiv(a, b):
    if b.lo <= 0 <= b.hi:
        raise ZeroDivisionError('interval division by an interval containing 0')
    qs = []
    for x in (a.lo, a.hi):
        for y in (b.lo, b.hi):
            q = x / y
            err = -math.fma(q, y, -x)                    # x - q y, exactly; true x/y = q + err/y
            e = err / y
            qs.append((q if e >= 0 else _dn(q), q if e <= 0 else _up(q)))
    return I(min(q[0] for q in qs), max(q[1] for q in qs))


def ddiv(x, y):
    """x / y for D numbers (y's value interval must exclude 0)."""
    if not isinstance(y, D):
        y = D.const(y, len(x.g))
    if not isinstance(x, D):
        x = D.const(x, len(y.g))
    xv, yv = x.enc(), y.enc()
    q = idiv(xv, yv)
    y2 = yv * yv
    return D(q, [idiv(a * yv - xv * b, y2) for a, b in zip(x.g, y.g)], idiv(x.v0, y.v0))
