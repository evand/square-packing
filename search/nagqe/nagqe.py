#!/usr/bin/env python3
"""Nagamochi's per-square scoring lemma (EJC 2005, Lemma 1) as an exact real-arithmetic decision problem.

Resource on the container [0,m]^2 (chelokot's formalisation, `NagamochiResource.measure`), with the layout made
parametric:
  * area of the inner box [1, m-1]^2, weight 1;
  * the four lines x, y in {1, m-1}, each over [1-e, m-1+e], density w;
  * Q points at the eight segment ends (1-e, 1), (1, 1-e), ...      weight q each;
  * P points (i, 1-d), (i, m-1+d), (1-d, i), (m-1+d, i), i = 2..m-2, weight p each.
Nagamochi: e = d = 1/10, w = 1/2, q = 9/20, p = 1/2 (total m^2 - 2).
A square of side lam > 1 scores the resource in its OPEN interior.  Lemma 1 claims score > 1 whenever the
closed square fits in [0,m]^2 and 1 < lam (chelokot's statement also has lam <= 101/100).

Pose: centre (X, Y), direction u = (c, s), v = (-s, c), c^2 + s^2 = 1.  By the square's own symmetry and the
container's D4 we may take 0 <= s <= c (angle in [0, 45deg]) and X, Y <= m/2: mirroring x -> m - x maps angle
a -> -a, and the quarter turns of the container keep the angle mod 90deg and move the centre to any quadrant.

`score_exact`: Fractions only (rational c, s), area by polygon clipping -- independent of the SMT encoding.
`counterexample_query`: z3 formula "exists a fitting pose with score <= 1", exact (NLSAT), area by Green.
"""
from fractions import Fraction as F
import sys

NAGAMOCHI = dict(e=F(1, 10), d=F(1, 10), w=F(1, 2), q=F(9, 20), p=F(1, 2))


def resources(m, e, d, w, q, p):
    """(points, lines): points = [(x, y, weight)], lines = [(axis, h, a, b, density)] with axis 'y' meaning y = h."""
    pts = []
    for (x, y) in [(1 - e, 1), (m - 1 + e, 1), (1 - e, m - 1), (m - 1 + e, m - 1),
                   (1, 1 - e), (1, m - 1 + e), (m - 1, 1 - e), (m - 1, m - 1 + e)]:
        pts.append((F(x), F(y), q))
    for i in range(2, m - 1):
        for (x, y) in [(i, 1 - d), (i, m - 1 + d), (1 - d, i), (m - 1 + d, i)]:
            pts.append((F(x), F(y), p))
    lines = [(ax, F(h), 1 - e, m - 1 + e, w) for ax in ('y', 'x') for h in (1, m - 1)]
    return pts, lines


def total_mass(m, e, d, w, q, p):
    return (m - 2) ** 2 + 4 * w * (m - 2 + 2 * e) + 8 * q + 4 * (m - 3) * p


# ---------------------------------------------------------------- exact evaluator (Fractions)

def _interval(k, r, half):
    """{t : |k t + r| < half} as (lo, hi) open interval; None = all of R; () = empty."""
    if k == 0:
        return None if abs(r) < half else ()
    a, b = (-half - r) / k, (half - r) / k
    return (min(a, b), max(a, b))


def _chord(X, Y, c, s, lam, ax, h, a, b):
    half = lam / 2
    if ax == 'y':   # points (t, h): u-coord (t-X)c + (h-Y)s, v-coord -(t-X)s + (h-Y)c
        ivs = [_interval(c, -X * c + (h - Y) * s, half), _interval(-s, X * s + (h - Y) * c, half)]
    else:           # points (h, t): u-coord (h-X)c + (t-Y)s, v-coord -(h-X)s + (t-Y)c
        ivs = [_interval(s, (h - X) * c - Y * s, half), _interval(c, -(h - X) * s - Y * c, half)]
    lo, hi = a, b
    for iv in ivs:
        if iv == ():
            return F(0)
        if iv is not None:
            lo, hi = max(lo, iv[0]), min(hi, iv[1])
    return max(F(0), hi - lo)


def vertices(X, Y, c, s, lam):
    h = lam / 2
    return [(X + sx * h * c - sy * h * s, Y + sx * h * s + sy * h * c)
            for (sx, sy) in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]   # counter-clockwise


def _clip(poly, inside, inter):
    out = []
    for i in range(len(poly)):
        P, Q = poly[i - 1], poly[i]
        if inside(Q):
            if not inside(P):
                out.append(inter(P, Q))
            out.append(Q)
        elif inside(P):
            out.append(inter(P, Q))
    return out


def _box_area(poly, lo, hi):
    for coord in (0, 1):
        for bound, sign in ((lo, 1), (hi, -1)):
            def inside(P, coord=coord, bound=bound, sign=sign):
                return sign * (P[coord] - bound) >= 0

            def inter(P, Q, coord=coord, bound=bound):
                t = (bound - P[coord]) / (Q[coord] - P[coord])
                return (P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1]))
            poly = _clip(poly, inside, inter)
            if not poly:
                return F(0)
    return abs(sum(poly[i - 1][0] * poly[i][1] - poly[i][0] * poly[i - 1][1] for i in range(len(poly)))) / 2


def fits(m, X, Y, c, s, lam):
    return all(0 <= x <= m and 0 <= y <= m for (x, y) in vertices(X, Y, c, s, lam))


def score_exact(m, X, Y, c, s, lam, e, d, w, q, p, parts=False):
    assert c * c + s * s == 1
    pts, lines = resources(m, e, d, w, q, p)
    half = lam / 2
    pt = sum((wt for (x, y, wt) in pts
              if abs((x - X) * c + (y - Y) * s) < half and abs(-(x - X) * s + (y - Y) * c) < half), F(0))
    ln = sum((dens * _chord(X, Y, c, s, lam, ax, h, a, b) for (ax, h, a, b, dens) in lines), F(0))
    ar = _box_area(vertices(X, Y, c, s, lam), F(1), F(m - 1))
    return (ar, ln, pt) if parts else ar + ln + pt


# ---------------------------------------------------------------- z3 encoding

def counterexample_query(m, e, d, w, q, p, tilted=True, lam_max=F(101, 100), quadrant=True, bound=1):
    """z3 Solver whose satisfiability = existence of a fitting pose (angle in (0,45deg] if tilted, else 0)
    with 1 < lam <= lam_max and score <= 1.  Returns (solver, vars)."""
    import z3
    R = lambda x: z3.RealVal(str(F(x)))
    X, Y, c, s, lam = z3.Reals('X Y c s lam')
    S = z3.Solver()
    S.add(lam > 1, lam <= R(lam_max))
    if tilted:
        S.add(c * c + s * s == 1, s > 0, s <= c)
    else:
        S.add(c == 1, s == 0)
    hx = lam / 2 * (c + s)
    S.add(X - hx >= 0, X + hx <= m, Y - hx >= 0, Y + hx <= m)
    if quadrant:
        S.add(X <= R(F(m, 2)), Y <= R(F(m, 2)))
    half = lam / 2

    def Max(*a):
        r = a[0]
        for x in a[1:]:
            r = z3.If(x > r, x, r)
        return r

    def Min(*a):
        r = a[0]
        for x in a[1:]:
            r = z3.If(x < r, x, r)
        return r

    pts, lines = resources(m, e, d, w, q, p)
    terms = []
    for (x, y, wt) in pts:
        x, y = R(x), R(y)
        inside = z3.And((x - X) * c + (y - Y) * s < half, (x - X) * c + (y - Y) * s > -half,
                        -(x - X) * s + (y - Y) * c < half, -(x - X) * s + (y - Y) * c > -half)
        terms.append(z3.If(inside, R(wt), R(0)))

    def chord(ax, h, a, b):
        h, a, b = R(h), R(a), R(b)
        if tilted:   # both coefficients nonzero; interval bounds as in _interval
            if ax == 'y':
                k1, r1, k2, r2 = c, -X * c + (h - Y) * s, -s, X * s + (h - Y) * c
            else:
                k1, r1, k2, r2 = s, (h - X) * c - Y * s, c, -(h - X) * s - Y * c
            # k1 > 0 always; k2 < 0 for 'y', > 0 for 'x'
            lo1, hi1 = (-half - r1) / k1, (half - r1) / k1
            if ax == 'y':
                lo2, hi2 = (half - r2) / k2, (-half - r2) / k2
            else:
                lo2, hi2 = (-half - r2) / k2, (half - r2) / k2
            L, U = Max(a, lo1, lo2), Min(b, hi1, hi2)
            return z3.If(U > L, U - L, R(0))
        if ax == 'y':
            ok = z3.And(h - Y < half, Y - h < half)
            L, U = Max(a, X - half), Min(b, X + half)
        else:
            ok = z3.And(h - X < half, X - h < half)
            L, U = Max(a, Y - half), Min(b, Y + half)
        return z3.If(z3.And(ok, U > L), U - L, R(0))

    for (ax, h, a, b, dens) in lines:
        terms.append(R(dens) * chord(ax, h, a, b))

    lo, hi = R(1), R(m - 1)
    if tilted:   # Green: area = sum over clipped square edges of int x dy + boundary chords of x = m-1 (up), x = 1 (down)
        hh = lam / 2
        V = [(X - hh * c + hh * s, Y - hh * s - hh * c), (X + hh * c + hh * s, Y + hh * s - hh * c),
             (X + hh * c - hh * s, Y + hh * s + hh * c), (X - hh * c - hh * s, Y - hh * s + hh * c)]
        D = [(c, s), (-s, c), (-c, -s), (s, -c)]
        sx = [1, -1, -1, 1]   # sign of dx per edge
        sy = [1, 1, -1, -1]
        area = R(0)
        for k in range(4):
            (vx, vy), (dx, dy) = V[k], D[k]
            ax1, ax2 = (lo - vx) / dx, (hi - vx) / dx
            ay1, ay2 = (lo - vy) / dy, (hi - vy) / dy
            lx, ux = (ax1, ax2) if sx[k] > 0 else (ax2, ax1)
            ly, uy = (ay1, ay2) if sy[k] > 0 else (ay2, ay1)
            t0, t1 = Max(R(0), lx, ly), Min(lam, ux, uy)
            area = area + z3.If(t1 > t0, dy * (vx * (t1 - t0) + dx * (t1 * t1 - t0 * t0) / 2), R(0))
        area = area + R(m - 1) * chord('x', m - 1, 1, m - 1) - chord('x', 1, 1, m - 1)
    else:
        ox = Min(hi, X + half) - Max(lo, X - half)
        oy = Min(hi, Y + half) - Max(lo, Y - half)
        area = z3.If(z3.And(ox > 0, oy > 0), ox * oy, R(0))
    terms.append(area)
    score = z3.Sum(terms)
    if bound is not None:
        S.add(score <= R(bound))
    return S, (X, Y, c, s, lam), score


def counterexample_query_df(m, e, d, w, q, p, lam_max=F(101, 100), quadrant=True, bound=1):
    """Division-free tilted version (0 < s <= c): the inequality score <= bound multiplied by K = c^2 s^2 > 0.
    Chord bounds are scaled by c*s, edge parameters by c*s; every atom is a polynomial in (X, Y, c, s, lam).
    Returns (solver, vars, scaled_score, K) with the assertion  scaled_score <= bound*K  added unless bound is None."""
    import z3
    R = lambda x: z3.RealVal(str(F(x)))
    X, Y, c, s, lam = z3.Reals('X Y c s lam')
    S = z3.Solver()
    S.add(lam > 1, lam <= R(lam_max), c * c + s * s == 1, s > 0, s <= c)
    hx = lam / 2 * (c + s)
    S.add(X - hx >= 0, X + hx <= m, Y - hx >= 0, Y + hx <= m)
    if quadrant:
        S.add(X <= R(F(m, 2)), Y <= R(F(m, 2)))
    half = lam / 2
    cs = c * s
    K = cs * cs

    def Max(*a):
        r = a[0]
        for x in a[1:]:
            r = z3.If(x > r, x, r)
        return r

    def Min(*a):
        r = a[0]
        for x in a[1:]:
            r = z3.If(x < r, x, r)
        return r

    pts, lines = resources(m, e, d, w, q, p)
    terms = []
    for (x, y, wt) in pts:
        x, y = R(x), R(y)
        uu, vv = (x - X) * c + (y - Y) * s, -(x - X) * s + (y - Y) * c
        terms.append(z3.If(z3.And(uu < half, uu > -half, vv < half, vv > -half), R(wt) * K, R(0)))

    def chord_cs(ax, h, a, b):
        """c*s*chord, polynomial pieces.  Bounds t = num/k scaled: t*cs = num*(cs/k)."""
        h, a, b = R(h), R(a), R(b)
        if ax == 'y':   # k1 = c, k2 = -s
            r1, r2 = -X * c + (h - Y) * s, X * s + (h - Y) * c
            lo1, hi1 = (-half - r1) * s, (half - r1) * s          # /c * cs = *s
            lo2, hi2 = -(half - r2) * c, -(-half - r2) * c        # /(-s) * cs = *(-c)
        else:           # k1 = s, k2 = c
            r1, r2 = (h - X) * c - Y * s, -(h - X) * s - Y * c
            lo1, hi1 = (-half - r1) * c, (half - r1) * c          # /s * cs = *c
            lo2, hi2 = (-half - r2) * s, (half - r2) * s          # /c * cs = *s
        L, U = Max(a * cs, lo1, lo2), Min(b * cs, hi1, hi2)
        return z3.If(U > L, U - L, R(0))

    for (ax, h, a, b, dens) in lines:
        terms.append(R(dens) * cs * chord_cs(ax, h, a, b))

    lo, hi = R(1), R(m - 1)
    hh = lam / 2
    V = [(X - hh * c + hh * s, Y - hh * s - hh * c), (X + hh * c + hh * s, Y + hh * s - hh * c),
         (X + hh * c - hh * s, Y + hh * s + hh * c), (X - hh * c - hh * s, Y - hh * s + hh * c)]
    D = [(c, s), (-s, c), (-c, -s), (s, -c)]
    # cs/dx and cs/dy as polynomials, per edge
    cs_over = {0: (s, c), 1: (-c, s), 2: (-s, -c), 3: (c, -s)}
    sx = [1, -1, -1, 1]
    sy = [1, 1, -1, -1]
    area_K = R(0)
    for k in range(4):
        (vx, vy), (dx, dy) = V[k], D[k]
        fx, fy = cs_over[k]
        ax1, ax2 = (lo - vx) * fx, (hi - vx) * fx
        ay1, ay2 = (lo - vy) * fy, (hi - vy) * fy
        lx, ux = (ax1, ax2) if sx[k] > 0 else (ax2, ax1)
        ly, uy = (ay1, ay2) if sy[k] > 0 else (ay2, ay1)
        T0, T1 = Max(R(0), lx, ly), Min(lam * cs, ux, uy)
        # K * dy*(vx*(t1-t0) + dx*(t1^2-t0^2)/2) with t = T/cs:  dy*(vx*cs*(T1-T0) + dx*(T1^2-T0^2)/2)
        area_K = area_K + z3.If(T1 > T0, dy * (vx * cs * (T1 - T0) + dx * (T1 * T1 - T0 * T0) / 2), R(0))
    area_K = area_K + cs * (R(m - 1) * chord_cs('x', m - 1, 1, m - 1) - chord_cs('x', 1, 1, m - 1))
    terms.append(area_K)
    score_K = z3.Sum(terms)
    if bound is not None:
        S.add(score_K <= R(bound) * K)
    return S, (X, Y, c, s, lam), score_K, K


def selftest():
    # chelokot's counterexample (docs/nagamochi-score-counterexample.md): bottom vertex A, side lam, u = (c, s)
    lam, c, s = F(10001, 10000), F(80, 1601), F(1599, 1601)
    A = (F(10419467, 5330000), F(0))
    X, Y = A[0] + lam / 2 * (c - s), A[1] + lam / 2 * (s + c)
    assert fits(4, X, Y, c, s, lam)
    sc = score_exact(4, X, Y, c, s, lam, **NAGAMOCHI)
    assert sc == F(25009470849041, 25584000000000), sc
    assert total_mass(4, **NAGAMOCHI) == 14
    print('selftest ok: chelokot counterexample scores', sc, '=', float(sc))


if __name__ == '__main__':
    selftest()
