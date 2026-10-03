"""Independent exact geometry (Fractions) for the Lemma E review.  Nothing imported from public/."""
from fractions import Fraction as F
import random

HALF = F(1, 2)


def trig(u):
    N = 1 + u * u
    return (1 - u * u) / N, 2 * u / N          # cos, sin


def square(cx, cy, u):
    """vertices of Q(c,u), CCW: c + R_theta(+-1/2, +-1/2)"""
    c, s = trig(u)
    out = []
    for (x, y) in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)):
        out.append((cx + x * c - y * s, cy + x * s + y * c))
    return out


def clip(poly, a, b, k):
    """keep a x + b y + k <= 0"""
    out = []
    n = len(poly)
    for i in range(n):
        P = poly[i]; Qp = poly[(i + 1) % n]
        fp = a * P[0] + b * P[1] + k; fq = a * Qp[0] + b * Qp[1] + k
        if fp <= 0: out.append(P)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Qp[0] - P[0]), P[1] + t * (Qp[1] - P[1])))
    return out


def area(poly):
    if len(poly) < 3: return F(0)
    s = F(0)
    for i in range(len(poly)):
        x0, y0 = poly[i]; x1, y1 = poly[(i + 1) % len(poly)]
        s += x0 * y1 - x1 * y0
    return abs(s) / 2


def chord(poly, orient, pos):
    """[lo, hi] of the intersection of the convex polygon with the line y = pos ('H', param x) or x = pos ('V'),
    computed from the polygon edges (independent of the option table), or None"""
    pts = []
    n = len(poly)
    for i in range(n):
        P = poly[i]; Qp = poly[(i + 1) % n]
        if orient == 'H':
            p0, p1, q0, q1 = P[1], Qp[1], P[0], Qp[0]
        else:
            p0, p1, q0, q1 = P[0], Qp[0], P[1], Qp[1]
        if p0 == pos: pts.append(q0)
        if (p0 - pos) * (p1 - pos) < 0:
            t = (pos - p0) / (p1 - p0)
            pts.append(q0 + t * (q1 - q0))
    if not pts: return None
    return min(pts), max(pts)


class Prof:
    """piecewise-uniform measure on a line: breakpoints b[0..n-1], densities rho[0..n] (rho[0] = rho[n] = 0)"""

    def __init__(self, b, rho):
        self.b = list(b); self.rho = list(rho)
        assert len(self.rho) == len(self.b) + 1 and self.rho[0] == 0 and self.rho[-1] == 0
        self.Gb = [F(0)]
        for i in range(1, len(self.b)):
            self.Gb.append(self.Gb[-1] + self.rho[i] * (self.b[i] - self.b[i - 1]))

    def G(self, t):
        if t <= self.b[0]: return F(0)
        for i in range(1, len(self.b)):
            if t <= self.b[i]:
                return self.Gb[i - 1] + self.rho[i] * (t - self.b[i - 1])
        return self.Gb[-1]


def mass_line(poly, orient, pos, prof):
    ch = chord(poly, orient, pos)
    if ch is None: return F(0)
    return max(F(0), prof.G(ch[1]) - prof.G(ch[0]))


def rand_prof(rng, lo, hi, nb, den=20):
    pts = sorted(set(F(rng.randint(int(lo * den), int(hi * den)), den) for _ in range(nb)))
    if len(pts) < 2: pts = [F(int(lo * den), den), F(int(hi * den), den)]
    rho = [F(0)] + [F(rng.choice([0, 0, 1, 2, 3, 5]), 4) for _ in range(len(pts) - 1)] + [F(0)]
    return Prof(pts, rho)
