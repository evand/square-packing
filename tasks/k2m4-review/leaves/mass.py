"""Own exact mass evaluator for the k2m4 box cover (no public/ code imported).
mu(Q) for the closed unit square Q(c, u): centre c, theta = 2 atan u (u rational, cos = (1-u^2)/(1+u^2), sin = 2u/(1+u^2)).
Lines: chord of Q (exact, closed) intersected with the line's piecewise-uniform profile.  U: exact polygon clipping."""
from fractions import Fraction as F
import bisect, math
import sys
sys.path.insert(0, '/home/evand/math/square-packing/private/s12/tasks/k2m4-review/leaves')
from cover_check import parse, line_profiles, canon, BOX

M, _pts, SEGS, POLYS = parse(BOX)
A_, B_ = F(14, 5), F(31, 5)
PROF = {k: canon(v) for k, v in line_profiles(SEGS).items()}


class G:
    def __init__(self, prof):
        self.lo = [a for a, b, r in prof]; self.hi = [b for a, b, r in prof]; self.r = [r for a, b, r in prof]
        self.cum = [F(0)]
        for a, b, r in prof: self.cum.append(self.cum[-1] + r * (b - a))
        self.lof = [float(x) for x in self.lo]; self.hif = [float(x) for x in self.hi]; self.rf = [float(x) for x in self.r]
        self.cumf = [float(x) for x in self.cum]

    def __call__(self, t):  # mass on (-inf, t]
        i = bisect.bisect_right(self.lo, t) - 1
        if i < 0: return F(0)
        if t >= self.hi[i]: return self.cum[i + 1]
        return self.cum[i] + self.r[i] * (t - self.lo[i])

    def f(self, t):
        i = bisect.bisect_right(self.lof, t) - 1
        if i < 0: return 0.0
        if t >= self.hif[i]: return self.cumf[i + 1]
        return self.cumf[i] + self.rf[i] * (t - self.lof[i])


GS = {k: G(v) for k, v in PROF.items()}
HL = sorted((k[1], g) for k, g in GS.items() if k[0] == 'H')
VL = sorted((k[1], g) for k, g in GS.items() if k[0] == 'V')
HLf = [(float(p), g) for p, g in HL]; VLf = [(float(p), g) for p, g in VL]


def trig(u):
    N = 1 + u * u
    return (1 - u * u) / N, 2 * u / N


def square(cx, cy, u):
    c, s = trig(u)
    out = []
    for X, Y in ((-1, -1), (1, -1), (1, 1), (-1, 1)):   # ccw
        X = F(X, 2); Y = F(Y, 2)
        out.append((cx + c * X - s * Y, cy + s * X + c * Y))
    return out


def chord(V, axis, p):
    """closed chord of the convex polygon V with the line {coord_axis = p}, as (t0, t1) in the other coord, or None"""
    o = 1 - axis
    ts = []
    n = len(V)
    for i in range(n):
        P, Q = V[i], V[(i + 1) % n]
        a, b = P[axis], Q[axis]
        if a == b:
            if a == p: ts += [P[o], Q[o]]
            continue
        if min(a, b) <= p <= max(a, b):
            ts.append(P[o] + (p - a) * (Q[o] - P[o]) / (b - a))
    if not ts: return None
    return min(ts), max(ts)


def clip(V, axis, val, keep_ge):
    out = []
    n = len(V)
    inside = lambda P: (P[axis] >= val) if keep_ge else (P[axis] <= val)
    for i in range(n):
        P, Q = V[i], V[(i + 1) % n]
        ip, iq = inside(P), inside(Q)
        if ip: out.append(P)
        if ip != iq:
            t = (val - P[axis]) / (Q[axis] - P[axis])
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    return out


def area(V):
    if len(V) < 3: return F(0)
    return abs(sum(V[i][0] * V[(i + 1) % len(V)][1] - V[(i + 1) % len(V)][0] * V[i][1] for i in range(len(V)))) / 2


def mass_poly(V):
    """mu of a convex polygon region (closed)."""
    xs = [v[0] for v in V]; ys = [v[1] for v in V]
    xlo, xhi, ylo, yhi = min(xs), max(xs), min(ys), max(ys)
    tot = F(0)
    for p, g in HL:
        if ylo <= p <= yhi:
            ch = chord(V, 1, p)
            if ch: tot += g(ch[1]) - g(ch[0]) + (point_mass_closed(g, ch[0]))
    for p, g in VL:
        if xlo <= p <= xhi:
            ch = chord(V, 0, p)
            if ch: tot += g(ch[1]) - g(ch[0]) + (point_mass_closed(g, ch[0]))
    W = V
    for ax, val, ge in ((0, A_, True), (0, B_, False), (1, A_, True), (1, B_, False)):
        W = clip(W, ax, val, ge)
        if not W: break
    tot += area(W)
    return tot


def point_mass_closed(g, t):
    return F(0)   # line measures have no atoms


def mass(cx, cy, u):
    return mass_poly(square(cx, cy, u))


# ---------------------------------------------------------------- floats (for local search only)
def mass_f(cx, cy, th):
    c, s = math.cos(th), math.sin(th)
    V = [(cx + c * X - s * Y, cy + s * X + c * Y) for X, Y in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]
    xs = [v[0] for v in V]; ys = [v[1] for v in V]
    tot = 0.0
    for lines, axis, lo, hi in ((HLf, 1, min(ys), max(ys)), (VLf, 0, min(xs), max(xs))):
        o = 1 - axis
        for p, g in lines:
            if not (lo <= p <= hi): continue
            ts = []
            for i in range(4):
                P, Q = V[i], V[(i + 1) % 4]
                a, b = P[axis], Q[axis]
                if a == b: continue
                if min(a, b) <= p <= max(a, b): ts.append(P[o] + (p - a) * (Q[o] - P[o]) / (b - a))
            if ts: tot += g.f(max(ts)) - g.f(min(ts))
    W = V
    for ax, val, ge in ((0, 2.8, True), (0, 6.2, False), (1, 2.8, True), (1, 6.2, False)):
        W = clip(W, ax, val, ge)
        if not W: break
    tot += float(area(W)) if False else (abs(sum(W[i][0] * W[(i + 1) % len(W)][1] - W[(i + 1) % len(W)][0] * W[i][1]
                                                for i in range(len(W)))) / 2 if len(W) >= 3 else 0.0)
    return tot


def admissible(cx, cy, u):
    c, s = trig(u)
    h = (abs(c) + abs(s)) / 2
    return h <= cx <= M - h and h <= cy <= M - h
