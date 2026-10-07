"""Independent evaluator for the zmx2 audit (written from FORMAT.md only, not from zmx2_tools.py).

mu(Q(c, theta)) for a mixed cover: exact with Fractions at rational (cx, cy, u), u = tan(theta/2);
float version (numpy) for searches.
"""
import sys
from fractions import Fraction as Fr
import math
import numpy as np


def load(path):
    t = []
    for l in open(path):
        t += l.split('#')[0].split()
    i = 0
    if t[0] == 'mixed':
        i = 2
    def nx():
        nonlocal i
        i += 1
        return int(t[i - 1])
    s_num, s_den, D, W = nx(), nx(), nx(), nx()
    pts = [(nx(), nx(), nx()) for _ in range(nx())]
    segs = [(nx(), nx(), nx(), nx(), nx()) for _ in range(nx())] if t[0] == 'mixed' else []
    return dict(s=Fr(s_num, s_den), D=D, W=W, pts=pts, segs=segs)


def write(path, cv):
    with open(path, 'w') as f:
        f.write('mixed 1\n%d %d\n%d\n%d\n%d\n' % (cv['s'].numerator, cv['s'].denominator, cv['D'], cv['W'], len(cv['pts'])))
        for p in cv['pts']:
            f.write('%d %d %d\n' % p)
        f.write('%d\n' % len(cv['segs']))
        for s in cv['segs']:
            f.write('%d %d %d %d %d\n' % s)
        f.write('0\n')


def inside(cv, x, y, u):
    """exact: is Q(x,y,u) inside [0,s]^2"""
    n = 1 + u * u
    w = (1 - u * u + 2 * u) / n
    s = cv['s']
    return w / 2 <= x <= s - w / 2 and w / 2 <= y <= s - w / 2


def mu(cv, x, y, u):
    """exact mu of the closed square centred (x,y), angle 2 atan u (0 <= u < 1).
    Membership: q in Q  iff  |<q-c, e1>| <= 1/2 and |<q-c, e2>| <= 1/2,
    e1 = (C, S), e2 = (-S, C)."""
    n = 1 + u * u
    C, S = (1 - u * u) / n, 2 * u / n
    D, W = cv['D'], cv['W']
    h = Fr(1, 2)
    tot = Fr(0)
    for (X, Y, w) in cv['pts']:
        if w == 0 or abs(X / D - float(x)) > 0.75 or abs(Y / D - float(y)) > 0.75:
            continue
        a, b = Fr(X, D) - x, Fr(Y, D) - y
        if abs(a * C + b * S) <= h and abs(-a * S + b * C) <= h:
            tot += w
    for (X0, Y0, X1, Y1, w) in cv['segs']:
        if w == 0:
            continue
        # parametrise p(t) = P0 + t (P1 - P0), t in [0,1]; each slab constraint is linear in t
        lo, hi = Fr(0), Fr(1)
        P0 = (Fr(X0, D) - x, Fr(Y0, D) - y)
        V = (Fr(X1 - X0, D), Fr(Y1 - Y0, D))
        for (e0, e1) in ((C, S), (-S, C)):
            a = P0[0] * e0 + P0[1] * e1
            b = V[0] * e0 + V[1] * e1
            if b == 0:
                if abs(a) > h:
                    lo, hi = Fr(1), Fr(0)
            else:
                t1, t2 = (-h - a) / b, (h - a) / b
                if t1 > t2:
                    t1, t2 = t2, t1
                lo, hi = max(lo, t1), min(hi, t2)
        if hi > lo:
            tot += w * (hi - lo)
    return tot / W


class FCover:
    """float vectorised evaluator"""
    def __init__(self, cv):
        D, W = cv['D'], cv['W']
        self.s = float(cv['s'])
        p = np.array(cv['pts'], dtype=float).reshape(-1, 3)
        self.px, self.py, self.pw = p[:, 0] / D, p[:, 1] / D, p[:, 2] / W
        g = np.array(cv['segs'], dtype=float).reshape(-1, 5)
        self.a0, self.b0 = g[:, 0] / D, g[:, 1] / D
        self.da, self.db = (g[:, 2] - g[:, 0]) / D, (g[:, 3] - g[:, 1]) / D
        self.sw = g[:, 4] / W

    def mu(self, x, y, th):
        C, S = math.cos(th), math.sin(th)
        a, b = self.px - x, self.py - y
        e = 1e-12
        m = self.pw[(np.abs(a * C + b * S) <= .5 + e) & (np.abs(-a * S + b * C) <= .5 + e)].sum()
        a0, b0 = self.a0 - x, self.b0 - y
        lo = np.zeros_like(a0)
        hi = np.ones_like(a0)
        for (e0, e1) in ((C, S), (-S, C)):
            A = a0 * e0 + b0 * e1
            B = self.da * e0 + self.db * e1
            with np.errstate(divide='ignore', invalid='ignore'):
                t1 = (-.5 - A) / B
                t2 = (.5 - A) / B
            z = np.abs(B) < 1e-15
            tl = np.where(z, np.where(np.abs(A) <= .5 + e, -np.inf, np.inf), np.minimum(t1, t2))
            th_ = np.where(z, np.where(np.abs(A) <= .5 + e, np.inf, -np.inf), np.maximum(t1, t2))
            lo = np.maximum(lo, tl)
            hi = np.minimum(hi, th_)
        return m + (self.sw * np.clip(hi - lo, 0, None)).sum()

    def adm(self, x, y, th):
        w = math.cos(th) + math.sin(th)
        return w / 2 <= x <= self.s - w / 2 and w / 2 <= y <= self.s - w / 2


def fr(v, den=10 ** 7):
    return Fr(v).limit_denominator(den)


if __name__ == '__main__':
    cv = load(sys.argv[1])
    x, y, u = (Fr(t) for t in sys.argv[2:5])
    m = mu(cv, x, y, u)
    print('mu', float(m), 'admissible', inside(cv, x, y, u))
