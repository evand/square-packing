"""Independent mass evaluator for the qx2 covers (lemmaE-code review).  Exact (Fraction) and float versions.
mu(Q) = area(Q ∩ [a,b]^2) + sum over axis lines of the segment mass inside the closed chord Q ∩ line.
Q(c,u) = c + R_theta [-1/2,1/2]^2, cos = (1-u^2)/(1+u^2), sin = 2u/(1+u^2).  Written from scratch (polygon clipping,
edge/line intersection); does not use qx2_zm's option formulas or Profile."""
import sys, os
from fractions import Fraction as F
sys.path.insert(0, os.path.expanduser('~/math/square-packing/public/s12/search'))
import zm_mixed as ZM, mixed_cover as MC

H = F(1, 2)


def verts(cx, cy, u, one=None):
    if one is None:
        c = (1 - u * u) / (1 + u * u); s = 2 * u / (1 + u * u); h = H
    else:
        c = (1 - u * u) / (1 + u * u); s = 2 * u / (1 + u * u); h = 0.5
    out = []
    for (x, y) in ((-h, -h), (h, -h), (h, h), (-h, h)):
        out.append((cx + x * c - y * s, cy + x * s + y * c))
    return out


def clip(P, axis, val, keep_ge):
    out = []
    n = len(P)
    for i in range(n):
        p, q = P[i], P[(i + 1) % n]
        pin = (p[axis] >= val) if keep_ge else (p[axis] <= val)
        qin = (q[axis] >= val) if keep_ge else (q[axis] <= val)
        if pin: out.append(p)
        if pin != qin:
            t = (val - p[axis]) / (q[axis] - p[axis])
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def area(P):
    if len(P) < 3: return 0
    s = 0
    for i in range(len(P)):
        x0, y0 = P[i]; x1, y1 = P[(i + 1) % len(P)]
        s += x0 * y1 - x1 * y0
    return abs(s) / 2


def chord(P, axis, val):
    """range of the other coordinate over P ∩ {coord[axis] = val} (closed), or None"""
    o = 1 - axis
    pts = []
    n = len(P)
    for i in range(n):
        p, q = P[i], P[(i + 1) % n]
        if p[axis] == val: pts.append(p[o])
        if (p[axis] - val) * (q[axis] - val) < 0:
            t = (val - p[axis]) / (q[axis] - p[axis])
            pts.append(p[o] + t * (q[o] - p[o]))
    if not pts: return None
    return min(pts), max(pts)


class Mass:
    def __init__(self, cov, a, b):
        self.a, self.b = a, b
        self.lines = []
        for key, L in cov.lines.items():
            self.lines.append((key[0], key[1], [(t0, t1, d) for t0, t1, d in L['segs']],
                               float(key[1]), [(float(t0), float(t1), float(d)) for t0, t1, d in L['segs']]))

    def exact(self, cx, cy, u):
        P = verts(cx, cy, u)
        Q = P
        for ax, v, ge in ((0, self.a, True), (0, self.b, False), (1, self.a, True), (1, self.b, False)):
            Q = clip(Q, ax, v, ge)
            if not Q: break
        m = area(Q) if Q else F(0)
        for o, p, segs, _, _ in self.lines:
            ch = chord(P, 1 if o == 'H' else 0, p)
            if ch is None: continue
            lo, hi = ch
            for t0, t1, d in segs:
                ov = min(t1, hi) - max(t0, lo)
                if ov > 0: m += d * ov
        return m

    def flt(self, cx, cy, u):
        P = verts(cx, cy, u, one=True)
        Q = P
        for ax, v, ge in ((0, float(self.a), True), (0, float(self.b), False), (1, float(self.a), True), (1, float(self.b), False)):
            Q = clip(Q, ax, v, ge)
            if not Q: break
        m = area(Q) if Q else 0.0
        xs = [p[0] for p in P]; ys = [p[1] for p in P]
        for o, _, _, pf, segs in self.lines:
            if o == 'H':
                if pf < min(ys) or pf > max(ys): continue
                ch = chord(P, 1, pf)
            else:
                if pf < min(xs) or pf > max(xs): continue
                ch = chord(P, 0, pf)
            if ch is None: continue
            lo, hi = ch
            for t0, t1, d in segs:
                ov = min(t1, hi) - max(t0, lo)
                if ov > 0: m += d * ov
        return m


def load(path=os.path.expanduser('~/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')):
    cov = ZM.Cover(MC.load(path))
    import qx2_zm as QZ
    a, b = QZ.u_square(cov)
    return cov, a, b
