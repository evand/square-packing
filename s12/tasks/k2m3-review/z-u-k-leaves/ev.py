# Independent exact evaluator of mu(Q) for rotated closed unit squares on the box cover (Fractions).
from fractions import Fraction as F
HALF = F(1, 2)
def load(fn):
    toks = []
    for line in open(fn):
        toks += line.split('#')[0].split()
    it = iter(toks[2:]); nx = lambda: int(next(it))
    s = F(nx(), nx()); D = nx(); W = nx(); assert nx() == 0
    H, V = {}, {}
    for _ in range(nx()):
        X0, Y0, X1, Y1, w = [nx() for _ in range(5)]
        if Y0 == Y1:
            H.setdefault(F(Y0, D), []).append((F(min(X0, X1), D), F(max(X0, X1), D), F(w, W)))
        else:
            assert X0 == X1
            V.setdefault(F(X0, D), []).append((F(min(Y0, Y1), D), F(max(Y0, Y1), D), F(w, W)))
    assert nx() == 1
    k = nx(); w = F(nx(), W); pts = [(F(nx(), D), F(nx(), D)) for _ in range(k)]
    a = min(p[0] for p in pts); b = max(p[0] for p in pts)
    return s, H, V, a, b

def clip(poly, f):  # keep f(p) >= 0, f affine
    out = []
    n = len(poly)
    for i in range(n):
        P, Qp = poly[i], poly[(i + 1) % n]
        fp, fq = f(P), f(Qp)
        if fp >= 0: out.append(P)
        if (fp >= 0) != (fq >= 0) and fp != fq:
            t = fp / (fp - fq)
            R = (P[0] + t * (Qp[0] - P[0]), P[1] + t * (Qp[1] - P[1]))
            if fq != 0 and fp != 0: out.append(R)
            elif fp == 0 or fq == 0: pass
    return out

def area(poly):
    s = F(0)
    for i in range(len(poly)):
        x0, y0 = poly[i]; x1, y1 = poly[(i + 1) % len(poly)]
        s += x0 * y1 - x1 * y0
    return abs(s) / 2

class Cov:
    def __init__(self, fn):
        self.s, self.H, self.V, self.a, self.b = load(fn)
    def verts(self, cx, cy, u):
        C, S, N = 1 - u * u, 2 * u, 1 + u * u
        c, s = C / N, S / N
        return [(cx + x * c - y * s, cy + x * s + y * c) for x, y in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF))]
    def chord_h(self, cx, cy, u, eta):
        C, S, N = 1 - u * u, 2 * u, 1 + u * u
        dy = eta - cy
        lo = cx + (-N / 2 - dy * S) / C; hi = cx + (N / 2 - dy * S) / C
        if S == 0:
            if abs(dy) > HALF: return None
        else:
            if S > 0:
                lo = max(lo, cx + (dy * C - N / 2) / S); hi = min(hi, cx + (dy * C + N / 2) / S)
            else:
                lo = max(lo, cx + (dy * C + N / 2) / S); hi = min(hi, cx + (dy * C - N / 2) / S)
        return (lo, hi) if lo <= hi else None
    def chord_v(self, cx, cy, u, xi):
        C, S, N = 1 - u * u, 2 * u, 1 + u * u
        dx = xi - cx
        lo = cy + (dx * S - N / 2) / C; hi = cy + (dx * S + N / 2) / C
        if S == 0:
            if abs(dx) > HALF: return None
        else:
            if S > 0:
                lo = max(lo, cy + (-N / 2 - dx * C) / S); hi = min(hi, cy + (N / 2 - dx * C) / S)
            else:
                lo = max(lo, cy + (N / 2 - dx * C) / S); hi = min(hi, cy + (-N / 2 - dx * C) / S)
        return (lo, hi) if lo <= hi else None
    def parts(self, cx, cy, u):
        vs = self.verts(cx, cy, u)
        ys = [v[1] for v in vs]; xs = [v[0] for v in vs]
        a, b = self.a, self.b
        P = vs
        for f in (lambda p: p[0] - a, lambda p: b - p[0], lambda p: p[1] - a, lambda p: b - p[1]):
            P = clip(P, f)
            if not P: break
        lam = area(P) if len(P) >= 3 else F(0)
        segm = F(0)
        for eta, segs in self.H.items():
            if eta < min(ys) or eta > max(ys): continue
            ch = self.chord_h(cx, cy, u, eta)
            if ch is None: continue
            for p, q, w in segs:
                o = min(q, ch[1]) - max(p, ch[0])
                if o > 0: segm += w * o / (q - p)
        for xi, segs in self.V.items():
            if xi < min(xs) or xi > max(xs): continue
            ch = self.chord_v(cx, cy, u, xi)
            if ch is None: continue
            for p, q, w in segs:
                o = min(q, ch[1]) - max(p, ch[0])
                if o > 0: segm += w * o / (q - p)
        return lam, segm
    def mass(self, cx, cy, u):
        l, s = self.parts(cx, cy, u); return l + s
    def admissible(self, cx, cy, u):
        vs = self.verts(cx, cy, u)
        return all(0 <= v[0] <= self.s and 0 <= v[1] <= self.s for v in vs)
