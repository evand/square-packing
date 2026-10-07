#!/usr/bin/env python3
"""Independent exact mass evaluator for mixed v1 covers (segments + convex polygons + points).

Written from private/s12/tasks/line-cover/FORMAT.md only.  Pose: centre (cx, cy), angle theta with
u = tan(theta/2) rational, so c = (1-u^2)/(1+u^2), s = 2u/(1+u^2) are rational.
Q = { (cx,cy) + R(theta)(a,b) : |a|,|b| <= 1/2 }  (closed).
"""
from fractions import Fraction as F
import numpy as np

HALF = F(1, 2)


def load(path):
    tok = []
    for ln in open(path):
        tok += ln.split('#', 1)[0].split()
    it = iter(tok)
    nx = lambda: next(it)
    assert nx() == 'mixed' and int(nx()) == 1
    sn, sd = int(nx()), int(nx())
    D = int(nx()); W = int(nx())
    pts = []
    for _ in range(int(nx())):
        X, Y, w = int(nx()), int(nx()), int(nx())
        pts.append((F(X, D), F(Y, D), F(w, W)))
    segs = []
    for _ in range(int(nx())):
        X0, Y0, X1, Y1, w = (int(nx()) for _ in range(5))
        segs.append((F(X0, D), F(Y0, D), F(X1, D), F(Y1, D), F(w, W)))
    polys = []
    for _ in range(int(nx())):
        k = int(nx()); w = int(nx())
        V = [(F(int(nx()), D), F(int(nx()), D)) for _ in range(k)]
        polys.append((F(w, W), V))
    rest = list(it)
    assert not rest, rest
    return dict(s=F(sn, sd), pts=pts, segs=segs, polys=polys)


def area(V):
    n = len(V)
    return sum(V[i][0] * V[(i + 1) % n][1] - V[(i + 1) % n][0] * V[i][1] for i in range(n)) / 2


def clip(V, a, b, c):
    """keep a*x + b*y <= c"""
    out = []
    n = len(V)
    for i in range(n):
        P, Qp = V[i], V[(i + 1) % n]
        fp = a * P[0] + b * P[1] - c
        fq = a * Qp[0] + b * Qp[1] - c
        if fp <= 0: out.append(P)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((P[0] + t * (Qp[0] - P[0]), P[1] + t * (Qp[1] - P[1])))
    return out


def cs(u):
    u = F(u)
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d


def square(cx, cy, u):
    c, s = cs(u)
    V = []
    for a, b in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)):
        V.append((cx + c * a - s * b, cy + s * a + c * b))
    return V  # CCW


class Cover:
    def __init__(self, path):
        d = load(path)
        self.s = d['s']; self.pts = d['pts']; self.segs = d['segs']; self.polys = d['polys']
        self.polyinfo = []
        for w, V in self.polys:
            A = area(V)
            assert A > 0
            n = len(V)
            hp = []  # CCW polygon: interior is left of each edge: cross(e, p-P) >= 0
            for i in range(n):
                P, Qp = V[i], V[(i + 1) % n]
                ex, ey = Qp[0] - P[0], Qp[1] - P[1]
                # ex*(y-Py) - ey*(x-Px) >= 0  <=>  ey*x - ex*y <= ey*Px - ex*Py
                hp.append((ey, -ex, ey * P[0] - ex * P[1]))
            self.polyinfo.append((w / A, hp))
        self.total = sum(p[2] for p in self.pts) + sum(g[4] for g in self.segs) + sum(w for w, _ in self.polys)
        # float arrays
        S = np.array([[float(v) for v in g] for g in self.segs])
        self.fP0 = S[:, 0:2]; self.fd = S[:, 2:4] - S[:, 0:2]; self.fw = S[:, 4]

    def admissible(self, cx, cy, u):
        V = square(cx, cy, u)
        return all(0 <= x <= self.s and 0 <= y <= self.s for x, y in V)

    def mass(self, cx, cy, u, detail=False):
        cx, cy = F(cx), F(cy)
        c, s = cs(u)
        e = (abs(c) + abs(s)) / 2
        xl, xh, yl, yh = cx - e, cx + e, cy - e, cy + e
        tot = F(0)
        parts = []
        for (x0, y0, x1, y1, w) in self.segs:
            if max(x0, x1) < xl or min(x0, x1) > xh or max(y0, y1) < yl or min(y0, y1) > yh:
                continue
            px, py = x0 - cx, y0 - cy
            dx, dy = x1 - x0, y1 - y0
            a0 = c * px + s * py; da = c * dx + s * dy
            b0 = -s * px + c * py; db = -s * dx + c * dy
            lo, hi = F(0), F(1)
            ok = True
            for v0, dv in ((a0, da), (b0, db)):
                if dv == 0:
                    if abs(v0) > HALF: ok = False; break
                else:
                    t1 = (-HALF - v0) / dv; t2 = (HALF - v0) / dv
                    if t1 > t2: t1, t2 = t2, t1
                    lo = max(lo, t1); hi = min(hi, t2)
            if ok and hi > lo:
                tot += w * (hi - lo)
                if detail: parts.append(((x0, y0, x1, y1), w * (hi - lo)))
        for (X, Y, w) in self.pts:
            px, py = X - cx, Y - cy
            if abs(c * px + s * py) <= HALF and abs(-s * px + c * py) <= HALF:
                tot += w
        Vq = square(cx, cy, u)
        for dens, hp in self.polyinfo:
            V = Vq
            for a, b, cc in hp:
                V = clip(V, a, b, cc)
                if not V: break
            if len(V) >= 3:
                ar = area(V)
                tot += dens * ar
                if detail: parts.append(('poly', dens * ar))
        return (tot, parts) if detail else tot

    def fmass(self, cx, cy, th):
        """float mass for arrays cx, cy, th (same shape (n,)); returns (n,) array.  Heuristic."""
        cx = np.asarray(cx, float)[:, None]; cy = np.asarray(cy, float)[:, None]; th = np.asarray(th, float)[:, None]
        c, s = np.cos(th), np.sin(th)
        px = self.fP0[None, :, 0] - cx; py = self.fP0[None, :, 1] - cy
        dx = self.fd[None, :, 0]; dy = self.fd[None, :, 1]
        lo = np.zeros(px.shape); hi = np.ones(px.shape)
        for v0, dv in ((c * px + s * py, c * dx + s * dy), (-s * px + c * py, -s * dx + c * dy)):
            small = np.abs(dv) < 1e-300
            with np.errstate(divide='ignore', invalid='ignore'):
                t1 = (-0.5 - v0) / dv; t2 = (0.5 - v0) / dv
            tl = np.where(small, np.where(np.abs(v0) <= 0.5 + 1e-13, -np.inf, np.inf), np.minimum(t1, t2))
            th_ = np.where(small, np.where(np.abs(v0) <= 0.5 + 1e-13, np.inf, -np.inf), np.maximum(t1, t2))
            lo = np.maximum(lo, tl); hi = np.minimum(hi, th_)
        segm = (np.clip(hi - lo, 0, None) * self.fw[None, :]).sum(1)
        # polygon (assume the single axis box polygon U) : area via sampling-free clipping in float
        out = segm
        for dens, hp in self.polyinfo:
            fa = np.array([[float(a), float(b), float(cc)] for a, b, cc in hp])
            areas = np.array([_farea_clip(cx[i, 0], cy[i, 0], th[i, 0], fa) for i in range(cx.shape[0])])
            out = out + float(dens) * areas
        return out


def _farea_clip(cx, cy, th, fa):
    c, s = np.cos(th), np.sin(th)
    V = [(cx + c * a - s * b, cy + s * a + c * b) for a, b in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]
    for a, b, cc in fa:
        out = []
        n = len(V)
        for i in range(n):
            P, Qp = V[i], V[(i + 1) % n]
            fp = a * P[0] + b * P[1] - cc; fq = a * Qp[0] + b * Qp[1] - cc
            if fp <= 0: out.append(P)
            if (fp < 0 < fq) or (fq < 0 < fp):
                t = fp / (fp - fq)
                out.append((P[0] + t * (Qp[0] - P[0]), P[1] + t * (Qp[1] - P[1])))
        V = out
        if not V: return 0.0
    n = len(V)
    if n < 3: return 0.0
    return 0.5 * sum(V[i][0] * V[(i + 1) % n][1] - V[(i + 1) % n][0] * V[i][1] for i in range(n))


D4 = [lambda x, y, m: (x, y), lambda x, y, m: (m - x, y), lambda x, y, m: (x, m - y), lambda x, y, m: (m - x, m - y),
      lambda x, y, m: (y, x), lambda x, y, m: (m - y, x), lambda x, y, m: (y, m - x), lambda x, y, m: (m - y, m - x)]
