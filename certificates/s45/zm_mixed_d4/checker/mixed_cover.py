#!/usr/bin/env python3
"""Reader for mixed covers (points + segments + polygons), format v1.

Format: tasks/line-cover/FORMAT.md (copied in search/LINE_COVER.md).  Contract (shared with zm_mixed.py):

    load(path)  -> dict(s=Fraction, s_num, s_den, D, W, points=[(X,Y,w)], segments=[(X0,Y0,X1,Y1,w)],
                        polygons=[(w, [(X1,Y1),...,(Xk,Yk)])])          all exact ints
    total(cover) -> Fraction                                           sum of all masses / W
    validate(cover)                                                    raises ValueError if ill-formed
    mass_in_square_float(cover, cx, cy, theta, tol=1e-12) -> float     mu(Q), Q the closed unit square with
                                                                       centre (cx,cy), angle theta (radians)

Closed semantics in the float evaluator: a point / segment piece within `tol` of Q counts (a segment lying on
an edge of Q counts in full).  This is a heuristic evaluator, not a proof.

Also: write(path, cover, comment=None), and a vectorised helper `MixedEval` for many poses.
"""
from fractions import Fraction
import math
import sys

import numpy as np


def _tokens(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.split('#', 1)[0]
            out.extend(line.split())
    return out


def load(path):
    t = _tokens(path)
    pos = 0

    def nxt():
        nonlocal pos
        if pos >= len(t):
            raise ValueError('unexpected end of file')
        v = t[pos]
        pos += 1
        return v

    def ni():
        v = nxt()
        try:
            return int(v)
        except ValueError:
            raise ValueError('expected integer, got %r' % v)

    if t and t[0] == 'mixed':
        nxt()
        ver = ni()
        if ver != 1:
            raise ValueError('unsupported mixed version %d' % ver)
        plain = False
    else:
        plain = True
    s_num, s_den = ni(), ni()
    D = ni()
    W = ni()
    npnt = ni()
    points = [(ni(), ni(), ni()) for _ in range(npnt)]
    segments, polygons = [], []
    if not plain:
        ns = ni()
        segments = [(ni(), ni(), ni(), ni(), ni()) for _ in range(ns)]
        npg = ni()
        for _ in range(npg):
            k = ni()
            if k < 3:
                raise ValueError('polygon with k=%d < 3 vertices' % k)
            w = ni()
            polygons.append((w, [(ni(), ni()) for _ in range(k)]))
    if pos != len(t):
        raise ValueError('trailing data (%d tokens)' % (len(t) - pos))
    cover = dict(s_num=s_num, s_den=s_den, D=D, W=W, points=points, segments=segments, polygons=polygons)
    if s_den > 0:
        cover['s'] = Fraction(s_num, s_den)
    return cover


def total(cover):
    tw = sum(p[2] for p in cover['points']) + sum(q[4] for q in cover['segments']) \
        + sum(g[0] for g in cover['polygons'])
    return Fraction(tw, cover['W'])


def _area2(vs):
    a = 0
    for i in range(len(vs)):
        x0, y0 = vs[i]
        x1, y1 = vs[(i + 1) % len(vs)]
        a += x0 * y1 - x1 * y0
    return a


def validate(cover):
    """Raise ValueError unless the cover is well formed (exact integer checks)."""
    for k in ('s_num', 's_den', 'D', 'W'):
        if not isinstance(cover[k], int) or cover[k] <= 0:
            raise ValueError('%s must be a positive integer' % k)
    s_num, s_den, D = cover['s_num'], cover['s_den'], cover['D']
    if (s_num * D) % s_den:
        raise ValueError('s_den must divide s_num*D')
    S = s_num * D // s_den            # container side in coordinate units

    def inb(X, Y):
        return 0 <= X <= S and 0 <= Y <= S

    for (X, Y, w) in cover['points']:
        if not inb(X, Y):
            raise ValueError('point (%d,%d) outside container' % (X, Y))
        if w < 0:
            raise ValueError('negative weight')
    for (X0, Y0, X1, Y1, w) in cover['segments']:
        if not (inb(X0, Y0) and inb(X1, Y1)):
            raise ValueError('segment endpoint outside container')
        if (X0, Y0) == (X1, Y1):
            raise ValueError('degenerate (zero-length) segment')
        if w < 0:
            raise ValueError('negative weight')
    for (w, vs) in cover['polygons']:
        if len(vs) < 3:
            raise ValueError('polygon with < 3 vertices')
        if w < 0:
            raise ValueError('negative weight')
        for (X, Y) in vs:
            if not inb(X, Y):
                raise ValueError('polygon vertex outside container')
        k = len(vs)
        for i in range(k):
            (x0, y0), (x1, y1), (x2, y2) = vs[i], vs[(i + 1) % k], vs[(i + 2) % k]
            if (x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1) < 0:
                raise ValueError('polygon not convex CCW')
        if _area2(vs) <= 0:
            raise ValueError('degenerate or clockwise polygon')
    return True


# ----------------------------------------------------------------------------------------------- float evaluation

def _square_frame(cx, cy, theta):
    c, s = math.cos(theta), math.sin(theta)
    return c, s


def _seg_frac_in_square(ax, ay, bx, by, cx, cy, c, s, tol):
    """Fraction (by length) of segment a-b inside closed unit square centre (cx,cy), axes (c,s),(-s,c).
    Liang-Barsky in the square's frame, with slack tol."""
    # local coords
    ux0 = c * (ax - cx) + s * (ay - cy)
    uy0 = -s * (ax - cx) + c * (ay - cy)
    ux1 = c * (bx - cx) + s * (by - cy)
    uy1 = -s * (bx - cx) + c * (by - cy)
    t0, t1 = 0.0, 1.0
    h = 0.5 + tol
    for p0, d in ((ux0, ux1 - ux0), (uy0, uy1 - uy0)):
        # -h <= p0 + t d <= h
        if abs(d) < 1e-300:
            if p0 < -h or p0 > h:
                return 0.0
            continue
        ta = (-h - p0) / d
        tb = (h - p0) / d
        if ta > tb:
            ta, tb = tb, ta
        if ta > t0:
            t0 = ta
        if tb < t1:
            t1 = tb
        if t0 >= t1:
            return 0.0
    return t1 - t0


def _clip_poly(vs, c, s, cx, cy, tol):
    """Area of convex polygon vs (float list, world) intersected with the closed unit square."""
    # to local frame
    loc = [(c * (x - cx) + s * (y - cy), -s * (x - cx) + c * (y - cy)) for (x, y) in vs]
    h = 0.5 + tol
    for axis in (0, 1):
        for sign in (1.0, -1.0):
            out = []
            n = len(loc)
            if n == 0:
                return 0.0
            for i in range(n):
                P = loc[i]
                Qp = loc[(i + 1) % n]
                fp = h - sign * P[axis]
                fq = h - sign * Qp[axis]
                if fp >= 0:
                    out.append(P)
                if (fp >= 0) != (fq >= 0):
                    tt = fp / (fp - fq)
                    out.append((P[0] + tt * (Qp[0] - P[0]), P[1] + tt * (Qp[1] - P[1])))
            loc = out
    if len(loc) < 3:
        return 0.0
    a = 0.0
    for i in range(len(loc)):
        x0, y0 = loc[i]
        x1, y1 = loc[(i + 1) % len(loc)]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def mass_in_square_float(cover, cx, cy, theta, tol=1e-12):
    """mu(Q) for the closed unit square Q with centre (cx, cy) and angle theta (float, heuristic)."""
    D, W = float(cover['D']), float(cover['W'])
    c, s = math.cos(theta), math.sin(theta)
    h = 0.5 + tol
    tot = 0.0
    for (X, Y, w) in cover['points']:
        x, y = X / D - cx, Y / D - cy
        if abs(c * x + s * y) <= h and abs(-s * x + c * y) <= h:
            tot += w
    for (X0, Y0, X1, Y1, w) in cover['segments']:
        if w:
            tot += w * _seg_frac_in_square(X0 / D, Y0 / D, X1 / D, Y1 / D, cx, cy, c, s, tol)
    for (w, vs) in cover['polygons']:
        if w:
            fv = [(X / D, Y / D) for (X, Y) in vs]
            A = abs(_area2(vs)) / (2.0 * D * D)
            tot += w * _clip_poly(fv, c, s, cx, cy, tol) / A
    return tot / W


class MixedEval:
    """Vectorised float evaluation of mu(Q) over many poses.  Pieces given in floats:
    pts (n,2), pw (n,), segs (m,4), sw (m,), polys list of (k,2) arrays, gw (g,) (masses as floats)."""

    def __init__(self, pts=None, pw=None, segs=None, sw=None, polys=None, gw=None, tol=1e-12):
        self.pts = np.zeros((0, 2)) if pts is None else np.asarray(pts, float).reshape(-1, 2)
        self.pw = np.zeros(0) if pw is None else np.asarray(pw, float)
        self.segs = np.zeros((0, 4)) if segs is None else np.asarray(segs, float).reshape(-1, 4)
        self.sw = np.zeros(0) if sw is None else np.asarray(sw, float)
        self.polys = [] if polys is None else [np.asarray(p, float) for p in polys]
        self.gw = np.zeros(0) if gw is None else np.asarray(gw, float)
        self.tol = tol

    @classmethod
    def from_cover(cls, cover, tol=1e-12):
        D, W = float(cover['D']), float(cover['W'])
        P = cover['points']
        S = cover['segments']
        G = cover['polygons']
        return cls(np.array([[x / D, y / D] for x, y, _ in P]).reshape(-1, 2), np.array([w / W for *_, w in P]),
                   np.array([[a / D, b / D, c / D, d / D] for a, b, c, d, _ in S]).reshape(-1, 4),
                   np.array([w / W for *_, w in S]),
                   [np.array(vs, float) / D for _, vs in G], np.array([w / W for w, _ in G]), tol)

    def seg_fracs(self, cx, cy, theta, idx=None):
        """(len(poses), m) matrix of length fractions of segments in the closed square (for LP rows)."""
        segs = self.segs if idx is None else self.segs[idx]
        cx = np.atleast_1d(np.asarray(cx, float))[:, None]
        cy = np.atleast_1d(np.asarray(cy, float))[:, None]
        th = np.atleast_1d(np.asarray(theta, float))[:, None]
        c, s = np.cos(th), np.sin(th)
        ax, ay, bx, by = segs[:, 0][None], segs[:, 1][None], segs[:, 2][None], segs[:, 3][None]
        u0 = c * (ax - cx) + s * (ay - cy)
        v0 = -s * (ax - cx) + c * (ay - cy)
        u1 = c * (bx - cx) + s * (by - cy)
        v1 = -s * (bx - cx) + c * (by - cy)
        h = 0.5 + self.tol
        t0 = np.zeros(np.broadcast(u0, cx).shape)
        t1 = np.ones_like(t0)
        for p0, p1 in ((u0, u1), (v0, v1)):
            d = p1 - p0
            small = np.abs(d) < 1e-15
            with np.errstate(divide='ignore', invalid='ignore'):
                ta = np.where(small, -np.inf, (-h - p0) / np.where(small, 1.0, d))
                tb = np.where(small, np.inf, (h - p0) / np.where(small, 1.0, d))
            lo = np.minimum(ta, tb)
            hi = np.maximum(ta, tb)
            inside = (p0 >= -h) & (p0 <= h)
            lo = np.where(small & ~inside, np.inf, lo)
            hi = np.where(small & ~inside, -np.inf, hi)
            t0 = np.maximum(t0, lo)
            t1 = np.minimum(t1, hi)
        return np.clip(t1 - t0, 0.0, 1.0)

    def point_in(self, cx, cy, theta, idx=None):
        pts = self.pts if idx is None else self.pts[idx]
        cx = np.atleast_1d(np.asarray(cx, float))[:, None]
        cy = np.atleast_1d(np.asarray(cy, float))[:, None]
        th = np.atleast_1d(np.asarray(theta, float))[:, None]
        c, s = np.cos(th), np.sin(th)
        x = pts[:, 0][None] - cx
        y = pts[:, 1][None] - cy
        h = 0.5 + self.tol
        return ((np.abs(c * x + s * y) <= h) & (np.abs(-s * x + c * y) <= h)).astype(float)

    def mass(self, cx, cy, theta, chunk=4096):
        cx = np.atleast_1d(np.asarray(cx, float))
        cy = np.atleast_1d(np.asarray(cy, float))
        th = np.broadcast_to(np.atleast_1d(np.asarray(theta, float)), cx.shape)
        out = np.zeros(len(cx))
        for a in range(0, len(cx), chunk):
            b = min(len(cx), a + chunk)
            if len(self.pw):
                out[a:b] += self.point_in(cx[a:b], cy[a:b], th[a:b]) @ self.pw
            if len(self.sw):
                out[a:b] += self.seg_fracs(cx[a:b], cy[a:b], th[a:b]) @ self.sw
            for P, w in zip(self.polys, self.gw):
                A = 0.5 * abs(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))
                for i in range(a, b):
                    c, s = math.cos(th[i]), math.sin(th[i])
                    out[i] += w * _clip_poly([tuple(v) for v in P], c, s, cx[i], cy[i], self.tol) / A
        return out


def write(path, cover, comment=None):
    with open(path, 'w') as f:
        f.write('mixed 1\n')
        if comment:
            for line in comment.splitlines():
                f.write('# %s\n' % line)
        f.write('%d %d\n%d\n%d\n' % (cover['s_num'], cover['s_den'], cover['D'], cover['W']))
        f.write('%d\n' % len(cover['points']))
        for p in cover['points']:
            f.write('%d %d %d\n' % tuple(p))
        f.write('%d\n' % len(cover['segments']))
        for q in cover['segments']:
            f.write('%d %d %d %d %d\n' % tuple(q))
        f.write('%d\n' % len(cover['polygons']))
        for w, vs in cover['polygons']:
            f.write('%d %d %s\n' % (len(vs), w, ' '.join('%d %d' % v for v in vs)))


if __name__ == '__main__':
    for p in sys.argv[1:]:
        cv = load(p)
        validate(cv)
        print('%s: s=%s D=%d W=%d points=%d segments=%d polygons=%d total=%s = %.9f' % (
            p, cv['s'], cv['D'], cv['W'], len(cv['points']), len(cv['segments']), len(cv['polygons']),
            total(cv), float(total(cv))))
