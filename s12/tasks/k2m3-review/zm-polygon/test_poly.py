#!/usr/bin/env python3
"""Independent exact tests of zm_mixed's polygon path (Lemma S(b)) as used by qx2_zm.QXChecker.

Run:  taskset -c 6,7 python3 test_poly.py [N]
Everything that could claim a violation is Fraction arithmetic; geometry (square vertices, clipping, area,
segment length in Q) is re-implemented here, not taken from zm_mixed.
"""
import sys, os, random, math, copy
from fractions import Fraction as F
S = os.path.expanduser('~/math/square-packing/public/s12/search')
sys.path.insert(0, S)
import zm_mixed as ZM, zeromargin as zm, mixed_cover as MC

HALF = F(1, 2)
rng = random.Random(int(os.environ.get('SEED', '1')))


# ------------------------------------------------------------------ independent exact geometry
def trig(u):
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d


def square(cx, cy, u):
    """CCW vertices of the closed unit square with centre (cx,cy), angle 2 atan u."""
    c, s = trig(u)
    out = []
    for a, b in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)):
        out.append((cx + a * c - b * s, cy + a * s + b * c))
    return out


def in_square(p, cx, cy, u):
    c, s = trig(u)
    dx, dy = p[0] - cx, p[1] - cy
    X = dx * c + dy * s; Y = -dx * s + dy * c
    return -HALF <= X <= HALF and -HALF <= Y <= HALF


def area(V):
    return sum(V[i][0] * V[(i + 1) % len(V)][1] - V[(i + 1) % len(V)][0] * V[i][1] for i in range(len(V))) / 2


def clip_rect(V, x0, x1, y0, y1):
    """polygon V ∩ [x0,x1]x[y0,y1], my own Sutherland-Hodgman on the 4 axis half-planes."""
    def clip(V, f):
        out = []
        for i in range(len(V)):
            P, Q = V[i], V[(i + 1) % len(V)]
            a, b = f(P), f(Q)
            if a >= 0: out.append(P)
            if (a > 0 > b) or (a < 0 < b):
                t = a / (a - b)
                out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
        return out
    for f in (lambda p: p[0] - x0, lambda p: x1 - p[0], lambda p: p[1] - y0, lambda p: y1 - p[1]):
        V = clip(V, f)
        if not V: return []
    return V


def area_QU(cx, cy, u, a, b):
    V = clip_rect(square(cx, cy, u), a, b, a, b)
    return area(V) if len(V) >= 3 else F(0)


def seg_len_in_Q(P0, P1, cx, cy, u):
    """length of the (closed) segment P0P1 inside the closed square, exact up to the sqrt of the direction norm;
    returns the parameter fraction in [0,1] of the segment inside Q."""
    c, s = trig(u)
    lo, hi = F(0), F(1)
    def coords(p):
        dx, dy = p[0] - cx, p[1] - cy
        return dx * c + dy * s, -dx * s + dy * c
    A = coords(P0); B = coords(P1)
    for k in (0, 1):
        a0, d = A[k], B[k] - A[k]
        for sg in (1, -1):          # sg*(a0 + t d) <= 1/2
            aa, dd = sg * a0, sg * d
            if dd == 0:
                if aa > HALF: return F(0)
            elif dd > 0: hi = min(hi, (HALF - aa) / dd)
            else: lo = max(lo, (HALF - aa) / dd)
    return max(F(0), hi - lo)


def exact_mass_indep(cv, cx, cy, u, a, b):
    D, W = cv['D'], cv['W']
    tot = F(0)
    for X0, Y0, X1, Y1, w in cv['segments']:
        tot += F(w, W) * seg_len_in_Q((F(X0, D), F(Y0, D)), (F(X1, D), F(Y1, D)), cx, cy, u)
    tot += area_QU(cx, cy, u, a, b)
    return tot


def admissible(m, cx, cy, u):
    c, s = trig(u)
    h = (c + s) / 2
    return h <= cx <= m - h and h <= cy <= m - h


def adm_rect(box, m, u):
    cx0, cx1, cy0, cy1 = box[:4]
    c, s = trig(u); h = (c + s) / 2
    x0, x1, y0, y1 = max(cx0, h), min(cx1, m - h), max(cy0, h), min(cy1, m - h)
    if x0 > x1 or y0 > y1: return None
    return x0, x1, y0, y1


def rand_frac(lo, hi, den=10 ** 6):
    if lo == hi: return lo
    return lo + (hi - lo) * F(rng.randrange(den + 1), den)


def sample_poses(box, m, nu=12, nc=4):
    cx0, cx1, cy0, cy1, u0, u1 = box
    us = [u0, u1] + [rand_frac(u0, u1) for _ in range(nu)]
    out = []
    for u in us:
        r = adm_rect(box, m, u)
        if r is None: continue
        x0, x1, y0, y1 = r
        for cx in (x0, x1):
            for cy in (y0, y1): out.append((cx, cy, u))
        for _ in range(nc): out.append((rand_frac(x0, x1), rand_frac(y0, y1), u))
    return out


# ------------------------------------------------------------------ the polygon region, as piece_bound builds it
def region_K(chk, box, B):
    """replicates zm_mixed.piece_bound l.1109-1117 by calling its own functions (for the vertex test)."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    specs = chk.zc._adm_specs(box, B)
    r = ZM.REACH
    frame = [(cx0 - r, cy0 - r), (cx1 + r, cy0 - r), (cx1 + r, cy1 + r), (cx0 - r, cy1 + r)]
    K = frame
    for c in range(4):
        Kc = ZM.cond_region(specs, c, u0, u1 - u0, chk.m, (cx0, cy0), frame)
        K = ZM.clip_convex(K, Kc) if Kc else []
        if not K: break
    return K, specs


def rand_box(m, a, b, scale):
    """a box whose centre range straddles / hugs a side or corner of U or a container wall."""
    kind = rng.random()
    def coord():
        t = rng.random()
        if t < 0.35: base = a + rng.choice([-1, 1]) * F(1, 2) + rand_frac(F(-1, 5), F(1, 5), 1000)
        elif t < 0.55: base = a + rand_frac(F(-1, 5), F(1, 5), 1000)
        elif t < 0.75: base = rand_frac(F(1, 2), F(3, 4), 1000)       # near the wall
        else: base = rand_frac(F(1, 2), m / 2, 1000)
        w = scale * F(rng.randrange(1, 101), 100)
        if rng.random() < 0.1: w = F(0)                                # zero-width centre range
        return base, base + w
    x0, x1 = coord(); y0, y1 = coord()
    t = rng.random()
    if t < 0.2: u0 = F(0)
    else: u0 = rand_frac(F(0), F(9, 20), 1000)
    uw = scale * F(rng.randrange(0, 101), 100)
    u1 = min(u0 + uw, F(1, 2))
    return (x0, x1, y0, y1, u0, u1)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    path = os.path.join(S, 'qx2_data/L4_k02_box7.txt')
    cv = MC.load(path)
    cov = ZM.Cover(cv)
    import qx2_zm as QX
    a, b = QX.u_square(cov)
    m = cov.m
    chk = ZM.MixedChecker(cov, max_depth=18, use_chain=False, cert_mode=False)
    cv_nopoly = dict(cv, polygons=[])
    chk0 = ZM.MixedChecker(ZM.Cover(cv_nopoly), max_depth=18, use_chain=False, cert_mode=True)
    stats = dict(boxes=0, poses=0, pos_poly=0, viol_poly=0, viol_vertex=0, viol_total=0, minratio=None,
                 tightest=None)
    for it in range(N):
        scale = rng.choice([F(1, 4), F(1, 20), F(1, 200), F(1, 5000)])
        box = rand_box(m, a, b, scale)
        if rng.random() < 0.5:   # mirror to the far side of U (b) sometimes
            x0, x1, y0, y1, u0, u1 = box
            if rng.random() < 0.5: x0, x1 = m - x1, m - x0
            if rng.random() < 0.5: y0, y1 = m - y1, m - y0
            box = (x0, x1, y0, y1, u0, u1)
        cu1 = zm.clip_bin(box, m) if box[5] > box[4] else box[5]
        box = box[:5] + (cu1,)
        poses = [p for p in sample_poses(box, m) if admissible(m, *p)]
        if not poses: continue
        B = zm.bin_data(box[4], box[5])
        L, _ = chk.piece_bound(box, B)
        L0, _ = chk0.piece_bound(box, B)
        polyL = L - L0
        stats['boxes'] += 1
        if polyL > 0: stats['pos_poly'] += 1
        K, specs = region_K(chk, box, B)
        # polygon contribution recomputed independently from K
        if K:
            Ki = clip_rect(K, a, b, a, b)
            ak = area(Ki) if len(Ki) >= 3 else F(0)
            assert ak == polyL, ('polygon part mismatch', box, ak, polyL)
        else:
            assert polyL == 0
        for (cx, cy, u) in poses:
            stats['poses'] += 1
            A = area_QU(cx, cy, u, a, b)
            if polyL > A:
                stats['viol_poly'] += 1
                print('POLY VIOLATION', box, (cx, cy, u), polyL, A, flush=True)
            if K:
                for v in K:
                    if not in_square(v, cx, cy, u):
                        stats['viol_vertex'] += 1
                        print('VERTEX OUTSIDE Q', box, (cx, cy, u), v, flush=True)
                        break
            if polyL > 0:
                gap = A - polyL
                if stats['tightest'] is None or gap < stats['tightest'][0]:
                    stats['tightest'] = (gap, box, (cx, cy, u))
            if rng.random() < 0.15:
                M = exact_mass_indep(cv, cx, cy, u, a, b)
                if L > M:
                    stats['viol_total'] += 1
                    print('TOTAL VIOLATION', box, (cx, cy, u), L, M, flush=True)
    t = stats.pop('tightest')
    print(stats)
    if t: print('smallest area(Q∩U) - polygon bound over samples: %.3e at box %s pose %s' % (
        float(t[0]), [str(v) for v in t[1]], [str(v) for v in t[2]]))


if __name__ == '__main__':
    main()
