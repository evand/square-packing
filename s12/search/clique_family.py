#!/usr/bin/env python3
"""Anchor cliques: certifiable pose families that pairwise closed-intersect  (task G, part 1).

An *anchor family* is  Q(A) = { S admissible closed unit square in [0,t]^2 : A subset of S }
for a compact convex A.  Q({p}) is the point clique of p, which is exactly the family the point
constraint  coverage(p) <= 1  imposes.  An *anchor clique* is a union

    K(A_1, ..., A_m) = Q(A_1) u ... u Q(A_m)

which is a clique (pairwise closed-intersecting) as soon as every S in Q(A_i) meets every
S' in Q(A_j) for all i, j.  Two sufficient conditions used here:

  (same j)   two squares containing A_j both contain A_j, so they meet   (A_j nonempty);
  (cross)    if every admissible S with A_i subset S meets A_j, then S meets every S' in Q(A_j)
             (because S n S' contains S n A_j).

The interesting case is m = 2 with A_1 = {p} a single point and A_2 = A a segment slightly
"inside" p:  K(p, A) = {S : p in S} u {S : A subset S}.  Its extra content over the point
constraint at p is  mu(Q(A) \ P_p) -- the squares that contain A but *not* p, i.e. the squares
whose boundary passes between p and A.  Without a container wall this family is empty
(Lemma 1 below); near a wall it is not.

Routines
    rho_star(p, d, eps, t)   the least half-length rho for which the segment
                             A = {p + eps*d + s*d_perp : |s| <= rho} is a transversal of
                             {admissible S : p in S};  +inf if no rho works.
    clique_mass(...)         mass of K(p, A) under a finitely supported measure.

Modes
    python3 search/clique_family.py rho      -- rho*(eps) tables and the analytic bound
    python3 search/clique_family.py measure  -- anchor cliques on the certified 3.99 measure
    python3 search/clique_family.py stress   -- randomised stress test of the transversal claim
"""
import sys, os, math, argparse, json, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
MAIN_RUNS = '/home/evand/math/square-packing/s12/runs'

INF = float('inf')


# ------------------------------------------------------------------ elementary pose geometry
def wid(th):
    return abs(math.cos(th)) + abs(math.sin(th))


def frame(th):
    ct, st = math.cos(th), math.sin(th)
    return (ct, st), (-st, ct)


def in_square(pt, c, th):
    """closed unit square at (c, th) contains pt?"""
    dx, dy = pt[0] - c[0], pt[1] - c[1]
    ct, st = math.cos(th), math.sin(th)
    return abs(dx * ct + dy * st) <= 0.5 and abs(-dx * st + dy * ct) <= 0.5


def admissible(c, th, t, tol=0.0):
    w2 = wid(th) / 2
    return (w2 - tol <= c[0] <= t - w2 + tol) and (w2 - tol <= c[1] <= t - w2 + tol)


# ------------------------------------------------------------------ convex polygon utilities
def clip(poly, a, b, cc):
    """keep the part of the polygon with a*x + b*y <= cc (Sutherland-Hodgman)"""
    if not poly:
        return poly
    out = []
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        d1 = a * x1 + b * y1 - cc
        d2 = a * x2 + b * y2 - cc
        if d1 <= 0:
            out.append((x1, y1))
        if (d1 < 0 < d2) or (d2 < 0 < d1):
            u = d1 / (d1 - d2)
            out.append((x1 + u * (x2 - x1), y1 + u * (y2 - y1)))
    return out


def centre_polygon(p, th, t):
    """the set of centres c with S(c, th) admissible and p in S(c, th)  (a convex polygon)"""
    w2 = wid(th) / 2
    poly = [(w2, w2), (t - w2, w2), (t - w2, t - w2), (w2, t - w2)]
    if t - w2 < w2:
        return []
    (e1x, e1y), (e2x, e2y) = frame(th)
    for (ex, ey) in ((e1x, e1y), (e2x, e2y)):
        pe = p[0] * ex + p[1] * ey
        poly = clip(poly, -ex, -ey, -(pe - 0.5))     # <c,e> >= <p,e> - 1/2
        poly = clip(poly, ex, ey, pe + 0.5)          # <c,e> <= <p,e> + 1/2
    return poly


# ------------------------------------------------------------------ the transversal threshold
def rho_star_theta(p, d, eps, t, th, want_witness=False):
    """sup over admissible S at angle th with p in S of the least rho making A a transversal.

    A = { p + eps*d + s*d_perp : |s| <= rho }, a segment perpendicular to d at offset eps.
    A square S and the segment A are disjoint iff one of THREE separating axes works
    (SAT: edge normals of both polygons; the segment's only edge normal is d):

      (0)  the axis d:  max_{x in S} <x, d>  <  <m, d>   (S does not reach the line of A)
           -- no rho can repair this, so it caps eps;
      (i, sigma), i = 1, 2:  sigma*<a - c, e_i> > 1/2 for BOTH endpoints a of A, i.e.
           sigma*<m - c, e_i> - rho*|<e_i, d_perp>| > 1/2.

    So S misses A iff (0) holds, or rho < (sigma*<m - c, e_i> - 1/2)/|<e_i, d_perp>| for some
    (i, sigma).  Returns the max of the latter over admissible centres (+inf if (0) is reachable,
    or if a separating edge is parallel to A)."""
    poly = centre_polygon(p, th, t)
    if not poly:
        return (0.0, None) if want_witness else 0.0
    dp = (-d[1], d[0])
    m = (p[0] + eps * d[0], p[1] + eps * d[1])
    es = frame(th)
    md = m[0] * d[0] + m[1] * d[1]
    # --- axis (0): the least reach of S in direction d over admissible centres
    wd = abs(es[0][0] * d[0] + es[0][1] * d[1]) + abs(es[1][0] * d[0] + es[1][1] * d[1])
    reach = min(cx * d[0] + cy * d[1] for (cx, cy) in poly) + wd / 2
    if reach < md - 1e-15:
        return (INF, ('reach', th, reach, md)) if want_witness else INF
    best = 0.0
    bw = None
    for e in es:
        den = abs(e[0] * dp[0] + e[1] * dp[1])
        for sg in (1.0, -1.0):
            # objective sigma*<m - c, e> - 1/2, linear in c: maximise over the polygon
            val = max(sg * ((m[0] - cx) * e[0] + (m[1] - cy) * e[1]) - 0.5 for (cx, cy) in poly)
            if val <= 0:
                continue
            if den < 1e-12:                      # edge parallel to A: no rho can help
                return (INF, ('parallel', th)) if want_witness else INF
            r = val / den
            if r > best:
                best = r
                if want_witness:
                    cbest = max(poly, key=lambda c: sg * ((m[0] - c[0]) * e[0] + (m[1] - c[1]) * e[1]))
                    bw = (th, cbest, e, sg)
    return (best, bw) if want_witness else best


def eps_max(p, d, t, nth=4096):
    """the largest offset eps for which the axis-d test (0) never separates:
    eps_max = min over theta of [ min_c <c,d> + w_d(theta)/2 ] - <p,d>."""
    pd = p[0] * d[0] + p[1] * d[1]
    best = INF
    bth = 0.0
    for k in range(nth + 1):
        th = (math.pi / 2) * k / nth
        poly = centre_polygon(p, th, t)
        if not poly:
            continue
        es = frame(th)
        wd = abs(es[0][0] * d[0] + es[0][1] * d[1]) + abs(es[1][0] * d[0] + es[1][1] * d[1])
        reach = min(cx * d[0] + cy * d[1] for (cx, cy) in poly) + wd / 2
        if reach - pd < best:
            best, bth = reach - pd, th
    return best, bth


def rho_star(p, d, eps, t, nth=2048, refine=3):
    """rho*(p, d, eps) = sup over all angles.  Grid + local refinement (heuristic but tight)."""
    best = 0.0
    bth = 0.0
    step = (math.pi / 2) / nth
    for k in range(nth):
        th = k * step
        r = rho_star_theta(p, d, eps, t, th)
        if r == INF:
            return INF, th
        if r > best:
            best, bth = r, th
    for _ in range(refine):
        lo, hi = bth - step, bth + step
        step = (hi - lo) / 64
        for k in range(65):
            th = lo + k * step
            r = rho_star_theta(p, d, eps, t, th)
            if r == INF:
                return INF, th
            if r > best:
                best, bth = r, th
    return best, bth


# ---------------------------------------------------------------- general convex anchors
#
# The complete test.  A (compact convex) is a transversal of U(p) = {admissible S : p in S}
# iff  for every angle th:   centre_polygon(p, th, t)  subset  A (+) R_th Q,
# because S(c, th) meets A iff c lies in the Minkowski sum A (+) R_th Q (Q the unit square,
# which is centrally symmetric).  Both sets are convex, so it is enough to test the VERTICES
# of centre_polygon, and "v in A (+) R_th Q" is just "the unit square at (v, th) meets A".

def poly_sq(c, th):
    ct, st = math.cos(th), math.sin(th)
    out = []
    for (sx, sy) in ((0.5, 0.5), (-0.5, 0.5), (-0.5, -0.5), (0.5, -0.5)):
        out.append((c[0] + sx * ct - sy * st, c[1] + sx * st + sy * ct))
    return out


def polys_meet(P, R, tol=0.0):
    """do two convex polygons (CCW or CW, any order) have a common point?  SAT, closed."""
    for poly in (P, R):
        n = len(poly)
        for i in range(n):
            x1, y1 = poly[i]
            x2, y2 = poly[(i + 1) % n]
            ax, ay = -(y2 - y1), (x2 - x1)
            nn = math.hypot(ax, ay)
            if nn < 1e-14:
                continue
            ax /= nn
            ay /= nn
            p0 = min(ax * x + ay * y for (x, y) in P)
            p1 = max(ax * x + ay * y for (x, y) in P)
            r0 = min(ax * x + ay * y for (x, y) in R)
            r1 = max(ax * x + ay * y for (x, y) in R)
            if p1 < r0 - tol or r1 < p0 - tol:
                return False
    return True


def is_transversal(p, A, t, nth=2048, refine=2, tol=0.0):
    """A (list of vertices of a compact convex set, or 1-2 points) a transversal of U(p)?
    Returns (ok, worst_witness).  Grid over angles with local refinement."""
    def bad_at(th):
        poly = centre_polygon(p, th, t)
        if not poly:
            return None
        for v in poly:
            if not polys_meet(poly_sq(v, th), A, tol=tol):
                return (th, v)
        return None
    step = (math.pi / 2) / nth
    hits = []
    for k in range(nth + 1):
        b = bad_at(k * step)
        if b is not None:
            return False, b
        hits.append(k * step)
    # refinement: re-test between grid points (the vertices move continuously except at
    # combinatorial changes of the polygon, so a second finer pass is a cheap guard)
    for r in range(refine):
        step /= 8
        for k in range(8 * nth + 1):
            b = bad_at(k * step)
            if b is not None:
                return False, b
        break
    return True, None


def _pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def transversal_margin(p, A, t, thetas, box=None):
    """vectorised: for each angle, min over the vertices of centre_polygon(p, th, t) of the
    SAT margin between the unit square at that centre and the convex set A.

    margin >= 0 at every angle  <=>  A is a transversal of U(p) = {admissible S : p in S}
    (the centres where S meets A form the convex set A (+) R_th Q, so it is enough to check
    the vertices of the convex centre polygon).  Returns (margin[T], witness centres[T,2]).

    A: (m, 2) array of vertices of a compact convex set (a point m=1, a segment m=2, a
    polygon m>=3, given in convex position, any order)."""
    import numpy as np
    th = np.asarray(thetas, dtype=float)
    T = th.size
    ct = np.cos(th)
    st = np.sin(th)
    w = np.abs(ct) + np.abs(st)
    px, py = p
    nc = 8 if box is None else 12
    # constraints  a*cx + b*cy <= d
    a = np.zeros((T, nc))
    b = np.zeros((T, nc))
    d = np.zeros((T, nc))
    if box is not None:
        bx0, bx1, by0, by1 = box
        a[:, 8], b[:, 8], d[:, 8] = -1.0, 0.0, -bx0
        a[:, 9], b[:, 9], d[:, 9] = 1.0, 0.0, bx1
        a[:, 10], b[:, 10], d[:, 10] = 0.0, -1.0, -by0
        a[:, 11], b[:, 11], d[:, 11] = 0.0, 1.0, by1
    a[:, 0], b[:, 0], d[:, 0] = -1.0, 0.0, -w / 2                     # cx >= w/2
    a[:, 1], b[:, 1], d[:, 1] = 1.0, 0.0, t - w / 2                   # cx <= t - w/2
    a[:, 2], b[:, 2], d[:, 2] = 0.0, -1.0, -w / 2
    a[:, 3], b[:, 3], d[:, 3] = 0.0, 1.0, t - w / 2
    pe1 = px * ct + py * st
    pe2 = -px * st + py * ct
    a[:, 4], b[:, 4], d[:, 4] = -ct, -st, -(pe1 - 0.5)                # <c,e1> >= <p,e1> - 1/2
    a[:, 5], b[:, 5], d[:, 5] = ct, st, pe1 + 0.5
    a[:, 6], b[:, 6], d[:, 6] = st, -ct, -(pe2 - 0.5)
    a[:, 7], b[:, 7], d[:, 7] = -st, ct, pe2 + 0.5
    P = _pairs(nc)
    ii = np.array([q[0] for q in P])
    jj = np.array([q[1] for q in P])
    det = a[:, ii] * b[:, jj] - a[:, jj] * b[:, ii]                   # (T, 28)
    ok = np.abs(det) > 1e-12
    dets = np.where(ok, det, 1.0)
    vx = (d[:, ii] * b[:, jj] - d[:, jj] * b[:, ii]) / dets
    vy = (a[:, ii] * d[:, jj] - a[:, jj] * d[:, ii]) / dets
    # feasibility of each candidate vertex
    feas = ok.copy()
    tolc = 1e-9
    for k in range(nc):
        feas &= (a[:, k][:, None] * vx + b[:, k][:, None] * vy <= d[:, k][:, None] + tolc)
    Av = np.asarray(A, dtype=float).reshape(-1, 2)
    # SAT axes: A's edge normals (fixed) + the square's two normals (angle dependent)
    fixed_axes = []
    m = Av.shape[0]
    if m >= 2:
        for i in range(m):
            j = (i + 1) % m
            ex, ey = -(Av[j, 1] - Av[i, 1]), (Av[j, 0] - Av[i, 0])
            nn = math.hypot(ex, ey)
            if nn > 1e-14:
                fixed_axes.append((ex / nn, ey / nn))
            if m == 2:
                break
    gapmax = np.full((T, len(P)), -np.inf)
    for (ex, ey) in fixed_axes:
        hA1 = float(np.max(Av[:, 0] * ex + Av[:, 1] * ey))
        hA0 = float(np.min(Av[:, 0] * ex + Av[:, 1] * ey))
        r = 0.5 * (np.abs(ct * ex + st * ey) + np.abs(-st * ex + ct * ey))     # (T,)
        cen = vx * ex + vy * ey
        s0 = cen - r[:, None]
        s1 = cen + r[:, None]
        gap = np.maximum(s0 - hA1, hA0 - s1)
        gapmax = np.maximum(gapmax, gap)
    for (ex, ey) in ((ct, st), (-st, ct)):
        hA1 = np.max(Av[:, 0] * ex[:, None] + Av[:, 1] * ey[:, None], axis=1)   # (T,)
        hA0 = np.min(Av[:, 0] * ex[:, None] + Av[:, 1] * ey[:, None], axis=1)
        cen = vx * ex[:, None] + vy * ey[:, None]
        s0 = cen - 0.5
        s1 = cen + 0.5
        gap = np.maximum(s0 - hA1[:, None], hA0[:, None] - s1)
        gapmax = np.maximum(gapmax, gap)
    marg = np.where(feas, -gapmax, np.inf)
    idx = np.argmin(marg, axis=1)
    r = np.arange(T)
    out = marg[r, idx]
    out = np.where(np.isfinite(out), out, np.inf)     # angle with no admissible pose through p
    wit = np.stack([vx[r, idx], vy[r, idx]], axis=1)
    return out, wit


def is_transversal_np(p, A, t, nth=4096, refine=True, box=None, thlo=0.0, thhi=math.pi / 2):
    """(ok, worst margin, worst angle, witness centre).  With box = (x0, x1, y0, y1) and an
    angle range, the claim is restricted to the poses of that box -- the family

        K(p, B, A) = {S : p in S, pose in B}  u  {S : A subset S}

    is a clique as soon as every admissible S in B with p in S meets A."""
    import numpy as np
    th = np.linspace(thlo, thhi, nth + 1)
    marg, wit = transversal_margin(p, A, t, th, box=box)
    k = int(np.argmin(marg))
    if refine:
        step = (thhi - thlo) / nth
        lo = max(thlo, th[k] - step)
        hi = min(thhi, th[k] + step)
        th2 = np.linspace(lo, hi, 2049)
        m2, w2 = transversal_margin(p, A, t, th2, box=box)
        k2 = int(np.argmin(m2))
        if m2[k2] < marg[k]:
            return bool(m2[k2] >= 0), float(m2[k2]), float(th2[k2]), tuple(w2[k2])
    return bool(marg[k] >= 0), float(marg[k]), float(th[k]), tuple(wit[k])


# ---------------------------------------------------------------- the maximal set K(p)
#
# K(p) = { S' admissible : S' meets every admissible S with p in S }.  For a FIXED angle th'
# of S', K(p) is convex in the centre c':  S' meets S iff c' lies in the convex set
# S (+) R_th' Q, so K_th'(p) = intersection over S in U(p) of (S (+) R_th' Q), intersected
# with admissibility.  U(p) is approximated from OUTSIDE by its extreme members (the vertices
# of the centre polygons over an angle grid), which makes K_th'(p) an OUTER approximation --
# the honest direction for "how big can the extra family be".

def extreme_squares(p, t, nth=512):
    """the vertex poses of centre_polygon(p, th, t) over an angle grid: (cx, cy, th) array"""
    import numpy as np
    out = []
    for k in range(nth + 1):
        th = (math.pi / 2) * k / nth
        poly = centre_polygon(p, th, t)
        for (cx, cy) in poly:
            out.append((cx, cy, th))
    return np.array(out) if out else np.zeros((0, 3))


def k_normals(E, nth):
    """the candidate outward normals: +-e1, +-e2 of every angle occurring in E"""
    import numpy as np
    ths = np.unique(np.round(E[:, 2], 12))
    ct, st = np.cos(ths), np.sin(ths)
    nx = np.concatenate([ct, -ct, -st, st])
    ny = np.concatenate([st, -st, ct, -ct])
    return nx, ny


def k_region(p, t, thp, E, NN=None):
    """K_th'(p) as a convex polygon of centres c' (outer approximation from E)"""
    import numpy as np
    ctp, stp = math.cos(thp), math.sin(thp)
    cS = E[:, :2]
    thS = E[:, 2]
    ctS, stS = np.cos(thS), np.sin(thS)
    nx0, ny0 = NN if NN is not None else k_normals(E, 0)
    nx = np.concatenate([nx0, [ctp, -ctp, -stp, stp]])
    ny = np.concatenate([ny0, [stp, -stp, ctp, -ctp]])
    # for each normal n:  H(n) = min over S of [ <c_S, n> + w_n(th_S)/2 ] + w_n(th')/2
    d1 = np.abs(np.outer(nx, ctS) + np.outer(ny, stS))
    d2 = np.abs(np.outer(-nx, stS) + np.outer(ny, ctS))
    wS = 0.5 * (d1 + d2)                                   # (|N|, |E|)
    wP = 0.5 * (np.abs(nx * ctp + ny * stp) + np.abs(-nx * stp + ny * ctp))    # (|N|,)
    proj = np.outer(nx, cS[:, 0]) + np.outer(ny, cS[:, 1])
    H = np.min(proj + wS, axis=1) + wP
    # admissibility of S'
    w2 = (abs(ctp) + abs(stp)) / 2
    poly = [(w2, w2), (t - w2, w2), (t - w2, t - w2), (w2, t - w2)]
    # clip by the tightest constraints first
    order = np.argsort(H - np.array([max(nx[i] * q[0] + ny[i] * q[1] for q in poly)
                                     for i in range(nx.size)]))
    for i in order:
        poly = clip(poly, float(nx[i]), float(ny[i]), float(H[i]))
        if not poly:
            return []
    return poly


def poly_area(poly):
    n = len(poly)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2


def p_region(p, t, thp):
    """{c' admissible : p in S'(c', th')} -- the point clique's slice"""
    return centre_polygon(p, thp, t)


# ---------------------------------------------------------------- graph cliques (float)
def closed_adj(imgs, tol=1e-9):
    """adjacency of the closed-intersection graph of poses (n, 3).  Touching counts."""
    import numpy as np
    cx, cy = imgs[:, 0], imgs[:, 1]
    ct, st = np.cos(imgs[:, 2]), np.sin(imgs[:, 2])
    n = len(imgs)
    adj = np.ones((n, n), dtype=bool)
    for (ex, ey) in ((ct, st), (-st, ct)):
        proj = np.outer(ex, cx) + np.outer(ey, cy)
        r = 0.5 * (np.abs(np.outer(ex, ct) + np.outer(ey, st))
                   + np.abs(np.outer(-ex, st) + np.outer(ey, ct)))
        lo, hi = proj - r, proj + r
        own = np.diagonal(lo).copy()
        ownh = np.diagonal(hi).copy()
        adj &= ~((lo > ownh[:, None] + tol) | (hi < own[:, None] - tol))
    np.fill_diagonal(adj, True)
    return adj & adj.T


def greedy_clique(adj, w, order):
    import numpy as np
    best = []
    bw = 0.0
    for v in order:
        cl = [v]
        cand = np.nonzero(adj[v])[0]
        cand = cand[cand != v]
        cand = cand[np.argsort(-w[cand])]
        for u in cand:
            if all(adj[u, x] for x in cl):
                cl.append(u)
        s = float(sum(w[x] for x in cl))
        if s > bw:
            bw, best = s, cl
    return bw, best


def max_weight_clique(adj, w, tlimit=30.0):
    """branch and bound with a greedy-colouring weight bound"""
    import numpy as np
    n = len(w)
    order = np.argsort(-w)
    bw, best = greedy_clique(adj, w, order[:min(n, 80)])
    t0 = time.time()

    def colour_bound(cand):
        rem = list(cand[np.argsort(-w[cand])])
        tot = 0.0
        cls = []
        while rem:
            cl = []
            rest = []
            for v in rem:
                if all(not adj[v, u] for u in cl):
                    cl.append(v)
                else:
                    rest.append(v)
            tot += max(w[v] for v in cl)
            cls.append(cl)
            rem = rest
        return tot, cls

    def expand(cur, curw, cand):
        nonlocal best, bw
        if time.time() - t0 > tlimit:
            return
        if len(cand) == 0:
            if curw > bw:
                bw, best = curw, list(cur)
            return
        ub, cls = colour_bound(cand)
        if curw + ub <= bw + 1e-12:
            return
        seq = [v for cl in cls for v in cl][::-1]
        for i, v in enumerate(seq):
            rest = np.array(seq[i + 1:], dtype=int)
            if curw + w[v] + (colour_bound(rest)[0] if len(rest) else 0.0) <= bw + 1e-12:
                continue
            nxt = rest[adj[v, rest]] if len(rest) else rest
            cur.append(v)
            expand(cur, curw + w[v], nxt)
            cur.pop()
            if time.time() - t0 > tlimit:
                return

    expand([], 0.0, np.arange(n))
    return bw, best


def poly_intersect(polys):
    """intersection of convex polygons (each a vertex list), by successive clipping"""
    if not polys:
        return []
    cur = polys[0]
    for R in polys[1:]:
        n = len(R)
        # orient R
        area = sum(R[i][0] * R[(i + 1) % n][1] - R[(i + 1) % n][0] * R[i][1] for i in range(n))
        RR = R if area > 0 else R[::-1]
        for i in range(n):
            x1, y1 = RR[i]
            x2, y2 = RR[(i + 1) % n]
            a, b = (y2 - y1), -(x2 - x1)          # inward: a*x + b*y <= a*x1 + b*y1
            cur = clip(cur, a, b, a * x1 + b * y1)
            if not cur:
                return []
    return cur


def seg_meets_square(a0, a1, c, th):
    """does the closed segment [a0, a1] meet the closed unit square at (c, th)?
    Independent of the SAT formulation above: Liang-Barsky clipping in the square's frame."""
    ct, st = math.cos(th), math.sin(th)
    u0 = (a0[0] - c[0]) * ct + (a0[1] - c[1]) * st
    v0 = -(a0[0] - c[0]) * st + (a0[1] - c[1]) * ct
    u1 = (a1[0] - c[0]) * ct + (a1[1] - c[1]) * st
    v1 = -(a1[0] - c[0]) * st + (a1[1] - c[1]) * ct
    lo, hi = 0.0, 1.0
    for (q0, q1) in ((u0, u1), (v0, v1)):
        dq = q1 - q0
        for (sgn, bnd) in ((1.0, 0.5), (-1.0, 0.5)):
            # sgn*q <= bnd
            num = bnd - sgn * q0
            den = sgn * dq
            if abs(den) < 1e-300:
                if num < 0:
                    return False
                continue
            r = num / den
            if den > 0:
                hi = min(hi, r)
            else:
                lo = max(lo, r)
            if lo > hi:
                return False
    return True


def brute_transversal(p, d, eps, rho, t, n=200000, seed=1):
    """random admissible poses through p: the worst (most separated) one found.
    Returns (nfail, worst_witness).  Independent numeric check of the transversal claim."""
    import random
    rng = random.Random(seed)
    dp = (-d[1], d[0])
    m = (p[0] + eps * d[0], p[1] + eps * d[1])
    a0 = (m[0] - rho * dp[0], m[1] - rho * dp[1])
    a1 = (m[0] + rho * dp[0], m[1] + rho * dp[1])
    nfail = 0
    worst = None
    for _ in range(n):
        th = rng.uniform(0, math.pi / 2)
        poly = centre_polygon(p, th, t)
        if not poly:
            continue
        # random point of the polygon by rejection in its bounding box
        xs = [q[0] for q in poly]
        ys = [q[1] for q in poly]
        for _try in range(40):
            cx = rng.uniform(min(xs), max(xs))
            cy = rng.uniform(min(ys), max(ys))
            if in_square(p, (cx, cy), th) and admissible((cx, cy), th, t, tol=-1e-12):
                break
        else:
            continue
        if not seg_meets_square(a0, a1, (cx, cy), th):
            nfail += 1
            worst = (th, cx, cy)
    return nfail, worst


def sat_meet(c1, th1, c2, th2, tol=0.0):
    """closed unit squares meet?  (separating-axis test with <=, touching counts)"""
    for (c, th) in ((c1, th1), (c2, th2)):
        ct, st = math.cos(th), math.sin(th)
        for (ex, ey) in ((ct, st), (-st, ct)):
            h1 = c1[0] * ex + c1[1] * ey
            h2 = c2[0] * ex + c2[1] * ey
            ct1, st1 = math.cos(th1), math.sin(th1)
            ct2, st2 = math.cos(th2), math.sin(th2)
            r1 = 0.5 * (abs(ct1 * ex + st1 * ey) + abs(-st1 * ex + ct1 * ey))
            r2 = 0.5 * (abs(ct2 * ex + st2 * ey) + abs(-st2 * ex + ct2 * ey))
            if h1 + r1 < h2 - r2 - tol or h2 + r2 < h1 - r1 - tol:
                return False
    return True


def transversal_violation(p, d, eps, rho, t, nth=4096):
    """max over admissible S ∋ p of (separation slack): > 0 means A is NOT a transversal."""
    worst = -INF
    ww = None
    for k in range(nth + 1):
        th = (math.pi / 2) * k / nth
        r, w = rho_star_theta(p, d, eps, t, th, want_witness=True)
        if r == INF:
            return INF, w
        if r - rho > worst:
            worst, ww = r - rho, w
    return worst, ww


# ------------------------------------------------------------------ measures
def read_exact_support(path):
    """(p, q, cx, cy, mass) with exact Fractions"""
    out = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        w = line.split()
        if w[0] == 'pose':
            out.append((int(w[1]), int(w[2]), Fr(w[3]), Fr(w[4]), Fr(w[5])))
    return out


def read_float_support(path):
    out = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        w = line.split()
        if w[0] == 'pose':
            out.append((float(w[1]), float(w[2]), math.radians(float(w[3])), float(w[4])))
    return out


def d4_images(cx, cy, th, t):
    return [(cx, cy, th), (t - cx, cy, -th), (cx, t - cy, -th), (t - cx, t - cy, th),
            (cy, cx, -th), (t - cy, cx, th), (cy, t - cx, th), (t - cy, t - cx, -th)]


def expand(poses, t):
    """poses = [(cx, cy, th, mu)] -> list of (cx, cy, th, mu/8) images"""
    imgs = []
    for (cx, cy, th, mu) in poses:
        for (x, y, a) in d4_images(cx, cy, th, t):
            imgs.append((x, y, a, mu / 8.0))
    return imgs


def load_measure(path, t):
    """returns list of images (cx, cy, th, mass) from an exact or float support file"""
    txt = open(path).readline()
    if 'p q cx cy' in open(path).read(2000):
        poses = []
        for (p, q, cx, cy, m) in read_exact_support(path):
            th = 2 * math.atan2(p, q)
            poses.append((float(cx), float(cy), th, float(m)))
    else:
        poses = read_float_support(path)
    return expand(poses, t), poses


def contains_segment(a0, a1, c, th):
    return in_square(a0, c, th) and in_square(a1, c, th)


def clique_mass(imgs, p, d, eps, rho):
    """mass of K(p, A) = {S ∋ p} u {S ⊇ A}; returns (total, cov_p, extra, n_extra)"""
    a0 = (p[0] + eps * d[0] - rho * (-d[1]), p[1] + eps * d[1] - rho * d[0])
    a1 = (p[0] + eps * d[0] + rho * (-d[1]), p[1] + eps * d[1] + rho * d[0])
    cov = 0.0
    extra = 0.0
    nex = 0
    for (cx, cy, th, m) in imgs:
        if in_square(p, (cx, cy), th):
            cov += m
        elif contains_segment(a0, a1, (cx, cy), th):
            extra += m
            nex += 1
    return cov + extra, cov, extra, nex


# ------------------------------------------------------------------ modes
def cmd_rho(args):
    t = args.t
    print(f"# rho*(p, d=+x, eps) for p = (px, py), container t = {t}")
    print(f"# analytic bound (Lemma 2):  rho* <= eps / (1 - px)     [needs px + eps < 1]")
    print()
    print(f"{'px':>7} {'py':>6} {'eps':>9} {'rho*':>10} {'eps/(1-px)':>12} {'theta*(deg)':>12} "
          f"{'ratio rho*/eps':>15}")
    for px in (0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5):
        for py in (1.275, 2.0):
            for frac in (0.1, 0.25, 0.5, 0.75, 0.9):
                eps = frac * (1 - px)
                r, th = rho_star((px, py), (1.0, 0.0), eps, t)
                bound = eps / (1 - px)
                print(f"{px:7.3f} {py:6.3f} {eps:9.5f} {r:10.5f} {bound:12.5f} "
                      f"{math.degrees(th):12.4f} {r/eps if eps else 0:15.4f}")
        print()


DIRS = [('+x', (1.0, 0.0)), ('-x', (-1.0, 0.0)), ('+y', (0.0, 1.0)), ('-y', (0.0, -1.0)),
        ('++', (0.7071067811865476, 0.7071067811865476)), ('+-', (0.7071067811865476, -0.7071067811865476)),
        ('-+', (-0.7071067811865476, 0.7071067811865476)), ('--', (-0.7071067811865476, -0.7071067811865476))]


def coverage_at(imgs, p):
    return sum(m for (cx, cy, th, m) in imgs if in_square(p, (cx, cy), th))


def tight_points(imgs, t, pitch, thr):
    """grid points with coverage >= thr (a coarse locator of the tight set)"""
    import numpy as np
    n = int(round(t / pitch))
    xs = (np.arange(n) + 0.5) * pitch
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    cov = np.zeros_like(X)
    for (cx, cy, th, m) in imgs:
        ct, st = math.cos(th), math.sin(th)
        dx = X - cx
        dy = Y - cy
        u = np.abs(dx * ct + dy * st)
        v = np.abs(-dx * st + dy * ct)
        cov += m * ((u <= 0.5) & (v <= 0.5))
    idx = np.argwhere(cov >= thr)
    return [((float(xs[i]), float(xs[j])), float(cov[i, j])) for i, j in idx]


def anchor_scan(imgs, p, t, nfrac, dirs=DIRS, nth=1024):
    """best anchor clique K(p, A) over directions and offsets; returns list of records"""
    out = []
    for dname, d in dirs:
        em, _ = eps_max(p, d, t, nth=nth)
        if em <= 0:
            continue
        for k in range(1, nfrac + 1):
            eps = em * k / (nfrac + 1.0)
            r, th = rho_star(p, d, eps, t, nth=nth, refine=2)
            if r == INF or r >= 0.5:
                continue
            rho = r * (1 + 1e-7) + 1e-12
            tot, cov, extra, nex = clique_mass(imgs, p, d, eps, rho)
            out.append(dict(mass=tot, cov=cov, extra=extra, nex=nex, p=p, dir=dname,
                            eps=eps, eps_max=em, rho=rho, theta=math.degrees(th)))
    return out


def cmd_measure(args):
    t = args.t
    imgs, poses = load_measure(args.support, t)
    tot = sum(m for (_, _, _, m) in imgs)
    print(f"# {args.support}: {len(poses)} poses, {len(imgs)} images, mass {tot:.9f}, t = {t}")
    t0 = time.time()
    cand = tight_points(imgs, t, args.pitch, args.thr)
    print(f"# {len(cand)} grid points at pitch {args.pitch} with coverage >= {args.thr} "
          f"({time.time()-t0:.1f} s)")
    if args.max_pts and len(cand) > args.max_pts:
        cand.sort(key=lambda z: -z[1])
        cand = cand[:args.max_pts]
    recs = []
    for (p, cv) in cand:
        recs += anchor_scan(imgs, p, t, args.nfrac)
    recs.sort(key=lambda r: -r['mass'])
    print(f"{'K mass':>10} {'cov(p)':>10} {'extra':>10} {'#ex':>4} {'p':>22} {'dir':>4} "
          f"{'eps':>9} {'epsmax':>9} {'rho':>9} {'th*':>8} {'chk':>10}")
    shown = 0
    for r in recs:
        if shown >= args.top:
            break
        chk = ''
        if args.verify and r['mass'] > 1.0:
            nf, w = brute_transversal(r['p'], dict(DIRS)[r['dir']], r['eps'], r['rho'], t,
                                      n=args.nbrute)
            npair = clique_pair_check(imgs, r['p'], dict(DIRS)[r['dir']], r['eps'], r['rho'])
            chk = f"{'OK' if nf == 0 else 'FAIL'}/{npair}"
            r['brute_fail'] = nf
            r['pair_fail'] = npair
        print(f"{r['mass']:10.6f} {r['cov']:10.6f} {r['extra']:10.6f} {r['nex']:4d} "
              f"({r['p'][0]:9.5f},{r['p'][1]:9.5f}) {r['dir']:>4} {r['eps']:9.5f} "
              f"{r['eps_max']:9.5f} {r['rho']:9.5f} {r['theta']:8.3f} {chk:>10}")
        shown += 1
    best = recs[0] if recs else None
    if best:
        print(f"\n# best anchor-clique mass {best['mass']:.9f}  (violation "
              f"{best['mass']-1:+.6f}); measure mass {tot:.6f}")
    json.dump(recs[:400], open(os.path.join(RUNS, f'clique_family_{args.tag}.json'), 'w'), indent=1)


def clique_pair_check(imgs, p, d, eps, rho, tol=1e-9):
    """number of member pairs of K(p, A) that do NOT closed-intersect (should be 0)"""
    dp = (-d[1], d[0])
    m = (p[0] + eps * d[0], p[1] + eps * d[1])
    a0 = (m[0] - rho * dp[0], m[1] - rho * dp[1])
    a1 = (m[0] + rho * dp[0], m[1] + rho * dp[1])
    mem = []
    for (cx, cy, th, mm) in imgs:
        if mm <= 0:
            continue
        if in_square(p, (cx, cy), th) or (in_square(a0, (cx, cy), th) and in_square(a1, (cx, cy), th)):
            mem.append((cx, cy, th))
    bad = 0
    for i in range(len(mem)):
        for j in range(i + 1, len(mem)):
            if not sat_meet(mem[i][:2], mem[i][2], mem[j][:2], mem[j][2], tol=tol):
                bad += 1
    return bad


def cmd_kset(args):
    """volume of K(p) and of K(p) \\ P_p in pose space (cx, cy, theta)"""
    import numpy as np
    t = args.t
    print(f"# K(p) = poses meeting every admissible square through p; container t = {t}")
    print(f"# outer approximation from {args.nth} angle samples of U(p); "
          f"theta' grid {args.nthp} over [0, 90)")
    print(f"{'p':>22} {'vol P_p':>12} {'vol K(p)':>12} {'extra vol':>12} {'extra/P':>10} "
          f"{'th-range of extra (deg)':>26}")
    for p in args.points:
        E = extreme_squares(p, t, nth=args.nth)
        NN = k_normals(E, args.nth)
        volP = volK = 0.0
        ths = []
        dthp = (math.pi / 2) / args.nthp
        for k in range(args.nthp):
            thp = (k + 0.5) * dthp
            aP = poly_area(p_region(p, t, thp))
            aK = poly_area(k_region(p, t, thp, E, NN))
            volP += aP * dthp
            volK += aK * dthp
            if aK - aP > 1e-9:
                ths.append(math.degrees(thp))
        rng = f"[{min(ths):.2f}, {max(ths):.2f}]" if ths else "(empty)"
        print(f"({p[0]:9.5f},{p[1]:9.5f}) {volP:12.6e} {volK:12.6e} {volK-volP:12.6e} "
              f"{(volK-volP)/volP if volP else 0:10.5f} {rng:>26}")


def q_region(A, t, thp):
    """{c' admissible : A subset S'(c', th')} as a convex polygon"""
    w2 = (abs(math.cos(thp)) + abs(math.sin(thp))) / 2
    poly = [(w2, w2), (t - w2, w2), (t - w2, t - w2), (w2, t - w2)]
    ct, st = math.cos(thp), math.sin(thp)
    for (ax, ay) in A:
        for (ex, ey) in ((ct, st), (-st, ct)):
            ae = ax * ex + ay * ey
            poly = clip(poly, ex, ey, ae + 0.5)
            poly = clip(poly, -ex, -ey, -(ae - 0.5))
            if not poly:
                return []
    return poly


def anchor_volume(p, A, t, nthp=180):
    """(vol Q_A, vol Q_A \\ P_p) in pose space"""
    volQ = volX = 0.0
    dthp = (math.pi / 2) / nthp
    for k in range(nthp):
        thp = (k + 0.5) * dthp
        q = q_region(A, t, thp)
        if not q:
            continue
        aq = poly_area(q)
        pr = p_region(p, t, thp)
        inter = poly_intersect([q, pr]) if pr else []
        volQ += aq * dthp
        volX += (aq - poly_area(inter)) * dthp
    return volQ, volX


def cmd_anchorvol(args):
    """certifiable anchor family vs the maximal K(p): how much of K(p) is reachable"""
    import numpy as np
    t = args.t
    print(f"# vertical-segment anchors A = {{px+eps}} x [py-rho, py+rho], rho = (1+s) kappa eps,")
    print(f"#   kappa = px/sqrt(1-px^2);  vol = Lebesgue measure in pose space (cx, cy, theta)")
    print(f"{'p':>20} {'eps':>9} {'rho':>9} {'vol Q_A':>11} {'vol Q_A\\P_p':>12} "
          f"{'vol K\\P_p':>11} {'frac of K':>10}")
    for p in args.points:
        px = p[0]
        kappa = px / math.sqrt(1 - px * px)
        E = extreme_squares(p, t, nth=args.nth)
        NN = k_normals(E, args.nth)
        volP = volK = 0.0
        dthp = (math.pi / 2) / args.nthp
        for k in range(args.nthp):
            thp = (k + 0.5) * dthp
            volP += poly_area(p_region(p, t, thp)) * dthp
            volK += poly_area(k_region(p, t, thp, E, NN)) * dthp
        extraK = volK - volP
        best = None
        for frac in np.linspace(0.05, 0.98, args.neps):
            eps = frac * (1 - px)
            rho = kappa * eps * (1 + args.slack)
            if rho >= 0.5:
                continue
            A = [(px + eps, p[1] - rho), (px + eps, p[1] + rho)]
            vq, vx = anchor_volume(p, A, t, nthp=args.nthp)
            if best is None or vx > best[3]:
                best = (eps, rho, vq, vx)
        if best is None:
            continue
        eps, rho, vq, vx = best
        print(f"({p[0]:8.5f},{p[1]:8.5f}) {eps:9.5f} {rho:9.5f} {vq:11.5e} {vx:12.5e} "
              f"{extraK:11.5e} {vx/extraK if extraK > 0 else 0:10.4f}")


def cmd_kclique(args):
    """is K(p) itself pairwise closed-intersecting?  (random pairs from its interior)"""
    import numpy as np
    import random
    t = args.t
    rng = random.Random(11)
    for p in args.points:
        E = extreme_squares(p, t, nth=args.nth)
        NN = k_normals(E, args.nth)
        mem = []
        for _ in range(args.n):
            thp = rng.uniform(0, math.pi / 2)
            poly = k_region(p, t, thp, E, NN)
            if not poly:
                continue
            xs = [q[0] for q in poly]
            ys = [q[1] for q in poly]
            for _try in range(30):
                cx = rng.uniform(min(xs), max(xs))
                cy = rng.uniform(min(ys), max(ys))
                if point_in_poly((cx, cy), poly):
                    mem.append((cx, cy, thp))
                    break
        bad = 0
        npair = 0
        for i in range(len(mem)):
            for j in range(i + 1, len(mem)):
                npair += 1
                if not sat_meet(mem[i][:2], mem[i][2], mem[j][:2], mem[j][2], tol=-1e-12):
                    bad += 1
        print(f"p = ({p[0]:.5f}, {p[1]:.5f}): {len(mem)} sampled members of K(p), {npair} pairs, "
              f"{bad} NON-intersecting")


def point_in_poly(q, poly):
    n = len(poly)
    area = sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1] for i in range(n))
    sgn = 1.0 if area > 0 else -1.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if sgn * ((x2 - x1) * (q[1] - y1) - (y2 - y1) * (q[0] - x1)) < -1e-12:
            return False
    return True


# ---------------------------------------------------------------- regularising a graph clique
def square_poly(c, th):
    return poly_sq(c, th)


def arrangement_points(members, ngrid=400, t=3.99):
    """candidate common points: corners of the members and pairwise edge crossings (float)"""
    import numpy as np
    pts = []
    for (cx, cy, th) in members:
        pts += poly_sq((cx, cy), th)
    return np.array(pts)


def contains_np(IM, pts, tol=1e-12):
    """(n,) mask: which of the squares IM (n,4: cx,cy,th,mass) contain ALL the given points"""
    import numpy as np
    ct, st = np.cos(IM[:, 2]), np.sin(IM[:, 2])
    m = np.ones(len(IM), dtype=bool)
    for (qx, qy) in pts:
        dx = qx - IM[:, 0]
        dy = qy - IM[:, 1]
        m &= (np.abs(dx * ct + dy * st) <= 0.5 + tol) & (np.abs(-dx * st + dy * ct) <= 0.5 + tol)
    return m


def meets_np(IM, A, tol=1e-12):
    """(n,) mask: which of the squares IM meet the convex polygon A (SAT, closed)"""
    import numpy as np
    Av = np.asarray(A, dtype=float).reshape(-1, 2)
    ct, st = np.cos(IM[:, 2]), np.sin(IM[:, 2])
    sep = np.zeros(len(IM), dtype=bool)
    # axes: the square's two normals
    for (ex, ey) in ((ct, st), (-st, ct)):
        cen = IM[:, 0] * ex + IM[:, 1] * ey
        hA1 = np.max(Av[:, 0] * ex[:, None] + Av[:, 1] * ey[:, None], axis=1)
        hA0 = np.min(Av[:, 0] * ex[:, None] + Av[:, 1] * ey[:, None], axis=1)
        sep |= (cen - 0.5 > hA1 + tol) | (cen + 0.5 < hA0 - tol)
    # axes: A's edge normals
    m = Av.shape[0]
    if m >= 2:
        rng = range(m) if m > 2 else range(1)
        for i in rng:
            j = (i + 1) % m
            ex, ey = -(Av[j, 1] - Av[i, 1]), (Av[j, 0] - Av[i, 0])
            nn = math.hypot(ex, ey)
            if nn < 1e-14:
                continue
            ex /= nn
            ey /= nn
            hA1 = float(np.max(Av[:, 0] * ex + Av[:, 1] * ey))
            hA0 = float(np.min(Av[:, 0] * ex + Av[:, 1] * ey))
            r = 0.5 * (np.abs(ct * ex + st * ey) + np.abs(-st * ex + ct * ey))
            cen = IM[:, 0] * ex + IM[:, 1] * ey
            sep |= (cen - r > hA1 + tol) | (cen + r < hA0 - tol)
    return ~sep


def kmass_np(IM, thru_p, A):
    """mass of K(p, A) with thru_p the precomputed mask {p in S}"""
    import numpy as np
    part1 = thru_p & meets_np(IM, A)
    part2 = (~thru_p) & contains_np(IM, A)
    sel = part1 | part2
    return float(IM[sel, 3].sum()), int(part1.sum()), int(part2.sum())


def kmass(imgs, p, A, tol=0.0):
    """mass, on a finitely supported measure, of the anchor clique

        K(p, A) = { S : p in S and S meets A }  u  { S : A subset S }

    which is a clique for EVERY point p and every nonempty compact convex A:
    two members of the first part share p; two of the second share A; and a mixed pair
    S (first), S' (second) has  S n S' contains S n A, nonempty."""
    tot = 0.0
    n1 = n2 = 0
    for (cx, cy, th, m) in imgs:
        if in_square(p, (cx, cy), th):
            if polys_meet(poly_sq((cx, cy), th), A, tol=tol):
                tot += m
                n1 += 1
        elif all(in_square(v, (cx, cy), th) for v in A):
            tot += m
            n2 += 1
    return tot, n1, n2


def regularise(imgs, members, mass, t, nth=1024, shrink=(1.0, 0.5, 0.25, 0.1, 0.0)):
    """turn a graph clique into an anchor clique K(p, A).

    p = a point in as much of the clique's mass as possible; E = the members that miss p;
    A = the intersection of a greedy sub-family of E with nonempty intersection (optionally
    shrunk towards its centroid, which enlarges {S : A subset S} at the cost of the filter
    'S meets A' on the point side).  K(p, A) is a clique whatever A is."""
    import numpy as np
    M = np.asarray(members, dtype=float)
    w = np.asarray(mass, dtype=float)
    cand = arrangement_points(M, t=t)
    best = None
    for q in cand:
        s = 0.0
        for i in range(len(M)):
            if in_square(q, M[i, :2], M[i, 2]):
                s += w[i]
        if best is None or s > best[0]:
            best = (s, q)
    p = (float(best[1][0]), float(best[1][1]))
    inC = np.array([in_square(p, M[i, :2], M[i, 2]) for i in range(len(M))])
    C = M[inC]
    E = M[~inC]
    wE = w[~inC]
    if len(C) == 0 or len(E) == 0:
        return None
    # greedy: add members of E (heaviest first) while the intersection stays nonempty, and
    # evaluate the anchor clique after every addition -- the trade-off is that a SMALLER A
    # is contained in more squares but is met by fewer of the squares through p
    order = np.argsort(-wE)
    A = None
    used = []
    out = []
    for i in order:
        cand_poly = poly_sq(E[i, :2], E[i, 2]) if A is None else \
            poly_intersect([A, poly_sq(E[i, :2], E[i, 2])])
        if not cand_poly or poly_area(cand_poly) <= 1e-13:
            continue
        A = cand_poly
        used.append(int(i))
        tot, n1, n2 = kmass(imgs, p, A)
        out.append(dict(p=p, A=[tuple(v) for v in A], step=len(used), mass=tot,
                        n_point=n1, n_anchor=n2, nC=int(len(C)), nE=int(len(E))))
    if not out:
        return None
    for r in out:
        ok, marg, thw, wit = is_transversal_np(p, r['A'], t, nth=nth)
        r['transversal'] = bool(ok)
        r['margin'] = float(marg)
    out.sort(key=lambda r: -r['mass'])
    return out


def anchor_separate(IM, p, nseed=8, maxstep=40, tol=1e-12):
    """best anchor clique K(p, A) for the finite measure IM (n, 4).

    For a fixed subset T of the images that miss p, the best anchor is A = intersection of T
    (the largest one: it maximises the 'S meets A' part).  So the search is a greedy over T,
    started from each of the `nseed` heaviest images that miss p, evaluating K after every
    addition.  Returns the best record."""
    import numpy as np
    thru = contains_np(IM, [p], tol=tol)
    E = np.nonzero(~thru)[0]
    if len(E) == 0:
        return None
    E = E[np.argsort(-IM[E, 3])]
    best = None
    for s in range(min(nseed, len(E))):
        A = poly_sq(IM[E[s], :2], IM[E[s], 2])
        used = [int(E[s])]
        m, n1, n2 = kmass_np(IM, thru, A)
        if best is None or m > best['mass']:
            best = dict(p=p, A=[tuple(v) for v in A], mass=m, n_point=n1, n_anchor=n2,
                        step=1, seed=int(E[s]))
        for step in range(2, maxstep + 1):
            # add the heaviest image not yet used whose square still meets A
            gain = None
            for i in E:
                if int(i) in used:
                    continue
                nxt = poly_intersect([A, poly_sq(IM[i, :2], IM[i, 2])])
                if not nxt or poly_area(nxt) <= 1e-14:
                    continue
                gain = (int(i), nxt)
                break
            if gain is None:
                break
            used.append(gain[0])
            A = gain[1]
            m, n1, n2 = kmass_np(IM, thru, A)
            if m > best['mass']:
                best = dict(p=p, A=[tuple(v) for v in A], mass=m, n_point=n1, n_anchor=n2,
                            step=step, seed=int(E[s]))
    return best


def cmd_separate(args):
    """scan candidate anchor points p and report the best anchor clique of the measure"""
    import numpy as np
    t = args.t
    imgs, poses = load_measure(args.support, t)
    IM = np.array([(a, b, c, d) for (a, b, c, d) in imgs if d > 1e-12])
    tot = float(IM[:, 3].sum())
    print(f"# {args.support}: {len(IM)} images, mass {tot:.9f}, t = {t}")
    cand = tight_points(imgs, t, args.pitch, args.thr)
    print(f"# {len(cand)} candidate points with coverage >= {args.thr} at pitch {args.pitch}")
    cand.sort(key=lambda z: -z[1])
    if args.max_pts and len(cand) > args.max_pts:
        step = max(1, len(cand) // args.max_pts)
        cand = cand[::step][:args.max_pts]
    recs = []
    t0 = time.time()
    for k, (p, cv) in enumerate(cand):
        r = anchor_separate(IM, p, nseed=args.nseed, maxstep=args.maxstep)
        if r is not None:
            r['cov'] = cv
            recs.append(r)
        if k % 20 == 0:
            print(f"  ... {k}/{len(cand)}  best so far "
                  f"{max((z['mass'] for z in recs), default=0):.6f}  ({time.time()-t0:.0f}s)",
                  flush=True)
    recs.sort(key=lambda r: -r['mass'])
    print(f"{'K mass':>10} {'cov(p)':>9} {'#thru p':>8} {'#anchor':>8} {'|A|':>4} {'diamA':>9} "
          f"{'p':>22}")
    for r in recs[:args.top]:
        A = r['A']
        dia = max(math.hypot(u[0] - v[0], u[1] - v[1]) for u in A for v in A) if len(A) > 1 else 0.0
        print(f"{r['mass']:10.6f} {r['cov']:9.6f} {r['n_point']:8d} {r['n_anchor']:8d} "
              f"{len(A):4d} {dia:9.5f} ({r['p'][0]:9.5f},{r['p'][1]:9.5f})")
    if recs:
        print(f"\n# BEST anchor-clique mass on this measure: {recs[0]['mass']:.9f} "
              f"(violation {recs[0]['mass']-1:+.6f})")
    json.dump(recs[:100], open(os.path.join(RUNS, f'clique_sep_{args.tag}.json'), 'w'), indent=1)


def cmd_regular(args):
    import numpy as np
    t = args.t
    imgs, poses = load_measure(args.support, t)
    tot = sum(m for (_, _, _, m) in imgs)
    print(f"# {args.support}: {len(poses)} poses, {len(imgs)} images, mass {tot:.9f}, t = {t}")
    im = np.array([(a, b, c, d) for (a, b, c, d) in imgs if d > 1e-12])
    print(f"# {len(im)} images with positive mass; building the closed-intersection graph")
    t0 = time.time()
    adj = closed_adj(im[:, :3])
    bw, cl = max_weight_clique(adj, im[:, 3], tlimit=args.tl)
    print(f"# max-mass clique of the SUPPORT graph: {bw:.6f} on {len(cl)} images "
          f"({time.time()-t0:.1f} s)")
    res = regularise(imgs, im[cl, :3], im[cl, 3], t, nth=args.nth)
    if res is None:
        print("# regularisation failed (no anchor clique found)")
        return
    r0 = res[0]
    print(f"# anchor p = ({r0['p'][0]:.6f}, {r0['p'][1]:.6f});  the graph clique splits into "
          f"{r0['nC']} images through p and {r0['nE']} others")
    print(f"{'step':>5} {'|A| vtx':>8} {'diam A':>9} {'K mass':>10} {'#thru p':>8} "
          f"{'#anchor':>8} {'transversal':>12} {'margin':>11}")
    for r in res[:20]:
        A = r['A']
        dia = max(math.hypot(u[0] - v[0], u[1] - v[1]) for u in A for v in A) if len(A) > 1 else 0.0
        print(f"{r['step']:5d} {len(A):8d} {dia:9.5f} {r['mass']:10.6f} {r['n_point']:8d} "
              f"{r['n_anchor']:8d} {str(r['transversal']):>12} {r['margin']:11.3e}")
    print(f"# best anchor-clique mass {res[0]['mass']:.6f} vs graph clique {bw:.6f}; "
          f"violation {res[0]['mass']-1:+.6f}")
    json.dump([{k: v for k, v in r.items() if k != 'A'} | {'A': [list(v) for v in r['A']]}
               for r in res], open(os.path.join(RUNS, 'clique_family_regular.json'), 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('rho')
    a.add_argument('--t', type=float, default=3.99)
    b = sub.add_parser('measure')
    b.add_argument('--t', type=float, default=3.99)
    b.add_argument('--support', default=os.path.join(MAIN_RUNS, 'dual_exact_3.99_support.txt'))
    b.add_argument('--top', type=int, default=30)
    b.add_argument('--pitch', type=float, default=0.005)
    b.add_argument('--thr', type=float, default=0.97)
    b.add_argument('--nfrac', type=int, default=6)
    b.add_argument('--max-pts', type=int, default=300)
    b.add_argument('--tag', default='PA2')
    b.add_argument('--verify', action='store_true')
    b.add_argument('--nbrute', type=int, default=30000)
    c = sub.add_parser('kset')
    c.add_argument('--t', type=float, default=3.99)
    c.add_argument('--nth', type=int, default=512)
    c.add_argument('--nthp', type=int, default=180)
    c.add_argument('--pts', default='')
    e = sub.add_parser('anchorvol')
    e.add_argument('--t', type=float, default=3.99)
    e.add_argument('--nth', type=int, default=180)
    e.add_argument('--nthp', type=int, default=180)
    e.add_argument('--neps', type=int, default=20)
    e.add_argument('--slack', type=float, default=0.05)
    e.add_argument('--pts', default='')
    f = sub.add_parser('kclique')
    f.add_argument('--t', type=float, default=3.99)
    f.add_argument('--nth', type=int, default=180)
    f.add_argument('--n', type=int, default=120)
    f.add_argument('--pts', default='')
    g = sub.add_parser('regular')
    g.add_argument('--t', type=float, default=3.99)
    g.add_argument('--support', default=os.path.join(MAIN_RUNS, 'dual_exact_3.99_support.txt'))
    g.add_argument('--tl', type=float, default=120.0)
    g.add_argument('--inflate', type=float, default=0.0)
    g.add_argument('--nth', type=int, default=1024)
    s2 = sub.add_parser('separate')
    s2.add_argument('--t', type=float, default=3.99)
    s2.add_argument('--support', default=os.path.join(MAIN_RUNS, 'dual_exact_3.99_support.txt'))
    s2.add_argument('--pitch', type=float, default=0.01)
    s2.add_argument('--thr', type=float, default=0.95)
    s2.add_argument('--max-pts', type=int, default=200)
    s2.add_argument('--nseed', type=int, default=8)
    s2.add_argument('--maxstep', type=int, default=30)
    s2.add_argument('--top', type=int, default=25)
    s2.add_argument('--tag', default='PA2')
    args = ap.parse_args()
    if args.cmd == 'separate':
        cmd_separate(args)
        return
    if args.cmd == 'regular':
        cmd_regular(args)
        return
    DEFPTS = [(0.99, 1.995), (0.95, 1.995), (0.9, 1.995), (0.8, 1.995), (0.6, 1.995),
              (0.5, 1.995), (0.99, 0.99), (0.95, 0.95), (0.9, 0.9)]
    if args.cmd in ('kset', 'anchorvol', 'kclique'):
        if getattr(args, 'pts', ''):
            args.points = [tuple(float(z) for z in s.split(',')) for s in args.pts.split(';')]
        else:
            args.points = DEFPTS
    if args.cmd == 'rho':
        cmd_rho(args)
        return
    if args.cmd == 'measure':
        cmd_measure(args)
        return
    if args.cmd == 'anchorvol':
        cmd_anchorvol(args)
        return
    if args.cmd == 'kclique':
        cmd_kclique(args)
        return
    if args.cmd == 'kset':
        if args.pts:
            args.points = [tuple(float(z) for z in s.split(',')) for s in args.pts.split(';')]
        else:
            args.points = [(0.99, 1.995), (0.99, 1.275), (0.95, 1.995), (0.9, 1.995),
                           (0.8, 1.995), (0.6, 1.995), (0.5, 1.995), (1.2, 1.995),
                           (1.0, 1.995), (0.99, 0.99), (0.95, 0.95), (0.9, 0.9),
                           (0.99, 1.0), (0.98, 0.98), (1.5, 1.995), (1.995, 1.995)]
        cmd_kset(args)


if __name__ == '__main__':
    main()
