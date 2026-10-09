"""Morphs between s(110) minima for the landscape explainer (10-08).

match(A, B): square assignment minimising total corner movement (each square's 4 corners may be relabelled, i.e. B's
angle is taken modulo 90 deg nearest to A's).  frames(A, B, perm, K): centres and angles interpolated linearly.
side_needed(frame): smallest side for that frame to be feasible after (a) uniform spreading only (slp2.repair) and
(b) a short tethered squeeze (slp2 with a capped number of small trust-region steps), an upper bound on the expansion
a path through nearby configurations needs.  Packings: (s, [(x, y, theta_rad)]), the slp2 convention.
"""
import math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
from scipy.optimize import linear_sum_assignment
import slp2

Q = math.pi / 2


def load(path):
    s, sq = slp2.load(path)
    return s, [(x, y, t) for x, y, t in sq]


def wrapq(d):
    """angle difference reduced to (-pi/4, pi/4] (a square is invariant under quarter turns)."""
    return d - Q * round(d / Q)


def corners(sq):
    P = np.array(sq)
    c, s = np.cos(P[:, 2]), np.sin(P[:, 2])
    V = np.array([[.5, .5], [-.5, .5], [-.5, -.5], [.5, -.5]])
    return P[:, None, :2] + np.stack([V[:, 0] * c[:, None] - V[:, 1] * s[:, None],
                                      V[:, 0] * s[:, None] + V[:, 1] * c[:, None]], -1)


def match(A, B):
    """perm[i] = index in B of the square matched to A's square i; total corner movement."""
    CA, CB = corners(A[1]), corners(B[1])
    best = None
    for k in range(4):
        d = np.linalg.norm(CA[:, None] - np.roll(CB, k, axis=1)[None], axis=-1).sum(-1)
        best = d if best is None else np.minimum(best, d)
    r, c = linear_sum_assignment(best)
    return c, float(best[r, c].sum()), best[r, c]


def frame(A, B, perm, t):
    sA, a = A
    sB, b = B
    out = []
    for i, (x, y, th) in enumerate(a):
        x2, y2, th2 = b[perm[i]]
        out.append((x + t * (x2 - x), y + t * (y2 - y), th + t * wrapq(th2 - th)))
    return sA + t * (sB - sA), out


def side_repair(s, sq):
    return slp2.repair(s, sq)


def side_relaxed(s, sq, steps=6, R=2e-3):
    """repair, then <= `steps` slp2 iterations of trust radius <= R: moves each coordinate by at most steps*R."""
    s1, sq1 = slp2.repair(s, sq)
    if not math.isfinite(s1):
        return s1, sq1
    s2, sq2, _ = slp2.slp2(s1, sq1, R=R, iters=steps, milp=False)
    return s2, sq2


# ------------------------------------------------------------------------------------- feasible walk at fixed side
import inc
from scipy.optimize import linprog


def _correct(s, sq, r):
    try:                    # inc.correct takes g.min() of the contact rows: empty when nothing is within reach
        return inc.correct(s, sq, r)
    except ValueError:
        return sq
from scipy.sparse import csr_matrix, hstack as shstack, vstack as svstack, identity

ANG_W = 0.5          # objective weight of an angle (rad) vs a centre coordinate: a corner moves ~0.7 rad^-1


def match2(A, B):
    """assignment by squared corner movement (favours many small moves over one long one: chains, not jumps)."""
    CA, CB = corners(A[1]), corners(B[1])
    best = None
    for k in range(4):
        d = (np.linalg.norm(CA[:, None] - np.roll(CB, k, axis=1)[None], axis=-1) ** 2).sum(-1)
        best = d if best is None else np.minimum(best, d)
    r, c = linear_sum_assignment(best)
    return c


def target(A, B, perm):
    a, b = A[1], B[1]
    return [(b[perm[i]][0], b[perm[i]][1], th + wrapq(b[perm[i]][2] - th)) for i, (_, _, th) in enumerate(a)]


def dist(sq, tg):
    return sum(abs(x - u) + abs(y - v) + ANG_W * abs(t - w) for (x, y, t), (u, v, w) in zip(sq, tg))


def step(s, sq, tg, R):
    """LP: min sum |sq + z - tg| (weighted)  s.t. linearised non-overlap at side s (fixed), |z| <= R."""
    n = len(sq)
    nz = 3 * n + 1
    rows, disj = inc.model_fast(s, sq, 3 * R)
    A, g = inc._mat(rows, nz)
    blocks, gs = [A], [g]
    d = np.array([v for (x, y, t), (u, w, h) in zip(sq, tg) for v in (u - x, w - y, h - t)])
    zdes = np.append(np.clip(d, -R, R), 0.0)
    for _, alts in disj:
        # corner-corner pair: the separating line that the step toward the target violates least
        q = max(range(len(alts)), key=lambda q: inc._alt_slack(alts[q], zdes))
        Aq, gq = inc._mat(alts[q], nz)
        blocks.append(Aq); gs.append(gq)
    A = svstack(blocks).tocsr()[:, :3 * n]
    g = np.maximum(np.concatenate(gs), 0.0)
    wt = np.tile([1.0, 1.0, ANG_W], n)
    m = 3 * n
    # variables [z (m), e (m)]: e >= |d - z|;  min wt.e
    I = identity(m, format='csr')
    Aub = svstack([shstack([-A, csr_matrix((A.shape[0], m))]),     # -A z <= g
                   shstack([-I, -I]),                               # d - z <= e
                   shstack([I, -I])]).tocsr()                       # z - d <= e
    bub = np.concatenate([g, -d, d])
    c = np.concatenate([np.zeros(m), wt])
    bnd = [(-R, R)] * m + [(0, None)] * m
    res = linprog(c, A_ub=Aub, b_ub=bub, bounds=bnd, method='highs')
    if res.status != 0:
        return None
    z = res.x[:m]
    return [(x + z[3 * i], y + z[3 * i + 1], t + z[3 * i + 2]) for i, (x, y, t) in enumerate(sq)]


def walk(s, A, B, perm, R=0.02, max_steps=600, tol=1e-9, keep=None):
    """Feasible path from A toward B at side s.  Returns (reached, final distance, path) -- path = list of sq."""
    tg = target(A, B, perm)
    sq = list(A[1])
    path = [sq]
    d0 = dist(sq, tg)
    stall = 0
    r = R
    for it in range(max_steps):
        nxt = step(s, sq, tg, r)
        if nxt is not None:
            c = _correct(s, nxt, r)
            nxt = c if c is not None else nxt
        if nxt is None or slp2.worst(s, nxt) > tol:
            r /= 2
            if r < 1e-5:
                break
            continue
        d1 = dist(nxt, tg)
        stall = stall + 1 if d1 > d0 - 1e-4 else 0
        sq, d0 = nxt, d1
        path.append(sq)
        r = min(R, r * 1.5)
        if d0 < 1e-3 or stall > 25:
            break
    reached = d0 < 1e-3
    if reached:
        path.append(list(B[1]) if False else tg)
    return reached, d0, path


def barrier(A, B, perm=None, lo=None, hi=0.3, res=2e-4, log=print, **kw):
    """smallest side s (to res) at which walk() connects A to B: an upper bound on the minimal expansion."""
    perm = match2(A, B) if perm is None else perm
    base = max(A[0], B[0])
    lo = base if lo is None else lo
    ok_hi = None
    h = 0.01
    while h <= hi:
        r, d, p = walk(base + h, A, B, perm, **kw)
        log(f'  probe s = {base + h:.5f}: reached {r} (dist {d:.4f}, {len(p)} steps)')
        if r:
            ok_hi = (base + h, p)
            break
        lo = base + h
        h *= 2
    if ok_hi is None:
        return None, None, perm
    shi, path = ok_hi
    while shi - lo > res:
        mid = (lo + shi) / 2
        r, d, p = walk(mid, A, B, perm, **kw)
        log(f'  bisect s = {mid:.5f}: reached {r} (dist {d:.4f})')
        if r:
            shi, path = mid, p
        else:
            lo = mid
    return shi, path, perm


def walk_carrot(s, A, B, perm, K=40, sub=8, R=0.02, tol=1e-9, finish=200):
    """Feasible path at side s that follows the straight-line morph: for t_k = k/K, up to `sub` feasible LP steps toward
    the interpolated frame (a moving target pushes blocking squares aside), then up to `finish` steps toward B.
    Returns (reached, final distance to B, path)."""
    fa, fb = s / A[0], s / B[0]
    tg = [(x * fb, y * fb, t) for x, y, t in target(A, B, perm)]
    a = [(x * fa, y * fa, t) for x, y, t in A[1]]       # spread A to fill the box: a feasible continuous motion
    sq = list(a)
    path = [sq]
    r = R

    def go(ref, nmax):
        nonlocal sq, r
        for _ in range(nmax):
            nxt = step(s, sq, ref, r)
            if nxt is not None:
                c = _correct(s, nxt, r)
                nxt = c if c is not None else nxt
            if nxt is None or slp2.worst(s, nxt) > tol:
                r /= 2
                if r < 1e-5:
                    r = 1e-5
                    return False
                continue
            moved = dist(nxt, sq)
            sq = nxt
            path.append(sq)
            r = min(R, r * 1.5)
            if dist(sq, ref) < 1e-4 or moved < 1e-7:
                return True
        return True

    for k in range(1, K + 1):
        t = k / K
        ref = [(x + t * (u - x), y + t * (v - y), th + t * (w - th)) for (x, y, th), (u, v, w) in zip(a, tg)]
        go(ref, sub)
    stall, d0 = 0, dist(sq, tg)
    for _ in range(finish):
        go(tg, 1)
        d1 = dist(sq, tg)
        stall = stall + 1 if d1 > d0 - 1e-5 else 0
        d0 = d1
        if d0 < 1e-3 or stall > 15:
            break
    return d0 < 1e-3, d0, path
