#!/usr/bin/env python3
"""Float exploration tool (NOT a proof): largest square avoiding a point set.

smax(P, t, theta) = sup of side s such that some square of side s at angle theta lies in [0,t]^2
with no point of P in its interior.  P is unavoidable for closed unit squares in [0,t]^2
iff max_theta smax <= 1.

At fixed theta the side available at centre c is
    min( min_p 2*max_k |<p - c, e_k>| ,  2*dist(c, walls)/w ),  w = |cos| + |sin|,
a piecewise-linear function built from the 4|P| + 4 linear functions  +-2<p - c, e_k>  and the
4 wall functions.  Its maximum is at a vertex where three of them agree.  We enumerate triples of
these functions whose owners (points / walls) are pairwise "close" (a square of side <= SMAX_CAP
cannot touch two owners farther apart than its diameter), solve, and evaluate the true value.
"""
import numpy as np, itertools, math

SMAX_CAP = 1.3          # only sides up to this are searched reliably (enough near s = 1)

_cache = {}

def _triples(P, t):
    """Triples of function indices with pairwise-close owners.  Functions 0..4n-1 belong to point
    i = f // 4; 4n..4n+3 are walls left, right, bottom, top."""
    n = len(P)
    D = SMAX_CAP * math.sqrt(2) + 1e-9
    owners = []
    for i in range(n): owners += [('p', i)] * 4
    owners += [('w', 0), ('w', 1), ('w', 2), ('w', 3)]
    def close(o1, o2):
        if o1 == o2: return True
        if o1[0] == 'p' and o2[0] == 'p':
            (x1, y1), (x2, y2) = P[o1[1]], P[o2[1]]
            return math.hypot(x1 - x2, y1 - y2) <= D
        if o1[0] == 'w' and o2[0] == 'w':
            return {o1[1], o2[1]} not in ({0, 1}, {2, 3}) or t <= D
        p, w = (o1, o2) if o1[0] == 'p' else (o2, o1)
        x, y = P[p[1]]
        d = (x, t - x, y, t - y)[w[1]]
        return d <= D
    m = len(owners)
    T = []
    for i, j, k in itertools.combinations(range(m), 3):
        if close(owners[i], owners[j]) and close(owners[i], owners[k]) and close(owners[j], owners[k]):
            T.append((i, j, k))
    return np.array(T)

def smax_theta(P, t, th, return_c=False, T=None):
    P = np.asarray(P, float)
    if T is None: T = _triples([tuple(p) for p in P], t)
    c, s = math.cos(th), math.sin(th)
    w = abs(c) + abs(s)
    E = np.array([[c, s], [-s, c]])
    L = []  # rows (a1, a2, b): f(c) = a1 cx + a2 cy + b
    for p in P:
        for e in E:
            for sg in (1, -1):
                a = -2 * sg * e
                L.append((a[0], a[1], 2 * sg * (p @ e)))
    L += [(2 / w, 0, 0), (-2 / w, 0, 2 * t / w), (0, 2 / w, 0), (0, -2 / w, 2 * t / w)]
    L = np.array(L)
    Ai = L[T]
    M = np.concatenate([Ai[:, :, :2], -np.ones((len(T), 3, 1))], axis=2)
    rhs = -Ai[:, :, 2]
    det = np.linalg.det(M)
    ok = np.abs(det) > 1e-12
    sol = np.linalg.solve(M[ok], rhs[ok][..., None])[..., 0]
    C = sol[:, :2]
    val = smax_at(P, t, E, w, C)
    i = np.argmax(val)
    return (val[i], C[i]) if return_c else val[i]

def smax_at(P, t, E, w, C):
    D = P[None, :, :] - C[:, None, :]
    proj = np.abs(D @ E.T)
    d = proj.max(axis=2).min(axis=1) * 2
    wall = np.minimum.reduce([C[:, 0], t - C[:, 0], C[:, 1], t - C[:, 1]]) * 2 / w
    return np.minimum(d, wall)

def smax(P, t, nth=180, refine=True, lo=0.0, hi=math.pi / 2, nloc=3):
    """max over theta in [lo, hi] (grid + local refinement of the nloc best grid angles)."""
    T = _triples([tuple(map(float, p)) for p in P], t)
    ths = np.linspace(lo, hi, nth + 1)
    vals = np.array([smax_theta(P, t, th, T=T) for th in ths])
    best_v, best_th = -1, None
    order = np.argsort(-vals)[:nloc] if refine else [int(np.argmax(vals))]
    for b in order:
        th, v = ths[b], vals[b]
        if refine:
            h = (hi - lo) / nth
            for _ in range(30):
                for cth in (th - h, th + h):
                    if lo <= cth <= hi:
                        vc = smax_theta(P, t, cth, T=T)
                        if vc > v: v, th = vc, cth
                h /= 2
        if v > best_v: best_v, best_th = v, th
    return best_v, best_th, vals
