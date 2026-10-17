#!/usr/bin/env python3
"""Same basin, proved: a continuous path from one packing to another along which the side never exceeds a level L.
(10-17; Evan: "'can you move from one to the other without increasing side length' is the real test of a local minimum
and separate basins … Can we make something provably correct there and use it?")

A path is a chain of configurations (n unit squares, centre (x, y), angle θ; box side s).  Between neighbouring nodes
A, B every square moves linearly in (x, y, θ) and the box side linearly, λ in [0, 1]; the configuration at λ is then
scaled about the origin by 1 + η(λ) (η continuous), so the path is continuous.  The step is certified at level L when
the scaled configuration is a packing with side <= L for every λ:

  pairs   for a pair (i, j) whose centres come within √2 during the step, one separating axis valid along it: the
          normal n of an edge of the owner o (4 choices each for o = i, j) and
            g(λ) = n(θ_o) · (c_t - c_o) - 1/2 - h(θ_t - θ_o),    h(φ) = (|cos φ| + |sin φ|) / 2  (t's support on n),
          g >= 0 means separated.  g'' <= q = |w|max Δθ_o² + 2 |Δw| |Δθ_o| + (√2/2) Δφ² on every smooth piece (w = c_t -
          c_o, linear in λ; |n'| = |n''| = 1; -h'' = h Δφ² <= (√2/2) Δφ²) and h's kinks are convex, so -h's are concave
          and only help.  Hence g(λ) >= (1-λ) g(A) + λ g(B) - λ(1-λ) q/2  (f'' <= q  =>  f >= chord - λ(1-λ) q/2).
          Pairs whose centres stay >= √2 apart can't touch (circumscribed discs).
  walls   x - e(θ) >= 0 and s - x - e(θ) >= 0 (and y), e = h(θ): the same bound with q = (√2/2) Δθ².
  scaling scaling the centres by 1 + η adds η n·w >= η (1 + g) to a pair's g (n·w = g + 1/2 + h, h >= 1/2) and
          η (x) >= η (1/2 + m) to a wall margin m.  With shortfall v = max(0, -G(λ)) for the lower bound G, the needed
          side is lerp(s) / (1 - v) (pairs) or lerp(s) / (1 - 2v) (walls).  So the step is at level L iff for every
          constraint, for all λ in [0, 1]:
              pairs  G(λ) + 1 - lerp(s)/L >= 0        walls  G(λ) + (1 - lerp(s)/L)/2 >= 0,
          each a convex quadratic in λ whose minimum on [0, 1] is evaluated in closed form (mpmath, 50 digits; a margin
          of 1e-30 absorbs its rounding).
A step that fails is bisected (midpoint node = linear interpolation; its own overlaps are covered by the scaling).

The chain from a posting P to a minimum M: P, the descent's accepted iterates (slp2.TRACE; every accepted step lowers
the side), the descended packing, M's exact solution.  Two postings P (higher) and Q: P's chain to M, then M to Q's
minimum M' (Q's chain matched to P's labels and box symmetry; a step within M's flat family), then Q's chain reversed.
L = P's level (its side, or slightly above if the posted floats overlap) times 1 + 1e-25.

  pathcert.py posting P.txt [--to Q.txt] [--exact M.exact.txt] [--exact-to M2.exact.txt]
"""
import argparse, json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
sys.path.insert(0, os.path.join(os.path.dirname(PK), 'exact'))
import mpmath
from mpmath import mp, mpf
import numpy as np

mp.dps = 50
SQ2 = mp.sqrt(2) / 2
MARGIN = mpf('1e-30')


# ------------------------------------------------------------------------------------------------ nodes
def read_node(path):
    """(s, [(x, y, θ)]) as mpf from an 'ours' file (fractions or decimals; angles in degrees)"""
    from fractions import Fraction
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n = int(L[0][0])
    f = lambda v: mpf(Fraction(v).numerator) / Fraction(v).denominator if '/' in v else mpf(v)
    return f(L[0][1]), [(f(r[0]), f(r[1]), f(r[2]) * mp.pi / 180) for r in L[1:n + 1]]


def to_mp(node):
    s, sq = node
    return mpf(s), [(mpf(x), mpf(y), mpf(t)) for x, y, t in sq]


def unwrap(prev, node):
    """choose each angle mod 90° nearest the previous node's (a square is invariant under quarter turns)"""
    s, sq = node
    q = mp.pi / 2
    out = []
    for (_, _, t0), (x, y, t) in zip(prev[1], sq):
        out.append((x, y, t - q * mpmath.nint((t - t0) / q)))
    return s, out


def transform(node, k):
    """box symmetry k (bit 4: transpose, 1: flip x, 2: flip y) applied to a node"""
    s, sq = node
    out = []
    for x, y, t in sq:
        if k & 4:
            x, y, t = y, x, mp.pi / 2 - t
        if k & 1:
            x, t = s - x, -t
        if k & 2:
            y, t = s - y, -t
        out.append((x, y, t))
    return s, out


def match(A, B):
    """relabel / transform B onto A: (k, perm) with B' = transform(B, k)[perm[i]] for A's square i"""
    from scipy.optimize import linear_sum_assignment
    a = np.array([(float(x), float(y)) for x, y, _ in A[1]])
    best = None
    for k in range(8):
        b = np.array([(float(x), float(y)) for x, y, _ in transform(B, k)[1]])
        C = np.abs(a[:, None, :] - b[None, :, :]).max(axis=2)
        r, c = linear_sum_assignment(C)
        cost = C[r, c].max()
        if best is None or cost < best[0]:
            best = (cost, k, [int(c[list(r).index(i)]) for i in range(len(a))])
    return best[1], best[2], best[0]


def apply(node, k, perm):
    s, sq = transform(node, k)
    return s, [sq[p] for p in perm]


# ------------------------------------------------------------------------------------------------ one step
def _h(phi):
    return (abs(mpmath.cos(phi)) + abs(mpmath.sin(phi))) / 2


def _qmin(c0, c1, c2):
    """min over [0, 1] of c0 + c1 λ + c2 λ², c2 >= 0"""
    v = min(c0, c0 + c1 + c2)
    if c2 > 0:
        lam = -c1 / (2 * c2)
        if 0 < lam < 1:
            v = min(v, c0 + c1 * lam + c2 * lam * lam)
    return v


def _cond(gA, gB, q, sA, sB, L, wall):
    """min over λ of  G(λ) + k (1 - lerp(s)/L),  G = (1-λ) gA + λ gB - λ(1-λ) q/2,  k = 1 (pair) or 1/2 (wall)"""
    k = mpf(1) / 2 if wall else mpf(1)
    # G = gA + λ (gB - gA) - (q/2) λ + (q/2) λ²;  k (1 - sA/L - λ (sB - sA)/L)
    c0 = gA + k * (1 - sA / L)
    c1 = (gB - gA) - q / 2 - k * (sB - sA) / L
    c2 = q / 2
    return _qmin(c0, c1, c2)


def _gap(o, t, k, sq):
    """g for the axis of owner o's edge normal k (k quarter turns from θ_o) against square t"""
    xo, yo, to = sq[o]
    xt, yt, tt = sq[t]
    a = to + k * mp.pi / 2
    return mpmath.cos(a) * (xt - xo) + mpmath.sin(a) * (yt - yo) - mpf(1) / 2 - _h(tt - to)


def _gapf(o, t, k, sq):
    xo, yo, to = (float(v) for v in sq[o])
    xt, yt, tt = (float(v) for v in sq[t])
    a = to + k * math.pi / 2
    d = tt - to
    return math.cos(a) * (xt - xo) + math.sin(a) * (yt - yo) - 0.5 - 0.5 * (abs(math.cos(d)) + abs(math.sin(d)))


def step_ok(A, B, L):
    """certify the step A -> B at level L.  Returns (ok, worst margin, what)"""
    sA, a = A
    sB, b = B
    n = len(a)
    worst = (mpf('inf'), None)
    # walls
    for i in range(n):
        dth = abs(b[i][2] - a[i][2])
        q = SQ2 * dth * dth
        eA, eB = _h(a[i][2]), _h(b[i][2])
        for c in (0, 1):
            for lo in (True, False):
                mA = (a[i][c] - eA) if lo else (sA - a[i][c] - eA)
                mB = (b[i][c] - eB) if lo else (sB - b[i][c] - eB)
                v = _cond(mA, mB, q, sA, sB, L, True)
                if v < worst[0]:
                    worst = (v, ('wall', i, c, lo))
    # pairs: centres within √2 at some λ (|w(λ)|² is a convex quadratic)
    ca = np.array([(float(x), float(y)) for x, y, _ in a])
    cb = np.array([(float(x), float(y)) for x, y, _ in b])
    for i in range(n):
        wa = ca - ca[i]
        wb = cb - cb[i]
        dw = wb - wa
        A2 = (dw * dw).sum(1); B2 = 2 * (wa * dw).sum(1); C2 = (wa * wa).sum(1)
        lam = np.clip(np.where(A2 > 0, -B2 / (2 * np.where(A2 > 0, A2, 1)), 0), 0, 1)
        dmin = np.minimum(np.minimum(C2, A2 + B2 + C2), A2 * lam * lam + B2 * lam + C2)
        for j in np.nonzero(dmin < 2.0 + 1e-6)[0]:
            j = int(j)
            if j <= i:
                continue
            # pick the axis by a float estimate, then certify it in mp
            cands = []
            for o, t in ((i, j), (j, i)):
                for k in range(4):
                    cands.append((min(_gapf(o, t, k, a), _gapf(o, t, k, b)), o, t, k))
            cands.sort(reverse=True)
            best = None
            for _, o, t, k in cands[:3]:
                wA = mpmath.matrix([a[t][0] - a[o][0], a[t][1] - a[o][1]])
                wB = mpmath.matrix([b[t][0] - b[o][0], b[t][1] - b[o][1]])
                dW = mpmath.norm(wB - wA)
                W = max(mpmath.norm(wA), mpmath.norm(wB))
                dto = abs(b[o][2] - a[o][2])
                dphi = abs((b[t][2] - b[o][2]) - (a[t][2] - a[o][2]))
                q = W * dto * dto + 2 * dW * dto + SQ2 * dphi * dphi
                v = _cond(_gap(o, t, k, a), _gap(o, t, k, b), q, sA, sB, L, False)
                if best is None or v > best[0]:
                    best = (v, ('pair', i, j, o, k))
                if v > MARGIN:
                    break
            if best[0] < worst[0]:
                worst = best
    return worst[0] > MARGIN, worst[0], worst[1]


def lerp(A, B, lam):
    sA, a = A
    sB, b = B
    return sA + lam * (sB - sA), [tuple(p + lam * (q - p) for p, q in zip(u, v)) for u, v in zip(a, b)]


def _axis_row(o, t, k, F, nz):
    """g and its gradient (float) for owner o's edge-normal k against t; F = [(x, y, θ)]"""
    xo, yo, to = F[o]; xt, yt, tt = F[t]
    a = to + k * math.pi / 2
    ca, sa = math.cos(a), math.sin(a)
    phi = tt - to
    h = 0.5 * (abs(math.cos(phi)) + abs(math.sin(phi)))
    dh = 0.5 * (-math.copysign(1, math.cos(phi)) * math.sin(phi) + math.copysign(1, math.sin(phi)) * math.cos(phi))
    g = ca * (xt - xo) + sa * (yt - yo) - 0.5 - h
    row = {3 * t: ca, 3 * t + 1: sa, 3 * o: -ca, 3 * o + 1: -sa,
           3 * o + 2: -sa * (xt - xo) + ca * (yt - yo) + dh, 3 * t + 2: -dh}
    return g, row


def midpoint_lp(A, B, free_side=False):
    """the bisection node by LP (10-17, 261's corner switch): from the linear midpoint, the smallest move (max norm) at
    fixed side such that every near pair satisfies (linearised) the separating line valid at both ends, or, if no line
    is valid at both ends (a corner-corner switch), both ends' lines: then the node is a corner touch, and each half-step
    has a common line.  Walls likewise.  None if the LP fails."""
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    Fa = [(float(x), float(y), float(t)) for x, y, t in A[1]]
    Fb = [(float(x), float(y), float(t)) for x, y, t in B[1]]
    s = (float(A[0]) + float(B[0])) / 2
    F = [tuple((u + v) / 2 for u, v in zip(p, q)) for p, q in zip(Fa, Fb)]
    n = len(F); nz = 3 * n + 2 if free_side else 3 * n + 1   # last: the max-norm bound; free_side: ds before it
    DS = 3 * n if free_side else None
    rows, rhs = [], []
    c = np.array([(x, y) for x, y, _ in F])
    for i in range(n):
        x, y, t = F[i]
        e = 0.5 * (abs(math.cos(t)) + abs(math.sin(t)))
        de = 0.5 * (-math.copysign(1, math.cos(t)) * math.sin(t) + math.copysign(1, math.sin(t)) * math.cos(t))
        hi = {DS: 1.0} if free_side else {}
        for (g, r) in ((x - e, {3 * i: 1, 3 * i + 2: -de}), (s - x - e, {3 * i: -1, 3 * i + 2: -de, **hi}),
                       (y - e, {3 * i + 1: 1, 3 * i + 2: -de}), (s - y - e, {3 * i + 1: -1, 3 * i + 2: -de, **hi})):
            if g < 1e-3:
                rows.append(r); rhs.append(g)
        d2 = ((c - c[i]) ** 2).sum(1)
        for j in np.nonzero(d2 < 2.0 + 1e-3)[0]:
            j = int(j)
            if j <= i:
                continue
            ax = [(o, t2, k) for o, t2 in ((i, j), (j, i)) for k in range(4)]
            ga = {x_: _gapf(*x_, Fa) for x_ in ax}
            gb = {x_: _gapf(*x_, Fb) for x_ in ax}
            com = max(ax, key=lambda x_: min(ga[x_], gb[x_]))
            if min(ga[com], gb[com]) > -1e-13:
                use = [com]
            else:                                       # switch: the best line at each end, both at the node
                use = list({max(ax, key=lambda x_: ga[x_]), max(ax, key=lambda x_: gb[x_])})
            for o, t2, k in use:
                g, r = _axis_row(o, t2, k, F, nz)
                if g < 1e-3:
                    rows.append(r); rhs.append(g)
    m = len(rows)
    I, J, V = [], [], []
    for r_, row in enumerate(rows):
        for jj, v in row.items():
            I.append(r_); J.append(jj); V.append(v)
    Aub = coo_matrix(([-v for v in V], (I, J)), shape=(m, nz)).tocsr()   # g + grad·z >= 0  ->  -grad·z <= g
    # |z_j| <= w for every position / angle j:  z_j - w <= 0, -z_j - w <= 0
    nv = 3 * n
    eye_r = list(range(2 * nv)); eye_c = [j for j in range(nv)] * 2
    eye_v = [1.0] * nv + [-1.0] * nv
    W = coo_matrix((eye_v + [-1.0] * (2 * nv), (eye_r + eye_r, eye_c + [nz - 1] * (2 * nv))),
                   shape=(2 * nv, nz)).tocsr()
    from scipy.sparse import vstack as svstack
    obj = np.zeros(nz); obj[-1] = 1.0
    if free_side:
        obj[DS] = 1e3                                   # the side's increase first, then the move
    # scale by the largest violation: HiGHS' feasibility tolerance (1e-7 by default) would accept the overlaps as is
    rhs = np.array(rhs)
    tau = max(1e-15, -min(rhs.min(), 0.0))
    res = linprog(obj, A_ub=svstack([Aub, W]), b_ub=np.concatenate([rhs / tau, np.zeros(2 * nv)]),
                  bounds=[(None, None)] * nv + ([(0, None)] if free_side else []) + [(0, None)], method='highs',
                  options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
    if res.status != 0:
        return None
    z = res.x * tau
    return s + (z[DS] if free_side else 0.0), [(x + z[3 * i], y + z[3 * i + 1], t + z[3 * i + 2]) for i, (x, y, t) in enumerate(F)]


def midpoint(A, B, L):
    """the bisection node: the linear midpoint, put back on its contacts at fixed side (inc.correct: the step's
    curvature leaves overlaps O(Δ²), corrected to O(Δ⁴)) and lifted to an exact packing; the plain midpoint if that
    fails or lands above L.  Any node is valid (the step check covers its overlaps); this one keeps steps short of
    second-order loss."""
    import inc
    M = lerp(A, B, mpf(1) / 2)
    try:
        lp = midpoint_lp(A, B)
    except Exception:
        lp = None
    if lp is not None:
        nd = unwrap(M, lift(to_mp(lp)))
        if nd[0] <= L:
            return nd
    sF = float(M[0]); sqF = [(float(x), float(y), float(t)) for x, y, t in M[1]]
    R = max(max(abs(float(p - q)) for p, q in zip(u, v)) for u, v in zip(A[1], B[1])) or 1e-12
    try:
        c = inc.correct(sF, sqF, max(R, 1e-12))
    except Exception:
        c = None
    if c is None:
        return M
    nd = lift(to_mp((sF, c)))
    nd = unwrap(M, nd)
    return nd if nd[0] <= L else M


def certify(chain, L, max_depth=40, log=print):
    """certify every step of chain at level L, bisecting failures.  Returns dict(ok, steps, bisections, worst)"""
    stats = dict(steps=0, bisections=0, worst=mpf('inf'), fail=None)

    def rec(A, B, depth):
        ok, m, what = step_ok(A, B, L)
        stats['steps'] += 1
        if ok:
            stats['worst'] = min(stats['worst'], m)
            return True
        if depth >= max_depth:
            stats['fail'] = (float(m), what, depth)
            return False
        stats['bisections'] += 1
        M = midpoint(A, B, L)
        return rec(A, M, depth + 1) and rec(M, B, depth + 1)

    t0 = time.time()
    for idx in range(len(chain) - 1):
        if not rec(chain[idx], chain[idx + 1], 0):
            log(f'  step {idx} fails: {stats["fail"]}')
            return dict(ok=False, at=idx, **{k: (float(v) if isinstance(v, mpf) else v) for k, v in stats.items()},
                        sec=time.time() - t0)
    return dict(ok=True, **{k: (float(v) if isinstance(v, mpf) else v) for k, v in stats.items()}, sec=time.time() - t0)


def level(node):
    """the smallest L at which the node itself is a packing after scaling about the origin (a zero-length step):
    max(s, s/(1+g) over pairs (best axis gap g), s/(1+2m) over wall margins m), in closed form"""
    s, sq = node
    n = len(sq)
    L = s
    for i in range(n):
        e = _h(sq[i][2])
        for c in (0, 1):
            for m in (sq[i][c] - e, s - sq[i][c] - e):
                if m < 0:
                    L = max(L, s / (1 + 2 * m))
    ca = np.array([(float(x), float(y)) for x, y, _ in sq])
    for i in range(n):
        d2 = ((ca - ca[i]) ** 2).sum(1)
        for j in np.nonzero(d2 < 2.0 + 1e-6)[0]:
            j = int(j)
            if j <= i:
                continue
            g = max(_gap(o, t, k, sq) for o, t in ((i, j), (j, i)) for k in range(4))
            if g < 0:
                L = max(L, s / (1 + g))
    return L


def lift(node):
    """the node scaled about the origin to its level (plus 1e-28 relative): a packing in exact terms"""
    s, sq = node
    L = level(node) * (1 + mpf('1e-28'))
    f = L / s
    return L, [(x * f, y * f, t) for x, y, t in sq]


def regularize(chain, L):
    """lift every node after the first to an exact packing; drop nodes whose lifted side exceeds L"""
    out = [chain[0]]
    for node in chain[1:]:
        nd = lift(node)
        if nd[0] <= L:
            out.append(nd)
    return out


# ------------------------------------------------------------------------------------------------ chains
def descent_chain(P, log=print):
    """P's descent (exact/descend.py: MILP kicks + slp2) with every accepted iterate: [P, iterates ...]"""
    import slp2, descend as D
    sF = float(P[0]); sqF = [(float(x), float(y), float(t)) for x, y, t in P[1]]
    slp2.TRACE = []
    try:
        D.descend(sF, sqF, log=lambda *a: None)
        tr = slp2.TRACE
    finally:
        slp2.TRACE = None
    chain = [P] + [to_mp((s, sq)) for s, sq in tr]
    log(f'  descent: {len(tr)} iterates, {float(P[0]):.15f} -> {float(chain[-1][0]):.15f}')
    return chain


def unwrap_chain(chain):
    out = [chain[0]]
    for node in chain[1:]:
        out.append(unwrap(out[-1], node))
    return out


def center(sq, s, i):
    """the max-clearance pose of square i with the others fixed (Nelder-Mead on the smallest best-axis gap / wall
    margin, from its current pose): a waypoint for rattlers, which exact solutions park within ~1e-6 of a neighbour"""
    from scipy.optimize import minimize
    F = [(float(x), float(y), float(t)) for x, y, t in sq]
    sF = float(s)
    near = [j for j in range(len(F)) if j != i and (F[j][0] - F[i][0]) ** 2 + (F[j][1] - F[i][1]) ** 2 < 6.0]

    def clear(v):
        G = list(F); G[i] = tuple(v)
        e = 0.5 * (abs(math.cos(v[2])) + abs(math.sin(v[2])))
        c = min(v[0] - e, sF - v[0] - e, v[1] - e, sF - v[1] - e)
        for j in near:
            c = min(c, max(_gapf(o, t, k, G) for o, t in ((i, j), (j, i)) for k in range(4)))
        return -c

    r = minimize(clear, F[i], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
    return tuple(mpf(v) for v in r.x), -r.fun


def junction_orders(A, B, groups, rattlers=()):
    """candidate node sequences from A to B (same minimum, different points of its flat family / rattler positions):
    move one group at a time, in a few orders; then one square at a time, largest move first"""
    n = len(A[1])
    moved = [i for i in range(n) if max(abs(p - q) for p, q in zip(A[1][i], B[1][i])) > mpf('1e-25')]
    rest = [i for i in moved if not any(i in g for g in groups)]
    gs = [sorted(set(g) & set(moved)) for g in groups] + [rest]
    gs = [g for g in gs if g]
    orders = []
    # flat squares (and anything else) first, then each rattler via its cage's max-clearance pose
    rat = [i for i in rattlers if i in moved]
    if rat:
        others = [i for i in moved if i not in rat]
        nodes, cur = [A], list(A[1])
        for i in others:
            cur[i] = B[1][i]
        if others:
            nodes.append((max(A[0], B[0]), list(cur)))
        for i in rat:
            c, _ = center(cur, max(A[0], B[0]), i)
            c = (c[0], c[1], c[2] - mp.pi / 2 * mpmath.nint((c[2] - cur[i][2]) / (mp.pi / 2)))
            cur[i] = c
            nodes.append((max(A[0], B[0]), list(cur)))
            cur[i] = B[1][i]
            nodes.append((max(A[0], B[0]), list(cur)))
        nodes[-1] = B
        yield nodes
    orders += [gs, gs[::-1]]
    big = sorted(moved, key=lambda i: -max(abs(p - q) for p, q in zip(A[1][i], B[1][i])))
    orders.append([[i] for i in big])
    for order in orders:
        nodes, cur = [A], list(A[1])
        for g in order:
            for i in g:
                cur[i] = B[1][i]
            nodes.append((max(A[0], B[0]), list(cur)))
        nodes[-1] = B
        yield nodes


def junction(A, B, L, groups, rattlers=(), log=print):
    """certify a path from A to B through group-wise moves.  Returns ('full', nodes) when every square follows it;
    else ('rattlers', nodes) when it is certified for the squares that are rattlers in neither minimum (rattlers
    deleted: both ends are full packings, but a rattler may sit in a different cage); else (False, None)"""
    def attempt(A, B, groups, rattlers):
        for order in junction_orders(A, B, groups, rattlers):
            nodes = [order[0]] + [lift(nd) for nd in order[1:-1]] + [order[-1]]
            if any(nd[0] > L for nd in nodes):
                continue
            if certify(nodes, L, log=lambda *a: None)['ok']:
                return nodes
        return None
    J = attempt(A, B, groups, rattlers)
    if J:
        return 'full', J
    keep = [i for i in range(len(A[1])) if i not in set(rattlers)]
    if len(keep) < len(A[1]):
        re = {i: k for k, i in enumerate(keep)}
        sub = lambda nd: (nd[0], [nd[1][i] for i in keep])
        g2 = [[re[i] for i in g if i in re] for g in groups]
        J = attempt(sub(A), sub(B), g2, ())
        if J:
            return 'rattlers', J
    return False, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('posting'); ap.add_argument('--exact', help="P's minimum (exactsolve .exact.txt); default: no exact end")
    ap.add_argument('--no-descent', action='store_true', help='step straight from P to its minimum')
    ap.add_argument('--to', help='a second posting Q (lower side)'); ap.add_argument('--exact-to', help="Q's minimum")
    ap.add_argument('--out')
    a = ap.parse_args()
    P = read_node(a.posting)
    L = level(P) * (1 + mpf('1e-25'))       # touching contacts have margin exactly 0 at P; 1e-25 covers the 1e-30 rounding margin
    print(f'P: n={len(P[1])}, side {mpmath.nstr(P[0], 20)}, level L = side + {mpmath.nstr(L - P[0], 3)}')
    chain = [P] if a.no_descent else descent_chain(P)
    if a.exact:
        M = read_node(a.exact)
        k, perm, cost = match(chain[-1], M)
        permP = perm
        chain.append(apply(M, k, perm))
        print(f'  + exact minimum {mpmath.nstr(M[0], 20)} (match cost {cost:.1e})')
    if a.to:
        Q = read_node(a.to)
        chainQ = [Q] if a.no_descent else descent_chain(Q)
        if a.exact_to:
            chainQ.append(read_node(a.exact_to))
        k, perm, cost = match(chain[-1], chainQ[-1])
        print(f"  Q's chain matched onto P's end (symmetry {k}, cost {cost:.1e})")
        back = unwrap_chain([chain[-1]] + [apply(node, k, perm) for node in reversed(chainQ)])[1:]
        chain = regularize(unwrap_chain(chain), L)
        back = unwrap_chain([chain[-1]] + back)[1:]                  # quarter turns against the final endpoint
        groups = []
        for f in (a.exact, a.exact_to):
            j = f and f.replace('.exact.txt', '.json')
            if j and os.path.exists(j):
                rep = json.load(open(j))
                fr = rep.get('free_squares'); fr = json.loads(fr) if isinstance(fr, str) else (fr or [])
                fl = (rep.get('second_order') or {}).get('flat_squares') or []
                pm = perm if f == a.exact_to else permP              # the report's labels -> the chain's
                inv = {p: i for i, p in enumerate(pm)}
                fr = [inv[x] for x in fr]; fl = [inv[x] for x in fl]
                groups += [fr, fl]
        ok, J = junction(chain[-1], back[0], L, groups, rattlers=sorted(set(groups[0]) | set(groups[2] if len(groups) > 2 else [])))
        print(f"  junction (P's minimum -> Q's, within the flat family): {(ok + ', ' + str(len(J) - 2) + ' intermediate nodes') if ok else 'NOT certified'}")
        r0 = certify(chain, L)
        rQ = certify(regularize(back, L), L) if len(back) > 1 else dict(ok=True, steps=0, bisections=0)
        r = dict(ok=bool(r0['ok'] and ok == 'full' and rQ['ok']), ok_mod_rattlers=bool(r0['ok'] and ok and rQ['ok']),
                 P_chain=r0, junction=ok, Q_chain=rQ)
    else:
        chain = regularize(unwrap_chain(chain), L)
        r = certify(chain, L)
    r.update(posting=a.posting, to=a.to, exact=a.exact, exact_to=a.exact_to, L=mpmath.nstr(L, 30),
             S_min=mpmath.nstr(chain[-1][0], 30))
    print(json.dumps(r, indent=None, default=str))
    if a.out:
        json.dump(r, open(a.out, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
