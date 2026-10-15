#!/usr/bin/env python3
"""Fast upright fill (FP1, 10-15): the same problem as gen.mis_fill2 (max set of non-overlapping axis-aligned unit
squares among gen.candidates2's candidates around fixed tilted squares), solved by safe reductions + a small kernel.

Corner rule (safe): a candidate v whose conflicting live neighbours all lie weakly right of and above it (or one of the
3 mirror images) has a clique neighbourhood (every such neighbour contains the point just inside v's far corner), so
v is simplicial and some maximum independent set contains it: take v, delete its neighbours, repeat (worklist).
Isolated candidates are the degenerate case.  What survives (the kernel: lattice seams, staircase bands at the tilted
squares) splits into connected components, each solved exactly (MILP) or by min-degree greedy.

  fill(s, T, mode='exact' | 'greedy', tau=0.0) -> np.array of T + the chosen uprights (x, y, deg)
  stats of the last call: fillfast.LAST (candidates, kernel size, components, largest component, times)
tau > 0: tolerant fill, candidates may overlap the tilted squares and each other by up to tau (overcounts; for starts
that go straight to an optimiser only).
"""
import math, time
from collections import deque
import numpy as np

LAST = {}
EPS = 1e-12


def candidates_slow(s, T, tau=0.0):
    import gen
    return gen.candidates2(s, [tuple(q) for q in T], tol=max(1e-6, tau))


def _strip_extent(P, lo, hi, axis):
    """vectorised gen._clip_extent: P (k, 4, 2) convex quads, strips [lo, hi] (m,) in coordinate `axis`; -> (min, max) of
    the other coordinate over quad ∩ strip, shape (k, m) each (nan where empty)"""
    a = P[:, :, axis][:, None, :]                       # (k, 1, 4)
    o = P[:, :, 1 - axis][:, None, :]
    lo = lo[None, :, None]; hi = hi[None, :, None]
    inside = (a >= lo - 1e-9) & (a <= hi + 1e-9)         # vertices inside the strip (inclusive: exact touches count)
    big = 1e18
    mn = np.where(inside, o, big).min(-1); mx = np.where(inside, o, -big).max(-1)
    a2 = np.roll(a, -1, axis=-1); o2 = np.roll(o, -1, axis=-1)
    for b in (lo, hi):                                   # edge crossings of the strip's two boundary lines
        cross = (a - b) * (a2 - b) < 0
        t = np.where(cross, (b - a) / np.where(a2 != a, a2 - a, 1.0), 0.0)
        v = o + t * (o2 - o)
        mn = np.minimum(mn, np.where(cross, v, big).min(-1)); mx = np.maximum(mx, np.where(cross, v, -big).max(-1))
    empty = mn > mx
    return np.where(empty, np.nan, mn), np.where(empty, np.nan, mx)


def candidates(s, T, tau=0.0, reach=4):
    """vectorised gen.candidates2 (same candidate set): corner-lattice points + staircase points on every wall-anchored
    row / column touching a tilted square from either side (up to `reach` steps), minus those overlapping a tilted
    square by more than tol (separating-axis test, as gen.pen)"""
    tol = max(1e-6, tau)
    T = np.asarray(T, float).reshape(-1, 3)
    m = int(math.floor(s))
    lat = np.r_[0.5 + np.arange(m), s - 0.5 - np.arange(m)]
    gx, gy = np.meshgrid(lat, lat, indexing='ij')
    pts = [np.c_[gx.ravel(), gy.ravel()]]
    if len(T):
        th = np.radians(T[:, 2]); c, si = np.cos(th), np.sin(th)
        off = np.array([[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]])
        P = T[:, None, :2] + np.stack([off[:, 0] * c[:, None] - off[:, 1] * si[:, None],
                                       off[:, 0] * si[:, None] + off[:, 1] * c[:, None]], -1)
        steps = np.arange(reach)
        for axis in (1, 0):                              # axis 1: rows (fixed y), axis 0: columns (fixed x)
            mn, mx = _strip_extent(P, lat - 0.5, lat + 0.5, axis)
            for start, d in ((mn - 0.5, -1.0), (mx + 0.5, 1.0)):
                u = start[..., None] + d * steps         # (k, m2, reach)
                cc = np.broadcast_to(lat[None, :, None], u.shape)
                ok = ~np.isnan(u) & (u >= 0.5 - 1e-9) & (u <= s - 0.5 + 1e-9)
                uu = np.clip(u[ok], 0.5, s - 0.5); c2 = cc[ok]
                pts.append(np.c_[uu, c2] if axis == 1 else np.c_[c2, uu])
    Q = np.unique(np.round(np.concatenate(pts), 12), axis=0)
    if len(T) and len(Q):
        keep = np.ones(len(Q), bool)
        h = 0.5 * (np.abs(np.cos(th)) + np.abs(np.sin(th)))
        for k0 in range(0, len(Q), 2048):
            q = Q[k0:k0 + 2048]
            dx = T[None, :, 0] - q[:, None, 0]; dy = T[None, :, 1] - q[:, None, 1]
            near = dx * dx + dy * dy < 2.0
            r = 0.5 + h[None, :]
            ox = r - np.abs(dx); oy = r - np.abs(dy)
            ou = r - np.abs(dx * c[None, :] + dy * si[None, :]); ov = r - np.abs(-dx * si[None, :] + dy * c[None, :])
            pen = np.minimum(np.minimum(ox, oy), np.minimum(ou, ov))
            keep[k0:k0 + 2048] = ~((pen > tol) & near).any(1)
        Q = Q[keep]
    return [(float(x), float(y), 0.0) for x, y in Q]


def conflict_graph(C, axtol=1e-6):
    """CSR-style adjacency via one vectorised pairwise test (row chunks)"""
    N = len(C)
    X = np.array([c[0] for c in C]); Y = np.array([c[1] for c in C])
    adj = [None] * N
    for k0 in range(0, N, 1024):
        B = (np.abs(X[k0:k0 + 1024, None] - X[None, :]) < 1 - axtol) & (np.abs(Y[k0:k0 + 1024, None] - Y[None, :]) < 1 - axtol)
        for i, row in enumerate(B):
            row[k0 + i] = False
            adj[k0 + i] = np.nonzero(row)[0]
    return X, Y, adj


def conflict_graph_slow(C, axtol=1e-6):
    """adjacency (list of int arrays): axis squares conflict iff |dx| < 1 - axtol and |dy| < 1 - axtol"""
    N = len(C)
    X = np.array([c[0] for c in C]); Y = np.array([c[1] for c in C])
    buckets = {}
    for k in range(N):
        buckets.setdefault((int(math.floor(X[k])), int(math.floor(Y[k]))), []).append(k)
    adj = [[] for _ in range(N)]
    for k in range(N):
        bx, by = int(math.floor(X[k])), int(math.floor(Y[k]))
        for ux in (bx - 1, bx, bx + 1):
            for uy in (by - 1, by, by + 1):
                for l in buckets.get((ux, uy), ()):
                    if l > k and abs(X[l] - X[k]) < 1 - axtol and abs(Y[l] - Y[k]) < 1 - axtol:
                        adj[k].append(l); adj[l].append(k)
    return X, Y, [np.array(a, dtype=np.int64) for a in adj]


def corner_reduce(X, Y, adj):
    """counter version: blk[v, o] = live neighbours that block v in orientation o (LL / LR / UL / UR); v is takeable
    when some blk[v, o] == 0; deaths decrement their neighbours' counters"""
    N = len(X)
    alive = np.ones(N, bool); taken = []
    if N == 0:
        return taken, alive
    src = np.repeat(np.arange(N), [len(a) for a in adj]); dst = np.concatenate(adj) if N else np.zeros(0, int)
    dx = X[dst] - X[src]; dy = Y[dst] - Y[src]
    # u blocks v in orientation LL if u is strictly left or strictly below v, etc.
    M = np.stack([(dx < -EPS) | (dy < -EPS), (dx > EPS) | (dy < -EPS), (dx < -EPS) | (dy > EPS), (dx > EPS) | (dy > EPS)], 1)
    blk = np.zeros((N, 4), np.int64)
    np.add.at(blk, src, M.astype(np.int64))
    start = np.r_[0, np.cumsum([len(a) for a in adj])]
    Ml = M.tolist(); dstl = dst.tolist(); blkl = blk.tolist()
    work = deque(v for v in range(N) if min(blkl[v]) == 0); queued = np.zeros(N, bool); queued[list(work)] = True
    alivel = [True] * N
    while work:
        v = work.popleft(); queued[v] = False
        if not alivel[v] or min(blkl[v]) != 0:
            continue
        taken.append(v); dead = [v] + [u for u in adj[v].tolist() if alivel[u]]
        for d in dead:
            alivel[d] = False
        for d in dead:                                    # d no longer blocks its live neighbours
            for e in range(start[d], start[d + 1]):
                w = dstl[e]
                if not alivel[w]:
                    continue
                # the edge stored at d -> w says how w blocks d; we need how d blocks w: mirror the masks
                md = Ml[e]                                # blocks of w as seen from d (w rel. d)
                bw = blkl[w]
                # d rel. w is the negation: LL<->UR, LR<->UL
                bw[3] -= md[0]; bw[2] -= md[1]; bw[1] -= md[2]; bw[0] -= md[3]
                if min(bw) == 0 and not queued[w]:
                    work.append(w); queued[w] = True
    alive[:] = alivel
    return taken, alive


def corner_reduce_slow(X, Y, adj):
    """-> (taken list, alive bool array): repeated corner-rule takes (safe), worklist-driven"""
    N = len(X)
    alive = np.ones(N, bool)
    taken = []
    work = deque(range(N)); queued = np.ones(N, bool)
    while work:
        v = work.popleft(); queued[v] = False
        if not alive[v]:
            continue
        nb = adj[v][alive[adj[v]]]
        if len(nb):
            dx = X[nb] - X[v]; dy = Y[nb] - Y[v]
            ok = ((dx >= -EPS) & (dy >= -EPS)).all() or ((dx <= EPS) & (dy >= -EPS)).all() or \
                 ((dx >= -EPS) & (dy <= EPS)).all() or ((dx <= EPS) & (dy <= EPS)).all()
            if not ok:
                continue
        taken.append(v); alive[v] = False; alive[nb] = False
        for u in nb:                                  # their live neighbours may have become corners
            for w in adj[u]:
                if alive[w] and not queued[w]:
                    work.append(w); queued[w] = True
    return taken, alive


def components(alive, adj):
    seen = np.zeros(len(alive), bool); comps = []
    for v in np.nonzero(alive)[0]:
        if seen[v]:
            continue
        comp, dq = [v], deque([v]); seen[v] = True
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if alive[w] and not seen[w]:
                    seen[w] = True; comp.append(w); dq.append(w)
        comps.append(comp)
    return comps


def _milp(comp, adj, time_limit=10):
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    idx = {v: i for i, v in enumerate(comp)}
    rows, cols, r = [], [], 0
    for v in comp:
        for w in adj[v]:
            if w > v and w in idx:
                rows += [r, r]; cols += [idx[v], idx[w]]; r += 1
    if r == 0:
        return list(comp)
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(r, len(comp))).tocsr()
    res = milp(c=-np.ones(len(comp)), constraints=LinearConstraint(A, -np.inf, 1), integrality=np.ones(len(comp)),
               bounds=Bounds(0, 1), options={'time_limit': time_limit})
    x = res.x if res.x is not None else np.zeros(len(comp))
    return [comp[i] for i in range(len(comp)) if x[i] > 0.5]


def _greedy(comp, adj):
    """min-degree greedy inside one component (degrees updated as vertices are removed)"""
    live = set(comp)
    deg = {v: sum(1 for w in adj[v] if w in live) for v in comp}
    out = []
    while live:
        v = min(live, key=lambda u: (deg[u], u))
        out.append(v)
        gone = [v] + [w for w in adj[v] if w in live]
        for g in gone:
            live.discard(g)
        for g in gone:
            for w in adj[g]:
                if w in live:
                    deg[w] -= 1
    return out


def fill(s, T, mode='exact', tau=0.0, time_limit=10):
    t0 = time.process_time()
    T = np.asarray(T, float).reshape(-1, 3)
    C = candidates(s, T, tau)
    t1 = time.process_time()
    X, Y, adj = conflict_graph(C, axtol=max(1e-6, tau))
    taken, alive = corner_reduce(X, Y, adj)
    t2 = time.process_time()
    comps = components(alive, adj)
    chosen = list(taken)
    for comp in comps:
        if len(comp) == 1:
            chosen.append(comp[0])
        else:
            chosen += _milp(comp, adj, time_limit) if mode == 'exact' else _greedy(comp, adj)
    t3 = time.process_time()
    LAST.clear()
    LAST.update(cand=len(C), taken=len(taken), kernel=int(alive.sum()), comps=len(comps),
                maxcomp=max((len(c) for c in comps), default=0), t_cand=t1 - t0, t_reduce=t2 - t1, t_kernel=t3 - t2)
    U = np.array([C[i] for i in chosen], float).reshape(-1, 3)
    return np.r_[T, U]
