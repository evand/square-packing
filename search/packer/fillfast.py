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


def candidates(s, T, tau=0.0):
    import gen
    return gen.candidates2(s, [tuple(q) for q in T], tol=max(1e-6, tau))


def conflict_graph(C, axtol=1e-6):
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
