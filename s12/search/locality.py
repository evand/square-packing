#!/usr/bin/env python3
"""locality.py -- LOC_D, the WINDOW-CONSISTENCY relaxation of a level-2 leaf on a fixed pose set.

    LOC_D(P) = max  sum mu
               s.t. cov(x) <= 1                 at every arrangement vertex (as `cliquelever`)
                    mu(R)  =  k_R               the leaf's pinned regions (boundary poses free)
                    mu(wall strip) <= 3         chord rows
                    mu|_{P_j} in Pi_j           for every window j

where `P_j = {S in P : centre(S) in W_j}` for a family of windows `W_j` (closed boxes of side `D`
or closed discs of diameter `D`, centres on a pitch `D/4`, clipped to the container) and
`Pi_j = conv{ 1_I : I independent in G[P_j] }` is the packing (stable-set) polytope of the exact
closed-intersection graph restricted to the window.  `tasks/locality-ceiling/README.md`.

`mu|_{P_j} in Pi_j` implies EVERY valid inequality of every local family supported in `W_j`, so
`LOC_D` is the ceiling of local reasoning at window diameter `D`:  `LOC_D >= alpha(P) = 11` always,
`LOC_D <= QSTAB(P)` as soon as `D - pitch >= sqrt(2)` (every clique of `G` has centre-extent
`<= sqrt 2`), and `LOC_4 = alpha(P)`.

Dantzig-Wolfe.  Variables `lam_{j,I} >= 0` over the independent sets `I` of `G[P_j]`:

    conv_j:      sum_I lam_{j,I} <= 1                                (one row per window)
    link_{j,i}:  mu_i <= sum_{I ni i} lam_{j,I}                      (one row per (j, S))

(the linking rows are inequalities because the stable-set polytope is down-monotone: subsets of
independent sets are independent).  The loop runs this master in its equivalent CUTTING-PLANE form,
which is the same LP read through its dual and converges from ABOVE instead of from below:

  * the outer LP is in the `mu` variables alone -- coverage rows, the region pins, the chord rows
    and the window cuts found so far -- so its value is an UPPER bound on `LOC_D(P)` at every
    iteration, every row being valid;
  * the window subproblem, for the current `mu`, is the Dantzig-Wolfe restricted master of window
    `j` restricted to the SUPPORT of `mu` (a pose with `mu_i = 0` constrains nothing, and an
    independent set of `G[P_j]` meets the support in an independent set of the support subgraph):

        z_j = min sum_I lam_I   s.t.   sum_{I ni i} lam_I >= mu_i,   lam >= 0

    over ALL independent sets of the support subgraph, which are enumerated exactly (Bron-Kerbosch
    on the complement: a window support is 5-60 poses).  `mu|_{P_j} in Pi_j` iff `z_j <= 1`, and
    the optimal dual `pi >= 0` has `sum_{i in I} pi_i <= 1` for every independent set, so
    `sum_i pi_i mu_i <= 1` is a valid inequality for the packing polytope violated by exactly
    `z_j - 1`.  That cut is then EXTENDED over the whole window -- a pose outside the support may
    carry `1 - (heaviest independent set of the support it can join)`, by exact max-weight
    independent set (`leaf_ceiling.max_clique_int` on the complement), and the extended vector is
    renormalised by its own exact max-weight independent set so that it is valid whatever happened;
  * convergence: no arrangement vertex violated, no window with `z_j > 1`.  The LP value is then
    `LOC_D(P)` and the solution is LOC_D-feasible.

Directions (as `T4SCREEN.md` 1.2 / `CLIQUELEVER.md` 1).  Restricting the poses LOWERS the value, so
every number here is a lower bound on the continuum `LOC_D`.  On the pose set itself every
intermediate LP value is an UPPER bound on `LOC_D(P)` (rows only relax), and at convergence it is
the value.  The exactly certified final measure -- round down to `/10^9`, top up the pinned regions
on exact window slack, then verify `M <= 1`, the region counts and a RATIONAL independent-set
decomposition of every window (exact simplex over the enumerated maximal independent sets, every
set re-verified pairwise disjoint by the exact integer SAT) -- is a rigorous LOWER bound on
`LOC_D(P)` whatever the loop did.

    python3 search/locality.py run TAG --poses runs/cl_A0101L2_poses.txt \\
        --start runs/cl_A0101L2_exact.txt --corners 1111 --patterns 01010101 --chord \\
        --D 2 --shape box
    python3 search/locality.py selftest
"""
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import leaf_ceiling as lc                                            # noqa: E402
from cliquelever import Poses, LABS, FTOL                            # noqa: E402

RUNS = os.path.join(REPO, 'runs')
DM = 10 ** 9                     # denominator of the exact final measure


# ====================================================================== exact closed graph
def build_adj(ps, log=None):
    """dense boolean closed-intersection matrix over the poses (exact: float SAT with the exact
    integer test on the |margin| <= FTOL band, centre-distance prefilter)"""
    t0 = time.time()
    n = ps.n
    A = np.zeros((n, n), dtype=bool)
    F = ps.F
    nex = 0
    for i in range(n):
        d2 = (F[:, 0] - F[i, 0]) ** 2 + (F[:, 1] - F[i, 1]) ** 2
        J = np.nonzero(d2 <= 2.0 + 1e-9)[0]
        m = ps.meet_margin(i, J)
        ok = m > FTOL
        for k in np.nonzero(np.abs(m) <= FTOL)[0]:
            nex += 1
            ok[k] = lc.sq_meets_sq(ps.sq[i], ps.sq[int(J[k])])
        A[i, J[ok]] = True
    np.fill_diagonal(A, False)
    assert (A == A.T).all(), 'closed-intersection matrix is not symmetric'
    if log:
        log(f'   closed-intersection graph: {n} poses, {int(A.sum()) // 2} edges '
            f'({nex} exact border tests, {time.time() - t0:.1f}s)')
    return A


# ====================================================================== windows
def make_windows(ps, D, shape, pitch_frac, log=None):
    """the window family: closed boxes of side D (or closed discs of diameter D) with centres on a
    pitch `pitch_frac * D` grid, clipped to the container, reduced to the inclusion-MAXIMAL pose
    sets (two windows with the same poses impose the same constraint; a window whose poses are
    contained in another's imposes a weaker one)."""
    t = float(ps.t)
    p = pitch_frac * D
    ks = list(np.arange(0.0, t + 1e-9, p))
    if ks[-1] < t - 1e-9:
        ks.append(t)
    C = ps.F[:, :2]
    cand = []
    for cx in ks:
        for cy in ks:
            if shape == 'box':
                x0, x1 = max(0.0, cx - D / 2), min(t, cx + D / 2)
                y0, y1 = max(0.0, cy - D / 2), min(t, cy + D / 2)
                sel = np.nonzero((C[:, 0] >= x0 - 1e-12) & (C[:, 0] <= x1 + 1e-12) &
                                 (C[:, 1] >= y0 - 1e-12) & (C[:, 1] <= y1 + 1e-12))[0]
                geo = ('box', round(x0, 9), round(x1, 9), round(y0, 9), round(y1, 9))
            else:
                d2 = (C[:, 0] - cx) ** 2 + (C[:, 1] - cy) ** 2
                sel = np.nonzero(d2 <= (D / 2) ** 2 + 1e-12)[0]
                geo = ('disc', round(cx, 9), round(cy, 9), round(D / 2, 9))
            if len(sel):
                cand.append((geo, frozenset(int(v) for v in sel)))
    # inclusion-maximal, deduplicated
    seen = {}
    for geo, S in cand:
        seen.setdefault(S, geo)
    sets = sorted(seen, key=lambda S: -len(S))
    keep = []
    for S in sets:
        if not any(S <= T for T in keep):
            keep.append(S)
    wins = [(seen[S], np.array(sorted(S), dtype=np.int64)) for S in keep]
    wins.sort(key=lambda w: (w[0][1], w[0][2]))
    if log:
        log(f'   windows: shape={shape} D={D} pitch={p:.4f} -> {len(cand)} candidates, '
            f'{len(wins)} inclusion-maximal; sizes {min(len(w[1]) for w in wins)}..'
            f'{max(len(w[1]) for w in wins)}, {sum(len(w[1]) for w in wins)} linking rows')
    return wins


# ====================================================================== independent sets
def greedy_is(sub, order, alive=None):
    """a maximal independent set of the dense-bool subgraph, taking vertices in `order`"""
    n = sub.shape[0]
    alive = np.ones(n, bool) if alive is None else alive.copy()
    out = []
    for v in order:
        v = int(v)
        if alive[v]:
            out.append(v)
            alive &= ~sub[v]
            alive[v] = False
            if not alive.any():
                break
    return out


def comp_nbr(sub):
    """complement adjacency as lists of sets (independent sets of G = cliques of the complement)"""
    comp = ~sub
    np.fill_diagonal(comp, False)
    return [set(np.nonzero(comp[v])[0].tolist()) for v in range(comp.shape[0])]


def mwis_exact(sub, w_int, time_limit, lb):
    """maximum-weight independent set of the subgraph `sub` (dense bool, LOCAL indices) with
    integer weights -> (weight, members, complete).  Independent sets of G are cliques of the
    complement, so this is `leaf_ceiling.max_clique_int` on the complement (which is sparse inside
    a window: two squares are non-adjacent only if their centres are more than ~1 apart)."""
    return lc.max_clique_int(comp_nbr(sub), list(w_int), time_limit, lb=lb)


# ====================================================================== exact rational simplex
def exact_cover_lp(cols, mu, nrow, cap=None):
    """EXACT rational LP:   min sum lam   s.t.  sum_{I ni i} lam_I >= mu_i (i < nrow),  lam >= 0,
    with `cols` a list of frozensets of row indices and `mu` a list of Fractions.
    Returns (value, {col index: Fraction}).  Tableau simplex in Fractions, Bland's rule (no
    cycling, no tolerance).  `cap`: stop early once the phase-2 objective is <= cap."""
    m = nrow
    ncol = len(cols)
    # variables: lam (ncol) | surplus s (m) | artificial a (m)
    N = ncol + 2 * m
    Z = Fr(0)
    T = [[Z] * (N + 1) for _ in range(m)]
    for j, I in enumerate(cols):
        for i in I:
            T[i][j] = Fr(1)
    for i in range(m):
        T[i][ncol + i] = Fr(-1)
        T[i][ncol + m + i] = Fr(1)
        T[i][N] = Fr(mu[i])
    basis = [ncol + m + i for i in range(m)]

    def run(cost, stop=None):
        # reduced-cost row  z_j - c_j  with the current basis
        while True:
            y = [cost[basis[i]] for i in range(m)]
            obj = sum(y[i] * T[i][N] for i in range(m))
            if stop is not None and obj <= stop:
                return obj
            enter = -1
            for j in range(N):
                if j in basis:
                    continue
                zj = sum(y[i] * T[i][j] for i in range(m))
                if zj - cost[j] > 0:
                    enter = j
                    break                       # Bland: smallest index
            if enter < 0:
                return obj
            ratio = None
            leave = -1
            for i in range(m):
                if T[i][enter] > 0:
                    r = T[i][N] / T[i][enter]
                    if ratio is None or r < ratio or (r == ratio and basis[i] < basis[leave]):
                        ratio, leave = r, i
            if leave < 0:
                return None                     # unbounded (cannot happen here)
            piv = T[leave][enter]
            T[leave] = [v / piv for v in T[leave]]
            for i in range(m):
                if i != leave and T[i][enter] != 0:
                    f = T[i][enter]
                    Ti, Tl = T[i], T[leave]
                    T[i] = [Ti[k] - f * Tl[k] for k in range(N + 1)]
            basis[leave] = enter

    c1 = [Z] * N
    for i in range(m):
        c1[ncol + m + i] = Fr(1)
    v1 = run(c1)
    if v1 is None or v1 > 0:
        return None, None                       # infeasible (impossible: mu >= 0)
    # drive artificials out of the basis where possible, then forbid them
    for i in range(m):
        if basis[i] >= ncol + m:
            for j in range(ncol + m):
                if T[i][j] != 0:
                    piv = T[i][j]
                    T[i] = [v / piv for v in T[i]]
                    for k in range(m):
                        if k != i and T[k][j] != 0:
                            f = T[k][j]
                            Tk, Ti = T[k], T[i]
                            T[k] = [Tk[q] - f * Ti[q] for q in range(N + 1)]
                    basis[i] = j
                    break
    c2 = [Z] * N
    for j in range(ncol):
        c2[j] = Fr(1)
    for i in range(m):
        c2[ncol + m + i] = Fr(10 ** 6)          # keep artificials at zero
    v2 = run(c2, stop=cap)
    sol = {}
    for i in range(m):
        if basis[i] < ncol and T[i][N] != 0:
            sol[basis[i]] = T[i][N]
    return v2, sol


def maximal_is_enum(sub, cap=200000):
    """every MAXIMAL independent set of the dense-bool subgraph `sub` (Bron-Kerbosch with a pivot
    on the complement).  Returns (list of frozensets, complete)."""
    n = sub.shape[0]
    comp = ~sub
    np.fill_diagonal(comp, False)
    nb = [set(np.nonzero(comp[v])[0].tolist()) for v in range(n)]
    out = []
    full = [True]

    def bk(R, Pv, X):
        if len(out) >= cap:
            full[0] = False
            return
        if not Pv and not X:
            out.append(frozenset(R))
            return
        u = max(Pv | X, key=lambda v: len(nb[v] & Pv))
        for v in list(Pv - nb[u]):
            bk(R | {v}, Pv & nb[v], X & nb[v])
            Pv = Pv - {v}
            X = X | {v}
            if not full[0]:
                return

    sys.setrecursionlimit(10000)
    bk(set(), set(range(n)), set())
    return out, full[0]


# ====================================================================== the LP
class Model:
    """a HiGHS model with two column blocks (the pose columns x, then the lam columns) and four
    row blocks (region equalities, chord rows, window convexity, window linking) plus the coverage
    rows appended by separation."""

    def __init__(self, threads, log, crossover=False):
        self.crossover = bool(crossover)
        import highspy
        self.hp = highspy
        self.INF = highspy.kHighsInf
        self.threads = threads
        self.log = log
        self.h = None
        self.ncol = 0

    def start(self, ncol):
        h = self.hp.Highs()
        h.setOptionValue('output_flag', False)
        h.setOptionValue('threads', self.threads)
        h.setOptionValue('presolve', 'on')
        h.addVars(ncol, np.zeros(ncol), np.full(ncol, self.INF))
        h.changeColsCost(ncol, np.arange(ncol, dtype=np.int32), -np.ones(ncol))
        self.h, self.ncol = h, ncol
        self.nrow = 0
        self.basis = False

    def add_rows(self, lo, hi, M):
        M = sp.csr_matrix(M)
        assert M.shape[1] == self.ncol, (M.shape, self.ncol)
        self.h.addRows(M.shape[0], np.asarray(lo, float), np.asarray(hi, float), M.nnz,
                       M.indptr.astype(np.int32), M.indices.astype(np.int32),
                       M.data.astype(float))
        self.nrow += M.shape[0]

    def add_cols(self, cost, M):
        """M: csc matrix (nrow x k) of the new columns"""
        M = sp.csc_matrix(M)
        assert M.shape[0] == self.nrow, (M.shape, self.nrow)
        k = M.shape[1]
        self.h.addCols(k, np.asarray(cost, float), np.zeros(k), np.full(k, self.INF), M.nnz,
                       M.indptr.astype(np.int32), M.indices.astype(np.int32),
                       M.data.astype(float))
        self.ncol += k

    def run(self, tlim=0.0):
        """ipm + crossover on the first solve (a fresh model has no basis, and cold dual simplex on
        a 26k x 10k coverage LP is far slower); afterwards a warm dual simplex capped at `tlim`
        seconds, falling back to ipm -- `t4leaf.Hi.run`'s rule.  `tlim <= 0` always uses ipm."""
        h = self.h
        warm = self.basis and tlim > 0
        h.setOptionValue('solver', 'simplex' if warm else 'ipm')
        # The cutting-plane loop never reads the outer LP's duals -- separation works on `mu`
        # alone -- so the crossover is pure cost on a 6k x 10k coverage LP with a large degenerate
        # optimal face.  Off by default; `--crossover` restores it.
        h.setOptionValue('run_crossover', 'on' if self.crossover else 'off')
        h.setOptionValue('time_limit', tlim if warm else 1e30)
        h.run()
        ms = h.getModelStatus()
        if ms != self.hp.HighsModelStatus.kOptimal and warm:
            self.log(f'   [HiGHS simplex -> {ms} after {tlim:.0f}s; ipm+crossover]')
            h.setOptionValue('solver', 'ipm')
            h.setOptionValue('time_limit', 1e30)
            h.run()
            ms = h.getModelStatus()
        if ms != self.hp.HighsModelStatus.kOptimal:
            self.log(f'   [HiGHS -> {ms}]')
            return None
        self.basis = True
        sol = h.getSolution()
        info = h.getInfo()
        self.iters = int(info.simplex_iteration_count)
        return (np.maximum(np.asarray(sol.col_value, float), 0.0),
                np.asarray(sol.row_dual, float), -float(info.objective_function_value))



# ====================================================================== window separation
def cover_lp_cg(sub, mu, tol=1e-9, max_iter=60, time_limit=3.0, cols=None, budget=30.0):
    """the window subproblem, by column generation on the SUPPORT.

        z = min sum lam   s.t.  sum_{I ni i} lam_I >= mu_i,  lam >= 0,  I independent in G[S]

    (`S` = the poses of the window carrying mass; a pose with `mu_i = 0` constrains nothing, and an
    independent set of `G[P_j]` meets `S` in an independent set of `G[S]`, so nothing is lost).
    `mu|_{P_j} in Pi_j` iff `z <= 1`, and `z` IS the Dantzig-Wolfe convexity multiplier the window
    needs.  The optimal dual `pi >= 0` satisfies `sum_{i in I} pi_i <= 1` for every independent set
    -- which is exactly what the pricing max-weight independent set proves on the iteration that
    ends the generation -- so `sum_i pi_i mu_i <= 1` is a VALID inequality for the packing polytope,
    violated by `z - 1`.  Returns (z, pi, columns, complete)."""
    from scipy.optimize import linprog
    n = len(mu)
    mu = np.asarray(mu, float)
    if cols is None:
        cols = []
    cols = list(cols)
    res = mu.copy()
    for _ in range(min(4 * n + 8, 400)):            # greedy peel: a feasible starting basis
        if res.max() <= 1e-12:
            break
        I = greedy_is(sub, np.argsort(-res))
        cols.append(frozenset(I))
        v = min(res[k] for k in I)
        for k in I:
            res[k] -= v
    seen = set(cols)
    for i in range(n):                              # every pose in at least one column
        if not any(i in I for I in cols):
            I = frozenset(greedy_is(sub, [i] + list(np.argsort(-mu))))
            if I not in seen:
                cols.append(I)
                seen.add(I)
    z, pi, comp = None, None, True
    t0 = time.time()
    for _ in range(max_iter):
        if time.time() - t0 > budget:
            comp = False
            break
        A = np.zeros((n, len(cols)))
        for k, I in enumerate(cols):
            for i in I:
                A[i, k] = 1.0
        r = linprog(np.ones(len(cols)), A_ub=-A, b_ub=-mu, bounds=(0, None), method='highs')
        if not r.success:
            return None, None, cols, False
        z = float(r.fun)
        pi = np.maximum(-r.ineqlin.marginals, 0.0)
        w = np.round(pi * 10 ** 9).astype(np.int64)
        v, mem, c = mwis_exact(sub, list(w), time_limit, lb=10 ** 9)
        comp = bool(c)
        if not c or v <= 10 ** 9 * (1 + tol) or not mem:
            comp = bool(c)
            break
        I = frozenset(int(k) for k in mem)
        if I in seen:
            break
        cols.append(I)
        seen.add(I)
    return z, pi, cols, comp


def extend_cut(sub_win, S_local, pi_S, a):
    """extend a cut `pi` supported on `S` (local indices into the window) over the WHOLE window.

    A pose `i` outside `S` may carry `pi_i = 1 - q_i`, where `q_i` is the weight of the heaviest
    independent set of `G[S]` that `i` can join: adding one such pose leaves every independent set
    at weight `<= 1`.  The great majority of poses have `q_i >= 1` -- a greedy independent set among
    the ones they do not meet already reaches 1 -- and are dismissed without any search; only the
    rest pay an exact max-weight independent set (`leaf_ceiling.max_clique_int` on the complement),
    memoised by the pose's adjacency pattern to `S`.  Several added poses at once need not keep
    validity, so the extended vector is renormalised by its OWN exact max-weight independent set,
    which leaves a valid cut whatever happened (and `None` if that B&B does not complete)."""
    nw = sub_win.shape[0]
    Sarr = np.array(S_local, dtype=np.int64)
    piS = np.asarray(pi_S, float)
    pi = np.zeros(nw)
    pi[Sarr] = piS
    subS = sub_win[np.ix_(Sarr, Sarr)]
    nbrS = comp_nbr(subS)
    wS = np.round(piS * 10 ** 9).astype(np.int64)
    adjS = sub_win[:, Sarr]
    ordS = np.argsort(-piS)
    memo = {}
    nadd = 0
    # the valuable extensions are the poses that meet MOST of the support (a pose meeting all of it
    # has q_i = 0 and carries the full weight 1 -- that is exactly `cliquelever`'s greedy maximal
    # clique extension when pi is the indicator of a clique), so take them in that order and cap.
    for i in np.argsort(-adjS.sum(1)):
        i = int(i)
        if pi[i] > 0 or nadd >= a.cut_extend:
            continue
        free = ~adjS[i]
        if not free.any():
            pi[i] = 1.0
            nadd += 1
            continue
        if adjS[i].sum() < a.cut_adj * len(Sarr):
            break
        alive = free.copy()
        tot = 0.0
        for v in ordS:
            if alive[v]:
                tot += piS[v]
                if tot >= 1.0 - 1e-12:
                    break
                alive &= ~subS[v]
                alive[v] = False
        if tot >= 1.0 - 1e-12:
            continue                                    # a greedy independent set already pays 1
        key = free.tobytes()
        if key in memo:
            q = memo[key]
        else:
            fr = np.nonzero(free)[0]
            pos = {int(v): k for k, v in enumerate(fr)}
            sub2 = [set(pos[u] for u in nbrS[int(v)] if u in pos) for v in fr]
            q, _m, _c = lc.max_clique_int(sub2, [int(wS[k]) for k in fr], 5.0, lb=0)
            memo[key] = q
        v = 1.0 - q / 10 ** 9
        if v > a.cut_min:
            pi[i] = v
            nadd += 1
    T = np.nonzero(pi > 1e-12)[0]
    wT = np.round(pi[T] * 10 ** 9).astype(np.int64)
    v, _m, comp = mwis_exact(sub_win[np.ix_(T, T)], list(wT), 20.0, lb=0)
    if not comp:
        return None, 0, False
    scale = max(1.0, v / 10 ** 9)
    return pi / scale, len(T), True


class Loc:
    """the LOC_D LP in the MU space: coverage rows, region pins, chord rows and the window cuts
    produced by `cover_lp_float`.  Every row is valid for LOC_D, so the LP value is an UPPER bound
    on `LOC_D(P)` at every iteration; when no window separates and no arrangement vertex is
    violated the solution is LOC_D-feasible and the value IS `LOC_D(P)`."""

    def __init__(self, ps, wins, ADJ, a, log):
        self.ps, self.wins, self.ADJ, self.a, self.log = ps, wins, ADJ, a, log
        self.t = ps.t
        self.target = lc.region_targets(self.t, a)
        self.nwin = len(wins)
        self.wpos = [{int(v): k for k, v in enumerate(P)} for _, P in wins]
        self.subs = [None] * self.nwin            # lazily built dense window subgraphs
        self.wcols = [set() for _ in range(self.nwin)]   # generated independent sets, by pose id
        self.rows = []                            # coverage rows (X, Y, D)
        self.rowkey = set()
        self.cuts = []                            # (window, pose array, pi array)
        self.cutkey = set()
        self.implied = bool(getattr(a, 'implied', False))
        self.m = Model(a.threads, log, getattr(a, 'crossover', False))
        self.build()

    def sub(self, j):
        if self.subs[j] is None:
            P = self.wins[j][1]
            self.subs[j] = self.ADJ[np.ix_(P, P)]
        return self.subs[j]

    def build(self):
        ps, a = self.ps, self.a
        m = self.m
        m.start(ps.ncol)
        INF = m.INF
        self.nE = 0
        if self.target:
            R, v = [], []
            reg = ps.col_reg_a
            for lab, kv in sorted(self.target.items()):
                R.append((reg == LABS.index(lab)).astype(float))
                v.append(float(kv))
            m.add_rows(np.array(v), np.array(v), sp.csr_matrix(np.array(R)))
            self.nE = len(v)
        self.nD = 0
        if a.chord:
            t, r = float(ps.t), float(ps.r)
            M = np.zeros((4, ps.ncol))
            for c, i in enumerate(ps.col_pose):
                _, _, cx, cy = ps.P[i]
                cx, cy = float(cx), float(cy)
                if cy <= r:
                    M[0, c] = 1
                if cy >= t - r:
                    M[1, c] = 1
                if cx <= r:
                    M[2, c] = 1
                if cx >= t - r:
                    M[3, c] = 1
            m.add_rows(np.full(4, -INF), np.full(4, 3.0), sp.csr_matrix(M))
            self.nD = 4
        self.log(f'   model: {m.nrow} rows ({self.nE} region, {self.nD} chord), {m.ncol} columns')

    def add_cov_rows(self, pts):
        new = [p for p in pts if p not in self.rowkey]
        if not new:
            return 0
        for p in new:
            self.rowkey.add(p)
        A, _ = self.ps.contains_rows(new)
        self.m.add_rows(np.full(len(new), -self.m.INF), np.ones(len(new)), A)
        self.rows += new
        return len(new)

    def add_cuts(self, cuts):
        if not cuts:
            return 0
        rr, cc, vv = [], [], []
        keep = []
        for (j, T, pi) in cuts:
            key = (j, tuple(int(v) for v in T), tuple(round(float(v), 6) for v in pi))
            if key in self.cutkey:
                continue
            self.cutkey.add(key)
            k = len(keep)
            keep.append((j, T, pi))
            for i, v in zip(T, pi):
                for c in self.ps.cols_of[int(i)]:
                    rr.append(k)
                    cc.append(c)
                    vv.append(float(v))
        if not keep:
            return 0
        A = sp.csr_matrix((vv, (rr, cc)), shape=(len(keep), self.ps.ncol))
        self.m.add_rows(np.full(len(keep), -self.m.INF), np.ones(len(keep)), A)
        self.cuts += keep
        return len(keep)

    def resume_rows(self, path, tight=0.0, mu0=None):
        """exact coverage rows `X Y D` from a `cl_*_rows.txt` checkpoint.  A coverage row is valid
        for LOC_D (it is part of its definition), so seeding with the rows a `cliquelever` run
        accumulated only starts the descent lower -- it never changes the value.

        `tight`: keep only the rows the SEED measure `mu0` already saturates (coverage
        `>= 1 - tight`).  A converged `cliquelever` measure saturates a few thousand of its 25,012
        rows and the rest carry no dual; the LP is quadratic in the row count and those rows are the
        densest in the model, so the tight subset is the same warm start at a fraction of the cost.
        Nothing about validity changes -- the dropped rows are still implied by the windows
        whenever `D(1-pitch) >= sqrt 2`, and the coverage separation would re-add any that bite."""
        pts = []
        for line in open(path):
            q = line.split()
            if len(q) == 3 and not q[0].startswith('#'):
                pts.append((int(q[0]), int(q[1]), int(q[2])))
        if tight > 0 and mu0 is not None and mu0.sum() > 0:
            A, _ = self.ps.contains_rows(pts)
            w = np.zeros(self.ps.ncol)
            for i in range(self.ps.n):
                cs = self.ps.cols_of[i]
                if cs:
                    w[cs[0]] = mu0[i]
            cov = A.dot(w)
            keep = np.nonzero(cov >= 1.0 - tight)[0]
            self.log(f'   [resume: {len(keep)} of {len(pts)} coverage rows are tight for the seed '
                     f'measure (coverage >= {1 - tight:.6f}); the rest are dropped]')
            pts = [pts[i] for i in keep]
        return self.add_cov_rows(pts)

    def resume_cliques(self, path):
        """clique rows from a `cl_*_cliques.txt` checkpoint: pose indices into the pose file the
        checkpoint was written with, which must be this run's `--poses`.  `mu(K) <= 1` is the
        `pi = 1_K` case of a window cut and is implied by LOC_D whenever the clique fits in a window
        (`D(1-pitch) >= sqrt 2`), so these rows are valid and only accelerate the descent.  Every
        row is re-verified pairwise before it enters the LP."""
        n, bad = 0, 0
        rows = []
        for line in open(path):
            q = line.split()
            if not q or q[0].startswith('#'):
                continue
            mem = [int(v) for v in q]
            if max(mem) >= self.ps.n:
                raise SystemExit(f'{path}: pose index {max(mem)} beyond the {self.ps.n} poses')
            if self.ps.verify_clique(mem):
                bad += 1
                continue
            rows.append((-1, np.array(mem, dtype=np.int64), np.ones(len(mem))))
            n += 1
        if bad:
            self.log(f'   [resume: {bad} checkpoint cliques FAILED the pairwise check, dropped]')
        self.add_cuts(rows)
        return n

    def pose_mass(self, x):
        mu = np.zeros(self.ps.n)
        np.add.at(mu, self.ps.col_pose_a, x[:self.ps.ncol])
        return mu

    def region_split(self, x):
        out = np.zeros(13)
        np.add.at(out, self.ps.col_reg_a, x[:self.ps.ncol])
        return out

    def certify_cov(self, mu, procs):
        sup = np.nonzero(mu > 1e-12)[0]
        squares = [self.ps.sq[i][:8] + (int(round(mu[i] * DM)),) for i in sup]
        V = lc.enumerate_vertices(squares, self.t, verbose=False, procs=procs)
        Inc = lc.incidences(squares, V, verbose=False, procs=procs)
        w = np.array([s[8] for s in squares], dtype=np.int64)
        cov = Inc.cov(w)
        M = float(Fr(int(cov.max()), DM)) if len(cov) else 0.0
        return sup, V, cov, M

    # ---------------------------------------------------------------- separation over the windows
    def separate(self, mu):
        """one exact window subproblem per window; returns (cuts, max z, all complete, z list)"""
        a = self.a
        cuts, zmax, allc, zs = [], 0.0, True, []
        capped = False
        for j, (_, P) in enumerate(self.wins):
            loc = np.nonzero(mu[P] > a.mu_tol)[0]
            if len(loc) < 2:
                zs.append(float(mu[P].sum()))
                continue
            if len(loc) > a.sep_cap:
                # early iterations spread the mass over thousands of poses; the subproblem on the
                # heaviest `--sep-cap` of them still yields a VALID cut (it is a valid inequality
                # for the sub-support, hence for the window), only the "nothing separates" test is
                # then not a proof -- which is what `capped` records.
                loc = loc[np.argsort(-mu[P][loc])[:a.sep_cap]]
                capped = True
            sub = self.sub(j)
            subS = sub[np.ix_(loc, loc)]
            inv = {int(P[k]): kk for kk, k in enumerate(loc)}
            seedc = set()
            for I in self.wcols[j]:
                J = frozenset(inv[i] for i in I if i in inv)
                if len(J) > 1:
                    seedc.add(J)
            z, pi, ncols, comp = cover_lp_cg(subS, mu[P][loc], tol=a.cut_tol,
                                             cols=sorted(seedc, key=sorted), budget=a.sep_time)
            keep = set()
            for I in ncols:
                keep.add(frozenset(int(P[loc[k]]) for k in I))
            self.wcols[j] = set(sorted(keep, key=sorted)[-a.warm_cols:]) if len(keep) > a.warm_cols else keep
            allc &= bool(comp)
            if z is None:
                allc = False
                zs.append(float('nan'))
                continue
            zs.append(z)
            zmax = max(zmax, z)
            if z > 1 + a.cut_tol:
                T = loc[pi > 1e-12]
                cuts.append((j, P[T], pi[pi > 1e-12]))      # the raw cut: valid, violated by z-1
                if a.cut_extend > 0:
                    pe, nT, okx = extend_cut(sub, loc, pi, a)
                    if okx:
                        Te = np.nonzero(pe > 1e-12)[0]
                        cuts.append((j, P[Te], pe[Te]))     # and its maximal extension
        return cuts, zmax, allc and not capped, zs

    # ---------------------------------------------------------------- the loop
    def run(self, tag):
        a = self.a
        T0 = time.time()
        hist = []
        x = None
        for it in range(a.iters):
            res = self.m.run(a.lp_tlim)
            if res is None:
                return hist, x
            x, _y, val = res
            mu = self.pose_mass(x)
            nrow, M, sup, V = 0, float('nan'), np.nonzero(mu > 1e-12)[0], []
            # The coverage rows are IMPLIED by the window constraints whenever D(1-pitch) >= sqrt 2
            # (every point clique fits in a window), so separating them during the descent is pure
            # cost; the final measure is checked for M <= 1 exactly regardless.  When they are NOT
            # implied they are separated at the arrangement vertices of the support, as
            # `cliquelever` does -- but only while the support is small enough for the exact
            # arrangement to be affordable.
            if self.a.cov_sep and len(sup) <= a.cov_sup:
                sup, V, cov, M = self.certify_cov(mu, a.procs)
                bad = np.nonzero(cov > DM * (1.0 + 1e-9))[0]
                if len(bad):
                    o = bad[np.argsort(-cov[bad])[:a.row_cap]]
                    nrow = self.add_cov_rows([V[i] for i in o])
            t1 = time.time()
            cuts, zmax, allc, zs = self.separate(mu)
            ncut = self.add_cuts(cuts)
            R = self.region_split(x)
            done = (nrow == 0 and ncut == 0 and zmax <= 1 + a.cut_tol
                    and (self.implied or (M == M and M <= 1 + 1e-9)))
            rec = dict(tag=tag, it=it, LP=float(val), M=M, zmax=float(zmax), complete=bool(allc),
                       rows=len(self.rows), cuts=len(self.cuts), nsup=int((mu > 1e-12).sum()),
                       nverts=len(V), interior=float(R[12]), corners=R[:4].tolist(),
                       slots=R[4:12].tolist(), z=[float(v) for v in zs],
                       secs=time.time() - T0, converged=bool(done and allc))
            hist.append(rec)
            self.log(f'   [{tag}.{it}] LP={val:.6f} M={M:.9f} zmax={zmax:.6f}{"" if allc else "*"} '
                     f'rows={len(self.rows)}(+{nrow}) cuts={len(self.cuts)}(+{ncut}) '
                     f'sup={int((mu > 1e-12).sum())} int={R[12]:.6f} '
                     f'(sep {time.time() - t1:.0f}s, {time.time() - T0:.0f}s)')
            if done:
                if not allc:
                    self.log('   [nothing separates, but an enumeration/B&B was capped: NOT a proof]')
                break
            if time.time() - T0 > a.time:
                self.log('   [time limit inside the loop]')
                break
            if a.ckpt and it % a.ckpt == a.ckpt - 1:
                json.dump(hist, open(os.path.join(RUNS, f'loc_{a.TAG}_hist.json'), 'w'), indent=0)
        self.last_x = x
        return hist, x


# ====================================================================== the shape of a window cut
def cut_anatomy(ps, ADJ, P, S, mi, log, tag=''):
    """the separating inequality of one window, named.

    `S` are the support poses of the window; the window subproblem's optimal dual `pi >= 0` has
    `sum_{i in I} pi_i <= 1` for every independent set, so `sum pi_i mu_i <= 1` is the cut.  Its
    SHAPE is read off the values `pi` takes: `pi = 1` on a clique is a clique row; `pi = 1/2` on
    five classes with `alpha = 2` is `RANKDIAG.md`'s pentagon (`alpha(C5) = 2`); `1/3` on seven
    classes with `alpha = 3` is the heptagon; anything else is new."""
    sub = ADJ[np.ix_(S, S)]
    muf = np.array([mi[i] / DM for i in S])
    z, pi, cols, comp = cover_lp_cg(sub, muf, tol=1e-9)
    if pi is None:
        return None
    T = np.nonzero(pi > 1e-9)[0]
    aw, am, ac = mwis_exact(sub[np.ix_(T, T)], [1] * len(T), 60.0, lb=0)
    vals = {}
    for k in T:
        vals.setdefault(round(float(pi[k]), 6), []).append(int(k))
    out = dict(z=float(z), complete=bool(comp), support=len(T), alpha=int(aw),
               alpha_complete=bool(ac), classes=[], pi_values=sorted(vals, reverse=True))
    log(f'   {tag}cut: z = {z:.9f}, {len(T)} poses carry pi > 0, alpha(cut support) = {aw} '
        f'(complete {ac}), {len(vals)} distinct pi values')
    for v in sorted(vals, reverse=True):
        ks = vals[v]
        cx = np.mean([float(ps.P[S[k]][2]) for k in ks])
        cy = np.mean([float(ps.P[S[k]][3]) for k in ks])
        th = sorted(round(math.degrees(math.atan2(ps.F[S[k], 3], ps.F[S[k], 2])) % 90, 2)
                    for k in ks)
        m = sum(mi[S[k]] for k in ks) / DM
        regs = {}
        for k in ks:
            lab = lc.regions_of(ps.P[S[k]][2], ps.P[S[k]][3], ps.boxes)[0]
            regs[lab] = regs.get(lab, 0) + 1
        out['classes'].append(dict(pi=v, poses=len(ks), mu=m, centre=[round(cx, 4), round(cy, 4)],
                                   theta=[th[0], th[-1]], regions=regs))
        log(f'     pi = {v:.6f} x {len(ks):3d} poses, mu = {m:.6f}, centroid '
            f'({cx:.3f}, {cy:.3f}), theta {th[0]:.2f}..{th[-1]:.2f} deg, regions {regs}')
    log(f'     -> pi . mu = {sum(pi[k] * muf[k] for k in T):.9f}; the inequality is '
        f'"{describe_cut(out)}"')
    out['name'] = describe_cut(out)
    return out


def describe_cut(out):
    """name the family a cut belongs to from its dual values and the independence number of its
    support: 1/r on 2r+1 classes with alpha = r is the odd (2r+1)-gon rank inequality."""
    vs = out['pi_values']
    a = out['alpha']
    if len(vs) == 1 and abs(vs[0] - 1.0) < 1e-6 and a == 1:
        return 'a clique row mu(K) <= 1'
    if len(vs) == 1 and a >= 2 and abs(vs[0] - 1.0 / a) < 1e-6:
        n = len(out['classes'][0]['poses'] and out['classes']) and out['support']
        return (f'a rank inequality at level 1/{a}: mu(X) <= {a} over {n} poses with alpha = {a} '
                f'-- the odd {2 * a + 1}-gon shape of RANKDIAG.md if the classes number {2 * a + 1}')
    return (f'mixed: pi takes {len(vs)} values {["%.4f" % v for v in vs[:6]]} on a support with '
            f'alpha = {a}')


# ====================================================================== exact certification
def window_certificates(ps, ADJ, wins, mi, log, cut_tol=1e-9):
    """for every window, an EXACT rational independent-set decomposition of the restriction of the
    integer measure `mi` (numerators over DM).  The columns come from `cover_lp_cg` (whose last
    pricing max-weight independent set PROVES that no column is missing when it is complete), and
    the decomposition itself is then recomputed by the exact rational simplex `exact_cover_lp`, so
    the `lam` in the certificate are rationals and the test `sum lam <= 1` is exact.  Every set in
    the certificate's support is re-verified pairwise disjoint by the exact integer SAT."""
    out = []
    okall = True
    for j, (geo, P) in enumerate(wins):
        S = [int(i) for i in P if mi[i] > 0]
        if not S:
            out.append(dict(window=j, geo=list(geo), poses=0, nsets=0, value='0', value_float=0.0,
                            ok=True, complete=True, mass='0', mass_float=0.0, support=0,
                            bad_pairs=0))
            continue
        sub = ADJ[np.ix_(S, S)]
        muf = np.array([mi[i] / DM for i in S])
        _z, _pi, cols, comp = cover_lp_cg(sub, muf, tol=cut_tol)
        cols = [frozenset(int(v) for v in I) for I in cols]
        mu = [Fr(int(mi[i]), DM) for i in S]
        val, sol = exact_cover_lp(cols, mu, len(S))
        ok = val is not None and val <= 1
        badp = 0
        for cj in (sol or {}):
            mem = sorted(cols[cj])
            for x in range(len(mem)):
                for y in range(x + 1, len(mem)):
                    if lc.sq_meets_sq(ps.sq[S[mem[x]]], ps.sq[S[mem[y]]]):
                        badp += 1
        ok = bool(ok and badp == 0)
        okall &= ok
        out.append(dict(window=j, geo=list(geo), poses=len(S), nsets=len(cols),
                        complete=bool(comp), value=str(val), value_float=float(val),
                        support=len(sol or {}), bad_pairs=badp, ok=ok,
                        mass=str(sum(mu)), mass_float=float(sum(mu))))
        log(f'     window {j:3d} {geo}: {len(S)} support poses, {len(cols)} independent-set '
            f'columns{"" if comp else " (pricing B&B INCOMPLETE)"}, mass {float(sum(mu)):.6f}, '
            f'sum lam = {val} = {float(val):.9f} {"OK" if ok else "FAIL"}'
            + (f' [{badp} meeting pairs!]' if badp else ''))
    return out, okall


def finalize(loc, x, log, procs):
    """round down to /DM, top up the pinned regions on exact WINDOW slack, then certify M, the
    region counts and every window's independent-set decomposition exactly."""
    ps = loc.ps
    xc = [int(math.floor(v * DM)) for v in x[:ps.ncol]]

    def pose_int():
        mi = [0] * ps.n
        for c, v in enumerate(xc):
            mi[ps.col_pose[c]] += v
        return mi

    def arrangement(mi):
        sup = [i for i in range(ps.n) if mi[i] > 0]
        squares = [ps.sq[i][:8] + (mi[i],) for i in sup]
        V = lc.enumerate_vertices(squares, loc.t, verbose=False, procs=procs)
        Inc = lc.incidences(squares, V, verbose=False, procs=procs)
        w = np.array([s[8] for s in squares], dtype=np.int64)
        return sup, squares, V, Inc, w, Inc.cov(w)

    mi = pose_int()
    sup, squares, V, Inc, w, cov = arrangement(mi)
    M = Fr(int(cov.max()), DM)
    log(f'   final: rounded down: mass {sum(mi) / DM:.9f}, exact M = {float(M):.12f}')
    if M > 1:
        log('   final: M > 1 after rounding down; scaling by 1/M')
        xc = [(v * M.denominator) // M.numerator for v in xc]
        mi = pose_int()
        sup, squares, V, Inc, w, cov = arrangement(mi)
        M = Fr(int(cov.max()), DM)
    # ---- region top-up on exact window slack
    if loc.target:
        cert, _ = window_certificates(ps, loc.ADJ, loc.wins, mi, lambda s: None)
        slack = {}
        for j, c in enumerate(cert):
            v = Fr(c['value']) if c['value'] not in (None, 'None') else Fr(2)
            slack[j] = Fr(1) - v
        wof = [[] for _ in range(ps.n)]
        for j, (_, P) in enumerate(loc.wins):
            for i in P:
                wof[int(i)].append(j)
        for lab, kv in sorted(loc.target.items()):
            li = LABS.index(lab)
            cols = [c for c in range(ps.ncol) if ps.col_reg[c] == li]
            d = int(kv * DM) - sum(xc[c] for c in cols)
            if d <= 0:
                log(f'   final: region {lab} residual {d / DM:+.3e}')
                continue
            for c in sorted(cols, key=lambda c: -xc[c]):
                if d <= 0:
                    break
                i = ps.col_pose[c]
                if mi[i] <= 0:
                    continue
                room = min([slack[j] for j in wof[i]] + [Fr(1)])
                k = int(cov[Inc.cols[sup.index(i)]].max()) if i in sup else 0
                room = min(room, Fr(DM - k, DM))
                add = min(d, int(room * DM))
                if add <= 0:
                    continue
                xc[c] += add
                d -= add
                for j in wof[i]:
                    slack[j] -= Fr(add, DM)
                mi = pose_int()
                sup, squares, V, Inc, w, cov = arrangement(mi)
            log(f'   final: region {lab} topped up, residual {d / DM:+.3e}')
    mi = pose_int()
    sup, squares, V, Inc, w, cov = arrangement(mi)
    M = Fr(int(cov.max()), DM)
    mass = Fr(sum(mi), DM)
    meta = [(Fr(s[0], s[2]), Fr(s[1], s[2]), s[8]) for s in squares]
    ok, assigned, amb, note = lc.region_report(meta, ps.boxes, loc.target, DM)
    nbr = lc.closed_graph(squares, verbose=False)
    Kw, Km, comp = lc.max_clique_int(nbr, [int(v) for v in w], 120.0, lb=0)
    Kq = Fr(Kw, DM)
    aw, am, acomp = mwis_exact(np.array([[b in nbr[a] for b in range(len(sup))]
                                         for a in range(len(sup))]),
                               [1] * len(sup), 300.0, lb=0)
    log(f'   FINAL EXACT: mass = {mass} = {float(mass):.9f}; M = {M} = {float(M):.12f} '
        f'{"OK" if M <= 1 else "FAIL"}; max clique = {float(Kq):.9f} (complete {comp}); '
        f'alpha(support) = {aw} (complete {acomp}); regions {"OK" if ok else "FAIL"} ({note}); '
        f'{len(sup)} poses')
    log('   window certificates:')
    cert, okw = window_certificates(ps, loc.ADJ, loc.wins, mi, log)
    res = dict(mass=str(mass), mass_float=float(mass), M=str(M), M_float=float(M),
               clique=str(Kq), clique_float=float(Kq), clique_complete=bool(comp),
               alpha=int(aw), alpha_complete=bool(acomp),
               regions_ok=bool(ok), poses=len(sup), vertices=len(V),
               regions={lab: str(Fr(assigned.get(lab, 0), DM)) for lab in LABS},
               windows=cert, windows_ok=bool(okw),
               certified=bool(M <= 1 and ok and okw))
    return res, mi


# ====================================================================== anatomy
def anatomy(loc, x, mu, log, res, hist):
    """where the windows disagree: which convexity rows are tight, which poses have
    mu_i < sum_{I ni i} lam, and the region/geometry shape of the excess over alpha"""
    ps = loc.ps
    out = {}
    zs = (hist[-1]['z'] if hist else [])
    out['window_mass'] = []
    for j, (geo, P) in enumerate(loc.wins):
        out['window_mass'].append(dict(window=j, geo=list(geo), poses=int(len(P)),
                                       mu=float(mu[P].sum()), nsup=int((mu[P] > 1e-9).sum()),
                                       z=(float(zs[j]) if j < len(zs) else None)))
    log('   window anatomy (z = the convexity multiplier the window needs; 1 = saturated):')
    for d in out['window_mass']:
        log(f'     {d["window"]:3d} {tuple(d["geo"])}: {d["poses"]:5d} poses '
            f'({d["nsup"]} in the support), mu = {d["mu"]:8.5f}, z = '
            + ('n/a' if d['z'] is None else f'{d["z"]:.6f}'))
    sup = np.nonzero(mu > 1e-9)[0]
    reg = {}
    for i in sup:
        lab = lc.regions_of(ps.P[i][2], ps.P[i][3], ps.boxes)[0]
        reg[lab] = reg.get(lab, 0.0) + float(mu[i])
    out['regions'] = {k: round(v, 6) for k, v in sorted(reg.items())}
    out['support'] = int(len(sup))
    out['alpha'] = res.get('alpha')
    out['excess'] = float(mu.sum()) - (res.get('alpha') or 0)
    return out


# ====================================================================== driver
def cmd_run(a):
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'loc_{a.TAG}.log'), 'a')

    def log(m):
        print(m, flush=True)
        logf.write(str(m) + '\n')
        logf.flush()

    T0 = time.time()
    ps = Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    n = ps.add_exact_file(a.poses)
    log(f'# locality {a.TAG}: +{n} poses from {os.path.basename(a.poses)}')
    ps.finish()
    log(f'# t={ps.t} r={ps.r} corners={a.corners} patterns={a.patterns} chord={a.chord} '
        f'D={a.D} shape={a.shape} pitch={a.pitch}; {ps.n} poses, {ps.ncol} columns, '
        f'initial mass {sum(ps.mu0):.6f}; args={vars(a)}')
    ADJ = build_adj(ps, log)
    wins = make_windows(ps, a.D, a.shape, a.pitch, log)
    # Every point clique {S : x in S} has its centres inside a box of side sqrt(2) (a unit square
    # containing x has its centre within sqrt(2)/2 of x), so as soon as D - pitch*D >= sqrt(2) that
    # clique lies inside SOME window of the grid and the coverage row  cov(x) <= 1  is IMPLIED by
    # mu|_{P_j} in Pi_j (an independent set meets a clique at most once).  The rows are then dead
    # weight -- and they are by far the densest rows in the model -- so they are not seeded; the
    # separation below still runs every iteration and would add any that were violated.
    a.implied = bool(a.D * (1.0 - a.pitch) >= math.sqrt(2) - 1e-12 and a.shape == 'box')
    if a.cov_sep is None:
        a.cov_sep = not a.implied
    log(f'   coverage rows implied by the windows: {a.implied} '
        f'(D(1-pitch) = {a.D * (1 - a.pitch):.4f} vs sqrt 2 = {math.sqrt(2):.4f}); '
        f'coverage separation {"on" if a.cov_sep else "off"}')
    loc = Loc(ps, wins, ADJ, a, log)
    mu0 = np.array(ps.mu0)
    # seed rows and columns
    if a.row_pitch > 0:
        Dg = 1000
        step = int(round(a.row_pitch * Dg))
        tn = int(float(ps.t) * Dg)
        pts = [lc.reduce3(X, Y, Dg) for X in range(step // 2, tn, step)
               for Y in range(step // 2, tn, step)]
        log(f'   +{loc.add_cov_rows(pts)} grid coverage rows (pitch {a.row_pitch})')
    if mu0.sum() > 0 and a.rows0 > 0:
        sup, V, cov, M0 = loc.certify_cov(mu0, a.procs)
        o = np.argsort(-cov)[:a.rows0]
        log(f'   +{loc.add_cov_rows([V[i] for i in o])} rows at the heaviest arrangement vertices '
            f'of the seed support ({len(sup)} poses, {len(V)} vertices, M0 = {M0:.6f})')
    for f in a.resume_rows:
        log(f'   +{loc.resume_rows(f, a.tight, mu0)} exact coverage rows resumed from '
            f'{os.path.basename(f)}')
    for f in a.resume_cliques:
        log(f'   +{loc.resume_cliques(f)} clique rows resumed from {os.path.basename(f)} '
            f'(each re-verified pairwise; valid for LOC_D, and implied by it when '
            f'D(1-pitch) >= sqrt 2)')
    # A seed of coverage rows is needed even when the windows imply them: without any row the LP
    # is unbounded on the free interior.  The `--row-pitch` grid does it; separation adds the rest.
    if not loc.rows:
        Dg = 1000
        step = int(round(0.08 * Dg))
        tn = int(float(ps.t) * Dg)
        pts = [lc.reduce3(X, Y, Dg) for X in range(step // 2, tn, step)
               for Y in range(step // 2, tn, step)]
        log(f'   +{loc.add_cov_rows(pts)} grid coverage rows (pitch 0.08, bounding rows)')
    log(f'   setup done, {time.time() - T0:.0f}s')

    hist, x = loc.run('s0')
    if x is None:
        raise SystemExit('LP failed')
    mu = loc.pose_mass(x)
    log(f'   loop done: LP = {hist[-1]["LP"]:.9f}, converged = {hist[-1]["converged"]}, '
        f'{time.time() - T0:.0f}s')
    res, mi = finalize(loc, x, log, a.procs)
    path = os.path.join(RUNS, f'loc_{a.TAG}_exact.txt')
    with open(path, 'w') as f:
        f.write(f'# locality LOC_D measure; t = {loc.t}; sym 1; D = {a.D} shape = {a.shape}; '
                f'exactly certified: M = {res["M"]}, windows_ok = {res["windows_ok"]}\n')
        f.write(f'# mass = {res["mass"]} = {res["mass_float"]:.12f}\n')
        f.write('# pose p q cx cy mass : theta = 2 arctan(p/q)\n')
        for i, m in enumerate(mi):
            if m > 0:
                p, q, cx, cy = ps.P[i]
                f.write(f'pose {p} {q} {cx} {cy} {Fr(int(m), DM)}\n')
    res['file'] = path
    ana = anatomy(loc, x, mu, log, res, hist)
    out = dict(tag=a.TAG, args={k: (str(v) if isinstance(v, Fr) else v) for k, v in vars(a).items()},
               windows=[[list(g), int(len(P))] for g, P in wins], hist=hist, final=res,
               anatomy=ana, secs=time.time() - T0)
    json.dump(out, open(os.path.join(RUNS, f'loc_{a.TAG}.json'), 'w'), indent=1)
    log(f'# DONE {a.TAG}: LOC_{a.D}({a.shape}) = {hist[-1]["LP"]:.9f} '
        f'(converged {hist[-1]["converged"]}), certified {res["mass_float"]:.9f}, '
        f'alpha = {res["alpha"]}, {time.time() - T0:.0f}s')


def cmd_check(a):
    """the window anatomy of a GIVEN measure: the exact convexity multiplier `z_j` every window
    needs (`mu|_{P_j} in Pi_j` iff `z_j <= 1`), by column generation plus the exact rational
    simplex.  `max_j z_j > 1` PROVES the measure is outside the LOC_D polytope."""
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'loc_{a.TAG}_check.log'), 'a')

    def log(m):
        print(m, flush=True)
        logf.write(str(m) + chr(10))
        logf.flush()

    ps = Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    n = ps.add_exact_file(a.poses)
    if os.path.abspath(a.measure) != os.path.abspath(a.poses):
        ps.add_exact_file(a.measure)
    ps.finish()
    t, sym, poses = lc.read_measure(a.measure)
    mi = [0] * ps.n
    for (p_, q_, cx, cy, m) in poses:
        i = ps.key.get((p_, q_, ps.bsnap(cx), ps.bsnap(cy)))
        assert i is not None, f'pose {(p_, q_, cx, cy)} of the measure is not in the pose set'
        assert Fr(int(m * DM), DM) == m, 'measure mass is not a multiple of 1e-9'
        mi[i] += int(m * DM)
    log(f'# loc check {a.TAG}: {n} poses, measure {os.path.basename(a.measure)} mass '
        f'{Fr(sum(mi), DM)} = {sum(mi) / DM:.9f} on {sum(1 for v in mi if v)} poses')
    ADJ = build_adj(ps, log)
    out = {}
    for shape in a.shape:
        for D in a.D:
            wins = make_windows(ps, D, shape, a.pitch, log)
            cert, ok = window_certificates(ps, ADJ, wins, mi, log)
            zs = [Fr(c['value']) for c in cert if c['poses']]
            zmax = max(zs) if zs else Fr(0)
            worst = max(cert, key=lambda c: float(c['value_float']))
            if float(worst['value_float']) > 1 + 1e-9:
                Pw = wins[worst['window']][1]
                Sw = [int(i) for i in Pw if mi[i] > 0]
                worst['anatomy'] = cut_anatomy(ps, ADJ, Pw, Sw, mi, log,
                                               tag=f'D={D} {shape} window {worst["window"]} ')
            log(f'## D = {D} {shape}: {len(wins)} windows, max_j z_j = {zmax} = '
                f'{float(zmax):.9f} at window {worst["window"]} {tuple(worst["geo"])} '
                f'(mass {worst["mass_float"]:.6f}, {worst["poses"]} support poses) -> the measure '
                f'is {"IN" if ok else "NOT in"} the LOC_{D} polytope; '
                f'mass/zmax = {sum(mi) / DM / max(float(zmax), 1.0):.6f}')
            out[f'{D}_{shape}'] = dict(D=D, shape=shape, windows=len(wins), zmax=str(zmax),
                                       zmax_float=float(zmax), ok=bool(ok), cert=cert)
    json.dump(out, open(os.path.join(RUNS, f'loc_{a.TAG}_check.json'), 'w'), indent=1)


# ====================================================================== selftest
def selftest():
    ok = [0, 0]

    def chk(name, c, extra=''):
        ok[1] += 1
        ok[0] += bool(c)
        print(f'  {"ok  " if c else "FAIL"} {name} {extra}')

    # exact cover LP: three pairwise-meeting squares with mass 1/3 each -> sum lam = 1
    cols = [frozenset([0]), frozenset([1]), frozenset([2])]
    v, sol = exact_cover_lp(cols, [Fr(1, 3)] * 3, 3)
    chk('exact_cover_lp clique of 3', v == 1, f'{v}')
    # two disjoint squares, mass 1/2 each, one independent set {0,1} -> sum lam = 1/2
    v, sol = exact_cover_lp([frozenset([0, 1])], [Fr(1, 2)] * 2, 2)
    chk('exact_cover_lp one pair', v == Fr(1, 2), f'{v}')
    # 5-cycle with mass 1/2 each: fractional chromatic number of C5 is 5/2, mass 5/2 -> 5/4
    C5 = [frozenset([0, 2]), frozenset([1, 3]), frozenset([2, 4]), frozenset([3, 0]),
          frozenset([4, 1])]
    v, sol = exact_cover_lp(C5, [Fr(1, 2)] * 5, 5)
    chk('exact_cover_lp C5 at 1/2', v == Fr(5, 4), f'{v}')
    # maximal independent set enumeration on C5 (5 maximal independent sets, each a pair)
    A = np.zeros((5, 5), bool)
    for i in range(5):
        A[i, (i + 1) % 5] = A[(i + 1) % 5, i] = True
    S, comp = maximal_is_enum(A)
    chk('maximal_is_enum C5', len(S) == 5 and comp and all(len(s) == 2 for s in S), f'{len(S)}')
    # mwis on C5 with weights 1,2,3,4,5 -> {1,3} (0-based) = 2+4 = 6? best is {0,2}=1+3=4,
    # {1,3}=2+4=6, {1,4}=2+5=7, {0,3}=1+4=5, {2,4}=3+5=8
    v, mem, comp = mwis_exact(A, [1, 2, 3, 4, 5], 10.0, lb=0)
    chk('mwis_exact C5', v == 8 and comp and set(mem) == {2, 4}, f'{v} {sorted(mem)}')
    # windows: a box of side 4 on a 4x4 container is one window containing everything
    print(f'  {ok[0]}/{ok[1]}')
    return 0 if ok[0] == ok[1] else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('run')
    c.add_argument('TAG')
    c.add_argument('--poses', required=True)
    c.add_argument('--t', default='4')
    c.add_argument('--r', default='1')
    c.add_argument('--bnd-delta', type=float, default=1e-6)
    c.add_argument('--corners', default='1111')
    c.add_argument('--patterns', default='01010101')
    c.add_argument('--chord', action='store_true')
    c.add_argument('--D', type=float, required=True)
    c.add_argument('--shape', default='box', choices=['box', 'disc'])
    c.add_argument('--pitch', type=float, default=0.25, help='window pitch as a fraction of D')
    c.add_argument('--resume-rows', action='append', default=[],
                   help='cl_*_rows.txt checkpoint (X Y D): valid coverage rows, seeds the descent')
    c.add_argument('--resume-cliques', action='append', default=[],
                   help='cl_*_cliques.txt checkpoint; --poses must be the pose file it indexes')
    c.add_argument('--crossover', action='store_true',
                   help='run the ipm crossover (off by default: the loop never reads the duals)')
    c.add_argument('--lp-tlim', type=float, default=0.0,
                   help='warm dual-simplex budget per solve after the first (0 = always ipm)')
    c.add_argument('--tight', type=float, default=0.0,
                   help='keep only resumed coverage rows the seed measure saturates to this slack')
    c.add_argument('--slack', type=float, default=0.0,
                   help='convexity rows are sum lam <= 1 - slack (room for the exact top-up)')
    c.add_argument('--threads', type=int, default=2)
    c.add_argument('--procs', type=int, default=2)
    c.add_argument('--iters', type=int, default=200)
    c.add_argument('--time', type=float, default=100000)
    c.add_argument('--row-pitch', type=float, default=0.0)
    c.add_argument('--rows0', type=int, default=0)
    c.add_argument('--row-cap', type=int, default=400)
    c.add_argument('--mu-tol', type=float, default=1e-9,
                   help='a pose below this carries no mass for the window subproblem')
    c.add_argument('--cut-tol', type=float, default=1e-7)
    c.add_argument('--cap-enum', type=int, default=200000,
                   help='cap on the maximal-independent-set enumeration of a window support')
    c.add_argument('--cut-adj', type=float, default=0.5,
                   help='only extend onto poses meeting at least this fraction of the cut support')
    c.add_argument('--cut-extend', type=int, default=400,
                   help='poses added to a cut by the extension pass (0 = no extension)')
    c.add_argument('--cut-min', type=float, default=0.02,
                   help='smallest weight an extended pose is given')
    c.add_argument('--ckpt', type=int, default=5)
    c.add_argument('--cov-sep', dest='cov_sep', action='store_true', default=None)
    c.add_argument('--no-cov-sep', dest='cov_sep', action='store_false')
    c.add_argument('--sep-time', type=float, default=25.0,
                   help='seconds per window subproblem')
    c.add_argument('--warm-cols', type=int, default=4000,
                   help='independent-set columns kept per window between iterations')
    c.add_argument('--sep-cap', type=int, default=140,
                   help='heaviest poses per window fed to the window subproblem')
    c.add_argument('--cov-sup', type=int, default=600,
                   help='skip the exact arrangement when the support is larger than this')
    c = sub.add_parser('check')
    c.add_argument('TAG')
    c.add_argument('--poses', required=True)
    c.add_argument('--measure', required=True)
    c.add_argument('--t', default='4')
    c.add_argument('--r', default='1')
    c.add_argument('--bnd-delta', type=float, default=1e-6)
    c.add_argument('--D', type=float, nargs='+', default=[1.5, 2.0, 2.5, 3.0, 3.5])
    c.add_argument('--shape', nargs='+', default=['box'])
    c.add_argument('--pitch', type=float, default=0.25)
    c.add_argument('--cut-tol', type=float, default=1e-9)
    sub.add_parser('selftest')
    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    if a.cmd == 'run':
        cmd_run(a)
    elif a.cmd == 'check':
        cmd_check(a)
    else:
        sys.exit(selftest())


if __name__ == '__main__':
    main()
