#!/usr/bin/env python3
"""cliquelever.py -- the QSTAB relaxation of a level-2 leaf on a FIXED pose set.

    max  sum mu
    s.t. cov(x)  <= 1     for every point x of [0,t]^2            (coverage)
         mu(K)   <= 1     for EVERY clique K of the closed-intersection graph of the poses
         mu(R)   =  k_R   for every pinned region R (corner boxes / wall slots, boundary poses
                          free to declare either side)
         [mu(wall strip) <= 3, the chord rows, optional]
         mu >= 0

over a fixed finite set of poses of closed unit squares in the closed container [0,t]^2, t = 4.
`search/T4LEAF.md` values the leaf 01010101 at 11.745 with the ANCHOR cliques K(p, A) only, and
`search/LEAF_CEILING.md` 5 measures that its certified measures violate GENERAL cliques by 0.54-0.65
while the anchor family sees at most half of that.  This tool answers: if every clique inequality
were enforced, where does the leaf LP go?  (`tasks/clique-lever/README.md`.)

Method: row generation on a fixed pose set, no pricing (one optional pricing stage at the end).
  * POSES are exact: `leaf_ceiling.py`'s rational poses (theta = 2 arctan(p/q), rational centres),
    read from an exact measure file and/or snapped from a `t4screen` float checkpoint.  A pose
    whose centre is within `--bnd-delta` of a region boundary is snapped ONTO the boundary and gets
    one column per adjacent region (`t4leaf.py`'s boundary semantics, exactly).
  * COVERAGE rows are exact arrangement vertices (`leaf_ceiling.enumerate_vertices`) of the
    SUPPORT of the current solution: complete for the support measure (LEAF_CEILING.md Lemma 1),
    so `M <= 1` certified here is `M <= 1` over the whole container.  Row coefficients over the
    whole pose set: float SAT with the exact integer test on every borderline pair.
  * CLIQUE rows: the exact closed-intersection graph of the support (`leaf_ceiling.closed_graph`,
    integer SAT), max-weight clique by exact branch and bound on integer weights
    (`leaf_ceiling.max_clique_int`), plus the max-weight clique THROUGH each of the heaviest
    support squares.  Each violated support clique is then EXTENDED greedily to a maximal clique
    over all poses -- a pose joins iff it closed-intersects every current member (float SAT with
    exact fallback) -- and every added row is re-verified pairwise before it enters the LP.
    Sound whether or not a B&B completes: an incomplete B&B only weakens the STOPPING test.
  * LP: HiGHS (`t4leaf.Hi`), warm simplex with an ipm fallback.

Directions.  The LP value with a finite subset of the rows is an UPPER bound on the QSTAB value
of this pose set (rows only relax), at every iteration and regardless of B&B completeness, as
long as every row is valid (coverage rows are; clique rows are verified).  When the loop converges
(`M <= 1`, max clique `<= 1 + ktol` with a COMPLETE B&B) the solution is QSTAB-feasible on the
pose set, hence the value is also a LOWER bound on the true QSTAB value over ALL poses; the final
exact certification (masses rounded down, regions topped up on exact slack, M and the max clique
recomputed in integers) turns that into an exact statement.

    python3 search/cliquelever.py run TAG --exact runs/lc_A0101.txt --corners 1111 --patterns 01010101
    python3 search/cliquelever.py run TAG --exact runs/lc_A0101.txt --load-poses runs/tl_A01010101_poses.txt \
        --load-rows runs/tl_A01010101_dual.txt --corners 1111 --patterns 01010101 --chord --price 1
    python3 search/cliquelever.py selftest
"""
import argparse
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import leaf_ceiling as lc                                            # noqa: E402
import rankfamily                                                    # noqa: E402

RUNS = os.path.join(REPO, 'runs')
LABS = ['C0', 'C1', 'C2', 'C3', 'W0', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'W7', 'I']
FTOL = 1e-7          # float SAT band: |margin| <= FTOL is decided by the exact integer test


# ====================================================================== the pose set
class Poses:
    """exact poses, their squares, float shadows, and the region columns"""

    def __init__(self, t, r, delta):
        self.t = Fr(t)
        self.r = Fr(r)
        self.delta = float(delta)
        self.boxes = lc.region_boxes(self.t, self.r)
        self.P = []            # (p, q, cx, cy) exact
        self.key = {}
        self.sq = []           # exact square records
        self.mu0 = []          # initial (float) mass per pose, from the source files
        self.cols_of = []      # pose -> list of column indices
        self.col_pose = []     # column -> pose
        self.col_reg = []      # column -> region index (0..12)
        self.F = np.zeros((0, 4))   # float (cx, cy, cos, sin)

    # ---------------------------------------------------------------- input
    def bsnap(self, x):
        """snap a coordinate within delta of a region boundary onto it (exactly)"""
        for b in (self.r, self.t / 2, self.t - self.r):
            if abs(float(x - b)) <= self.delta:
                return b
        return x

    def add(self, p, q, cx, cy, mass=0.0):
        cx, cy = self.bsnap(Fr(cx)), self.bsnap(Fr(cy))
        a, b, rr = lc.cos_sin(p, q)
        w2 = lc.half_width(a, b, rr)
        cx = min(max(cx, w2), self.t - w2)
        cy = min(max(cy, w2), self.t - w2)
        k = (p, q, cx, cy)
        if k in self.key:
            i = self.key[k]
            self.mu0[i] += float(mass)
            return i
        i = len(self.P)
        self.key[k] = i
        self.P.append(k)
        self.sq.append(lc.make_square(cx, cy, p, q, self.t, 0, tag=i))
        self.mu0.append(float(mass))
        regs = lc.regions_of(cx, cy, self.boxes)
        self.cols_of.append([])
        for lab in regs:
            self.cols_of[i].append(len(self.col_pose))
            self.col_pose.append(i)
            self.col_reg.append(LABS.index(lab))
        return i

    def add_exact_file(self, path):
        t, sym, poses = lc.read_measure(path)
        assert t == self.t, (t, self.t)
        assert sym == 1
        n0 = len(self.P)
        for (p, q, cx, cy, m) in poses:
            self.add(p, q, cx, cy, float(m))
        return len(self.P) - n0

    def add_float_file(self, path, Q, Dc):
        tf, src = lc.read_float_poses(path)
        n0 = len(self.P)
        for (cx, cy, th, mu) in src:
            p, q, x, y = lc.snap_pose(cx, cy, th, self.t, Q, Dc)
            self.add(p, q, x, y, mu)
        return len(self.P) - n0

    def add_float_poses(self, cand, Q, Dc):
        n0 = len(self.P)
        new = []
        for (cx, cy, th) in cand:
            p, q, x, y = lc.snap_pose(cx, cy, th, self.t, Q, Dc)
            k = len(self.P)
            i = self.add(p, q, x, y, 0.0)
            if i == k:
                new.append(i)
        return new

    def finish(self):
        self.F = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2])),
                            s[3] / s[5], s[4] / s[5]] for s in self.sq], dtype=float)
        self.ncol = len(self.col_pose)
        self.n = len(self.P)
        self.col_pose_a = np.array(self.col_pose, dtype=np.int64)
        self.col_reg_a = np.array(self.col_reg, dtype=np.int64)

    # ---------------------------------------------------------------- float geometry
    def contains_rows(self, pts, pose_lo=0, col_lo=0):
        """coverage-row coefficients for exact points [(X, Y, D)]: R x ncol sparse 0/1.
        Float test with the exact integer test on |margin| <= FTOL.

        `pose_lo`/`col_lo` (both 0 = the whole pose set, the only thing the row/stage code uses)
        restrict the scan to poses `>= pose_lo` and return the `R x (ncol - col_lo)` BLOCK of those
        columns: poses are only ever appended and never move, so after a pricing pass the old
        block of `A` is unchanged and only the new one has to be computed (`--lattice-every`)."""
        R = len(pts)
        rows, cols = [], []
        nex = 0
        F = self.F[pose_lo:]
        for lo in range(0, R, 256):
            chunk = pts[lo:lo + 256]
            x = np.array([X / D for (X, Y, D) in chunk])
            y = np.array([Y / D for (X, Y, D) in chunk])
            dx = x[:, None] - F[None, :, 0]
            dy = y[:, None] - F[None, :, 1]
            u = np.abs(dx * F[None, :, 2] + dy * F[None, :, 3])
            v = np.abs(-dx * F[None, :, 3] + dy * F[None, :, 2])
            m = 0.5 - np.maximum(u, v)
            inside = m > FTOL
            border = np.abs(m) <= FTOL
            for (ri, j) in zip(*np.nonzero(border)):
                X, Y, D = chunk[ri]
                nex += 1
                if lc.sq_contains(self.sq[pose_lo + j], X, Y, D):
                    inside[ri, j] = True
            ri, j = np.nonzero(inside)
            for a, b in zip(ri, j):
                for c in self.cols_of[pose_lo + b]:
                    rows.append(lo + a)
                    cols.append(c - col_lo)
        A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(R, self.ncol - col_lo))
        return A, nex

    def meet_margin(self, i, J):
        """float closed-intersection margin of square i against squares J (array): >= 0 iff they
        meet (up to float error; FTOL band is decided exactly by the caller)"""
        F = self.F
        ci, si = F[i, 2], F[i, 3]
        cj, sj = F[J, 2], F[J, 3]
        ddx = F[J, 0] - F[i, 0]
        ddy = F[J, 1] - F[i, 1]
        cr = ci * cj + si * sj            # cos(theta_j - theta_i)
        sr = ci * sj - si * cj            # sin(theta_j - theta_i)
        w = 0.5 * (np.abs(cr) + np.abs(sr))
        m1 = 0.5 + w - np.abs(ddx * ci + ddy * si)
        m2 = 0.5 + w - np.abs(-ddx * si + ddy * ci)
        m3 = 0.5 + w - np.abs(ddx * cj + ddy * sj)
        m4 = 0.5 + w - np.abs(-ddx * sj + ddy * cj)
        return np.minimum(np.minimum(m1, m2), np.minimum(m3, m4))

    def meets_all(self, i, J):
        """boolean: does square i closed-meet each of J?  (float, exact on the band)"""
        J = np.asarray(J, dtype=np.int64)
        if not len(J):
            return np.zeros(0, bool)
        m = self.meet_margin(i, J)
        out = m > FTOL
        for k in np.nonzero(np.abs(m) <= FTOL)[0]:
            out[k] = lc.sq_meets_sq(self.sq[i], self.sq[int(J[k])])
        return out

    def extend_clique(self, members, mu, order_by=None):
        """greedy maximal extension of a clique (pose indices) over ALL poses.  A pose joins iff
        it closed-meets every CURRENT member (the ones already added included), so the result is
        a clique at every step.  Returns the extended member list (original members first)."""
        members = list(members)
        mem = set(members)
        cand = np.array([j for j in range(self.n) if j not in mem], dtype=np.int64)
        # prefilter: centre within sqrt(2) of every member's centre
        for k in members:
            d2 = (self.F[cand, 0] - self.F[k, 0]) ** 2 + (self.F[cand, 1] - self.F[k, 1]) ** 2
            cand = cand[d2 <= 2.0 + 1e-6]
            if not len(cand):
                break
        for k in members:
            if not len(cand):
                break
            cand = cand[self.meets_all(k, cand)]
        if not len(cand):
            return members
        cxm = self.F[members, 0].mean()
        cym = self.F[members, 1].mean()
        d2 = (self.F[cand, 0] - cxm) ** 2 + (self.F[cand, 1] - cym) ** 2
        if order_by is None:
            order = np.lexsort((d2, -mu[cand]))
        else:
            order = np.argsort(-order_by[cand], kind='stable')
        added = []
        for j in cand[order]:
            j = int(j)
            if added and not self.meets_all(j, np.array(added)).all():
                continue
            added.append(j)
        return members + added

    def joins_all(self, cand, members):
        """the subset of `cand` (pose indices, ORDER PRESERVED) that closed-meets every member.

        Same predicate as `[j for j in cand if self.meets_all(j, members).all()]` and, since the
        closed-intersection margin is symmetric in the two squares, the same float comparisons --
        but vectorised over the candidates with an early exit, so a clique row is extended in a
        handful of numpy calls instead of one per candidate.  (2,645 clique rows x 1,200 new poses
        is 3M calls the other way round: minutes per pricing pass.)"""
        cand = np.asarray(cand, dtype=np.int64)
        for k in members:
            if not len(cand):
                break
            cand = cand[self.meets_all(int(k), cand)]
        return cand

    def verify_clique(self, members, exact=False):
        """pairwise closed intersection of every unordered pair; float with exact fallback, or
        fully exact.  Returns the number of non-intersecting pairs (0 = a clique)."""
        bad = 0
        M = np.array(members, dtype=np.int64)
        for a in range(len(M)):
            J = M[a + 1:]
            if not len(J):
                continue
            if exact:
                for j in J:
                    if not lc.sq_meets_sq(self.sq[int(M[a])], self.sq[int(j)]):
                        bad += 1
            else:
                bad += int((~self.meets_all(int(M[a]), J)).sum())
        return bad


# ====================================================================== exact max clique through j
def max_clique_through_members(nbr, w, j, time_limit, lb=0):
    """as `leaf_ceiling.max_clique_through` but returns (weight, members, complete)"""
    nb = sorted(nbr[j])
    pos = {v: i for i, v in enumerate(nb)}
    sub = [set() for _ in nb]
    for a, v in enumerate(nb):
        for u in nbr[v]:
            if u in pos:
                sub[a].add(pos[u])
    ww = [w[v] for v in nb]
    kw, km, comp = lc.max_clique_int(sub, ww, time_limit, lb=max(0, lb - w[j]))
    return w[j] + kw, [j] + [nb[i] for i in km], comp


# ====================================================================== the loop
class Lever:
    def __init__(self, ps, args, log):
        self.ps = ps
        self.a = args
        self.log = log
        self.t = ps.t
        self.DM = 10 ** 9
        self.target = lc.region_targets(self.t, args)     # {label: Fraction} or None
        self.rows = []            # exact points (X, Y, D)
        self.rowkey = set()
        self.A = sp.csr_matrix((0, ps.ncol))
        self.cliques = []         # list of dict(members=[pose idx], sup=[...], mass0, ext, it)
        self.ckeys = set()
        self.K = sp.csr_matrix((0, ps.ncol))
        self.cage = np.zeros(0, dtype=int)
        self.z = np.zeros(0)
        self.chord = bool(args.chord)
        from t4leaf import Hi
        self.hi = Hi(args.threads, log, 'simplex', args.lp_tlim)
        self.hi_nrows = 0
        self.hi_cq = []
        self.nE = 0
        self.lam = {}
        self.y = np.zeros(0)
        self.chord_dual = np.zeros(4)
        self.found = []           # every violated support clique found, with its anatomy
        # ---- rank-family (odd-polygon) rows, right-hand side (k-1)/2; see search/rankfamily.py
        self.pgons = []           # dict(members, rhs, pts, k, ...)
        self.pgkeys = set()
        self.PG = sp.csr_matrix((0, ps.ncol))
        self.pgz = np.zeros(0)
        self.hi_pg = []
        self.pent_seeds = {}
        self.pent_rng = random.Random(getattr(args, 'pent_seed', 20260911))
        # ---- restricted master (--master, search/CLMASTER.md).  All OFF by default.
        self.master = bool(getattr(args, 'master', False))
        self.act = None           # bool over the COLUMNS: the active (master) set
        self.colage = None        # consecutive zero-mass solves per active column
        self.rc = None            # reduced cost of every loaded column under the last duals
        self.rc_npos = 0          # inactive loaded columns with rc > price_tol
        self.rc_max = 0.0
        self.mcols = None         # column index array currently in the HiGHS model (None = all)
        self.master_n = 0
        self.master_grew = 0
        self.master_passes_used = 0
        self.force_full_master = False   # set once the separation is done: close the master exactly
        self._Acsc = None         # column-sliceable cache of self.A
        self._Kcsc = None         # column-sliceable cache of self.K
        self._PGcsc = None        # column-sliceable cache of self.PG
        self._chordM = None

    # ---------------------------------------------------------------- fixed rows
    def eq_matrix(self):
        rows, vals = [], []
        if self.target:
            reg = self.ps.col_reg_a
            for lab, kv in sorted(self.target.items()):
                li = LABS.index(lab)
                rows.append((reg == li).astype(float))
                vals.append(float(kv))
        if not rows:
            return None, np.zeros(0)
        return sp.csr_matrix(np.array(rows)), np.array(vals)

    def chord_matrix_cached(self):
        if self._chordM is None or self._chordM.shape[1] != self.ps.ncol:
            self._chordM = self.chord_matrix()
        return self._chordM

    def chord_matrix(self):
        ps = self.ps
        t, r = ps.t, ps.r
        m = np.zeros((4, ps.ncol))
        for c, i in enumerate(ps.col_pose):
            p, q, cx, cy = ps.P[i]
            if cy <= r:
                m[0, c] = 1
            if cy >= t - r:
                m[1, c] = 1
            if cx <= r:
                m[2, c] = 1
            if cx >= t - r:
                m[3, c] = 1
        return sp.csr_matrix(m)

    def add_rows(self, pts):
        new = [p for p in pts if p not in self.rowkey]
        if not new:
            return 0
        for p in new:
            self.rowkey.add(p)
        An, nex = self.ps.contains_rows(new)
        self.A = sp.vstack([self.A, An], format='csr') if self.A.shape[0] else An
        self._Acsc = None
        self.rows += new
        return len(new)

    def add_grid_rows(self, pitch):
        D = 1000
        step = int(round(pitch * D))
        pts = []
        tn = int(self.t * D)
        for X in range(step // 2, tn, step):
            for Y in range(step // 2, tn, step):
                pts.append(lc.reduce3(X, Y, D))
        return self.add_rows(pts)

    def resume_rows(self, path):
        """exact rows `X Y D` from a `cl_*_rows.txt` checkpoint"""
        pts = []
        for line in open(path):
            q = line.split()
            if len(q) == 3 and not q[0].startswith('#'):
                pts.append((int(q[0]), int(q[1]), int(q[2])))
        return self.add_rows(pts)

    def resume_cliques(self, path):
        """clique rows from a `cl_*_cliques.txt` checkpoint: pose indices into the pose file the
        checkpoint was written with (which must be the ONLY pose source of this run, in order).
        Every row is re-verified pairwise before it enters the LP."""
        n = 0
        bad = 0
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
            if self.add_clique(mem, dict(it=-1, mass0=None, size0=None, size=len(mem),
                                         complete=None, sup_members=None)):
                n += 1
        if bad:
            self.log(f'   [resume: {bad} checkpoint cliques FAILED the pairwise check and were dropped]')
        return n

    def add_row_file(self, path, D=10 ** 6):
        pts = []
        for line in open(path):
            q = line.split()
            if len(q) >= 2 and not q[0].startswith('#'):
                X = int(round(float(q[0]) * D))
                Y = int(round(float(q[1]) * D))
                X = min(max(X, 0), int(self.t * D))
                Y = min(max(Y, 0), int(self.t * D))
                pts.append(lc.reduce3(X, Y, D))
        return self.add_rows(pts)

    # ---------------------------------------------------------------- cliques
    def add_clique(self, members, info):
        key = frozenset(members)
        if key in self.ckeys:
            return False
        self.ckeys.add(key)
        cols = [c for i in members for c in self.ps.cols_of[i]]
        row = sp.csr_matrix((np.ones(len(cols)), (np.zeros(len(cols), dtype=int), cols)),
                            shape=(1, self.ps.ncol))
        self.K = sp.vstack([self.K, row], format='csr') if self.K.shape[0] else row
        info['members'] = list(members)
        self.cliques.append(info)
        self.cage = np.concatenate([self.cage, [0]])
        self._Kcsc = None
        return True

    # ---------------------------------------------------------------- rank-family rows
    def add_pgon(self, members, rhs, info, pts):
        """an odd-polygon row `mu(members) <= rhs`, rhs = (k-1)/2 (see search/rankfamily.py).
        Same column expansion as a clique row; never aged out."""
        key = frozenset(members)
        if key in self.pgkeys:
            return False
        self.pgkeys.add(key)
        cols = [c for i in members for c in self.ps.cols_of[i]]
        row = sp.csr_matrix((np.ones(len(cols)), (np.zeros(len(cols), dtype=int), cols)),
                            shape=(1, self.ps.ncol))
        self.PG = sp.vstack([self.PG, row], format='csr') if self.PG.shape[0] else row
        self.pgons.append(dict(info, members=list(members), rhs=float(rhs), pts=list(pts)))
        return True

    def sift_cliques(self, age, keep_min=0):
        if not self.cliques or age <= 0:
            return 0
        keep = self.cage < age
        if keep.all() or keep.sum() < keep_min:
            return 0
        idx = np.nonzero(keep)[0]
        self.cliques = [self.cliques[i] for i in idx]
        self.ckeys = set(frozenset(c['members']) for c in self.cliques)
        self.K = self.K[idx]
        self._Kcsc = None
        self.cage = self.cage[idx]
        self.z = self.z[idx] if len(self.z) == len(keep) else np.zeros(len(idx))
        return int((~keep).sum())

    # ---------------------------------------------------------------- LP
    def rebuild(self, cols=None):
        """build the HiGHS model from scratch.  `cols` = None is the full model (the default path,
        unchanged); an index array restricts it to those columns (the restricted master)."""
        ps = self.ps
        INF = self.hi.INF
        self.mcols = cols
        ncol = ps.ncol if cols is None else len(cols)
        self.hi.reset(ncol)

        def sub(M):
            return M if cols is None else M[:, cols]

        E, ev = self.eq_matrix()
        self.nE = 0
        if E is not None:
            self.hi.append([('E', i) for i in range(E.shape[0])], ev, ev, sub(E))
            self.nE = E.shape[0]
        if self.chord:
            self.hi.append([('D', i) for i in range(4)], np.full(4, -INF), np.full(4, 3.0),
                           sub(self.chord_matrix_cached() if cols is not None else self.chord_matrix()))
        if len(self.rows):
            if self._Acsc is None or self._Acsc.shape[1] != ps.ncol:
                self._Acsc = self.A.tocsc() if cols is not None else None
            A = self.A if cols is None else self._Acsc[:, cols].tocsr()
            self.hi.append([('P', i) for i in range(len(self.rows))],
                           np.full(len(self.rows), -INF), np.ones(len(self.rows)), A)
        self.hi_nrows = len(self.rows)
        self.hi_cq = []
        self.hi_pg = []
        if cols is not None and len(self.cliques):
            if self._Kcsc is None or self._Kcsc.shape != self.K.shape or self._Kcsc.nnz != self.K.nnz:
                self._Kcsc = self.K.tocsc()
            K = self._Kcsc[:, cols].tocsr()
            cur = [frozenset(c['members']) for c in self.cliques]
            self.hi.append([('C', k) for k in cur], np.full(len(cur), -INF), np.ones(len(cur)), K)
            self.hi_cq = list(cur)
        if len(self.pgons):                       # rank-family rows (search/rankfamily.py)
            if self._PGcsc is None or self._PGcsc.shape != self.PG.shape \
                    or self._PGcsc.nnz != self.PG.nnz:
                self._PGcsc = self.PG.tocsc()
            PG = self.PG if cols is None else self._PGcsc[:, cols].tocsr()
            gi = list(range(len(self.pgons)))
            self.hi.append([('G', i) for i in gi], np.full(len(gi), -INF),
                           np.array([self.pgons[i]['rhs'] for i in gi]), PG)
            self.hi_pg = gi

    def solve(self):
        if self.master:
            return self.solve_master()
        return self.solve_full()

    # ---------------------------------------------------------------- the restricted master
    @staticmethod
    def _grow(arr, n, dtype, fill):
        if arr is None:
            return np.full(n, fill, dtype=dtype)
        if len(arr) == n:
            return arr
        out = np.full(n, fill, dtype=dtype)
        out[:len(arr)] = arr
        return out

    def cover_regions(self, nper=64):
        """the master is infeasible iff a pinned region cannot reach its count on the active
        columns; the cheap necessary condition is that it has NO active column at all.  Cover such
        a region greedily with its best-priced columns (`_solve_restricted` in search/branch.py)."""
        if not self.target:
            return 0
        reg = self.ps.col_reg_a
        n = 0
        for lab, kv in sorted(self.target.items()):
            if float(kv) <= 0:
                continue
            idx = np.nonzero(reg == LABS.index(lab))[0]
            if not len(idx) or self.act[idx].any():
                continue
            sc = self.rc[idx] if (self.rc is not None and len(self.rc) == len(reg)) else np.zeros(len(idx))
            take = idx[np.argsort(-sc)[:nper]]
            self.act[take] = True
            self.colage[take] = 0
            n += len(take)
        return n

    def price_columns(self):
        """reduced cost of EVERY loaded column under the master's duals,

            rc_j = 1 - sum_r y_r A[r,j] - sum_k z_k K[k,j] - lam(region j) - sum_d chd_d D[d,j],

        i.e. exactly the quantity `price()` evaluates on new lattice candidates, but on the loaded
        columns and with their own region label and their own clique memberships.  The LP is a
        MAXIMISATION, so `rc_j > 0` is the improving direction.  Two sparse mat-vecs."""
        ps = self.ps
        rc = np.ones(ps.ncol)
        if self.A.shape[0] and len(self.y) == self.A.shape[0]:
            rc -= self.A.T.dot(self.y)
        if self.K.shape[0] and len(self.z) == self.K.shape[0]:
            rc -= self.K.T.dot(self.z)
        if self.PG.shape[0] and len(self.pgz) == self.PG.shape[0]:
            rc -= self.PG.T.dot(self.pgz)          # rank-family rows (search/rankfamily.py)
        if self.lam:
            lamv = np.zeros(len(LABS))
            for li, v in self.lam.items():
                lamv[li] = v
            rc -= lamv[ps.col_reg_a]
        if self.chord and np.any(self.chord_dual):
            rc -= np.asarray(self.chord_matrix_cached().T.dot(self.chord_dual)).ravel()
        return rc

    def solve_master(self):
        """the restricted master.  Active columns = the current support, the columns that entered
        this iteration, and the `--master-add` loaded columns with the most positive reduced cost;
        all rows kept; ipm+crossover on the small model.  After each solve EVERY loaded column is
        priced by its reduced cost (two sparse mat-vecs) and the improving ones are added, for up
        to `--master-passes` passes (0 = always to convergence).  The cap is LIFTED on any
        iteration whose separation found nothing (`force_full_master`), so the loop can only stop
        at the exact full-LP optimum of the loaded column set, while the iterations that still
        have rows to separate pay for at most `--master-passes` small ipm solves instead of the
        8-12 the dual chase needs to close a capped master exactly.
        Columns at zero mass for `--master-trim` consecutive solves and with NON-positive reduced
        cost leave the master; they stay LOADED and are priced again every solve.  (Trimming on
        mass alone cycles: the LP is dual-degenerate, so a zero-mass column can carry rc = +0.4 and
        must stay in -- measured on `lc_A0101.txt`, which never terminated.)"""
        ps = self.ps
        a = self.a
        if self.act is None:
            # the first master is the SUPPORT of the loaded masses (a checkpoint carries a few
            # hundred positive-mass poses among its 8k), which is what "current support" means at
            # iteration 0; pricing brings in whatever else is needed.  No masses -> everything.
            m0 = np.array(ps.mu0, dtype=float)
            seed = (m0[ps.col_pose_a] > 1e-12) if m0.size else np.zeros(ps.ncol, bool)
            self.act = seed if seed.any() else np.ones(ps.ncol, bool)
        self.act = self._grow(self.act, ps.ncol, bool, True)
        self.colage = self._grow(self.colage, ps.ncol, np.int64, 0)
        cap = 200 if (a.master_passes <= 0 or self.force_full_master) else a.master_passes
        grew = 0
        passes = 0
        lp_secs = 0.0
        build_secs = 0.0
        res = None
        cols = None
        x = None
        while True:
            if self.rc is not None and len(self.rc) == ps.ncol and a.master_add:
                cand = np.nonzero((~self.act) & (self.rc > a.price_tol))[0]
                if len(cand):
                    # closing the master exactly (nothing left to separate): take EVERY improving
                    # column at once.  One ipm solve at n columns costs what n/2000 of them cost
                    # (t ~ n^1.5 with a small constant: 3 s at 750 columns, 23 s at 2,149, 260 s at
                    # 9,669 on the B40KL2 resume), so wide jumps beat many narrow ones.
                    lim = len(cand) if self.force_full_master else a.master_add
                    pick = cand[np.argsort(-self.rc[cand])[:lim]]
                    self.act[pick] = True
                    self.colage[pick] = 0
                    grew += len(pick)
            grew += self.cover_regions()
            res = None
            for attempt in range(3):
                cols = np.nonzero(self.act)[0]
                if not len(cols):
                    self.act[:] = True
                    continue
                tb = time.time()
                self.rebuild(cols)
                t0 = time.time()
                res = self.hi.run()
                lp_secs += time.time() - t0
                build_secs += t0 - tb
                if os.environ.get('CLDBG'):
                    self.log(f'      [master pass {passes + 1}] {len(cols)} cols: '
                             f'build {t0 - tb:.1f}s solve {time.time() - t0:.1f}s')
                if res is not None:
                    break
                if attempt == 0:
                    self.log('   [master: LP failed; activating every pinned-region column '
                             'and the 5000 best-priced]')
                    if self.target:
                        reg = ps.col_reg_a
                        for lab in self.target:
                            self.act[reg == LABS.index(lab)] = True
                    if self.rc is not None and len(self.rc) == ps.ncol:
                        self.act[np.argsort(-self.rc)[:5000]] = True
                else:
                    self.log('   [master: LP failed twice; falling back to the FULL column set]')
                    self.act[:] = True
            if res is None:
                self.lp_secs = lp_secs
                return None
            passes += 1
            xs, d, obj = res
            x = np.zeros(ps.ncol)
            x[cols] = xs
            self.read_duals(d)
            self.rc = self.price_columns()
            inact = ~self.act
            self.rc_npos = int((inact & (self.rc > a.price_tol)).sum())
            self.rc_max = float(self.rc[inact].max()) if inact.any() else 0.0
            if self.rc_npos == 0 or passes >= cap:
                break
        self.lp_secs = lp_secs + build_secs
        self.master_build_secs = build_secs
        self.master_passes_used = passes
        self.master_grew = grew
        zero = x <= 1e-12
        self.colage[self.act & zero] += 1
        self.colage[~zero] = 0
        if a.master_trim > 0:
            self.act &= ~(zero & (self.colage > a.master_trim) & (self.rc <= a.price_tol))
        self.master_n = int(len(cols))
        return x, obj

    def read_duals(self, d):
        """scatter the HiGHS row duals into y (coverage), z (cliques), lam (regions), chord"""
        cur = [frozenset(c['members']) for c in self.cliques]
        y = np.zeros(len(self.rows))
        z = np.zeros(len(self.cliques))
        pgz = np.zeros(len(self.pgons))
        lam = np.zeros(self.nE)
        chd = np.zeros(4)
        ci = {k: i for i, k in enumerate(cur)}
        for mi, k in enumerate(self.hi.rk):
            v = -d[mi]
            if k[0] == 'P':
                y[k[1]] = max(v, 0.0)
            elif k[0] == 'E':
                lam[k[1]] = v
            elif k[0] == 'D':
                chd[k[1]] = max(v, 0.0)
            elif k[0] == 'C':
                z[ci[k[1]]] = max(v, 0.0)
            elif k[0] == 'G':
                pgz[k[1]] = max(v, 0.0)
        self.pgz = pgz
        self.z = z
        self.cage = np.where(self.z > 1e-9, 0, self.cage + 1)
        self.y = y
        self.chord_dual = chd
        self.lam = {}
        if self.target:
            for i, (lab, kv) in enumerate(sorted(self.target.items())):
                self.lam[LABS.index(lab)] = float(lam[i])

    # ---------------------------------------------------------------- lattice pricing in the loop
    def lattice_price(self, x, mu):
        """`--lattice-every`: one pass of the stage pricer (the `price_pitch`/`price_dth` lattice
        with clique charging) INSIDE the loop -- the new poses join the LOADED column set, the
        coverage rows get their new block and every clique row is extended to the new poses that
        closed-meet all of its members, exactly as the stage code does at a stage boundary."""
        ps = self.ps
        a = self.a
        n0, c0 = ps.n, ps.ncol
        # inside the loop the cliques separated after the last solve have no dual yet (the stage
        # code only ever priced a CONVERGED iterate, where the two always matched): give them 0
        if len(self.z) < len(self.cliques):
            self.z = np.concatenate([self.z, np.zeros(len(self.cliques) - len(self.z))])
        if len(self.pgz) < len(self.pgons):          # ditto for the rank-family rows
            self.pgz = np.concatenate([self.pgz, np.zeros(len(self.pgons) - len(self.pgz))])
        new, gap = price(self, x, mu, a, self.log)
        idx = ps.add_float_poses(new, a.Q, a.Dc)
        ps.finish()
        if not idx:
            return 0, gap, 0, 0
        An, _ = ps.contains_rows(self.rows, pose_lo=n0, col_lo=c0)
        self.A = sp.hstack([self.A, An], format='csr') if self.A.shape[0] else An
        self._Acsc = None
        self._chordM = None
        K = []
        nmem = 0
        for c in self.cliques:
            mem = c['members']
            add = ps.joins_all(idx, mem)
            kept = []
            for j in add:
                j = int(j)
                if not kept or ps.meets_all(j, np.array(kept)).all():
                    kept.append(j)
            c['members'] = mem + kept
            nmem += len(kept)
            cols = [cc for i in c['members'] for cc in ps.cols_of[i]]
            K.append(sp.csr_matrix((np.ones(len(cols)), (np.zeros(len(cols), dtype=int), cols)),
                                   shape=(1, ps.ncol)))
        self.ckeys = set(frozenset(c['members']) for c in self.cliques)
        self.K = sp.vstack(K, format='csr') if K else sp.csr_matrix((0, ps.ncol))
        self._Kcsc = None
        pmem = rankfamily.regrow(self) if self.pgons else 0   # keep polygon rows maximal too
        self.act = self._grow(self.act, ps.ncol, bool, True)        # new columns enter the master
        self.colage = self._grow(self.colage, ps.ncol, np.int64, 0)
        self.rc = None
        self.hi.h = None
        return len(idx), gap, nmem, pmem

    def solve_full(self):
        ps = self.ps
        INF = self.hi.INF
        if self.hi.h is None or self.hi.ncols != ps.ncol:
            self.rebuild()
        if len(self.rows) > self.hi_nrows:
            k0 = self.hi_nrows
            self.hi.append([('P', i) for i in range(k0, len(self.rows))],
                           np.full(len(self.rows) - k0, -INF), np.ones(len(self.rows) - k0),
                           self.A[k0:])
            self.hi_nrows = len(self.rows)
        cur = [frozenset(c['members']) for c in self.cliques]
        curs = set(cur)
        gone = [k for k in self.hi_cq if k not in curs]
        if gone:
            self.hi.delete([('C', k) for k in gone])
            self.hi_cq = [k for k in self.hi_cq if k in curs]
        have = set(self.hi_cq)
        addi = [i for i, k in enumerate(cur) if k not in have]
        if addi:
            self.hi.append([('C', cur[i]) for i in addi], np.full(len(addi), -INF),
                           np.ones(len(addi)), self.K[addi])
            self.hi_cq += [cur[i] for i in addi]
        if len(self.pgons) > len(self.hi_pg):                 # rank-family rows, never removed
            g0 = len(self.hi_pg)
            gi = list(range(g0, len(self.pgons)))
            self.hi.append([('G', i) for i in gi], np.full(len(gi), -INF),
                           np.array([self.pgons[i]['rhs'] for i in gi]), self.PG[g0:])
            self.hi_pg = list(range(len(self.pgons)))
        t0 = time.time()
        res = self.hi.run()
        self.lp_secs = time.time() - t0
        if res is None:
            return None
        x, d, obj = res
        self.read_duals(d)
        return x, obj

    # ---------------------------------------------------------------- per-pose masses
    def pose_mass(self, x):
        mu = np.zeros(self.ps.n)
        np.add.at(mu, self.ps.col_pose_a, x)
        return mu

    def region_split(self, x):
        out = np.zeros(13)
        np.add.at(out, self.ps.col_reg_a, x)
        return out

    # ---------------------------------------------------------------- certification of one solve
    def certify(self, mu, procs):
        """exact arrangement of the support; returns (sup, V, Inc, M, cov, w_int)"""
        sup = np.nonzero(mu > 1e-12)[0]
        squares = [self.ps.sq[i][:8] + (int(round(mu[i] * self.DM)),) for i in sup]
        V = lc.enumerate_vertices(squares, self.t, verbose=False, procs=procs)
        Inc = lc.incidences(squares, V, verbose=False, procs=procs)
        w = np.array([s[8] for s in squares], dtype=np.int64)
        cov = Inc.cov(w)
        M = Fr(int(cov.max()), self.DM) if len(cov) else Fr(0)
        return sup, squares, V, Inc, float(M), cov, w

    def separate(self, sup, squares, w, mu, it, want, time_limit, top, ktol, Inc=None):
        """exact max-weight cliques on the support; returns (kmax, all_complete, new_rows, best)"""
        DM = self.DM
        nbr = lc.closed_graph(squares, verbose=False)
        wl = [int(v) for v in w]
        lb = int(DM * (1.0 + ktol))
        found = {}
        allc = True
        t0 = time.time()
        Kw, Km, comp = lc.max_clique_int(nbr, wl, time_limit, lb=0)
        allc &= comp
        kmax = Kw / DM
        best = None
        if Kw > lb:
            found[frozenset(Km)] = (Kw, comp)
        # through the heaviest support squares (and the members of the best clique)
        heavy = list(np.argsort(-w)[:top])
        seen = set()
        for j in heavy + list(Km):
            j = int(j)
            if j in seen or time.time() - t0 > 4 * time_limit:
                continue
            seen.add(j)
            Kw2, Km2, comp2 = max_clique_through_members(nbr, wl, j, time_limit / 4, lb=lb)
            allc &= comp2
            if Kw2 > lb:
                found[frozenset(Km2)] = (Kw2, comp2)
                kmax = max(kmax, Kw2 / DM)
        sep_secs = time.time() - t0
        # extension over all poses, verification, and the rows
        cl = sorted(found.items(), key=lambda kv: -kv[1][0])[:want]
        nnew = 0
        for key, (Kw2, comp2) in cl:
            mem_sup = sorted(key, key=lambda k: -w[k])
            members = [int(sup[k]) for k in mem_sup]
            ext = self.ps.extend_clique(members, mu) if self.a.extend else list(members)
            bad = self.ps.verify_clique(ext)
            if bad:
                self.log(f'   [clique DROPPED: {bad} non-intersecting pairs after extension]')
                continue
            info = dict(it=it, mass0=Kw2 / DM, size0=len(members), size=len(ext),
                        complete=bool(comp2), sup_members=members)
            if Inc is not None:
                # anatomy-lite at separation time: the heaviest point sub-clique of the members
                # (exact coverage of the member sub-measure over the support's arrangement)
                cv = Inc.cov(w, subset=list(mem_sup))
                ip = int(np.argmax(cv))
                info['point_mass'] = int(cv[ip]) / DM
                info['excess'] = (Kw2 - int(cv[ip])) / DM
                regs = {}
                for k in members:
                    lab = lc.regions_of(self.ps.P[k][2], self.ps.P[k][3], self.ps.boxes)[0]
                    regs[lab] = regs.get(lab, 0.0) + float(mu[k])
                info['regions'] = {k: round(v, 4) for k, v in regs.items()}
                th = np.degrees(np.arctan2(self.ps.F[members, 3], self.ps.F[members, 2])) % 90
                info['nangles'] = len(set(int(round(v / 5)) for v in th))
                info['centre'] = [round(float(self.ps.F[members, 0].mean()), 3),
                                  round(float(self.ps.F[members, 1].mean()), 3)]
            if self.add_clique(ext, info):
                nnew += 1
                self.found.append(dict(info, members=ext))
                if best is None:
                    best = info
        return kmax, allc, nnew, best, sep_secs, len(found)

    # ---------------------------------------------------------------- one stage
    def run_stage(self, tag, x0=None):
        a = self.a
        T0 = time.time()
        hist = []
        rec = None
        for it in range(a.iters):
            ndrop = self.sift_cliques(a.cq_age, keep_min=0)
            sol = self.solve()
            if sol is None:
                self.log(f'   [{tag}.{it}] LP failed')
                return None
            x, val = sol
            mu = self.pose_mass(x)
            sup, squares, V, Inc, M, cov, w = self.certify(mu, a.procs)
            bad = np.nonzero(cov > self.DM * (1.0 + 1e-9))[0]
            nrow = 0
            if len(bad):
                o = bad[np.argsort(-cov[bad])[:a.row_cap]]
                nrow = self.add_rows([V[i] for i in o])
            kmax, allc, ncq, best, sep_secs, nfound = self.separate(
                sup, squares, w, mu, it, a.cq_want, a.clique_time, a.cq_top, a.ktol, Inc=Inc)
            npg, pbest, pg_secs = 0, None, 0.0
            if getattr(a, "pent", None):    # rank-family (odd-polygon) rows, RHS (k-1)/2
                npg, pbest, pg_secs, _pf = rankfamily.separate(
                    self, sup, squares, V, Inc, w, mu, it)
            R = self.region_split(x)
            rec = dict(tag=tag, it=it, LP=float(val), M=M, kmax=float(kmax), complete=bool(allc),
                       rows=len(self.rows), cq=len(self.cliques), nsup=len(sup), nverts=len(V),
                       bad=int(len(bad)), interior=float(R[12]), corners=R[:4].tolist(),
                       slots=R[4:12].tolist(), chord_dual=self.chord_dual.tolist(),
                       lp_secs=self.lp_secs, secs=time.time() - T0,
                       best=(None if best is None else dict(mass0=best['mass0'], size0=best['size0'],
                                                             size=best['size'])))
            if self.master:
                rec.update(master_cols=self.master_n, master_new=self.master_grew,
                           master_passes=self.master_passes_used,
                           rc_npos=self.rc_npos, rc_max=self.rc_max)
            hist.append(rec)
            self.log(f'   [{tag}.{it}] LP={val:.6f} M={M:.9f} kmax={kmax:.6f}{"" if allc else "*"} '
                     f'rows={len(self.rows)}(+{nrow}) cq={len(self.cliques)}(+{ncq}/{nfound},-{ndrop}) '
                     f'sup={len(sup)} verts={len(V)} int={R[12]:.6f} '
                     + (f'master={self.master_n}/{self.ps.ncol}(+{self.master_grew}'
                        f'/{self.master_passes_used}p) rc+={self.rc_npos}/{self.rc_max:+.2e} '
                        if self.master else '')
                     + (f'best={best["mass0"]:.4f}/{best["size0"]}->{best["size"]} ' if best else '')
                     + (f'pg={len(self.pgons)}(+{npg}) ' if getattr(a, 'pent', None) else '')
                     + (f'pbest=k{pbest["k"]}:{pbest["mass0"]:.4f}/{pbest["rhs"]:.1f}'
                        f'(+{pbest["excess"]:.4f})/{pbest["size"]} ' if pbest else '')
                     + f'(lp {self.lp_secs:.1f}s it{self.hi.iters}, sep {sep_secs:.1f}s'
                     + (f', pg {pg_secs:.1f}s' if getattr(a, 'pent', None) else '')
                     + f', {time.time() - T0:.0f}s)')
            rec['pg'] = len(self.pgons)
            rec['npg'] = npg
            rec['pbest'] = pbest
            geo = (nrow == 0 and ncq == 0 and npg == 0 and M <= 1 + 1e-9 and kmax <= 1 + a.ktol)
            # nothing left to separate: the next master solve runs to full pricing convergence,
            # so the loop can only terminate on the exact full-LP optimum of the loaded columns
            self.force_full_master = bool(geo)
            done = geo and (not self.master or self.rc_npos == 0)
            # --lattice-every: price the full lattice every M iterations, and whenever the row/clique
            # separation has stalled (no new row of either kind), so that a run never spins
            stall = (nrow == 0 and ncq == 0)
            if a.lattice_every > 0 and (stall or it % a.lattice_every == a.lattice_every - 1):
                npos, gap, nmem, pmem = self.lattice_price(x, mu)
                self.log(f'   [{tag}.{it}] lattice: +{npos} poses -> {self.ps.n} poses, '
                         f'{self.ps.ncol} columns (gap {gap:+.6f}, {nmem} new clique memberships'
                         + (f', {pmem} new polygon memberships' if self.pgons else '')
                         + f', {time.time() - T0:.0f}s)')
                rec['lattice'] = dict(new=npos, gap=gap, memberships=nmem, pgon_memberships=pmem,
                                      poses=self.ps.n, cols=self.ps.ncol)
                if npos:
                    # the pose set just grew: pad this iterate onto the new columns (they carry no
                    # mass in it) so that the checkpoint, `region_split` and `finalize` -- all of
                    # which index by the CURRENT ps.ncol / ps.n -- stay in step with it
                    x = np.concatenate([x, np.zeros(self.ps.ncol - len(x))])
                    mu = np.concatenate([mu, np.zeros(self.ps.n - len(mu))])
                if npos:
                    done = False
            rec['converged'] = bool(done and allc)
            if done:
                if not allc:
                    self.log('   [no violated row found, but a B&B was incomplete: NOT a proof]')
                break
            if time.time() - T0 > a.time:
                self.log('   [time limit inside the loop]')
                break
            if a.ckpt and (it % a.ckpt == a.ckpt - 1):
                self.checkpoint(tag, x, mu, rec)
        self.last = (x, mu, rec)
        return hist

    # ---------------------------------------------------------------- output
    def write_measure(self, path, mu_int, note=''):
        ps = self.ps
        with open(path, 'w') as f:
            f.write(f"# cliquelever measure; t = {self.t}; sym 1; {note}\n")
            f.write(f"# mass = {Fr(int(sum(mu_int)), self.DM)} = {sum(mu_int) / self.DM:.12f}\n")
            f.write("# pose p q cx cy mass : theta = 2 arctan(p/q)\n")
            for i, m in enumerate(mu_int):
                if m > 0:
                    p, q, cx, cy = ps.P[i]
                    f.write(f"pose {p} {q} {cx} {cy} {Fr(int(m), self.DM)}\n")

    def checkpoint(self, tag, x, mu, rec):
        ps = self.ps
        base = os.path.join(RUNS, f'cl_{self.a.TAG}')
        mu_int = [int(round(v * self.DM)) for v in mu]
        self.write_measure(base + '_measure.txt', mu_int, note=f'{tag} LP={rec["LP"]:.9f} (float masses rounded to 1e-9)')
        with open(base + '_poses.txt', 'w') as f:
            f.write(f"# cliquelever pose set; t = {self.t}; sym 1; {ps.n} poses, {ps.ncol} columns\n")
            for i, (p, q, cx, cy) in enumerate(ps.P):
                f.write(f"pose {p} {q} {cx} {cy} {Fr(mu_int[i], self.DM)}\n")
        with open(base + '_cliques.txt', 'w') as f:
            f.write(f"# cliquelever cliques: pose indices into cl_{self.a.TAG}_poses.txt; "
                    f"one line per row in the LP\n")
            for c in self.cliques:
                f.write(' '.join(str(i) for i in c['members']) + '\n')
        with open(base + '_rows.txt', 'w') as f:
            f.write("# cliquelever coverage rows: X Y D\n")
            for (X, Y, D) in self.rows:
                f.write(f"{X} {Y} {D}\n")
        # every violated clique found so far, with its anatomy at the time it was separated
        json.dump([{k: v for k, v in c.items() if k != 'members'} for c in self.found],
                  open(base + '_found.json', 'w'), indent=0)
        if self.pgons:
            json.dump([dict({k: v for k, v in c.items() if k not in ('members', 'pts')},
                            members=c['members'],
                            dual=float(self.pgz[i]) if i < len(self.pgz) else 0.0,
                            mass_now=float(sum(mu[j] for j in c['members'])))
                       for i, c in enumerate(self.pgons)],
                      open(base + '_pgons.json', 'w'), indent=1)

    # ---------------------------------------------------------------- exact final certification
    def finalize(self, x, mu, procs):
        """round down, top up the pinned regions on exact coverage AND clique slack, and certify
        M, the max clique (complete B&B) and the region masses exactly"""
        ps = self.ps
        DM = self.DM
        os.makedirs(RUNS, exist_ok=True)
        # per-column integer masses (the LP's explicit region assignment), rounded DOWN
        xc = [int(math.floor(v * DM)) for v in x]
        res = {}

        def pose_int():
            m = [0] * ps.n
            for c, v in enumerate(xc):
                m[ps.col_pose[c]] += v
            return m

        def arrangement(mi):
            sup = [i for i in range(ps.n) if mi[i] > 0]
            squares = [ps.sq[i][:8] + (mi[i],) for i in sup]
            V = lc.enumerate_vertices(squares, self.t, verbose=False, procs=procs)
            Inc = lc.incidences(squares, V, verbose=False, procs=procs)
            w = np.array([s[8] for s in squares], dtype=np.int64)
            return sup, squares, V, Inc, w, Inc.cov(w)

        mi = pose_int()
        sup, squares, V, Inc, w, cov = arrangement(mi)
        M = Fr(int(cov.max()), DM)
        self.log(f'   final: rounded down: mass {sum(mi) / DM:.9f}, exact M = {float(M):.12f}')
        if M > 1:
            self.log('   final: M > 1 after rounding down (float LP infeasibility); scaling by 1/M')
            xc = [(v * M.denominator) // M.numerator for v in xc]
            mi = pose_int()
            sup, squares, V, Inc, w, cov = arrangement(mi)
            M = Fr(int(cov.max()), DM)
        nbr = lc.closed_graph(squares, verbose=False)
        wl = [int(v) for v in w]
        pos = {i: k for k, i in enumerate(sup)}
        # ---- region top-up on exact slack (coverage, cliques AND the rank-family rows)
        pg = rankfamily.TopUpGuard(self, mi)
        if self.target:
            for lab, kv in sorted(self.target.items()):
                li = LABS.index(lab)
                cols = [c for c in range(ps.ncol) if ps.col_reg[c] == li]
                want = int(kv * DM)
                have = sum(xc[c] for c in cols)
                d = want - have
                if d < 0:
                    self.log(f'   final: region {lab} OVER by {-d / DM:.3e} (should not happen)')
                    continue
                tries = 0
                while d > 0 and tries < 50:
                    tries += 1
                    cand = [c for c in cols if xc[c] > 0 and ps.col_pose[c] in pos]
                    bestc, bests = None, 0
                    for c in cand:
                        i = ps.col_pose[c]
                        k = pos[i]
                        cs = int(DM - cov[Inc.cols[k]].max()) if len(Inc.cols[k]) else DM
                        if cs <= bests:
                            continue
                        Kw, _, comp = max_clique_through_members(nbr, wl, k, 30.0, lb=0)
                        ks = DM - Kw if comp else 0
                        s = min(cs, ks, pg.slack(i))
                        if s > bests:
                            bestc, bests = c, s
                        if bests >= d:
                            break
                    if bestc is None or bests <= 0:
                        break
                    add = min(d, bests)
                    xc[bestc] += add
                    d -= add
                    i = ps.col_pose[bestc]
                    pg.credit(i, add)
                    k = pos[i]
                    wl[k] += add
                    w[k] += add
                    cov[Inc.cols[k]] += add
                self.log(f'   final: region {lab} topped up to {Fr(sum(xc[c] for c in cols), DM)} '
                         f'(target {kv}, residual {d / DM:+.3e})')
        mi = pose_int()
        sup, squares, V, Inc, w, cov = arrangement(mi)
        M = Fr(int(cov.max()), DM)
        nbr = lc.closed_graph(squares, verbose=False)
        Kw, Km, comp = lc.max_clique_int(nbr, [int(v) for v in w], self.a.clique_time * 4, lb=0)
        Kq = Fr(Kw, DM)
        mass = Fr(sum(mi), DM)
        meta = [(Fr(s[0], s[2]), Fr(s[1], s[2]), s[8]) for s in squares]
        ok, assigned, amb, note = lc.region_report(meta, ps.boxes, self.target, DM)
        pgok, pgworst, pgrep = rankfamily.check_final(self, mi)
        if self.pgons:
            tight = [r for r in pgrep if r['slack'] <= 1e-9]
            self.log(f'   FINAL EXACT: {len(self.pgons)} rank-family rows: '
                     f'{"OK" if pgok else "FAIL"} (worst violation {pgworst:+.9f}, '
                     f'{len(tight)} tight, all re-derived from their exact anchors: '
                     f'{all(r["rederived"] for r in pgrep)})')
        self.log(f'   FINAL EXACT: mass = {mass} = {float(mass):.9f}; M = {M} = {float(M):.12f} '
                 f'{"OK" if M <= 1 else "FAIL"}; max clique = {Kq} = {float(Kq):.9f} '
                 f'(size {len(Km)}, complete {comp}) {"OK" if Kq <= 1 else "FAIL"}; '
                 f'regions {"OK" if ok else "FAIL"} ({note}); {len(sup)} poses')
        res = dict(mass=str(mass), mass_float=float(mass), M=str(M), M_float=float(M),
                   clique=str(Kq), clique_float=float(Kq), clique_complete=bool(comp),
                   clique_size=len(Km), regions_ok=bool(ok), poses=len(sup), vertices=len(V),
                   regions={lab: str(Fr(assigned.get(lab, 0), DM)) for lab in LABS},
                   pgons=len(self.pgons), pgons_ok=bool(pgok), pgons_worst=pgworst,
                   pgon_rows=pgrep,
                   certified=bool(M <= 1 and Kq <= 1 and comp and ok and pgok))
        path = os.path.join(RUNS, f'cl_{self.a.TAG}_exact.txt')
        self.write_measure(path, mi, note=f'exactly certified: M = {M}, max clique = {Kq} (complete {comp})')
        res['file'] = path
        return res


# ====================================================================== anatomy of a clique
def anatomy(ps, members, mu, lever):
    """what a violated clique looks like: region split, best common point, angle spread, and how
    much of it the anchor family K(p, A) can see (exact scan + ascent on the sub-measure)"""
    DM = lever.DM
    mem = [i for i in members if mu[i] > 1e-12]
    squares = [ps.sq[i][:8] + (int(round(mu[i] * DM)),) for i in mem]
    w = np.array([s[8] for s in squares], dtype=np.int64)
    tot = int(w.sum())
    V = lc.enumerate_vertices(squares, ps.t, verbose=False, procs=1)
    Inc = lc.incidences(squares, V, verbose=False, procs=1)
    cov = Inc.cov(w)
    ipt = int(np.argmax(cov))
    X, Y, D = V[ipt]
    point = int(cov[ipt])
    regs = {}
    for i in mem:
        for lab in lc.regions_of(ps.P[i][2], ps.P[i][3], ps.boxes)[:1]:
            regs[lab] = regs.get(lab, 0) + float(mu[i])
    th = np.degrees(np.arctan2(ps.F[mem, 3], ps.F[mem, 2])) % 90
    out = dict(size=len(mem), mass=tot / DM, point_mass=point / DM, point=(X / D, Y / D),
               excess_over_point=(tot - point) / DM, regions={k: round(v, 4) for k, v in regs.items()},
               centre=(float(ps.F[mem, 0].mean()), float(ps.F[mem, 1].mean())),
               spread=(float(ps.F[mem, 0].std()), float(ps.F[mem, 1].std())),
               angles=sorted(set(int(round(v / 5) * 5) for v in th)),
               heaviest=[(round(mu[i], 4), round(ps.F[i, 0], 3), round(ps.F[i, 1], 3),
                          round(float(th[k]), 1)) for k, i in
                         sorted(enumerate(mem), key=lambda kv: -mu[kv[1]])[:6]])
    # anchor family on the sub-measure (exact scan + ascent, LEAF_CEILING.md 2.5 `scan`)
    try:
        lc.MCAP[0] = int(cov.max())          # the sound cap of `best_on_line`, for THIS sub-measure
        F = lc.float_squares(squares, DM)
        Vf = np.array([[X / D, Y / D] for (X, Y, D) in V])
        cand = lc.scan_candidates(squares, F, ps.t, Vf, cov.astype(float) / DM, 60, 24, 12, 10,
                                  0.9, covfrac=0.3, verbose=False)
        best = (Fr(point, DM), None)
        for (val, p, P, Q) in cand[:30]:
            P0, P1, Ds = lc.rationalise_segment(P, Q, 10 ** 6)
            if P0 == P1:
                continue
            v = lc.anchor_value(squares, Inc, DM, P0, P1, Ds, w)
            if v[0] > best[0]:
                best = (v[0], (P0, P1, Ds))
        if best[1] is not None:
            bv, bseg, rounds = lc.local_ascent(squares, Inc, DM, w, V, best[1], 20.0, verbose=False)
            if bv[0] > best[0]:
                best = (bv[0], bseg)
        out['anchor_best'] = float(best[0])
    except Exception as e:                                # noqa: BLE001
        out['anchor_best'] = None
        out['anchor_err'] = repr(e)
    return out


# ====================================================================== pricing (stage 3)
def price(lever, x, mu, a, log):
    """one stage of `t4leaf.py`'s stratified pricer: rc = 1 - capture - lam(region) - clique
    duals, where a candidate is charged a clique's dual iff it closed-meets every member (so
    charging it is legitimate: the row extends to it).  Returns the new pose indices."""
    import packing_dual as pd
    import level2_lp as L2
    from level2_regions import classify
    ps = lever.ps
    t = float(ps.t)
    lib = pd.build_lib()
    y = lever.y
    sel = np.nonzero(y > 1e-9)[0]
    ax = np.ascontiguousarray(np.array([lever.rows[i][0] / lever.rows[i][2] for i in sel]))
    ay = np.ascontiguousarray(np.array([lever.rows[i][1] / lever.rows[i][2] for i in sel]))
    aw = np.ascontiguousarray(y[sel])
    o = np.argsort(ax, kind='stable')
    ax, ay, aw = np.ascontiguousarray(ax[o]), np.ascontiguousarray(ay[o]), np.ascontiguousarray(aw[o])

    zi = [i for i in np.argsort(-lever.z) if lever.z[i] > 1e-9][:500]

    def region_of(p):
        k, j = classify(p[0], p[1], t, float(ps.r))
        return j if k == 'C' else (4 + j if k == 'W' else 12)

    def clique_cost(cand):
        # identical arithmetic to the original mask version, but carried on the SHRINKING index
        # set of surviving candidates instead of the full 170k-long mask: the member loop is
        # 800-deep on a wide clique row and the survivors are a handful, so this is the difference
        # between ~20 min and a few seconds per pricing pass (search/CLMASTER.md)
        out = np.zeros(len(cand))
        if not zi or not len(cand):
            return out
        C = np.asarray(cand, dtype=float).reshape(-1, 3)
        G = np.column_stack([C[:, 0], C[:, 1], np.cos(C[:, 2]), np.sin(C[:, 2])])
        for ci in zi:
            mem = lever.cliques[ci]['members']
            # a member of K must be within sqrt(2) of EVERY member's centre: prefilter by the
            # bounding box of the member centres, expanded by sqrt(2)
            bx0, bx1 = ps.F[mem, 0].min() - 1.4143, ps.F[mem, 0].max() + 1.4143
            by0, by1 = ps.F[mem, 1].min() - 1.4143, ps.F[mem, 1].max() + 1.4143
            sel = np.nonzero((G[:, 0] >= bx0) & (G[:, 0] <= bx1)
                             & (G[:, 1] >= by0) & (G[:, 1] <= by1))[0]
            for k in mem:
                if not len(sel):
                    break
                Fk = ps.F[k]
                g0, g1, g2, g3 = G[sel, 0], G[sel, 1], G[sel, 2], G[sel, 3]
                cr = Fk[2] * g2 + Fk[3] * g3
                sr = Fk[2] * g3 - Fk[3] * g2
                wj = 0.5 * (np.abs(cr) + np.abs(sr))
                ddx = g0 - Fk[0]
                ddy = g1 - Fk[1]
                m1 = 0.5 + wj - np.abs(ddx * Fk[2] + ddy * Fk[3])
                m2 = 0.5 + wj - np.abs(-ddx * Fk[3] + ddy * Fk[2])
                m3 = 0.5 + wj - np.abs(ddx * g2 + ddy * g3)
                m4 = 0.5 + wj - np.abs(-ddx * g3 + ddy * g2)
                sel = sel[np.minimum(np.minimum(m1, m2), np.minimum(m3, m4)) >= 0]
            out[sel] += lever.z[ci]
        return out

    def rc_of(poses):
        cap = pd.capture(lib, ax, ay, aw, poses, a.threads)
        lam = np.array([lever.lam.get(region_of(p), 0.0) for p in poses])
        # the rank-family rows charge a candidate exactly when they extend to it (rankfamily.py)
        return 1.0 - cap - lam - clique_cost(poses) - rankfamily.pgon_cost(lever, poses)

    # rc check on the support COLUMNS (should be ~0 at an optimum): the column's own region label
    # and its actual clique rows, not the geometric proxies used for new candidates
    sc = np.nonzero(x > 1e-9)[0]
    scp = [tuple(map(float, (ps.F[ps.col_pose[c], 0], ps.F[ps.col_pose[c], 1],
                             math.atan2(ps.F[ps.col_pose[c], 3], ps.F[ps.col_pose[c], 2])))) for c in sc]
    capc = pd.capture(lib, ax, ay, aw, scp, a.threads)
    lamc = np.array([lever.lam.get(int(ps.col_reg[c]), 0.0) for c in sc])
    zc = np.asarray(lever.K[:, sc].T @ lever.z).ravel() if len(lever.z) else np.zeros(len(sc))
    gc = (np.asarray(lever.PG[:, sc].T @ lever.pgz).ravel()
          if len(lever.pgz) == lever.PG.shape[0] and lever.PG.shape[0] else np.zeros(len(sc)))
    rcs = 1.0 - capc - lamc - zc - gc
    log(f'   pricing: rc on the {len(sc)} support columns: max |rc| = {np.abs(rcs).max():.2e} '
        f'(mean {rcs.mean():+.2e})')
    supp = [tuple(map(float, (ps.F[i, 0], ps.F[i, 1], math.atan2(ps.F[i, 3], ps.F[i, 2]))))
            for i in np.nonzero(mu > 1e-9)[0]]
    cand = L2.seed_poses_full(t, a.price_pitch, a.price_dth)
    rc = rc_of(cand)
    creg = np.array([region_of(p) for p in cand])
    idx = []
    regs = sorted(set(int(v) for v in creg))
    quota = max(a.cg_want // max(len(regs), 1), 1)
    for rr in regs:
        w = np.nonzero(creg == rr)[0]
        o = w[np.argsort(-rc[w])[:quota]]
        idx += [int(i) for i in o if rc[i] > a.price_tol]
    extra = [int(i) for i in np.argsort(-rc) if rc[i] > a.price_tol and int(i) not in set(idx)]
    idx = (idx + extra)[:a.cg_want]
    best = [cand[i] for i in idx]
    gi = np.nonzero(creg == 12)[0]
    log(f'   pricing: {len(cand)} candidates, best rc {rc.max():+.6f} (interior '
        f'{rc[gi].max() if len(gi) else float("nan"):+.6f}), kept {len(best)} over {len(regs)} regions')
    ref = pd.refine(lib, ax, ay, aw, best, t, a.threads) if best else []
    ref = [(c[0], c[1], c[2]) for c in ref]
    rr = rc_of(ref) if ref else np.zeros(0)
    heavy = [supp[i] for i in np.argsort(-mu[mu > 1e-9])[:100]]
    nb = pd.neighbours(heavy, t)
    rn = rc_of(nb)
    keepnb = [nb[i] for i in np.argsort(-rn)[:a.cg_want] if rn[i] > a.price_tol]
    new = list(best) + [ref[i] for i in range(len(ref)) if rr[i] > a.price_tol] + keepnb
    log(f'   pricing: +{len(ref)} refined (best rc {rr.max() if len(rr) else 0:+.6f}), '
        f'+{len(keepnb)} neighbours of the support (best rc {rn.max():+.6f})')
    return new, float(rc.max())


# ====================================================================== main
def cmd_run(a):
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'cl_{a.TAG}.log'), 'a')

    def log(m):
        print(m, flush=True)
        logf.write(m + '\n')
        logf.flush()

    T0 = time.time()
    ps = Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    for f in a.exact:
        n = ps.add_exact_file(f)
        log(f'   +{n} exact poses from {os.path.basename(f)}')
    for f in a.load_poses:
        n = ps.add_float_file(f, a.Q, a.Dc)
        log(f'   +{n} snapped poses from {os.path.basename(f)} (Q={a.Q}, Dc={a.Dc})')
    ps.finish()
    nb = sum(1 for i in range(ps.n) if len(ps.cols_of[i]) > 1)
    log(f'# cliquelever t={ps.t} r={ps.r} corners={a.corners} patterns={a.patterns} chord={a.chord} '
        f'extend={a.extend}: {ps.n} poses, {ps.ncol} columns ({nb} boundary poses duplicated), '
        f'initial mass {sum(ps.mu0):.6f}; args={vars(a)}')
    lever = Lever(ps, a, log)
    if a.row_pitch > 0:
        log(f'   +{lever.add_grid_rows(a.row_pitch)} grid rows (pitch {a.row_pitch})')
    for f in a.load_rows:
        log(f'   +{lever.add_row_file(f)} rows from {os.path.basename(f)}')
    for f in a.resume_rows:
        log(f'   +{lever.resume_rows(f)} exact rows resumed from {os.path.basename(f)}')
    for f in a.resume_cliques:
        log(f'   +{lever.resume_cliques(f)} clique rows resumed from {os.path.basename(f)} '
            f'(each re-verified pairwise)')
    for f in getattr(a, 'resume_pgons', []) or []:
        nadd, nbad = rankfamily.resume(lever, f)
        log(f'   +{nadd} rank-family rows resumed from {os.path.basename(f)} (membership '
            f're-derived from the exact anchors over the current pose set; {nbad} dropped)')
    # rows at the arrangement vertices of the INITIAL support (from the source masses)
    mu0 = np.array(ps.mu0)
    if mu0.sum() > 0:
        sup, squares, V, Inc, M0, cov, w = lever.certify(mu0, a.procs)
        o = np.argsort(-cov)[:a.rows0]
        n = lever.add_rows([V[i] for i in o])
        log(f'   +{n} rows at the heaviest arrangement vertices of the initial support '
            f'({len(sup)} poses, {len(V)} vertices, initial M = {M0:.6f})')
    log(f'   {len(lever.rows)} rows, {time.time() - T0:.0f}s')

    hist = []
    results = []
    stage = 0
    while True:
        h = lever.run_stage(f's{stage}')
        if h is None:
            break
        hist += h
        x, mu, rec = lever.last
        lever.checkpoint(f's{stage}', x, mu, rec)
        R = lever.region_split(x)
        log(f'   STAGE {stage}: LP={rec["LP"]:.6f} M={rec["M"]:.9f} kmax={rec["kmax"]:.6f} '
            f'complete={rec["complete"]} conv={rec["converged"]} rows={rec["rows"]} cq={rec["cq"]} '
            f'poses={ps.n} int={R[12]:.6f} corners={np.round(R[:4], 6).tolist()} '
            f'slots={np.round(R[4:12], 6).tolist()} chord_dual={np.round(lever.chord_dual, 4).tolist()}')
        fin = lever.finalize(x, mu, a.procs) if a.finalize else None
        srec = dict(stage=stage, LP=rec['LP'], M=rec['M'], kmax=rec['kmax'], complete=rec['complete'],
                    converged=rec['converged'], rows=rec['rows'], cq=rec['cq'], poses=ps.n,
                    cols=ps.ncol, interior=float(R[12]), iters=len(h), secs=time.time() - T0,
                    exact=fin)
        results.append(srec)
        json.dump(dict(t=str(ps.t), tag=a.TAG, args=vars(a), stages=results, hist=hist),
                  open(os.path.join(RUNS, f'cl_{a.TAG}.json'), 'w'), indent=1)
        if stage >= a.price or time.time() - T0 > a.time:
            break
        # ---- one pricing stage
        new, gap = price(lever, x, mu, a, log)
        idx = ps.add_float_poses(new, a.Q, a.Dc)
        ps.finish()
        log(f'   stage {stage}: +{len(idx)} new poses -> {ps.n} poses, {ps.ncol} columns '
            f'(pricing gap {gap:+.6f})')
        if not idx:
            log('   [no improving column]')
            break
        # extend the coverage rows and the clique rows to the new columns; rebuild the model
        A2, _ = ps.contains_rows(lever.rows)
        lever.A = A2
        K = []
        next_ = 0
        for c in lever.cliques:
            mem = c['members']
            add = ps.joins_all(idx, mem)
            # the new members must also meet each other: greedy, as in extend_clique
            kept = []
            for j in add:
                j = int(j)
                if not kept or ps.meets_all(j, np.array(kept)).all():
                    kept.append(j)
            c['members'] = mem + kept
            next_ += len(kept)
            cols = [cc for i in c['members'] for cc in ps.cols_of[i]]
            K.append(sp.csr_matrix((np.ones(len(cols)), (np.zeros(len(cols), dtype=int), cols)),
                                   shape=(1, ps.ncol)))
        lever.ckeys = set(frozenset(c['members']) for c in lever.cliques)
        lever.K = sp.vstack(K, format='csr') if K else sp.csr_matrix((0, ps.ncol))
        log(f'   stage {stage}: clique rows extended by {next_} memberships')
        if lever.pgons:
            log(f'   stage {stage}: {len(lever.pgons)} polygon rows re-derived over the new pose '
                f'set, +{rankfamily.regrow(lever)} memberships')
        lever.hi.h = None
        stage += 1

    # ---- anatomy of the binding cliques and of the raw (first-iteration) cliques
    if a.anatomy and lever.found:
        x, mu, rec = lever.last
        log('   ANATOMY of the violated cliques (raw = first iteration, coverage-only optimum; '
            'binding = positive dual at the end)')
        out = dict(binding=[])
        binding = [lever.cliques[i] for i in np.argsort(-lever.z)[:a.anatomy] if lever.z[i] > 1e-9]
        for name, cl, mm in (('binding', binding, mu),):
            for c in cl:
                an = anatomy(ps, c['members'], mm, lever)
                an['dual'] = float(lever.z[lever.cliques.index(c)]) if c in lever.cliques else None
                an['size_row'] = len(c['members'])
                out[name].append(an)
                log(f'   [{name}] size {an["size"]}/{an["size_row"]} mass {an["mass"]:.4f} '
                    f'point {an["point_mass"]:.4f} at ({an["point"][0]:.3f},{an["point"][1]:.3f}) '
                    f'excess {an["excess_over_point"]:.4f} anchor {an["anchor_best"]} '
                    f'regs {an["regions"]} centre ({an["centre"][0]:.2f},{an["centre"][1]:.2f}) '
                    f'angles {an["angles"]} dual {an["dual"]}')
        json.dump(out, open(os.path.join(RUNS, f'cl_{a.TAG}_anatomy.json'), 'w'), indent=1)
    # ---- the raw cliques: recorded with their masses at the time they were found
    log('RESULT tag=%s ' % a.TAG + ' '.join(
        f's{r["stage"]}=LP{r["LP"]:.6f}/M{r["M"]:.6f}/K{r["kmax"]:.6f}/conv{int(r["converged"])}'
        + (f'/exact{r["exact"]["mass_float"]:.6f}{"c" if r["exact"]["certified"] else "u"}' if r['exact'] else '')
        for r in results) + f' secs={time.time() - T0:.0f}')


def cmd_raw(a):
    """the anatomy of the max-weight cliques of an exact measure file (no LP): what the
    coverage-only optimum violates"""
    ps = Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    ps.add_exact_file(a.FILE)
    ps.finish()
    mu = np.array(ps.mu0)
    a.chord = False
    a.threads = 1
    a.lp_tlim = 1
    lever = Lever(ps, a, print)
    sup, squares, V, Inc, M, cov, w = lever.certify(mu, a.procs)
    nbr = lc.closed_graph(squares, verbose=True)
    wl = [int(v) for v in w]
    DM = lever.DM
    print(f'{a.FILE}: mass {mu.sum():.9f}, M = {M:.9f}, {len(sup)} poses')
    found = {}
    lb = int(a.lb * DM)                  # report cliques of mass > lb (1 = violated; 1-1e-9 = tight)
    Kw, Km, comp = lc.max_clique_int(nbr, wl, a.clique_time, lb=0)
    found[frozenset(Km)] = (Kw, comp)
    for j in list(np.argsort(-w)[:a.cq_top]) + list(Km):
        Kw2, Km2, comp2 = max_clique_through_members(nbr, wl, int(j), a.clique_time / 4, lb=lb)
        if Kw2 > lb:
            found[frozenset(Km2)] = (Kw2, comp2)
    print(f'{len(found)} distinct cliques of mass > {a.lb} (global max {Kw / DM:.6f}, complete {comp})')
    out = []
    for key, (Kw2, comp2) in sorted(found.items(), key=lambda kv: -kv[1][0])[:a.anatomy]:
        members = [int(sup[k]) for k in key]
        an = anatomy(ps, members, mu, lever)
        an['complete'] = bool(comp2)
        out.append(an)
        print(f'   size {an["size"]} mass {an["mass"]:.6f} point {an["point_mass"]:.4f} at '
              f'({an["point"][0]:.3f},{an["point"][1]:.3f}) excess {an["excess_over_point"]:.4f} '
              f'anchor {an["anchor_best"]} regs {an["regions"]} centre '
              f'({an["centre"][0]:.2f},{an["centre"][1]:.2f}) spread ({an["spread"][0]:.2f},'
              f'{an["spread"][1]:.2f}) angles {an["angles"]}')
        print(f'      heaviest: {an["heaviest"]}')
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


# ====================================================================== selftest
def selftest():
    """the pinwheel of LEAF_CEILING.md 4: three squares, mass 1/2 each, coverage-feasible
    (M = 1), pairwise closed-intersecting with empty triple intersection.  The point LP says 3/2;
    QSTAB says 1, and the loop must find the 3-clique and return exactly 1."""
    import tempfile
    d = tempfile.mkdtemp()
    f = os.path.join(d, 'pin.txt')
    with open(f, 'w') as fh:
        fh.write("# t = 4/1 ; sym 1\npose 0 1 2 13/5 1/2\npose 26795 100000 1480385/1000000 17/10 1/2\n"
                 "pose 57735 100000 2519615/1000000 17/10 1/2\n")
    ns = argparse.Namespace(TAG='selftest', t='4', r='1', bnd_delta=1e-6, corners=None,
                            patterns=None, chord=False, extend=True, threads=1, lp_tlim=10.0,
                            iters=20, row_cap=1000, cq_want=10, clique_time=10.0, cq_top=5,
                            ktol=1e-6, cq_age=0, procs=1, time=60, ckpt=0, price=0,
                            finalize=True, anatomy=0, rows0=1000, row_pitch=0.5, load_rows=[],
                            exact=[f], load_poses=[], Q=10 ** 5, Dc=10 ** 6,
                            pent=[5], pent_want=4, pent_cands=60, pent_restarts=20,
                            pent_time=20.0, pent_seed=1,
                            master=False, master_add=2000, master_trim=3, master_passes=2,
                            lattice_every=0,
                            price_tol=1e-7, price_pitch=0.04, price_dth=2.5, cg_want=400)
    ps = Poses(Fr(4), Fr(1), 1e-6)
    ps.add_exact_file(f)
    ps.finish()
    lever = Lever(ps, ns, print)
    lever.add_grid_rows(0.5)
    sup, squares, V, Inc, M0, cov, w = lever.certify(np.array(ps.mu0), 1)
    lever.add_rows(V)                       # as `cmd_run`: the initial support's own vertices
    h = lever.run_stage('pin')
    x, mu, rec = lever.last
    ok1 = abs(h[0]['LP'] - 1.5) < 1e-6 and abs(h[0]['kmax'] - 1.5) < 1e-9
    ok2 = abs(rec['LP'] - 1.0) < 1e-6 and rec['converged']
    fin = lever.finalize(x, mu, 1)
    ok3 = fin['certified'] and abs(fin['mass_float'] - 1.0) < 1e-6
    # the restricted master must reach the same value, with the same exact certificate
    ns2 = argparse.Namespace(**dict(vars(ns), master=True, lp_tlim=0.0))
    ps2 = Poses(Fr(4), Fr(1), 1e-6)
    ps2.add_exact_file(f)
    ps2.finish()
    lv2 = Lever(ps2, ns2, print)
    lv2.add_grid_rows(0.5)
    s2, sq2, V2, I2, M2, c2, w2 = lv2.certify(np.array(ps2.mu0), 1)
    lv2.add_rows(V2)
    h2 = lv2.run_stage('pinM')
    x2, mu2, rec2 = lv2.last
    fin2 = lv2.finalize(x2, mu2, 1)
    ok5 = (abs(rec2['LP'] - 1.0) < 1e-6 and rec2['converged'] and fin2['certified']
           and abs(fin2['mass_float'] - 1.0) < 1e-6 and lv2.rc_npos == 0)
    print(f'selftest: master converged to {rec2["LP"]:.9f} in {len(h2)} it, '
          f'master cols {lv2.master_n}/{ps2.ncol}, rc+ {lv2.rc_npos}: {ok5}')
    # a pose that meets all three (a big axis-parallel square at the centre) must be picked up by
    # the extension; one far away must not
    j = ps.add(0, 1, Fr(2), Fr(2))
    k = ps.add(0, 1, Fr(1, 2), Fr(1, 2))
    ps.finish()
    ext = ps.extend_clique([0, 1, 2], np.zeros(ps.n))
    ok4 = (j in ext) and (k not in ext) and ps.verify_clique(ext, exact=True) == 0
    print(f'selftest: first LP 1.5 with clique 1.5: {ok1}; converged to 1.0: {ok2}; exact: {ok3}; '
          f'extension: {ok4}; master: {ok5}')
    return 0 if (ok1 and ok2 and ok3 and ok4 and ok5) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(c):
        c.add_argument('--t', default='4')
        c.add_argument('--r', default='1')
        c.add_argument('--bnd-delta', type=float, default=1e-6)
        c.add_argument('--clique-time', type=float, default=60.0, help='B&B budget per global call')
        c.add_argument('--cq-top', type=int, default=30, help='max clique THROUGH the heaviest k squares')
        c.add_argument('--anatomy', type=int, default=8)
        c.add_argument('--procs', type=int, default=4)

    c = sub.add_parser('run')
    c.add_argument('TAG')
    common(c)
    c.add_argument('--exact', action='append', default=[], help='exact measure file (pose p q cx cy mass)')
    c.add_argument('--load-poses', action='append', default=[], help='t4screen float checkpoint')
    c.add_argument('--load-rows', action='append', default=[], help='t4leaf dual checkpoint (x y dual)')
    c.add_argument('--resume-rows', action='append', default=[], help='cl_*_rows.txt checkpoint (X Y D)')
    c.add_argument('--resume-cliques', action='append', default=[],
                   help='cl_*_cliques.txt checkpoint; the run must have cl_*_poses.txt as its only --exact')
    c.add_argument('--Q', type=int, default=10 ** 5)
    c.add_argument('--Dc', type=int, default=10 ** 6)
    c.add_argument('--corners', default='1111')
    c.add_argument('--patterns', default='01010101')
    c.add_argument('--chord', action='store_true')
    c.add_argument('--no-extend', dest='extend', action='store_false')
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--lp-tlim', type=float, default=60.0)
    c.add_argument('--iters', type=int, default=400)
    c.add_argument('--row-cap', type=int, default=3000)
    c.add_argument('--rows0', type=int, default=6000)
    c.add_argument('--row-pitch', type=float, default=0.08)
    c.add_argument('--cq-want', type=int, default=40, help='clique rows added per iteration')
    c.add_argument('--cq-age', type=int, default=12, help='drop a clique row after this many dual-free solves (0 = never)')
    c.add_argument('--ktol', type=float, default=1e-6)
    c.add_argument('--time', type=float, default=20000)
    c.add_argument('--ckpt', type=int, default=5, help='checkpoint every N iterations')
    c.add_argument('--price', type=int, default=0, help='number of pricing stages after convergence')
    c.add_argument('--price-pitch', type=float, default=0.04)
    c.add_argument('--price-dth', type=float, default=2.5)
    c.add_argument('--price-tol', type=float, default=1e-7)
    c.add_argument('--cg-want', type=int, default=400)
    c.add_argument('--no-finalize', dest='finalize', action='store_false')
    # ---- rank family: odd-polygon rows mu(X) <= (k-1)/2 (search/rankfamily.py); off by default
    c.add_argument('--pent', default='', help='comma-separated odd k (e.g. 5,7,9); empty = off')
    c.add_argument('--pent-want', type=int, default=6, help='polygon rows added per iteration per k')
    c.add_argument('--pent-cands', type=int, default=240, help='candidate anchor points')
    c.add_argument('--pent-restarts', type=int, default=60, help='hill-climb restarts per k')
    c.add_argument('--pent-time', type=float, default=90.0, help='separation budget per iteration')
    c.add_argument('--pent-seed', type=int, default=20260911)
    c.add_argument('--resume-pgons', action='append', default=[],
                   help='cl_*_pgons.json checkpoint; each row is rebuilt from its exact rational '
                        'anchors over the current pose set and re-verified')
    # ---- restricted master and continuous lattice pricing (search/CLMASTER.md); all default OFF
    c.add_argument('--no-warm', action='store_true',
                   help='skip the warm dual-simplex attempt and go straight to ipm+crossover '
                        '(same as --lp-tlim 0)')
    c.add_argument('--master', action='store_true',
                   help='restricted master: solve over the support + the newly priced-in columns '
                        'only, pricing every loaded column by reduced cost each iteration')
    c.add_argument('--master-add', type=int, default=2000,
                   help='columns with the most positive reduced cost added to the master per pass '
                        '(no cap on the pass that closes the master exactly)')
    c.add_argument('--master-passes', type=int, default=2,
                   help='price-and-resolve passes per iteration (0 = always to convergence).  The '
                        'cap is lifted automatically on any iteration where the row and clique '
                        'separation found nothing, so the loop can only STOP at the exact full-LP '
                        'optimum of the loaded column set')
    c.add_argument('--master-trim', type=int, default=3,
                   help='drop a master column after this many consecutive zero-mass solves (0 = never)')
    c.add_argument('--lattice-every', type=int, default=0,
                   help='price the price-pitch/price-dth lattice every M iterations (0 = off) and '
                        'add the --cg-want best poses to the LOADED set, without ending the stage')

    c = sub.add_parser('raw')
    c.add_argument('FILE')
    common(c)
    c.add_argument('--json', default=None)
    c.add_argument('--lb', type=float, default=1.0,
                   help='report cliques of mass > this (1 = violated; 0.999999 = the TIGHT cliques of a converged measure)')

    sub.add_parser('selftest')
    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    if getattr(a, 'pent', None) is not None:
        a.pent = [int(v) for v in str(a.pent).replace(',', ' ').split()]
        if any(k < 5 or k % 2 == 0 for k in a.pent):
            sys.exit('--pent takes odd k >= 5 (k = 3 is a clique row)')
    if a.cmd == 'run':
        if a.no_warm:
            a.lp_tlim = 0.0
        cmd_run(a)
    elif a.cmd == 'raw':
        cmd_raw(a)
    else:
        sys.exit(selftest())


if __name__ == '__main__':
    main()
