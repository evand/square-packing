#!/usr/bin/env python3
"""t4leaf.py -- one level-2 leaf at t = 4, taken to honest convergence.

`search/t4screen.py` (which this imports and subclasses) screened the fifteen level-2 leaves under
the k = 4 corner leaf on starved pose sets and in 450 s chunks; `search/T4SCREEN.md` records why
none of its numbers is a bound.  This script changes three things and nothing else:

1. **BOUNDARY SEMANTICS.**  `t4screen.py` assigns a pose to exactly one region by
   `level2_regions.classify`, a deterministic tie-break on closed boxes, and the packing optimum
   parks whole units of mass exactly on a slot boundary (`T4SCREEN.md` 5) -- so its value is the
   value of an LP whose branch a cell-based verifier cannot implement.  A verifier's region test is
   "is this cell's bounding box inside the box?", and a cell straddling a boundary satisfies
   NEITHER box, hence must meet BOTH regions' thresholds.  The packing-side dual of that is: a pose
   within `delta` of a region boundary may be assigned to EITHER adjacent region, at the packing's
   choice.  Implemented by DUPLICATING such a pose, one column per adjacent region (identical
   coverage and clique coefficients, different region label).  This can only RAISE the value, so
   the result stays a lower bound on what a certificate must beat.

2. **THE CHORD INEQUALITY.**  `mu({centres within 1 of a given wall}) <= 3`, one row per wall
   (`notes/level2-design.md` 2.3, Nagamochi Lemma 7(i) / Stromquist).  Still a HYPOTHESIS in this
   repo -- it is being proved separately -- so it is a flag and every value is reported with and
   without it.  The row is built from the pose coordinates (`c_y <= r + 1e-9` etc.), NOT from the
   region labels, so it is not implied by the region equalities by construction.

3. **SPEED, so that one process can converge.**  A HiGHS backend (`highspy`) that keeps the model
   across the row/clique-generation loop and warm-starts the simplex, instead of scipy's
   `linprog`, which rebuilds and re-solves from scratch every time.  `T4SCREEN.md` 8: at 100-310 s
   per LP a 450 s chunk never completed a stage.  Same LP, same answer (`--backend scipy`
   reproduces `t4screen.py` bit for bit); only the time changes.

Everything else -- semantics of the container and the squares, the exact arrangement-vertex
certification, the anchor-clique separators, the stratified pricer -- is `t4screen.py`'s, unchanged.
Nothing in `verify/`, `xcheck.py`, `lean/`, `certificates/` or `branch.py` is touched.

    python3 search/t4leaf.py 4.0 A01010101 --corners 1111 --patterns 01010101 \
        --price-pattern 01010101 --cliques --cq-wall --cq-interior --chord --bnd-delta 1e-6
"""
import argparse, json, math, os, sys, time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import packing_dual as pd                                            # noqa: E402
import level2_lp as L2                                               # noqa: E402
import anchorclique as AC                                            # noqa: E402
import t4screen as T4                                                # noqa: E402
import linecuts as LC                                                # noqa: E402
from level2_regions import classify                                  # noqa: E402

RUNS = os.path.join(REPO, 'runs')
pd.RUNS = RUNS
L2.RUNS = RUNS
T4.RUNS = RUNS

CHORD_EPS = 1e-9        # a pose counts in a wall strip if c <= r + CHORD_EPS (see 2. above)


# ====================================================================== the HiGHS backend
class Hi:
    """A HiGHS model kept alive across the row/clique generation loop.

    Model rows are addressed by LOGICAL KEY -- ('E', i) region equality, ('D', i) chord row,
    ('P', i) coverage row at `leaf.pts[i]`, ('C', k) the clique with `anchorclique.key` == k --
    because HiGHS only ever APPENDS rows, so the model order is the order rows were generated and
    is not the order any of the python-side arrays are in.  `self.rk` is the model order and
    `self.pos` the inverse map; both are rebuilt after every deletion.
    """

    def __init__(self, threads, log, solver='simplex', tlim=300.0):
        self.threads, self.log, self.solver, self.tlim = threads, log, solver, float(tlim)
        self.h = None
        self.ncols = 0
        self.rk = []
        self.pos = {}
        self.iters = 0
        self.basis = False
        import highspy
        self.hp = highspy
        self.INF = highspy.kHighsInf

    def reset(self, ncols):
        h = self.hp.Highs()
        h.setOptionValue('output_flag', False)
        h.setOptionValue('threads', self.threads)
        h.setOptionValue('presolve', 'on')
        h.addVars(ncols, np.zeros(ncols), np.full(ncols, self.INF))
        h.changeColsCost(ncols, np.arange(ncols, dtype=np.int32), -np.ones(ncols))
        self.h, self.ncols = h, ncols
        self.rk, self.pos = [], {}
        self.basis = False

    def append(self, keys, lo, hi, mat):
        if not len(keys):
            return
        M = sp.csr_matrix(mat)
        self.h.addRows(M.shape[0], np.asarray(lo, dtype=float), np.asarray(hi, dtype=float),
                       M.nnz, M.indptr.astype(np.int32), M.indices.astype(np.int32),
                       M.data.astype(float))
        for k in keys:
            self.pos[k] = len(self.rk)
            self.rk.append(k)

    def delete(self, keys):
        idx = sorted(self.pos[k] for k in keys if k in self.pos)
        if not idx:
            return
        self.h.deleteRows(len(idx), np.array(idx, dtype=np.int32))
        drop = set(idx)
        self.rk = [k for i, k in enumerate(self.rk) if i not in drop]
        self.pos = {k: i for i, k in enumerate(self.rk)}

    def bounds(self, keys, lo, hi):
        idx = np.array([self.pos[k] for k in keys if k in self.pos], dtype=np.int32)
        if not len(idx):
            return
        self.h.changeRowsBounds(len(idx), idx, np.full(len(idx), float(lo)),
                                np.full(len(idx), float(hi)))

    def run(self):
        # a fresh model has no basis, and cold dual simplex on a 10k x 10k coverage LP is far
        # slower than interior point: use ipm+crossover for the first solve after a rebuild (it
        # leaves a basis behind), then warm-started simplex for every incremental solve
        if not self.basis:
            self.h.setOptionValue('solver', 'ipm')
            self.h.setOptionValue('run_crossover', 'on')
            self.h.setOptionValue('time_limit', 1e30)
        else:
            # A warm dual-simplex reoptimisation is usually seconds, but on this LP it sometimes
            # does not terminate at all: after a model REBUILD the basis comes from an interior
            # point crossover, and appending rows to it can send the simplex into a degenerate
            # crawl (measured: 25+ minutes against 4 s for the same model built incrementally).
            # Cap it and fall back to the interior point solver, which is slower but bounded.
            self.h.setOptionValue('solver', self.solver)
            self.h.setOptionValue('time_limit', self.tlim)
        st = self.h.run()
        ms = self.h.getModelStatus()
        if ms == self.hp.HighsModelStatus.kOptimal:
            self.basis = True
        if ms != self.hp.HighsModelStatus.kOptimal:
            self.log(f'   [HiGHS {self.solver} -> {ms} after {self.tlim:.0f}s; ipm+crossover]')
            self.h.setOptionValue('solver', 'ipm')
            self.h.setOptionValue('run_crossover', 'on')
            self.h.setOptionValue('time_limit', 1e30)
            self.h.run()
            self.h.setOptionValue('solver', self.solver)
            ms = self.h.getModelStatus()
            if ms != self.hp.HighsModelStatus.kOptimal:
                self.log(f'   [HiGHS ipm+crossover -> {ms}: no solution]')
                return None
            self.basis = True
        sol = self.h.getSolution()
        info = self.h.getInfo()
        self.iters = int(info.simplex_iteration_count)
        x = np.maximum(np.asarray(sol.col_value, dtype=float), 0.0)
        d = np.asarray(sol.row_dual, dtype=float)
        return x, d, -float(info.objective_function_value)


# ====================================================================== the leaf
class BLeaf(T4.CLeaf):
    """`t4screen.CLeaf` with boundary duplication, chord rows and the HiGHS backend."""

    def __init__(self, t, r, lib, threads, log, D=1000000, delta=1e-6, chord=True,
                 backend='highs', hs_solver='simplex', lp_tlim=300.0, lines=None,
                 line_eps=-1e-12, quad=None, quad_mode='eq'):
        super().__init__(t, r, lib, threads, log, D=D)
        self.delta = float(delta)
        self.chord = bool(chord)
        # --- interior quadrant counts (tasks/depth-test).  `quad` is None (OFF: every code path
        # below is then bit for bit the one without it) or a list of four entries, one per
        # quadrant Q0 = [r,t/2]^2, Q1 = [t/2,t-r]x[r,t/2], Q2 = [r,t/2]x[t/2,t-r], Q3 = [t/2,t-r]^2
        # of the interior (label 12 + i), each a count or None (free).  `quad_mode` is 'eq'
        # (mu(Q_i) = q_i, the branch's own semantics) or 'le' (mu(Q_i) <= q_i, the cap
        # relaxation; see search/DEPTH.md 1 for why both are needed).
        self.kq = None if quad is None else list(quad)
        self.quad_mode = quad_mode
        self.nq = 0                           # quadrant rows in the model
        self.lines = lines                    # a linecuts.LineFamily, or None
        self.line_eps = float(line_eps)
        self.line_dual = np.zeros(len(lines) if lines is not None else 0)
        self.gkey = []                        # pose key of each COLUMN (duplicates share one)
        self.backend = backend
        self.hi = Hi(threads, log, hs_solver, lp_tlim) if backend == 'highs' else None
        self.hi_pattern = None                # (tuple(kc), tuple(kw)) currently in the model
        self.hi_npts = 0
        self.hi_cq = []                       # clique keys currently in the model, in model order
        self.nvar = 0                         # number of equality rows in the model
        self.chord_dual = np.zeros(4)

    # ---------------------------------------------------------------- interior quadrants
    def quad_of(self, x, y):
        """0..3 for a centre strictly inside a quadrant of the interior (Q0 bottom-left, Q1
        bottom-right, Q2 top-left, Q3 top-right); the boundaries `x = t/2`, `y = t/2` are handled
        by `region_labels`, which never calls this on them."""
        h = self.t / 2.0
        return (1 if x > h else 0) + (2 if y > h else 0)

    def label_at(self, x, y):
        """region label of a centre: 0..3 corner box, 4..11 slot, 12 interior -- or, with the
        quadrant counts on, 12 + quadrant."""
        k, j = classify(x, y, self.t, self.r)
        lab = j if k == 'C' else (4 + j if k == 'W' else 12)
        if lab == 12 and self.kq is not None:
            lab = 12 + self.quad_of(x, y)
        return lab

    def region_of(self, p):
        """`level2_lp.Leaf.region_of` (used by the pricer and the reduced-cost check), made
        quadrant-aware; identical to it when the quadrant counts are off."""
        return self.label_at(p[0], p[1])

    # ---------------------------------------------------------------- boundary semantics
    def region_labels(self, q):
        """every region a cell containing this pose's centre could straddle into.

        The region boundaries in the centre plane are `c in {r, t/2, t-r}` in each coordinate
        (`level2_regions.classify`).  A centre within `delta` of one of them is offered BOTH sides;
        `classify` is then evaluated at points strictly inside each side, so the answer does not
        depend on `classify`'s own tie-break at all."""
        d = self.delta
        bs = (self.r, self.t / 2.0, self.t - self.r)
        xs, ys = [q[0]], [q[1]]
        for b in bs:
            if abs(q[0] - b) <= d:
                xs = [b - 2 * d, b + 2 * d]
            if abs(q[1] - b) <= d:
                ys = [b - 2 * d, b + 2 * d]
        out = []
        for x in xs:
            for y in ys:
                lab = self.label_at(x, y)
                if lab not in out:
                    out.append(lab)
        return out

    def add_poses(self, cand):
        """as `t4screen.CLeaf.add_poses`, but a pose on a region boundary contributes one column
        PER adjacent region (identical geometry, different region label)."""
        n0 = len(self.poses)
        new = []
        for c in cand:
            q = pd.clamp_pose(c[0], c[1], pd.snap_angle(c[2]), self.t)
            if q is None or not pd.admissible(*q, self.t):
                continue
            k = (round(q[0] * 1e7), round(q[1] * 1e7), round(q[2] / pd.ANG_UNIT))
            if k in self.key:
                continue
            self.key[k] = len(self.poses)
            for lab in self.region_labels(q):
                self.poses.append(q)
                self.reg.append(lab)
                self.gkey.append(k)
                new.append(q)
        if new and len(self.pts):
            An = self.incidence(self.pts, new)
            self.A = An if self.A is None else sp.hstack([self.A, An], format='csr')
        k = len(new)
        if self.cliques:
            if k:
                Kn = self.clique_matrix(self.cliques, self.poses[n0:])
                self.K = sp.hstack([self.K, Kn], format='csr') if self.K.shape[1] else Kn
            elif self.K.shape[1] != len(self.poses):
                self.K = self.clique_matrix(self.cliques, self.poses)
        return k

    def sift_poses(self, keep):
        """a pose survives if ANY of its region variants does: dropping one variant of a boundary
        pose would silently undo the boundary semantics for it (and the pose key would then block
        it from ever coming back)."""
        keep = np.asarray(keep, dtype=bool).copy()
        if len(self.gkey) == len(keep):
            live = set(self.gkey[i] for i in np.nonzero(keep)[0])
            for i, g in enumerate(self.gkey):
                if g in live:
                    keep[i] = True
        idx = np.nonzero(keep)[0]
        if len(idx) == len(self.poses):
            return 0
        n = len(self.poses) - len(idx)
        self.poses = [self.poses[i] for i in idx]
        self.reg = [self.reg[i] for i in idx]
        self.gkey = [self.gkey[i] for i in idx]
        self.key = {}
        for j, g in enumerate(self.gkey):
            self.key.setdefault(g, j)
        self.A = self.A[:, idx]
        if self.cliques:
            self.K = self.K[:, idx]
        return n

    def boundary_mass(self, mu, tol=1e-7):
        """mass on poses whose centre is within `tol` of a region boundary"""
        if not len(self.poses):
            return 0.0
        P = np.asarray(self.poses, dtype=float).reshape(-1, 3)
        onb = np.zeros(len(P), bool)
        for b in (self.r, self.t / 2.0, self.t - self.r):
            onb |= (np.abs(P[:, 0] - b) <= tol) | (np.abs(P[:, 1] - b) <= tol)
        return float(mu[onb].sum())

    # ---------------------------------------------------------------- the chord rows
    def chord_matrix(self):
        """mu(centres within 1 of wall w) <= 3, w = bottom, top, left, right.  Built from the pose
        COORDINATES, not the region labels, so a pose at c_y = r exactly (which `classify` puts in
        a slot or the interior) still counts in the bottom strip."""
        P = np.asarray(self.poses, dtype=float).reshape(-1, 3)
        e = CHORD_EPS
        m = [P[:, 1] <= self.r + e, P[:, 1] >= self.t - self.r - e,
             P[:, 0] <= self.r + e, P[:, 0] >= self.t - self.r - e]
        return sp.csr_matrix(np.array([q.astype(float) for q in m]))

    # ---------------------------------------------------------------- the line-cut rows
    def line_matrix(self, poses=None):
        """`mu(B_L) <= cap(L)` for every line of the family, where `B_L` is the set of poses whose
        chord on `L` has length >= 1 (`linecuts.member`, with the SUBSET boundary convention of
        `search/LINECUTS.md` 2).  Built from the pose coordinates, like the chord rows."""
        P = np.asarray(poses if poses is not None else self.poses, dtype=float).reshape(-1, 3)
        return sp.csr_matrix(LC.member(self.lines, P, self.line_eps).astype(float))

    def line_cost(self, cand):
        """the line part of a candidate pose's reduced cost: `sum_L v_L [pose in B_L]`.  Without
        this the pricer is blind to the cuts and prices poses that the cuts forbid."""
        n = len(cand)
        if self.lines is None or not len(self.line_dual) or n == 0:
            return np.zeros(n)
        idx = np.nonzero(self.line_dual > 1e-12)[0]
        if not len(idx):
            return np.zeros(n)
        P = np.asarray(cand, dtype=float).reshape(-1, 3)
        sub = LC.LineFamily(self.t)
        sub.nx, sub.ny, sub.c0 = self.lines.nx[idx], self.lines.ny[idx], self.lines.c0[idx]
        sub.phi, sub.cap = self.lines.phi[idx], self.lines.cap[idx]
        return LC.member(sub, P, self.line_eps).astype(float).T @ self.line_dual[idx]

    def clique_cost(self, cand, zmin=1e-9, topz=500):
        """`t4screen.price`, `pose_rc` and `rc_check` all take the non-coverage part of a pose's
        reduced cost from here, so folding the line duals in makes all three line-aware at once."""
        return super().clique_cost(cand, zmin, topz) + self.line_cost(cand)

    def line_report(self, mu):
        if self.lines is None:
            return None
        m = self.line_matrix() @ np.asarray(mu, dtype=float)
        v = m - self.lines.cap
        return m, v

    def eq_matrix(self, kc, kw):
        reg = np.array(self.reg)
        rows, vals = [], []
        for i in range(4):
            if kc[i] is not None:
                rows.append((reg == i).astype(float)); vals.append(kc[i])
        for j in range(8):
            if kw[j] is not None:
                rows.append((reg == 4 + j).astype(float)); vals.append(kw[j])
        return (sp.csr_matrix(np.array(rows)) if rows else None), np.array(vals)

    def quad_matrix(self):
        """the quadrant rows `mu({reg == 12 + i}) (= or <=) q_i`, in the order i = 0..3 over the
        pinned quadrants; (None, empty) when the quadrant counts are off."""
        if self.kq is None:
            return None, np.zeros(0)
        reg = np.array(self.reg)
        rows, vals = [], []
        for i in range(4):
            if self.kq[i] is not None:
                rows.append((reg == 12 + i).astype(float)); vals.append(float(self.kq[i]))
        if not rows:
            return None, np.zeros(0)
        self.nq = len(rows)
        return sp.csr_matrix(np.array(rows)), np.array(vals)

    # ---------------------------------------------------------------- LP
    def solve(self, kc, kw, method='highs', no_cliques=False, no_chord=False, no_lines=False):
        if self.backend != 'highs':
            return self._solve_scipy(kc, kw, method, no_cliques, no_chord, no_lines)
        t0 = time.time()
        out = self._solve_highs(kc, kw, no_cliques, no_chord, no_lines)
        self.lp_secs = time.time() - t0
        return out

    def _solve_fresh(self, kc, kw, no_cliques, no_chord, no_lines=False):
        """the matched variants (clique rows and/or chord rows switched off) on a ONE-SHOT model.

        Relaxing rows in place by widening their bounds to +-inf leaves the warm basis with
        nonbasic-at-bound logicals whose bounds no longer exist; HiGHS then spends longer repairing
        that than solving the model outright (26 min against 4 s for the base model on a
        12.8k x 22.8k instance).  A fresh model solved by ipm+crossover is predictable, and it
        leaves the main model's basis alone."""
        import highspy
        n = len(self.poses)
        INF = highspy.kHighsInf
        h = highspy.Highs()
        h.setOptionValue('output_flag', False)
        h.setOptionValue('threads', self.threads)
        h.setOptionValue('solver', 'ipm')
        h.setOptionValue('run_crossover', 'on')
        h.addVars(n, np.zeros(n), np.full(n, INF))
        h.changeColsCost(n, np.arange(n, dtype=np.int32), -np.ones(n))
        E, ev = self.eq_matrix(kc, kw)
        blocks, lo, hi = [], [], []
        if E is not None:
            blocks.append(E); lo.append(ev); hi.append(ev)
        Qm, qv = self.quad_matrix()
        nq = 0
        if Qm is not None:
            nq = Qm.shape[0]
            blocks.append(Qm)
            lo.append(qv.copy() if self.quad_mode == 'eq' else np.full(nq, -INF)); hi.append(qv)
        if self.chord and not no_chord:
            blocks.append(self.chord_matrix())
            lo.append(np.full(4, -INF)); hi.append(np.full(4, 3.0))
        nL = 0
        if self.lines is not None and not no_lines:
            nL = len(self.lines)
            blocks.append(self.line_matrix())
            lo.append(np.full(nL, -INF)); hi.append(self.lines.cap.copy())
        blocks.append(self.A)
        lo.append(np.full(len(self.pts), -INF)); hi.append(np.ones(len(self.pts)))
        if self.cliques and not no_cliques:
            blocks.append(self.K)
            lo.append(np.full(len(self.cliques), -INF)); hi.append(np.ones(len(self.cliques)))
        M = sp.vstack(blocks, format='csr')
        h.addRows(M.shape[0], np.concatenate(lo), np.concatenate(hi), M.nnz,
                  M.indptr.astype(np.int32), M.indices.astype(np.int32), M.data.astype(float))
        h.run()
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
            return None
        sol = h.getSolution()
        info = h.getInfo()
        mu = np.maximum(np.asarray(sol.col_value, dtype=float), 0.0)
        d = np.asarray(sol.row_dual, dtype=float)
        ne = (E.shape[0] if E is not None else 0) + nq
        o = ne + (4 if (self.chord and not no_chord) else 0) + nL
        y = np.maximum(-d[o:o + len(self.pts)], 0.0)
        lam = -d[:ne]
        return mu, y, -float(info.objective_function_value), lam

    def _solve_highs(self, kc, kw, no_cliques, no_chord, no_lines=False):
        if no_cliques or no_chord or no_lines:
            return self._solve_fresh(kc, kw, no_cliques, no_chord, no_lines)
        n = len(self.poses)
        pat = (tuple(kc), tuple(kw), None if self.kq is None else tuple(self.kq), self.quad_mode)
        use_cq = bool(self.cliques)
        rebuild = (self.hi.h is None or self.hi.ncols != n or self.hi_pattern != pat
                   or self.hi_npts > len(self.pts))
        INF = self.hi.INF
        if rebuild:
            self.hi.reset(n)
            E, ev = self.eq_matrix(kc, kw)
            if E is not None:
                self.hi.append([('E', i) for i in range(E.shape[0])], ev, ev, E)
                self.nvar = E.shape[0]
            else:
                self.nvar = 0
            Qm, qv = self.quad_matrix()
            self.nq = 0
            if Qm is not None:
                self.nq = Qm.shape[0]
                self.hi.append([('Q', i) for i in range(self.nq)],
                               qv.copy() if self.quad_mode == 'eq' else np.full(self.nq, -INF),
                               qv, Qm)
                self.nvar += self.nq
            if self.chord:
                C = self.chord_matrix()
                self.hi.append([('D', i) for i in range(4)], np.full(4, -INF), np.full(4, 3.0), C)
            if self.lines is not None:
                nL = len(self.lines)
                self.hi.append([('L', i) for i in range(nL)], np.full(nL, -INF),
                               self.lines.cap.copy(), self.line_matrix())
            self.hi.append([('P', i) for i in range(len(self.pts))],
                           np.full(len(self.pts), -INF), np.ones(len(self.pts)), self.A)
            self.hi_npts = len(self.pts)
            self.hi_cq = []
            self.hi_pattern = pat
        # --- rows added since the last solve
        if len(self.pts) > self.hi_npts:
            k0 = self.hi_npts
            self.hi.append([('P', i) for i in range(k0, len(self.pts))],
                           np.full(len(self.pts) - k0, -INF), np.ones(len(self.pts) - k0),
                           self.A[k0:])
            self.hi_npts = len(self.pts)
        # --- cliques: delete the sifted ones, append the new ones
        cur = [AC.key(cl) for cl in self.cliques]
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
        # --- variants: relax the clique / chord rows in place
        if no_cliques and self.hi_cq:
            self.hi.bounds([('C', k) for k in self.hi_cq], -INF, INF)
        if no_chord and self.chord:
            self.hi.bounds([('D', i) for i in range(4)], -INF, INF)
        res = self.hi.run()
        if no_cliques and self.hi_cq:
            self.hi.bounds([('C', k) for k in self.hi_cq], -INF, 1.0)
        if no_chord and self.chord:
            self.hi.bounds([('D', i) for i in range(4)], -INF, 3.0)
        if res is None:
            return None
        x, d, obj = res
        mu = x
        y = np.zeros(len(self.pts))
        z = np.zeros(len(self.cliques))
        lam = np.zeros(self.nvar)
        chd = np.zeros(4)
        lnd = np.zeros(len(self.lines) if self.lines is not None else 0)
        ci = {k: i for i, k in enumerate(cur)}
        for mi, k in enumerate(self.hi.rk):
            v = -d[mi]
            if k[0] == 'P':
                y[k[1]] = max(v, 0.0)
            elif k[0] == 'E':
                lam[k[1]] = v
            elif k[0] == 'Q':
                lam[self.nvar - self.nq + k[1]] = v
            elif k[0] == 'D':
                chd[k[1]] = max(v, 0.0)
            elif k[0] == 'L':
                lnd[k[1]] = max(v, 0.0)
            elif k[0] == 'C':
                z[ci[k[1]]] = max(v, 0.0)
        self.chord_dual = chd
        self.line_dual = lnd
        if not no_cliques:
            self.z = z if use_cq else np.zeros(0)
            if use_cq and len(self.z) == len(self.cage):
                self.cage = np.where(self.z > 1e-9, 0, self.cage + 1)
        return mu, y, obj, lam

    def _solve_scipy(self, kc, kw, method, no_cliques, no_chord, no_lines=False):
        from scipy.optimize import linprog
        n = len(self.poses)
        E, ev = self.eq_matrix(kc, kw)
        Qm, qv = self.quad_matrix()
        nq_le = 0
        if Qm is not None:
            if self.quad_mode == 'eq':
                E = Qm if E is None else sp.vstack([E, Qm], format='csr')
                ev = np.concatenate([ev, qv]) if len(ev) else qv
            else:
                nq_le = Qm.shape[0]
        blocks, ub = [], []
        if nq_le:
            blocks.append(Qm); ub.append(qv)
        if self.chord and not no_chord:
            blocks.append(self.chord_matrix()); ub.append(np.full(4, 3.0))
        nL = 0
        if self.lines is not None and not no_lines:
            nL = len(self.lines)
            blocks.append(self.line_matrix()); ub.append(self.lines.cap.copy())
        blocks.append(self.A); ub.append(np.ones(len(self.pts)))
        if self.cliques and not no_cliques:
            blocks.append(self.K); ub.append(np.ones(len(self.cliques)))
        A = sp.vstack(blocks, format='csr')
        b = np.concatenate(ub)
        t0 = time.time()
        res = None
        for meth in (method, 'highs-ds', 'highs-ipm'):
            res = linprog(c=-np.ones(n), A_ub=A, b_ub=b, A_eq=(E.toarray() if E is not None else None),
                          b_eq=(ev if E is not None else None), bounds=(0, None), method=meth)
            if res.success:
                break
        self.lp_secs = time.time() - t0
        if not res.success:
            return None
        mu = np.maximum(res.x, 0.0)
        yz = np.maximum(-res.ineqlin.marginals, 0.0)
        o = 0
        qlam = np.zeros(0)
        if nq_le:
            qlam = yz[:nq_le]; o = nq_le
        if self.chord and not no_chord:
            self.chord_dual = yz[o:o + 4]; o += 4
        if nL:
            self.line_dual = yz[o:o + nL]; o += nL
        y = yz[o:o + len(self.pts)]
        if not no_cliques:
            self.z = yz[o + len(self.pts):] if (self.cliques and not no_cliques) else np.zeros(0)
            if self.cliques and len(self.z) == len(self.cage):
                self.cage = np.where(self.z > 1e-9, 0, self.cage + 1)
        lam = -res.eqlin.marginals if E is not None else np.zeros(0)
        if nq_le:
            lam = np.concatenate([lam, qlam])
        return mu, y, -res.fun, lam


# ====================================================================== reporting
def split(leaf, mu):
    """region split of a measure: the four corner boxes, the eight slots, the interior"""
    reg = np.array(leaf.reg)
    mc = [float(mu[reg == i].sum()) for i in range(4)]
    mw = [float(mu[reg == 4 + j].sum()) for j in range(8)]
    mi = float(mu[reg >= 12].sum())            # 12 alone, or 12..15 with the quadrants on
    return mc, mw, mi


def quads(leaf, mu):
    """mass per interior quadrant (by column label), or None with the quadrant counts off"""
    if leaf.kq is None:
        return None
    reg = np.array(leaf.reg)
    return [float(mu[reg == 12 + i].sum()) for i in range(4)]


def strips(leaf, mu):
    P = np.asarray(leaf.poses, dtype=float).reshape(-1, 3)
    e = CHORD_EPS
    m = [P[:, 1] <= leaf.r + e, P[:, 1] >= leaf.t - leaf.r - e,
         P[:, 0] <= leaf.r + e, P[:, 0] >= leaf.t - leaf.r - e]
    return [float(mu[q].sum()) for q in m]


def save_dual(leaf, y, path):
    sel = np.nonzero(y > 1e-12)[0]
    with open(path + '.tmp', 'w') as f:
        f.write(f'# t4leaf point duals t={leaf.t} n={len(sel)} sum={float(y[sel].sum()):.9f}\n')
        for i in sel:
            f.write(f'{leaf.pts[i, 0]:.17g} {leaf.pts[i, 1]:.17g} {y[i]:.12g}\n')
    os.replace(path + '.tmp', path)


def save_measure(leaf, mu, path, t, note=''):
    sel = np.nonzero(mu > 1e-12)[0]
    with open(path + '.tmp', 'w') as f:
        f.write(f'# t4leaf measure t={t} n={len(sel)} mass={float(mu.sum()):.9f} {note}\n')
        for i in sel:
            cx, cy, th = leaf.poses[i]
            f.write(f'pose {cx:.17g} {cy:.17g} {math.degrees(th):.17g} {mu[i]:.12g}'
                    f'   # region {leaf.reg[i]}\n')
    os.replace(path + '.tmp', path)


# ====================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('T'); ap.add_argument('TAG')
    ap.add_argument('--r', type=float, default=1.0)
    ap.add_argument('--corners', default='1111')
    ap.add_argument('--patterns', default='........')
    ap.add_argument('--stages', type=int, default=200)
    ap.add_argument('--warm-stages', type=int, default=0)
    ap.add_argument('--rowloops', type=int, default=25)
    ap.add_argument('--seed-pitch', type=float, default=0.25)
    ap.add_argument('--seed-dth', type=float, default=15.0)
    ap.add_argument('--row-pitch', type=float, default=0.08)
    ap.add_argument('--price-pitch', type=float, default=0.04)
    ap.add_argument('--price-dth', type=float, default=2.5)
    ap.add_argument('--price-tol', type=float, default=1e-7)
    ap.add_argument('--cg-want', type=int, default=400)
    ap.add_argument('--pose-max', type=int, default=9000)
    ap.add_argument('--corner-rows', type=int, default=1200)
    ap.add_argument('--row-cap', type=int, default=8000)
    ap.add_argument('--row-age', type=int, default=0, help='OFF by default and meant to stay off')
    ap.add_argument('--warm', action='append', default=[])
    ap.add_argument('--load-poses', action='append', default=[])
    ap.add_argument('--price-pattern', default=None, help='default: the first --patterns entry')
    ap.add_argument('--method', default='highs-ipm', help='scipy backend only')
    ap.add_argument('--backend', default='highs', choices=('highs', 'scipy'))
    ap.add_argument('--hs-solver', default='simplex', choices=('simplex', 'ipm'))
    ap.add_argument('--lp-tlim', type=float, default=300.0,
                    help='seconds before a warm simplex solve is abandoned for ipm+crossover')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--time', type=float, default=21000)
    ap.add_argument('--ckpt-secs', type=float, default=600, help='checkpoint at least this often')
    # --- the two new semantics
    ap.add_argument('--bnd-delta', type=float, default=1e-6,
                    help='a pose within this of a region boundary gets one column per adjacent region')
    ap.add_argument('--chord', action='store_true', help='add mu(wall strip) <= 3 (a HYPOTHESIS)')
    # --- the interior quadrant counts (tasks/depth-test, search/DEPTH.md)
    ap.add_argument('--quad', default=None,
                    help='four digits, the count pinned in each interior quadrant Q0..Q3 '
                         '(bottom-left, bottom-right, top-left, top-right; "." = free).  A pose '
                         'within --bnd-delta of x = t/2 or y = t/2 gets one column per adjacent '
                         'quadrant, as on the slot boundaries.  Off (default): no change at all.')
    ap.add_argument('--quad-mode', default='eq', choices=('eq', 'le'),
                    help='eq: mu(Q_i) = q_i (the branch; with sum q = 4 the LP is then 12 or '
                         'infeasible).  le: mu(Q_i) <= q_i, the cap relaxation whose value is a '
                         'continuous margin (>= 12 iff the eq child is feasible).')
    ap.add_argument('--load-poses-c4', action='append', default=[],
                    help='like --load-poses, but only the poses carrying mass in the checkpoint, '
                         'together with their three images under the C4 rotations about the '
                         'centre (the symmetry of a C4-invariant slot pattern)')
    # --- the general-line chord cuts (search/LINECUTS.md)
    ap.add_argument('--lines', action='store_true',
                    help='add mu({chord of S on L >= 1}) <= ceil(Lam_L) - 1 for a family of lines')
    ap.add_argument('--line-pitch', type=float, default=0.05,
                    help='pitch of the horizontal/vertical line grid (0 = no axis-parallel lines)')
    ap.add_argument('--line-diag-pitch', type=float, default=0.0,
                    help='pitch of the two 45 deg line families (0 = off)')
    ap.add_argument('--line-exact', default='1,2,3',
                    help='extra axis-parallel lines at exactly these offsets')
    ap.add_argument('--line-eps', type=float, default=-1e-12,
                    help='<0: exclude poses on the boundary d = D(psi) (a SUBSET of the true set, '
                         'hence unconditionally valid).  0: the exact closed criterion.')
    ap.add_argument('--line-maxcap', type=int, default=3,
                    help='skip lines whose bound is above this (they are valid but weak)')
    ap.add_argument('--variants', action='store_true',
                    help='at each stage also value the same instance with cliques and/or the chord rows off')
    ap.add_argument('--variant-every', type=int, default=1,
                    help='with --variants: only on every Nth stage (the pure solve is expensive)')
    ap.add_argument('--load-rows', action='append', default=[],
                    help='seed the row set from a t4leaf dual checkpoint (`x y dual` lines)')
    # --- cliques (as t4screen.py)
    ap.add_argument('--cliques', action='store_true')
    ap.add_argument('--cq-wall', action='store_true')
    ap.add_argument('--cq-interior', action='store_true')
    ap.add_argument('--cq-load', action='append', default=[])
    ap.add_argument('--cq-want', type=int, default=200)
    ap.add_argument('--cq-max', type=int, default=900)
    ap.add_argument('--cq-age', type=int, default=6)
    ap.add_argument('--cq-pitch', type=float, default=0.02)
    ap.add_argument('--cq-top', type=int, default=250)
    ap.add_argument('--cq-int-pitch', type=float, default=0.04)
    ap.add_argument('--cq-int-top', type=int, default=60)
    ap.add_argument('--cq-ndir', type=int, default=12)
    ap.add_argument('--cq-neps', type=int, default=6)
    ap.add_argument('--cq-nrho', type=int, default=6)
    ap.add_argument('--ktol', type=float, default=1e-6)
    ap.add_argument('--D', type=int, default=1000000)
    a = ap.parse_args()
    if a.cliques and not (a.cq_wall or a.cq_interior):
        a.cq_interior = True
    if not a.cliques and not a.cq_load:
        a.cq_wall = a.cq_interior = False
    patterns = a.patterns.split(',')
    if a.price_pattern is None:
        a.price_pattern = patterns[0]

    t = float(Fr(a.T))
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'tl_{a.TAG}.log'), 'a')

    def log(m):
        print(m, flush=True); logf.write(m + '\n'); logf.flush()

    T0 = time.time()
    lib = pd.build_lib()
    F = None
    if a.lines:
        ex = [float(v) for v in a.line_exact.split(',') if v.strip()]
        F = LC.build(t, a.line_pitch, exact=ex, diag_pitch=a.line_diag_pitch,
                     maxcap=a.line_maxcap, hv=(a.line_pitch > 0))
    kq = None
    if a.quad is not None:
        if len(a.quad) != 4:
            raise SystemExit('--quad wants four characters, e.g. 1111 or 2..0')
        kq = [None if c == '.' else float(c) for c in a.quad]
    leaf = BLeaf(t, a.r, lib, a.threads, log, D=a.D, delta=a.bnd_delta, chord=a.chord,
                 backend=a.backend, hs_solver=a.hs_solver, lp_tlim=a.lp_tlim, lines=F,
                 line_eps=a.line_eps, quad=kq, quad_mode=a.quad_mode)
    kc = [None if c == '.' else float(c) for c in a.corners]
    log(f'# t4leaf t={t} r={a.r} corners={a.corners} patterns={a.patterns} chord={a.chord} '
        f'bnd_delta={a.bnd_delta} backend={a.backend} args={vars(a)}')
    if kq is not None:
        log(f'   interior quadrants: {a.quad} ({a.quad_mode}); labels 12..15 = Q0..Q3')
    if F is not None:
        caps = {}
        for i in range(len(F)):
            caps[int(F.cap[i])] = caps.get(int(F.cap[i]), 0) + 1
        log(f'   line family: {len(F)} lines, caps {sorted(caps.items())}, eps={a.line_eps:g} '
            f'({"subset, valid" if a.line_eps < 0 else "closed criterion"})')

    # ---- rows and columns
    g = np.arange(a.row_pitch / 2, t, a.row_pitch)
    GX, GY = np.meshgrid(g, g, indexing='ij')
    leaf.add_points(np.column_stack([GX.ravel(), GY.ravel()]))
    cand = L2.seed_poses_full(t, a.seed_pitch, a.seed_dth)
    for w in a.warm:
        cand += L2.read_measure_poses(w, t)
    for w in a.load_poses:
        cand += T4.load_poses(w)
    for w in a.load_poses_c4:
        sup = [(float(q[1]), float(q[2]), math.radians(float(q[3]))) for q in
               (l.split() for l in open(w)) if len(q) >= 5 and q[0] == 'pose' and float(q[4]) > 1e-9]
        imgs = []
        for (x, y, th) in sup:              # a square's angle is unchanged by a quarter turn
            imgs += [(x, y, th), (t - y, x, th), (t - x, t - y, th), (y, t - x, th)]
        log(f'   +{len(imgs)} C4 images of {len(sup)} support poses from {os.path.basename(w)}')
        cand += imgs
    n = leaf.add_poses(cand)
    leaf.add_points(L2.pose_corners(leaf.poses[:a.corner_rows]))
    for rf in a.load_rows:
        pts = [[float(q[0]), float(q[1])] for q in
               (l.split() for l in open(rf)) if len(q) >= 2 and not q[0].startswith('#')]
        log(f'   +{leaf.add_points(np.array(pts))} rows from {os.path.basename(rf)}')
    leaf.rfix = len(leaf.pts)
    ndup = len(leaf.poses) - len(set(leaf.gkey))
    log(f'   seed: {n} columns for {len(set(leaf.gkey))} poses ({ndup} boundary duplicates), '
        f'{len(leaf.pts)} rows ({time.time()-T0:.0f}s)')
    for cf in a.cq_load:
        leaf.load_pool(cf, log)

    hist = []
    _hj = os.path.join(RUNS, f'tl_{a.TAG}.json')
    if os.path.exists(_hj):
        try:
            hist = json.load(open(_hj)).get('hist', [])
        except Exception:
            hist = []
    last_ck = [time.time()]
    varflag = [True]

    def line_stats(mu):
        """(number of lines carrying dual, sum of the line duals, largest mu(B_L) - cap)"""
        if leaf.lines is None:
            return None
        m, v = leaf.line_report(mu)
        nz = int((leaf.line_dual > 1e-9).sum())
        return dict(nz=nz, dualsum=float(leaf.line_dual.sum()),
                    tight=int((v > -1e-9).sum()), maxexc=float(v.max()),
                    top=[[leaf.lines.tag[i], float(leaf.line_dual[i]), float(m[i])]
                         for i in np.argsort(-leaf.line_dual)[:8] if leaf.line_dual[i] > 1e-9])

    def checkpoint(mu, y):
        T4.save_poses(leaf, mu, os.path.join(RUNS, f'tl_{a.TAG}_poses.txt'), t)
        if leaf.cliques:
            T4.save_cliques(leaf, os.path.join(RUNS, f'tl_{a.TAG}_cliques.txt'))
        if y is not None:
            save_dual(leaf, y, os.path.join(RUNS, f'tl_{a.TAG}_dual.txt'))
        json.dump(dict(t=t, tag=a.TAG, args=vars(a), hist=hist),
                  open(os.path.join(RUNS, f'tl_{a.TAG}.json'), 'w'), indent=1)
        last_ck[0] = time.time()

    def value(kc, kw, tag):
        rec = None
        pmu = py = None
        for it in range(a.rowloops):
            ndrop = 0
            if a.cliques or leaf.pool:
                ndrop = leaf.sift_cliques(a.cq_age)
            sol = leaf.solve(kc, kw, a.method)
            if sol is None:
                log(f'   [{tag}.{it}] LP infeasible/failed')
                return None
            mu, y, val, lam = sol
            pmu, py = mu, y
            M, bad, nv, ovf = leaf.certify(mu)
            kmax, ncq, best = (0.0, 0, None)
            if a.cliques or leaf.pool:
                room = max(a.cq_max - len(leaf.cliques), 0)
                kmax, ncq, best = leaf.separate(mu, a, log, want=min(a.cq_want, room))
            mc, mw, mi = split(leaf, mu)
            L = val / max(M, kmax, 1.0)
            ls = line_stats(mu)
            log(f'   [{tag}.{it}] cols={len(leaf.poses)} rows={len(leaf.pts)} cq={len(leaf.cliques)} '
                f'LP={val:.6f} M={M:.9f} kmax={kmax:.6f} L={L:.6f} bad={len(bad)} +cq={ncq} '
                f'-cq={ndrop} int={mi:.6f} bnd={leaf.boundary_mass(mu):.6f} '
                f'chd={np.round(leaf.chord_dual, 4).tolist() if leaf.chord else "-"} '
                + (f'lin=nz{ls["nz"]}/dsum{ls["dualsum"]:.4f}/tight{ls["tight"]}/'
                   f'exc{ls["maxexc"]:+.2e} ' if ls else '')
                + f'(lp {leaf.lp_secs:.1f}s it{leaf.hi.iters if leaf.hi else 0}, {time.time()-T0:.0f}s)')
            rec = dict(tag=tag, it=it, LP=val, M=M, kmax=kmax, L=L, cols=len(leaf.poses),
                       poses=len(set(leaf.gkey)), rows=len(leaf.pts), cq=len(leaf.cliques),
                       bad=len(bad), interior=mi, corners=mc, slots=mw, quads=quads(leaf, mu),
                       strips=strips(leaf, mu), bnd=leaf.boundary_mass(mu),
                       chord_dual=[float(v) for v in leaf.chord_dual] if leaf.chord else None,
                       secs=time.time() - T0)
            nrow = 0
            if M > 1 + 1e-9 and len(bad):
                o = np.argsort(-bad[:, 2])[:a.row_cap]
                nrow = leaf.add_points(bad[o, :2])
            done = (nrow == 0 and ncq == 0 and M <= 1 + 1e-9 and kmax <= 1 + a.ktol)
            rec['converged'] = bool(done)
            if done or (nrow == 0 and ncq == 0):
                break
            if time.time() - T0 > a.time:
                log('   [time limit inside the inner loop]')
                break
            if time.time() - last_ck[0] > a.ckpt_secs:
                checkpoint(mu, y)
        # ---- settle: the loop's last `separate` call added clique rows AFTER the last solve, so
        # the record above is not a value of the model that is now in memory.  Re-solve the base
        # model, recompute `M` and `kmax` (with `want=0`, i.e. separate but add nothing) against
        # THAT measure, and report the variants against the same rows and columns.  Everything in
        # `rec` from here on is one self-consistent LP.
        sol = leaf.solve(kc, kw, a.method)
        if sol is None:
            return None
        mu, y, val, lam = sol
        M, bad, nv, ovf = leaf.certify(mu)
        kmax = 0.0
        if a.cliques or leaf.pool:
            kmax, _, _ = leaf.separate(mu, a, log, want=0)
        mc, mw, mi = split(leaf, mu)
        rec = dict(tag=tag, it='final', LP=val, M=M, kmax=kmax, L=val / max(M, kmax, 1.0),
                   cols=len(leaf.poses), poses=len(set(leaf.gkey)), rows=len(leaf.pts),
                   cq=len(leaf.cliques), bad=len(bad), interior=mi, corners=mc, slots=mw,
                   quads=quads(leaf, mu), quad_dual=(None if leaf.kq is None else
                                                     [float(v) for v in lam[len(lam) - leaf.nq:]]),
                   strips=strips(leaf, mu), bnd=leaf.boundary_mass(mu),
                   chord_dual=[float(v) for v in leaf.chord_dual] if leaf.chord else None,
                   lines=line_stats(mu),
                   secs=time.time() - T0,
                   converged=bool(M <= 1 + 1e-9 and kmax <= 1 + a.ktol))
        log(f'   [{tag}.final] cols={len(leaf.poses)} rows={len(leaf.pts)} cq={len(leaf.cliques)} '
            f'LP={val:.6f} M={M:.9f} kmax={kmax:.6f} int={mi:.6f} '
            f'bnd={rec["bnd"]:.6f} conv={rec["converged"]} (lp {leaf.lp_secs:.1f}s)'
            + (f' quads={np.round(rec["quads"], 6).tolist()} '
               f'qdual={np.round(rec["quad_dual"], 4).tolist()}' if rec['quads'] else ''))
        if rec['lines']:
            log(f'   [{tag}.lines] {rec["lines"]["nz"]} of {len(leaf.lines)} lines carry dual, '
                f'sum {rec["lines"]["dualsum"]:.6f}, {rec["lines"]["tight"]} tight, '
                f'max excess {rec["lines"]["maxexc"]:+.3e}; top '
                + ', '.join(f'{q[0]} v={q[1]:.4f} mu={q[2]:.4f}' for q in rec['lines']['top']))
        # matched variants on the SAME columns and the SAME rows.  Relaxing every clique row at
        # once throws the warm basis away, so this is by far the most expensive solve of a stage
        # (25 min against 4 s for the base model on a 12.8k x 22.8k instance); `--variant-every`
        # keeps it off the critical path.  The chord variant needs no solve at all when the chord
        # rows carry zero dual: that dual solution is feasible for the model without them and has
        # the same objective, so `LP_nochord = LP` exactly.
        if leaf.chord and float(np.max(leaf.chord_dual)) <= 1e-12:
            rec['LP_nochord'] = val
            rec['chord_worth'] = 0.0
            rec['nochord_by_dual'] = True
        # the same dual argument for the line rows: if every line dual is 0, the optimal dual
        # solution of this model is feasible for the model WITHOUT the line rows and has the same
        # objective, so `LP_nolines = LP` exactly -- no solve needed.
        if leaf.lines is not None and float(np.max(leaf.line_dual)) <= 1e-12:
            rec['LP_nolines'] = val
            rec['lines_worth'] = 0.0
            rec['nolines_by_dual'] = True
        if a.variants:
            for (nc, nd, nl, name) in ((True, False, False, 'pure'),
                                       (False, True, False, 'nochord'),
                                       (False, False, True, 'nolines'),
                                       (True, False, True, 'pure_nolines'),
                                       (True, True, True, 'pure_nochord')):
                if nc and not leaf.cliques:
                    continue
                if nd and not leaf.chord:
                    continue
                if nl and leaf.lines is None:
                    continue
                if nd and rec.get('nochord_by_dual') and not (nc or nl):
                    continue                     # already known exactly, from the duals
                if nl and rec.get('nolines_by_dual') and not (nc or nd):
                    continue
                if not varflag[0]:
                    continue
                s2 = leaf.solve(kc, kw, a.method, no_cliques=nc, no_chord=nd, no_lines=nl)
                if s2 is None:
                    continue
                mu2, _, val2, _ = s2
                M2, _, _, _ = leaf.certify(mu2)
                rec[f'LP_{name}'] = val2
                rec[f'M_{name}'] = M2
                mc2, mw2, mi2 = split(leaf, mu2)
                rec[f'int_{name}'] = mi2
                save_measure(leaf, mu2, os.path.join(RUNS, f'tl_{a.TAG}_measure_{name}.txt'), t,
                             note=f'{tag} variant={name} LP={val2:.9f} M={M2:.9f}')
                if name == 'pure':
                    rec['gain'] = val2 - rec['LP']
                if name == 'nochord':
                    rec['chord_worth'] = val2 - rec['LP']
                if name == 'nolines':
                    rec['lines_worth'] = val2 - rec['LP']
        elif leaf.cliques and varflag[0]:
            s2 = leaf.solve(kc, kw, a.method, no_cliques=True)
            if s2 is not None:
                mu2, _, val2, _ = s2
                M2, _, _, _ = leaf.certify(mu2)
                rec['LP_pure'] = val2; rec['M_pure'] = M2; rec['gain'] = val2 - rec['LP']
        return rec, mu, y, lam

    def lam_map(kw, lam):
        out = {}
        idx = 0
        for i in range(4):
            if kc[i] is not None:
                out[i] = float(lam[idx]); idx += 1
        for j in range(8):
            if kw[j] is not None:
                out[4 + j] = float(lam[idx]); idx += 1
        if leaf.kq is not None:
            for i in range(4):
                if leaf.kq[i] is not None:
                    out[12 + i] = float(lam[idx]); idx += 1
        return out

    results = {}
    ppat = a.price_pattern
    for stage in range(a.stages):
        pats = [ppat] if stage < a.warm_stages else ([ppat] + [p for p in patterns if p != ppat])
        pmu = py = plam = None
        varflag[0] = bool(a.variant_every and stage % a.variant_every == 0)
        for pat in pats:
            kw = [None if c == '.' else float(c) for c in pat]
            out = value(kc, kw, f's{stage}.{pat}')
            if out is None:
                results[pat] = None
                continue
            rec, mu, y, lam = out
            rec['stage'] = stage; rec['pattern'] = pat
            hist.append(rec)
            results[pat] = rec
            if pat == ppat:
                pmu, py, plam = mu, y, lam_map(kw, lam)
            extra = ''.join(f' | {k[3:]} LP={rec[k]:.6f}' for k in
                            ('LP_pure', 'LP_nochord', 'LP_nolines', 'LP_pure_nolines',
                             'LP_pure_nochord') if k in rec)
            log(f'   STAGE {stage} {pat}: LP={rec["LP"]:.6f} M={rec["M"]:.9f} kmax={rec["kmax"]:.6f} '
                f'cols={rec["cols"]} rows={rec["rows"]} cq={rec["cq"]} int={rec["interior"]:.6f} '
                f'strips={np.round(rec["strips"], 4).tolist()} conv={rec["converged"]}'
                + (f' quads={np.round(rec["quads"], 4).tolist()}' if rec.get('quads') else '')
                + extra)
            if time.time() - T0 > a.time:
                break
        if pmu is not None:
            checkpoint(pmu, py)
            save_measure(leaf, pmu, os.path.join(RUNS, f'tl_{a.TAG}_measure.txt'), t,
                         note=f'stage={stage} pattern={ppat}')
        if time.time() - T0 > a.time or stage == a.stages - 1 or pmu is None:
            log('   [stopping: budget or last stage]')
            break
        chk = T4.rc_check(leaf, pmu, py, plam, log)
        if chk is not None:
            log(f'   rc-check: {chk:.2e}')
        new, gap, _ = T4.price(leaf, py, plam, a, log)
        heavy = [leaf.poses[i] for i in np.argsort(-pmu)[:100]]
        nn = leaf.add_poses(new + pd.neighbours(heavy, t))
        ndp = 0
        if a.pose_max and len(leaf.poses) > a.pose_max:
            rc_all = T4.pose_rc(leaf, py, plam)
            keep = np.zeros(len(leaf.poses), bool)
            keep[:len(pmu)] = pmu > 1e-9
            room = a.pose_max - int(keep.sum())
            if room > 0:
                preg = np.array(leaf.reg)
                regs = sorted(set(int(v) for v in preg))
                q = max(room // max(len(regs), 1), 1)
                for rr in regs:
                    w = np.nonzero((preg == rr) & ~keep)[0]
                    if len(w):
                        keep[w[np.argsort(-rc_all[w])[:q]]] = True
                room2 = a.pose_max - int(keep.sum())
                if room2 > 0:
                    w = np.nonzero(~keep)[0]
                    keep[w[np.argsort(-rc_all[w])[:room2]]] = True
            ndp = leaf.sift_poses(keep)
        log(f'   stage {stage}: +{nn} columns, -{ndp} sifted -> {len(leaf.poses)} columns for '
            f'{len(set(leaf.gkey))} poses (pricing gap {gap:+.6f})')
        if nn == 0:
            log('   [no improving column]')
            break

    log('RESULT tag=%s t=%s corners=%s chord=%d ' % (a.TAG, t, a.corners, int(a.chord))
        + ' '.join((f'{k}=LP{v["LP"]:.6f}/M{v["M"]:.6f}/K{v["kmax"]:.4f}'
                    + (f'/pure{v["LP_pure"]:.6f}' if 'LP_pure' in v else '')
                    + (f'/nochord{v["LP_nochord"]:.6f}' if 'LP_nochord' in v else '')
                    + (f'/nolines{v["LP_nolines"]:.6f}' if 'LP_nolines' in v else '')
                    + f'/int{v["interior"]:.6f}') if v else f'{k}=FAIL'
                   for k, v in results.items())
        + f' cols={len(leaf.poses)} rows={len(leaf.pts)} cq={len(leaf.cliques)} '
          f'secs={time.time()-T0:.0f}')


if __name__ == '__main__':
    main()
