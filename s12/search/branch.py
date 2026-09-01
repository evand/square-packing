#!/usr/bin/env python3
"""branch.py -- corner-branch certificates: the LP of tighten.py with Lagrange multipliers.

The pure cover LP cannot prove s(12) >= t for t >= 3.99 (search/DUAL_EXACT.md: a fractional
packing of mass 12.008 exists there).  That packing puts mass 0.85 in each corner; a real packing
puts 0 or 1 square there.  Branch on it.

Region: the four corner boxes  B_1 = [0,r]^2,  B_2 = [s-r,s]x[0,r],  B_3 = [0,r]x[s-r,s],
B_4 = [s-r,s]^2  (closed; a pose belongs to B_j iff its CENTRE does), r < 1/2 + 1/sqrt2 so that at
most one square of a packing is centred in each box.  Two forms:

  single (--k K, K in 0..4): one multiplier lambda, K = number of squares centred in the union;
        w(S) >= 1 + lambda (centre in a box),  w(S) >= 1 (otherwise)  =>  n + lambda*K <= W.
  per box (--k K1K2K3K4, K_j in {0,1}): multipliers lambda_1..lambda_4, K_j = 1 iff a square is
        centred in B_j;   w(S) >= 1 + lambda_j (centre in B_j)  =>  n + sum_j lambda_j K_j <= W.

So W - sum lambda*K < n refutes the leaf; the 5 (single) or 16 (per box; 6 up to symmetry) leaves
together refute every packing.  The verifier (verify/) checks a certificate with the trailer
`region corner r_num r_den / lambda L [L2 L3 L4] / k K [K2 K3 K4]` exactly as a plain one.

LP = tighten.py's cutting-plane loop (exact verifier as separation oracle, column generation on
the dual) with the multipliers as free columns:

        min  sum_k |orbit_k| x_k  -  sum_j K_j lambda_j
        s.t. A x - lambda_{box(r)} >= 1 + margin  (rows in a box),   A x >= 1 + margin  (others),  x >= 0.

Rows whose box has K_j = 0 are dropped and lambda_j is fixed at -2 (the box is free); otherwise
-1 <= lambda_j <= 4.  Row flags come from the verifier's witnesses (box index) or, for the lattice
warm start, from the centre.  Single-lambda leaves keep D4-symmetric point sets (orbit columns,
angles [0,45deg]); per-box leaves are not symmetric in general, so every point is its own
column, pricing runs over the whole container and the verifier sweeps [0,90deg).

Usage
    python3 search/branch.py CERT TAG --k 4    [--r 1] [--Dp 784000 --mul 199] [--N 6000] [--colgen 30]
    python3 search/branch.py CERT TAG --k 1100 ...          -> runs/branch_TAG.txt (+ .log, .json)
"""
import sys, os, math, time, argparse, json
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tighten as T
import nu_f as NF
import anchorclique as AC
import anchorsep as ASEP
import floatsep as FS

LAM_FREE = -2.0          # multiplier of a box with K_j = 0: its rows are vacuous (A x >= -1 + margin)
SOLVER = os.environ.get('BRANCH_SOLVER', 'highs-ipm')   # scipy method: 'highs', 'highs-ds', 'highs-ipm'; 'warm' = highspy dual simplex with a carried basis (measured: a 40k-iteration phase 1 per round, no gain -- search/LPSPEED.md); 'restricted' = column-sifted master (see solve_restricted; the measured win)
RESTRICTED_ADD = int(os.environ.get('BRANCH_RESTRICTED_ADD', '1000'))      # columns priced into the master per pass
RESTRICTED_PASSES = int(os.environ.get('BRANCH_RESTRICTED_PASSES', '2'))  # passes per round (0 = to convergence: exact full-LP optimum)
RESTRICTED_TOL = float(os.environ.get('BRANCH_RESTRICTED_TOL', '1e-7'))
try:
    import highspy
except ImportError:
    highspy = None
    if SOLVER == 'warm': SOLVER = 'highs'
LAM_LO, LAM_HI = -1.0, 1.5   # the LP has a free direction (corner-only atoms up, lambda up); the cap keeps certificates sane (--lam-hi overrides; the k = 4 leaf sits at the cap)


class BModel(T.Model):
    """tighten.Model plus a box flag per row (0 = none, j = box j) and the multiplier columns."""

    def __init__(self, K, D, orbits_int, r, kreg, sym):
        super().__init__(K, D, orbits_int)
        # anchor-clique columns (search/anchorclique.py), appended after the point columns:
        # column n_orb + j is the D4 orbit of clique j, of cost = number of distinct images and
        # entry (number of images containing the row's pose) -- exactly like a point orbit.
        self.cliques = []; self.cparams = []; self.ckey = {}
        self.csizes = np.zeros(0); self.CR = []; self.CC = []; self.CV = []
        self.r = float(r); self.rfrac = Fraction(r); self.rflag = []; self.rowkeys = []; self.sym = sym
        self.row_grid = (1e-7, 1e-8)     # (position, angle) rounding of the row key: rows closer than this are one row
        if isinstance(kreg, (tuple, list)):
            self.perbox = True; self.kvec = [int(v) for v in kreg]; assert len(self.kvec) == 4 and all(v in (0, 1) for v in self.kvec)
        else:
            self.perbox = False; self.kvec = [int(kreg)] * 4; assert 0 <= int(kreg) <= 4
        self.nl = 4 if self.perbox else 1

    def ncols(self):
        return len(self.orbits) + len(self.cliques)

    def costv(self):
        """objective coefficients of all columns: orbit sizes, then clique-orbit sizes"""
        return np.concatenate([self.sizes, self.csizes])

    def matrix(self):
        A = super().matrix()
        if not self.cliques: return A
        nr = len(self.rows); nc = len(self.cliques)
        if self.CR:
            R = np.concatenate(self.CR); C = np.concatenate(self.CC); V = np.concatenate(self.CV)
            B = sp.coo_matrix((V, (R, C)), shape=(nr, nc)).tocsr()
        else:
            B = sp.csr_matrix((nr, nc))
        return sp.hstack([A, B]).tocsr()

    def add_clique(self, cl, params=None):
        """add the D4 orbit of an anchor clique as one column (single clique if the model is not
        symmetric); returns False if it is already there"""
        k = AC.key(cl)
        if k in self.ckey: return False
        sfr = Fraction(self.K, self.D)
        imgs = AC.images(sfr, cl) if self.sym else [cl]
        j = len(self.cliques); self.ckey[k] = j
        self.cliques.append(imgs); self.cparams.append(params)
        self.csizes = np.append(self.csizes, float(len(imgs)))
        if self.rows:
            Q = np.array(self.rows); v = AC.coeff(imgs, Q)
            hit = np.nonzero(v)[0]
            if hit.size:
                self.CR.append(hit.astype(np.int32)); self.CC.append(np.full(hit.size, j, dtype=np.int32)); self.CV.append(v[hit])
        return True

    def box_of(self, cx, cy, eps=1e-9):
        s, r = self.s, self.r
        lo_x, hi_x = cx <= r + eps, cx >= s - r - eps
        lo_y, hi_y = cy <= r + eps, cy >= s - r - eps
        if lo_x and lo_y: return 1
        if hi_x and lo_y: return 2
        if lo_x and hi_y: return 3
        if hi_x and hi_y: return 4
        return 0

    def row_active(self, flag):
        if flag == 0: return True
        return self.kvec[flag - 1] != 0 if self.perbox else self.kvec[0] != 0

    def add_rows(self, pl):
        """pl: iterable of (cx, cy, theta, h) or (cx, cy, theta, h, flag).  Returns number added."""
        P, own = self.P, self.own; added = 0
        for row in pl:
            cx, cy, tm, h = row[:4]
            flag = int(row[4]) if len(row) > 4 else self.box_of(cx, cy)
            if not self.row_active(flag): continue
            gp, ga = self.row_grid
            key = (round(cx / gp), round(cy / gp), round(tm / ga), round(h * 1e7), flag)
            if key in self.rkey: continue
            self.rkey.add(key)
            ct, st = math.cos(tm), math.sin(tm)
            u0 = cx * ct + cy * st; u1 = -cx * st + cy * ct
            q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
            ok = np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= h - T.TOL
            cnt = np.bincount(own[ok], minlength=len(self.orbits)); idx = np.nonzero(cnt)[0]
            r = len(self.rows); self.rows.append((cx, cy, tm, h)); self.rflag.append(flag); self.rowkeys.append(key)
            if idx.size:
                self.R.append(np.full(idx.size, r, dtype=np.int32)); self.C.append(idx.astype(np.int32)); self.V.append(cnt[idx].astype(float))
            added += 1
        if self.cliques and added:
            n0 = len(self.rows) - added
            Q = np.array(self.rows[n0:]); rot = AC._rot(Q)
            for j, imgs in enumerate(self.cliques):
                v = AC.coeff(imgs, Q, rot); hit = np.nonzero(v)[0]
                if hit.size:
                    self.CR.append((hit + n0).astype(np.int32)); self.CC.append(np.full(hit.size, j, dtype=np.int32)); self.CV.append(v[hit])
        return added

    def add_orbit(self, X, Y):
        if self.sym: return super().add_orbit(X, Y)
        K = self.K; X = min(max(int(X), 0), K); Y = min(max(int(Y), 0), K)
        key = (X, Y)
        if key in self.okey: return False
        k = len(self.orbits); self.okey[key] = k
        oi = np.array([[X, Y]], dtype=np.int64); self.orb_int.append(oi); of = oi / self.D; self.orbits.append(of)
        self.sizes = np.append(self.sizes, 1.0)
        self.P = np.concatenate([self.P, of]); self.own = np.concatenate([self.own, np.full(1, k)])
        if self.rows:
            Q = np.array(self.rows); ct = np.cos(Q[:, 2]); st = np.sin(Q[:, 2])
            u0 = Q[:, 0] * ct + Q[:, 1] * st; u1 = -Q[:, 0] * st + Q[:, 1] * ct; h = Q[:, 3]
            px, py = of[0]
            q0 = px * ct + py * st; q1 = -px * st + py * ct
            hit = np.nonzero(np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= h - T.TOL)[0]
            if hit.size:
                self.R.append(hit.astype(np.int32)); self.C.append(np.full(hit.size, k, dtype=np.int32)); self.V.append(np.ones(hit.size))
        return True

    def prune(self, keep):
        if self.CR:
            idx = np.full(len(self.rows), -1, dtype=np.int64); idx[np.nonzero(keep)[0]] = np.arange(int(keep.sum()))
            R = np.concatenate(self.CR); C = np.concatenate(self.CC); V = np.concatenate(self.CV)
            nr = idx[R]; sel = nr >= 0
            self.CR = [nr[sel].astype(np.int32)]; self.CC = [C[sel]]; self.CV = [V[sel]]
        super().prune(keep)
        self.rflag = [f for f, k in zip(self.rflag, keep) if k]
        self.rowkeys = [rk for rk, k in zip(self.rowkeys, keep) if k]
        self.rkey = set(self.rowkeys)

    def flagmat(self):
        """rows x nl matrix F with F[r, j] = 1 iff row r is charged multiplier j"""
        nr = len(self.rows); fl = np.array(self.rflag, dtype=int)
        if self.perbox:
            rr = np.nonzero(fl)[0]
            return sp.csr_matrix((np.ones(len(rr)), (rr, fl[rr] - 1)), shape=(nr, 4))
        return sp.csr_matrix((fl > 0).astype(float).reshape(-1, 1))

    def kcoef(self):
        return np.array(self.kvec, dtype=float) if self.perbox else np.array([float(self.kvec[0])])

    def solve(self, margin, cost=None, budget=None, fixed_zero=None):
        """min cost.x - K.lam  s.t.  A x - F lam >= 1+margin, x >= 0 [, sizes.x <= budget] [, x_k = 0].
        Returns (obj, x, lam, y)."""
        plain = cost is None and budget is None and not (fixed_zero is not None and len(fixed_zero))
        if SOLVER == 'warm' and plain:
            return self.solve_warm(margin)
        if SOLVER == 'restricted' and plain:
            return self.solve_restricted(margin)
        A = self.matrix(); n = self.ncols(); nr = len(self.rows); nl = self.nl
        c = np.concatenate([self.costv() if cost is None else np.asarray(cost, dtype=float), -self.kcoef()])
        Aub = sp.hstack([-A, self.flagmat()]).tocsr(); bub = -np.full(nr, 1.0 + margin)
        extra_A = []; extra_b = []
        if budget is not None:
            extra_A.append(np.concatenate([self.costv(), np.zeros(nl)]).reshape(1, -1)); extra_b.append(budget)
        if fixed_zero is not None and len(fixed_zero):
            Z = sp.csr_matrix((np.ones(len(fixed_zero)), (np.arange(len(fixed_zero)), np.asarray(fixed_zero))), shape=(len(fixed_zero), n + nl))
            extra_A.append(Z); extra_b += [0.0] * len(fixed_zero)
        if extra_A:
            Aub = sp.vstack([Aub] + [sp.csr_matrix(a) for a in extra_A]).tocsr(); bub = np.concatenate([bub, np.array(extra_b)])
        kc = self.kcoef()
        bounds = [(0, None)] * n + [((LAM_FREE, LAM_FREE) if kc[j] == 0 else (LAM_LO, LAM_HI)) for j in range(nl)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c=c, A_ub=Aub, b_ub=bub, bounds=bounds, method=SOLVER, options=dict(T.HIGHS))
        if not res.success: return None
        y = -res.ineqlin.marginals[:nr]
        return res.fun, res.x[:n], np.array(res.x[n:n + nl]), np.maximum(y, 0.0)


def _solve_warm(self, margin):
    """the same LP through highspy, dual simplex, with the basis of the previous solve carried over:
    rows are only ever appended (new rows start basic, i.e. their slack is in the basis) or deleted
    (prune, handled by index bookkeeping), columns only appended (new columns start nonbasic at 0),
    so the previous basis is a valid, near-optimal starting basis and the re-solve is cheap."""
    A = self.matrix(); n = self.ncols(); nr = len(self.rows); nl = self.nl
    M = sp.hstack([A, -self.flagmat()]).tocsc()          # A x - F lam >= 1 + margin
    kc = self.kcoef()
    lp = highspy.HighsLp()
    lp.num_col_ = n + nl; lp.num_row_ = nr
    lp.col_cost_ = np.concatenate([self.costv(), -kc])
    lo = np.concatenate([np.zeros(n), [LAM_FREE if kc[j] == 0 else LAM_LO for j in range(nl)]])
    hi = np.concatenate([np.full(n, highspy.kHighsInf), [LAM_FREE if kc[j] == 0 else LAM_HI for j in range(nl)]])
    lp.col_lower_ = lo; lp.col_upper_ = hi
    lp.row_lower_ = np.full(nr, 1.0 + margin); lp.row_upper_ = np.full(nr, highspy.kHighsInf)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = M.indptr.astype(np.int32); lp.a_matrix_.index_ = M.indices.astype(np.int32); lp.a_matrix_.value_ = M.data.astype(float)
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'simplex'); h.setOptionValue('simplex_strategy', 1)
    h.setOptionValue('random_seed', int(T.HIGHS.get('random_seed', 0)))
    h.passModel(lp)
    # basis from the previous solve: row statuses are keyed by row key so that pruning is handled
    prev = getattr(self, '_basis', None)
    if prev is not None:
        pcol, prow, pn = prev
        b = highspy.HighsBasis(); b.valid = True
        cs = [highspy.HighsBasisStatus.kLower] * (n + nl)
        for j in range(min(pn, n)): cs[j] = pcol[j]
        for j in range(nl): cs[n + j] = pcol[pn + j] if pn + j < len(pcol) else highspy.HighsBasisStatus.kLower
        rs = [prow.get(self.rowkeys[i], highspy.HighsBasisStatus.kBasic) for i in range(nr)]
        b.col_status = cs; b.row_status = rs
        h.setBasis(b)
    h.run()
    st = h.getModelStatus()
    if st != highspy.HighsModelStatus.kOptimal:
        return None
    sol = h.getSolution(); bas = h.getBasis()
    x = np.array(sol.col_value); y = np.array(sol.row_dual)
    self._basis = (list(bas.col_status), {self.rowkeys[i]: bas.row_status[i] for i in range(nr)}, n)
    # HiGHS: for a minimisation the dual of a row at its lower bound is >= 0
    return float(h.getInfo().objective_function_value), x[:n], x[n:n + nl], np.maximum(y, 0.0)
BModel.solve_warm = _solve_warm


def _solve_restricted(self, margin, log=None, no_cliques=False):
    """column-sifted master (search/LPSPEED.md).  The IPM's cost per iteration is ~ (columns)^1.5 here
    and only ~10 % of the columns ever carry weight, so: solve the LP over an ACTIVE column set (the
    previous support, this round's new columns, and whatever was priced in) with scipy highs-ipm, price
    every column of the model by its reduced cost c_j - y.A_j under the master's dual, add the most
    negative RESTRICTED_ADD, repeat for RESTRICTED_PASSES passes (0 = until none is below
    -RESTRICTED_TOL, which is the exact full-LP optimum).  Rows are all kept.  The active set persists
    across rounds in self._active (points) and self._active_cl (cliques); columns that leave the
    support stay active until the set is trimmed (trim: drop zero-weight columns not added in the
    last 3 rounds).  Clique columns are sifted like point columns since 2026-09-01: they are DENSE
    (a clique holds thousands of poses) and a master carrying 1,500 of them took an hour per solve
    (runs/branch_J16h.log), so only the ones in the support or with negative reduced cost are in;
    a newly separated clique is active for its first solves.  The dual y and the multipliers are
    those of the master, which is what the loop's pricing and verifier need."""
    A = self.matrix(); n = self.ncols(); npt = len(self.orbits); ncl = len(self.cliques)
    nr = len(self.rows); nl = self.nl
    kc = self.kcoef(); c_full = np.concatenate([self.costv(), -kc])
    lo = np.concatenate([np.zeros(n), [LAM_FREE if kc[j] == 0 else LAM_LO for j in range(nl)]])
    hi = np.concatenate([np.full(n, np.inf), [LAM_FREE if kc[j] == 0 else LAM_HI for j in range(nl)]])
    M = sp.hstack([A, -self.flagmat()]).tocsr(); Mt = M.T.tocsr(); b = np.full(nr, 1.0 + margin)
    def grow(mask, m, default_new, default_all):
        if mask is None or len(mask) != m:
            new = np.full(m, default_all if mask is None else default_new, dtype=bool)
            if mask is not None: new[:len(mask)] = mask
            return new
        return mask
    def grow_age(age, m):
        if age is None or len(age) != m:
            a2 = np.zeros(m, dtype=np.int64)
            if age is not None: a2[:len(age)] = age
            return a2
        return age
    act = grow(getattr(self, '_active', None), npt, True, n <= 3000)   # start empty if large: greedy cover + pricing fill it
    age = grow_age(getattr(self, '_age', None), npt)
    actc = grow(getattr(self, '_active_cl', None), ncl, True, True)     # new cliques are in for their first solves
    agec = grow_age(getattr(self, '_age_cl', None), ncl)
    # feasibility: every row needs an active column with a positive entry
    Pb = A[:, :npt].copy(); Pb.data = (Pb.data > 0).astype(float)
    fixed = (A[:, npt:] @ actc.astype(float)) if ncl else np.zeros(nr)
    unc = (Pb @ act.astype(float) + fixed) <= 0
    while unc.any():
        hits = Pb[unc].sum(axis=0).A1; j = int(np.argmax(hits))
        if hits[j] <= 0: break
        act[j] = True; unc &= (Pb[:, j].toarray().ravel() <= 0)
    passes = 0; t0 = time.time(); res = None
    if no_cliques: act = act.copy()          # the matched pure solve must not disturb the master
    while True:
        cols = np.concatenate([np.nonzero(act)[0], np.zeros(0, dtype=np.int64) if no_cliques else npt + np.nonzero(actc)[0],
                               np.arange(n, n + nl)])
        Mc = M[:, cols].tocsr()
        bounds = [(l, None if not np.isfinite(h) else h) for l, h in zip(lo[cols], hi[cols])]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c=c_full[cols], A_ub=-Mc, b_ub=-b, bounds=bounds, method='highs-ipm', options=dict(T.HIGHS))
        passes += 1
        if not res.success: return None
        # the matched pure solve must use exactly the master's own column set, or it prices in
        # columns the clique solve never saw and comes out LOWER (measured: -0.09)
        if no_cliques:
            if os.environ.get('CQDBG'): print(f"    [matched] pure cols={len(cols)-nl} obj={res.fun:.6f}", flush=True)
            break
        y = np.maximum(-res.ineqlin.marginals, 0.0)
        rc = c_full - (Mt @ y); rc[cols] = np.inf
        neg = np.nonzero(rc[:npt] < -RESTRICTED_TOL)[0]
        negc = np.nonzero(rc[npt:n] < -RESTRICTED_TOL)[0]
        if os.environ.get('CQDBG') and not no_cliques: print(f"    [main] pass {passes} cols={len(cols)-nl} obj={res.fun:.6f}", flush=True)
        if log: log(f"    restricted pass {passes}: {len(cols)-nl} cols ({int(actc.sum())} cliques) obj={res.fun:.7f} neg_rc={len(neg)}+{len(negc)}cl min_rc={rc[:n].min() if (len(neg) or len(negc)) else 0:+.2e} t={time.time()-t0:.0f}s")
        if (len(neg) == 0 and len(negc) == 0) or (RESTRICTED_PASSES and passes >= RESTRICTED_PASSES): break
        pick = neg[np.argsort(rc[neg])[:RESTRICTED_ADD]]; act[pick] = True; age[pick] = 0
        actc[negc] = True; agec[negc] = 0                  # cliques are few: every improving one comes in
    x = np.zeros(n); x[cols[:-nl]] = res.x[:-nl]; lam = np.array(res.x[-nl:])
    if no_cliques: return res.fun, x, lam, np.maximum(-res.ineqlin.marginals, 0.0)
    # trim: zero-weight columns that have been idle for 3 rounds leave the master (they stay in the model and are priced every round)
    xp = x[:npt]; xc = x[npt:n]
    age[act & (xp <= 1e-12)] += 1; age[xp > 1e-12] = 0
    act &= ~((xp <= 1e-12) & (age > 3))
    agec[actc & (xc <= 1e-12)] += 1; agec[xc > 1e-12] = 0
    actc &= ~((xc <= 1e-12) & (agec > 3))
    self._active = act; self._age = age; self._active_cl = actc; self._age_cl = agec
    self._restricted_info = dict(passes=passes, cols=int(act.sum()) + int(actc.sum()), neg_rc=int(len(neg) + len(negc)), t=time.time() - t0)
    return res.fun, x, lam, y
BModel.solve_restricted = _solve_restricted


def parse_k(txt):
    """'4' -> 4 (single), '1100' -> (1,1,0,0) (per box)"""
    txt = str(txt).strip()
    if len(txt) == 4 and set(txt) <= set('01'): return tuple(int(ch) for ch in txt)
    return int(txt)


def build_model(cert, r, kreg, Dp=None, mul=1, raw_points=False):
    """raw_points (asymmetric leaves only): the certificate's points are the columns, no D4 images --
    for restarting a leaf from its own probe/checkpoint, whose points already are the column set."""
    sn, sd, D, WD, rows = T.read_cert(cert)
    K = Fraction(sn, sd) * D; assert K.denominator == 1; K = int(K)
    if Dp is None: Dp = D * mul
    K2 = K * mul
    sym = not isinstance(kreg, tuple)
    okey = {}; orbs = []; w0 = []
    for X, Y, W in rows:
        X = int(X) * mul; Y = int(Y) * mul
        if raw_points and not sym:
            p = (X, Y)
            if p not in okey: okey[p] = len(orbs); orbs.append(np.array([p], dtype=np.int64)); w0.append(0.0)
            w0[okey[p]] += W / WD
            continue
        imgs = sorted(set([(X, Y), (K2 - X, Y), (X, K2 - Y), (K2 - X, K2 - Y), (Y, X), (K2 - Y, X), (Y, K2 - X), (K2 - Y, K2 - X)]))
        if sym:
            key = imgs[0]
            if key not in okey:
                okey[key] = len(orbs); orbs.append(np.array(imgs, dtype=np.int64)); w0.append(W / WD)
        else:
            # every image is its own column; a symmetric input is symmetrised so that all images exist
            # as columns, but only the point actually listed carries its weight (images start at 0)
            for p in imgs:
                if p not in okey:
                    okey[p] = len(orbs); orbs.append(np.array([p], dtype=np.int64)); w0.append(0.0)
            w0[okey[(X, Y)]] += W / WD
    m = BModel(K2, Dp, orbs, r, kreg, sym)
    m.okey = okey
    return m, np.array(w0), Fraction(K2, Dp)


def export(m, x, lam, path, WD=10 ** 7, factor=1.0, up=True):
    """certificate at container K/D with the region trailer; weights x_k*factor over WD rounded up
    (final) or down (probe); lambda*factor rounded the other way (down for the final file, up for
    the probe -- the probe must be at least as demanding as the true (x, lambda))."""
    s = Fraction(m.K, m.D)
    lines = []; tot = 0
    for k, o in enumerate(m.orb_int):
        if x[k] <= 0: continue
        v = x[k] * factor * WD
        w = int(math.ceil(v - 1e-9)) if up else int(math.floor(v))
        if w <= 0: continue
        for X, Y in o: lines.append(f"{X} {Y} {w}"); tot += w
    cls = []; clw = []                                    # one file clique per image, same weight
    n_orb = len(m.orb_int)
    for j, imgs in enumerate(getattr(m, 'cliques', [])):
        v = x[n_orb + j] * factor * WD
        w = int(math.ceil(v - 1e-9)) if up else int(math.floor(v))
        if w <= 0: continue
        for im in imgs: cls.append(im); clw.append(w); tot += w
    # Every anchor point is emitted as a ZERO-WEIGHT atom.  It carries no weight, but the sweep's
    # cells are the atoms' breakpoints and a cell is credited only if it lies WHOLLY inside a piece:
    # without them the cells straddle the boundary of {S : p in S} and of {S : A subseteq S} and a
    # band around each loses the credit -- which is exactly where the covering is tight.
    extra = set()
    for cl in cls:
        for anc in cl[0]:
            for (vx, vy) in ([(anc[1], anc[2])] if anc[0] == 'P' else [(anc[1], anc[2]), (anc[3], anc[4])]):
                gx = Fraction(vx) * m.D; gy = Fraction(vy) * m.D
                if gx.denominator == 1 and gy.denominator == 1: extra.add((int(gx), int(gy)))
    for (X, Y) in sorted(extra): lines.append(f"{X} {Y} 0")
    Ls = [int(math.floor(l * factor * WD)) if up else int(math.ceil(l * factor * WD)) for l in lam]
    with open(path, 'w') as f:
        f.write(f"{s.numerator} {s.denominator}\n{m.D}\n{WD}\n{len(lines)}\n" + "\n".join(lines) + "\n")
        if cls: f.write(AC.block(cls, clw, m.D))
        f.write(f"region corner {m.rfrac.numerator} {m.rfrac.denominator}\nlambda {' '.join(str(L) for L in Ls)}\nk {' '.join(str(v) for v in (m.kvec if m.perbox else [m.kvec[0]]))}\n")
    return tot / WD, len(lines), np.array(Ls) / WD


def read_witnesses(path, N):
    """verifier witness file (v theta cx cy flag) -> rows (cx, cy, theta_k, sigma_k/2, flag)"""
    out = []
    for line in open(path):
        q = line.split()
        if len(q) not in (4, 5): continue
        v, th, cx, cy = map(float, q[:4]); fl = int(q[4]) if len(q) == 5 else 0
        k = int(round(N * math.tan(th / 2)))
        out.append((v, th, cx, cy, k, fl))
    out.sort()
    return [(cx, cy, th, T.sigma_k(N, k) / 2, fl) for v, th, cx, cy, k, fl in out]


def lattice_rows(m, w, eta=0.01, dt=0.01, thr=1.05, perang=400, B=6, nproc=8):
    """closed unit squares on an (eta, dt) lattice capturing < thr under w; angles [0,45deg] for a
    symmetric model, [0,90deg) otherwise (tighten.lattice_rows covers [0,45deg] only)."""
    import multiprocessing as mp
    top = math.pi / 4 if m.sym else math.pi / 2
    thetas = np.linspace(0.0, top, int(round(top / dt)) + 1)
    if not m.sym: thetas = thetas[:-1]
    pool = mp.get_context('fork').Pool(nproc)
    gmin, nviol, pend = NF.separate(m.P, w, m.s, thetas, eta, perang=perang, B=B, pool=pool, thr=thr)
    pool.close(); pool.join()
    s = m.s; out = []
    for cx, cy, tm in pend:
        wid = abs(math.cos(tm)) + abs(math.sin(tm))
        if wid / 2 - 1e-12 <= cx <= s - wid / 2 + 1e-12 and wid / 2 - 1e-12 <= cy <= s - wid / 2 + 1e-12:
            out.append((cx, cy, tm, 0.5))
    return out


def corner_rows(m, eta=0.05, dt=0.05):
    """Poses whose CENTRE lies in a corner box, on an (eta, dt) lattice, whatever they capture
    (task K deliverable 4a).  The lattice warm start of `lattice_rows` keeps only the poses that
    capture least under the input weights, and on the k = 4 leaves that is never a corner pose:
    round 0 therefore has no row charged to any multiplier, lambda goes straight to its cap and
    the first LP value is meaningless.  These rows are ordinary unit-square constraints
    (h = 1/2, admissible centre), so they are valid whatever the LP does with them."""
    s = m.s; r = m.r; out = []
    boxes = [(0.0, 0.0), (s - r, 0.0), (0.0, s - r), (s - r, s - r)]
    # A symmetric model has D4 orbits as columns, so a pose and its D4 image have the SAME
    # coefficient vector: box 0 with theta in [0, 45deg] already generates every distinct row
    # (the map (x,y) -> (y,x) fixes box 0 and sends theta to 90deg - theta).
    top = math.pi / 4 if m.sym else math.pi / 2
    thetas = np.arange(0.0, top + (dt if m.sym else 0.0) * 0.5, dt)
    if m.sym: boxes = boxes[:1]
    for tm in thetas:
        wid = abs(math.cos(tm)) + abs(math.sin(tm)); lo = wid / 2; hi = s - wid / 2
        for j, (bx, by) in enumerate(boxes):
            if not m.row_active(j + 1): continue
            x0 = max(lo, bx); x1 = min(hi, bx + r); y0 = max(lo, by); y1 = min(hi, by + r)
            if x1 < x0 or y1 < y0: continue
            gx = np.arange(x0, x1 + 1e-12, eta); gy = np.arange(y0, y1 + 1e-12, eta)
            if len(gx) == 0 or len(gy) == 0: continue
            if gx[-1] < x1 - 1e-12: gx = np.append(gx, x1)
            if gy[-1] < y1 - 1e-12: gy = np.append(gy, y1)
            for cx in gx:
                for cy in gy: out.append((float(cx), float(cy), float(tm), 0.5, j + 1))
    return out


def seed_rows(path):
    """rows of a sibling leaf's dual checkpoint (runs/branch_TAG_dual.txt: `cx cy theta h flag y`,
    written every round) -- the binding row set of a converged leaf, which is what makes a sibling
    start where the sibling finished instead of at the lattice warm start."""
    out = []
    for line in open(path):
        if line.startswith('#'): continue
        q = line.split()
        if len(q) < 5: continue
        out.append((float(q[0]), float(q[1]), float(q[2]), float(q[3]), int(q[4])))
    return out


def price_asym(m, y, pitch=0.02, want=200, ysup=1e-9):
    """reduced-cost pricing without symmetry: the dual y is a packing measure on the rows; the
    reduced cost of the point p is 1 - cov(p), cov(p) = sum_r y_r [p in Q_r].  Evaluates cov on a
    pitch-grid of the whole container, refines the best candidates on the integer grid, returns
    integer points with cov > 1, best first."""
    Q = np.array(m.rows); sup = y > ysup; Q = Q[sup]; yy = y[sup]
    if len(yy) == 0: return []
    s = m.s; cx0, cy0, t0, h0 = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
    ct_, st_ = np.cos(t0), np.sin(t0)
    R = np.c_[cx0 * ct_ + cy0 * st_, -cx0 * st_ + cy0 * ct_, ct_, st_, h0, yy]
    def cov_at(G):
        out = np.zeros(len(G))
        for i in range(0, len(R), 2000):
            B = R[i:i + 2000]
            q0 = G[:, 0:1] * B[:, 2] + G[:, 1:2] * B[:, 3]; q1 = -G[:, 0:1] * B[:, 3] + G[:, 1:2] * B[:, 2]
            ins = np.maximum(np.abs(q0 - B[:, 0]), np.abs(q1 - B[:, 1])) <= B[:, 4] + 1e-9
            out += ins @ B[:, 5]
        return out
    n0 = int(math.ceil(s / pitch)); gx = np.arange(n0 + 1) * (s / n0)
    G0, G1 = np.meshgrid(gx, gx, indexing='ij'); G = np.c_[G0.ravel(), G1.ravel()]
    c = np.zeros(len(G))
    for i in range(0, len(G), 20000): c[i:i + 20000] = cov_at(G[i:i + 20000])
    o = np.argsort(c)[::-1][:want]
    cand = []; D = m.D; step = s / n0
    for i in o:
        if c[i] <= 1.0 + 1e-7: break
        px, py = G[i]
        X0, Y0 = int(round(px * D)), int(round(py * D)); r_ = int(math.ceil(step * D))
        xs = np.arange(max(0, X0 - r_), min(m.K, X0 + r_) + 1); ys_ = np.arange(max(0, Y0 - r_), min(m.K, Y0 + r_) + 1)
        if len(xs) * len(ys_) > 4000:
            xs = xs[::max(1, len(xs) // 60)]; ys_ = ys_[::max(1, len(ys_) // 60)]
        GX, GY = np.meshgrid(xs, ys_, indexing='ij'); Gf = np.c_[GX.ravel(), GY.ravel()] / D
        cf = cov_at(Gf); j = int(np.argmax(cf))
        if cf[j] > 1.0 + 1e-7: cand.append((float(cf[j]), int(GX.ravel()[j]), int(GY.ravel()[j])))
    cand.sort(reverse=True)
    return cand


CQ_MAX = int(os.environ.get('BRANCH_CQ_MAX', '400'))     # cap on the number of clique-orbit columns
INT_TOP = int(os.environ.get('BRANCH_CQ_INT_TOP', '150'))  # anchor points for the interior separator


def price_cliques(m, y, x, want=40, thr=1e-7, pitch=0.01, top=600, fracs=(0.95, 0.85, 0.7, 0.55, 0.4, 0.25), log=None, interior=False):
    """Price anchor-clique columns `K(p, A)` on the current dual (search/anchorsep.py).

    The dual `y` is a fractional packing on the rows (D4-averaged for a symmetric model), so the
    reduced cost of the clique orbit is `|orbit| (1 - ybar(K))` and a column is worth adding iff
    `ybar(K) > 1` -- which is the same statement as "the packing violates the clique constraint
    `mu(K) <= 1`".  Candidates are wall points `p` (Lemma 1: nowhere else can a clique beat the
    coverage row at `p`) on a grid, ranked by their coverage, each with the Lemma-2 anchors for
    several `eps`; `K(p, A)` then contains the whole point clique of `p`, so the column dominates
    the point column of `p` and what it adds is the dual mass of `{S : A subseteq S}` outside
    `P_p`.  Returns (number added, best ybar seen, best ybar of the corresponding point row)."""
    want = min(want, CQ_MAX - len(m.cliques))
    if want <= 0 or not len(m.rows): return 0, None, None
    Q = np.array(m.rows); sup = y > 1e-9
    if not sup.any(): return 0, None, None
    best = ASEP.separate(Q[sup], y[sup], float(Fraction(m.K, m.D)), m.D, pitch=pitch, top=top,
                         fracs=fracs, sym=m.sym)
    if interior:
        # the interior family as well (search/RECONCILE.md): the largest violations of these
        # packings are at wall distance > 1, which `separate` never looks at
        best = best + ASEP.separate_interior(Q[sup], y[sup], float(Fraction(m.K, m.D)), m.D,
                                             pitch=max(pitch, 0.02), top=INT_TOP, sym=m.sym)
        best.sort(key=lambda t: -t[0])
    added = 0
    for (mk, mp, params, cl) in best:
        if mk <= 1.0 + thr or added >= want: break
        if AC.key(cl) in m.ckey: continue
        if m.add_clique(cl, params): added += 1
    return added, (best[0][0] if best else None), (best[0][1] if best else None)


def save_cliques(m, path):
    """checkpoint the anchor-clique columns as their generating parameters (exact integers)"""
    with open(path + ".tmp", 'w') as f:
        f.write(f"{m.K} {m.D} {int(m.sym)}\n")
        for p in m.cparams:
            if p is not None: f.write(" ".join(str(int(v)) for v in p) + "\n")
    os.replace(path + ".tmp", path)


def load_cliques(m, path):
    n = 0
    with open(path) as f:
        K, D, sym = (int(v) for v in f.readline().split())
        assert K == m.K and D == m.D, f"clique file {path} is for container {K}/{D}, model is {m.K}/{m.D}"
        sfr = Fraction(m.K, m.D)
        for line in f:
            q = [int(v) for v in line.split()]
            if len(q) == 5:                       # wall-perpendicular Lemma-2 anchor
                cl = AC.kpa(sfr, m.D, q[0], q[1], q[2], q[3], q[4])
            elif len(q) == 6:                     # general segment anchor (search/RECONCILE.md)
                cl = AC.kseg(sfr, m.D, *q)
            else:
                continue
            if cl is not None and m.add_clique(cl, tuple(q)): n += 1
    return n


def dump_dual(m, y, lam, val, path, thr=1e-9):
    """the LP dual at this round as a fractional packing: one line per row with y > thr,
    `cx cy theta h flag y`.  For a symmetric model the measure is y/8 on each D4 image of the
    pose (total mass sum(y) = LP value), its coverage is <= 1 at every column point."""
    with open(path + ".tmp", 'w') as f:
        f.write(f"# s={m.K}/{m.D} value={val:.9f} lambda={' '.join(f'{l:.7f}' for l in lam)} k={m.kvec} sym={int(m.sym)} mass={float(np.sum(y)):.9f} mass_in_R={float(sum(v for v, fl in zip(y, m.rflag) if fl)):.9f}\n")
        for (cx, cy, tm, h), fl, v in zip(m.rows, m.rflag, y):
            if v > thr: f.write(f"{cx:.9f} {cy:.9f} {tm:.9f} {h:.7f} {fl} {v:.9e}\n")
    os.replace(path + ".tmp", path)


def save_cols(m, path):
    """all integer points of the model (one representative per orbit; the D4 images are implied for
    symmetric models), so that a later run can start from the same column set (--cols)"""
    with open(path + ".tmp", 'w') as f:
        f.write(f"{m.K} {m.D} {int(m.sym)}\n")
        for o in m.orb_int: f.write(f"{int(o[0][0])} {int(o[0][1])}\n")
    os.replace(path + ".tmp", path)


def write_lp_dump(m, margin, path, val=None, x=None, lam=None, y=None):
    """the LP of this round, exactly as solve() builds it, as `path.npz` + `path_A.npz`:
    min c.z  s.t.  M z >= b,  lo <= z <= hi,  z = (x, lambda);  M = [A, -F] (csr), plus the row keys
    and flags, and the solution if given.  Read back with load_lp_dump()."""
    A = m.matrix(); n = m.ncols(); nr = len(m.rows); nl = m.nl
    M = sp.hstack([A, -m.flagmat()]).tocsr()
    kc = m.kcoef()
    c = np.concatenate([m.costv(), -kc])
    lo = np.concatenate([np.zeros(n), [LAM_FREE if kc[j] == 0 else LAM_LO for j in range(nl)]])
    hi = np.concatenate([np.full(n, np.inf), [LAM_FREE if kc[j] == 0 else LAM_HI for j in range(nl)]])
    b = np.full(nr, 1.0 + margin)
    sp.save_npz(path + "_A.npz", M, compressed=False)
    np.savez(path + ".npz", c=c, b=b, lo=lo, hi=hi, n=n, nl=nl, nr=nr, K=m.K, D=m.D, sym=int(m.sym), kvec=np.array(m.kvec), r=m.r,
             rows=np.array(m.rows), rflag=np.array(m.rflag), rowkeys=np.array([list(k) for k in m.rowkeys], dtype=np.int64),
             val=np.nan if val is None else val, x=np.zeros(0) if x is None else x, lam=np.zeros(0) if lam is None else lam, y=np.zeros(0) if y is None else y)


def load_lp_dump(path):
    d = dict(np.load(path + ".npz", allow_pickle=False)); d['M'] = sp.load_npz(path + "_A.npz").tocsr()
    return d


def load_cols(m, path, raw=False):
    """add the points of a column file (written for the same K, D) as columns; returns number added.
    Asymmetric models get all 8 images of each point unless raw=True (then the file's points only:
    the checkpoint of an asymmetric run already lists every column, and symmetrising it makes the LP
    up to 8x wider than the run it came from)."""
    n = 0
    with open(path) as f:
        K, D, sym = (int(v) for v in f.readline().split())
        assert K == m.K and D == m.D, f"column file {path} is for container {K}/{D}, model is {m.K}/{m.D}"
        for line in f:
            X, Y = (int(v) for v in line.split())
            if m.sym or raw: n += m.add_orbit(X, Y)
            else:
                for p in set([(X, Y), (K - X, Y), (X, K - Y), (K - X, K - Y), (Y, X), (K - Y, X), (Y, K - X), (K - Y, K - X)]):
                    n += m.add_orbit(*p)
    return n


def float_rows(m, x, lam, N, topk, probe_margin, pitch, dtheta, pool, nproc, seeds=None):
    """the float separation oracle of search/floatsep.py on the current LP point, in exactly the
    units the probe file would be written in.  Returns (min value, rows, info).  Only which rows
    are added depends on this; validity of a row does not (see floatsep's "Rows are valid")."""
    n_orb = len(m.orb_int); f = 1.0 / (1.0 + probe_margin)
    w = x[:n_orb][m.own] * f
    cq = [(float(x[n_orb + j]) * f, imgs) for j, imgs in enumerate(m.cliques) if x[n_orb + j] > 1e-12]
    lamv = [float(l) * f for l in (lam if m.perbox else [lam[0]] * 4)]
    ks = FS.bin_list(N, m.sym and not m.cliques, math.radians(dtheta))
    return FS.separate(m.P, w, m.s, N, ks, topk=topk, pitch=pitch, cliques=cq, lam=lamv, r=m.r,
                       pool=pool, nproc=nproc, seeds=seeds)


def float_sep_two_tier(m, x, lam, N, topk, probe_margin, pitch, dtheta, dtheta_full, pool, nproc, log=None, seeds=None):
    """The float oracle at the coarse angle pitch, CONFIRMED at the full one before it is allowed
    to report `probe_min >= 1` (task K follow-up; see search/LOOPSPEED.md 6).

    A coarse scan looks at a tenth of the angle net, and on a converged leaf the worst poses live
    in the bins it skips: measured on `runs/branch_L1110f_probe.txt`, whose exact minimum is
    0.636114, the scan reports 0.950 at 0.4 deg, 0.688 at 0.2 deg and 0.636 at every bin.  Taking
    the coarse number as "clean" is what made the loop spend an exact round to be told the deficit
    was 0.64, not 0.99.  So: separate coarsely every round (that is only about choosing rows), and
    the moment the coarse scan says it is clean, re-scan every bin and believe THAT."""
    mv, rows, info = float_rows(m, x, lam, N, topk, probe_margin, pitch, dtheta, pool, nproc, seeds)
    info['coarse_min'] = mv; info['full'] = False
    if mv >= 1.0 and dtheta_full != dtheta:
        mv2, rows2, info2 = float_rows(m, x, lam, N, topk, probe_margin, pitch, dtheta_full, pool, nproc, seeds)
        info = dict(info2); info['coarse_min'] = mv; info['full'] = True
        if log: log(f"    float: coarse ({info2['bins']} of the net) said {mv:.7f}, full net says {mv2:.7f}")
        return mv2, rows2, info
    return mv, rows, info


def loop(m, tag, margin=2e-6, N=6000, topk=None, max_iters=400, log=print, probe_margin=None, x0=None,
         prune_at=None, colgen=0, cg_want=200, cg_pitch=0.01, n=12, threads=None, dump_lp=None,
         cliques=0, cq_want=40, cq_pitch=0.01, matched=False, cq_interior=False,
         sep_mode='exact', exact_every=5, sep_pitch=0.004, sep_dtheta=0.4, sep_dtheta_full=0.0, stab=0.0):
    """cutting-plane loop; returns (x, lam, obj, info).  Asymmetric models sweep [0,90deg), i.e.
    2.4x the bins: fewer witnesses per bin and a higher prune threshold keep the row set stable
    (pruning every round makes the LP vertex jump and the dropped rows come straight back).

    `sep_mode='float'` separates with search/floatsep.py instead and runs the exact verifier only
    every `exact_every` rounds and at convergence (the float oracle clean twice running).  The
    loop still only ever STOPS on an exact verdict, and nothing is claimed until the final file
    passes `verify` + `xcheck.py`."""
    if probe_margin is None: probe_margin = margin / 2
    if topk is None: topk = 6 if m.sym else 3
    if prune_at is None: prune_at = 40000 if m.sym else 60000       # solve time is superlinear in rows: ~20 s at 30k, hours at 120k
    tmp = f"runs/branch_{tag}_probe.txt"; sep = f"runs/branch_{tag}_sep.txt"
    t0 = time.time(); x = x0; lam = None; val = None; it = 0
    pool = None; fclean = 0; xstab = None; lamstab = None; seeds = []
    if sep_mode == 'float':
        import multiprocessing as mp
        pool = mp.get_context('fork').Pool(threads or 8)
    while it < max_iters:
        t1 = time.time(); out = m.solve(margin); t_lp = time.time() - t1
        if out is None:
            log(f"  it{it} LP infeasible/failed (rows={len(m.rows)})")
            if pool is not None: pool.close(); pool.join()
            return None, None, None, dict(status='infeasible')
        val, x, lam, y = out; tw = float(m.costv() @ x)
        # the matched pure value: the SAME master with the clique columns removed.  It has to run
        # here, before this round's column generation: a column added in between is marked active
        # by the next solve, and the "pure" LP then has MORE columns than the clique one and comes
        # out lower (measured: 328 columns against 168, and a gain of -0.02).
        pureval = None
        if matched and m.cliques and SOLVER == 'restricted':
            o2 = m.solve_restricted(margin, no_cliques=True)
            if o2 is not None: pureval = o2[0]
        dump_dual(m, y, lam, val, f"runs/branch_{tag}_dual.txt")
        if dump_lp: write_lp_dump(m, margin, f"{dump_lp}_it{it}", val=val, x=x, lam=lam, y=y)
        # ---- separation ------------------------------------------------------------------
        # `sep_mode='exact'` (default) is the original: export the probe and run the verifier.
        # `sep_mode='float'` runs the numpy oracle every round and the verifier only every
        # `exact_every` rounds or once the float oracle has been clean twice running; the loop
        # only ever STOPS on an exact verdict.
        mv = None; mvf = None; t_ver = 0.0; t_fsep = 0.0; nfviol = 0; wit = []
        want_exact = (sep_mode != 'float') or (it % exact_every == 0)
        if sep_mode == 'float':
            # in-out stabilisation (--stab): separate at a convex combination of this round's LP
            # point and the previous separation point.  Rows are valid at any point, and the
            # exact rounds always separate at the true LP point, so the stopping rule is unchanged.
            if stab > 0 and xstab is not None and not want_exact:
                nx0 = min(len(x), len(xstab)); xs = x.copy(); xs[:nx0] = (1 - stab) * x[:nx0] + stab * xstab[:nx0]
                lams = (1 - stab) * np.asarray(lam) + stab * np.asarray(lamstab)
            else:
                xs = x; lams = lam
            xstab = xs.copy(); lamstab = np.asarray(lams).copy()
            t1 = time.time()
            mvf, fwit, finfo = float_sep_two_tier(m, xs, lams, N, topk, probe_margin, sep_pitch,
                                                  sep_dtheta, sep_dtheta_full, pool, threads or 8, log=log,
                                                  seeds=seeds)
            t_fsep = time.time() - t1; nfviol = len(fwit)
            fclean = fclean + 1 if mvf >= 1.0 else 0
            if fclean >= 2: want_exact = True
            wit = fwit
        if want_exact:
            export(m, x, lam, tmp, WD=10 ** 12, factor=1.0 / (1.0 + probe_margin), up=False)
            t1 = time.time(); mv, ok = T.run_verifier(tmp, N, topk=topk, sep=sep, n=n, threads=threads or T.NPROC); t_ver = time.time() - t1
            wit = read_witnesses(sep, N) if os.path.exists(sep) else []
            if os.path.exists(sep): os.remove(sep)
            if mv is not None and mv < 1: fclean = 0
        elif sep_mode == 'float':
            # keep the probe checkpoint fresh even in float rounds (a restart reads it)
            export(m, x, lam, tmp, WD=10 ** 12, factor=1.0 / (1.0 + probe_margin), up=False)
        n_orb = len(m.orb_int); xo = x[:n_orb]
        nz = int((xo > 1e-12).sum()); natoms = int(m.sizes[xo > 1e-12].sum())
        ncol = 0; Mcov = None; t_cg = 0.0; cand = []
        if it < colgen:
            t1 = time.time()
            cand = T.price(m, y, pitch=cg_pitch, want=cg_want) if m.sym else price_asym(m, y, pitch=max(cg_pitch, 0.02), want=cg_want)
            Mcov = cand[0][0] if cand else 1.0
            for (cv, X, Y) in cand: ncol += m.add_orbit(X, Y)
            t_cg = time.time() - t1
        wcl = float(m.csizes[:len(x) - n_orb] @ x[n_orb:]) if len(x) > n_orb else 0.0
        ncl_used = int((x[n_orb:] > 1e-12).sum()) if len(x) > n_orb else 0
        ncq = 0; cqmass = None; cqpt = None; t_cq = 0.0
        if it < cliques:
            t1 = time.time(); ncq, cqmass, cqpt = price_cliques(m, y, x, want=cq_want, pitch=cq_pitch, interior=cq_interior); t_cq = time.time() - t1
        # the poses that were worst this round seed next round's local descent: they barely move,
        # and a grid alone will not land on them again (search/LOOPSPEED.md 6)
        if sep_mode == 'float': seeds = list(wit[:4000])
        t1 = time.time(); added = m.add_rows(wit); t_rows = time.time() - t1
        save_cols(m, f"runs/branch_{tag}_cols.txt")          # checkpoint: every column, weight or not, so a restart loses nothing
        if m.cliques: save_cliques(m, f"runs/branch_{tag}_cliques.txt")
        nR = sum(1 for f in m.rflag if f)
        lamtxt = ",".join(f"{l:+.4f}" for l in lam)
        ri = getattr(m, '_restricted_info', None) if SOLVER == 'restricted' else None
        pm = f"{float(mv):.7f}" if mv is not None else f"{mvf:.7f}(float)"
        log(f"  it{it} total={tw:.6f} lam=[{lamtxt}] obj={val:.6f} probe_min={pm} viol={len(wit)} new_rows={added} rows={len(m.rows)} (inR {nR}) support={nz} cols/{natoms} atoms"
            + (f" dualcov={Mcov:.4f} new_cols={ncol} cols={len(m.orbits)}" if Mcov is not None else "")
            + (f" pure={pureval:.6f} gain={pureval - val:+.6f}" if pureval is not None else "")
            + (f" cliques={len(m.cliques)} (+{ncq}, used {ncl_used}, weight {wcl:.6f}" + (f", ybar {cqmass:.4f} vs point {cqpt:.4f})" if cqmass is not None else ")") if (m.cliques or ncq) else "")
            + f" t={time.time()-t0:.0f}s [lp {t_lp:.0f}s ver {t_ver:.0f}s fsep {t_fsep:.0f}s cg {t_cg:.0f}s cq {t_cq:.0f}s rows {t_rows:.0f}s]"
            + (f" [sep float_min={mvf:.7f} fviol={nfviol} clean={fclean}{' full' if finfo.get('full') else ''}{' exact' if mv is not None else ''}]" if sep_mode == 'float' else "")
            + (f" [restricted: {ri['cols']} active cols, {ri['passes']} passes, {ri['neg_rc']} neg rc left]" if ri else ""))
        # only an EXACT verdict may stop the loop
        if added == 0 and mv is not None and mv >= 1 and ncol == 0 and ncq == 0: break
        if len(m.rows) > prune_at:
            # keep rows with positive dual, rows that are nearly tight (zero dual but no slack --
            # pruning those re-violates them next round and the loop cycles), the rows added this
            # iteration (y does not cover them; x does not cover this iteration's new columns), and
            # the most recent 30000
            xp = np.zeros(m.ncols()); xp[:len(x)] = x
            A = m.matrix(); ny = len(y)
            slack = A[:ny] @ xp - (m.flagmat()[:ny] @ lam) - (1.0 + margin)
            keep = np.ones(len(m.rows), dtype=bool); keep[:ny] = (y > 1e-12) | (slack < 2e-2); keep[-15000:] = True
            m.prune(keep); log(f"   pruned to {len(m.rows)} rows")
        it += 1
    if pool is not None: pool.close(); pool.join()
    return x, lam, val, dict(iters=it, rows=len(m.rows), t=time.time() - t0)


def finalize(m, x, lam, tag, out_path, N=6000, log=print, WD=10 ** 7, margin=2e-6, n=12, threads=None):
    for rnd in range(20):
        tw, npts, L = export(m, x, lam, out_path, WD=WD, up=True)
        mv1, ok1 = T.run_verifier(out_path, N, n=n, threads=threads or T.NPROC)
        sep = f"runs/branch_{tag}_sep2.txt"
        mv2, ok2 = T.run_verifier(out_path, 2 * N, topk=6, sep=sep, n=n, threads=threads or T.NPROC)
        wit = read_witnesses(sep, 2 * N) if os.path.exists(sep) else []
        if os.path.exists(sep): os.remove(sep)
        eff = tw - float(m.kcoef() @ L)
        log(f"  final[{rnd}] {out_path}: points={npts} total={tw:.7f} lambda=[{','.join(f'{l:+.7f}' for l in L)}] k={m.kvec if m.perbox else m.kvec[0]} total-lambda.k={eff:.7f} min@{N}={float(mv1):.7f} {'OK' if ok1 else 'FAIL'}  min@{2*N}={float(mv2):.7f} {'OK' if ok2 else 'FAIL'} (viol {len(wit)})")
        if ok1 and ok2: return tw, npts, L, mv1, mv2
        if eff >= n and mv1 >= 1 and mv2 >= 1:
            log(f"  final: covering holds at both nets but total - lambda.k = {eff:.6f} >= {n}: not a certificate")
            return tw, npts, L, mv1, mv2
        m.add_rows(wit)
        x, lam, val, info = loop(m, tag, margin=margin, N=N, log=log, x0=x, n=n, threads=threads)
        if x is None: return None
    return None


def run(cert, tag, kreg, r, Dp=None, mul=1, N=6000, topk=6, margin=2e-6, log=print, warm=True, out=None,
        colgen=0, cg_want=200, cg_pitch=0.01, n=12, threads=None, cols=None, warm_thr=1.05, prune_at=None, row_grid=None, max_iters=None, dump_lp=None, cols_raw=False,
        cliques=0, cq_want=40, cq_load=None, cq_pitch=0.01, matched=False, cq_interior=False,
        sep_mode='exact', exact_every=5, sep_pitch=0.004, sep_dtheta=0.4, sep_dtheta_full=0.0, stab=0.0, seed_from=None, corner_rows_on=False):
    m, w0, s = build_model(cert, r, kreg, Dp, mul, raw_points=cols_raw)
    if row_grid: m.row_grid = tuple(row_grid)
    log(f"[{tag}] {cert}: {len(m.orbits)} columns / {len(m.P)} atoms ({'D4 orbits' if m.sym else 'single points'}), container {s} = {float(s):.7f}, r = {m.rfrac} = {m.r:.4f}, k = {kreg}, input total {float(m.sizes @ w0):.6f}")
    t0 = time.time()
    if warm:
        lr = lattice_rows(m, w0[m.own], thr=warm_thr)
        m.add_rows(lr); log(f"[{tag}] warm start: {len(m.rows)} lattice rows (capture < {warm_thr} under the input weights), {sum(1 for f in m.rflag if f)} in boxes ({time.time()-t0:.0f}s)")
    if corner_rows_on:
        nb = m.add_rows(corner_rows(m))
        log(f"[{tag}] corner rows: +{nb} poses centred in the active corner boxes -> {len(m.rows)} rows, {sum(1 for f in m.rflag if f)} in boxes")
    for sf in (seed_from or []):
        base = sf if os.path.sep in sf else f"runs/branch_{sf}"
        for suf, fn, what in (("_cols.txt", lambda p: load_cols(m, p, raw=cols_raw), "columns"),
                              ("_cliques.txt", lambda p: load_cliques(m, p), "clique columns"),
                              ("_dual.txt", lambda p: m.add_rows(seed_rows(p)), "rows")):
            p = base + suf
            if os.path.exists(p):
                k = fn(p); log(f"[{tag}] seed-from {p}: +{k} {what}")
            else:
                log(f"[{tag}] seed-from: {p} does not exist, skipped")
    for cf in (cols or []):
        nc = load_cols(m, cf, raw=cols_raw); log(f"[{tag}] columns from {cf}{' (raw)' if cols_raw else ''}: +{nc} -> {len(m.orbits)} columns")
    for cf in (cq_load or []):
        nc = load_cliques(m, cf); log(f"[{tag}] anchor cliques from {cf}: +{nc} -> {len(m.cliques)} clique columns")
    if SOLVER == 'restricted':
        # the restricted master starts from the input certificate's support; checkpoint columns (mostly
        # dead) are in the model but priced in only when their reduced cost says so
        m._active = np.concatenate([w0 > 0, np.zeros(len(m.orbits) - len(w0), dtype=bool)])
        log(f"[{tag}] restricted master: {int(m._active.sum())} of {len(m.orbits)} columns active at start")
    x, lam, val, info = loop(m, tag, margin=margin, N=N, topk=topk, log=log, colgen=colgen, cg_want=cg_want, cg_pitch=cg_pitch, n=n, threads=threads, prune_at=prune_at,
                             max_iters=max_iters if max_iters is not None else 400, dump_lp=dump_lp, cliques=cliques, cq_want=cq_want, cq_pitch=cq_pitch, matched=matched, cq_interior=cq_interior,
                             sep_mode=sep_mode, exact_every=exact_every, sep_pitch=sep_pitch, sep_dtheta=sep_dtheta,
                             sep_dtheta_full=sep_dtheta_full, stab=stab)
    if x is None: log(f"[{tag}] FAILED"); return None
    if max_iters is not None:
        log(f"[{tag}] stopped after {info.get('iters')} rounds (--max-iters): LP value {val:.7f}, dual in runs/branch_{tag}_dual.txt; no finalize")
        return None
    out = out or f"runs/branch_{tag}.txt"
    res = finalize(m, x, lam, tag, out, N=N, log=log, margin=margin, n=n, threads=threads)
    if res is None: log(f"[{tag}] finalize FAILED"); return None
    tw, npts, L, mv1, mv2 = res
    eff = tw - float(m.kcoef() @ L)
    log(f"[{tag}] RESULT s={s} ({float(s):.7f}) k={kreg} r={m.rfrac} points={npts} total={tw:.7f} lambda=[{','.join(f'{l:+.7f}' for l in L)}] total-lambda.k={eff:.7f} {'< ' + str(n) + ': LEAF CLOSED' if eff < n else '>= ' + str(n) + ': leaf open'} min@{N}={mv1} min@{2*N}={mv2} t={time.time()-t0:.0f}s")
    return dict(cert=out, s=str(s), sf=float(s), k=str(kreg), r=str(m.rfrac), points=npts, total=tw, lam=[float(l) for l in L], effective=eff, closed=bool(eff < n),
                min1=str(mv1), min2=str(mv2), t=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cert'); ap.add_argument('tag')
    ap.add_argument('--k', required=True, help='single: number of squares centred in the corner boxes (0..4); per box: a pattern like 1100')
    ap.add_argument('--r', default='1', help='corner box side (fraction), needs (2r-1)^2 < 2')
    ap.add_argument('--Dp', type=int, default=None); ap.add_argument('--mul', type=int, default=1)
    ap.add_argument('--N', type=int, default=6000); ap.add_argument('--topk', type=int, default=None, help='witnesses per angle bin (default 6, asymmetric 3)')
    ap.add_argument('--margin', type=float, default=2e-6)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--out', default=None)
    ap.add_argument('--colgen', type=int, default=0); ap.add_argument('--cg-want', type=int, default=200); ap.add_argument('--cg-pitch', type=float, default=0.01)
    ap.add_argument('--n', type=int, default=12); ap.add_argument('--threads', type=int, default=None)
    ap.add_argument('--no-warm', action='store_true')
    ap.add_argument('--cols', action='append', default=None, help='column checkpoint file(s) (runs/branch_TAG_cols.txt) whose points are added as columns')
    ap.add_argument('--warm-thr', type=float, default=1.05, help='warm start with lattice poses capturing < this under the input weights (use ~3 when restarting from a converged probe)')
    ap.add_argument('--prune-at', type=int, default=None, help='prune the row set when it exceeds this (default 40000 symmetric / 60000 asymmetric; asymmetric leaves converge better unpruned)')
    ap.add_argument('--lam-hi', type=float, default=None, help='upper bound on the multipliers (default 1.5; the k = 4 leaf sits at the cap, so a higher cap can only lower its value)')
    ap.add_argument('--lam-lo', type=float, default=None, help='lower bound on the multipliers (default -1); --lam-lo 0 --lam-hi 0 pins every multiplier at 0, i.e. runs the PURE cover LP with a vacuous trailer')
    ap.add_argument('--max-iters', type=int, default=None, help='stop the cutting-plane loop after this many rounds and skip finalize (diagnostic runs: the dual of every round is in runs/branch_TAG_dual.txt)')
    ap.add_argument('--row-grid', type=float, nargs=2, default=None, metavar=('POS', 'ANG'), help='merge rows whose centre/angle agree to this resolution (default 1e-7 1e-8; e.g. 0.002 0.005 keeps the LP small)')
    ap.add_argument('--cols-raw', action='store_true', help='asymmetric leaves: take the certificate\'s and --cols\' points as they are (no D4 images); use when restarting a leaf from its own probe + checkpoint, which already list every column')
    ap.add_argument('--cliques', type=int, default=0, metavar='ROUNDS', help='price anchor-clique columns K(p,A) (search/anchorclique.py, certificates/FORMAT.md) for this many rounds')
    ap.add_argument('--cq-want', type=int, default=40, help='anchor-clique columns added per round (default 40)')
    ap.add_argument('--matched', action='store_true', help='every round, also solve the SAME LP with the clique columns removed and log the pure value and the gain (search/CLIQUE_CONTINUUM.md calls this a matched pair)')
    ap.add_argument('--cq-pitch', type=float, default=0.01, help='grid pitch of the anchor-clique separator (default 0.01)')
    ap.add_argument('--cq-interior', action='store_true', help='also separate anchor cliques at INTERIOR anchor points and arbitrary anchor directions (search/anchorsep.py separate_interior; search/RECONCILE.md -- the largest violations are at wall distance > 1)')
    ap.add_argument('--cq-load', action='append', default=None, help='anchor-clique checkpoint file(s) (runs/branch_TAG_cliques.txt) to start from')
    ap.add_argument('--dump-lp', default=None, metavar='PREFIX', help='write the LP of every round of the cutting-plane loop as PREFIX_itN.npz + PREFIX_itN_A.npz (see write_lp_dump / load_lp_dump; used by search/lp_bench.py)')
    # ---- task K (search/LOOPSPEED.md): the separation oracle, and fewer rounds.  Every default
    # here reproduces the old behaviour exactly; nothing changes unless a flag is given.
    ap.add_argument('--sep', dest='sep_mode', choices=['exact', 'float'], default='exact',
                    help="separation oracle: 'exact' (default) runs the Rust verifier every round, as before; 'float' runs the numpy oracle of search/floatsep.py every round and the verifier only every --exact-every rounds and at convergence.  The loop still only STOPS on an exact verdict, and the final file is verified unchanged")
    ap.add_argument('--exact-every', type=int, default=5, help='--sep float: run the exact verifier every this many rounds (default 5; round 0 is always exact)')
    ap.add_argument('--sep-pitch', type=float, default=0.004, help='--sep float: centre pitch of the dense scan (default 0.004)')
    ap.add_argument('--sep-dtheta', type=float, default=0.4, help='--sep float: angle pitch of the dense scan, in degrees (default 0.4).  Only which rows are added depends on it')
    ap.add_argument('--sep-dtheta-full', type=float, default=0.0, help="--sep float: angle pitch of the CONFIRMING scan, run whenever the coarse one reports probe_min >= 1 (default 0.0 = every bin of the net).  Set equal to --sep-dtheta to switch the confirmation off")
    ap.add_argument('--stab', type=float, default=0.0, help='in-out stabilisation: separate at (1-stab)*x + stab*(previous separation point) in the float rounds, to damp the probe sawtooth (0 = off, the default; 0.3 is a reasonable value)')
    ap.add_argument('--seed-from', action='append', default=None, metavar='LEAF', help='start from a sibling leaf\'s checkpoint: its columns, clique columns and the rows of its dual (runs/branch_LEAF_{cols,cliques,dual}.txt, or a path prefix).  The per-slot leaves differ only in the multipliers, so a converged sibling\'s binding rows are almost the right row set')
    ap.add_argument('--corner-rows', action='store_true', help='seed the warm start with poses centred in the ACTIVE corner boxes (lattice, h = 1/2).  Without them round 0 has no row charged to any multiplier and lambda goes straight to its cap')
    a = ap.parse_args()
    global LAM_HI, LAM_LO
    if a.lam_hi is not None: LAM_HI = float(a.lam_hi)
    if a.lam_lo is not None: LAM_LO = float(a.lam_lo)
    r = Fraction(a.r); assert (2 * r - 1) ** 2 < 2, "r too large: (2r-1)^2 < 2 needed"
    kreg = parse_k(a.k)
    T.HIGHS['random_seed'] = a.seed
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/branch_{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"branch.py {' '.join(sys.argv[1:])}")
    res = run(a.cert, a.tag, kreg, r, Dp=a.Dp, mul=a.mul, N=a.N, topk=a.topk, margin=a.margin, log=log, warm=not a.no_warm, out=a.out,
              colgen=a.colgen, cg_want=a.cg_want, cg_pitch=a.cg_pitch, n=a.n, threads=a.threads, cols=a.cols, warm_thr=a.warm_thr, prune_at=a.prune_at, row_grid=a.row_grid, max_iters=a.max_iters, dump_lp=a.dump_lp, cols_raw=a.cols_raw,
              cliques=a.cliques, cq_want=a.cq_want, cq_load=a.cq_load, cq_pitch=a.cq_pitch, matched=a.matched, cq_interior=a.cq_interior,
              sep_mode=a.sep_mode, exact_every=a.exact_every, sep_pitch=a.sep_pitch, sep_dtheta=a.sep_dtheta,
              sep_dtheta_full=a.sep_dtheta_full, stab=a.stab,
              seed_from=a.seed_from, corner_rows_on=a.corner_rows)
    if res: json.dump(res, open(f"runs/branch_{a.tag}.json", 'w'), indent=1)


if __name__ == '__main__':
    main()
