#!/usr/bin/env python3
"""lp_bench.py -- benchmark LP backends on a dumped branch.py / tighten.py LP (tasks/lp-speed).

    min c.z   s.t.  M z >= b,   lo <= z <= hi        (z = (x, lambda);  M = [A, -F])

Dumps come from `branch.py --dump-lp PREFIX` (PREFIX_itN.npz + PREFIX_itN_A.npz, one per round).

    python3 search/lp_bench.py runs/lp1110_it2 --methods ipm,ds,hs-ds,hs-ipm,pdlp,sift --reps 2
    python3 search/lp_bench.py runs/lp1110_it2 --warm-from runs/lp1110_it1 --methods warm,inc,sift-warm

Methods
  ipm      scipy linprog method='highs-ipm' (what branch.py runs today)       -- baseline
  ds       scipy linprog method='highs-ds', cold
  hs-ds    highspy dual simplex, cold (model passed once)
  hs-ipm   highspy IPM, cold
  pdlp     highspy solver='pdlp' (first-order; kkt_tolerance --pdlp-tol)
  warm     highspy dual simplex on THIS dump, starting from the optimal basis of --warm-from
           (rows matched by key, new rows basic; columns by index, new columns nonbasic at lower)
  inc      highspy: load --warm-from, solve, then addCols/addRows/changeColsBounds to reach THIS dump
           on the same Highs object (hot start: internal basis + factorization kept), solve again;
           the second solve is what is timed
  inc2     two-stage warm start on one Highs object: zero-change re-run (basis cleanliness), then
           new rows only + DUAL simplex (basis stays dual feasible), then new columns + PRIMAL simplex
           (basis stays primal feasible) -- no phase 1 in either stage
  csift    column sifting / restricted master: IPM over the columns in --warm-from's support (+ the
           last --csift-extra columns), price all columns by reduced cost, add --csift-add per pass;
           csift-ds = the same with dual simplex sub-solves
  sift     row sifting: solve on a subset of rows, scan all rows for violation (sparse mat-vec), add
           the most violated (--sift-add per pass), repeat until none violated.  Inner solves are
           highspy dual simplex on one Highs object (rows appended, basis kept).  Starts from the rows
           of --warm-from that were tight (y > 0 or slack < --sift-slack) if given, else the first
           --sift-init rows.
  sift-warm  sift with the starting set = tight rows of --warm-from (same as sift when --warm-from given)

Every solve reports objective, seconds, iterations; objectives are compared to the ipm value of the
same dump (or the recorded value in the dump).  Pin with `taskset` and `--threads`.
"""
import sys, os, time, argparse, json, warnings
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import highspy

INF = highspy.kHighsInf


def load(path):
    d = dict(np.load(path + ".npz", allow_pickle=False)); d['M'] = sp.load_npz(path + "_A.npz").tocsr()
    d['hi'] = np.where(np.isinf(d['hi']), INF, d['hi'])
    d['keys'] = [tuple(int(v) for v in k) for k in d['rowkeys']]
    return d


def describe(d):
    M = d['M']; return f"rows={M.shape[0]} cols={M.shape[1]} nnz={M.nnz} (n={int(d['n'])} nl={int(d['nl'])}) recorded_val={float(d['val']):.9f}"


# ---------------------------------------------------------------- scipy backends
def run_scipy(d, method, threads):
    M, c, b, lo, hi = d['M'], d['c'], d['b'], d['lo'], d['hi']
    bounds = [(l, None if h >= INF else h) for l, h in zip(lo, hi)]
    t = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = linprog(c=c, A_ub=-M, b_ub=-b, bounds=bounds, method=method, options=dict(random_seed=0))
    return dict(obj=float(res.fun) if res.success else np.nan, t=time.time() - t, it=int(getattr(res, 'nit', -1)), status=int(res.status),
                z=res.x if res.success else None, y=(-res.ineqlin.marginals if res.success else None))


# ---------------------------------------------------------------- highspy helpers
def hs_model(d, rows=None):
    """HighsLp for the dump, optionally restricted to a row subset (index array)"""
    M = d['M'] if rows is None else d['M'][rows]
    Mc = M.tocsc()
    lp = highspy.HighsLp()
    lp.num_col_ = M.shape[1]; lp.num_row_ = M.shape[0]
    lp.col_cost_ = np.asarray(d['c'], dtype=float); lp.col_lower_ = np.asarray(d['lo'], dtype=float); lp.col_upper_ = np.asarray(d['hi'], dtype=float)
    bb = d['b'] if rows is None else d['b'][rows]
    lp.row_lower_ = np.asarray(bb, dtype=float); lp.row_upper_ = np.full(M.shape[0], INF)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = Mc.indptr.astype(np.int32); lp.a_matrix_.index_ = Mc.indices.astype(np.int32); lp.a_matrix_.value_ = Mc.data.astype(float)
    return lp


def hs_new(threads, solver, log=False, pdlp_tol=1e-7, seed=0):
    h = highspy.Highs(); h.setOptionValue('output_flag', bool(log)); h.setOptionValue('threads', int(threads)); h.setOptionValue('random_seed', int(seed))
    if solver == 'simplex':
        h.setOptionValue('solver', 'simplex'); h.setOptionValue('simplex_strategy', 1)     # dual
    elif solver == 'ipm':
        h.setOptionValue('solver', 'ipm')
    elif solver == 'pdlp':
        h.setOptionValue('solver', 'pdlp'); h.setOptionValue('kkt_tolerance', float(pdlp_tol))
        h.setOptionValue('run_crossover', 'off')
    return h


def hs_result(h, t0):
    st = h.getModelStatus(); info = h.getInfo(); sol = h.getSolution()
    ok = st == highspy.HighsModelStatus.kOptimal
    it = max(int(info.simplex_iteration_count), int(info.ipm_iteration_count), int(getattr(info, 'pdlp_iteration_count', 0)))
    return dict(obj=float(info.objective_function_value) if ok else np.nan, t=time.time() - t0, it=it, status=str(st).split('.')[-1],
                z=np.array(sol.col_value) if ok else None, y=np.array(sol.row_dual) if ok else None)


def run_hs(d, solver, threads, pdlp_tol, log=False):
    t0 = time.time()
    h = hs_new(threads, solver, log=log, pdlp_tol=pdlp_tol); h.passModel(hs_model(d)); t_build = time.time() - t0
    h.run(); r = hs_result(h, t0); r['t_build'] = t_build
    return r


def basis_from(hprev, dprev, d):
    """basis for dump d from the optimal basis of hprev (solved on dprev): rows by key, cols by index"""
    bp = hprev.getBasis()
    n = d['M'].shape[1]; nr = d['M'].shape[0]
    b = highspy.HighsBasis(); b.valid = True
    cs = [highspy.HighsBasisStatus.kLower] * n
    pc = list(bp.col_status); np_ = min(len(pc), n)
    cs[:np_] = pc[:np_]
    prow = dict(zip(dprev['keys'], bp.row_status))
    rs = [prow.get(k, highspy.HighsBasisStatus.kBasic) for k in d['keys']]
    b.col_status = cs; b.row_status = rs
    nb = sum(1 for s in cs if s == highspy.HighsBasisStatus.kBasic) + sum(1 for s in rs if s == highspy.HighsBasisStatus.kBasic)
    return b, nb, nr


PREV_SOLVER = 'simplex'      # how the --warm-from dump is solved to obtain the starting basis ('simplex' or 'ipm' = IPM + crossover)


HLOG = False     # HiGHS console output for the warm/inc solves (--log)


def run_warm(d, dprev, threads):
    """solve dprev cold (not timed), then d with the carried basis (timed)"""
    h0 = hs_new(threads, PREV_SOLVER, log=HLOG); h0.passModel(hs_model(dprev)); t = time.time(); h0.run(); t_prev = time.time() - t
    st0 = h0.getModelStatus()
    t0 = time.time()
    h = hs_new(threads, 'simplex', log=HLOG); h.passModel(hs_model(d))
    b, nb, nr = basis_from(h0, dprev, d)
    h.setBasis(b)
    t_build = time.time() - t0
    h.run(); r = hs_result(h, t0); r['t_build'] = t_build; r['t_prev_cold'] = t_prev; r['basic_in_basis'] = nb; r['nr'] = nr; r['prev_status'] = str(st0).split('.')[-1]
    return r


def run_inc(d, dprev, threads):
    """solve dprev cold, then modify the same Highs object into d (append cols, append rows, bounds), re-solve (timed)"""
    h = hs_new(threads, PREV_SOLVER, log=HLOG); h.passModel(hs_model(dprev)); t = time.time(); h.run(); t_prev = time.time() - t
    prev_obj = float(h.getInfo().objective_function_value)
    print(f"    [inc] prev solved by {PREV_SOLVER}: obj={prev_obj:.9f} t={t_prev:.0f}s status={str(h.getModelStatus()).split('.')[-1]}", flush=True)
    if PREV_SOLVER != 'simplex':          # continue with dual simplex from the crossover basis
        h.setOptionValue('solver', 'simplex'); h.setOptionValue('simplex_strategy', 1)
    t0 = time.time()
    n0, nr0 = dprev['M'].shape[1], dprev['M'].shape[0]; n1, nr1 = d['M'].shape[1], d['M'].shape[0]
    # column index map: prev columns 0..n0-nl-1 are the same points; lambda columns are the last nl in both
    nl = int(d['nl']); assert nl == int(dprev['nl'])
    # rebuild as: [old point cols][new point cols][lambda] -- so permute d's columns accordingly
    perm = np.concatenate([np.arange(n0 - nl), np.arange(n0 - nl, n1 - nl), np.arange(n1 - nl, n1)])   # identity: d already has old points first
    M = d['M']
    # check the old-row/old-col block is unchanged (row keys of prev must be a prefix-set of d's)
    kprev = dprev['keys']; kd = d['keys']; pos = {k: i for i, k in enumerate(kd)}
    old_rows = np.array([pos[k] for k in kprev if k in pos]); assert len(old_rows) == len(kprev), "prev rows missing in d (pruned?) -- inc needs monotone rows"
    new_rows = np.array(sorted(set(range(nr1)) - set(old_rows.tolist())), dtype=np.int64)
    # new columns (point columns n0-nl .. n1-nl-1): append with entries against the old rows only (in prev row order)
    newc = np.arange(n0 - nl, n1 - nl)
    if len(newc):
        Mc = M[old_rows][:, newc].tocsc()
        h.addCols(len(newc), np.asarray(d['c'][newc], dtype=float), np.asarray(d['lo'][newc], dtype=float), np.asarray(d['hi'][newc], dtype=float),
                  Mc.nnz, Mc.indptr.astype(np.int32), Mc.indices.astype(np.int32), Mc.data.astype(float))
        # lambda columns now sit before the new columns in h; move nothing -- instead permute: h's cols = [old points][lambda][new points]
    # rows: new rows, with entries mapped to h's column order [old points][lambda][new points]
    if len(new_rows):
        Mr = M[new_rows].tocsr()
        colmap = np.empty(n1, dtype=np.int64); colmap[:n0 - nl] = np.arange(n0 - nl); colmap[n1 - nl:] = np.arange(n0 - nl, n0); colmap[n0 - nl:n1 - nl] = np.arange(n0, n1)
        idx = colmap[Mr.indices]
        h.addRows(len(new_rows), np.asarray(d['b'][new_rows], dtype=float), np.full(len(new_rows), INF), Mr.nnz, Mr.indptr.astype(np.int32), idx.astype(np.int32), Mr.data.astype(float))
    # lambda bounds may differ (they don't in branch.py, but be safe)
    lam_cols = np.arange(n0 - nl, n0)
    h.changeColsBounds(nl, lam_cols.astype(np.int32), np.asarray(d['lo'][n1 - nl:], dtype=float), np.asarray(d['hi'][n1 - nl:], dtype=float))
    t_build = time.time() - t0
    h.run(); r = hs_result(h, t0); r['t_build'] = t_build; r['t_prev_cold'] = t_prev; r['prev_obj'] = prev_obj; r['prev_solver'] = PREV_SOLVER; r['new_rows'] = int(len(new_rows)); r['new_cols'] = int(len(newc))
    # objective is invariant under the column permutation
    return r


def run_inc2(d, dprev, threads):
    """two-stage warm start on one Highs object: (0) re-run on the unmodified model to measure how clean
    the carried basis is; (1) append the new rows only and run DUAL simplex (the old basis with the new
    logicals basic is dual feasible, so no phase 1); (2) append the new columns and run PRIMAL simplex
    (the basis from (1) is primal feasible for the extended model).  Timed: (1) + (2); (0) reported."""
    h = hs_new(threads, PREV_SOLVER, log=HLOG); h.passModel(hs_model(dprev)); t = time.time(); h.run(); t_prev = time.time() - t
    prev_obj = float(h.getInfo().objective_function_value)
    print(f"    [inc2] prev solved by {PREV_SOLVER}: obj={prev_obj:.9f} t={t_prev:.0f}s status={str(h.getModelStatus()).split('.')[-1]}", flush=True)
    h.setOptionValue('solver', 'simplex'); h.setOptionValue('simplex_strategy', 1)
    it_base = int(h.getInfo().simplex_iteration_count)
    t = time.time(); h.run(); t_zero = time.time() - t; it_zero = int(h.getInfo().simplex_iteration_count) - it_base
    print(f"    [inc2] zero-change re-run: {it_zero} simplex iterations, {t_zero:.1f}s, obj={float(h.getInfo().objective_function_value):.9f}", flush=True)
    n0 = dprev['M'].shape[1]; n1 = d['M'].shape[1]; nr1 = d['M'].shape[0]; nl = int(d['nl']); assert nl == int(dprev['nl'])
    kprev = dprev['keys']; kd = d['keys']; pos = {k: i for i, k in enumerate(kd)}
    old_rows = np.array([pos[k] for k in kprev]); assert len(old_rows) == len(kprev)
    new_rows = np.array(sorted(set(range(nr1)) - set(old_rows.tolist())), dtype=np.int64)
    newc = np.arange(n0 - nl, n1 - nl); M = d['M']
    t0 = time.time()
    # stage 1: rows only (entries restricted to the old columns, in h's order [old points][lambda])
    if len(new_rows):
        Mr = M[new_rows].tocsr()
        colmap = np.full(n1, -1, dtype=np.int64); colmap[:n0 - nl] = np.arange(n0 - nl); colmap[n1 - nl:] = np.arange(n0 - nl, n0)
        keep = colmap[Mr.indices] >= 0
        # rebuild csr without the new-column entries
        import scipy.sparse as sp_
        Mr2 = sp_.csr_matrix((Mr.data[keep], colmap[Mr.indices][keep], np.concatenate([[0], np.cumsum(np.add.reduceat(keep.astype(np.int64), Mr.indptr[:-1]) if Mr.nnz else [])])), shape=(len(new_rows), n0)) if Mr.nnz else sp_.csr_matrix((len(new_rows), n0))
        h.addRows(len(new_rows), np.asarray(d['b'][new_rows], dtype=float), np.full(len(new_rows), INF), Mr2.nnz, Mr2.indptr.astype(np.int32), Mr2.indices.astype(np.int32), Mr2.data.astype(float))
    it1 = int(h.getInfo().simplex_iteration_count); t = time.time(); h.run(); t_rows = time.time() - t; it_rows = int(h.getInfo().simplex_iteration_count) - it1
    obj_rows = float(h.getInfo().objective_function_value)
    print(f"    [inc2] stage 1 (+{len(new_rows)} rows, dual simplex): {it_rows} iterations, {t_rows:.1f}s, obj={obj_rows:.9f}, {str(h.getModelStatus()).split('.')[-1]}", flush=True)
    # stage 2: columns (entries against all rows of h: old rows in prev order, then new rows in new_rows order)
    if len(newc):
        rows_h = np.concatenate([old_rows, new_rows])
        Mc = M[rows_h][:, newc].tocsc()
        h.addCols(len(newc), np.asarray(d['c'][newc], dtype=float), np.asarray(d['lo'][newc], dtype=float), np.asarray(d['hi'][newc], dtype=float),
                  Mc.nnz, Mc.indptr.astype(np.int32), Mc.indices.astype(np.int32), Mc.data.astype(float))
    h.setOptionValue('simplex_strategy', 4)          # primal
    it2 = int(h.getInfo().simplex_iteration_count); t = time.time(); h.run(); t_cols = time.time() - t; it_cols = int(h.getInfo().simplex_iteration_count) - it2
    print(f"    [inc2] stage 2 (+{len(newc)} cols, primal simplex): {it_cols} iterations, {t_cols:.1f}s, obj={float(h.getInfo().objective_function_value):.9f}, {str(h.getModelStatus()).split('.')[-1]}", flush=True)
    r = hs_result(h, t0); r['it'] = it_rows + it_cols; r['t_prev'] = t_prev; r['prev_obj'] = prev_obj; r['it_zero'] = it_zero; r['t_zero'] = t_zero
    r['it_rows'] = it_rows; r['t_rows'] = t_rows; r['obj_rows'] = obj_rows; r['it_cols'] = it_cols; r['t_cols'] = t_cols; r['new_rows'] = int(len(new_rows)); r['new_cols'] = int(len(newc))
    return r


def run_csift(d, threads, init_cols, add_per_pass, tol=1e-7, max_pass=50, solver='ipm', log=lambda s: print(s, flush=True), hlog=False):
    """column sifting (restricted master): solve the LP over a column subset (all rows), price every
    column by its reduced cost c_j - y.A_j under the restricted dual, add the most negative ones,
    repeat until none is below -tol.  Sub-solves are fresh HiGHS runs (IPM by default: the cost of
    an IPM iteration is ~ (columns)^2 here, so a 3k-column master is ~20x cheaper than the 14k one).
    The lambda columns are always in."""
    M, c, b, lo, hi = d['M'], d['c'], d['b'], d['lo'], d['hi']; nr, nc = M.shape; nl = int(d['nl'])
    t0 = time.time()
    active = np.zeros(nc, dtype=bool); active[init_cols] = True; active[nc - nl:] = True
    Mt = M.T.tocsr()
    # feasibility: every row needs at least one active point column with a positive entry (then A x >= b
    # is satisfiable by scaling); cover the uncovered rows greedily by the column hitting most of them
    P = M[:, :nc - nl].tocsr(); Pb = P.copy(); Pb.data = (Pb.data > 0).astype(float)
    cov = Pb @ active[:nc - nl].astype(float); unc = cov <= 0; added_cover = 0
    while unc.any():
        hits = Pb[unc].sum(axis=0).A1; j = int(np.argmax(hits))
        if hits[j] <= 0: break                       # a row with no positive entry at all: the full LP is infeasible too
        active[j] = True; added_cover += 1; unc &= (Pb[:, j].toarray().ravel() <= 0)
    log(f"      csift start: {int(active.sum())} columns ({added_cover} added to cover rows)")
    passes = 0; t_solve = 0.0; t_price = 0.0; its = 0; hist = []
    while True:
        cols = np.nonzero(active)[0]
        e = dict(M=M[:, cols].tocsr(), c=c[cols], b=b, lo=lo[cols], hi=hi[cols])
        h = hs_new(threads, solver, log=hlog); h.passModel(hs_model(e))
        t = time.time(); h.run(); t_solve += time.time() - t; passes += 1
        st = h.getModelStatus()
        if st != highspy.HighsModelStatus.kOptimal:
            return dict(obj=np.nan, t=time.time() - t0, it=its, status=str(st).split('.')[-1], passes=passes)
        info = h.getInfo(); its += max(int(info.simplex_iteration_count), int(info.ipm_iteration_count))
        sol = h.getSolution(); y = np.array(sol.row_dual); obj = float(info.objective_function_value)
        t = time.time(); rc = c - (Mt @ y); t_price += time.time() - t
        rc[active] = np.inf
        viol = np.nonzero(rc < -tol)[0]
        hist.append((int(active.sum()), obj, int(len(viol)), float(rc[viol].min()) if len(viol) else 0.0))
        log(f"      csift pass {passes}: cols={hist[-1][0]} obj={obj:.7f} neg_rc={hist[-1][2]} min_rc={hist[-1][3]:.2e} it={its} t_solve={t_solve:.0f}s t_price={t_price:.1f}s")
        if len(viol) == 0 or passes >= max_pass: break
        pick = viol[np.argsort(rc[viol])[:add_per_pass]]; active[pick] = True
    r = hs_result(h, t0); r['it'] = its; r['passes'] = passes; r['cols_final'] = int(active.sum()); r['t_solve'] = t_solve; r['t_price'] = t_price; r['hist'] = hist
    if r['z'] is not None:
        z = np.zeros(nc); z[np.nonzero(active)[0]] = r['z']; r['z'] = z
    return r


def support_cols(dprev, d, thr=1e-12):
    """columns of d that carried weight in dprev's recorded solution (by index: columns are appended, never reordered)"""
    if len(dprev['x']) == 0: return np.zeros(0, dtype=np.int64)
    nl = int(dprev['nl']); x = dprev['x']
    return np.nonzero(x > thr)[0]


def run_sift(d, threads, init_rows, add_per_pass, tol=1e-9, max_pass=200, log=lambda s: print(s, flush=True), hlog=False):
    """row sifting on one Highs object: rows appended, basis kept between passes"""
    M, b = d['M'], d['b']; nr = M.shape[0]
    t0 = time.time()
    active = np.zeros(nr, dtype=bool); active[init_rows] = True
    order = np.nonzero(active)[0]                     # rows in h, in order
    h = hs_new(threads, 'simplex', log=hlog); h.passModel(hs_model(d, order))
    passes = 0; t_solve = 0.0; t_scan = 0.0; its = 0; hist = []
    while True:
        t = time.time(); h.run(); t_solve += time.time() - t; passes += 1
        st = h.getModelStatus()
        if st != highspy.HighsModelStatus.kOptimal:
            return dict(obj=np.nan, t=time.time() - t0, it=its, status=str(st).split('.')[-1], passes=passes)
        its += int(h.getInfo().simplex_iteration_count)
        z = np.array(h.getSolution().col_value)
        t = time.time(); slack = M @ z - b; t_scan += time.time() - t
        viol = np.nonzero((slack < -tol) & ~active)[0]
        hist.append((int(active.sum()), float(h.getInfo().objective_function_value), int(len(viol)), float(slack[viol].min()) if len(viol) else 0.0))
        log(f"      sift pass {passes}: rows={hist[-1][0]} obj={hist[-1][1]:.7f} viol={hist[-1][2]} worst={hist[-1][3]:.2e} simplex_it={int(h.getInfo().simplex_iteration_count)} t_solve={t_solve:.0f}s t_scan={t_scan:.0f}s")
        if len(viol) == 0 or passes >= max_pass: break
        pick = viol[np.argsort(slack[viol])[:add_per_pass]]
        Mr = M[pick].tocsr()
        h.addRows(len(pick), np.asarray(b[pick], dtype=float), np.full(len(pick), INF), Mr.nnz, Mr.indptr.astype(np.int32), Mr.indices.astype(np.int32), Mr.data.astype(float))
        active[pick] = True; order = np.concatenate([order, pick])
    r = hs_result(h, t0); r['it'] = its; r['passes'] = passes; r['rows_final'] = int(active.sum()); r['t_solve'] = t_solve; r['t_scan'] = t_scan; r['hist'] = hist
    # duals mapped back to the full row set (inactive rows: 0)
    if r['y'] is not None:
        y = np.zeros(nr); y[order] = r['y']; r['y'] = y
    return r


def tight_rows(dprev, d, slack_thr):
    """rows of d that were tight in dprev's recorded solution (by key); falls back to y > 0"""
    if len(dprev['y']) == 0: return np.zeros(0, dtype=np.int64)
    z = np.concatenate([dprev['x'], dprev['lam']]) if len(dprev['x']) else None
    tight = dprev['y'] > 1e-12
    if z is not None and len(z) == dprev['M'].shape[1]:
        tight |= (dprev['M'] @ z - dprev['b']) < slack_thr
    keys = set(k for k, tg in zip(dprev['keys'], tight) if tg)
    return np.array([i for i, k in enumerate(d['keys']) if k in keys], dtype=np.int64)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dump'); ap.add_argument('--warm-from', default=None)
    ap.add_argument('--methods', default='ipm,ds,hs-ds,hs-ipm,pdlp,sift'); ap.add_argument('--reps', type=int, default=2)
    ap.add_argument('--threads', type=int, default=12); ap.add_argument('--pdlp-tol', type=float, default=1e-7)
    ap.add_argument('--sift-init', type=int, default=20000); ap.add_argument('--sift-add', type=int, default=5000); ap.add_argument('--sift-slack', type=float, default=2e-2)
    ap.add_argument('--csift-init', type=int, default=2000); ap.add_argument('--csift-add', type=int, default=1000); ap.add_argument('--csift-extra', type=int, default=0, help='also start with the last N point columns (the round\'s new columns)')
    ap.add_argument('--ref', type=float, default=None, help='reference objective (default: ipm result if run, else recorded)')
    ap.add_argument('--json', default=None); ap.add_argument('--log', action='store_true')
    ap.add_argument('--note', default='', help='what else was running (recorded with every timing row)')
    ap.add_argument('--prev-solver', default='simplex', choices=['simplex', 'ipm'], help='warm/inc: how --warm-from is solved for the starting basis (ipm = IPM + crossover)')
    a = ap.parse_args()
    global PREV_SOLVER, HLOG; PREV_SOLVER = a.prev_solver; HLOG = a.log
    d = load(a.dump); dprev = load(a.warm_from) if a.warm_from else None
    print(f"dump {a.dump}: {describe(d)}");
    if dprev is not None: print(f"warm-from {a.warm_from}: {describe(dprev)}")
    ref = a.ref if a.ref is not None else float(d['val'])
    out = []
    for meth in a.methods.split(','):
        for rep in range(a.reps):
            if meth == 'ipm': r = run_scipy(d, 'highs-ipm', a.threads)
            elif meth == 'ds': r = run_scipy(d, 'highs-ds', a.threads)
            elif meth == 'hs-ds': r = run_hs(d, 'simplex', a.threads, a.pdlp_tol, log=a.log)
            elif meth == 'hs-ipm': r = run_hs(d, 'ipm', a.threads, a.pdlp_tol, log=a.log)
            elif meth == 'pdlp': r = run_hs(d, 'pdlp', a.threads, a.pdlp_tol, log=a.log)
            elif meth == 'warm': r = run_warm(d, dprev, a.threads)
            elif meth == 'inc': r = run_inc(d, dprev, a.threads)
            elif meth == 'inc2': r = run_inc2(d, dprev, a.threads)
            elif meth in ('csift', 'csift-ds'):
                init = support_cols(dprev, d) if dprev is not None else np.arange(min(a.csift_init, d['M'].shape[1]))
                if len(init) == 0: init = np.arange(min(a.csift_init, d['M'].shape[1]))
                if a.csift_extra: init = np.unique(np.concatenate([init, np.arange(max(0, d['M'].shape[1] - int(d['nl']) - a.csift_extra), d['M'].shape[1] - int(d['nl']))]))
                r = run_csift(d, a.threads, init, a.csift_add, solver='ipm' if meth == 'csift' else 'simplex', hlog=a.log)
            elif meth in ('sift', 'sift-warm'):
                init = tight_rows(dprev, d, a.sift_slack) if dprev is not None else np.arange(min(a.sift_init, d['M'].shape[0]))
                if len(init) == 0: init = np.arange(min(a.sift_init, d['M'].shape[0]))
                r = run_sift(d, a.threads, init, a.sift_add, hlog=a.log)
            else: raise SystemExit(f"unknown method {meth}")
            if meth == 'ipm' and a.ref is None and rep == 0 and np.isfinite(r['obj']): ref = r['obj']
            extra = {k: v for k, v in r.items() if k not in ('z', 'y', 'obj', 't', 'it', 'status', 'hist')}
            dev = (r['obj'] - ref) / max(1.0, abs(ref)) if np.isfinite(r['obj']) else np.nan
            ldavg = os.getloadavg()
            print(f"  {meth:9s} rep{rep}  obj={r['obj']:.9f}  dev={dev:+.2e}  t={r['t']:.1f}s  it={r['it']}  {r['status']}  {extra}  load={ldavg[0]:.1f}/{ldavg[1]:.1f} threads={a.threads} affinity={sorted(os.sched_getaffinity(0))} note='{a.note}'", flush=True)
            if 'hist' in r and a.log:
                for hrow in r['hist']: print(f"      sift pass: rows={hrow[0]} obj={hrow[1]:.7f} viol={hrow[2]} worst={hrow[3]:.2e}")
            out.append(dict(method=meth, rep=rep, obj=r['obj'], dev=dev, t=r['t'], it=r['it'], status=r['status'], load=ldavg, threads=a.threads, affinity=sorted(os.sched_getaffinity(0)), note=a.note, **{k: (v if not isinstance(v, np.generic) else v.item()) for k, v in extra.items()}, hist=r.get('hist')))
    if a.json: json.dump(dict(dump=a.dump, warm_from=a.warm_from, ref=ref, rows=int(d['M'].shape[0]), cols=int(d['M'].shape[1]), nnz=int(d['M'].nnz), results=out), open(a.json, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
