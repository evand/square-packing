"""Load-bearing set and multiplier check without tiny-multiplier LPs (10-15, for exactsolve's fallback).

At a KKT point of min S with closed contacts A (g >= 0), the critical cone is C = {d : J_A d >= 0, d_S = 0}.  Its implicit
equalities E (rows with J d = 0 on all of C) are exactly the contacts that carry force in some multiplier; the rest can
open.  Computed by homogeneous LPs (scale-free, unlike max_support / full_lambda_lp, whose tolerances miss contacts with
lambda ~1e-9 .. 1e-11, as in 47 s(110) census points):
    repeat: max sum_{k undecided} t_k  s.t.  J d >= t (undecided rows), J d >= 0 (opened rows), 0 <= t <= 1, |d| <= 1;
            rows with t_k > 1e-9 are opened;  until none opens.
Rows with a positive exact (mp) multiplier are imposed as equalities through an orthonormal null-space basis
(smaller, better-conditioned LPs; HiGHS 1.8 returns status Unknown on the full LP for some s(110) points).
mp_nnls: min |J_E^T lambda - e_S| over lambda >= 0, Lawson-Hanson with least squares by iterative refinement (f64 solves,
mp residuals), so the KKT residual is exact to ~1e-70 even when some lambda are ~1e-11.
"""
import numpy as np
import scipy.linalg as sla
from scipy.optimize import linprog, nnls
from mpmath import mpf


def _lp(B, R, opened, rng=None):
    """max sum t over rows R of B y >= t (others in `opened`: B y >= 0), |y| <= 1, 0 <= t <= 1.  Returns t or None."""
    n = B.shape[1]; k = len(R)
    rows_t = B[R]; rows_0 = B[opened] if len(opened) else np.zeros((0, n))
    w = np.ones(k) if rng is None else rng.uniform(0.5, 1.0, k)
    A_ub = np.vstack([np.hstack([-rows_t, np.eye(k)]), np.hstack([-rows_0, np.zeros((len(rows_0), k))])])
    b_ub = np.zeros(len(A_ub))
    tol = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}
    for meth, extra in (('highs', {}), ('highs-ds', {}), ('highs-ipm', {}), ('highs-ds', {'presolve': False}),
                        ('highs-ipm', {'presolve': False})):       # HiGHS 1.8 postsolve sometimes ends 'Unknown'
        r = linprog(np.r_[np.zeros(n), -w], A_ub=A_ub, b_ub=b_ub, bounds=[(-1, 1)] * n + [(0, 1)] * k, method=meth,
                    options=dict(tol, **extra))
        if r.x is not None:
            return r.x[n:]
    return None


def cone_support(JV, seed, open_tol=1e-9):
    """JV: closed contacts x variables, S column last; seed: mask of rows known to be implicit equalities (positive exact
    multiplier).  Returns (E mask, min opening of opened rows) or (None, msg)."""
    K, m = JV.shape
    Jd = JV[:, :-1].copy(); Jd[np.abs(Jd) < 1e-13] = 0.0
    tries = [np.asarray(seed, bool)]
    for P in tries:
        N = sla.null_space(Jd[P]) if P.any() else np.eye(m - 1)
        if N.shape[1] == 0:
            return P.copy(), None
        B = Jd @ N; B[np.abs(B) < 1e-13] = 0.0
        undecided = np.where(~P)[0].tolist(); opened = []; tmin = None; ok = True
        rng = None
        for it in range(50):
            if not undecided:
                break
            t = _lp(B, undecided, opened, rng)
            if t is None and rng is None:
                rng = np.random.default_rng(it); t = _lp(B, undecided, opened, rng)
            if t is None:                                  # last resort: one small LP per undecided row
                t = np.zeros(len(undecided))
                for q, k in enumerate(undecided):
                    tq = _lp(B, [k], opened + [j for j in undecided if j != k])
                    if tq is None:
                        ok = False; break
                    t[q] = tq[0]
                if not ok:
                    break
            new = [undecided[q] for q in range(len(undecided)) if t[q] > open_tol]
            if not new:
                break
            tmin = min([tmin] + [t[q] for q in range(len(undecided)) if t[q] > open_tol] if tmin is not None
                       else [t[q] for q in range(len(undecided)) if t[q] > open_tol])
            opened += new
            undecided = [k for k in undecided if k not in set(new)]
        if ok:
            E = np.ones(K, bool); E[opened] = False
            return E, tmin
    return None, 'cone LP failed'


def jac_mp(P, A, V, eval_contact):
    """Contact gradients at the exact point in mp: one dict {column in V: value} per contact."""
    C, Sn = P.trig(); H = mpf(1) / 2
    col = {v: k for k, v in enumerate(V)}
    rows = []
    for c in A:
        _, gr, _ = eval_contact(c, P.X, P.Y, C, Sn, P.S, H, order=1)
        row = {}
        for v, val in gr:
            if v in col:
                row[col[v]] = row.get(col[v], mpf(0)) + val
        rows.append(row)
    return rows


def mp_nnls(JV, Jm, maxit=400, tol=mpf(10) ** -60):
    """Lawson-Hanson NNLS for JV^T lambda = e_S (S = last column), mp residuals.  Returns (lam list mp, residual max mp)."""
    K, m = JV.shape
    e = np.zeros(m); e[-1] = 1.0
    lam0, _ = nnls(JV.T, e, maxiter=100 * K)
    Ps = set(np.where(lam0 > 0)[0].tolist()) or {0}
    em = [mpf(0)] * m; em[-1] = mpf(1)
    lm, Pi, rm = [], [], em
    for _ in range(maxit):
        Pl = sorted(Ps)
        _, Rq, piv = sla.qr(JV[Pl].T, pivoting=True, mode='economic')
        dg = np.abs(np.diag(Rq)); rk = int((dg > 1e-10 * dg[0]).sum())
        Pi = sorted(Pl[k] for k in piv[:rk])
        G = JV[Pi] @ JV[Pi].T
        lm = [mpf(0)] * len(Pi)
        for _ in range(14):
            rm = list(em)
            for q, k in enumerate(Pi):
                for c, val in Jm[k].items():
                    rm[c] -= lm[q] * val
            gm = [sum((val * rm[c] for c, val in Jm[k].items()), mpf(0)) for k in Pi]
            ng = max(abs(x) for x in gm)
            if ng < mpf(10) ** -70:
                break
            dl = np.linalg.solve(G, np.array([float(x / ng) for x in gm]))
            lm = [l + mpf(float(x)) * ng for l, x in zip(lm, dl)]
        Ps = set(Pi)
        if min(lm) < 0:
            Ps.discard(Pi[min(range(len(Pi)), key=lambda q: lm[q])]); continue
        Jr = [sum((val * rm[c] for c, val in Jm[k].items()), mpf(0)) for k in range(K)]
        cand = [k for k in range(K) if k not in Ps and Jr[k] > tol]
        if not cand:
            break
        Ps.add(max(cand, key=lambda k: Jr[k]))
    lam = [mpf(0)] * K
    for q, k in enumerate(Pi):
        lam[k] = lm[q]
    return lam, max(abs(x) for x in rm)
