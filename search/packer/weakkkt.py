#!/usr/bin/env python3
"""Redo exactsolve's multiplier + second-order step for a weak census minimum, with the multiplier computed in mp (10-15).

Why: in the s(110) 'PSD with k zero modes' weak rows, exactsolve's max_support counts a contact as load-bearing only if
its lambda > 0.01 tau0 = 1e-8; contacts that need lambda ~1e-9 (typically a pair closed by the active-set step) land in the
weakly active list, full_lambda_lp over the rest is infeasible, and second_order runs with lambda = 0 (so 'PSD', every
null direction a zero mode).  Here: all closed contacts at the exact point (load-bearing + weakly active, as captured by
weakmodes.py), Lawson-Hanson NNLS for J^T lambda = e_S with least squares by iterative refinement in mp (exact to ~1e-70).
Residual ~0 => KKT; then exactsolve.second_order with this lambda (positive part = equalities, zero part = weakly active,
critical cone) and probe_flat on any zero modes.  Residual > 0 => no multiplier: see weakfarkas.py (descent certificate).

  python3 weakkkt.py DIR        (DIR/modes.pkl; prints one JSON line)
"""
import contextlib, io, json, os, pickle, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exact')); sys.path.insert(0, HERE)
import numpy as np
import scipy.linalg as sla
from scipy.optimize import nnls, linprog
from mpmath import mp, mpf
import exactsolve as E
from weakfarkas import jac_mp, refine


def mp_nnls(P, A, V, maxit=400):
    """min |J_A^T lambda - e_S| over lambda >= 0, active set in mp.  Returns (lam dict k -> mp, residual vector rm, J r)."""
    nS = 3 * P.n
    J, _ = E.jacobian(P, A); JV = J[:, V]
    e = np.zeros(len(V)); e[V.index(nS)] = 1.0
    lam0, _ = nnls(JV.T, e, maxiter=50 * len(A))
    _, lam0, _ = refine(JV, e, lam0)
    Ps = set(np.where(lam0 > 0)[0].tolist())
    Jm = jac_mp(P, A, V)
    em = [mpf(0)] * len(V); em[V.index(nS)] = mpf(1)
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
        Jr = [sum((val * rm[c] for c, val in Jm[k].items()), mpf(0)) for k in range(len(A))]
        cand = [k for k in range(len(A)) if k not in Ps and Jr[k] > mpf(10) ** -60]
        if not cand:
            break
        Ps.add(max(cand, key=lambda k: Jr[k]))
    return dict(zip(Pi, lm)), rm, Jr


def main(d):
    mp.dps = 80
    D = pickle.load(open(f'{d}/modes.pkl', 'rb'))
    P, V = D['P'], D['V']
    Vset = set(V)
    inV = lambda c: all(v in Vset for v in [3 * c[1]] + ([3 * c[3]] if c[0] == 'C' else []))
    A = [c for c in D['A'] + [w for w in D['weak'] if w not in set(D['A'])] if inV(c)]
    # (1) strictly positive multiplier over every closed contact (max t, lambda_k >= t), f64 LP, then an exact (mp)
    #     min-norm correction of J^T lambda = e_S; positivity must survive the correction.
    J0, _ = E.jacobian(P, A)
    _, supp, _ = E.max_support(J0, tau0=1e-10)           # contacts that can carry force (default tau0 1e-6 misses 1e-9 ones)
    weak0 = [A[k] for k in range(len(A)) if not supp[k]]
    A = [A[k] for k in range(len(A)) if supp[k]]
    J, _ = E.jacobian(P, A); JV = J[:, V]; K = len(A)
    e = np.zeros(len(V)); e[V.index(3 * P.n)] = 1.0
    r = linprog(np.r_[np.zeros(K), -1.0], A_ub=np.hstack([-np.eye(K), np.ones((K, 1))]), b_ub=np.zeros(K),
                A_eq=np.hstack([JV.T, np.zeros((len(V), 1))]), b_eq=e, bounds=[(0, None)] * K + [(None, 1)],
                method='highs', options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
    out = dict(dir=os.path.basename(d), S=str(D['S'])[:22], n_closed=K, lp_status=int(r.status))
    if r.status == 0 and r.x[-1] > 1e-13:
        Jm = jac_mp(P, A, V)
        lm = [mpf(float(x)) for x in r.x[:K]]
        Gp = np.linalg.pinv(JV.T @ JV)
        for _ in range(8):
            rm = [mpf(0)] * len(V); rm[V.index(3 * P.n)] = mpf(1)
            for k in range(K):
                for c, val in Jm[k].items():
                    rm[c] -= lm[k] * val
            nr = max(abs(x) for x in rm)
            if nr < mpf(10) ** -70:
                break
            y = Gp @ np.array([float(x / nr) for x in rm])     # lambda += J y (min-norm)
            dl = JV @ y
            lm = [l + mpf(float(x)) * nr for l, x in zip(lm, dl)]
        out.update(t_lp=float(r.x[-1]), residual=float(nr), lam_min_mp=float(min(lm)), n_weak0=len(weak0))
        if nr < mpf(10) ** -50 and min(lm) > 0:
            Ap, lp, wk = A, [float(x) for x in lm], weak0
        else:
            r = None
    else:
        r = None
    if r is None:
        A = A + weak0
        lam, rm, Jr = mp_nnls(P, A, V)
        res = float(max(abs(x) for x in rm)); out['residual_nnls'] = res
        if res > 1e-50:
            out['verdict'] = 'no multiplier (not KKT): run weakfarkas'
            print(json.dumps(out)); return
        Ap = [A[k] for k in sorted(lam) if lam[k] > 0]
        lp = [float(lam[k]) for k in sorted(lam) if lam[k] > 0]
        wk = [A[k] for k in range(len(A)) if k not in lam or lam[k] <= 0]
    squares = sorted({c[1] for c in Ap} | {c[3] for c in Ap if c[0] == 'C'})
    Vp = [3 * i + k for i in squares for k in range(3)] + [3 * P.n]
    with contextlib.redirect_stdout(io.StringIO()):
        so = E.second_order(P, Ap, lp, Vp, squares, weak=wk)
    zm = so.pop('_zero_modes', None)
    out.update(n_support=len(Ap), lam_min=min(lp), lam_max=max(lp), n_weak=len(wk), dim_null=so['dim_null'],
               n_zero=so['n_zero'], n_neg=so['n_neg'], eig=[float(f'{x:.3g}') for x in (so.get('eig') or [])[:4]],
               eig_min_pos=so.get('eig_min_positive'), status=so['status'][:90],
               descent=[{k: v for k, v in x.items() if k != 'squares'} for x in so.get('descent', [])][:3])
    if zm is not None and zm.shape[1]:
        with contextlib.redirect_stdout(io.StringIO()):
            nre, wr, wc, wds = E.probe_flat(P, Ap, Vp, zm, log_=lambda *a: None)
        out.update(realised=nre, probe_res=wr, probe_corr=wc, probe_dS=wds)
        if so['n_neg'] == 0 and nre == zm.shape[1]:
            out['status'] = f'local minimum of S, strict modulo {zm.shape[1]} exact flat motions'
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
