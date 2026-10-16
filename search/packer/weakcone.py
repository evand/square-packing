#!/usr/bin/env python3
"""Second-order test of a weak census minimum with the critical cone found without a multiplier (10-15).

At a KKT point of min S (all closed contacts A, g >= 0), the critical cone is  C = {d : J_A d >= 0, d_S = 0}  (since
e_S . d = lambda^T J_A d >= 0, with equality iff complementary).  Implicit equalities of C (rows with J d = 0 on all of
C) = the contacts that carry force in some multiplier (E); the rest can open (weakly active).  One homogeneous LP finds E
(scale-free, unlike the multiplier LPs: the s(110) weak rows need lambda ~1e-9 .. 1e-11, below HiGHS tolerances):
    max sum t   s.t.  J d >= t,  0 <= t <= 1,  d_S = 0,  |d|_inf <= 1.
Squares with no E contact (wedged rattlers) are dropped together with their contacts: a relaxation, so a strict local
minimum of the relaxed problem is one of the full problem.  Then exactsolve.second_order(P, E, lambda_E, V, squares,
weak = A \\ E among kept squares) with lambda >= 0 on E from NNLS (f64; only its curvature contribution is used), and
exactsolve.probe_flat on the zero modes.  KKT itself is checked: NNLS residual of J_E^T lambda = e_S.

  python3 weakcone.py DIR            (DIR/modes.pkl from weakmodes.py; one JSON line)
"""
import contextlib, io, json, os, pickle, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exact')); sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import linprog, nnls
from mpmath import mp
import exactsolve as E


def sqs(c):
    return [c[1]] + ([c[3]] if c[0] == 'C' else [])


def main(d):
    mp.dps = 80
    D = pickle.load(open(f'{d}/modes.pkl', 'rb'))
    P = D['P']; nS = 3 * P.n
    A = D['A'] + [w for w in D['weak'] if w not in set(D['A'])]
    V = sorted({3 * i + k for c in A for i in sqs(c) for k in range(3)}) + [nS]
    J, _ = E.jacobian(P, A); JV = J[:, V]; K = len(A); m = len(V)
    # implicit equalities of C = {J d >= 0, d_S = 0}
    c = np.r_[np.zeros(m), -np.ones(K)]
    Aub = np.hstack([-JV, np.eye(K)])
    bnd = [(-1, 1)] * (m - 1) + [(0, 0)] + [(0, 1)] * K
    import conelp
    JVc = np.delete(JV, m - 1, axis=1); JVc[np.abs(JVc) < 1e-12] = 0.0
    t = conelp.solve(JVc)
    if t is None:                                          # HiGHS 1.8 'Unknown': retry under the nix scipy (newer HiGHS)
        import subprocess, tempfile
        with tempfile.TemporaryDirectory() as td:
            np.savez(f'{td}/in.npz', JV=JVc)
            subprocess.run([f'{HERE}/ipenv.sh', 'python3', f'{HERE}/conelp.py', f'{td}/in.npz', f'{td}/t.npy'],
                           check=True, capture_output=True)
            t = np.load(f'{td}/t.npy')
        if not len(t):
            print(json.dumps(dict(dir=os.path.basename(d), status='cone LP failed'))); return
        out_lp = 'nix'
    else:
        out_lp = 'system'
    Eset = [k for k in range(K) if t[k] < 1e-7]
    out = dict(dir=os.path.basename(d), S=str(D['S'])[:22], n_closed=K, n_E=len(Eset), t_gap=float(np.sort(t)[len(Eset)]) if len(Eset) < K else None,
               t_max_on_E=float(t[Eset].max()) if Eset else None, lp=out_lp)
    # multiplier on E
    e = np.zeros(m); e[-1] = 1.0
    lam, rn = nnls(JV[Eset].T, e, maxiter=100 * len(Eset))
    out['kkt_resid'] = float(rn)
    keep = sorted({i for k in Eset for i in sqs(A[k])})
    ks = set(keep)
    Ep = [A[k] for k in Eset]
    weak = [A[k] for k in range(K) if k not in set(Eset) and all(i in ks for i in sqs(A[k]))]
    Vp = [3 * i + q for i in keep for q in range(3)] + [nS]
    with contextlib.redirect_stdout(io.StringIO()):
        so = E.second_order(P, Ep, [float(x) for x in lam], Vp, keep, weak=weak)
    zm = so.pop('_zero_modes', None)
    out.update(dropped_squares=sorted(set(range(P.n)) - ks), n_weak=len(weak), dim_null=so['dim_null'], n_zero=so['n_zero'],
               n_neg=so['n_neg'], eig=[float(f'{x:.3g}') for x in (so.get('eig') or [])[:4]],
               eig_min_pos=so.get('eig_min_positive'), status=so['status'][:90],
               descent=[{k: v for k, v in x.items() if k != 'squares'} for x in so.get('descent', [])][:3])
    if zm is not None and zm.shape[1]:
        with contextlib.redirect_stdout(io.StringIO()):
            nre, wr, wc, wds = E.probe_flat(P, Ep, Vp, zm, log_=lambda *a: None)
        out.update(realised=nre, probe_res=wr, probe_corr=wc, probe_dS=wds)
        if so['n_neg'] == 0 and nre == zm.shape[1]:
            out['status'] = f'local minimum of S, strict modulo {zm.shape[1]} exact flat motions'
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
