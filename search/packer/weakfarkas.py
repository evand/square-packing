#!/usr/bin/env python3
"""Exact non-minimality certificate for a census point whose contact set has no multiplier lambda >= 0 (10-15).

exactsolve's 'PSD with k zero modes' weak verdicts (s(110), 47 rows) come from full_lambda_lp failing (no lambda >= 0
balances the load for the closed contact set A at the exact point) and the second-order test then running with
lambda = 0.  By Farkas, the NNLS residual r = e_S - J_A^T lambda* (lambda* = argmin_{lambda >= 0} |J_A^T lambda - e_S|)
gives a first-order descent d = -r:  J_A d >= 0, dS = -|r|^2 < 0 (slopes ~1e-9: below the f64 MILP's tolerance).

For each t: start at v* + t d^ (d^ = d / max|d|), hold S at S* + t d^_S < S*, Gauss-Newton (mp residuals) the contacts of A
with J d^ < 1e-6 back to zero over the other variables, then check every pair (SAT) and wall in mp and write an exact
rational certificate (exactsolve.certificate: scale by 1 + eps about the origin, round, verify in exact arithmetic).
A verified certificate at S' < S* within distance ~t of the point, for t -> 0, shows the point is not a local minimum.

  python3 weakfarkas.py DIR [--ts 1e-4,1e-6,1e-8,1e-10]       (DIR/modes.pkl from weakmodes.py; writes DIR/farkas_*.cert)
"""
import argparse, contextlib, io, json, os, pickle, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exact')); sys.path.insert(0, HERE)
import numpy as np
import scipy.linalg as sla
from scipy.optimize import nnls
import mpmath
from mpmath import mp, mpf
import exactsolve as E


def min_clearance(P):
    C, Sn = P.trig(); H = mpf(1) / 2
    m, w = None, None
    for i, j in E.near_pairs(P, 0.01):
        g = E.sat_gap(i, j, P.X, P.Y, C, Sn, H)
        if m is None or g < m: m, w = g, (i, j)
    for i in range(P.n):
        for k, g in enumerate(E.wall_gaps(i, P.X, P.Y, C, Sn, P.S, H)):
            if g < m: m, w = g, (i, 'LRBT'[k])
    min_clearance.where = w
    return m


def project(Q, rows, Vf, maxit=60):
    """Gauss-Newton onto {g_rows = 0} over variables Vf (S fixed), mp residuals, independent subset by pivoted QR."""
    H = mpf(1) / 2
    J0, _ = E.jacobian(Q, rows)
    _, Rq, piv = sla.qr(J0[:, Vf].T, pivoting=True, mode='economic')
    dg = np.abs(np.diag(Rq)); rk = int((dg > 1e-9 * dg[0]).sum()) if len(dg) else 0
    B = [rows[k] for k in sorted(piv[:rk].tolist())]
    hist = []
    for it in range(maxit):
        C, Sn = Q.trig()
        g = [E.eval_contact(c, Q.X, Q.Y, C, Sn, Q.S, H, order=1)[0] for c in B]
        nr = max(abs(x) for x in g); hist.append(float(nr))
        if nr < mpf(10) ** (-(mp.dps - 12)) or (it > 6 and hist[-1] > 0.5 * hist[-5]):
            break
        sc = nr if nr < mpf(10) ** -200 else 1
        J, _ = E.jacobian(Q, B)
        d, *_ = np.linalg.lstsq(J[:, Vf], -np.array([float(x / sc) for x in g]), rcond=1e-12)
        for k, v in enumerate(Vf):
            Q.add(v, mpf(float(d[k])) * sc)
    # redundant rows (not in B): residual after the projection
    C, Sn = Q.trig()
    allr = max(abs(E.eval_contact(c, Q.X, Q.Y, C, Sn, Q.S, H, order=1)[0]) for c in rows)
    return hist, float(allr), len(B)


def jac_mp(P, A, V):
    """Contact gradients at the exact point in mp: one dict {column in V: value} per contact."""
    C, Sn = P.trig(); H = mpf(1) / 2
    col = {v: k for k, v in enumerate(V)}
    rows = []
    for c in A:
        _, gr, _ = E.eval_contact(c, P.X, P.Y, C, Sn, P.S, H, order=1)
        row = {}
        for v, val in gr:
            if v in col:
                row[col[v]] = row.get(col[v], mpf(0)) + val
        rows.append(row)
    return rows


def refine(JV, e, lam, tol=1e-14, maxit=2000):
    """Lawson-Hanson active-set polish of an NNLS solution with exact least-squares solves on the support: returns
    r = e - JV^T lam with lam >= 0, (JV r)_k ~ 0 on the support and <= tol elsewhere (Farkas: d = -r is a descent)."""
    Ps = set(np.where(lam > 0)[0].tolist())
    for it in range(maxit):
        idx = sorted(Ps)
        lP, *_ = np.linalg.lstsq(JV[idx].T, e, rcond=1e-13) if idx else (np.zeros(0),)
        if idx and lP.min() < 0:
            Ps.discard(idx[int(np.argmin(lP))]); continue
        lam = np.zeros(JV.shape[0]); lam[idx] = lP
        r = e - JV.T @ lam
        w = JV @ r
        out = [k for k in range(JV.shape[0]) if k not in Ps]
        if not out:
            break
        k = max(out, key=lambda k: w[k])
        if w[k] <= tol * max(1.0, np.linalg.norm(r)):
            break
        Ps.add(k)
    return r, lam, it


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('dir'); ap.add_argument('--ts', default='1e-4,1e-6,1e-8,1e-10')
    a = ap.parse_args()
    mp.dps = 80
    D = pickle.load(open(f'{a.dir}/modes.pkl', 'rb'))
    P, A, V = D['P'], D['A'] + [c for c in D['weak'] if c not in set(D['A'])], D['V']
    Vset = set(V)
    A = [c for c in A if all(v in Vset for v in [3 * c[1], 3 * c[1] + 1, 3 * c[1] + 2] + ([3 * c[3], 3 * c[3] + 1, 3 * c[3] + 2] if c[0] == 'C' else []))] or A
    nS = 3 * P.n
    J, _ = E.jacobian(P, A); JV = J[:, V]
    e = np.zeros(len(V)); e[V.index(nS)] = 1.0
    lam, rn = nnls(JV.T, e, maxiter=50 * len(A))
    r, lam, it = refine(JV, e, lam)
    # mp stage: support rows (independent subset), lambda by iterative refinement (f64 solves, mp residuals)
    Ps = set(np.where(lam > 0)[0].tolist())
    Jm = jac_mp(P, A, V)                                   # list of mp rows (dicts v-index -> value)
    em = [mpf(0)] * len(V); em[V.index(nS)] = mpf(1)
    lh = 0
    for lh in range(400):                                  # Lawson-Hanson in mp (least squares by iterative refinement)
        Pl = sorted(Ps)
        _, Rq, piv = sla.qr(JV[Pl].T, pivoting=True, mode='economic')
        dg = np.abs(np.diag(Rq)); rk = int((dg > 1e-10 * dg[0]).sum())
        Pi = sorted(Pl[k] for k in piv[:rk])
        JPf = JV[Pi]; G = JPf @ JPf.T
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
    rnorm2 = sum(x * x for x in rm)
    rmax = max(abs(x) for x in rm)
    dm = [-x / rmax for x in rm]                            # descent direction, max-norm 1 (mp)
    Jd = np.array([float(-x / rmax) for x in Jr])
    out = dict(dir=os.path.basename(a.dir), S=str(D['S'])[:22], resid=float(mpmath.sqrt(rnorm2)), rS=float(rm[V.index(nS)]),
               rnorm2=float(rnorm2), lam_min=float(min(lm)), dS_unit=float(dm[V.index(nS)]), Jd_min=float(Jd.min()),
               normal_eq=float(ng), support=len(Ps), indep=len(Pi), lh_iters=lh, runs=[])
    if dm[V.index(nS)] >= 0 or min(lm) < 0:
        out['verdict'] = 'no descent direction'
        print(json.dumps(out)); return
    tight = [A[k] for k in range(len(A)) if Jd[k] < 1e-9]
    Vf = [v for v in V if v != nS]
    for t in [float(x) for x in a.ts.split(',')]:
        Q = P.copy()
        for k, v in enumerate(V):
            Q.add(v, mpf(t) * dm[k])
        hist, allr, nB = project(Q, tight, Vf)
        mc = min_clearance(Q)
        gain = P.S - Q.S
        row = dict(t=t, gain=float(gain), proj=hist[-1], redundant=allr, min_clear=float(mc), where=str(min_clearance.where), nB=nB)
        if mc > -mpf(10) ** -60:
            with contextlib.redirect_stdout(io.StringIO()):
                ok, Sp, wall, pair = E.certificate(Q, mpf(10) ** -30, f'{a.dir}/farkas_{t:.0e}.cert')
            row.update(cert_ok=bool(ok), cert_below=bool(Sp < P.S), cert_gain=float(P.S - mpf(Sp.numerator) / Sp.denominator))
        out['runs'].append(row)
    good = [x for x in out['runs'] if x.get('cert_ok') and x.get('cert_below')]
    out['verdict'] = (f'NOT a local minimum: exact certificates below S* at t = {[x["t"] for x in good]}'
                      if good else 'no certificate')
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
