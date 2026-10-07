"""Numerical test of grade-C hypothesis (v): along the flat position family M = z* + F, the force balance
L(m)^T lam = e_S still has a solution with lam > 0.  Rows as in localmin.py (load-bearing incidences; aligned
side-end pairs as midpoint rows); only the load-bearing squares are variables."""
import json, sys, numpy as np, geom
from mpmath import mp, mpf
from scipy.optimize import linprog
import scipy.linalg as sla
mp.dps = 40
def rows_and_L(n, X, Y, C, Sn, S, A):
    N = 3 * n + 1
    rows, ends = [], {}
    for a in A:
        if a[0] == 'W': rows.append(('W', a)); continue
        tau = geom.tangential(a, X, Y, C, Sn, 0.5)
        if abs(abs(tau) - 0.5) < 1e-9: ends.setdefault((a[1], a[3], a[4]), []).append(a)
        else: rows.append(('C', a))
    groups = {}
    for key, inc in ends.items(): groups.setdefault(frozenset((key[0], key[1])), []).append(inc)
    def grad(ct):
        r = np.zeros(N)
        for v, val in geom.eval_contact(ct, X, Y, C, Sn, S, 0.5, order=1)[1]: r[v] += val
        return r
    Ls = [grad(a) for _, a in rows]
    bad = 0
    for g, incs in groups.items():
        for inc in incs[:1]:                       # one line per pair (the other line has the same rows to 1st order)
            if len(inc) != 2: bad += 1
            for a in inc: Ls.append(grad(a))
    return np.array(Ls), bad
def maxmin(L, N, n):
    m = L.shape[0]; es = np.zeros(N); es[3 * n] = 1
    res = linprog(np.r_[np.zeros(m), -1], A_eq=np.c_[L.T, np.zeros(N)], b_eq=es, A_ub=np.c_[-np.eye(m), np.ones(m)],
                  b_ub=np.zeros(m), bounds=[(None, None)] * (m + 1), method='highs')
    return res.x[-1] if res.status == 0 else None, (np.linalg.lstsq(L.T, es, rcond=None)[1] if res.status != 0 else 0)
for n in map(int, sys.argv[1:]):
    Lf = [l.split() for l in open(f'minpoly/solve/n-{n}.exact.txt') if l.strip()]
    X = [float(mpf(l[0])) for l in Lf[1:]]; Y = [float(mpf(l[1])) for l in Lf[1:]]; TH = [float(mpf(l[2])) for l in Lf[1:]]
    S = float(mpf(Lf[0][1]))
    C = [np.cos(np.radians(t)) for t in TH]; Sn = [np.sin(np.radians(t)) for t in TH]
    ct = json.load(open(f'minpoly/solve/n-{n}.contacts.json'))
    A = [tuple(a) for a in ct['load_bearing']]
    N = 3 * n + 1
    L, bad = rows_and_L(n, X, Y, C, Sn, S, A)
    lb = sorted({a[1] for a in A} | {a[3] for a in A if a[0] == 'C'})
    pos = [3 * i + c for i in lb for c in (0, 1)]
    # flat position directions: null space of L on the position columns of load-bearing squares
    Lp = L[:, pos]
    U, sv, Vt = np.linalg.svd(Lp)
    rk = int((sv > 1e-9 * sv[0]).sum()); F = Vt[rk:].T
    mm0, _ = maxmin(L, N, n)
    out = []
    for k in range(F.shape[1]):
        for tau in (1e-3, -1e-3, 1e-2):
            Xm, Ym = list(X), list(Y)
            for idx, col in enumerate(pos):
                i, c = divmod(col, 3)
                if c == 0: Xm[i] += tau * F[idx, k]
                else: Ym[i] += tau * F[idx, k]
            Lm, _ = rows_and_L(n, Xm, Ym, C, Sn, S, A)
            mmv, res = maxmin(Lm, N, n)
            out.append(mmv)
    ok = [x for x in out if x is not None and x > 0]
    print(f'n={n}: rows {L.shape[0]}, lone side-end incidences {bad}, flat position dims {F.shape[1]}, max-min lambda at z* {mm0}; '
          f'along flat dirs: {len(ok)}/{len(out)} feasible with lambda>0, min {min(ok) if ok else None}')
