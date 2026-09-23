#!/usr/bin/env python3
"""s6exact.py -- closed forms for the fixed-assignment margin along a one-parameter tilt (2026-09-20).

For a FIXED separating-axis assignment the margin is an LP whose data are trigonometric in the
angles.  Any dual-feasible w >= 0 gives, for that assignment and for every t where w(t) >= 0,

        delta  <=  - sum_r w_r(t) b_r(t)                                  (exact, by weak duality)

where row r reads  a_r(t).x + d_r*delta >= b_r(t)  (d_r = -1 on pair rows, 0 on containment rows),
and dual feasibility is  sum_r w_r a_r = 0,  sum_r w_r d_r = -1.  This script

  1. finds the best configuration at one numerical t (s6local), its assignment and the support
     of the HiGHS dual;
  2. solves the dual on that support SYMBOLICALLY in t (sympy), and
  3. prints the closed form, its series at t = 0, and a numerical cross-check.

    python3 search/s6exact.py --n 6  --T 3 --pattern uniform --tdeg 1.6
    python3 search/s6exact.py --n 6  --T 3 --pattern 0,0,0,1,1,1 --tdeg 1.6
    python3 search/s6exact.py --n 12 --T 4 --pattern uniform --tdeg 1.6
"""
import argparse
import math
import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6local                                                       # noqa: E402


def rows_numeric(F, asg):
    """all rows (a, d, b, tag): a.x + d*delta >= b.  Containment is hard (d = 0), as in
    s6local.Fixed (its wall rows are implied by the bounds once delta <= 0)."""
    n = F.n
    rows = []
    for i in range(n):
        for ax, off in (('x', 0), ('y', n)):
            a = np.zeros(2 * n); a[off + i] = 1.0
            rows.append((a, 0.0, F.P[i], ('lo', ax, i)))
            a = np.zeros(2 * n); a[off + i] = -1.0
            rows.append((a, 0.0, F.P[i] - F.T, ('hi', ax, i)))
    for (i, j), (o, k, sg) in zip(F.pairs, asg):
        cc, ss = F.cs[o]
        d0, d1 = ((cc, ss), (-ss, cc))[k]
        a = np.zeros(2 * n)
        a[j] = sg * d0; a[i] = -sg * d0; a[n + j] = sg * d1; a[n + i] = -sg * d1
        rows.append((a, -1.0, F.m[(i, j)], ('pair', i, j, o, k, sg)))
    return rows


def solve_dual_numeric(rows, n):
    """min sum w b  s.t.  sum w a = 0, sum w d = -1, w >= 0   (value = -delta)."""
    A = np.array([np.append(a, d) for (a, d, _b, _t) in rows]).T
    rhs = np.zeros(2 * n + 1); rhs[-1] = -1.0
    b = np.array([r[2] for r in rows])
    res = linprog(-b, A_eq=A, b_eq=rhs, bounds=[(0, None)] * len(rows), method='highs',
                  options={'primal_feasibility_tolerance': 1e-10,
                           'dual_feasibility_tolerance': 1e-10})
    assert res.success, res.message
    return res.x, res.fun          # delta = res.fun  (since we minimise -sum w b ... see below)


def sym_row(tag, mult, base, T, t, n):
    """the same row with symbolic t."""
    th = [sp.nsimplify(b0) + sp.nsimplify(m) * t for b0, m in zip(base, mult)]
    if tag[0] in ('lo', 'hi'):
        _k, ax, i = tag
        P = (sp.cos(th[i]) + sp.sin(th[i])) / 2          # tilts in [0, 90 deg)
        a = [0] * (2 * n)
        off = 0 if ax == 'x' else n
        if tag[0] == 'lo':
            a[off + i] = 1
            return a, 0, P
        a[off + i] = -1
        return a, 0, P - sp.nsimplify(T)
    _k, i, j, o, k, sg = tag
    cc, ss = sp.cos(th[o]), sp.sin(th[o])
    d0, d1 = ((cc, ss), (-ss, cc))[k]
    a = [0] * (2 * n)
    a[j] = sg * d0; a[i] = -sg * d0; a[n + j] = sg * d1; a[n + i] = -sg * d1
    D = th[j] - th[i]
    # |cos D| + |sin D| with the sign of sin D fixed by the ordering of the tilts at the sample
    return a, -1, sp.Rational(1, 2) + (sp.cos(D) + sp.sin(sp.Abs(D))) / 2   # |D| < 90 deg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--T', type=float, required=True)
    ap.add_argument('--pattern', default='uniform')
    ap.add_argument('--tdeg', type=float, default=1.6)
    ap.add_argument('--limit', type=int, default=2000)
    ap.add_argument('--jitters', type=int, default=4)
    ap.add_argument('--nproc', type=int, default=14)
    ap.add_argument('--order', type=int, default=5)
    a = ap.parse_args()

    n, T = a.n, a.T
    mult = s6local.pattern_vec(a.pattern, n)
    base = [0.0] * n
    t0 = math.radians(a.tdeg)
    theta = [m * t0 for m in mult]
    rng = np.random.default_rng(0)
    starts = s6local.structured_starts(n, T, mult, rng, a.limit, a.jitters, 0.08)
    top, _all = s6local.best_at(theta, T, starts, a.nproc, keep=1)
    v, X, Y = top[0]
    F = s6local.Fixed(theta, T)
    v2, X, Y, asg = F.ascent(X, Y)
    print(f"# n={n} T={T} pattern={a.pattern}  sample t={a.tdeg} deg: geometric value {v2:+.12e}")

    rows = rows_numeric(F, asg)
    w, obj = solve_dual_numeric(rows, n)
    delta_num = obj            # min(-b.w) = -(max b.w);  delta <= -sum w b, tight at optimum
    print(f"# fixed-assignment LP value (dual)  delta = {delta_num:+.12e}")
    supp = [r for r in range(len(rows)) if w[r] > 1e-9]
    print(f"# dual support: {len(supp)} rows of {len(rows)}")
    for r in supp:
        print(f"    w={w[r]:.6f}  {rows[r][3]}")

    t = sp.symbols('t', positive=True)
    srows = [sym_row(rows[r][3], mult, base, T, t, n) for r in supp]
    ws = sp.symbols(f'w0:{len(supp)}')
    eqs = []
    for col in range(2 * n):
        e = sum(ws[k] * srows[k][0][col] for k in range(len(supp)))
        if e != 0:
            eqs.append(e)
    eqs.append(sum(ws[k] * srows[k][1] for k in range(len(supp))) + 1)
    sol = sp.solve(eqs, ws, dict=True)
    if not sol:
        print("# symbolic dual: no solution on this support"); return
    sol = sol[0]
    free = [x for x in ws if x not in sol]
    if free:
        print(f"# symbolic dual has {len(free)} free multipliers; pinning them to the numeric values")
        sub = {x: sp.nsimplify(w[supp[ws.index(x)]], rational=True, tolerance=1e-9) for x in free}
        sol = {k: sp.simplify(val.subs(sub)) for k, val in sol.items()}
        sol.update(sub)
    bound = -sum(sol[ws[k]] * srows[k][2] for k in range(len(supp)))
    bound = sp.simplify(sp.trigsimp(bound))
    print("# closed form   delta(t) <=", bound)
    print("# series        ", sp.series(bound, t, 0, a.order))
    print(f"# check at sample: closed form = {float(bound.subs(t, t0)):+.12e}")
    print("# multipliers (series):")
    for k in range(len(supp)):
        print(f"    {rows[supp[k]][3]}:  {sp.series(sp.simplify(sol[ws[k]]), t, 0, 3)}")


if __name__ == '__main__':
    main()
