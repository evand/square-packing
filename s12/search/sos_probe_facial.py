#!/usr/bin/env python3
"""sos_probe_facial.py -- the EXACT bound T >= 3 (no epsilon) on a family, by facial reduction (2026-10-03).

    T - 3 = sigma_0 + sum_g sigma_g g + sigma_q (3 - T) + ...     (angle multipliers, mod the ideal)

is degenerate: at the zero-margin point (all tilts 0, T = 3) every term vanishes, so each Gram matrix
whose generator is not identically 0 on the plateau must kill the monomial vector v0 = v(theta = 0)
(entries 1 on s-free monomials, 0 otherwise).  The SDP at gamma = 3 then has no interior point and
float solutions cannot be rounded.  Facial reduction step used here:

  1. solve the bound SDP (sos_probe_sdp --bound) in floats;
  2. drop blocks with negligible trace; for every remaining block with v0' Q v0 ~ 0 replace its
     monomial basis by a rational basis of v0-perp:  {m : m contains some s_i} u {m - 1 : m s-free};
  3. re-solve with gamma = 3 EXACTLY, maximising the interior margin mu;
  4. round + exact projection + exact LDL (sos_probe_exact.exact_round).

    python3 search/sos_probe_facial.py --cfg PINWHEEL --sub 0,1,2,3,5 --orthant=+++++ --prod --bound \
        --classes a,a,a,a,a --order 2
"""
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sos_probe_enc as E                                               # noqa: E402
import sos_probe_sdp as S                                               # noqa: E402
import sos_probe_exact as X                                             # noqa: E402


def main():
    ap = S.parser()
    ap.add_argument('--bits', type=int, default=30)
    ap.add_argument('--droptol', type=float, default=1e-7)
    ap.add_argument('--v0tol', type=float, default=1e-6)
    a = ap.parse_args()
    assert a.bound
    t0 = time.time()
    n, typ, G, names, idx = S.setup(a)
    P = S.make_program(n, G, names, a.order, a.mult, a.prod, a.pdeg, a.csp, a.csp0)
    info1, vals = S.solve_infeas(P, a.solver, bound=True)
    out = dict(cfg=a.cfg, sub=idx, classes=a.classes, orthant=a.orthant, order=a.order,
               gamma_float=info1['lam'], status1=info1['status'], n_blocks=len(P.blocks), n_scal=len(P.scal))
    psd_vals = vals[:-1]
    # plateau points: vertices of the polytope {theta = 0, T = 3, all rows >= 0} (random LP objectives)
    from scipy.optimize import linprog
    th0 = {}
    for c, sv in zip(P.cvars, P.svars):
        th0[c], th0[sv] = 1.0, 0.0
    def at(g, xy):
        z = np.zeros(P.nv); z[0] = 3.0
        for i in range(n):
            z[1 + 4 * i], z[2 + 4 * i] = xy[i], xy[n + i]
            z[3 + 4 * i], z[4 + 4 * i] = 1.0, 0.0
        return g.ev(z)
    rowsA, rowsb = [], []
    for g in G:
        # g is affine in the centres at theta = 0, T = 3: g = a.xy + b
        b = at(g, np.zeros(2 * n))
        a_ = np.array([at(g, np.eye(2 * n)[k]) - b for k in range(2 * n)])
        rowsA.append(-a_); rowsb.append(b)
    rng = np.random.default_rng(3)
    plateau = []
    for _ in range(40):
        r = linprog(rng.normal(size=2 * n), A_ub=np.array(rowsA), b_ub=np.array(rowsb),
                    bounds=[(0, 3)] * (2 * n), method='highs')
        if r.status == 0:
            plateau.append(r.x)
    out['plateau_points'] = len(plateau)
    slack = [max((at(g, xy) for xy in plateau), default=0.0) > 1e-9 for g in G]
    gslack = {str(nm): sl for nm, sl in zip(names, slack)}
    gslack['sigma0'] = True
    gslack["('T<=3',)"] = True     # 3 - T = 0 on the plateau but the block also sees T > 3 points
    new_blocks, n_drop, n_red = [], 0, 0
    for (name, B, g), Q in zip(P.blocks, psd_vals):
        if gslack.get(name, True):
            n_red += 1
            one = E.Poly.const(P.nv, 1)
            nb = []
            for m in B:
                pm = E.Poly(P.nv, {m: Fr(1)})
                if any(m[s] for s in P.svars):
                    nb.append(pm)
                elif any(m):
                    nb.append(pm - one)
            if nb:
                new_blocks.append((name, nb, g))
            else:
                n_drop += 1
        else:
            new_blocks.append((name, [E.Poly(P.nv, {m: Fr(1)}) for m in B], g))
    # a scalar multiplier of g_a g_b is forced to 0 if g_a g_b > 0 somewhere on the plateau
    keep_sc = []
    for nm, g in P.scal:
        if max((at(g, xy) for xy in plateau), default=0.0) <= 1e-9:
            keep_sc.append((nm, g))
    out.update(dropped=n_drop, reduced=n_red, kept=len(new_blocks), scal_kept=len(keep_sc))
    P2 = S.Prog(P.nv, P.n, P.mult)
    P2.active = P.active
    P2.idx((0,) * P.nv)
    P2.blocks = new_blocks
    P2.scal = keep_sc
    for e in [tuple(1 if i == 0 else 0 for i in range(P.nv))]:
        P2.idx(e)
    info2, vals2 = S.solve_infeas(P2, a.solver, lamfix=3.0, bound=True)
    out.update(mu=info2['mu'], status2=info2['status'], resid2=info2.get('resid'))
    out['t_sdp'] = round(time.time() - t0, 1)
    if info2['mu'] is None or info2['mu'] <= 0:
        out['verdict'] = 'facial reduction (v0 step) left no interior: mu <= 0'
        print(json.dumps(out)); return
    X.exact_round(P2, vals2, Fr(3), a.bits, True, out)
    out['verdict'] = ('EXACT CERTIFICATE T >= 3' if out['identity_exact'] and out['blocks_not_psd'] == 0
                      else 'rounding failed')
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
