#!/usr/bin/env python3
"""sos_probe_sdp.py -- Positivstellensatz infeasibility certificates for one separation type (2026-10-03).

Task `tasks/sos-probe/README.md`, write-up `search/SOS_PROBE.md`, encoding `search/sos_probe_enc.py`.

At a fixed container side T = t (rational) we look for

    -lam  =  sigma_0  +  sum_g sigma_g g  +  sum_{a,b} tau_ab g_a g_b          (mod  c_i^2 + s_i^2 - 1)

with sigma SOS, tau SOS (or nonnegative scalars), lam > 0: an algebraic proof that the type has no
configuration in [0, t]^2.  (Feasible sets are monotone in T, so this proves T > t on the type.)
The equality constraints are handled exactly by reducing every monomial to normal form
(s_i^2 -> 1 - c_i^2), so there are no free multipliers.

--mult angle  (default): every multiplier is a polynomial in the ANGLE variables (c_i, s_i) only.
   This is the natural form: for fixed angles the problem is an LP in the centres, and an
   angle-multiplier certificate is an LP dual w(theta) >= 0 that is polynomial/SOS in theta.
--mult full: multipliers in all variables (Lasserre/Putinar proper; much larger).
--order d: total degree budget 2d, counted in the angle variables (angle mode) or all variables.
--prod: products g_a g_b with g_a angle-only (and g_b any generator), with SOS multiplier of the
   remaining degree (--pdeg) or scalar (--pdeg 0).

Normalisation: sum of Gram traces + scalars <= 1; then lam is a margin (scale-free comparisons
across t are meaningful only within one program size).
"""
import argparse
import itertools
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sos_probe_enc as E                                               # noqa: E402


def subsystem(cells, typ, idx):
    m = {o: k for k, o in enumerate(idx)}
    t2 = {}
    for (i, j), (o, k, sg) in typ.items():
        if i in m and j in m:
            t2[(m[i], m[j])] = (m[o], k, sg)
    return [cells[i] for i in idx], t2


class Prog:
    def __init__(self, nv, n, mult='angle'):
        self.nv, self.n, self.mult = nv, n, mult
        self.cvars = [3 + 4 * i for i in range(n)]
        self.svars = [4 + 4 * i for i in range(n)]
        self.avars = sorted(self.cvars + self.svars)
        self.mono_index = {}
        self.blocks = []         # (name, basis, g)  -> PSD Gram
        self.scal = []           # (name, g)         -> scalar >= 0
        self._red = {}

    def basis(self, d, squares=None):
        """standard monomials (s-degree <= 1 per square) of degree <= d in the multiplier variables
        (restricted to the angle variables of `squares` when given: correlative sparsity)."""
        vs = self.avars if self.mult == 'angle' else list(range(self.nv))
        if squares is not None:
            vs = [v for v in vs if v > 0 and (v - 1) // 4 in squares]
        if getattr(self, 'active', None) is not None:
            vs = [v for v in vs if v in self.active]
        out = []
        for deg in range(d + 1):
            for comb in itertools.combinations_with_replacement(vs, deg):
                e = [0] * self.nv
                for i in comb:
                    e[i] += 1
                if all(e[s] <= 1 for s in self.svars):
                    out.append(tuple(e))
        return out

    def adeg(self, g):
        """degree of g counted in the multiplier-relevant variables."""
        if self.mult == 'angle':
            return max(sum(e[v] for v in self.avars) for e in g.t)
        return g.deg()

    def reduce(self, e):
        """normal form of x^e modulo s_i^2 = 1 - c_i^2: dict e' -> Fr."""
        r = self._red.get(e)
        if r is not None:
            return r
        res = {e: Fr(1)}
        for c, s in zip(self.cvars, self.svars):
            new = {}
            for f, cf in res.items():
                k = f[s]
                if k <= 1:
                    new[f] = new.get(f, 0) + cf
                    continue
                q, rmd = divmod(k, 2)
                for j in range(q + 1):             # (1 - c^2)^q = sum C(q,j) (-1)^j c^{2j}
                    g = list(f); g[s] = rmd; g[c] += 2 * j
                    g = tuple(g)
                    new[g] = new.get(g, 0) + cf * math.comb(q, j) * (-1) ** j
            res = {k: v for k, v in new.items() if v != 0}
        self._red[e] = res
        return res

    def idx(self, e):
        k = self.mono_index.get(e)
        if k is None:
            k = self.mono_index[e] = len(self.mono_index)
        return k

    def _add_terms(self, rows, cols, vals, col, e, cf, exact):
        for f, c2 in self.reduce(e).items():
            rows.append(self.idx(f)); cols.append(col)
            vals.append(cf * c2 if exact else float(cf * c2))

    def build(self, exact=False):
        mats = []
        for name, B, g in self.blocks:
            k = len(B)
            rows, cols, vals = [], [], []
            polyb = k and isinstance(B[0], E.Poly)          # basis of polynomials (facial reduction)
            for a in range(k):
                for b in range(k):
                    if polyb:
                        for e, cf in (B[a] * B[b] * g).t.items():
                            self._add_terms(rows, cols, vals, a + b * k, e, cf, exact)
                        continue
                    ab = tuple(x + y for x, y in zip(B[a], B[b]))
                    for e, cf in g.t.items():
                        self._add_terms(rows, cols, vals, a + b * k, tuple(x + y for x, y in zip(ab, e)), cf, exact)
            mats.append(('psd', k, rows, cols, vals, name))
        rows, cols, vals = [], [], []
        for col, (name, g) in enumerate(self.scal):
            for e, cf in g.t.items():
                self._add_terms(rows, cols, vals, col, e, cf, exact)
        mats.append(('scal', len(self.scal), rows, cols, vals, 'scalars'))
        if exact:
            return mats
        M = len(self.mono_index)
        out = []
        for kind, k, r, c, v, name in mats:
            ncol = k * k if kind == 'psd' else k
            out.append((kind, k, sps.csr_matrix((v, (r, c)), shape=(M, ncol)), name))
        return out


def squares_of(g):
    return frozenset((v - 1) // 4 for v in g.vars() if v > 0)


def make_program(n, G, names, order, mult='angle', prod=False, pdeg=0, csp=0, csp0=0):
    """csp = 0: dense multipliers.  csp = k >= 2: each multiplier only in the angles of the squares its
    generator touches, and sigma_0 a sum of SOS over all k-subsets of squares.  csp0 = k: only sigma_0
    is split over k-subsets (multipliers stay dense in all angles)."""
    nv = 1 + 4 * n
    P = Prog(nv, n, mult)
    P.active = set(v for g in G for v in g.vars())
    one = E.Poly.const(nv, 1)
    P.idx((0,) * nv)
    if csp or csp0:
        for sq in itertools.combinations(range(n), csp or csp0):
            P.blocks.append((f'sigma0{sq}', P.basis(order, set(sq)), one))
    else:
        P.blocks.append(('sigma0', P.basis(order), one))
    for g, nm in zip(G, names):
        dm = (2 * order - P.adeg(g)) // 2
        if dm >= 0:
            P.blocks.append((str(nm), P.basis(dm, squares_of(g) if csp else None), g))
    if prod:
        cfree = [all((v - 1) % 4 >= 2 for v in g.vars() if v > 0) and 0 not in g.vars() for g in G]
        for a in range(len(G)):
            if not cfree[a]:
                continue
            for b in range(len(G)):
                if cfree[b] and b <= a:
                    continue
                gg = G[a] * G[b]
                if not gg.t:
                    continue
                room = 2 * order - P.adeg(gg)
                if room < 0:
                    continue
                dm = min(pdeg, room // 2)
                if dm == 0:
                    P.scal.append((f'{names[a]}*{names[b]}', gg))
                else:
                    P.blocks.append((f'{names[a]}*{names[b]}', P.basis(dm, squares_of(gg) if csp else None), gg))
    return P


def solve_infeas(P, solver='CLARABEL', verbose=False, lamfix=None, bound=False):
    """max lam s.t. sum(...) + lam = 0 identically, traces + scalars <= 1.
    With lamfix: lam fixed, maximise the common margin mu (Gram - mu I >= 0, scalars >= mu)."""
    import cvxpy as cp
    t0 = time.time()
    mats = P.build()
    M = len(P.mono_index)
    e0 = np.zeros(M); e0[P.idx((0,) * P.nv)] = 1.0
    eT = np.zeros(M)
    tkey = tuple(1 if i == 0 else 0 for i in range(P.nv))
    if tkey in P.mono_index:
        eT[P.mono_index[tkey]] = 1.0
    # drop linearly dependent monomial rows (Clarabel's KKT fails on them); consistency of the
    # right-hand side is checked after the solve on ALL rows
    Afull = sps.hstack([m[2] for m in mats]).tocsr()
    rows = np.arange(M)
    if M < 25000:
        import scipy.linalg as sla
        AAT = (Afull @ Afull.T).toarray()          # column j of AAT <-> row j of A; same rank pattern
        _q, rr, piv = sla.qr(AAT, mode='economic', pivoting=True)
        dg = np.abs(np.diag(rr))
        rank = int((dg > 1e-12 * dg[0]).sum())
        rows = np.sort(piv[:rank])
    full_mats = mats
    mats = [(kind, k, A[rows], name) for kind, k, A, name in mats]
    e0f, eTf = e0, eT
    e0, eT = e0[rows], eT[rows]
    M = len(rows)
    lam = cp.Variable()
    mu = cp.Variable()
    expr = 0
    tr = 0
    var = []
    for kind, k, A, name in mats:
        if kind == 'psd':
            R = cp.Variable((k, k), PSD=True)
            expr = expr + A @ cp.vec(R, order='F')
            tr = tr + cp.trace(R)
            if lamfix is not None:
                expr = expr + mu * (A @ np.eye(k).flatten(order='F'))
                tr = tr + mu * k
        else:
            if k == 0:
                var.append((kind, name, None)); continue
            R = cp.Variable(k, nonneg=True)
            expr = expr + A @ R
            tr = tr + cp.sum(R)
            if lamfix is not None:
                expr = expr + mu * np.asarray(A.sum(axis=1)).ravel()
                tr = tr + mu * k
        var.append((kind, name, R))
    if bound and lamfix is not None:
        # gamma fixed (e.g. exactly 3 after facial reduction): maximise the interior margin mu
        prob = cp.Problem(cp.Maximize(mu), [expr + lamfix * e0 == eT, mu <= 1])
    elif bound:
        # T - gamma = sum(...):  lam plays gamma; no normalisation needed
        prob = cp.Problem(cp.Maximize(lam), [expr + lam * e0 == eT])
    elif lamfix is None and os.environ.get('SOS_NORM', 'trace') == 'lam1':
        # equivalent scaling: lam = 1, minimise the trace; reported lam = 1 / trace
        prob = cp.Problem(cp.Minimize(tr), [expr + e0 == 0])
    elif lamfix is None:
        prob = cp.Problem(cp.Maximize(lam), [expr + lam * e0 == 0, tr <= 1])
    else:
        prob = cp.Problem(cp.Maximize(mu), [expr + lamfix * e0 == 0, tr <= 1])
    t1 = time.time()
    kw = {}
    if solver == 'CLARABEL':
        kw = dict(tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=300,
                  direct_solve_method=os.environ.get('SOS_KKT', 'faer'),
                  equilibrate_enable=os.environ.get('SOS_EQ', '1') == '1')
    elif solver == 'SCS':
        kw = dict(eps_abs=1e-9, eps_rel=1e-9, max_iters=200000, acceleration_lookback=20)
    try:
        prob.solve(solver=solver, verbose=verbose, **kw)
    except cp.error.SolverError:
        pass
    t2 = time.time()
    vals = []
    for kind, name, R in var:
        if R is None:
            vals.append(np.zeros(0)); continue
        v = R.value
        if v is None:
            vals.append(None); continue
        if lamfix is not None and v is not None:
            v = v + (mu.value * np.eye(v.shape[0]) if kind == 'psd' else mu.value)
        vals.append(v)
    lamv = lamfix
    if lamfix is None and lam.value is not None:
        lamv = float(lam.value)
    elif lamfix is None and os.environ.get('SOS_NORM') == 'lam1' and prob.value is not None and not bound:
        lamv = 1.0 / float(prob.value) if prob.status.startswith('optimal') else 0.0
    resid = None
    if all(v is not None for v in vals) and lamv is not None:
        tot = np.zeros(len(e0f))
        for (kind, k, A, name), v in zip(full_mats, vals):
            if kind == 'psd':
                tot += A @ np.asarray(v).flatten(order='F')
            elif k:
                tot += A @ np.asarray(v)
        tot += (lamv if lamfix is None and os.environ.get('SOS_NORM') != 'lam1' else (1.0 if lamfix is None else lamfix)) * e0f
        if bound:
            tot -= eTf
        resid = float(np.abs(tot).max())
    # support: blocks (generators) carrying weight, by trace
    supp = []
    if all(v is not None for v in vals):
        wts = []
        for (kind, k, A, name), v in zip(full_mats, vals):
            if kind == 'psd':
                wts.append((float(np.trace(v)), name))
            elif k:
                for (nm, _g), x in zip(P.scal, np.asarray(v).ravel()):
                    wts.append((float(x), nm))
        mx = max(w for w, _ in wts) if wts else 0
        supp = [(round(w, 8), nm) for w, nm in sorted(wts, reverse=True) if w > 1e-5 * mx]
    info = dict(status=prob.status, lam=lamv, resid=resid, rank_rows=int(M), support=supp[:400],
                mu=(float(mu.value) if lamfix is not None and mu.value is not None else None),
                n_mono=M, gram=sorted({k for kind, k, A, n in mats if kind == 'psd'}, reverse=True),
                n_gram=sum(1 for kind, k, A, n in mats if kind == 'psd'),
                n_scal=len(P.scal), t_build=round(t1 - t0, 2), t_solve=round(t2 - t1, 2))
    return info, vals


def setup(a):
    cells = getattr(E, a.cfg)
    over = None
    if a.override:
        over = {}
        for item in a.override.split('/'):
            ij, t = item.split('=')
            i, j = map(int, ij.split(','))
            over[(i, j)] = None if t == 'none' else tuple(int(x) for x in t.split(','))
    typ = E.make_type(cells, diag=a.diag, owner=a.owner, far=not a.nofar, overrides=over)
    idx = [int(x) for x in a.sub.split(',')] if a.sub else list(range(len(cells)))
    cells, typ = subsystem(cells, typ, idx)
    n = len(cells)
    Tfix = None if a.bound else Fr(a.T).limit_denominator(10 ** 6)
    assert not a.orthant or len(a.orthant) == n, 'orthant string must have one sign per square'
    nv, G, H, names = E.system(n, typ, Tfix=Tfix, extra=not a.noextra, orth=a.orthant or None)
    if a.maxtilt is not None:
        s2 = Fr(math.sin(math.radians(a.maxtilt)) ** 2).limit_denominator(10 ** 6)
        for i in range(n):
            S = E.Poly.var(nv, 4 + 4 * i)
            G.append(s2 - S * S); names.append(('tilt', i))
    if a.bound:     # prove T >= gamma on {T <= 3}: the (3 - T) generator lets the T-weight vary with angles
        G.append(3 - E.Poly.var(nv, 0)); names.append(('T<=3',))
    if a.classes:   # angle classes: '0' = exactly axis-parallel, equal labels share one angle
        lab = a.classes.split(',')
        assert len(lab) == n
        rep, m = {}, {}
        for i, l in enumerate(lab):
            if l == '0':
                m[3 + 4 * i] = E.Poly.const(nv, 1); m[4 + 4 * i] = E.Poly.const(nv, 0)
            elif l in rep:
                m[3 + 4 * i] = E.Poly.var(nv, 3 + 4 * rep[l]); m[4 + 4 * i] = E.Poly.var(nv, 4 + 4 * rep[l])
            else:
                rep[l] = i
        G2, N2 = [], []
        for g, nm in zip(G, names):
            g2 = g.subs(m)
            if any(sum(k) > 0 for k in g2.t):               # drop generators that became constants
                G2.append(g2); N2.append(nm)
            else:
                assert all(v >= 0 for v in g2.t.values()), ('class substitution made a row infeasible', nm)
        G, names = G2, N2
    for i in [int(x) for x in a.axis.split(',') if x]:    # square i exactly axis-parallel
        S = E.Poly.var(nv, 4 + 4 * i)
        G.append(-S * S); names.append(('axis', i))
    return n, typ, G, names, idx


def parser():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cfg', default='ROW_A')
    ap.add_argument('--sub', default='')
    ap.add_argument('--diag', default='x')
    ap.add_argument('--owner', default='low')
    ap.add_argument('--override', default='', help="'i,j=o,k,sg/i,j=none/...' (original cell indices)")
    ap.add_argument('--nofar', action='store_true')
    ap.add_argument('--order', type=int, default=2)
    ap.add_argument('--T', type=float, default=2.9)
    ap.add_argument('--mult', default='angle')
    ap.add_argument('--prod', action='store_true')
    ap.add_argument('--pdeg', type=int, default=0)
    ap.add_argument('--csp', type=int, default=0)
    ap.add_argument('--csp0', type=int, default=0)
    ap.add_argument('--noextra', action='store_true')
    ap.add_argument('--orthant', default='', help="sign of s_i per square, e.g. '++-+0+' ('0' = free)")
    ap.add_argument('--maxtilt', type=float, default=None)
    ap.add_argument('--axis', default='', help='squares forced axis-parallel (s_i = 0)')
    ap.add_argument('--classes', default='', help="angle classes, e.g. '0,0,0,a,a' (0 = axis-parallel)")
    ap.add_argument('--bound', action='store_true', help='maximise gamma with T - gamma in the module (T free, T <= 3)')
    ap.add_argument('--solver', default='CLARABEL')
    ap.add_argument('--lamfix', type=float, default=None)
    ap.add_argument('--save', default='')
    ap.add_argument('--verbose', action='store_true')
    return ap


def main():
    a = parser().parse_args()
    n, typ, G, names, idx = setup(a)
    P = make_program(n, G, names, a.order, a.mult, a.prod, a.pdeg, a.csp, a.csp0)
    info, vals = solve_infeas(P, a.solver, a.verbose, a.lamfix, a.bound)
    info.update(cfg=a.cfg, sub=idx, diag=a.diag, owner=a.owner, override=a.override, far=not a.nofar,
                order=a.order, T=a.T, mult=a.mult, prod=a.prod, pdeg=a.pdeg, csp=a.csp, csp0=a.csp0, orthant=a.orthant,
                maxtilt=a.maxtilt, axis=a.axis, classes=a.classes, bound=a.bound, n_gen=len(G), n_pairs=len(typ))
    print(json.dumps(info), flush=True)
    if a.save and info['status'] in ('optimal', 'optimal_inaccurate'):
        np.savez_compressed(a.save, *vals)


if __name__ == '__main__':
    main()
