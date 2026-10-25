#!/usr/bin/env python3
"""Certificates for the band leaves (2026-10-25).

Each sub-leaf: 0 <= P(C, S, H, r) under generators g >= 0 (frame: C, S, C + S - 1, 1 - C, 1 - S; band: H - H9, H87 - H;
r: r - rlo, rhi - r; leaf hypotheses) and C^2 + S^2 = 1.  LP (HiGHS): P = sum lam_k prod(g) + sum mu_m m (C^2 + S^2 - 1),
lam >= 0, products of total degree <= deg.  Lean: `linarith` with the used products as hints plus the needed
`m * (C^2 + S^2) = m` facts recomputes exact coefficients (float LP only proposes the support).

Run: python3 band_cert.py [deg ...]  -> cache runs/rectfront/band/certs.pkl (key: sub-leaf id)."""
import math, sys, os, pickle, itertools, time
import numpy as np
import sympy as sp
from scipy.optimize import linprog
import band_gen as G

C, S, H, r = G.C, G.S, G.H, G.r
VARS = (C, S, H, r)
RLO = sp.Rational(141421356, 100000000); RHI = sp.Rational(141421357, 100000000)
CACHE = os.path.join(os.path.dirname(__file__), '../../../runs/rectfront/band/certs.pkl')

def base_gens():
    return [('C', C), ('S', S), ('p1', C + S - 1), ('C1', 1 - C), ('S1', 1 - S),
            ('H9', H - G.H9), ('H87', G.H87 - H), ('rlo', r - RLO), ('rhi', RHI - r)]

def todict(e):
    return {m: float(c) for m, c in sp.Poly(sp.expand(e), *VARS).terms()}

def dmul(x, y):
    out = {}
    for ma, ca in x.items():
        for mb, cb in y.items():
            m = tuple(i + j for i, j in zip(ma, mb)); out[m] = out.get(m, 0.0) + ca * cb
    return out

def gdeg(e):
    return max(1, sp.Poly(sp.expand(e), *VARS).total_degree())

def certify(P, gens, deg, eqdeg=None):
    """returns (used products [(tuple of gen indices), weight], used eq monomials) or None"""
    eqdeg = deg if eqdeg is None else eqdeg
    gd = [todict(e) for _, e in gens]; gdg = [gdeg(e) for _, e in gens]
    cols, labels = [], []
    # multisets of generators with total degree <= deg (r generators at most once)
    def rec(start, budget, cur, poly):
        cols.append(poly); labels.append(tuple(cur))
        for k in range(start, len(gens)):
            if gdg[k] <= budget:
                if gens[k][0] in ('rlo', 'rhi') and any(gens[c][0] in ('rlo', 'rhi') for c in cur): continue
                rec(k, budget - gdg[k], cur + [k], dmul(poly, gd[k]))
    rec(0, deg, [], {(0, 0, 0, 0): 1.0})
    neq = len(cols)
    eqm = []
    circle = todict(C**2 + S**2 - 1)
    for i in range(eqdeg - 1):
        for j in range(eqdeg - 1 - i):
            for k in range(3):
                for l in range(2):
                    m = (i, j, k, l)
                    cols.append(dmul({m: 1.0}, circle)); eqm.append(m)
    target = todict(P)
    monos = sorted({m for c in cols for m in c} | set(target)); idx = {m: n for n, m in enumerate(monos)}
    A = np.zeros((len(monos), len(cols)))
    for jj, c in enumerate(cols):
        for m, v in c.items(): A[idx[m], jj] = v
    bvec = np.zeros(len(monos))
    for m, v in target.items(): bvec[idx[m]] = v
    bounds = [(0, None)] * neq + [(None, None)] * len(eqm)
    cost = np.concatenate([np.ones(neq), np.zeros(len(eqm))])
    res = linprog(cost, A_eq=A, b_eq=bvec, bounds=bounds, method='highs',
                  options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
    if res.status != 0: return None
    used = [(labels[jj], res.x[jj]) for jj in range(neq) if res.x[jj] > 1e-11]
    ueq = [eqm[jj - neq] for jj in range(neq, len(cols)) if abs(res.x[jj]) > 1e-11]
    return used, ueq, float(np.abs(A @ res.x - bvec).max())

def leaf_gens(hyps):
    gens = base_gens()
    for n, (k, h) in enumerate(hyps):
        gens.append((f'h{n}', sp.expand(h)))
    return gens

def leaf_id(lf, ch):
    return f"{lf.cond}_{lf.cell[0]}{lf.cell[1]}_{lf.tag}_{''.join(map(str, ch))}"

TIGHT = {'Sa_10_xk_yk_', 'Sb_01_xk_yk_', 'C0_11_xk_yk_', 'C1_00_xk_yk_'}

if __name__ == '__main__':
    degs = [int(d) for d in sys.argv[1:]] or [2, 3, 4, 5]
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    todo = []
    for lf in G.all_leaves():
        for ch, e, hy in G.subcases(lf):
            lid = leaf_id(lf, ch)
            if lid in TIGHT or lid in cache: continue
            todo.append((lid, e, hy))
    print(len(todo), 'to certify', flush=True)
    fails = []
    for lid, e, hy in todo:
        t0 = time.time(); got = None
        for d in degs:
            gens = leaf_gens(hy)
            res = certify(e, gens, d)
            if res and res[2] < 1e-7:
                got = (d, [(tuple(gens[k][0] for k in lab), w) for lab, w in res[0]], res[1], [gname for gname, _ in gens])
                break
        if got:
            cache[lid] = got
            print(f"  {lid}: deg {got[0]}, {len(got[1])} products, {len(got[2])} eq ({time.time() - t0:.1f}s)", flush=True)
        else:
            fails.append(lid); print(f"  {lid}: FAIL ({time.time() - t0:.1f}s)", flush=True)
        pickle.dump(cache, open(CACHE, 'wb'))
    print("failed:", fails)


# ---------------------------------------------------------------- subdivision for near-tight leaves (36.87 deg)
SPLITS = [[], [('Cs', C - sp.Rational(4, 5))], [('Hs', H - sp.Rational(1858, 1000))],
          [('Cs', C - sp.Rational(4, 5)), ('Hs', H - sp.Rational(1858, 1000))]]

def certify_split(P, hyps, degs):
    """try split sets in order; for a split set, every sign combination must certify.  Returns (splits, plans)."""
    for sp_set in SPLITS[1:]:
        plans = []
        for signs in itertools.product((1, -1), repeat=len(sp_set)):
            extra = [(nm + ('+' if sg > 0 else '-'), sp.expand(sg * e)) for (nm, e), sg in zip(sp_set, signs)]
            gens = leaf_gens(hyps) + extra
            got = None
            for d in degs:
                res = certify(P, gens, d)
                if res and res[2] < 1e-7:
                    got = (d, [(tuple(gens[k][0] for k in lab), w) for lab, w in res[0]], res[1], [g for g, _ in gens], extra)
                    break
            if not got: break
            plans.append(got)
        if len(plans) == 2 ** len(sp_set):
            return [nm for nm, _ in sp_set], plans
    return None

def run_splits(degs=(4, 5, 6)):
    cache = pickle.load(open(CACHE, 'rb'))
    for lf in G.all_leaves():
        for ch, e, hy in G.subcases(lf):
            lid = leaf_id(lf, ch)
            if lid in TIGHT or lid in cache: continue
            t0 = time.time()
            got = certify_split(e, hy, degs)
            print(f"  {lid}: {'split ' + str(got[0]) if got else 'FAIL'} ({time.time() - t0:.1f}s)", flush=True)
            if got: cache[lid] = ('split', got); pickle.dump(cache, open(CACHE, 'wb'))
