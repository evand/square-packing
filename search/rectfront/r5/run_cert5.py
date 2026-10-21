#!/usr/bin/env python3
"""Enumerate all Final5 leaves, certify each (2026-10-21).  Output JSON lines to runs/rectfront/r5cert/leaves.jsonl:
case, Hint, region, D, sign, status (empty / cert / FAIL), certificate labels, the polynomial P (as a string)."""
import itertools, json, math, os, sys
import numpy as np
import sympy as sp
from concurrent.futures import ProcessPoolExecutor
from cert5 import parts, leaf_P, gens_for, certify, u, v, H, r, H0

OUT = '/home/evand/math/square-packing/public/runs/rectfront/r5cert'

def feasible(region, Hint, n=4000):
    """region (polynomials >= 0 in u, v) nonempty with interior?  Returns a sample point or None."""
    fs = [sp.lambdify((u, v), g) for g in region]
    rng = np.random.default_rng(1)
    for _ in range(n):
        a, b = rng.uniform(1e-3, 1 - 1e-3, 2)
        if all(f(a, b) > 1e-9 for f in fs): return (a, b)
    return None

def one(args):
    case, Hint = args
    try:
        N, D, region = leaf_P(case, Hint)
    except Exception as ex:
        return dict(case=case, Hint=Hint, status='ERR ' + repr(ex))
    pt = feasible(region, Hint)
    rec = dict(case=case, Hint=Hint, region=[str(sp.factor(g)) for g in region], D=str(D))
    if pt is None:
        rec['status'] = 'empty'
        return rec
    Dval = sp.lambdify((u, v, r, H), sp.expand(D))(pt[0], pt[1], math.sqrt(2), 1.5)
    P = N if Dval > 0 else -N
    rec['Dsign'] = 1 if Dval > 0 else -1
    rec['P'] = str(sp.expand(P))
    G = gens_for(region, Hint)
    for deg in (6, 8, 10):
        out = certify(P, G, deg)
        if out:
            used, resid = out
            rec.update(status='cert', deg=deg, resid=resid, used=[[list(l[0]), l[1], l[2], float(w)] for l, w in used],
                       gens=[[g[0], str(g[1]), g[2]] for g in G])
            return rec
    rec['status'] = 'FAIL'
    return rec

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for e1, e2, b1, b2 in itertools.product((False, True), (False, True), 'UP', 'UP'):
        for cand in ('up', 'a1', 'a2'):
            for mp in (True, False):
                for kp in ((True, False) if cand != 'up' else (True,)):
                    for Hint in ('lo', 'hi'):
                        jobs.append(((e1, e2, b1, b2, cand, mp, kp), Hint))
    with ProcessPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 12) as ex, open(f'{OUT}/leaves.jsonl', 'w') as f:
        for rec in ex.map(one, jobs):
            f.write(json.dumps(rec) + '\n'); f.flush()
            print(rec['case'], rec['Hint'], rec['status'], rec.get('deg', ''), flush=True)
