"""Contact forces of certified minima (f64 multipliers at the exact point) for the explainer figures.
forces(path) -> dict(S, sq=[(x, y, th_rad)], contacts=[(kind, i, j_or_wall, px, py, lam)], free=[...], flat=[...])
Contacts: the load-bearing set from exactsolve's contacts.json; lambda = full_lambda_lp (max-min multipliers)."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), '..', 'exact'))
import numpy as np
from mpmath import mp
import exactsolve as E
from geom import offset


def forces(path):
    base = path.replace('/b', '/exact/b')[:-4]
    mp.dps = 40
    P = E.Pack.load(base + '.exact.txt')
    c = json.load(open(base + '.contacts.json'))
    rep = json.load(open(base + '.json'))
    A = [tuple(x) for x in c['load_bearing']]
    sqs = sorted({x[1] for x in A} | {x[3] for x in A if x[0] == 'C'})
    V = [3 * i + k for i in sqs for k in range(3)] + [3 * P.n]
    t, lam = E.full_lambda_lp(P, A, V)
    if lam is None:                                   # maxmin LP fails on some healthy points: max-support multipliers
        J, _ = E.jacobian(P, A)
        lam, _, _ = E.max_support(J)
    X, Y, C, Sn, S = P.flt()
    out = []
    for ct, l in zip(A, lam):
        if ct[0] == 'W':
            _, i, a, w = ct
            ox, oy = offset(C[i], Sn[i], a, 0.5)
            out.append(('W', i, w, X[i] + ox, Y[i] + oy, float(l)))
        else:
            _, j, a, i, k = ct
            ox, oy = offset(C[j], Sn[j], a, 0.5)
            out.append(('C', j, i, X[j] + ox, Y[j] + oy, float(l)))
    so = rep.get('second_order') or {}
    return dict(S=float(P.S), sq=[(X[i], Y[i], float(P.T[i])) for i in range(P.n)], contacts=out,
                free=c.get('free') or [], flat=so.get('flat_squares') or [] if isinstance(so, dict) else [])
