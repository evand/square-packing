#!/usr/bin/env python3
"""Float feasibility of sub-configurations of an n = 17 closed-cell state at cap U (penalty minimisation).

Squares are placed with centres in their own exp-247 cells (all 24 cells are convex polygons), inside [0, U]^2,
pairwise disjoint (separating-axis penetration depth).  Reports the best total squared violation over multistarts:
~0 means a placement exists (to float accuracy), clearly positive suggests (does not prove) infeasibility.
Usage: subsets.py MASK [--only cell,cell,...] [--drop cell,...] [--starts N] [--cap U]"""
import sys, os, json, argparse
import numpy as np
from fractions import Fraction as Fr
from scipy.optimize import minimize

COVER = os.path.expanduser('~/math/_untrusted-third-party/jlevy-squares/packing/campaign/series/'
                           'series-000-smoke-and-calibration/results/exp-247-n17-unique-state-cover/run-001/receipt.json')


def load():
    r = json.load(open(COVER))
    cells = {}
    names = []
    for c in r['cells']:
        V = np.array([[float(Fr(a)), float(Fr(b))] for a, b in c['vertices']])
        # half-planes a.x <= b for a counter-clockwise polygon
        if np.cross(V[1] - V[0], V[2] - V[1]) < 0:
            V = V[::-1]
        A, B = [], []
        for i in range(len(V)):
            P, Q = V[i], V[(i + 1) % len(V)]
            e = Q - P
            n = np.array([e[1], -e[0]])          # outward for CCW
            n /= np.linalg.norm(n)
            A.append(n); B.append(n @ P)
        cells[c['name']] = (np.array(A), np.array(B), V)
        names.append(c['name'])
    return cells, names


def violation(z, A_list, B_list, U):
    k = len(A_list)
    X, Y, T = z[:k], z[k:2 * k], z[2 * k:]
    v = 0.0
    hw = 0.5 * (np.abs(np.cos(T)) + np.abs(np.sin(T)))
    v += np.sum(np.maximum(0, hw - X) ** 2 + np.maximum(0, X + hw - U) ** 2
                + np.maximum(0, hw - Y) ** 2 + np.maximum(0, Y + hw - U) ** 2)
    for i in range(k):
        c = A_list[i] @ np.array([X[i], Y[i]]) - B_list[i]
        v += np.sum(np.maximum(0, c) ** 2)
    i, j = np.triu_indices(k, 1)
    dx, dy = X[j] - X[i], Y[j] - Y[i]
    pen = np.full(len(i), np.inf)
    for phi in (T[i], T[i] + np.pi / 2, T[j], T[j] + np.pi / 2):
        nx, ny = np.cos(phi), np.sin(phi)
        H = 0.5 * (np.abs(np.cos(phi - T[i])) + np.abs(np.sin(phi - T[i]))) + \
            0.5 * (np.abs(np.cos(phi - T[j])) + np.abs(np.sin(phi - T[j])))
        pen = np.minimum(pen, H - np.abs(nx * dx + ny * dy))
    v += np.sum(np.maximum(0, pen) ** 2)
    return v


def solve(sel, cells, U, starts=40, seed=0):
    rng = np.random.default_rng(seed)
    A_list = [cells[c][0] for c in sel]
    B_list = [cells[c][1] for c in sel]
    best = (np.inf, None)
    k = len(sel)
    for s in range(starts):
        X0 = [cells[c][2][rng.integers(len(cells[c][2]))] * 0 + cells[c][2].mean(0) +
              rng.normal(0, 0.12, 2) for c in sel]
        z0 = np.concatenate([[p[0] for p in X0], [p[1] for p in X0], rng.uniform(0, np.pi / 2, k)])
        f = lambda z: violation(z, A_list, B_list, U)
        r = minimize(f, z0, method='L-BFGS-B', options=dict(maxiter=20000, maxfun=10 ** 6, ftol=1e-16, gtol=1e-12))
        for _ in range(3):
            r = minimize(f, r.x, method='L-BFGS-B', options=dict(maxiter=20000, maxfun=10 ** 6, ftol=1e-16, gtol=1e-12))
        if r.fun < best[0]:
            best = (r.fun, r.x)
        if best[0] < 1e-12:
            break
    return best


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('mask', type=int)
    a.add_argument('--only')
    a.add_argument('--drop')
    a.add_argument('--starts', type=int, default=40)
    a.add_argument('--cap', default='1169/250')
    o = a.parse_args()
    cells, names = load()
    U = float(Fr(o.cap))
    sel = [names[i] for i in range(24) if o.mask >> i & 1]
    if o.only:
        sel = [c for c in sel if c in o.only.split(',')]
    if o.drop:
        sel = [c for c in sel if c not in o.drop.split(',')]
    v, z = solve(sel, cells, U, o.starts)
    print(json.dumps(dict(mask=o.mask, cap=o.cap, k=len(sel), cells=sel, best_violation=v,
                          sqrt=float(np.sqrt(v)))), flush=True)
