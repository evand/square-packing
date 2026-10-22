#!/usr/bin/env python3
"""Two-phase certificate for the slope walk, odd n (2026-10-22).
Find p0, p1: slopes -> [0, 1] (line-index parity) with, for every step (m at parity j) -> (m' at parity 1 - j),
   K(m, m') - lam >= a min(p, p') + g |p - p'|        (both orders; = max of two linear forms since 2g >= a)
and wall steps  wall -> (m at parity 1),  (m at parity 0) -> wall  with wall p = 1.
lam = 2 sqrt2 - H, a = 1 - lam, 2g = 4 - H - 2 lam.  Then every odd-n walk costs >= min(n, n lam + 2g) = min(n, F_n).
Reports the largest uniform slack eps (constraints with - eps), so eps >= 0 means a certificate exists on this grid."""
import math, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from kwalk import Kmat

def solve(H, phis_deg, parity=True):
    ms = np.tan(np.radians(phis_deg)); ms[np.abs(ms) < 1e-12] = 0
    K, _, _ = Kmat(H, ms)
    M = len(ms); r2 = math.sqrt(2)
    lam = 2 * r2 - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    nv = 2 * M + 1                                 # p0, p1, eps
    P = lambda j, i: (j if parity else 0) * M + i
    rows, cols, vals, rhs = [], [], [], []
    r = 0
    def add(terms, b):
        nonlocal r
        for c, v in terms: rows.append(r); cols.append(c); vals.append(v)
        rows.append(r); cols.append(nv - 1); vals.append(1.0)     # + eps <= b
        rhs.append(b); r += 1
    I, J = np.meshgrid(np.arange(M), np.arange(M), indexing='ij')
    for j in (0, 1):
        for i, k in zip(I.ravel(), J.ravel()):
            b = K[i, k] - lam
            u, v = P(j, i), P(1 - j, k)
            # max(g p_u + (a-g) p_v, g p_v + (a-g) p_u) <= b
            add([(u, g), (v, a - g)], b) if u != v else add([(u, a)], b)
            if u != v: add([(v, g), (u, a - g)], b)
    z = int(np.argmin(np.abs(ms)))
    for i in range(M):                              # wall (p = 1) -> parity-1 line, parity-0 line -> wall
        for j in ((1, 0) if parity else (0,)):
            b = K[z, i] - lam
            u = P(j, i)
            # max(g*1 + (a-g) p, g p + (a-g)*1) <= b
            add([(u, a - g)], b - g); add([(u, g)], b - (a - g))
    A = coo_matrix((vals, (rows, cols)), shape=(r, nv)).tocsr()
    c = np.zeros(nv); c[-1] = -1
    bounds = [(0, 1)] * (nv - 1) + [(None, 1)]
    res = linprog(c, A_ub=A, b_ub=np.array(rhs), bounds=bounds, method='highs')
    return res, ms

if __name__ == '__main__':
    phis = np.linspace(-84, 84, 337)    # 0.5 deg, includes 0, +-45
    for H in [float(x) for x in sys.argv[1:]] or [1.83, 1.85, 1.88, 1.9, 1.95, 1.99]:
        for par in (False, True):
            res, ms = solve(H, phis, par)
            print(f"H={H:.3f} parity={par}: status {res.status} max slack eps = {res.x[-1]:+.3e}", flush=True)
