#!/usr/bin/env python3
"""Two-phase certificate with a transport potential, odd n (2026-10-22).
Unknowns on a slope grid, per line-index parity j: p_j(m) in [0, 1] (phase) and psi_j(m) (free).
Step (m, parity j) -> (m', parity 1 - j):
    K(m, m') - lam >= a min(p, p') + g |p - p'| + psi_{1-j}(m') - psi_j(m)
Walls: index 0 (parity 0) and index n (parity 1 for odd n), p = 1, psi = 0.
Summing along an odd-n walk, psi telescopes to 0, so cost >= n lam + a sum min + g TV >= min(n, n lam + 2g) = min(n, F_n).
Reports the best uniform slack eps (>= 0: certificate on this grid)."""
import math, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from kwalk import Kmat

def solve(H, phis_deg, K=None, ms=None):
    if ms is None:
        ms = np.tan(np.radians(phis_deg)); ms[np.abs(ms) < 1e-12] = 0
    if K is None:
        K, _, _ = Kmat(H, ms)
    M = len(ms); r2 = math.sqrt(2)
    lam = 2 * r2 - H; a = 1 - lam; g = (4 - H - 2 * lam) / 2
    z = int(np.argmin(np.abs(ms)))
    # variable layout: p0, p1, s0, s1 (psi), eps
    p = lambda j, i: j * M + i
    s = lambda j, i: 2 * M + j * M + i
    nv = 4 * M + 1; E = nv - 1
    R, C, V, B = [], [], [], []
    r = [0]
    def add(terms, b):
        for c_, v_ in terms: R.append(r[0]); C.append(c_); V.append(v_)
        R.append(r[0]); C.append(E); V.append(1.0); B.append(b); r[0] += 1
    ii, kk = np.meshgrid(np.arange(M), np.arange(M), indexing='ij'); ii = ii.ravel(); kk = kk.ravel()
    for j in (0, 1):
        for i, k in zip(ii, kk):
            b = K[i, k] - lam
            u, v = p(j, i), p(1 - j, k)
            ps = [(s(1 - j, k), 1.0), (s(j, i), -1.0)]
            add([(u, g), (v, a - g)] + ps, b)
            add([(v, g), (u, a - g)] + ps, b)
    for i in range(M):
        b = K[z, i] - lam
        # start: wall (p=1, psi=0, parity 0) -> (m_i, parity 1):  f(1, p) + psi_1(m) <= b
        add([(p(1, i), a - g), (s(1, i), 1.0)], b - g)
        add([(p(1, i), g), (s(1, i), 1.0)], b - (a - g))
        # end: (m_i, parity 0) -> wall (parity 1):  f(p, 1) + 0 - psi_0(m) <= b
        add([(p(0, i), a - g), (s(0, i), -1.0)], b - g)
        add([(p(0, i), g), (s(0, i), -1.0)], b - (a - g))
    A = coo_matrix((V, (R, C)), shape=(r[0], nv)).tocsr()
    c = np.zeros(nv); c[E] = -1
    bounds = [(0, 1)] * (2 * M) + [(-10, 10)] * (2 * M) + [(None, 1)]
    res = linprog(c, A_ub=A, b_ub=np.array(B), bounds=bounds, method='highs')
    return res, ms

if __name__ == '__main__':
    phis = np.linspace(-84, 84, 337)
    for H in [float(x) for x in sys.argv[1:]] or [1.83, 1.85, 1.88, 1.9, 1.95, 1.99]:
        res, ms = solve(H, phis)
        print(f"H={H:.3f}: status {res.status} max slack eps = {res.x[-1]:+.3e}", flush=True)
