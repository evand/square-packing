#!/usr/bin/env python3
"""n unit squares in [0,w0] x [0,h] at FIXED tilts: min h over separation types (float LP explorer, 2026-10-20).

Separating-axis theorem for squares: two disjoint squares are separated along an edge normal of one of them, so
for fixed tilts the packing problem is a union over types (per pair: which square's axis, which of its two axes,
which sign) of LPs in the centres and h:
    walls  k_i <= x_i <= w0 - k_i,  k_i <= y_i <= h - k_i,   k_i = (|cos t_i| + |sin t_i|)/2
    pair   u . (c_j - c_i) >= s_i(u) + s_j(u),  s(u) = (|u.a| + |u.b|)/2 (support of a unit square, axes a, b)
min h over types = the least height at width w0 with these tilts (inf = infeasible).
"""
import itertools, math, sys
import numpy as np
from scipy.optimize import linprog


def axes(t):
    c, s = math.cos(t), math.sin(t)
    return [(c, s), (-s, c)]


def supp(t, u):
    a, b = axes(t)
    return (abs(u[0] * a[0] + u[1] * a[1]) + abs(u[0] * b[0] + u[1] * b[1])) / 2


def min_h(tilts_deg, w0, hmax=10.0, want_type=False):
    n = len(tilts_deg); T = [math.radians(t) for t in tilts_deg]
    k = [supp(t, (1, 0)) for t in T]
    if max(2 * x for x in k) > w0 + 1e-15:
        return math.inf
    pairs = list(itertools.combinations(range(n), 2))
    opts = [(o, ax, sg) for o in (0, 1) for ax in (0, 1) for sg in (1, -1)]
    nv = 2 * n + 1                             # x_i, y_i, h
    best, bt = math.inf, None
    # bounds: x in [k, w0 - k]; y >= k; y <= h - k as rows
    bounds = [(k[i], w0 - k[i]) for i in range(n)] + [(k[i], None) for i in range(n)] + [(0, hmax)]
    base_A, base_b = [], []
    for i in range(n):
        r = [0.0] * nv; r[n + i] = 1; r[2 * n] = -1; base_A.append(r); base_b.append(-k[i])
    cobj = [0.0] * nv; cobj[2 * n] = 1
    for combo in itertools.product(opts, repeat=len(pairs)):
        A, b = list(base_A), list(base_b)
        for (i, j), (o, ax, sg) in zip(pairs, combo):
            u = axes(T[i] if o == 0 else T[j])[ax]; u = (sg * u[0], sg * u[1])
            need = supp(T[i], u) + supp(T[j], u)
            r = [0.0] * nv                       # -(u.(c_j - c_i)) <= -need
            r[i] += u[0]; r[n + i] += u[1]; r[j] -= u[0]; r[n + j] -= u[1]
            A.append(r); b.append(-need)
        res = linprog(cobj, A_ub=A, b_ub=b, bounds=bounds, method='highs')
        if res.status == 0 and res.fun < best:
            best, bt = res.fun, combo
    return (best, bt) if want_type else best


if __name__ == '__main__':
    w0 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.95
    print('target 4+2sqrt2-2w0 =', 4 + 2 * math.sqrt(2) - 2 * w0)
    for t in [(0, 45, 0), (0, 40, 0), (0, 50, 0), (5, 45, 0), (0, 45, 5), (0, 0, 0), (45, 45, 45), (0, 30, 0)]:
        print(t, min_h(t, w0))
