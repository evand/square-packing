#!/usr/bin/env python3
"""Chain relaxation for n unit squares in a strip of height H < 2, all n at once (2026-10-22).
Relaxation (as r4/relax.py, r5/relax5.py): squares sorted by centre x; end squares replaced by upright shelves
(left shelf on square 2, right shelf on square n-1, bands free); middle squares 2..n-1 with tilt t and height y;
only consecutive middle pairs kept disjoint (x order).  L >= S(s_2) + sum sigma(s_i, s_{i+1}) + E(s_{n-1}).
Computed as a min-plus DP on a (tilt, height) grid: an upper estimate of the relaxation's minimum (grid), then
polished by local optimisation (polish.py)."""
import math, sys
import numpy as np

def grid(H, nt=61, ny=33):
    ts = np.linspace(0, math.pi / 2, nt, endpoint=False)
    T, Y = [], []
    for t in ts:
        k = (math.cos(t) + math.sin(t)) / 2
        if 2 * k > H: continue
        for u in np.linspace(0, 1, ny):
            T.append(t); Y.append(k + u * (H - 2 * k))
    return np.array(T), np.array(Y)

def supp(t, nx, ny):
    C, S = np.cos(t), np.sin(t)
    return np.abs(nx * C + ny * S) + np.abs(-nx * S + ny * C)

def sigma_mat(T, Y):
    """sigma[i, j] = least x_j - x_i (>= 0) with squares i, j disjoint, x_i <= x_j."""
    tb, tc = T[:, None], T[None, :]
    dy = Y[None, :] - Y[:, None]
    best = np.full(dy.shape, np.inf)
    for t in (tb, tc):
        C, S = np.cos(t), np.sin(t)
        for nx, ny in ((C, S), (-C, -S), (-S, C), (S, -C)):
            ok = nx > 1e-9
            val = ((supp(tb, nx, ny) + supp(tc, nx, ny)) / 2 - ny * dy) / np.where(ok, nx, 1)
            best = np.where(ok, np.minimum(best, val), best)
    return np.maximum(best, 0)

def minx_window(t, y0, y1):
    C, S = math.cos(t), math.sin(t)
    e, f = np.array([C, S]) / 2, np.array([-S, C]) / 2
    V = [e + f, -e + f, -e - f, e - f]
    best = None
    for i in range(4):
        p, q = V[i], V[(i + 1) % 4]
        cands = []
        if y0 <= p[1] <= y1: cands.append(p[0])
        for yy in (y0, y1):
            if (p[1] - yy) * (q[1] - yy) < 0:
                s = (yy - p[1]) / (q[1] - p[1]); cands.append(p[0] + s * (q[0] - p[0]))
        for x in cands:
            best = x if best is None else min(best, x)
    return best

def start_cost(t, y, H, nlo=81):
    """least x of a square (tilt t, height y) right of an upright shelf [0,1] x [lo, lo+1], lo free."""
    k = (abs(math.cos(t)) + abs(math.sin(t))) / 2
    best = math.inf
    for lo in np.linspace(0, H - 1, nlo):
        m = minx_window(t, lo - y, lo + 1 - y)
        best = min(best, k if m is None else max(k, 1 - m))
    return best

def dp(H, nmax, nt=61, ny=33):
    T, Y = grid(H, nt, ny)
    Sg = sigma_mat(T, Y)
    St = np.array([start_cost(t, y, H) for t, y in zip(T, Y)])
    En = np.array([start_cost(-t, y, H) for t, y in zip(T, Y)])   # mirror x -> -x
    out = {}
    V = St.copy()                    # value with the current last middle square
    back = []
    for n in range(3, nmax + 1):     # n squares = n - 2 middle squares
        out[n] = float(np.min(V + En))
        back.append(V.copy())
        V = np.min(V[:, None] + Sg, axis=0)
    return out, (T, Y, Sg, St, En, back)

if __name__ == '__main__':
    sys.path.insert(0, '.')
    from targets import target
    nmax = 16
    for H in [float(a) for a in sys.argv[1:]] or [1.80, 1.84, 1.86, 1.88, 1.90, 1.93, 1.96, 1.99]:
        out, _ = dp(H, nmax)
        print(f"H={H:.3f} " + " ".join(f"{n}:{out[n] - target(n, H):+.4f}" for n in range(3, nmax + 1)), flush=True)
