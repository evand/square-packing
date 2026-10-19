#!/usr/bin/env python3
"""Fixed-angle LP polish (10-18).  With every square's angle fixed, non-overlap of a pair is linear in the centres once
a separating axis is chosen (an edge normal of either square): sigma * d.(p_j - p_i) >= h_i(d) + h_j(d), with h the
constant half-width of a square projected on d.  Containment is linear too.  So min s is an LP per separator choice;
re-pick separators (largest SAT gap at the current point) and repeat until s stops falling.

At Pythagorean angles (cos, sin rational) the optimum of each LP is rational: the source of the rational records
(53/7 with 3-4-5, 723/41 with 20-21-29).

  polish(sq, s, iters=30) -> (sq, s)       sq: (n, 3) array of x, y, deg
"""
import math
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog


def axes(th):
    r = np.radians(th)
    return np.stack([np.cos(r), np.sin(r)], -1), np.stack([-np.sin(r), np.cos(r)], -1)


def halfw(th, d):
    """half extent of squares at angles th (k,) projected on unit directions d (k, 2)"""
    e1, e2 = axes(th)
    return 0.5 * (np.abs((e1 * d).sum(-1)) + np.abs((e2 * d).sum(-1)))


def pairs(sq, cut=1.5):
    P = sq[:, :2]
    D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1))
    i, j = np.nonzero(np.triu(D < cut, 1))
    return i, j


def separators(sq, i, j):
    """per pair: (direction d pointing i -> j side, required gap h_i + h_j, current slack) maximizing the SAT gap"""
    th = sq[:, 2]
    e1, e2 = axes(th)
    cand = np.stack([e1[i], e2[i], e1[j], e2[j]], 1)              # (m, 4, 2)
    dp = sq[j, :2] - sq[i, :2]
    proj = (cand * dp[:, None]).sum(-1)                           # (m, 4)
    sgn = np.where(proj >= 0, 1.0, -1.0)
    cand = cand * sgn[..., None]
    need = np.stack([halfw(th[i], cand[:, k]) + halfw(th[j], cand[:, k]) for k in range(4)], 1)
    gap = np.abs(proj) - need
    k = gap.argmax(1)
    r = np.arange(len(i))
    return cand[r, k], need[r, k], gap[r, k]


def lp(sq, i, j, d, need, trust=None, want_stress=False):
    n = len(sq); nv = 2 * n + 1                                   # x0, y0, x1, y1, ..., s
    th = sq[:, 2]
    hx = halfw(th, np.tile([1.0, 0.0], (n, 1))); hy = halfw(th, np.tile([0.0, 1.0], (n, 1)))
    rows, cols, vals, b = [], [], [], []
    m = len(i); r = np.arange(m)
    # -(d.(p_j - p_i)) <= -need
    for (idx, sign) in ((i, 1.0), (j, -1.0)):
        rows += [r, r]; cols += [2 * idx, 2 * idx + 1]; vals += [sign * d[:, 0], sign * d[:, 1]]
    b.append(-need)
    off = m
    k = np.arange(n)
    # -x <= -hx ; x - s <= -hx ; same y
    for c, h in ((0, hx), (1, hy)):
        rr = off + k; rows += [rr]; cols += [2 * k + c]; vals += [-np.ones(n)]; b.append(-h); off += n
        rr = off + k; rows += [rr, rr]; cols += [2 * k + c, np.full(n, 2 * n)]; vals += [np.ones(n), -np.ones(n)]
        b.append(-h); off += n
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(off, nv))
    cobj = np.zeros(nv); cobj[-1] = 1.0
    x0 = sq[:, :2].ravel()
    bounds = [(None, None)] * nv
    if trust is not None:
        bounds = [(x0[q] - trust, x0[q] + trust) for q in range(2 * n)] + [(0, None)]
    res = linprog(cobj, A_ub=A, b_ub=np.concatenate(b), bounds=bounds, method='highs')
    if res.status != 0:
        return None
    if want_stress:
        mu = np.abs(res.ineqlin.marginals)                        # contact forces of the side-pinning network
        st = np.zeros(n)
        np.add.at(st, i, mu[:m]); np.add.at(st, j, mu[:m])
        st += sum(mu[m + q * n: m + (q + 1) * n] for q in range(4))
        lp.stress = st
    return res.x


def polish(sq, s, iters=40, trust=0.25, cut=1.6, tol=1e-12, verbose=False):
    """-> (sq, s); polish.stress = per-square contact force sum at the last LP (the squares that pin the side)"""
    sq = np.array(sq, float).copy()
    best = s
    polish.stress = np.zeros(len(sq))
    for it in range(iters):
        i, j = pairs(sq, cut)
        d, need, gap = separators(sq, i, j)
        x = lp(sq, i, j, d, need, trust, want_stress=True)
        if x is not None:
            polish.stress = lp.stress
        if x is None:
            if verbose: print('LP infeasible at iter', it)
            break
        sq[:, :2] = x[:-1].reshape(-1, 2); snew = x[-1]
        if verbose: print(it, snew)
        if snew > best - tol:
            best = min(best, snew)
            if it > 2: break
        best = min(best, snew)
    return sq, best


def max_pen(sq, s):
    i, j = pairs(sq, 1.5)
    _, _, gap = separators(sq, i, j)
    th = sq[:, 2]
    hx = halfw(th, np.tile([1.0, 0.0], (len(sq), 1))); hy = halfw(th, np.tile([0.0, 1.0], (len(sq), 1)))
    wall = max((hx - sq[:, 0]).max(), (sq[:, 0] + hx - s).max(), (hy - sq[:, 1]).max(), (sq[:, 1] + hy - s).max())
    return max(-gap.min() if len(gap) else 0.0, wall)
