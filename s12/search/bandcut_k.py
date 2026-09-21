#!/usr/bin/env python3
"""bandcut_k.py -- delta* with exactly k near-axis squares, with and without a chain (2026-09-20).

THE QUESTION (`tasks/bandcut-k/README.md`).  Fix eps > 0 and call a square *near-axis* if its tilt
(distance of its angle to 0 mod 90 deg) is < eps.  Split Theta_T = [0,90)^n by the number k of
near-axis squares.  In the hole  T <= k <= (T-1)^2  chain counting does NOT force a wall-to-wall
chain of T among the near-axis squares (`notes/proof-architecture.md` sec 0a items 1-2).  Measure:

    (a)  max delta* over angle vectors with exactly k near-axis squares            -- unrestricted
    (b)  the same, restricted to configurations whose near-axis squares contain
         NO chain of T                                                             -- c(k, eps)

WHAT A CHAIN IS.  For two squares i, j the separating-axis gap along one of the four edge normals
is  |n.(c_j - c_i)| - m_ij,  m_ij = 1/2 + (|cos D| + |sin D|)/2.  A normal is *x-type* if it is
closer to the global x axis than to y (for tilt < eps < 45 deg this is unambiguous).  Put an edge
i -> j in the x-DAG when some x-type normal separates the pair (gap >= -tol) with j on the + side.
A *chain of T* is a path on T squares in the x-DAG (or in the y-DAG); it need not be a straight row
(`search/RANK8.md` sec 3.2: one of the T = 4 certificate chains is a staircase).

WHY A CHAIN OF T FORCES delta <= 0 (the easy direction, exact).  Along the chain the T-1 link rows
give  sum n_s.(c_{s+1} - c_s) >= sum m >= T - 1  (m >= 1 always), the first and last square are
inside the box, and for a COMMON tilt t this telescopes to the identity of `search/ARCH_TLEDGER.md`
sec 2.1,  delta <= [(T - u) cos t + R sin t]/(T-1) - 1  with u = cos t + sin t and R the rise.

HOW CHAIN-FREENESS IS ENFORCED (complete, by Mirsky).  x-height <= T-1 and y-height <= T-1 among
the near-axis squares iff there are level maps L, M : near -> {0..T-2} with L(i) < L(j) for every
x-edge i -> j and M(i) < M(j) for every y-edge.  Any two near-axis squares are separated along an
x-type or a y-type normal, so (L, M) is injective: a chain-free configuration is exactly one that
admits an injective *pattern* P : near -> {0..T-2}^2 with

    L(i) < L(j)  =>  no x-edge j -> i,     L(i) = L(j)  =>  no x-edge either way,

and the same in y.  "No x-edge j -> i" is the linear-in-the-centres row  m_ij - n_{o,0}.(c_i - c_j)
>= eta  for both owners o.  So the chain-free maximum is a maximum over the C((T-1)^2, k) patterns
(deduplicated by the dihedral group of the level grid) of a smooth program, and NOTHING is assumed
about the far squares -- they are free obstacles at any angle with tilt >= eps.

EVERY NUMBER OUT OF THIS FILE IS A FEASIBLE POINT, i.e. a LOWER bound on the relevant sup; an
outlier above a trend convicts the trend.  Reported values are re-verified independently of the
pattern rows: the tilt band is checked, and chain-freeness is checked by walking the DAG.

    python3 search/bandcut_k.py scan --n 6 --T 3 --eps 1,5,10 --starts 600 --out runs/bandcut_k_T3.jsonl
    python3 search/bandcut_k.py report runs/bandcut_k_T3.jsonl
"""
import argparse
import itertools
import json
import math
import os
import sys

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6skel                                                        # noqa: E402

SGN = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def tilt(th):
    return abs(s6skel.norm_tilt(th))


# ===================================================================== chains (verification)
def normal_gaps(z, n):
    """-> gx, dirx, gy, diry.  gx[i,j] = best separating gap along an x-TYPE edge normal of i or j;
    dirx[i,j] = +1 if that normal puts j on the + (global +x) side, -1 otherwise.  Likewise y."""
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    C, S = np.cos(TH), np.sin(TH)
    gx = np.full((n, n), -1e18)
    gy = np.full((n, n), -1e18)
    dx_ = np.zeros((n, n), dtype=int)
    dy_ = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            ddx, ddy = X[j] - X[i], Y[j] - Y[i]
            D = TH[j] - TH[i]
            m = 0.5 + 0.5 * (abs(math.cos(D)) + abs(math.sin(D)))
            for o in (i, j):
                for kind, (d0, d1) in enumerate(((C[o], S[o]), (-S[o], C[o]))):
                    p = d0 * ddx + d1 * ddy
                    g = abs(p) - m
                    if abs(d0) > abs(d1):               # x-type normal
                        if g > gx[i, j]:
                            gx[i, j] = gx[j, i] = g
                            s = 1 if p * d0 > 0 else -1
                            dx_[i, j], dx_[j, i] = s, -s
                    else:                               # y-type normal
                        if g > gy[i, j]:
                            gy[i, j] = gy[j, i] = g
                            s = 1 if p * d1 > 0 else -1
                            dy_[i, j], dy_[j, i] = s, -s
    return gx, dx_, gy, dy_


def longest_path(idx, g, dr, coord, lvl):
    """longest path (in vertices) of the DAG on `idx` with edge i->j iff g[i,j] >= lvl and
    dr[i,j] = +1.  Topological order = increasing `coord` (valid for near-axis squares)."""
    order = sorted(idx, key=lambda i: coord[i])
    best = {i: 1 for i in order}
    pred = {i: None for i in order}
    for b, j in enumerate(order):
        for i in order[:b]:
            if g[i, j] >= lvl and dr[i, j] > 0 and best[i] + 1 > best[j]:
                best[j] = best[i] + 1
                pred[j] = i
    if not order:
        return 0, []
    end = max(order, key=lambda i: best[i])
    path, v = [], end
    while v is not None:
        path.append(v)
        v = pred[v]
    return best[end], list(reversed(path))


def chain_report(z, n, near, lvl, tol=1e-7):
    """-> dict with the x/y heights among `near`, among all n, and the longest paths.

    THE LEVEL.  An edge counts when its separating gap is at least `lvl`, and `lvl` is the
    configuration's own margin delta (minus a tolerance).  That is the definition under which the
    chain lemma is a theorem: the T-1 links of a chain each give  n.(c_{s+1} - c_s) >= m + delta
    and the two end walls give delta each, so a chain of T of squares at a COMMON tilt t forces
    delta <= B_T(R,t) (ARCH_TLEDGER sec 2.1), which is <= 0 for R <= T tan(t/2) + cos t - sin t.
    Defining edges at level 0 instead would be wrong: a chain whose links are at gap exactly
    delta < 0 would not be seen, and the optimiser would 'break' chains for free by overlapping
    one link by 1e-6."""
    X, Y = z[1::3], z[2::3]
    gx, dxr, gy, dyr = normal_gaps(z, n)
    L = lvl - tol
    hxn, px = longest_path(list(near), gx, dxr, X, L)
    hyn, py = longest_path(list(near), gy, dyr, Y, L)
    hxa, pxa = longest_path(list(range(n)), gx, dxr, X, L)
    hya, pya = longest_path(list(range(n)), gy, dyr, Y, L)
    out = dict(hx_near=hxn, hy_near=hyn, hx_all=hxa, hy_all=hya,
               path_x=px, path_y=py, path_x_all=pxa, path_y_all=pya)
    for nm, pth, co in (('x', px, X), ('y', py, Y)):
        if len(pth) >= 2:
            perp = (Y if nm == 'x' else X)
            out['rise_' + nm] = float(abs(perp[pth[-1]] - perp[pth[0]]))
            out['span_' + nm] = float(co[pth[-1]] - co[pth[0]])
    return out


# ===================================================================== the band problem
class BandProblem(s6skel.Problem):
    """s6skel.Problem (max-min gap over centres AND angles) plus HARD rows (no delta in them):

       near i:  |th_i| <= eps                  far i:  eps <= sgn_i * th_i <= 45 deg
       pattern: m_ij - s * n_{o,kind}.(c_j - c_i) >= eta     for the forbidden separations.
    """

    def __init__(self, n, T, k, eps, signs, pattern=None, eta=1e-6, maxtilt=math.pi / 4):
        super().__init__(n, T)
        self.k, self.eps, self.eta, self.maxtilt = k, eps, eta, maxtilt
        self.near = list(range(k))
        self.far = list(range(k, n))
        self.signs = np.array(signs, dtype=float)          # length n - k, +-1
        self.pattern = pattern
        rows = []            # (i, j, o, kind, s)  ->  m_ij - s * n_{o,kind}.(c_j - c_i) >= eta
        if pattern is not None:
            for a in range(k):
                for b in range(a + 1, k):
                    (la, ma), (lb, mb) = pattern[a], pattern[b]
                    for (kind, ca, cb) in ((0, la, lb), (1, ma, mb)):
                        for o in (a, b):
                            if ca <= cb:                   # forbid b -> a  (a on + side of b)
                                rows.append((a, b, o, kind, -1))
                            if ca >= cb:                   # forbid a -> b
                                rows.append((a, b, o, kind, +1))
        self.Ri = np.array([r[0] for r in rows], dtype=int)
        self.Rj = np.array([r[1] for r in rows], dtype=int)
        self.Ro = np.array([r[2] for r in rows], dtype=int)
        self.Rk = np.array([r[3] for r in rows], dtype=int)
        self.Rs = np.array([r[4] for r in rows], dtype=float)
        self.Rsd = np.ones(len(rows))                      # sign of th_j - th_i, fixed per round
        self.nextra = 2 * n + len(rows)

    def set_branch(self, z):
        TH = z[3::3]
        if len(self.Ri):
            d = TH[self.Rj] - TH[self.Ri]
            self.Rsd = np.where(d >= 0.0, 1.0, -1.0)

    # ---- hard rows
    def extra(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        k, n = self.k, self.n
        out = [self.eps - TH[:k], self.eps + TH[:k],
               self.signs * TH[k:] - self.eps, self.maxtilt - self.signs * TH[k:]]
        if len(self.Ri):
            C, S = np.cos(TH), np.sin(TH)
            co, so = C[self.Ro], S[self.Ro]
            d0 = np.where(self.Rk == 0, co, -so)
            d1 = np.where(self.Rk == 0, so, co)
            ddx = X[self.Rj] - X[self.Ri]
            ddy = Y[self.Rj] - Y[self.Ri]
            D = TH[self.Rj] - TH[self.Ri]
            m = 0.5 + 0.5 * (np.cos(D) + self.Rsd * np.sin(D))
            # the forbidden separation must be WORSE than the margin:  gap <= g - eta
            out.append(m - self.Rs * (d0 * ddx + d1 * ddy) + z[0] - self.eta)
        return np.concatenate(out)

    def extra_jac(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        k, n, nv = self.k, self.n, self.nv
        J = np.zeros((self.nextra, nv))
        for a in range(k):
            J[a, 3 + 3 * a] = -1.0
            J[k + a, 3 + 3 * a] = +1.0
        for b, i in enumerate(range(k, n)):
            J[2 * k + b, 3 + 3 * i] = self.signs[b]
            J[2 * k + (n - k) + b, 3 + 3 * i] = -self.signs[b]
        if len(self.Ri):
            r0 = 2 * n
            C, S = np.cos(TH), np.sin(TH)
            co, so = C[self.Ro], S[self.Ro]
            d0 = np.where(self.Rk == 0, co, -so)
            d1 = np.where(self.Rk == 0, so, co)
            dp0 = np.where(self.Rk == 0, -so, -co)
            dp1 = np.where(self.Rk == 0, co, -so)
            ddx = X[self.Rj] - X[self.Ri]
            ddy = Y[self.Rj] - Y[self.Ri]
            D = TH[self.Rj] - TH[self.Ri]
            rows = r0 + np.arange(len(self.Ri))
            J[rows, 0] = 1.0
            np.add.at(J, (rows, 1 + 3 * self.Rj), -self.Rs * d0)
            np.add.at(J, (rows, 2 + 3 * self.Rj), -self.Rs * d1)
            np.add.at(J, (rows, 1 + 3 * self.Ri), self.Rs * d0)
            np.add.at(J, (rows, 2 + 3 * self.Ri), self.Rs * d1)
            np.add.at(J, (rows, 3 + 3 * self.Ro), -self.Rs * (dp0 * ddx + dp1 * ddy))
            dm = 0.5 * (-np.sin(D) + self.Rsd * np.cos(D))
            np.add.at(J, (rows, 3 + 3 * self.Rj), dm)
            np.add.at(J, (rows, 3 + 3 * self.Ri), -dm)
        return J

    def cons_all(self, z):
        return np.concatenate([self.cons(z), self.extra(z)])

    def jac_all(self, z):
        return np.vstack([self.jac(z), self.extra_jac(z)])


class FixedBand:
    """the disjunctive LP in the CENTRES at a fixed angle vector (s6local.Fixed) plus the chain
    rows, which are also linear in the centres.  HiGHS solves this to optimality for a fixed
    assignment, so it is a far more reliable polish than SLSQP; the angles do not move."""

    def __init__(self, theta, T, k, pattern, eta):
        from scipy.optimize import linprog                     # noqa: F401  (import cost once)
        self.theta, self.T, self.n = list(theta), T, len(theta)
        n = self.n
        self.NV = 1 + 2 * n
        self.cs = [(math.cos(t), math.sin(t)) for t in theta]
        self.P = [0.5 * (abs(c) + abs(s)) for (c, s) in self.cs]
        self.pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        self.m = {(i, j): 0.5 + 0.5 * (abs(math.cos(theta[j] - theta[i]))
                                       + abs(math.sin(theta[j] - theta[i])))
                  for (i, j) in self.pairs}
        A, b, tags = [], [], []
        for i in range(n):
            for (off, ax) in ((1, 'x'), (1 + n, 'y')):
                r = np.zeros(self.NV); r[off + i] = 1.0; r[0] = -1.0
                A.append(r); b.append(self.P[i]); tags.append(('lo', ax, i))
                r = np.zeros(self.NV); r[off + i] = -1.0; r[0] = -1.0
                A.append(r); b.append(self.P[i] - T); tags.append(('hi', ax, i))
        if pattern is not None:
            for a in range(k):
                for bb in range(a + 1, k):
                    (la, ma), (lb, mb) = pattern[a], pattern[bb]
                    for (kind, ca, cb) in ((0, la, lb), (1, ma, mb)):
                        for o in (a, bb):
                            for s in ((-1,) if ca < cb else (1,) if ca > cb else (-1, 1)):
                                cc, ss = self.cs[o]
                                d0, d1 = ((cc, ss), (-ss, cc))[kind]
                                r = np.zeros(self.NV); r[0] = 1.0
                                r[1 + bb] -= s * d0; r[1 + a] += s * d0
                                r[1 + n + bb] -= s * d1; r[1 + n + a] += s * d1
                                A.append(r); b.append(eta - self.m[(a, bb)])
                                tags.append(('chain', 'xy'[kind], a, bb, o, s))
        self.baseA, self.baseb, self.basetags = A, b, tags
        # delta is carried by the WALL rows too (S6_SKELETON sec 3.1), unlike s6local.Fixed, which
        # uses hard containment; the two differ exactly where delta < 0, which is this whole file.
        self.bounds = [(-2.0, 0.5)] + [(-1.0, T + 1.0)] * (2 * n)
        self.c = np.zeros(self.NV); self.c[0] = -1.0

    def assign(self, X, Y):
        out = []
        for (i, j) in self.pairs:
            dx, dy = X[j] - X[i], Y[j] - Y[i]
            bv, ba = -1e18, None
            for o in (i, j):
                cc, ss = self.cs[o]
                for kd, (d0, d1) in enumerate(((cc, ss), (-ss, cc))):
                    pr = d0 * dx + d1 * dy
                    for sg in (1, -1):
                        if sg * pr > bv:
                            bv, ba = sg * pr, (o, kd, sg)
            out.append(ba)
        return tuple(out)

    def rows(self, asg):
        n = self.n
        A, b, tags = list(self.baseA), list(self.baseb), list(self.basetags)
        for (i, j), (o, kd, sg) in zip(self.pairs, asg):
            cc, ss = self.cs[o]
            d0, d1 = ((cc, ss), (-ss, cc))[kd]
            r = np.zeros(self.NV)
            r[1 + j] = sg * d0; r[1 + i] = -sg * d0
            r[1 + n + j] = sg * d1; r[1 + n + i] = -sg * d1
            r[0] = -1.0
            A.append(r); b.append(self.m[(i, j)]); tags.append(('pair', i, j, o, kd, sg))
        return A, b, tags

    def lp(self, asg, full=False):
        from scipy.optimize import linprog
        A, b, tags = self.rows(asg)
        res = linprog(self.c, A_ub=-np.array(A), b_ub=-np.array(b), bounds=self.bounds,
                      method='highs',
                      options={'primal_feasibility_tolerance': 1e-10,
                               'dual_feasibility_tolerance': 1e-10})
        if not res.success:
            return (None, None, None) if full else None
        return (res.x, res, tags) if full else res.x

    def ascent(self, X, Y, iters=20):
        n = self.n
        z = np.zeros(1 + 3 * n)
        z[3::3] = self.theta
        best, bX, bY, prev = -1e18, None, None, None
        for _ in range(iters):
            asg = self.assign(X, Y)
            if asg == prev:
                break
            prev = asg
            x = self.lp(asg)
            if x is None:
                break
            X, Y = x[1:1 + n].copy(), x[1 + n:1 + 2 * n].copy()
            z[1::3], z[2::3] = X, Y
            v = s6skel.value(z, n, self.T)
            if v > best:
                best, bX, bY = v, X.copy(), Y.copy()
        return best, bX, bY


def lp_polish(z, n, T, k, pattern, eta):
    """re-optimise the centres at the angles of `z`.  -> new z (or z if it failed)."""
    F = FixedBand(list(z[3::3]), T, k, pattern, eta)
    v, X, Y = F.ascent(z[1::3].copy(), z[2::3].copy())
    if X is None:
        return z, s6skel.value(z, n, T)
    w = z.copy()
    w[1::3], w[2::3] = X, Y
    w[0] = v
    return w, v


def solve_band(pr, z0, rounds=8, maxiter=200):
    from scipy.optimize import minimize
    n, T = pr.n, pr.T
    z = np.array(z0, dtype=float)
    z[0] = s6skel.value(z, n, T)
    best = (-1e18, None)
    seen = set()
    for _it in range(rounds):
        pa = pr.pair_assign(z)
        key = tuple(pa)
        pr.build(pa)
        pr.set_branch(z)
        res = minimize(lambda w: -w[0], z,
                       jac=lambda w: np.eye(1, len(w), 0).ravel() * -1.0,
                       constraints=[{'type': 'ineq', 'fun': pr.cons_all, 'jac': pr.jac_all}],
                       method='SLSQP', options={'maxiter': maxiter, 'ftol': 1e-14})
        zn = res.x
        val = s6skel.value(zn, n, T)
        if val > best[0]:
            best = (val, zn.copy())
        if key in seen:
            break
        seen.add(key)
        z = zn
        z[0] = s6skel.value(z, n, T)
    return best


# ===================================================================== verification
def verify(z, n, T, k, eps, val, tol_band=1e-9, tol_edge=1e-7, need_chainfree=False):
    """-> (ok, info).  Checks the tilt band and (if asked) chain-freeness among the near squares,
    both from the geometry alone -- the pattern rows are not trusted."""
    TH = z[3::3]
    tl = np.array([tilt(t) for t in TH])
    if k > 0 and float(np.max(tl[:k])) > eps + tol_band:
        return False, dict(reason='near tilt %.3e > eps' % math.degrees(float(np.max(tl[:k]))))
    if k < n and float(np.min(tl[k:])) < eps - tol_band:
        return False, dict(reason='far tilt %.3e < eps' % math.degrees(float(np.min(tl[k:]))))
    ch = chain_report(z, n, range(k), val, tol=tol_edge)
    ch['nboundary'] = int(sum(1 for t in tl if abs(t - eps) <= 1e-7))
    if need_chainfree and max(ch['hx_near'], ch['hy_near']) >= T - 0.5:
        return False, dict(reason='chain of %d present' % max(ch['hx_near'], ch['hy_near']), **ch)
    return True, ch


# ===================================================================== patterns
def level_patterns(T, k):
    """all injective near -> {0..T-2}^2 patterns of size k, up to the dihedral group of the grid."""
    L = int(round(T)) - 1
    cells = [(a, b) for a in range(L) for b in range(L)]
    if k > len(cells):
        return []
    def syms(s):
        out = []
        for (fx, fy, tr) in itertools.product((0, 1), (0, 1), (0, 1)):
            t = []
            for (a, b) in s:
                aa, bb = (L - 1 - a if fx else a), (L - 1 - b if fy else b)
                t.append((bb, aa) if tr else (aa, bb))
            out.append(tuple(sorted(t)))
        return out
    seen, out = set(), []
    for s in itertools.combinations(cells, k):
        key = min(syms(s))
        if key in seen:
            continue
        seen.add(key)
        out.append(list(s))
    return out


# ===================================================================== starts
def make_start(rng, n, T, k, eps, pattern, kind):
    """`kind`:  spread/compact/jitgrid  put the near squares on the (T-1)^2 LEVEL grid at three
    different scales (pattern runs only);  tilingall puts all n squares on distinct cells of the
    TxT tiling;  random is uniform.  Far tilts get a random sign -- sign coherence matters
    (`search/S6_LOCAL.md` sec 5)."""
    z = np.zeros(1 + 3 * n)
    Ti = int(round(T))
    L = Ti - 1
    cells = [(a + 0.5, b + 0.5) for a in range(Ti) for b in range(Ti)]
    pos = []
    if pattern is not None and kind in ('spread', 'compact', 'jitgrid'):
        if kind == 'spread':
            sc = T / L
        elif kind == 'compact':
            sc = 1.0
        else:
            sc = 1.0 + rng.uniform(0.0, T / L - 1.0)
        bx = rng.uniform(0.0, max(0.0, T - sc * (L - 1) - 1.0))
        by = rng.uniform(0.0, max(0.0, T - sc * (L - 1) - 1.0))
        for (a, b) in pattern:
            pos.append((bx + 0.5 + sc * a, by + 0.5 + sc * b))
        for _ in range(n - k):
            pos.append((rng.uniform(0.7, T - 0.7), rng.uniform(0.7, T - 0.7)))
    elif kind == 'tilingall':
        sub = rng.permutation(len(cells))[:n]
        pos = [cells[c] for c in sub]
    elif kind == 'chainrow':
        # the near squares on one full row or column of the tiling (a wall-to-wall chain start),
        # or, if k < Ti, on a run of k consecutive cells of one; the rest anywhere else.
        vert = rng.random() < 0.5
        line = int(rng.integers(Ti))
        run = list(range(Ti)) if k >= Ti else list(range(k))
        near = [((line + 0.5, q + 0.5) if vert else (q + 0.5, line + 0.5)) for q in run]
        near = near[:k] + [(rng.uniform(0.6, T - 0.6), rng.uniform(0.6, T - 0.6))
                           for _ in range(max(0, k - len(near)))]
        rest = [c for c in cells if c not in near]
        rng.shuffle(rest)
        pos = near + [rest[q % len(rest)] for q in range(n - k)]
    else:
        pos = [(rng.uniform(0.6, T - 0.6), rng.uniform(0.6, T - 0.6)) for _ in range(n)]
    fartilt = eps if rng.random() < 0.6 else rng.uniform(eps, math.pi / 4)
    coherent = rng.random() < 0.55           # one sign for every far square (S6_LOCAL sec 5.1)
    gsign = 1.0 if rng.random() < 0.5 else -1.0
    rho = 0.10 if kind != 'random' else 0.22
    for q in range(n):
        x, y = pos[q]
        z[1 + 3 * q] = x + rng.normal(0.0, rho)
        z[2 + 3 * q] = y + rng.normal(0.0, rho)
        if q < k:
            z[3 + 3 * q] = rng.uniform(-1.0, 1.0) * eps * 0.3
        else:
            s = gsign if coherent else (1.0 if rng.random() < 0.5 else -1.0)
            tq = fartilt if rng.random() < 0.75 else rng.uniform(eps, math.pi / 4)
            z[3 + 3 * q] = s * tq
    return z


def nudge(z, rng, n, T, k, eps, rho, alpha):
    w = z.copy()
    for q in range(n):
        w[1 + 3 * q] += rng.normal(0.0, rho)
        w[2 + 3 * q] += rng.normal(0.0, rho)
        th = w[3 + 3 * q] + rng.normal(0.0, alpha)
        if q < k:
            th = min(max(th, -eps), eps)
        else:
            s = 1.0 if w[3 + 3 * q] >= 0 else -1.0
            th = s * min(max(s * th, eps), math.pi / 4)
        w[3 + 3 * q] = th
    return w


def start_signs(z, n, k):
    TH = z[3::3]
    return [1.0 if TH[i] >= 0 else -1.0 for i in range(k, n)]


# ===================================================================== the scan
_G = {}


def _init(a):
    _G['a'] = a


def clamp_band(z, n, k, eps):
    w = z.copy()
    for q in range(n):
        th = s6skel.norm_tilt(w[3 + 3 * q])
        if q < k:
            th = min(max(th, -eps), eps)
        else:
            s = 1.0 if th >= 0 else -1.0
            th = s * min(max(abs(th), eps), math.pi / 4)
        w[3 + 3 * q] = th
    return w


def _job(arg):
    """one (k, pattern); sweeps eps down and up carrying the best configurations (continuation is
    what makes the n = 12 multistart usable, `search/S6_LOCAL.md` sec 1)."""
    (k, pidx, pattern, seed, nstart) = arg
    a = _G['a']
    n, T = a.n, a.T
    rng = np.random.default_rng(seed)
    L2 = (int(round(T)) - 1) ** 2
    pats = level_patterns(T, k) if (pattern is None and 0 < k <= L2) else None
    kinds = (['spread', 'compact', 'jitgrid', 'tilingall', 'random']
             if pattern is not None else
             ['tilingall', 'chainrow', 'random', 'spread'] if pats else
             ['tilingall', 'chainrow', 'random'])
    need = pattern is not None
    epsdegs = [float(s) for s in a.eps.split(',')]

    def geom():
        return pattern if pattern is not None else (
            pats[int(rng.integers(len(pats)))] if pats else None)

    def run(z0, eps):
        signs = start_signs(z0, n, k)
        pr = BandProblem(n, T, k, eps, signs, pattern=pattern, eta=a.eta)
        best = (-1e18, None)
        z, v = lp_polish(z0, n, T, k, pattern, a.eta)     # LP first: exact centres, angles fixed
        if v > best[0]:
            best = (v, z.copy())
        for _cyc in range(a.cycles):                      # then SLSQP (angles) / LP alternately
            v, z = solve_band(pr, z, rounds=a.rounds, maxiter=a.maxiter)
            if z is None:
                break
            if v > best[0]:
                best = (v, z.copy())
            z, v = lp_polish(z, n, T, k, pattern, a.eta)
            if v > best[0]:
                best = (v, z.copy())
        v, z = best
        if z is None:
            return None
        z[0] = v
        ok, info = verify(z, n, T, k, eps, v, tol_edge=a.tol_edge, need_chainfree=need)
        return (v, z, info) if ok else None

    table = {}
    carried = []
    order = sorted(epsdegs, reverse=True) + sorted(epsdegs)
    for epsdeg in order:
        eps = math.radians(epsdeg)
        starts = [make_start(rng, n, T, k, eps, geom(), kinds[s % len(kinds)])
                  for s in range(nstart)]
        for (_v, z, _i) in carried + [r for r in table.get(epsdeg, ([], []))[0]]:
            w = clamp_band(z, n, k, eps)
            starts.append(w)
            for (rho, alpha) in ((0.0, 0.0), (0.03, 0.3 * eps), (0.12, 0.7 * eps)):
                for _ in range(3):
                    starts.append(nudge(w, rng, n, T, k, eps, rho, alpha))
        keep, vals = [], []
        for z0 in starts:
            r = run(z0, eps)
            if r is None:
                continue
            vals.append(r[0])
            keep.append(r)
        keep.sort(key=lambda r: -r[0])
        keep = keep[:6]
        for _rd in range(a.polish):
            newk = list(keep)
            for (v, z, _i) in keep:
                for (rho, alpha) in ((0.02, 0.2 * eps), (0.08, 0.5 * eps), (0.2, eps)):
                    for _ in range(3):
                        r = run(nudge(z, rng, n, T, k, eps, rho, alpha), eps)
                        if r is not None:
                            vals.append(r[0])
                            newk.append(r)
            newk.sort(key=lambda r: -r[0])
            keep = newk[:6]
        old = table.get(epsdeg)
        if old is None or (keep and keep[0][0] > old[0][0][0]):
            table[epsdeg] = (keep, vals + (old[1] if old else []))
        elif old:
            table[epsdeg] = (old[0], old[1] + vals)
        carried = keep

    out = []
    for epsdeg in epsdegs:
        keep, vals = table.get(epsdeg, ([], []))
        pat = None if pattern is None else [list(c) for c in pattern]
        if not keep:
            out.append(dict(k=k, eps=epsdeg, pidx=pidx, value=None, nfeas=0, nstart=len(vals),
                            nhit=0, pattern=pat, z=None, info=None))
            continue
        v, z, info = keep[0]
        out.append(dict(k=k, eps=epsdeg, pattern=pat, pidx=pidx, value=float(v),
                        nfeas=len(vals), nstart=len(vals),
                        nhit=sum(1 for u in vals if u >= v - 1e-9),
                        z=[float(q) for q in z],
                        info={kk: (vv if not isinstance(vv, list) else list(vv))
                              for kk, vv in info.items()}))
    return out


def cmd_scan(a):
    import multiprocessing as mp
    n, T = a.n, a.T
    ks = [int(s) for s in a.ks.split(',')] if a.ks else list(range(n + 1))
    epss = [float(s) for s in a.eps.split(',')]
    jobs = []
    seed = a.seed * 1000003
    for k in ks:
        if a.mode in ('both', 'free'):
            jobs.append((k, -1, None, seed, a.starts))
            seed += 7919
        # for k < T no chain of T can exist among k squares, so the chain-free cell IS the free
        # cell; --patmink skips those pattern runs.
        if a.mode in ('both', 'chainfree') and a.patmink <= k <= (int(round(T)) - 1) ** 2:
            for pidx, pat in enumerate(level_patterns(T, k)):
                jobs.append((k, pidx, pat, seed, a.pstarts))
                seed += 7919
    print(f"# n={n} T={T} eps={epss} eta={a.eta}: {len(jobs)} jobs "
          f"({sum(1 for j in jobs if j[2] is None)} unrestricted), each sweeping eps",
          flush=True)
    with mp.Pool(a.nproc, initializer=_init, initargs=(a,)) as pool, open(a.out, 'w') as f:
        for i, rs in enumerate(pool.imap_unordered(_job, jobs)):
            for r in rs:
                f.write(json.dumps(r) + '\n')
            f.flush()
            r0 = rs[0]
            tag = 'free    ' if r0['pidx'] < 0 else f"pat{r0['pidx']:<4d}"
            vs = ' '.join(f"{r['eps']:g}:{'--' if r['value'] is None else round(r['value'], 8)}"
                          f"[{r['nhit']}/{r['nstart']}]" for r in rs)
            print(f"  [{i + 1}/{len(jobs)}] k={r0['k']:<2d} {tag} {vs}", flush=True)


# ===================================================================== reporting
def cmd_report(a):
    recs = []
    for f in a.file:
        recs += [json.loads(l) for l in open(f)]
    Ts = {}
    for r in recs:
        Ts.setdefault((r['eps'], r['k']), []).append(r)
    epss = sorted({r['eps'] for r in recs})
    ks = sorted({r['k'] for r in recs})
    print("# every entry is a FEASIBLE POINT: a LOWER bound on the cell's true maximum.")
    print("# (hit/total) counts the local optimisations that reached the value.")
    for tagname, sel in (('(a) unrestricted: max delta* with exactly k near-axis squares  '
                          '(best over the free AND the chain-free runs)', lambda r: True),
                         ('(b) chain-free: c(k,eps) = max delta* with NO chain of T among the '
                          'near-axis squares', lambda r: r['pidx'] >= 0)):
        print(f"\n## {tagname}")
        print("  k  | " + " | ".join(f"{'eps=%g deg' % e:^30s}" for e in epss))
        for k in ks:
            cells = []
            for e in epss:
                rs = [r for r in Ts.get((e, k), []) if sel(r) and r['value'] is not None]
                if not rs:
                    cells.append(' ' * 30)
                    continue
                b = max(rs, key=lambda r: r['value'])
                q = b['value'] / math.radians(e) ** 2
                cells.append(f"{b['value']:+.6e} /e^2 {q:+.4f} ({b['nhit']}/{b['nstart']})")
            print(f" {k:2d}  | " + " | ".join(cells))
    print("\n## gamma = -value/eps^2, and its linear-in-eps extrapolation to eps = 0")
    for tagname, sel in (('unrestricted', lambda r: True),
                         ('chain-free  ', lambda r: r['pidx'] >= 0)):
        print(f"  ## {tagname}")
        for k in ks:
            qs = []
            for e in epss:
                rs = [r for r in Ts.get((e, k), []) if sel(r) and r['value'] is not None]
                qs.append(max(r['value'] for r in rs) / math.radians(e) ** 2 if rs else None)
            if all(q is None for q in qs):
                continue
            s = f"   k={k:<2d} " + ' '.join('   --   ' if q is None else f"{q:+8.4f}" for q in qs)
            if qs[0] is not None and qs[1] is not None:
                e1, e2 = math.radians(epss[0]), math.radians(epss[1])
                q0 = qs[0] - (qs[1] - qs[0]) / (e2 - e1) * e1
                s += f"   -> eps=0: {q0:+.4f}" + (f"  = -1/{-1 / q0:.2f}" if q0 < -1e-6 else '')
            print(s)
    print("\n## the best chain-free pattern per cell (level-grid cells of the (T-1)x(T-1) grid)")
    for k in ks:
        for e in epss:
            rs = [r for r in Ts.get((e, k), []) if r['pidx'] >= 0 and r['value'] is not None]
            if not rs:
                continue
            b = max(rs, key=lambda r: r['value'])
            i = b['info'] or {}
            print(f"  k={k:<2d} eps={e:<4g} c={b['value']:+.6e}  pattern={b['pattern']}  "
                  f"near heights x={i.get('hx_near')} y={i.get('hy_near')}, all-squares "
                  f"x={i.get('hx_all')} y={i.get('hy_all')}")
    print("\n## chain structure of the unrestricted optima (heights among near-axis / among all)")
    for k in ks:
        for e in epss:
            rs = [r for r in Ts.get((e, k), []) if r['value'] is not None]
            if not rs:
                continue
            b = max(rs, key=lambda r: r['value'])
            i = b['info'] or {}
            rx, ry = i.get('rise_x'), i.get('rise_y')
            print(f"  k={k:<2d} eps={e:<4g} value={b['value']:+.3e}  near x={i.get('hx_near')} "
                  f"y={i.get('hy_near')}   all x={i.get('hx_all')} y={i.get('hy_all')}   "
                  f"rise x={'--' if rx is None else round(rx, 5)} "
                  f"y={'--' if ry is None else round(ry, 5)}"
                  + ('   [chain-free run]' if b['pidx'] >= 0 else ''))


def cmd_show(a):
    recs = [json.loads(l) for l in open(a.file)]
    recs = [r for r in recs if r['value'] is not None and r['k'] == a.k and r['eps'] == a.eps
            and ((r['pidx'] >= 0) == (not a.free))]
    recs.sort(key=lambda r: -r['value'])
    for r in recs[:a.top]:
        z = np.array(r['z'])
        n = (len(z) - 1) // 3
        print(f"# k={r['k']} eps={r['eps']} pattern={r['pattern']} value={r['value']:+.9e} "
              f"feasible {r['nfeas']}/{r['nstart']} hit {r['nhit']}  info={r['info']}")
        order = sorted(range(n), key=lambda q: (round(z[2 + 3 * q] * 2) / 2, z[1 + 3 * q]))
        for q in order:
            print(f"    sq{q:2d}  x={z[1 + 3 * q]:10.6f}  y={z[2 + 3 * q]:10.6f}  "
                  f"tilt={math.degrees(s6skel.norm_tilt(z[3 + 3 * q])):+8.4f}"
                  f"{'  NEAR' if q < r['k'] else ''}")


def all_chains(idx, g, dr, coord, lvl, L):
    """every path on exactly L vertices of the DAG (edge i->j iff g >= lvl and dr = +1)."""
    order = sorted(idx, key=lambda i: coord[i])
    adj = {i: [j for j in order if g[i, j] >= lvl and dr[i, j] > 0 and coord[j] > coord[i]]
           for i in order}
    out = []

    def walk(path):
        if len(path) == L:
            out.append(list(path))
            return
        for j in adj[path[-1]]:
            if j not in path:
                walk(path + [j])
    for i in order:
        walk([i])
    return out


def chain_certificate(z, n, near, lvl, T, tol=1e-7):
    """the STRONGEST wall-to-wall chain bound available from a chain of T near-axis squares:
    B = [span cos t0 + R sin t0]/(T-1) - 1 with R the signed rise, minimised over chains
    (`search/arch_chain.py`).  -> (nchains, best B, best rise, axis, path) or (0, None, ...)."""
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    gx, dxr, gy, dyr = normal_gaps(z, n)
    best = None
    tot = 0
    for (g, dr, coord, perp, ax) in ((gx, dxr, X, Y, 'x'), (gy, dyr, Y, X, 'y')):
        for ch in all_chains(list(near), g, dr, coord, lvl - tol, int(round(T))):
            tot += 1
            # the SIGNED tilt matters: for a y-chain the rise enters with -sin t, and a chain of
            # squares tilted the other way is a different (and much weaker) certificate.
            t0 = float(np.mean([s6skel.norm_tilt(TH[q]) for q in ch]))
            R = (perp[ch[-1]] - perp[ch[0]]) * (1.0 if ax == 'x' else -1.0)
            span = coord[ch[-1]] - coord[ch[0]]
            u = math.cos(t0) + abs(math.sin(t0))
            B = (span * math.cos(t0) + R * math.sin(t0)) / (T - 1) - 1
            Bw = ((T - u) * math.cos(t0) + R * math.sin(t0)) / (T - 1) - 1
            if best is None or B < best[0]:
                best = (B, Bw, R, span, ax, ch, math.degrees(t0))
    if best is None:
        return dict(nchains=0)
    B, Bw, R, span, ax, ch, t0 = best
    return dict(nchains=tot, B=B, B_walltowall=Bw, rise=R, span=span, axis=ax, path=ch,
                tilt_deg=math.degrees(t0))


def cmd_chains(a):
    recs = []
    for f in a.file:
        recs += [json.loads(l) for l in open(f)]
    cells = {}
    for r in recs:
        if r['value'] is None:
            continue
        key = (r['k'], r['eps'], r['pidx'] >= 0)
        if key not in cells or r['value'] > cells[key]['value']:
            cells[key] = r
    print("# for the best configuration in each cell: chains of T among the NEAR-AXIS squares,")
    print("# at the level of the configuration's own margin, and the strongest chain bound B_T.")
    print("#  k eps   run         delta*        nchains   minrise   B_T (this chain)   "
          "B_T (wall-to-wall)")
    for key in sorted(cells):
        r = cells[key]
        z = np.array(r['z'])
        n = (len(z) - 1) // 3
        c = chain_certificate(z, n, range(r['k']), r['value'], a.T)
        tag = 'chainfree' if key[2] else 'free     '
        if c['nchains'] == 0:
            print(f"  {r['k']:2d} {r['eps']:<5g} {tag} {r['value']:+.6e}   0        "
                  f"--        --                 --")
        else:
            print(f"  {r['k']:2d} {r['eps']:<5g} {tag} {r['value']:+.6e}   {c['nchains']:<4d}     "
                  f"{c['rise']:+.5f}  {c['B']:+.6e}      {c['B_walltowall']:+.6e}  "
                  f"{c['axis']}-chain {c['path']}")


def cmd_cert(a):
    """the dual of the (chain-constrained) LP in the centres at the recorded optimum: which rows
    hold it shut, and do any of them belong to a TILTED square?"""
    recs = [json.loads(l) for l in open(a.file)]
    recs = [r for r in recs if r['value'] is not None and r['k'] == a.k and r['eps'] == a.eps
            and ((r['pidx'] >= 0) == (not a.free))]
    recs.sort(key=lambda r: -r['value'])
    for r in recs[:a.top]:
        z = np.array(r['z'])
        n = (len(z) - 1) // 3
        eps = math.radians(r['eps'])
        pat = None if r['pattern'] is None else [tuple(c) for c in r['pattern']]
        F = FixedBand(list(z[3::3]), a.T, r['k'], pat, a.eta)
        asg = F.assign(z[1::3], z[2::3])
        x, res, tags = F.lp(asg, full=True)
        if x is None:
            print('# LP failed'); continue
        w = -np.array(res.ineqlin.marginals)
        tl = [math.degrees(tilt(t)) for t in z[3::3]]
        print(f"\n# k={r['k']} eps={r['eps']} pattern={r['pattern']} value={r['value']:+.9e} "
              f"LP={x[0]:+.9e}")
        print(f"#   tilts (deg): " + ' '.join(f"{q}:{tl[q]:.3f}" for q in range(n)))
        tot = {'wall': 0.0, 'pair': 0.0, 'chain': 0.0}
        touched = set()
        for q, tg in enumerate(tags):
            if w[q] <= a.wtol:
                continue
            kind = 'chain' if tg[0] == 'chain' else ('pair' if tg[0] == 'pair' else 'wall')
            tot[kind] += float(w[q])
            sq = [tg[2]] if kind == 'wall' else ([tg[1], tg[2]] if kind == 'pair'
                                                 else [tg[2], tg[3]])
            touched |= set(sq)
            mx = max(tl[s] for s in sq)
            print(f"    w={w[q]:9.6f}  {str(tg):<30s} max tilt {mx:7.3f} deg"
                  + ('   <-- FAR square' if max(sq) >= r['k'] else ''))
        near_only = all(s < r['k'] for s in touched)
        print(f"#   support weight: wall {tot['wall']:.4f}  pair {tot['pair']:.4f}  "
              f"chain {tot['chain']:.4f};  squares in the support: {sorted(touched)};  "
              f"{'ALL NEAR-AXIS' if near_only else 'INCLUDES A TILTED SQUARE'}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('scan')
    p.add_argument('--n', type=int, required=True)
    p.add_argument('--T', type=float, required=True)
    p.add_argument('--eps', default='1,5,10', help='degrees')
    p.add_argument('--ks', default=None)
    p.add_argument('--mode', default='both', choices=('both', 'free', 'chainfree'))
    p.add_argument('--starts', type=int, default=800)
    p.add_argument('--pstarts', type=int, default=300)
    p.add_argument('--eta', type=float, default=1e-6)
    p.add_argument('--tol-edge', dest='tol_edge', type=float, default=1e-7)
    p.add_argument('--patmink', type=int, default=0)
    p.add_argument('--polish', type=int, default=2)
    p.add_argument('--cycles', type=int, default=2)
    p.add_argument('--rounds', type=int, default=8)
    p.add_argument('--maxiter', type=int, default=200)
    p.add_argument('--nproc', type=int, default=10)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_scan)
    p = sub.add_parser('report')
    p.add_argument('file', nargs='+')
    p.set_defaults(fn=cmd_report)
    p = sub.add_parser('show')
    p.add_argument('file')
    p.add_argument('--k', type=int, required=True)
    p.add_argument('--eps', type=float, required=True)
    p.add_argument('--free', action='store_true')
    p.add_argument('--top', type=int, default=1)
    p.set_defaults(fn=cmd_show)
    p = sub.add_parser('chains')
    p.add_argument('file', nargs='+')
    p.add_argument('--T', type=float, required=True)
    p.set_defaults(fn=cmd_chains)
    p = sub.add_parser('cert')
    p.add_argument('file')
    p.add_argument('--T', type=float, required=True)
    p.add_argument('--k', type=int, required=True)
    p.add_argument('--eps', type=float, required=True)
    p.add_argument('--free', action='store_true')
    p.add_argument('--top', type=int, default=1)
    p.add_argument('--eta', type=float, default=1e-6)
    p.add_argument('--wtol', type=float, default=1e-7)
    p.set_defaults(fn=cmd_cert)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
