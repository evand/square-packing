#!/usr/bin/env python3
"""rank8.py -- the margin of the rank-8 statement, and the minimal core of the missing unit.

Task: `tasks/rank8-margin/README.md`.  Note: `search/RANK8.md`.

Setting (`search/BENTZ.md` §0, §5; `notes/s13-casefree.md` §1): `t = 4`, closed unit squares,
closed containment, a packing is a family pairwise disjoint **as closed sets**.  `P0` is Bentz's
sixteen points (`search/bentz.py points`).  **Leaf A** is the fully pinned pattern leaf `(AB)^4`:
twelve *labelled* squares,

    T_i   (i = 0..3)  pattern exactly {a_i, b_i}      -- the four corner squares
    C_j   (j = 0..3)  pattern exactly {c_j}           -- the four mid-wall singletons
    D_j   (j = 0..3)  pattern exactly {d_j}           -- the four interior singletons

each square confined to its pattern region: it must contain its own points and must **not**
contain any of the other sixteen.

The measurement is the sup over leaf A of

    F(S_1..S_12) = min_{i<j} gap(S_i, S_j)

with `gap` the separating-axis signed distance.  For two closed unit squares at
`(c_i, th_i)`, `(c_j, th_j)`, with `Dc = c_j - c_i` and `Delta = th_j - th_i`,

    gap = max( |Dc . u_i|, |Dc . v_i|, |Dc . u_j|, |Dc . v_j| ) - W(Delta),
    W(Delta) = 1/2 + 1/2 (|cos Delta| + |sin Delta|),

because the support half-width of a unit square in ANY direction `d` is
`1/2 (|d.u| + |d.v|)`, which is `1/2` on its own two axes and `W - 1/2` on the other square's.
`gap > 0` iff the two closed squares are disjoint; `gap = 0` iff they touch; `gap < 0` iff they
overlap.  (SAT on the four edge normals is exact for two convex quadrilaterals.)

**Everything in this file is floating point.**  It measures a margin; it certifies nothing.  The
one exact ingredient is the re-derivation of a candidate's pattern and of its pairwise
disjointness through `leaf_ceiling.sq_contains` / `sq_meets_sq` after snapping to rationals
(`verify` sub-command), which turns a *positive* answer into a theorem about that configuration.
A *negative* answer (`sup <= 0`) is never certified here: it is a multistart local-optimisation
measurement.

Sub-commands
------------
  margin   : (1) sup of the min pairwise closed gap over leaf A, by multistart local optimisation
  subsets  : (2) which subsets of the eight singletons are jointly realisable with the four
             corner squares; the minimal infeasible ones up to D4
  linear   : (3) the first-order (linearised) disjointness system around the tiling, as an LP,
             plus the validity radius of the linearisation
  sa2      : (4) Sherali-Adams level 2 on the 162-pose support of a measure file
  verify   : exact re-check of a saved configuration (patterns + pairwise disjointness)
  selftest : the geometry primitives against `leaf_ceiling`'s exact ones
"""
import argparse
import itertools
import json
import math
import os
import sys
from fractions import Fraction as Fr

# tiny dense linear algebra: BLAS threading is pure overhead here, and the multistart is
# parallelised over processes instead (at most 8, `notes/status.md`).
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_ceiling as lc                                             # noqa: E402
import bentz                                                         # noqa: E402

T = 4.0

# ===================================================================== P0 and the twelve labels
_PTS = bentz.bentz_points()
PNAMES = [n for (n, _x, _y) in _PTS]
PXY = np.array([[float(x), float(y)] for (_n, x, y) in _PTS])
PIDX = {n: k for k, n in enumerate(PNAMES)}
IPTS = bentz.as_int_points(_PTS)

# label i -> (name, the point indices its pattern must contain)
LABELS = ([(f'T{i}', (PIDX[f'a{i}'], PIDX[f'b{i}'])) for i in range(4)]
          + [(f'C{j}', (PIDX[f'c{j}'],)) for j in range(4)]
          + [(f'D{j}', (PIDX[f'd{j}'],)) for j in range(4)])
LNAME = [n for (n, _s) in LABELS]
LIDX = {n: k for k, n in enumerate(LNAME)}
SINGLETONS = list(range(4, 12))          # the eight non-corner labels C0..C3, D0..D3

# the D4 action on the twelve labels, induced by its action on P0
def label_perms():
    perms = bentz.d4_perm(_PTS)
    out = []
    for pm in perms:
        img = []
        for (_n, pts) in LABELS:
            q = frozenset(pm[p] for p in pts)
            hit = [k for k, (_m, s) in enumerate(LABELS) if frozenset(s) == q]
            assert len(hit) == 1
            img.append(hit[0])
        out.append(tuple(img))
    return out


LPERM = label_perms()


# ===================================================================== float geometry
def axes(th):
    c, s = np.cos(th), np.sin(th)
    return c, s


def margins(z, labels):
    """(n, 16, 2) array: the two signed coordinates of each P0 point in each square's own frame.
    The point is in the closed square iff max(|.|) <= 1/2."""
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    c, s = np.cos(TH), np.sin(TH)
    dx = PXY[None, :, 0] - X[:, None]
    dy = PXY[None, :, 1] - Y[:, None]
    u = c[:, None] * dx + s[:, None] * dy
    v = -s[:, None] * dx + c[:, None] * dy
    return np.stack([u, v], axis=2)


def patterns_float(z, tol=0.0):
    m = np.max(np.abs(margins(z, None)), axis=2)
    return [frozenset(np.nonzero(row <= 0.5 + tol)[0].tolist()) for row in m]


def pair_gap(z, i, j):
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    dcx, dcy = X[j] - X[i], Y[j] - Y[i]
    best = -1e18
    for o in (i, j):
        c, s = math.cos(TH[o]), math.sin(TH[o])
        for (d0, d1) in ((c, s), (-s, c)):
            best = max(best, abs(d0 * dcx + d1 * dcy))
    d = TH[j] - TH[i]
    return best - 0.5 - 0.5 * (abs(math.cos(d)) + abs(math.sin(d)))


def all_gaps(z, n):
    return {(i, j): pair_gap(z, i, j) for i in range(n) for j in range(i + 1, n)}


def true_objective(z, n):
    g = all_gaps(z, n)
    return (min(g.values()), min(g, key=g.get)) if g else (float('inf'), None)


def admissibility_slack(z, n):
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    w = 0.5 * (np.abs(np.cos(TH)) + np.abs(np.sin(TH)))
    return float(min(np.min(X - w), np.min(T - X - w), np.min(Y - w), np.min(T - Y - w)))


def pattern_slack(z, labels):
    """(worst containment margin over required points, worst exclusion margin over the others).
    Containment wants >= 0, exclusion wants < 0 (we report the max, which wants < 0)."""
    m = np.max(np.abs(margins(z, labels)), axis=2)
    inc, exc = [], []
    for k, li in enumerate(labels):
        want = set(LABELS[li][1])
        for p in range(16):
            (inc if p in want else exc).append(0.5 - m[k, p])
    return min(inc), max(exc)


# ===================================================================== the smooth sub-problem
SGN = ((1, 1), (1, -1), (-1, 1), (-1, -1))


class Problem:
    """leaf A (or a sub-configuration): `labels` is a list of label ids, one square each."""

    def __init__(self, labels, eps_excl=0.0):
        self.labels = list(labels)
        self.n = len(self.labels)
        self.eps = eps_excl
        self.pairs = [(i, j) for i in range(self.n) for j in range(i + 1, self.n)]
        self.inc = []                     # (square, point) that must be contained
        self.exc = []                     # (square, point) that must not be
        for k, li in enumerate(self.labels):
            want = set(LABELS[li][1])
            for p in range(16):
                (self.inc if p in want else self.exc).append((k, p))
        self.nv = 1 + 3 * self.n

    # ---------- assignments (which axis separates a pair; which axis excludes a point) ----------
    def pair_assign(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        out = []
        for (i, j) in self.pairs:
            dcx, dcy = X[j] - X[i], Y[j] - Y[i]
            best, arg = -1e18, None
            for oi, o in enumerate((i, j)):
                c, s = math.cos(TH[o]), math.sin(TH[o])
                for kind, (d0, d1) in enumerate(((c, s), (-s, c))):
                    pr = d0 * dcx + d1 * dcy
                    for sg in (1, -1):
                        if sg * pr > best:
                            best, arg = sg * pr, (oi, kind, sg)
            out.append(arg)
        return out

    def exc_assign(self, z):
        m = margins(z, self.labels)
        out = []
        for (k, p) in self.exc:
            best, arg = -1e18, None
            for kind in (0, 1):
                for sg in (1, -1):
                    val = sg * m[k, p, kind]
                    if val > best:
                        best, arg = val, (kind, sg)
            out.append(arg)
        return out

    # ---------------------------------- constraint assembly ----------------------------------
    def build(self, pa, ea):
        """precompute the index arrays for the fixed-assignment smooth NLP"""
        P = []
        for (pi, (i, j)) in enumerate(self.pairs):
            oi, kind, sg = pa[pi]
            o = i if oi == 0 else j
            for (s1, s2) in SGN:
                P.append((i, j, o, kind, sg, s1, s2))
        self.P = np.array(P, dtype=float) if P else np.zeros((0, 7))
        self.Pi = np.array([p[0] for p in P], dtype=int)
        self.Pj = np.array([p[1] for p in P], dtype=int)
        self.Po = np.array([p[2] for p in P], dtype=int)
        self.Pk = np.array([p[3] for p in P], dtype=int)
        self.Ps = np.array([p[4] for p in P], dtype=float)
        self.P1 = np.array([p[5] for p in P], dtype=float)
        self.P2 = np.array([p[6] for p in P], dtype=float)

        I = [(k, p, kind, sg) for (k, p) in self.inc for kind in (0, 1) for sg in (1, -1)]
        self.Ik = np.array([a[0] for a in I], dtype=int)
        self.Ip = np.array([a[1] for a in I], dtype=int)
        self.Ii = np.array([a[2] for a in I], dtype=int)
        self.Is = np.array([a[3] for a in I], dtype=float)

        self.Ek = np.array([k for (k, _p) in self.exc], dtype=int)
        self.Ep = np.array([p for (_k, p) in self.exc], dtype=int)
        self.Ei = np.array([a[0] for a in ea], dtype=int)
        self.Es = np.array([a[1] for a in ea], dtype=float)

        A = [(k, side, s1, s2) for k in range(self.n) for side in range(4) for (s1, s2) in SGN]
        self.Ak = np.array([a[0] for a in A], dtype=int)
        self.Aside = np.array([a[1] for a in A], dtype=int)
        self.A1 = np.array([a[2] for a in A], dtype=float)
        self.A2 = np.array([a[3] for a in A], dtype=float)
        self.m = len(P) + len(I) + len(self.exc) + len(A)

    def cons(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        g = z[0]
        C, S = np.cos(TH), np.sin(TH)
        out = []
        # pairs
        if len(self.Pi):
            co, so = C[self.Po], S[self.Po]
            d0 = np.where(self.Pk == 0, co, -so)
            d1 = np.where(self.Pk == 0, so, co)
            dcx = X[self.Pj] - X[self.Pi]
            dcy = Y[self.Pj] - Y[self.Pi]
            dth = TH[self.Pj] - TH[self.Pi]
            out.append(self.Ps * (d0 * dcx + d1 * dcy) - 0.5
                       - 0.5 * (self.P1 * np.cos(dth) + self.P2 * np.sin(dth)) - g)
        # containment
        if len(self.Ik):
            c, s = C[self.Ik], S[self.Ik]
            dx = PXY[self.Ip, 0] - X[self.Ik]
            dy = PXY[self.Ip, 1] - Y[self.Ik]
            e0 = np.where(self.Ii == 0, c, -s)
            e1 = np.where(self.Ii == 0, s, c)
            out.append(0.5 - self.Is * (e0 * dx + e1 * dy))
        # exclusion
        if len(self.Ek):
            c, s = C[self.Ek], S[self.Ek]
            dx = PXY[self.Ep, 0] - X[self.Ek]
            dy = PXY[self.Ep, 1] - Y[self.Ek]
            e0 = np.where(self.Ei == 0, c, -s)
            e1 = np.where(self.Ei == 0, s, c)
            out.append(self.Es * (e0 * dx + e1 * dy) - 0.5 - self.eps)
        # admissibility
        c, s = C[self.Ak], S[self.Ak]
        h = 0.5 * (self.A1 * c + self.A2 * s)
        base = np.where(self.Aside == 0, X[self.Ak],
                        np.where(self.Aside == 1, T - X[self.Ak],
                                 np.where(self.Aside == 2, Y[self.Ak], T - Y[self.Ak])))
        out.append(base - h)
        return np.concatenate(out)

    def jac(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        C, S = np.cos(TH), np.sin(TH)
        J = np.zeros((self.m, self.nv))
        r = 0
        if len(self.Pi):
            k = len(self.Pi)
            co, so = C[self.Po], S[self.Po]
            d0 = np.where(self.Pk == 0, co, -so)
            d1 = np.where(self.Pk == 0, so, co)
            dp0 = np.where(self.Pk == 0, -so, -co)     # d/dth of d0
            dp1 = np.where(self.Pk == 0, co, -so)      # d/dth of d1
            dcx = X[self.Pj] - X[self.Pi]
            dcy = Y[self.Pj] - Y[self.Pi]
            dth = TH[self.Pj] - TH[self.Pi]
            rows = np.arange(k)
            J[rows, 0] = -1.0
            np.add.at(J, (rows, 1 + 3 * self.Pj), self.Ps * d0)
            np.add.at(J, (rows, 2 + 3 * self.Pj), self.Ps * d1)
            np.add.at(J, (rows, 1 + 3 * self.Pi), -self.Ps * d0)
            np.add.at(J, (rows, 2 + 3 * self.Pi), -self.Ps * d1)
            np.add.at(J, (rows, 3 + 3 * self.Po), self.Ps * (dp0 * dcx + dp1 * dcy))
            w = -0.5 * (-self.P1 * np.sin(dth) + self.P2 * np.cos(dth))
            np.add.at(J, (rows, 3 + 3 * self.Pj), w)
            np.add.at(J, (rows, 3 + 3 * self.Pi), -w)
            r += k
        if len(self.Ik):
            k = len(self.Ik)
            c, s = C[self.Ik], S[self.Ik]
            dx = PXY[self.Ip, 0] - X[self.Ik]
            dy = PXY[self.Ip, 1] - Y[self.Ik]
            e0 = np.where(self.Ii == 0, c, -s)
            e1 = np.where(self.Ii == 0, s, c)
            ep0 = np.where(self.Ii == 0, -s, -c)
            ep1 = np.where(self.Ii == 0, c, -s)
            rows = r + np.arange(k)
            np.add.at(J, (rows, 1 + 3 * self.Ik), self.Is * e0)
            np.add.at(J, (rows, 2 + 3 * self.Ik), self.Is * e1)
            np.add.at(J, (rows, 3 + 3 * self.Ik), -self.Is * (ep0 * dx + ep1 * dy))
            r += k
        if len(self.Ek):
            k = len(self.Ek)
            c, s = C[self.Ek], S[self.Ek]
            dx = PXY[self.Ep, 0] - X[self.Ek]
            dy = PXY[self.Ep, 1] - Y[self.Ek]
            e0 = np.where(self.Ei == 0, c, -s)
            e1 = np.where(self.Ei == 0, s, c)
            ep0 = np.where(self.Ei == 0, -s, -c)
            ep1 = np.where(self.Ei == 0, c, -s)
            rows = r + np.arange(k)
            np.add.at(J, (rows, 1 + 3 * self.Ek), -self.Es * e0)
            np.add.at(J, (rows, 2 + 3 * self.Ek), -self.Es * e1)
            np.add.at(J, (rows, 3 + 3 * self.Ek), self.Es * (ep0 * dx + ep1 * dy))
            r += k
        k = len(self.Ak)
        c, s = C[self.Ak], S[self.Ak]
        rows = r + np.arange(k)
        sx = np.where(self.Aside == 0, 1.0, np.where(self.Aside == 1, -1.0, 0.0))
        sy = np.where(self.Aside == 2, 1.0, np.where(self.Aside == 3, -1.0, 0.0))
        np.add.at(J, (rows, 1 + 3 * self.Ak), sx)
        np.add.at(J, (rows, 2 + 3 * self.Ak), sy)
        np.add.at(J, (rows, 3 + 3 * self.Ak), -0.5 * (-self.A1 * s + self.A2 * c))
        return J


def solve(labels, z0, eps_excl=0.0, rounds=8, maxiter=300):
    """max-min-gap by assignment-fixed SLSQP, iterated until the assignment stops changing.
    Returns (true objective at the returned point, z, info)."""
    from scipy.optimize import minimize
    pr = Problem(labels, eps_excl)
    z = np.array(z0, dtype=float)
    z[0] = true_objective(z, pr.n)[0]
    best = (true_objective(z, pr.n)[0], z.copy())
    seen = set()
    for it in range(rounds):
        pa, ea = pr.pair_assign(z), pr.exc_assign(z)
        key = (tuple(pa), tuple(ea))
        pr.build(pa, ea)
        res = minimize(lambda w: -w[0], z, jac=lambda w: np.eye(1, len(w), 0).ravel() * -1.0,
                       constraints=[{'type': 'ineq', 'fun': pr.cons, 'jac': pr.jac}],
                       method='SLSQP', options={'maxiter': maxiter, 'ftol': 1e-12})
        zn = res.x
        val = _feasible_value(zn, pr)
        if val > best[0]:
            best = (val, zn.copy())
        if key in seen:
            break
        seen.add(key)
        z = zn
        z[0] = true_objective(z, pr.n)[0]
    return best[0], best[1], pr


def _feasible_value(z, pr):
    """the honest value of a candidate: -inf unless it really is in the (relaxed) leaf."""
    inc, exc = pattern_slack(z, pr.labels)
    adm = admissibility_slack(z, pr.n)
    tol = 1e-9
    if inc < -tol or exc > pr.eps + tol or adm < -tol:
        return -float('inf')
    return true_objective(z, pr.n)[0]


# ===================================================================== starts
def tiling_starts():
    """the `4x4` tiling minus four wall squares: four corner tiles, four interior tiles, and one
    wall tile per wall nudged along the wall so its pattern collapses to the mid-wall singleton.
    Each wall has two tiles and either can be used, so there are 2^4 = 16 such families."""
    tile = {}
    for i in range(4):
        for j in range(4):
            tile[(i, j)] = (i + 0.5, j + 0.5)
    corner = {0: (0, 0), 1: (3, 0), 2: (0, 3), 3: (3, 3)}
    #  wall j -> the two tiles that can hold c_j, and the nudge that keeps only c_j
    #  c0 = (0.914, 2) left wall;  c1 = (2, 0.914) bottom;  c2 = (2, 3.086) top; c3 = (3.086, 2) right
    wall = {0: [((0, 1), (0.0, +1.0)), ((0, 2), (0.0, -1.0))],
            1: [((1, 0), (+1.0, 0.0)), ((2, 0), (-1.0, 0.0))],
            2: [((1, 3), (+1.0, 0.0)), ((2, 3), (-1.0, 0.0))],
            3: [((3, 1), (0.0, +1.0)), ((3, 2), (0.0, -1.0))]}
    interior = {0: (1, 1), 1: (1, 2), 2: (2, 1), 3: (2, 2)}
    out = []
    for pick in itertools.product((0, 1), repeat=4):
        used = set()
        z = np.zeros(1 + 36)
        ok = True
        for i in range(4):
            c = tile[corner[i]]
            z[1 + 3 * i], z[2 + 3 * i], z[3 + 3 * i] = c[0], c[1], 0.0
            used.add(corner[i])
        for j in range(4):
            t, nd = wall[j][pick[j]]
            if t in used:
                ok = False
            used.add(t)
            k = 4 + j
            z[1 + 3 * k] = tile[t][0] + 1e-3 * nd[0]
            z[2 + 3 * k] = tile[t][1] + 1e-3 * nd[1]
            z[3 + 3 * k] = 0.0
        for j in range(4):
            k = 8 + j
            z[1 + 3 * k], z[2 + 3 * k], z[3 + 3 * k] = tile[interior[j]][0], tile[interior[j]][1], 0.0
        if ok:
            out.append((pick, z))
    return out


def sample_pools(nsamp, rng, tol=0.0):
    """rejection-sample poses whose exact float pattern is one of the twelve labels"""
    pools = [[] for _ in range(12)]
    want = {frozenset(s): k for k, (_n, s) in enumerate(LABELS)}
    B = 200000
    done = 0
    while done < nsamp:
        b = min(B, nsamp - done)
        done += b
        th = rng.uniform(0, math.pi / 2, b)
        w = 0.5 * (np.abs(np.cos(th)) + np.abs(np.sin(th)))
        x = rng.uniform(0, 1, b) * (T - 2 * w) + w
        y = rng.uniform(0, 1, b) * (T - 2 * w) + w
        c, s = np.cos(th), np.sin(th)
        dx = PXY[None, :, 0] - x[:, None]
        dy = PXY[None, :, 1] - y[:, None]
        u = np.abs(c[:, None] * dx + s[:, None] * dy)
        v = np.abs(-s[:, None] * dx + c[:, None] * dy)
        inside = np.maximum(u, v) <= 0.5 + tol
        for r in np.nonzero(inside.sum(1) <= 2)[0]:
            f = frozenset(np.nonzero(inside[r])[0].tolist())
            k = want.get(f)
            if k is not None:
                pools[k].append((x[r], y[r], th[r]))
    return pools


def support_pools(path):
    """the 162 support poses of a measure file, grouped by which label they realise"""
    t, sym, poses = lc.read_measure(path)
    pools = [[] for _ in range(12)]
    want = {frozenset(s): k for k, (_n, s) in enumerate(LABELS)}
    for (p, q, cx, cy, m) in poses:
        pat = bentz.pattern_of(lc.make_square(cx, cy, p, q, t), IPTS)
        k = want.get(pat)
        if k is not None:
            pools[k].append((float(cx), float(cy), 2 * math.atan2(p, q)))
    return pools


def z_from(poses):
    z = np.zeros(1 + 3 * len(poses))
    for k, (x, y, th) in enumerate(poses):
        z[1 + 3 * k], z[2 + 3 * k], z[3 + 3 * k] = x, y, th
    return z


# ===================================================================== (1) the margin
_JOB = {}


def _job(arg):
    tag, z0 = arg
    val, z, _pr = solve(_JOB['labels'], z0, eps_excl=_JOB['eps'], rounds=_JOB['rounds'],
                        maxiter=_JOB['maxiter'])
    return (val, tag, z)


def _init(labels, eps, rounds, maxiter):
    _JOB.update(labels=labels, eps=eps, rounds=rounds, maxiter=maxiter)


def run_starts(starts, labels, a, nproc=8):
    import multiprocessing as mp
    if nproc <= 1:
        _init(labels, a.eps, a.rounds, a.maxiter)
        return [_job(s) for s in starts]
    with mp.Pool(nproc, initializer=_init,
                 initargs=(labels, a.eps, a.rounds, a.maxiter)) as pool:
        return pool.map(_job, starts, chunksize=4)


def cmd_margin(a):
    rng = np.random.default_rng(a.seed)
    labels = list(range(12))
    starts = []

    base = tiling_starts()
    print(f'# {len(base)} tiling-minus-4 families')
    for (pick, z) in base:
        starts.append((f'tiling{"".join(map(str, pick))}', z))

    print(f'# {a.nudge} nudged-tiling starts (sigma {a.sigma} centre, {a.sigma_th} rad angle)')
    for k in range(a.nudge):
        pick, z = base[rng.integers(len(base))]
        z = z.copy()
        z[1::3] += rng.normal(0, a.sigma, 12)
        z[2::3] += rng.normal(0, a.sigma, 12)
        z[3::3] += rng.normal(0, a.sigma_th, 12)
        starts.append((f'nudge{k}', z))

    sp = support_pools(a.support) if a.support else None
    if sp:
        print(f'# {a.nsupport} starts from the 162-pose support of the certified measure')
        print('  support poses per label: '
              + ', '.join(f'{LNAME[k]}:{len(sp[k])}' for k in range(12)))
        for k in range(a.nsupport):
            z = z_from([sp[i][rng.integers(len(sp[i]))] if sp[i] else (2.0, 2.0, 0.0)
                        for i in range(12)])
            starts.append((f'support{k}', z))

    print(f'# {a.nrandom} random starts from rejection-sampled pattern pools')
    pools = sample_pools(a.nsamp, rng)
    print('  pool sizes: ' + ', '.join(f'{LNAME[k]}:{len(pools[k])}' for k in range(12)))
    for k in range(a.nrandom):
        z = z_from([pools[i][rng.integers(len(pools[i]))] for i in range(12)])
        starts.append((f'rand{k}', z))

    print(f'# running {len(starts)} local optimisations on {a.nproc} processes', flush=True)
    results = run_starts(starts, labels, a, a.nproc)
    results.sort(key=lambda r: -r[0])
    fin = [r for r in results if r[0] > -float('inf')]
    print(f'\n{len(fin)} of {len(results)} starts ended inside the (relaxed) leaf')
    print(f'BEST min pairwise closed gap over leaf A: {results[0][0]:+.12e}   ({results[0][1]})')
    print('top 15:')
    for (v, tag, z) in results[:15]:
        _, exc = pattern_slack(z, labels)
        print(f'  {v:+.12e}  {tag:<14} max tilt {max_tilt_deg(z):8.4f} deg, '
              f'dist to tiling {_tiling_distance(z):8.6f}, worst exclusion {exc:+.2e}')
    near = [(v, tag, z) for (v, tag, z) in fin if v > -a.near]
    far = [(v, tag, z) for (v, tag, z) in near if _tiling_distance(z) > a.far]
    print(f'\n{len(near)} starts within {a.near} of 0; of those, {len(far)} are further than '
          f'{a.far} (max-norm, centres and angles mod 90 deg) from any tiling-minus-4 family')
    for (v, tag, z) in sorted(far, key=lambda r: -r[0])[:12]:
        _, exc = pattern_slack(z, labels)
        print(f'  {v:+.12e}  {tag:<14} max tilt {max_tilt_deg(z):8.4f} deg, '
              f'dist to tiling {_tiling_distance(z):8.6f}, worst exclusion {exc:+.2e}')
    tilted = [(v, tag, z) for (v, tag, z) in near if max_tilt_deg(z) > a.tilt_deg]
    print(f'\n{len(tilted)} of the near-0 configurations have a square tilted by more than '
          f'{a.tilt_deg} deg (mod 90):')
    for (v, tag, z) in sorted(tilted, key=lambda r: -r[0])[:12]:
        print(f'  {v:+.12e}  {tag:<14} max tilt {max_tilt_deg(z):8.4f} deg, '
              f'dist to tiling {_tiling_distance(z):8.6f}')
    strict = [r for r in near if pattern_slack(r[2], labels)[1] < -1e-12]
    print(f'\n{len(strict)} of the {len(near)} near-0 configurations satisfy the exclusions '
          f'STRICTLY (are in leaf A itself, not just its closed relaxation); the best of those is '
          f'{max([r[0] for r in strict], default=float("nan")):+.12e}')
    if a.out:
        save(a.out, results[:a.nsave], labels)
        print(f'\nwrote {min(len(results), a.nsave)} configurations to {a.out}')
    return results


def norm_theta(th):
    """a unit square is invariant under a quarter turn, so the pose angle only matters mod 90 deg;
    normalise into (-45, 45] degrees so that "tilt" means what it says"""
    t = np.remainder(np.asarray(th) + math.pi / 4, math.pi / 2) - math.pi / 4
    return t


def max_tilt_deg(z):
    return float(np.max(np.abs(norm_theta(z[3::3])))) * 180.0 / math.pi


_TILING_CACHE = None


def _tiling_distance(z, strict=False):
    """max-norm distance (centres and angles mod 90 deg) to the nearest tiling-minus-4 family.
    With `strict`, the four interior squares must match too; otherwise they are compared as well
    (they always are -- the flag is kept for symmetry with the report)."""
    global _TILING_CACHE
    if _TILING_CACHE is None:
        _TILING_CACHE = tiling_starts()
    best = 1e18
    th = norm_theta(z[3::3])
    for (_pick, t) in _TILING_CACHE:
        d = max(float(np.max(np.abs(z[1::3] - t[1::3]))),
                float(np.max(np.abs(z[2::3] - t[2::3]))),
                float(np.max(np.abs(th - norm_theta(t[3::3])))))
        best = min(best, d)
    return best


def save(path, results, labels):
    with open(path, 'w') as f:
        f.write('# rank8.py configurations: label cx cy theta_rad, then the value.  FLOAT.\n')
        for (v, tag, z) in results:
            f.write(f'config {tag} value {v:.17g} n {len(labels)}\n')
            for k, li in enumerate(labels):
                f.write(f'  {LNAME[li]} {z[1+3*k]:.17g} {z[2+3*k]:.17g} {z[3+3*k]:.17g}\n')


def load(path):
    out, cur = [], None
    for line in open(path):
        s = line.split()
        if not s or s[0] == '#':
            continue
        if s[0] == 'config':
            cur = (s[1], float(s[3]), [])
            out.append(cur)
        else:
            cur[2].append((s[0], float(s[1]), float(s[2]), float(s[3])))
    return out


# ===================================================================== (2) sub-configurations
def cmd_subsets(a):
    rng = np.random.default_rng(a.seed)
    pools = sample_pools(a.nsamp, rng)
    sp = support_pools(a.support) if a.support else None
    base = tiling_starts()

    def feasible(labels, nstart):
        """multistart: the best min-gap found, and the configuration"""
        best = (-float('inf'), None)
        starts = []
        for (_pick, z) in base:
            starts.append(np.concatenate([[0.0]] + [z[1 + 3 * li:4 + 3 * li] for li in labels]))
        for _ in range(nstart):
            src = pools if (sp is None or rng.random() < 0.7) else sp
            starts.append(z_from([src[li][rng.integers(len(src[li]))] for li in labels]))
        for z0 in starts:
            v, z, _pr = solve(labels, z0, eps_excl=a.eps, rounds=a.rounds, maxiter=a.maxiter)
            if v > best[0]:
                best = (v, z.copy())
            if best[0] > a.pos:
                break
        return best

    # D4 orbits of the 256 subsets of the eight singletons
    orb, seen = [], set()
    for mask in range(256):
        if mask in seen:
            continue
        o = set()
        for pm in LPERM:
            o.add(sum(1 << (pm[SINGLETONS[b]] - 4) for b in range(8) if mask >> b & 1))
        seen |= o
        orb.append((min(o), len(o)))
    orb.sort(key=lambda kv: (bin(kv[0]).count('1'), kv[0]))
    print(f'{len(orb)} D4 classes of subsets of the eight singletons')

    status = {}          # mask -> (value, z) ; feasible iff value > 0
    canon = {}
    for (rep, _sz) in orb:
        for pm in LPERM:
            canon[sum(1 << (pm[SINGLETONS[b]] - 4) for b in range(8) if rep >> b & 1)] = rep

    def subs(mask):
        bits = [b for b in range(8) if mask >> b & 1]
        return [mask ^ (1 << b) for b in bits]

    # top-down: a subset of a feasible set is feasible, so only test what is not already implied
    order = sorted({r for (r, _s) in orb}, key=lambda m: -bin(m).count('1'))
    implied = {}
    for rep in order:
        if rep in implied:
            continue
        labels = list(range(4)) + [4 + b for b in range(8) if rep >> b & 1]
        ns = a.nstart if bin(rep).count('1') >= 6 else max(6, a.nstart // 3)
        v, z = feasible(labels, ns)
        status[rep] = (v, z)
        tag = 'FEASIBLE' if v > a.pos else ('infeasible' if v < a.pos else '?')
        print(f'  {_mask_name(rep):<28} k={bin(rep).count("1")}  best min-gap {v:+.6e}   {tag}',
              flush=True)
        if v > a.pos:
            # every subset is feasible too
            stack = [rep]
            while stack:
                m = stack.pop()
                for s in subs(m):
                    c = canon[s]
                    if c not in implied:
                        implied[c] = rep
                        stack.append(s)
    print(f'\n{len(status)} classes tested, {len(implied)} implied feasible by a larger '
          f'feasible class')
    feas = {m for m, (v, _z) in status.items() if v > a.pos} | set(implied)
    infeas = {m for m in canon.values() if m not in feas}
    minimal = [m for m in sorted(infeas, key=lambda m: bin(m).count('1'))
               if all(canon[s] in feas for s in subs(m))]
    print('\nMINIMAL INFEASIBLE subsets of the eight singletons (with all four corner squares), '
          'up to D4:')
    for m in minimal:
        print(f'  {_mask_name(m)}   |T| = {bin(m).count("1")}   '
              f'best min-gap {status.get(m, (float("nan"),))[0]:+.6e}')
    maximal = [m for m in sorted(feas, key=lambda m: -bin(m).count('1'))
               if all(canon[m | (1 << b)] not in feas for b in range(8) if not m >> b & 1)]
    print('\nMAXIMAL FEASIBLE subsets up to D4:')
    for m in maximal:
        print(f'  {_mask_name(m)}   |T| = {bin(m).count("1")}   '
              f'best min-gap {status.get(m, (float("nan"),))[0]:+.6e}')
    if a.out:
        keep = [(v, _mask_name(m), z) for m, (v, z) in status.items() if z is not None]
        with open(a.out, 'w') as f:
            f.write('# rank8.py subsets: label cx cy theta_rad.  FLOAT.\n')
            for (v, tag, z) in keep:
                labels = list(range(4)) + [4 + b for b in range(8)
                                           if _name_mask(tag) >> b & 1]
                f.write(f'config {tag} value {v:.17g} n {len(labels)}\n')
                for k, li in enumerate(labels):
                    f.write(f'  {LNAME[li]} {z[1+3*k]:.17g} {z[2+3*k]:.17g} {z[3+3*k]:.17g}\n')
        print(f'\nwrote {len(keep)} configurations to {a.out}')
    return status


def _mask_name(mask):
    return '{' + ','.join(LNAME[4 + b] for b in range(8) if mask >> b & 1) + '}'


def _name_mask(name):
    s = set(name.strip('{}').split(','))
    return sum(1 << b for b in range(8) if LNAME[4 + b] in s)


# ===================================================================== (3) the linearised system
def cmd_linear(a):
    from scipy.optimize import linprog
    pick = tuple(int(c) for c in a.pick)
    base = dict(tiling_starts())[pick]
    n = 12
    labels = list(range(12))
    X0, Y0 = base[1::3].copy(), base[2::3].copy()
    # the un-nudged base configuration: the tiling itself (nudges are what the LP is asked to find)
    X0 = np.round(X0 - 0.5) + 0.5
    Y0 = np.round(Y0 - 0.5) + 0.5
    print(f'linearisation about the tiling-minus-4 family {pick}:')
    for k in range(12):
        print(f'  {LNAME[labels[k]]:<3} tile centre ({X0[k]:.1f}, {Y0[k]:.1f})')
    # variables: g, then (dx, dy, phi) per square
    nv = 1 + 3 * n
    A, b, rows = [], [], []

    def add(coef, rhs, tag):
        A.append(coef)
        b.append(rhs)
        rows.append(tag)

    # --- pattern constraints, linearised.  A point p on the boundary of tile k gives an active
    #     linear constraint; a point strictly inside / outside gives a constraint with a constant
    #     slack, which we record as the validity radius rather than as a row.
    slack = []
    for k in range(n):
        want = set(LABELS[labels[k]][1])
        for p in range(16):
            dx, dy = PXY[p, 0] - X0[k], PXY[p, 1] - Y0[k]
            # |dx - ddx + phi*dy| <= 1/2  and  |dy - ddy - phi*dx| <= 1/2   (first order)
            for (val, gx, gy, gphi) in ((dx, -1.0, 0.0, dy), (dy, 0.0, -1.0, -dx)):
                for sg in (1.0, -1.0):
                    lhs0 = sg * val                     # nominal value of the signed coordinate
                    c = np.zeros(nv)
                    c[1 + 3 * k] = sg * gx
                    c[2 + 3 * k] = sg * gy
                    c[3 + 3 * k] = sg * gphi
                    if p in want:
                        if abs(lhs0 - 0.5) < 1e-12:     # active: sg*coord <= 1/2
                            add(c, 0.5 - lhs0, f'{LNAME[labels[k]]} holds {PNAMES[p]}')
                        else:
                            slack.append((0.5 - lhs0, f'{LNAME[labels[k]]} holds {PNAMES[p]}'))
            if p not in want:
                # NOT contained: at least one of the four signed coordinates exceeds 1/2.
                # Nominally exactly one can be active (a point on the boundary of the tile);
                # then that branch is forced and gives a strict linear row.
                act = []
                for (val, gx, gy, gphi) in ((dx, -1.0, 0.0, dy), (dy, 0.0, -1.0, -dx)):
                    for sg in (1.0, -1.0):
                        if abs(sg * val - 0.5) < 1e-12:
                            act.append((sg, gx, gy, gphi))
                if act:
                    if len(act) > 1:
                        print(f'  note: {LNAME[labels[k]]} excludes {PNAMES[p]} with '
                              f'{len(act)} tight coordinates (a corner incidence)')
                    sg, gx, gy, gphi = act[0]
                    c = np.zeros(nv)
                    c[1 + 3 * k] = -sg * gx
                    c[2 + 3 * k] = -sg * gy
                    c[3 + 3 * k] = -sg * gphi
                    add(c, 0.0, f'{LNAME[labels[k]]} excludes {PNAMES[p]}')
                else:
                    m = max(abs(dx), abs(dy)) - 0.5
                    slack.append((m, f'{LNAME[labels[k]]} excludes {PNAMES[p]}'))
    print(f'  {len(A)} active first-order pattern rows; smallest inactive slack '
          f'{min(s for s, _t in slack):.6f} ({min(slack)[1]})')

    # --- pairwise disjointness, linearised.  For tiles whose nominal centres differ by
    #     (1,0), (0,1) or (1,1) the nominal gap is 0 and the pair contributes a row; tiles two
    #     apart have nominal gap >= 1 and contribute only to the validity radius.
    disj = []
    for i in range(n):
        for j in range(i + 1, n):
            ddx, ddy = X0[j] - X0[i], Y0[j] - Y0[i]
            if max(abs(ddx), abs(ddy)) > 1.5:
                continue
            disj.append((i, j, ddx, ddy))
    print(f'  {len(disj)} adjacent pairs (nominal gap 0) of {n*(n-1)//2}')

    # branch over which axis separates each diagonal pair
    diag = [k for k, (_i, _j, ddx, ddy) in enumerate(disj) if abs(ddx) > 0.5 and abs(ddy) > 0.5]
    print(f'  {len(diag)} of them are diagonal (two axes tie at first order): '
          f'{2**len(diag)} LP branches')
    best = (-1e18, None, None)
    A0, b0, rows0 = list(A), list(b), list(rows)
    for choice in itertools.product((0, 1), repeat=len(diag)):
        A, b, rows = list(A0), list(b0), list(rows0)
        ch = dict(zip(diag, choice))
        for k, (i, j, ddx, ddy) in enumerate(disj):
            ax = 0 if abs(ddx) > 0.5 else 1
            if k in ch:
                ax = ch[k]
            sg = math.copysign(1.0, ddx if ax == 0 else ddy)
            #  gap = sg*(dd + delta_j - delta_i) - 1/2 - 1/2(|cos|+|sin|) >= g
            #      = sg*(delta_j - delta_i) - 1/2|phi_j - phi_i| >= g     (first order)
            for s2 in (1.0, -1.0):
                c = np.zeros(nv)
                c[0] = 1.0
                c[(1 if ax == 0 else 2) + 3 * j] = -sg
                c[(1 if ax == 0 else 2) + 3 * i] = sg
                c[3 + 3 * j] = 0.5 * s2
                c[3 + 3 * i] = -0.5 * s2
                A.append(c)
                b.append(0.0)
                rows.append(f'gap({LNAME[labels[i]]},{LNAME[labels[j]]}) axis {"xy"[ax]} '
                            f's{int(s2):+d}')
        # normalisation so the homogeneous cone is bounded: |delta| <= 1, |phi| <= 1
        bounds = [(None, None)] + [(-1.0, 1.0)] * (3 * n)
        c = np.zeros(nv)
        c[0] = -1.0
        res = linprog(c, A_ub=np.array(A), b_ub=np.array(b), bounds=bounds, method='highs')
        if res.status == 0 and -res.fun > best[0]:
            best = (-res.fun, res, (list(A), list(b), list(rows), choice))
    val, res, (A, b, rows, choice) = best
    print(f'\nLP value (the best first-order min-gap direction, normalised to |delta|,|phi| <= 1):'
          f'  {val:+.12e}')
    print(f'  branch {choice}; the system is homogeneous, so value 0 means NO first-order '
          f'direction opens every gap, and any value > 0 would scale to infinity.')
    if val <= 1e-9:
        print('  => the linearised disjointness system around the tiling is INFEASIBLE '
              '(first-order rigid: no perturbation makes all twelve pairwise gaps positive).')
    y = res.ineqlin.marginals if hasattr(res, 'ineqlin') else None
    if y is not None:
        w = -np.asarray(y)
        idx = np.argsort(-w)
        print('\n  dual (Farkas) certificate, the rows carrying it:')
        for k in idx[:a.ndual]:
            if w[k] <= 1e-9:
                break
            print(f'    {w[k]:10.6f}  {rows[k]}')
        print(f'    ({int((w > 1e-9).sum())} rows with positive multiplier of {len(rows)})')
    return val, rows, res


# ===================================================================== (4) Sherali-Adams level 2
def cmd_sa2(a):
    from scipy.optimize import linprog
    import scipy.sparse as sp
    t, sym, poses = lc.read_measure(a.FILE)
    sq = [lc.make_square(cx, cy, p, q, t) for (p, q, cx, cy, _m) in poses]
    want = {frozenset(s): k for k, (_n, s) in enumerate(LABELS)}
    lab = []
    for s in sq:
        lab.append(want.get(bentz.pattern_of(s, IPTS), -1))
    keep = [i for i in range(len(sq)) if lab[i] >= 0]
    sq = [sq[i] for i in keep]
    lab = [lab[i] for i in keep]
    n = len(sq)
    print(f'{a.FILE}: {n} support poses carrying one of the twelve leaf-A patterns')
    cls = [[i for i in range(n) if lab[i] == k] for k in range(12)]
    print('  per label: ' + ', '.join(f'{LNAME[k]}:{len(cls[k])}' for k in range(12)))
    dis = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            if not lc.sq_meets_sq(sq[i], sq[j]):
                dis[i, j] = dis[j, i] = True
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if dis[i, j]]
    pidx = {p: k for k, p in enumerate(pairs)}
    print(f'  {len(pairs)} of {n*(n-1)//2} pairs disjoint as closed sets')

    # variables: y_i (n), y_ij (len(pairs))
    nv = n + len(pairs)

    def Y(i):
        return i

    def Yp(i, j):
        return n + pidx[(min(i, j), max(i, j))]

    R, C, V, rhs = [], [], [], []
    r = 0

    def row(entries, b):
        nonlocal r
        for (c, v) in entries:
            R.append(r)
            C.append(c)
            V.append(v)
        rhs.append(b)
        r += 1

    # clique rows (one square per pattern): sum_{s in cls[k]} y_s <= 1
    for k in range(12):
        row([(Y(s), 1.0) for s in cls[k]], 1.0)
    # y_ij <= y_i, y_ij <= y_j, y_i + y_j - y_ij <= 1
    for (i, j) in pairs:
        row([(Yp(i, j), 1.0), (Y(i), -1.0)], 0.0)
        row([(Yp(i, j), 1.0), (Y(j), -1.0)], 0.0)
        row([(Y(i), 1.0), (Y(j), 1.0), (Yp(i, j), -1.0)], 1.0)
    # SA level 2: multiply each clique row by x_u and by (1 - x_u)
    for k in range(12):
        for u in range(n):
            if u in cls[k]:
                continue
            ent = [(Yp(s, u), 1.0) for s in cls[k] if dis[s, u]]
            row(ent + [(Y(u), -1.0)], 0.0)
            ent2 = [(Y(s), 1.0) for s in cls[k]] + [(Yp(s, u), -1.0) for s in cls[k] if dis[s, u]]
            row(ent2 + [(Y(u), 1.0)], 1.0)
    Aub = sp.coo_matrix((V, (R, C)), shape=(r, nv)).tocsr()
    c = np.zeros(nv)
    c[:n] = -1.0
    res = linprog(c, A_ub=Aub, b_ub=np.array(rhs), bounds=[(0, 1)] * nv, method='highs')
    print(f'\nSA level-2 value on the 162-pose support: {-res.fun:.9f}   '
          f'(LP with the clique rows alone: 12; integer optimum alpha = 11)')
    if res.status != 0:
        print('  LP status', res.status, res.message)
        return
    y = res.x
    w = -np.asarray(res.ineqlin.marginals)
    # which rows carry the dual
    names = []
    r = 0
    for k in range(12):
        names.append(('clique', LNAME[k], None))
    for (i, j) in pairs:
        names.append(('y_ij<=y_i', i, j))
        names.append(('y_ij<=y_j', i, j))
        names.append(('y_i+y_j-y_ij<=1', i, j))
    for k in range(12):
        for u in range(n):
            if u in cls[k]:
                continue
            names.append(('SA*x_u', LNAME[k], u))
            names.append(('SA*(1-x_u)', LNAME[k], u))
    idx = np.argsort(-w)
    print(f'\n  {int((w > 1e-9).sum())} of {len(w)} rows carry a positive dual; the largest:')
    shown = 0
    bylab = {}
    for k in idx:
        if w[k] <= 1e-9 or shown >= a.ndual:
            break
        kind, x, u = names[k]
        if kind in ('SA*x_u', 'SA*(1-x_u)'):
            print(f'    {w[k]:10.6f}  {kind} label {x} pose {u} (label {LNAME[lab[u]]})')
            bylab[(x, LNAME[lab[u]])] = bylab.get((x, LNAME[lab[u]]), 0.0) + w[k]
        elif kind == 'clique':
            print(f'    {w[k]:10.6f}  clique row {x}')
        else:
            print(f'    {w[k]:10.6f}  {kind} poses {x},{u} '
                  f'(labels {LNAME[lab[x]]},{LNAME[lab[u]]})')
        shown += 1
    tot = {}
    r = 0
    for k in range(12):
        tot[('clique', LNAME[k])] = tot.get(('clique', LNAME[k]), 0.0) + w[r]
        r += 1
    for (i, j) in pairs:
        key = tuple(sorted((LNAME[lab[i]], LNAME[lab[j]])))
        tot[key] = tot.get(key, 0.0) + w[r] + w[r + 1] + w[r + 2]
        r += 3
    for k in range(12):
        for u in range(n):
            if u in cls[k]:
                continue
            key = tuple(sorted((LNAME[k], LNAME[lab[u]])))
            tot[key] = tot.get(key, 0.0) + w[r] + w[r + 1]
            r += 2
    print('\n  dual mass aggregated by LABEL PAIR (which pairs carry the certificate):')
    for key, v in sorted(tot.items(), key=lambda kv: -kv[1])[:a.npairs]:
        if v <= 1e-9:
            break
        print(f'    {v:10.6f}  {key}')
    return -res.fun


# ===================================================================== exact verification
def cmd_analyse(a):
    """anatomy of a saved configuration: poses, contact graph, which constraints are tight"""
    for (tag, val, sqs) in load(a.FILE)[:a.n]:
        labels = [LIDX[nm] for (nm, _x, _y, _t) in sqs]
        z = z_from([(x, y, th) for (_nm, x, y, th) in sqs])
        n = len(labels)
        inc, exc = pattern_slack(z, labels)
        print(f'\n=== {tag}: min pairwise closed gap {val:+.6e}, {n} squares, '
              f'worst containment {inc:+.3e}, worst exclusion {exc:+.3e}, '
              f'max tilt {max_tilt_deg(z):.4f} deg, dist to tiling {_tiling_distance(z):.6f}')
        th = norm_theta(z[3::3])
        for k in range(n):
            print(f'   {LNAME[labels[k]]:<3} c = ({z[1+3*k]:9.6f}, {z[2+3*k]:9.6f})  '
                  f'theta = {math.degrees(th[k]):+9.5f} deg')
        g = all_gaps(z, n)
        touch = sorted((v, i, j) for (i, j), v in g.items() if v < a.tol)
        print(f'   contact graph: {len(touch)} of {n*(n-1)//2} pairs with gap < {a.tol}')
        print('     ' + ', '.join(f'{LNAME[labels[i]]}-{LNAME[labels[j]]}:{v:+.1e}'
                                  for (v, i, j) in touch))
        deg = {}
        for (_v, i, j) in touch:
            deg[i] = deg.get(i, 0) + 1
            deg[j] = deg.get(j, 0) + 1
        print('     contacts per square: '
              + ', '.join(f'{LNAME[labels[k]]}:{deg.get(k,0)}' for k in range(n)))


def cmd_verify(a):
    cfgs = load(a.FILE)
    Q, Dc = a.Q, a.Dc
    for (tag, val, sqs) in cfgs[:a.n]:
        recs, names, ok = [], [], True
        for (nm, x, y, th) in sqs:
            p, q, cx, cy = lc.snap_pose(x, y, th, Fr(4), Q, Dc)
            recs.append(lc.make_square(cx, cy, p, q, Fr(4)))
            names.append(nm)
        pats = [bentz.pattern_of(s, IPTS) for s in recs]
        want = [frozenset(LABELS[LIDX[nm]][1]) for nm in names]
        badpat = [(names[k], bentz.pname(pats[k], IPTS)) for k in range(len(recs))
                  if pats[k] != want[k]]
        meets = [(names[i], names[j]) for i in range(len(recs)) for j in range(i + 1, len(recs))
                 if lc.sq_meets_sq(recs[i], recs[j])]
        ok = not badpat and not meets
        print(f'{tag}: float value {val:+.6e} -> snapped (Q={Q}, Dc={Dc}): '
              f'{"EXACT PACKING IN THE LEAF" if ok else "not a packing"}')
        if badpat:
            print('   wrong pattern after snapping: '
                  + ', '.join(f'{n} -> {p}' for (n, p) in badpat[:6]))
        if meets:
            print(f'   {len(meets)} meeting pairs: ' + ', '.join(f'{i}-{j}' for (i, j) in meets[:8]))


# ===================================================================== selftest
def selftest():
    bad = []

    def ck(name, cond, extra=''):
        print(f'  {"ok  " if cond else "FAIL"}  {name}' + (f'   {extra}' if extra else ''))
        if not cond:
            bad.append(name)

    rng = np.random.default_rng(7)
    # the float SAT gap agrees with leaf_ceiling's exact sq_meets_sq
    agree = 0
    for _ in range(4000):
        th = rng.uniform(0, math.pi / 2, 2)
        w = 0.5 * (np.abs(np.cos(th)) + np.abs(np.sin(th)))
        x = rng.uniform(w[0], T - w[0])
        y = rng.uniform(w[0], T - w[0])
        x2 = rng.uniform(w[1], T - w[1])
        y2 = rng.uniform(w[1], T - w[1])
        z = np.array([0.0, x, y, th[0], x2, y2, th[1]])
        g = pair_gap(z, 0, 1)
        if abs(g) < 1e-9:
            continue
        s = []
        for k in range(2):
            p, q, cx, cy = lc.snap_pose(z[1 + 3 * k], z[2 + 3 * k], z[3 + 3 * k],
                                        Fr(4), 10 ** 7, 10 ** 9)
            s.append(lc.make_square(cx, cy, p, q, Fr(4)))
        if (g > 0) == (not lc.sq_meets_sq(s[0], s[1])):
            agree += 1
        else:
            agree -= 10000
    ck('the float SAT gap agrees in sign with leaf_ceiling.sq_meets_sq', agree > 3500,
       f'{agree} agreements of ~4000 samples')
    # float pattern test vs the exact one
    ok = True
    for _ in range(2000):
        th = rng.uniform(0, math.pi / 2)
        w = 0.5 * (abs(math.cos(th)) + abs(math.sin(th)))
        x, y = rng.uniform(w, T - w), rng.uniform(w, T - w)
        z = np.array([0.0, x, y, th])
        pf = patterns_float(z)[0]
        p, q, cx, cy = lc.snap_pose(x, y, th, Fr(4), 10 ** 7, 10 ** 9)
        pe = bentz.pattern_of(lc.make_square(cx, cy, p, q, Fr(4)), IPTS)
        m = np.max(np.abs(margins(z, None)), axis=2)[0]
        if pf != pe and np.min(np.abs(m - 0.5)) > 1e-5:
            ok = False
    ck('the float pattern test agrees with the exact one away from the boundary', ok)
    # the tiling-minus-4 families are in the relaxed leaf with min gap exactly 0
    st = tiling_starts()
    ck('there are 16 tiling-minus-4 families', len(st) == 16, f'{len(st)}')
    vals = []
    for (_pick, z) in st:
        inc, exc = pattern_slack(z, list(range(12)))
        v = true_objective(z, 12)[0]
        vals.append(v)
        if inc < -1e-12 or exc >= 0:
            bad.append('tiling family not in the leaf')
    ck('every tiling-minus-4 family lies in leaf A (strict exclusion) ...',
       all(pattern_slack(z, list(range(12)))[1] < 0 for (_p, z) in st))
    ck('... with min pairwise closed gap exactly 0 (the tiles touch)',
       all(abs(v) < 1e-12 for v in vals), f'max |gap| = {max(abs(v) for v in vals):.3e}')
    # analytic jacobian against finite differences
    pr = Problem(list(range(12)))
    z = st[0][1].copy()
    z[0] = 0.0
    z[1:] += rng.normal(0, 0.01, 36)
    pr.build(pr.pair_assign(z), pr.exc_assign(z))
    J = pr.jac(z)
    Jn = np.zeros_like(J)
    for k in range(len(z)):
        e = np.zeros_like(z)
        e[k] = 1e-7
        Jn[:, k] = (pr.cons(z + e) - pr.cons(z - e)) / 2e-7
    err = float(np.max(np.abs(J - Jn)))
    ck('the analytic constraint jacobian matches finite differences', err < 1e-5,
       f'max error {err:.3e}')
    print('rank8 selftest:', 'PASS' if not bad else f'FAIL ({len(bad)})')
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('margin')
    c.add_argument('--seed', type=int, default=20260913)
    c.add_argument('--nudge', type=int, default=120)
    c.add_argument('--nrandom', type=int, default=200)
    c.add_argument('--nsupport', type=int, default=80)
    c.add_argument('--nsamp', type=int, default=2000000)
    c.add_argument('--sigma', type=float, default=0.04)
    c.add_argument('--sigma-th', type=float, default=0.05)
    c.add_argument('--eps', type=float, default=0.0)
    c.add_argument('--rounds', type=int, default=8)
    c.add_argument('--maxiter', type=int, default=300)
    c.add_argument('--near', type=float, default=1e-3)
    c.add_argument('--tilt-deg', type=float, default=1.0)
    c.add_argument('--far', type=float, default=0.05)
    c.add_argument('--support', default=None)
    c.add_argument('--out', default=None)
    c.add_argument('--nsave', type=int, default=25)
    c.add_argument('--nproc', type=int, default=8)
    c = sub.add_parser('subsets')
    c.add_argument('--seed', type=int, default=20260913)
    c.add_argument('--nstart', type=int, default=30)
    c.add_argument('--nsamp', type=int, default=2000000)
    c.add_argument('--eps', type=float, default=0.0)
    c.add_argument('--rounds', type=int, default=8)
    c.add_argument('--maxiter', type=int, default=300)
    c.add_argument('--pos', type=float, default=1e-7)
    c.add_argument('--support', default=None)
    c.add_argument('--out', default=None)
    c.add_argument('--nproc', type=int, default=8)
    c = sub.add_parser('linear')
    c.add_argument('--pick', default='0000')
    c.add_argument('--ndual', type=int, default=40)
    c = sub.add_parser('sa2')
    c.add_argument('FILE')
    c.add_argument('--ndual', type=int, default=25)
    c.add_argument('--npairs', type=int, default=25)
    c = sub.add_parser('analyse')
    c.add_argument('FILE')
    c.add_argument('-n', type=int, default=6)
    c.add_argument('--tol', type=float, default=1e-6)
    c = sub.add_parser('verify')
    c.add_argument('FILE')
    c.add_argument('-n', type=int, default=10)
    c.add_argument('--Q', type=int, default=10 ** 7)
    c.add_argument('--Dc', type=int, default=10 ** 9)
    sub.add_parser('selftest')
    a = ap.parse_args()
    if a.cmd == 'margin':
        cmd_margin(a)
    elif a.cmd == 'subsets':
        cmd_subsets(a)
    elif a.cmd == 'linear':
        cmd_linear(a)
    elif a.cmd == 'sa2':
        cmd_sa2(a)
    elif a.cmd == 'analyse':
        cmd_analyse(a)
    elif a.cmd == 'verify':
        cmd_verify(a)
    elif a.cmd == 'selftest':
        sys.exit(selftest())


if __name__ == '__main__':
    main()
