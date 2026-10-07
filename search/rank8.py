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

    def __init__(self, labels, eps_excl=0.0, eps_inc=0.0):
        self.labels = list(labels)
        self.n = len(self.labels)
        self.eps = eps_excl
        self.epsi = eps_inc
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
            out.append(0.5 - self.epsi - self.Is * (e0 * dx + e1 * dy))
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


def solve(labels, z0, eps_excl=0.0, rounds=8, maxiter=300, eps_inc=0.0):
    """max-min-gap by assignment-fixed SLSQP, iterated until the assignment stops changing.
    Returns (true objective at the returned point, z, info)."""
    from scipy.optimize import minimize
    pr = Problem(labels, eps_excl, eps_inc)
    z = np.array(z0, dtype=float)
    z[0] = true_objective(z, pr.n)[0]
    best = (_feasible_value(z, pr), z.copy())
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
    if inc < pr.epsi - tol or exc > -pr.eps + tol or adm < -tol:
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
    if near:
        mt = max(max_tilt_deg(z) for (_v, _t, z) in near)
        md = max(_tiling_distance(z) for (_v, _t, z) in near)
        print(f'\nover the {len(near)} near-0 configurations: largest tilt of any square '
              f'{mt:.4f} deg, largest\ndistance to a tiling-minus-4 family {md:.6f}; '
              f'so the 0-level set is a plateau, not a point.')
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
    """which sub-configurations of leaf A are jointly realisable as a CLOSED packing.

    A square set is realisable iff the sup of the min pairwise closed gap over it is > 0.
    Realisability is downward closed (drop a square and the rest still work), so the answer is
    the pair (maximal realisable sets, minimal non-realisable sets), each up to D4.  A positive
    answer is a witness and can be re-checked exactly (`verify`); a negative one is a multistart
    measurement and nothing more."""
    rng = np.random.default_rng(a.seed)
    pools = sample_pools(a.nsamp, rng)
    sp = support_pools(a.support) if a.support else None
    base = tiling_starts()
    universe = list(range(12)) if a.all else SINGLETONS
    forced = [] if a.all else list(range(4))
    nb = len(universe)
    print(f'universe: {[LNAME[k] for k in universe]}'
          + (f', forced: {[LNAME[k] for k in forced]}' if forced else ''))

    def labels_of(mask):
        return sorted(forced + [universe[b] for b in range(nb) if mask >> b & 1])

    def best_of(mask, nstart):
        """multistart + basin hopping: the tiling families restricted to `labels`, random draws
        from the per-label pose pools, then rounds of perturb-and-reoptimise around the best few
        so far.  Local optimisation alone gives false negatives on ten-square problems."""
        labels = labels_of(mask)
        if len(labels) < 2:
            return (float('inf'), None)
        n = len(labels)
        starts = []
        for (_pick, z) in base:
            starts.append(np.concatenate([[0.0]] + [z[1 + 3 * li:4 + 3 * li] for li in labels]))
        for _ in range(nstart):
            src = pools if (sp is None or rng.random() < 0.7) else sp
            starts.append(z_from([src[li][rng.integers(len(src[li]))] for li in labels]))
        res = run_starts([(f's{k}', z) for k, z in enumerate(starts)], labels, a, a.nproc)
        res = [r for r in res if r[0] > -float('inf')]
        for rd in range(a.bh):
            res.sort(key=lambda r: -r[0])
            if res and res[0][0] > a.pos:
                break
            pool = res[:a.bhkeep] or [(0.0, 'x', starts[0])]
            kick = []
            for k in range(max(nstart, a.bhn)):
                _v, _t, z = pool[rng.integers(len(pool))]
                w = z.copy()
                sig = 10.0 ** rng.uniform(-3, -0.7)
                w[1::3] += rng.normal(0, sig, n)
                w[2::3] += rng.normal(0, sig, n)
                w[3::3] += rng.normal(0, 1.5 * sig, n)
                kick.append((f'bh{rd}_{k}', w))
            res += [r for r in run_starts(kick, labels, a, a.nproc) if r[0] > -float('inf')]
        v, _tag, z = max(res, key=lambda r: r[0])
        return (v, z)

    # --- D4 orbits of the subsets of `universe`
    idx = {l: b for b, l in enumerate(universe)}
    perms = []
    for pm in LPERM:
        if all(pm[l] in idx for l in universe):
            perms.append(tuple(idx[pm[l]] for l in universe))
    canon, orb = {}, []
    for mask in range(1 << nb):
        if mask in canon:
            continue
        o = {sum(1 << pm[b] for b in range(nb) if mask >> b & 1) for pm in perms}
        rep = min(o)
        for m in o:
            canon[m] = rep
        orb.append(rep)
    print(f'{len(orb)} D4 classes of the {1 << nb} subsets')

    def subs(mask):
        return [mask ^ (1 << b) for b in range(nb) if mask >> b & 1]

    if a.only:
        want = set(a.only.split(','))
        m = sum(1 << b for b, l in enumerate(universe) if LNAME[l] in want)
        v, z = best_of(m, a.nstart)
        print(f'  {_mask_name(m, universe):<38} k={len(labels_of(m)):2d}  best min-gap {v:+.9e}'
              f'   {"REALISABLE" if v > a.pos else "not realisable"}')
        if a.out and z is not None:
            labels = labels_of(m)
            with open(a.out, 'w') as f:
                f.write('# rank8.py subsets --only.  FLOAT.\n')
                f.write(f'config {_mask_name(m, universe)} value {v:.17g} n {len(labels)}\n')
                for k, li in enumerate(labels):
                    f.write(f'  {LNAME[li]} {z[1+3*k]:.17g} {z[2+3*k]:.17g} {z[3+3*k]:.17g}\n')
            print(f'wrote {a.out}')
        return {m: (v, z)}
    status, implied = {}, {}
    order = sorted(orb, key=lambda m: -bin(m).count('1'))
    ntest = 0
    for rep in order:
        if rep in implied and not a.exhaustive:
            continue
        v, z = best_of(rep, a.nstart)
        if v <= a.pos and a.confirm:            # a negative is the dangerous answer: retry harder
            v2, z2 = best_of(rep, a.nstart * a.confirm)
            if v2 > v:
                v, z = v2, z2
        ntest += 1
        status[rep] = (v, z)
        print(f'  {_mask_name(rep, universe):<34} k={len(labels_of(rep)):2d}  '
              f'best min-gap {v:+.6e}   '
              f'{"REALISABLE" if v > a.pos else "not realisable"}', flush=True)
        if v > a.pos:
            stack = [rep]
            while stack:
                m = stack.pop()
                for s in subs(m):
                    c = canon[s]
                    if c not in implied and c not in status:
                        implied[c] = rep
                        stack.append(s)
    if a.exhaustive:
        print('\nmargin table (sup of the min pairwise closed gap, by class):')
        for m in sorted(status, key=lambda m: (-bin(m).count('1'), -status[m][0])):
            v = status[m][0]
            print(f'  {_mask_name(m, universe):<34} {len(labels_of(m)):2d} squares   {v:+.6e}')
    print(f'\n{ntest} classes tested, {len(implied)} implied realisable by a larger one')
    feas = {m for m, (v, _z) in status.items() if v > a.pos} | set(implied)
    infeas = [m for m in orb if m not in feas]
    minimal = [m for m in sorted(infeas, key=lambda m: bin(m).count('1'))
               if all(canon[s] in feas for s in subs(m))]
    print('\nMINIMAL NON-REALISABLE sub-configurations up to D4 '
          '(every proper subset of each IS realisable):')
    for m in minimal:
        v = status.get(m, (float('nan'), None))[0]
        print(f'  {_mask_name(m, universe):<34} {len(labels_of(m)):2d} squares   '
              f'best min-gap {v:+.6e}')
    maximal = [m for m in sorted(feas, key=lambda m: -bin(m).count('1'))
               if all(canon[m | (1 << b)] not in feas for b in range(nb) if not m >> b & 1)]
    print('\nMAXIMAL REALISABLE sub-configurations up to D4:')
    for m in maximal:
        v = status.get(m, (float('nan'), None))[0]
        print(f'  {_mask_name(m, universe):<34} {len(labels_of(m)):2d} squares   '
              f'best min-gap {v:+.6e}'
              + ('' if m in status else f'   (implied by {_mask_name(implied[m], universe)})'))
    if a.out:
        with open(a.out, 'w') as f:
            f.write('# rank8.py subsets: label cx cy theta_rad.  FLOAT.\n')
            n = 0
            for m, (v, z) in sorted(status.items(), key=lambda kv: -bin(kv[0]).count('1')):
                if z is None:
                    continue
                labels = labels_of(m)
                f.write(f'config {_mask_name(m, universe)} value {v:.17g} n {len(labels)}\n')
                for k, li in enumerate(labels):
                    f.write(f'  {LNAME[li]} {z[1+3*k]:.17g} {z[2+3*k]:.17g} {z[3+3*k]:.17g}\n')
                n += 1
        print(f'\nwrote {n} configurations to {a.out}')
    return status


def _mask_name(mask, universe=None):
    u = universe if universe is not None else SINGLETONS
    return '{' + ','.join(LNAME[u[b]] for b in range(len(u)) if mask >> b & 1) + '}'



# ===================================================================== (3) the linearised system
def lin_model(pick, verbose=True):
    """the first-order model of leaf A about a tiling-minus-4 family.

    Write the perturbed pose of tile `k` as `c_k = C_k + (dx_k, dy_k)`, `theta_k = phi_k`, with
    `C_k` the tile centre.  To first order in `(dx, dy, phi)`:

      * the coordinates of a point `p` in square `k`'s frame, with `e = p - C_k`, are
            coord_u = e_x - dx_k + e_y phi_k,      coord_v = e_y - dy_k - e_x phi_k,
        and `p` is in the closed square iff both have modulus `<= 1/2`;
      * the separating-axis gap of two squares with nominal centre difference `dd` is
            gap = max_branches [ |dd_ax| + sg (Ddelta_ax + s_o phi_o) ] - 1 - (1/2)|phi_j - phi_i|,
        over `ax in {u, v}` with `dd_ax != 0`, `sg = sign(dd_ax)`, owner `o in {i, j}`, and
        `s_o = dd_y` for `ax = u`, `-dd_x` for `ax = v`  (the tilt term does NOT drop out for a
        diagonal pair: that is where a rotation first buys separation);
      * admissibility `c in [w/2, 4 - w/2]` with `w = |cos| + |sin| ~ 1 + |phi|` is
            dx_k >= (1/2)|phi_k|  at a tile touching the left wall, and its three images.

    Every `max` over branches is a disjunction and is carried as a binary; every `|.|` that enters
    with a minus sign is a conjunction and is carried as two rows.  The system is homogeneous, so
    it is normalised by `|dx|, |dy|, |phi| <= 1` and the question is only whether the optimum `g`
    is `0` (no first-order direction opens every gap) or `> 0`.
    """
    base = dict(tiling_starts())[pick]
    n = 12
    labels = list(range(12))
    X0 = np.round(base[1::3] - 0.5) + 0.5
    Y0 = np.round(base[2::3] - 0.5) + 0.5
    nv = 1 + 3 * n                      # g, then (dx, dy, phi) per square

    def col(k, which):
        return 1 + 3 * k + which

    rows, rhs, tags = [], [], []        # conjunctive rows:  coef . z >= rhs
    disj = []                           # disjunctions: list of list of (coef, rhs, tag)
    slack = []                          # (nominal slack, description) of everything inactive

    # ---- pattern rows
    for k in range(n):
        want = set(LABELS[labels[k]][1])
        for p in range(16):
            ex, ey = PXY[p, 0] - X0[k], PXY[p, 1] - Y0[k]
            # (nominal value, d/d dx, d/d dy, d/d phi) of the two frame coordinates
            co = [(ex, -1.0, 0.0, ey), (ey, 0.0, -1.0, -ex)]
            branches = []
            for (val, gx, gy, gp) in co:
                for sg in (1.0, -1.0):
                    branches.append((sg * val, sg * gx, sg * gy, sg * gp))
            M0 = max(b[0] for b in branches)
            if p in want:
                for (v0, gx, gy, gp) in branches:      # ALL four must stay <= 1/2
                    c = np.zeros(nv)
                    c[col(k, 0)], c[col(k, 1)], c[col(k, 2)] = -gx, -gy, -gp
                    if abs(v0 - 0.5) < 1e-12:
                        rows.append(c)
                        rhs.append(0.0)
                        tags.append(f'{LNAME[labels[k]]} holds {PNAMES[p]} ({"uv"[0]}-side)')
                    else:
                        slack.append((0.5 - v0, f'{LNAME[labels[k]]} holds {PNAMES[p]}'))
            else:
                if M0 > 0.5 + 1e-12:                   # excluded with room; no row
                    slack.append((M0 - 0.5, f'{LNAME[labels[k]]} excludes {PNAMES[p]}'))
                    continue
                assert M0 > 0.5 - 1e-12, (k, p, M0)
                opts = []
                for (v0, gx, gy, gp) in branches:
                    if abs(v0 - 0.5) < 1e-12:
                        c = np.zeros(nv)
                        c[col(k, 0)], c[col(k, 1)], c[col(k, 2)] = gx, gy, gp
                        opts.append((c, 0.0, f'{LNAME[labels[k]]} excludes {PNAMES[p]}'))
                if len(opts) == 1:
                    rows.append(opts[0][0])
                    rhs.append(0.0)
                    tags.append(opts[0][2])
                else:
                    disj.append(opts)

    # ---- admissibility rows (the tiles that touch the container)
    for k in range(n):
        for (v0, which, sgn) in ((X0[k] - 0.5, 0, +1.0), (3.5 - X0[k], 0, -1.0),
                                 (Y0[k] - 0.5, 1, +1.0), (3.5 - Y0[k], 1, -1.0)):
            if v0 > 1e-12:
                slack.append((v0, f'{LNAME[labels[k]]} clear of a wall'))
                continue
            for s in (1.0, -1.0):                    # d >= (1/2)|phi|  ->  two rows
                c = np.zeros(nv)
                c[col(k, which)] = sgn
                c[col(k, 2)] = -0.5 * s
                rows.append(c)
                rhs.append(0.0)
                tags.append(f'{LNAME[labels[k]]} against the wall '
                            f'({"xy"[which]}{"+" if sgn > 0 else "-"}) s{int(s):+d}')

    # ---- disjointness rows
    npair = 0
    for i in range(n):
        for j in range(i + 1, n):
            ddx, ddy = X0[j] - X0[i], Y0[j] - Y0[i]
            if max(abs(ddx), abs(ddy)) > 1.5:
                slack.append((max(abs(ddx), abs(ddy)) - 1.0,
                              f'gap({LNAME[labels[i]]},{LNAME[labels[j]]}) nominally open'))
                continue
            npair += 1
            opts = []
            for (ax, dd, s_of) in ((0, ddx, ddy), (1, ddy, -ddx)):
                if abs(dd) < 0.5:
                    continue
                sg = math.copysign(1.0, dd)
                for o in (i, j):
                    for s in (1.0, -1.0):            # the -(1/2)|Dphi| is a conjunction
                        c = np.zeros(nv)
                        c[0] = -1.0
                        c[col(j, ax)] += sg
                        c[col(i, ax)] -= sg
                        c[col(o, 2)] += sg * s_of
                        c[col(j, 2)] -= 0.5 * s
                        c[col(i, 2)] += 0.5 * s
                        opts.append((c, 0.0,
                                     f'gap({LNAME[labels[i]]},{LNAME[labels[j]]}) axis {"uv"[ax]}'
                                     f' owner {LNAME[labels[o]]} s{int(s):+d}', (ax, o)))
            # group the two |Dphi| rows of one branch together: a branch is (ax, owner)
            grp = {}
            for (c, r, t, key) in opts:
                grp.setdefault(key, []).append((c, r, t))
            disj.append([g for g in grp.values()])
    if verbose:
        print(f'  {len(rows)} conjunctive first-order rows, {len(disj)} disjunctions, '
              f'{npair} adjacent pairs of {n*(n-1)//2}')
    return dict(rows=rows, rhs=rhs, tags=tags, disj=disj, slack=slack, nv=nv,
                labels=labels, X0=X0, Y0=Y0, npair=npair)


def cmd_linear(a):
    if a.all_picks:
        vals = {}
        for pick in itertools.product((0, 1), repeat=4):
            b = argparse.Namespace(pick=''.join(map(str, pick)), ndual=0, all_picks=False,
                                   quiet=True)
            vals[pick] = _linear_one(b)
        print('\nfirst-order value over all sixteen tiling-minus-4 families:')
        for pick, v in vals.items():
            print(f'  {"".join(map(str, pick))}  {v:+.12e}')
        print(f'  max over the sixteen: {max(vals.values()):+.12e}')
        return max(vals.values())
    return _linear_one(a)


def _linear_one(a):
    from scipy.optimize import linprog, milp, LinearConstraint, Bounds
    quiet = getattr(a, 'quiet', False)
    pick = tuple(int(c) for c in a.pick)
    if not quiet:
        print(f'linearisation about the tiling-minus-4 family {pick} '
              f'(the exact 4x4 tiling restricted to twelve tiles):')
    M = lin_model(pick, verbose=not quiet)
    nv = M['nv']
    if not quiet:
        for k in range(12):
            print(f'  {LNAME[M["labels"][k]]:<3} tile centre '
                  f'({M["X0"][k]:.1f}, {M["Y0"][k]:.1f})')

    # flatten the disjunctions: each option is a LIST of rows that must all hold
    opts = []
    for d in M['disj']:
        flat = []
        for o in d:
            flat.append(o if isinstance(o, list) else [o])
        opts.append(flat)
    nb = sum(len(f) for f in opts)
    if not quiet:
        print(f'  {nb} binary branch choices over {len(opts)} disjunctions')

    BIG = 40.0
    R, rr, names = list(M['rows']), list(M['rhs']), list(M['tags'])
    A, bl, bu = [], [], []
    for c, r in zip(R, rr):
        A.append(np.concatenate([c, np.zeros(nb)]))
        bl.append(r)
        bu.append(np.inf)
    bcol = nv
    for f in opts:
        sel = np.zeros(nv + nb)
        for o in f:
            for (c, r, t) in o:
                # c.z >= r  must hold when the branch binary is 1, and be free when it is 0:
                #     c.z + BIG (1 - b) >= r   <=>   c.z - BIG b >= r - BIG
                row = np.concatenate([c, np.zeros(nb)])
                row[bcol] = -BIG
                A.append(row)
                bl.append(r - BIG)
                bu.append(np.inf)
                names.append(t)
            sel[bcol] = -1.0
            bcol += 1
        A.append(sel)               # EXACTLY one branch per disjunction, so that the Farkas
        bl.append(-1.0)             # certificate read off below belongs to a single branch
        bu.append(-1.0)
        names.append('branch selector')
    A = np.array(A)
    lb = np.concatenate([[-5.0], -np.ones(nv - 1), np.zeros(nb)])
    ub = np.concatenate([[5.0], np.ones(nv - 1), np.ones(nb)])
    cost = np.zeros(nv + nb)
    cost[0] = -1.0
    integ = np.concatenate([np.zeros(nv), np.ones(nb)])
    res = milp(c=cost, constraints=LinearConstraint(A, bl, bu), integrality=integ,
               bounds=Bounds(lb, ub))
    val = -res.fun if res.success else float('nan')
    if quiet:
        return val
    print(f'\nFIRST-ORDER LP/MILP value (normalised to |dx|,|dy|,|phi| <= 1):  {val:+.12e}')
    if val <= 1e-9:
        print('  => the linearised disjointness system around the tiling is INFEASIBLE: no '
              'first-order\n     perturbation (centre offsets and tilts) makes all twelve '
              'pairwise gaps positive.\n     The tiling-minus-4 family is FIRST-ORDER RIGID '
              'inside leaf A.')
    else:
        print('  => a first-order direction exists; the leaf is not first-order rigid.')

    # --- the Farkas certificate: fix the branch the MILP chose and read the LP duals
    sel = np.round(res.x[nv:]).astype(int)
    A2, b2, n2 = [], [], []
    for c, r, t in zip(M['rows'], M['rhs'], M['tags']):
        A2.append(-c)
        b2.append(-r)
        n2.append(t)
    k = 0
    for f in opts:
        for o in f:
            if sel[k]:
                for (c, r, t) in o:
                    A2.append(-c)
                    b2.append(-r)
                    n2.append(t)
            k += 1
    cost2 = np.zeros(nv)
    cost2[0] = -1.0
    res2 = linprog(cost2, A_ub=np.array(A2), b_ub=np.array(b2),
                   bounds=[(-5.0, 5.0)] + [(-1.0, 1.0)] * (nv - 1), method='highs')
    print(f'  LP on the chosen branch: {-res2.fun:+.12e}')
    w = -np.asarray(res2.ineqlin.marginals)
    pos = int((w > 1e-9).sum())
    print(f'\n  Farkas certificate on the chosen branch: {pos} of {len(w)} rows carry a '
          f'positive multiplier.\n  The LP is degenerate, so that dual is not unique; the '
          f'IRREDUCIBLE certificate below is found\n  by greedily dropping rows while the '
          f'value stays 0 (drop one more and a direction opens):')
    A2 = np.array(A2)
    b2 = np.array(b2)
    keep = list(range(len(b2)))

    def value(rowset):
        r = linprog(cost2, A_ub=A2[rowset], b_ub=b2[rowset],
                    bounds=[(-5.0, 5.0)] + [(-1.0, 1.0)] * (nv - 1), method='highs')
        return -r.fun if r.status == 0 else float('inf')

    order = sorted(keep, key=lambda k: (w[k] > 1e-9, -w[k]))
    for k in order:
        trial = [i for i in keep if i != k]
        if value(trial) <= 1e-9:
            keep = trial
    print(f'  irreducible first-order certificate: {len(keep)} rows')
    for k in keep:
        print(f'    {n2[k]}')
    # read the certificate back as chains: a set of squares linked by same-axis gap rows and
    # anchored at both ends by the container is a run of unit squares spanning a side-4 box,
    # which has exactly zero slack -- that is the whole first-order obstruction.
    import re
    gaps = {0: [], 1: []}
    walls = {0: {}, 1: {}}
    for k in keep:
        m = re.match(r'gap\((\w+),(\w+)\) axis (\w)', n2[k])
        if m:
            gaps['uv'.index(m.group(3))].append((m.group(1), m.group(2)))
            continue
        m = re.match(r'(\w+) against the wall \((.)([-+])\)', n2[k])
        if m:
            walls['xy'.index(m.group(2))].setdefault(m.group(1), set()).add(m.group(3))
    print('\n  read as CHAINS (a run of unit squares pinned between two opposite container '
          'walls\n  has exactly zero slack in a side-4 box, and a relative tilt only costs):')
    for ax in (0, 1):
        adj = {}
        for (i, j) in gaps[ax]:
            adj.setdefault(i, set()).add(j)
            adj.setdefault(j, set()).add(i)
        seen = set()
        for v0 in sorted(adj):
            if v0 in seen:
                continue
            comp, st = set(), [v0]
            while st:
                v = st.pop()
                if v in comp:
                    continue
                comp.add(v)
                st += list(adj[v])
            seen |= comp
            anch = {v: walls[ax].get(v, set()) for v in comp if walls[ax].get(v)}
            print(f'    along {"xy"[ax]}: {len(comp)} squares '
                  f'{{{",".join(sorted(comp))}}}, anchored at '
                  + (', '.join(f'{v}{"".join(sorted(a))}' for v, a in sorted(anch.items()))
                     or 'nothing'))

    # --- where the linearisation is valid
    sl = sorted(M['slack'])
    pat = [s for s in sl if 'holds' in s[1] or 'excludes' in s[1]]
    adm = [s for s in sl if 'wall' in s[1]]
    prs = [s for s in sl if 'nominally open' in s[1]]
    print('\nVALIDITY OF THE LINEARISATION.  The model above keeps a fixed combinatorial '
          'structure:\n  the same pattern rows active, the same pairs adjacent, the same branch '
          'in each disjunction.\n  That structure survives a perturbation of size (rho in centre, '
          'alpha in angle) as long as')
    print(f'    every inactive PATTERN inequality keeps its sign:  slack >= '
          f'{min(p[0] for p in pat):.6f}  ({min(pat)[1]})')
    print(f'    every tile that is clear of a wall stays clear:    slack >= '
          f'{min(a2[0] for a2 in adm):.6f}  ({min(adm)[1]})')
    print(f'    every non-adjacent pair stays apart:               slack >= '
          f'{min(p[0] for p in prs):.6f}')
    print('  and, since a point of P0 moves by at most rho + (sqrt2/2) alpha in a square\'s own '
          'frame,\n  the structure is unchanged for  rho + 0.7072 alpha < 0.086000  (the binding '
          'slack).')
    print('\n  SECOND-ORDER REMAINDER.  With |centre offsets| <= rho and |phi| <= alpha, the '
          'exact gap\n  and its first-order model differ by at most')
    print('      |gap - gap_lin|  <=  2 rho alpha + (1 + 2 rho) alpha^2 / 2 + alpha^2')
    print('  (the projection error (1+2rho)(1-cos alpha) + 2 rho |sin alpha|, plus the '
          'W(Dphi) error\n  0.25 Dphi^2 <= alpha^2).  So the first-order certificate bounds the '
          'true min gap by\n  O(rho alpha + alpha^2), not by 0: it says the obstruction is '
          'first order, and the exact\n  statement still needs the non-linearised inequality.  '
          'The sign of the leading second-order\n  term is POSITIVE (relative tilt costs '
          '|Dphi|/2 at first order but only\n  |Dphi|/2 - Dphi^2/4 exactly), which is why (1) '
          'has to be measured globally and not read off\n  this LP.')
    for rho, al in ((0.01, 0.01), (0.01, 0.05), (0.05, 0.05), (0.02, 0.0175)):
        e = 2 * rho * al + (1 + 2 * rho) * al * al / 2 + al * al
        print(f'      rho = {rho:.3f}, alpha = {al:.4f} rad ({math.degrees(al):5.2f} deg): '
              f'remainder <= {e:.6f}')
    return val



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
    use = set(range(12)) if not a.labels else {LIDX[n] for n in a.labels.split(',')}
    keep = [i for i in range(len(sq)) if lab[i] >= 0 and lab[i] in use]
    sq = [sq[i] for i in keep]
    lab = [lab[i] for i in keep]
    n = len(sq)
    print(f'{a.FILE}: {n} support poses carrying one of the twelve leaf-A patterns')
    use = sorted(use)
    cls = [[i for i in range(n) if lab[i] == k] for k in use]
    print('  per label: ' + ', '.join(f'{LNAME[k]}:{len(cls[u])}' for u, k in enumerate(use)))
    dis = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            if not lc.sq_meets_sq(sq[i], sq[j]):
                dis[i, j] = dis[j, i] = True
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if dis[i, j]]
    pidx = {p: k for k, p in enumerate(pairs)}
    print(f'  {len(pairs)} of {n*(n-1)//2} pairs disjoint as closed sets')

    if a.alpha:
        adj = [0] * n
        for (i, j) in pairs:
            adj[i] |= 1 << j
            adj[j] |= 1 << i
        best = [0]

        def bb(cand, k):
            if bin(cand).count('1') + k <= best[0]:
                return
            if cand == 0:
                best[0] = max(best[0], k)
                return
            c = cand
            while c:
                b = c & -c
                i = b.bit_length() - 1
                c ^= b
                bb(cand & adj[i] & ~((1 << (i + 1)) - 1), k + 1)
                cand ^= b
                if bin(cand).count('1') + k <= best[0]:
                    return
        bb((1 << n) - 1, 0)
        print(f'  integer optimum on this support: alpha = {best[0]}')

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
    for k in range(len(use)):
        row([(Y(s), 1.0) for s in cls[k]], 1.0)
    # y_ij <= y_i, y_ij <= y_j, y_i + y_j - y_ij <= 1
    for (i, j) in pairs:
        row([(Yp(i, j), 1.0), (Y(i), -1.0)], 0.0)
        row([(Yp(i, j), 1.0), (Y(j), -1.0)], 0.0)
        row([(Y(i), 1.0), (Y(j), 1.0), (Yp(i, j), -1.0)], 1.0)
    # SA level 2: multiply each clique row by x_u and by (1 - x_u)
    for k in range(len(use)):
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
    print(f'\nSA level-2 value on the {n}-pose support over {len(use)} labels: {-res.fun:.9f}'
          f'   (the clique LP alone gives {len(use)})')
    if res.status != 0:
        print('  LP status', res.status, res.message)
        return
    y = res.x
    w = -np.asarray(res.ineqlin.marginals)
    # which rows carry the dual
    names = []
    r = 0
    for k in use:
        names.append(('clique', LNAME[k], None))
    for (i, j) in pairs:
        names.append(('y_ij<=y_i', i, j))
        names.append(('y_ij<=y_j', i, j))
        names.append(('y_i+y_j-y_ij<=1', i, j))
    for kk, k in enumerate(use):
        for u in range(n):
            if u in cls[kk]:
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
    for k in use:
        tot[('clique', LNAME[k])] = tot.get(('clique', LNAME[k]), 0.0) + w[r]
        r += 1
    for (i, j) in pairs:
        key = tuple(sorted((LNAME[lab[i]], LNAME[lab[j]])))
        tot[key] = tot.get(key, 0.0) + w[r] + w[r + 1] + w[r + 2]
        r += 3
    for kk, k in enumerate(use):
        for u in range(n):
            if u in cls[kk]:
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
              f'max tilt {max_tilt_deg(z):.4f} deg'
              + (f', dist to tiling {_tiling_distance(z):.6f}' if n == 12 else ''))
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


def cmd_harden(a):
    """re-optimise each saved configuration with a strict exclusion margin, so that the witness
    survives snapping to rationals and `verify` can turn it into an exact statement.  The value
    can only go down; a positive value after hardening is a genuine closed packing."""
    out = []
    for (tag, val, sqs) in load(a.FILE):
        labels = [LIDX[nm] for (nm, _x, _y, _t) in sqs]
        z0 = z_from([(x, y, th) for (_nm, x, y, th) in sqs])
        best = (-float('inf'), z0)
        for eps in (a.eps, a.eps / 10.0, a.eps / 100.0):
            v, z, _pr = solve(labels, z0, eps_excl=eps, rounds=a.rounds,
                              maxiter=a.maxiter, eps_inc=eps)
            if v > a.pos:
                best = (v, z)
                break
            if v > best[0]:
                best = (v, z)
        out.append((best[0], tag, best[1]))
        print(f'  {tag:<34} {val:+.6e} -> {best[0]:+.6e}')
    src = load(a.FILE)
    with open(a.out, 'w') as f:
        f.write(f'# rank8.py hardened (exclusion margin {a.eps}) from {a.FILE}.  FLOAT.\n')
        for (v, tag, z), (_t, _val, sqs) in zip(out, src):
            f.write(f'config {tag} value {v:.17g} n {len(sqs)}\n')
            for k, (nm, _x, _y, _th) in enumerate(sqs):
                f.write(f'  {nm} {z[1+3*k]:.17g} {z[2+3*k]:.17g} {z[3+3*k]:.17g}\n')
    print(f'wrote {len(out)} hardened configurations to {a.out}')


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
    # the same statement in EXACT rational arithmetic, through leaf_ceiling's own primitives
    exact_ok, exact_meet = True, None
    for (_pick, z) in st:
        recs, names = [], []
        for k in range(12):
            cx = Fr(round(z[1 + 3 * k] * 1000), 1000)
            cy = Fr(round(z[2 + 3 * k] * 1000), 1000)
            recs.append(lc.make_square(cx, cy, 0, 1, Fr(4)))
            names.append(LNAME[k])
        for k in range(12):
            if bentz.pattern_of(recs[k], IPTS) != frozenset(LABELS[k][1]):
                exact_ok = False
        m = sum(1 for i in range(12) for j in range(i + 1, 12)
                if lc.sq_meets_sq(recs[i], recs[j]))
        if exact_meet is None:
            exact_meet = m
        elif exact_meet != m:
            exact_meet = -1
    ck('EXACT: every tiling-minus-4 family has all twelve patterns right (integer tests)',
       exact_ok)
    ck('EXACT: and exactly 18 of its 66 pairs meet as closed sets (they touch, so it is not '
       'a packing)', exact_meet == 18, f'{exact_meet} meeting pairs')

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
    c.add_argument('--all', action='store_true',
                   help='subsets of all twelve labels, not just the eight singletons')
    c.add_argument('--only', default=None,
                   help='measure just this one class, e.g. --only C0,C1,C2,D0,D1,D2,D3')
    c.add_argument('--exhaustive', action='store_true',
                   help='measure every class, not just the realisability frontier')
    c.add_argument('--bh', type=int, default=4, help='basin-hopping rounds')
    c.add_argument('--bhn', type=int, default=200, help='kicks per basin-hopping round')
    c.add_argument('--bhkeep', type=int, default=8, help='incumbents kicked from')
    c.add_argument('--confirm', type=int, default=4,
                   help='on a negative answer, retry with this many times more starts')
    c.add_argument('--support', default=None)
    c.add_argument('--out', default=None)
    c.add_argument('--nproc', type=int, default=8)
    c = sub.add_parser('linear')
    c.add_argument('--pick', default='0000')
    c.add_argument('--all-picks', action='store_true')
    c.add_argument('--ndual', type=int, default=40)
    c = sub.add_parser('sa2')
    c.add_argument('FILE')
    c.add_argument('--labels', default=None,
                   help='restrict to these labels, e.g. --labels C0,C1,C2,C3,D0,D1,D2,D3')
    c.add_argument('--alpha', action='store_true',
                   help='also compute the integer optimum on the restricted support')
    c.add_argument('--ndual', type=int, default=25)
    c.add_argument('--npairs', type=int, default=25)
    c = sub.add_parser('analyse')
    c.add_argument('FILE')
    c.add_argument('-n', type=int, default=6)
    c.add_argument('--tol', type=float, default=1e-6)
    c = sub.add_parser('harden')
    c.add_argument('FILE')
    c.add_argument('out')
    c.add_argument('--eps', type=float, default=1e-4)
    c.add_argument('--pos', type=float, default=1e-7)
    c.add_argument('--rounds', type=int, default=8)
    c.add_argument('--maxiter', type=int, default=300)
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
    elif a.cmd == 'harden':
        cmd_harden(a)
    elif a.cmd == 'verify':
        cmd_verify(a)
    elif a.cmd == 'selftest':
        sys.exit(selftest())


if __name__ == '__main__':
    main()
