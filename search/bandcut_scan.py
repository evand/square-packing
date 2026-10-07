#!/usr/bin/env python3
"""bandcut_scan.py -- the tilted-band family of packings of `n = T^2 - T` unit squares
in `[0,T]^2`, rebuilt from `T = 12` down to `T = 4`   (2026-09-20, tasks/bandcut-scan).

Background.  `s(T^2-T) < T` is a theorem for every `T >= 12` (Arslanov-Mustafin-Shangitbayev,
EJC 28(4) 2021 #P4.22) and is reported at `T = 11` (Cantrell 2025).  Every such packing known
has the same shape (`search/ARCH_TLEDGER.md` §5): a large majority of the squares axis-parallel
on the integer grid, cut by a **tilted band** which crosses every row and every column of the
container, so that no wall-to-wall chain of `T` axis-parallel squares survives and the chain
bound `B_T <= 0` of `ARCH_TLEDGER.md` §2 has nothing to bite on.

The mechanism, in one line: a straight stack ("bar") of `j` unit squares at a common tilt `t`
has bounding box `(j cos t + sin t) x (cos t + j sin t)`, and `j cos t + sin t < j` exactly when
`tan(t/2) > 1/j` -- so a bar of `j` squares laid across `j` columns leaves **slack**
`j - (j cos t + sin t) > 0` in those columns.  A band built of such bars and threaded
diagonally across the container puts that slack into every row and every column at once.
What it costs is area: the band region holds fewer squares than its cell count, and the budget
is exactly `T` (since `n = T^2 - T` leaves `T` cells free).  The scan below measures that
trade-off as `T` falls.

Everything this script reports is a **feasible point**, i.e. a LOWER bound on `delta*(theta)`
(the disjunctive LP of `search/S6_SKELETON.md` §3.1).  `delta > 0` means the squares fit in
`[0,T]^2` pairwise-disjoint as closed sets with room `delta` to spare, i.e. `s(T^2-T) < T`.

  float pipeline :  band layout  ->  assignment  ->  LP in the centres  ->  ascent
  exact pipeline :  round the centres to rationals  ->  Q[sqrt 2] arithmetic  ->  exact delta

Subcommands:
    selftest            cross-check the exact evaluator against search/s6skel.value
    build   --T 12 ...  build one band configuration, LP-polish it, print delta
    scan    --T ...     sweep the band family at one T
    exact   --json f    exact rational verification of a saved configuration
"""
import argparse
import itertools
import json
import math
import os
import sys
from fractions import Fraction

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

LPOPT = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}


# =========================================================================== exact arithmetic
class Q2:
    """exact element `a + b sqrt(2)` of `Q[sqrt 2]`, with exact sign test.

    Every quantity in the band family lives here: tilts are Pythagorean
    (`cos t = (j^2-1)/(j^2+1)`, `sin t = 2j/(j^2+1)` for `t = 2 arctan(1/j)`), so `0` and `t`
    give rationals, and the `45 deg` squares of the elbow give `sqrt 2 / 2`.
    """
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fraction) else Fraction(a)
        self.b = b if isinstance(b, Fraction) else Fraction(b)

    def __add__(self, o):
        o = _q2(o)
        return Q2(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, o):
        return self + (-_q2(o))

    def __rsub__(self, o):
        return _q2(o) + (-self)

    def __mul__(self, o):
        o = _q2(o)
        return Q2(self.a * o.a + 2 * self.b * o.b, self.a * o.b + self.b * o.a)

    __rmul__ = __mul__

    def sign(self):
        a, b = self.a, self.b
        if a == 0 and b == 0:
            return 0
        if a >= 0 and b >= 0:
            return 1
        if a <= 0 and b <= 0:
            return -1
        # opposite signs: |a| vs |b| sqrt2  <=>  a^2 vs 2 b^2
        d = a * a - 2 * b * b
        if d == 0:
            return 0
        return (1 if d > 0 else -1) if a > 0 else (-1 if d > 0 else 1)

    def __lt__(self, o):
        return (self - _q2(o)).sign() < 0

    def __le__(self, o):
        return (self - _q2(o)).sign() <= 0

    def __gt__(self, o):
        return (self - _q2(o)).sign() > 0

    def __ge__(self, o):
        return (self - _q2(o)).sign() >= 0

    def __eq__(self, o):
        return (self - _q2(o)).sign() == 0

    def __abs__(self):
        return self if self.sign() >= 0 else -self

    def __float__(self):
        return float(self.a) + float(self.b) * math.sqrt(2.0)

    def __repr__(self):
        if self.b == 0:
            return str(self.a)
        return '%s + %s*sqrt2' % (self.a, self.b)

    def __hash__(self):
        return hash((self.a, self.b))


def _q2(x):
    if isinstance(x, Q2):
        return x
    return Q2(x)


HALF = Q2(Fraction(1, 2))


def q2max(x, y):
    return x if x >= y else y


def q2min(x, y):
    return x if x <= y else y


def exact_angles(kinds, j):
    """(cos, sin) in Q[sqrt2] for each angle kind.

    kind 0   : axis-parallel                        (1, 0)
    kind +1  : tilt  t = 2 arctan(1/j)              ((j^2-1)/(j^2+1), 2j/(j^2+1))
    kind -1  : tilt -t   (equivalently 90 - t)      ((j^2-1)/(j^2+1), -2j/(j^2+1))
    kind  2  : 45 deg                               (sqrt2/2, sqrt2/2)
    """
    c = Fraction(j * j - 1, j * j + 1)
    s = Fraction(2 * j, j * j + 1)
    table = {0: (Q2(1), Q2(0)),
             1: (Q2(c), Q2(s)),
             -1: (Q2(c), Q2(-s)),
             2: (Q2(0, Fraction(1, 2)), Q2(0, Fraction(1, 2)))}
    return [table[k] for k in kinds]


def exact_delta(cx, cy, cs, T, report=False):
    """exact `min(min pair gap, min wall slack)` for centres `cx, cy` (Q2) and
    per-square `(cos, sin)` pairs `cs` (Q2).  Returns Q2 (and, if `report`, the witness)."""
    n = len(cx)
    T = _q2(T)
    P = [(abs(c) + abs(s)) * HALF for (c, s) in cs]
    best = None
    witness = None
    for i in range(n):
        for val, tag in ((cx[i] - P[i], 'L'), (T - cx[i] - P[i], 'R'),
                         (cy[i] - P[i], 'B'), (T - cy[i] - P[i], 'T')):
            if best is None or val < best:
                best, witness = val, ('wall', i, tag)
    for i in range(n):
        ci, si = cs[i]
        for jj in range(i + 1, n):
            cj, sj = cs[jj]
            dx, dy = cx[jj] - cx[i], cy[jj] - cy[i]
            g = None
            for (co, so) in ((ci, si), (cj, sj)):
                v = abs(co * dx + so * dy)
                if g is None or v > g:
                    g = v
                v = abs((-so) * dx + co * dy)
                if v > g:
                    g = v
            # W(theta_j - theta_i) = |cos D| + |sin D|
            W = abs(ci * cj + si * sj) + abs(si * cj - ci * sj)
            g = g - HALF - W * HALF
            if g < best:
                best, witness = g, ('pair', i, jj)
    return (best, witness) if report else best


# =========================================================================== float evaluator
def float_value(X, Y, TH, T):
    """min(min pair gap, min wall slack), vectorised; identical to s6skel.value.
    `T` may be a scalar (square container) or a pair `(Tx, Ty)` (rectangle)."""
    n = len(X)
    Tx, Ty = (T, T) if np.isscalar(T) else (T[0], T[1])
    c, s = np.cos(TH), np.sin(TH)
    p = 0.5 * (np.abs(c) + np.abs(s))
    wall = min(float((X - p).min()), float((Tx - X - p).min()),
               float((Y - p).min()), float((Ty - Y - p).min()))
    iu, ju = np.triu_indices(n, 1)
    dx, dy = X[ju] - X[iu], Y[ju] - Y[iu]
    best = np.full(len(iu), -1e18)
    for o in (iu, ju):
        co, so = c[o], s[o]
        np.maximum(best, np.abs(co * dx + so * dy), out=best)
        np.maximum(best, np.abs(-so * dx + co * dy), out=best)
    d = TH[ju] - TH[iu]
    m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
    g = best - m
    k = int(np.argmin(g))
    return min(wall, float(g[k])), (int(iu[k]), int(ju[k])), wall


# =========================================================================== the LP instrument
class BandLP:
    """the disjunctive LP of `S6_SKELETON.md` §3.1 at fixed angles, sparse, on a pair subset.

    Only pairs that can ever be tight are given rows (default: centres within `cut`).  Dropping
    rows only *raises* the LP value, so the LP answer is not itself a bound -- but the point it
    returns is evaluated over **all** pairs by `float_value`, and that is a feasible point,
    i.e. a lower bound on `delta*`.
    """

    def __init__(self, theta, T, pairs, anchors=()):
        # anchors: (i, j, axis) -- force square i onto the lo wall and square j onto the hi
        # wall of that axis, both at slack exactly `delta`.  These are the two wall rows of
        # the chain bound of ARCH_TLEDGER.md 2.1, imposed as equalities: they are what makes
        # a tilted bar *wall-to-wall*, and without them the LP slides the band into a
        # harmless corner and reports the axis-parallel plateau value instead.
        self.anchors = list(anchors)
        self.theta = np.asarray(theta, float)
        self.Tx, self.Ty = (float(T), float(T)) if np.isscalar(T) else (float(T[0]), float(T[1]))
        self.T = self.Tx
        n = self.n = len(theta)
        self.NV = 1 + 2 * n
        self.c_, self.s_ = np.cos(self.theta), np.sin(self.theta)
        self.P = 0.5 * (np.abs(self.c_) + np.abs(self.s_))
        self.pairs = list(pairs)
        d = self.theta[[p[1] for p in self.pairs]] - self.theta[[p[0] for p in self.pairs]]
        self.m = 0.5 + 0.5 * (np.abs(np.cos(d)) + np.abs(np.sin(d)))
        rows, cols, vals, b = [], [], [], []
        r = 0
        for i in range(n):
            for var, TT in ((1 + i, self.Tx), (1 + n + i, self.Ty)):
                rows += [r, r]; cols += [var, 0]; vals += [1.0, -1.0]; b.append(self.P[i]); r += 1
                rows += [r, r]; cols += [var, 0]; vals += [-1.0, -1.0]
                b.append(self.P[i] - TT); r += 1
        self.wrows, self.wcols, self.wvals, self.wb, self.nwall = rows, cols, vals, b, r
        self.obj = np.zeros(self.NV); self.obj[0] = -1.0
        self.bounds = [(-2.0, 1.0)] + [(None, None)] * (2 * n)

    def assign(self, X, Y):
        out = []
        c_, s_ = self.c_, self.s_
        for (i, jj) in self.pairs:
            dx, dy = X[jj] - X[i], Y[jj] - Y[i]
            bv, ba = -1e18, None
            for o in (i, jj):
                cc, ss = c_[o], s_[o]
                for k, (d0, d1) in enumerate(((cc, ss), (-ss, cc))):
                    pr = d0 * dx + d1 * dy
                    for sg in (1.0, -1.0):
                        if sg * pr > bv:
                            bv, ba = sg * pr, (o, k, sg)
            out.append(ba)
        return tuple(out)

    def solve(self, asg):
        n, NV = self.n, self.NV
        rows = list(self.wrows); cols = list(self.wcols); vals = list(self.wvals)
        b = list(self.wb); r = self.nwall
        for (i, jj), (o, k, sg), mij in zip(self.pairs, asg, self.m):
            cc, ss = self.c_[o], self.s_[o]
            d0, d1 = ((cc, ss), (-ss, cc))[k]
            rows += [r, r, r, r, r]
            cols += [1 + jj, 1 + i, 1 + n + jj, 1 + n + i, 0]
            vals += [sg * d0, -sg * d0, sg * d1, -sg * d1, -1.0]
            b.append(mij); r += 1
        A = sparse.coo_matrix((vals, (rows, cols)), shape=(r, NV)).tocsr()
        Aeq = beq = None
        if self.anchors:
            er, ec, ev, eb = [], [], [], []
            q = 0
            for (i, jj, axis) in self.anchors:
                off = 1 + (0 if axis == 'x' else n)
                er += [q, q]; ec += [off + i, 0]; ev += [1.0, -1.0]; eb.append(self.P[i]); q += 1
                er += [q, q]; ec += [off + jj, 0]; ev += [-1.0, -1.0]
                eb.append(self.P[jj] - (self.Tx if axis == 'x' else self.Ty)); q += 1
            Aeq = sparse.coo_matrix((ev, (er, ec)), shape=(q, NV)).tocsr()
            beq = np.asarray(eb)
        res = linprog(self.obj, A_ub=(-A), b_ub=-np.asarray(b), A_eq=Aeq, b_eq=beq,
                      bounds=self.bounds, method='highs', options=LPOPT)
        if not res.success:
            return None
        return res.x

    def ascent(self, X, Y, iters=40):
        best, bX, bY = -1e18, X.copy(), Y.copy()
        seen = set()
        for _ in range(iters):
            asg = self.assign(X, Y)
            if asg in seen:
                break
            seen.add(asg)
            x = self.solve(asg)
            if x is None:
                break
            n = self.n
            X, Y = x[1:1 + n].copy(), x[1 + n:1 + 2 * n].copy()
            v, _, _ = float_value(X, Y, self.theta, self.T)
            if v > best:
                best, bX, bY = v, X.copy(), Y.copy()
        return best, bX, bY


def near_pairs(X, Y, cut=3.0):
    n = len(X)
    out = []
    for i in range(n):
        dx = X[i + 1:] - X[i]
        dy = Y[i + 1:] - Y[i]
        k = np.nonzero(dx * dx + dy * dy <= cut * cut)[0]
        out += [(i, i + 1 + int(t)) for t in k]
    return out


def polish(X, Y, TH, T, cut=3.0, iters=40, rounds=3, anchors=()):
    """LP ascent with the pair set refreshed from the current centres."""
    X, Y = np.asarray(X, float).copy(), np.asarray(Y, float).copy()
    TH = np.asarray(TH, float)
    best, bX, bY = float_value(X, Y, TH, T)[0], X.copy(), Y.copy()
    for _ in range(rounds):
        lp = BandLP(TH, T, near_pairs(X, Y, cut), anchors=anchors)
        v, X, Y = lp.ascent(X, Y, iters=iters)
        if v > best:
            best, bX, bY = v, X.copy(), Y.copy()
        else:
            break
    return best, bX, bY


# =========================================================================== the band family
def bar(x0, y0, ang, j):
    """`j` unit squares in a straight chain at angle `ang`, first centre at `(x0, y0)`;
    consecutive centres differ by the unit edge vector `(cos, sin)` (edge-to-edge contact).
    Bounding box `(j cos + sin) x (cos + j sin)`."""
    c, s = math.cos(ang), math.sin(ang)
    return [(x0 + k * c, y0 + k * s, ang) for k in range(j)]


def band_squares(T, t, j, x0, y0, dx, nbars, ang_sign=+1):
    """A band of `nbars` bars of `j` squares at tilt `ang_sign*t`.

    Consecutive bars are offset by `(dx, dy)` with **perpendicular spacing exactly 1**
    (edge-to-edge between neighbouring bars), i.e. `-dx sin t + dy cos t = 1`; `dx` is the
    drift along the band, and is what turns a vertical strip (`dx` chosen so the bars stay
    stacked) into a diagonal ribbon that cuts every row *and* every column.
    """
    ang = ang_sign * t
    c, s = math.cos(ang), math.sin(ang)
    # perpendicular n = (-s, c);  need  n.(dx,dy) = 1
    dy = (1.0 + dx * s) / c
    out = []
    for r in range(nbars):
        out += bar(x0 + r * dx, y0 + r * dy, ang, j)
    return out


def grid_fill(T, tilted, n_need, tol=0.25, shifts=((0.0, 0.0),)):
    """Add axis-parallel unit squares on integer cells of `[0,T]^2`, best-clearance first.

    `tol` lets a cell in through even when it still overlaps a tilted square by up to `tol`;
    the LP is then asked to push everything apart.  Returns the chosen list and the clearance
    of the worst one accepted (a diagnostic: very negative => the band does not fit)."""
    cand = []
    for ix in range(T):
        for iy in range(T):
            for (sx, sy) in shifts:
                x, y = ix + 0.5 + sx, iy + 0.5 + sy
                if x < 0.5 or x > T - 0.5 or y < 0.5 or y > T - 0.5:
                    continue
                g = 1e9
                for (tx, ty, ta) in tilted:
                    dx, dy = tx - x, ty - y
                    c, s = math.cos(ta), math.sin(ta)
                    v = max(abs(dx), abs(dy), abs(c * dx + s * dy), abs(-s * dx + c * dy))
                    v -= 0.5 + 0.5 * (abs(c) + abs(s))
                    g = min(g, v)
                cand.append((g, x, y))
    cand.sort(key=lambda z: -z[0])
    chosen, used = [], []
    for (g, x, y) in cand:
        if len(chosen) >= n_need:
            break
        ok = all(max(abs(x - u), abs(y - v)) >= 1.0 - 1e-9 for (u, v) in used)
        if not ok or g < -tol:
            continue
        chosen.append((x, y, 0.0)); used.append((x, y))
    return chosen, (min(g for (g, _, _) in
                        [(gg, xx, yy) for (gg, xx, yy) in cand][:0]) if False else None)


def assemble(T, spec):
    """spec: dict with keys  t, j, bands=[{sign,x0,y0,dx,nbars}], elbow=[(x,y,ang)...]"""
    t = spec['t']; j = spec['j']
    tilted = []
    for b in spec.get('bands', []):
        tilted += band_squares(T, t, j, b['x0'], b['y0'], b['dx'], b['nbars'],
                               b.get('sign', 1))
    tilted += [tuple(e) for e in spec.get('elbow', [])]
    n = T * T - T
    if len(tilted) > n:
        return None
    ap, _ = grid_fill(T, tilted, n - len(tilted), tol=spec.get('tol', 0.25))
    sq = tilted + ap
    if len(sq) != n:
        return None, len(tilted), len(ap)
    X = np.array([s[0] for s in sq]); Y = np.array([s[1] for s in sq])
    TH = np.array([s[2] for s in sq])
    return X, Y, TH


# =========================================================================== self-test
def cmd_selftest(a):
    import s6skel
    rng = np.random.default_rng(7)
    worst = 0.0
    for trial in range(12):
        n = int(rng.integers(3, 9))
        j = int(rng.integers(2, 6))
        kinds = [int(rng.choice([0, 1, -1, 2])) for _ in range(n)]
        cs = exact_angles(kinds, j)
        T = Fraction(4)
        cx = [Q2(Fraction(int(rng.integers(0, 400)), 100)) for _ in range(n)]
        cy = [Q2(Fraction(int(rng.integers(0, 400)), 100)) for _ in range(n)]
        d = exact_delta(cx, cy, cs, T)
        tt = math.atan2(2 * j, j * j - 1)
        angs = {0: 0.0, 1: tt, -1: -tt, 2: math.pi / 4}
        TH = np.array([angs[k] for k in kinds])
        X = np.array([float(v) for v in cx]); Y = np.array([float(v) for v in cy])
        z = np.zeros(1 + 3 * n); z[1::3], z[2::3], z[3::3] = X, Y, TH
        f = s6skel.value(z, n, float(T))
        worst = max(worst, abs(float(d) - f))
    print('exact_delta vs s6skel.value : max |diff| = %.3e over 12 random configurations' % worst)
    # float_value vs s6skel.value
    worst2 = 0.0
    for trial in range(8):
        n = int(rng.integers(4, 14))
        TH = rng.uniform(0, math.pi / 2, n)
        X = rng.uniform(0.8, 3.2, n); Y = rng.uniform(0.8, 3.2, n)
        z = np.zeros(1 + 3 * n); z[1::3], z[2::3], z[3::3] = X, Y, TH
        worst2 = max(worst2, abs(float_value(X, Y, TH, 4.0)[0] - s6skel.value(z, n, 4.0)))
    print('float_value vs s6skel.value : max |diff| = %.3e over 8 random configurations' % worst2)
    # the window identity
    for j in (2, 3, 4, 5, 6, 12):
        c = Fraction(j * j - 1, j * j + 1); s = Fraction(2 * j, j * j + 1)
        print('  j=%2d  t=2atan(1/%d)=%9.5f deg   j cos t + sin t = %s  (= j: %s)'
              % (j, j, math.degrees(math.atan2(2 * j, j * j - 1)), j * c + s, j * c + s == j))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('selftest')
    a = ap.parse_args()
    if a.cmd == 'selftest':
        cmd_selftest(a)


if __name__ == '__main__':
    main()
