#!/usr/bin/env python3
"""s6local.py -- fixed-angle margin delta*(theta) with STRUCTURED starts (2026-09-20).

Replaces s6skel.fixed_angle_value where it matters: that instrument starts from uniform random
centres, which at n = 12 lands on a sub-optimal family most of the time (it reported
delta*(t,..,t) = -(2/3) t^2 where a family at <= -(1/3) t^2 exists).  Here the starts are the
tiling sub-configurations (choose n of the T^2 cells, assign the angle classes to the cells)
plus jitter -- the jitter is what breaks the ties at diagonal neighbours, i.e. what chooses
between the row and the pinwheel assignments -- followed by continuation in t.

Everything reported is a LOWER bound on delta* (a feasible point of the disjunctive LP of
S6_SKELETON.md section 3.1), evaluated geometrically by s6skel.value.

    python3 search/s6local.py tilt --n 12 --T 4 --pattern uniform --tdeg 0.2,0.4,0.8,1.6
    python3 search/s6local.py tilt --n 6  --T 3 --pattern 0,0,0,1,1,1 --tdeg 1.6,3.2,6.4,12.8
"""
import argparse
import itertools
import math
import os
import sys

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6skel                                                        # noqa: E402

LPOPT = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}


class Fixed:
    """the disjunctive LP at one exact angle vector."""

    def __init__(self, theta, T):
        self.theta, self.T = list(theta), T
        n = self.n = len(theta)
        self.NV = NV = 1 + 2 * n
        self.cs = [(math.cos(t), math.sin(t)) for t in theta]
        self.P = [0.5 * (abs(c) + abs(s)) for (c, s) in self.cs]
        self.pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        self.m = {(i, j): 0.5 + 0.5 * (abs(math.cos(theta[j] - theta[i]))
                                       + abs(math.sin(theta[j] - theta[i])))
                  for (i, j) in self.pairs}
        A, b = [], []
        for i in range(n):
            for var in (1 + i, 1 + n + i):
                r = np.zeros(NV); r[var] = 1.0; r[0] = -1.0; A.append(r); b.append(self.P[i])
                r = np.zeros(NV); r[var] = -1.0; r[0] = -1.0; A.append(r); b.append(self.P[i] - T)
        self.wallA, self.wallb = A, b
        self.bounds = [(-1.0, 0.5)] + [(self.P[i], T - self.P[i]) for i in range(n)] * 2
        self.c = np.zeros(NV); self.c[0] = -1.0

    def assign(self, X, Y):
        """for every pair the most separated of the 8 edge-normal directions."""
        out = []
        for (i, j) in self.pairs:
            dx, dy = X[j] - X[i], Y[j] - Y[i]
            bv, ba = -1e18, None
            for o in (i, j):
                cc, ss = self.cs[o]
                for k, (d0, d1) in enumerate(((cc, ss), (-ss, cc))):
                    pr = d0 * dx + d1 * dy
                    for sg in (1, -1):
                        if sg * pr > bv:
                            bv, ba = sg * pr, (o, k, sg)
            out.append(ba)
        return tuple(out)

    def lp(self, asg):
        n, NV = self.n, self.NV
        A, b = list(self.wallA), list(self.wallb)
        for (i, j), (o, k, sg) in zip(self.pairs, asg):
            cc, ss = self.cs[o]
            d0, d1 = ((cc, ss), (-ss, cc))[k]
            r = np.zeros(NV)
            r[1 + j] = sg * d0; r[1 + i] = -sg * d0
            r[1 + n + j] = sg * d1; r[1 + n + i] = -sg * d1
            r[0] = -1.0
            A.append(r); b.append(self.m[(i, j)])
        res = linprog(self.c, A_ub=-np.array(A), b_ub=-np.array(b), bounds=self.bounds,
                      method='highs', options=LPOPT)
        if not res.success:
            return None
        return res.x

    def value(self, X, Y):
        n = self.n
        z = np.zeros(1 + 3 * n)
        z[1::3], z[2::3], z[3::3] = X, Y, self.theta
        return s6skel.value(z, n, self.T)

    def ascent(self, X, Y, iters=25):
        best, bX, bY, basg, prev = -1.0, None, None, None, None
        for _ in range(iters):
            asg = self.assign(X, Y)
            if asg == prev:
                break
            prev = asg
            x = self.lp(asg)
            if x is None:
                break
            n = self.n
            X, Y = x[1:1 + n].copy(), x[1 + n:1 + 2 * n].copy()
            v = self.value(X, Y)
            if v > best:
                best, bX, bY, basg = v, X.copy(), Y.copy(), asg
        return best, bX, bY, basg


# --------------------------------------------------------------------------- starts
def pattern_vec(pattern, n):
    if pattern == 'uniform':
        return [1.0] * n
    v = [float(s) for s in pattern.split(',')]
    assert len(v) == n, 'pattern length != n'
    return v


def sampled_starts(n, T, rng, limit, jitters, rho):
    """(X, Y) starts for an arbitrary angle vector: a random n-subset of the tiling cells with a
    random dealing of the labels, `limit` times, each with `jitters` jitters."""
    Ti = int(round(T))
    cells = [(a + 0.5, b + 0.5) for a in range(Ti) for b in range(Ti)]
    out = []
    for _ in range(limit):
        sub = rng.choice(len(cells), size=n, replace=False)
        for _j in range(jitters):
            X = np.array([cells[k][0] for k in sub]) + rng.uniform(-rho, rho, n)
            Y = np.array([cells[k][1] for k in sub]) + rng.uniform(-rho, rho, n)
            out.append((X, Y))
    return out


def structured_starts(n, T, mult, rng, limit, jitters, rho):
    """(X, Y) starts: n of the T^2 tiling cells, the angle classes dealt onto the cells in every
    distinct way (sampled down to `limit` cell/label choices), each with `jitters` jitters."""
    Ti = int(round(T))
    cells = [(a + 0.5, b + 0.5) for a in range(Ti) for b in range(Ti)]
    subsets = list(itertools.combinations(range(len(cells)), n))
    classes = sorted(set(mult))
    base = []
    if len(classes) == 1:
        for sub in subsets:
            base.append([cells[k] for k in sub])
    else:
        # labels matter only through the multiset of multipliers: deal class positions
        idx_by_class = {c: [i for i, m in enumerate(mult) if m == c] for c in classes}
        first = idx_by_class[classes[0]]
        rest = [i for i in range(n) if i not in first]
        assert len(classes) == 2, 'two angle classes supported'
        for sub in subsets:
            for pos in itertools.combinations(range(n), len(first)):
                other = [p for p in range(n) if p not in pos]
                cfg = [None] * n
                for lab, p in zip(first, pos):
                    cfg[lab] = cells[sub[p]]
                for lab, p in zip(rest, other):
                    cfg[lab] = cells[sub[p]]
                base.append(cfg)
    if len(base) > limit:
        pick = rng.choice(len(base), size=limit, replace=False)
        base = [base[k] for k in pick]
    out = []
    for cfg in base:
        for _ in range(jitters):
            X = np.array([c[0] for c in cfg]) + rng.uniform(-rho, rho, n)
            Y = np.array([c[1] for c in cfg]) + rng.uniform(-rho, rho, n)
            out.append((X, Y))
    return out


_F = None


def _init(theta, T):
    global _F
    _F = Fixed(theta, T)


def _job(start):
    X, Y = start
    v, bX, bY, asg = _F.ascent(X, Y)
    return v, bX, bY


def best_at(theta, T, starts, nproc, keep=8):
    from multiprocessing import Pool
    with Pool(nproc, initializer=_init, initargs=(theta, T)) as pool:
        res = pool.map(_job, starts, chunksize=max(1, len(starts) // (nproc * 8)))
    res = [r for r in res if r[1] is not None]
    res.sort(key=lambda r: -r[0])
    return res[:keep], [r[0] for r in res]


def grid_print(X, Y, theta, T):
    order = sorted(range(len(X)), key=lambda k: (-round(Y[k] * 2) / 2, X[k]))
    for k in order:
        print(f"      sq{k:2d}  x={X[k]:.6f}  y={Y[k]:.6f}  tilt={math.degrees(theta[k]):7.3f}")


def cmd_tilt(a):
    rng = np.random.default_rng(a.seed)
    mult = pattern_vec(a.pattern, a.n)
    nz = sum(1 for m in mult if m != 0)
    tdegs = [float(s) for s in a.tdeg.split(',')]
    base = [math.radians(float(s)) for s in a.base.split(',')] if a.base else [0.0] * a.n
    assert len(base) == a.n
    print(f"# n={a.n} T={a.T} pattern={a.pattern} base={a.base}  starts: limit={a.limit} x jitters={a.jitters} "
          f"rho={a.rho}, then continuation from the {a.keep} best of the neighbouring t")
    carried = []
    table = {}
    for sweep in ('down', 'up'):
        seq = sorted(tdegs, reverse=(sweep == 'down'))
        carried = []
        for td in seq:
            t = math.radians(td)
            theta = [b0 + m * t for b0, m in zip(base, mult)]
            if sweep != 'down':
                starts = []
            elif a.base:
                starts = sampled_starts(a.n, a.T, rng, a.limit, a.jitters, a.rho)
            else:
                starts = structured_starts(a.n, a.T, mult, rng, a.limit, a.jitters, a.rho)
            starts += [(X, Y) for (_v, X, Y) in carried]
            if td in table:
                starts += [(X, Y) for (_v, X, Y) in table[td][0]]
            top, allv = best_at(theta, a.T, starts, a.nproc, a.keep)
            if td not in table or top[0][0] > table[td][0][0][0]:
                table[td] = (top, allv)
            carried = table[td][0]
    for td in sorted(tdegs):
        top, allv = table[td]
        t = math.radians(td)
        v = top[0][0]
        nbest = sum(1 for u in allv if u >= v - 1e-9)
        tt = t if t > 0 else float('nan')
        print(f"  t={td:7.3f} deg  delta* >= {v:+.12f}   /t={v / tt:+.5f}   /t^2={v / tt ** 2:+.5f}   "
              f"/t^3={v / tt ** 3:+.5f}   (hit by {nbest}/{len(allv)} starts)")
    if a.show:
        td = float(a.show)
        top, _ = table[td]
        t = math.radians(td)
        print(f"# best configuration at t={td} deg (value {top[0][0]:+.3e})")
        grid_print(top[0][1], top[0][2], [b0 + m * t for b0, m in zip(base, mult)], a.T)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('tilt')
    p.add_argument('--n', type=int, required=True)
    p.add_argument('--T', type=float, required=True)
    p.add_argument('--pattern', default='uniform')
    p.add_argument('--tdeg', default='0.2,0.4,0.8,1.6')
    p.add_argument('--limit', type=int, default=2000)
    p.add_argument('--jitters', type=int, default=4)
    p.add_argument('--rho', type=float, default=0.08)
    p.add_argument('--keep', type=int, default=8)
    p.add_argument('--nproc', type=int, default=14)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--show', default=None)
    p.add_argument('--base', default=None, help='base angles in degrees; theta = base + pattern*t')
    p.set_defaults(fn=cmd_tilt)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
