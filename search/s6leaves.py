#!/usr/bin/env python3
"""s6leaves.py -- how many leaves of the fixed-angle decision tree are TIGHT? (2026-09-20)

s6skel.Decider proves delta*(theta0) <= 0 at an exact angle vector by exhausting the pair
disjunctions; every leaf is an LP (a partial assignment) whose value is <= 0.  For a local theorem
around theta0 the leaves split in two:

  * value < 0 strictly: the same leaf stays closed on a neighbourhood of theta0 (LP values are
    Lipschitz in the data), nothing to prove;
  * value = 0 exactly (TIGHT): the leaf needs an exact inequality in theta (s6exact.py is the
    one-parameter version).

This script re-runs the tree with delta free below, records every leaf value, and for the tight
leaves the support of the LP dual, reduced to a label-free signature (which walls, how many pair
rows, along which axes).

    python3 search/s6leaves.py --n 6 --T 3 --theta 0,0,0,0,0,0
    python3 search/s6leaves.py --n 6 --T 3 --theta 0,0,0,0,20,33
"""
import argparse
import collections
import math
import os
import sys

import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6skel                                                        # noqa: E402

OPT = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10}


class LeafDecider(s6skel.Decider):
    def _lp_free(self, g, rows):
        n, T = self.n, self.T
        NV = 1 + 2 * n
        A = np.array([r for (r, _m) in rows])
        b = np.array([m for (_r, m) in rows])
        bounds = [(-1.0, self.dcap)] + [(g['P'][i], T - g['P'][i]) for i in range(n)] * 2
        c = np.zeros(NV); c[0] = -1.0
        res = linprog(c, A_ub=-A, b_ub=-b, bounds=bounds, method='highs', options=OPT)
        if not res.success:
            return -1.0, None, None
        return -res.fun, res.x, res

    def leaves(self, theta, tol=1e-9):
        g = self.prep(theta, theta)
        nwall = len(self._walls(g))
        stack = [(self._walls(g), frozenset())]
        nodes, out = 0, []
        while stack:
            rows, asg = stack.pop()
            nodes += 1
            if nodes > self.node_cap:
                return None, nodes
            d, x, res = self._lp_free(g, rows)
            if d <= tol:
                sig = None
                if d >= -tol:
                    sig = self.signature(rows, res, nwall)
                out.append((d, len(asg), sig))
                continue
            worst, arg = self._worst(g, x, asg)
            if arg is None or worst >= d - 1e-12:
                return ('SURVIVE', d, x), nodes
            opts = (self._pair_rows(g, arg[1], arg[2]) if arg[0] == 'p'
                    else self._sub_rows(g, arg[1], arg[2], arg[3]))
            na = asg | {arg}
            for (r, m) in opts:
                stack.append((rows + [(r, m)], na))
        return out, nodes

    def signature(self, rows, res, nwall):
        """label-free shape of the dual support: the connected chains it is made of."""
        n = self.n
        y = -res.ineqlin.marginals
        lo = np.abs(res.lower.marginals[1:]); up = np.abs(res.upper.marginals[1:])
        sup_pairs = []
        for k in range(nwall, len(rows)):
            if y[k] > 1e-9:
                r = rows[k][0]
                xs = np.nonzero(np.abs(r[1:1 + n]) > 1e-9)[0]
                ys = np.nonzero(np.abs(r[1 + n:]) > 1e-9)[0]
                ax = ('x' if len(xs) else '') + ('y' if len(ys) else '')
                sq = tuple(sorted(set(xs) | set(ys)))
                sup_pairs.append((ax, sq))
        walls = []
        for k in range(nwall):
            if y[k] > 1e-9:
                r = rows[k][0]
                v = int(np.nonzero(np.abs(r[1:]) > 1e-9)[0][0])
                walls.append(('x' if v < n else 'y'))
        for v in range(2 * n):
            if lo[v] > 1e-9 or up[v] > 1e-9:
                walls.append('x' if v < n else 'y')
        nx = sum(1 for (ax, _s) in sup_pairs if ax == 'x')
        ny = sum(1 for (ax, _s) in sup_pairs if ax == 'y')
        nxy = sum(1 for (ax, _s) in sup_pairs if ax == 'xy')
        squares = set()
        for (_a, sq) in sup_pairs:
            squares |= set(sq)
        return (f"walls x{walls.count('x')} y{walls.count('y')} | pair rows x{nx} y{ny} tilted{nxy}"
                f" | squares {len(squares)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--T', type=float, required=True)
    ap.add_argument('--theta', required=True, help='degrees, comma separated')
    ap.add_argument('--kills', default='N')
    ap.add_argument('--node-cap', type=int, default=400000)
    a = ap.parse_args()
    theta = [math.radians(float(s)) for s in a.theta.split(',')]
    D = LeafDecider(a.n, a.T, kills=a.kills, node_cap=a.node_cap)
    out, nodes = D.leaves(theta)
    print(f"# n={a.n} T={a.T} theta={a.theta} kills={a.kills}: {nodes} LP nodes")
    if out is None:
        print("# node budget exhausted"); return
    if isinstance(out, tuple):
        print(f"# SURVIVES with delta = {out[1]:+.3e} (not in the closed region)"); return
    vals = np.array([d for (d, _k, _s) in out])
    tight = [o for o in out if o[2] is not None]
    neg = vals[vals < -1e-9]
    print(f"# leaves: {len(out)}   tight (|value| <= 1e-9): {len(tight)}   strictly negative: {len(neg)}")
    if len(neg):
        qs = np.quantile(-neg, [0, 0.01, 0.1, 0.5])
        print(f"# strictly negative leaves: smallest |value| = {qs[0]:.3e}, 1% = {qs[1]:.3e}, "
              f"10% = {qs[2]:.3e}, median = {qs[3]:.3e}   (infeasible leaves count as 1)")
        hist = collections.Counter(int(math.floor(math.log10(v))) for v in -neg)
        print("# decades of |value|:", dict(sorted(hist.items())))
    depth = collections.Counter(k for (_d, k, _s) in tight)
    print("# tight leaves by number of branched pairs:", dict(sorted(depth.items())))
    sigs = collections.Counter(s for (_d, _k, s) in tight)
    print(f"# distinct tight dual signatures: {len(sigs)}")
    for s, c in sigs.most_common(20):
        print(f"    {c:6d}  {s}")


if __name__ == '__main__':
    main()
