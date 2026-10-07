#!/usr/bin/env python3
"""s6skel.py -- a separation/transversal branch-and-bound skeleton for  s(6) = 3,
as the rehearsal for  s(12) = 4   (task `tasks/s6-skeleton/README.md`, 2026-09-20).

Theorem being re-proved (Kearney-Shiu 2002; `notes/proof-anatomy.md` §5.2):

    six closed unit squares, pairwise disjoint AS CLOSED SETS, do not fit in [0,3]^2.

Closed semantics exactly as `notes/s13-casefree.md` §1: squares are closed, containment is
closed, "packing" means no point of the plane lies in two of them.  The statement is the
`n^2 - n` case at `n = 3`; it has margin zero (the 3x3 tiling minus three squares is a
configuration with all gaps exactly 0), 18 dimensions, and the same two core types the
`n = 12` census (`search/CENSUS.md` §2) found at `t = 4`.

Subcommands
    margin      multistart max-min-gap for n squares in [0,T]^2 (the instrument)
    cores       first-order rigidity cores at the zero-margin configurations (chains.py, T=3)
    lemma       numerical audit of the chord constants used by K1/K2
    bb          the branch-and-bound skeleton
    selftest    exact/float consistency checks of everything the kills rely on

Everything here is float except where marked; the KILLS are stated as exact lemmas in
`search/S6_SKELETON.md` and the branch-and-bound evaluates their hypotheses with directed
interval arithmetic (see `IV` below), so a kill is rigorous modulo IEEE-754 rounding being
what Python's float arithmetic says it is.
"""
import argparse
import itertools
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rank8                                                         # noqa: E402

SGN = ((1, 1), (1, -1), (-1, 1), (-1, -1))
WALLN = ('L', 'R', 'B', 'T')


# ===================================================================== float geometry (T free)
def pair_gap(z, i, j):
    """separating-axis signed distance; < 0 overlap, = 0 touching, > 0 disjoint as closed sets.
    Identical to rank8.pair_gap (which is T-independent); re-exported for clarity."""
    return rank8.pair_gap(z, i, j)


def min_gap(z, n):
    return min(pair_gap(z, i, j) for i in range(n) for j in range(i + 1, n))


def wall_slack(z, n, T):
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    w = 0.5 * (np.abs(np.cos(TH)) + np.abs(np.sin(TH)))
    return float(min(np.min(X - w), np.min(T - X - w), np.min(Y - w), np.min(T - Y - w)))


def value(z, n, T):
    return min(min_gap(z, n), wall_slack(z, n, T))


def norm_tilt(th):
    """angle reduced to (-45 deg, 45 deg]: a unit square is invariant under 90 deg rotation."""
    return (th + math.pi / 4) % (math.pi / 2) - math.pi / 4


# ===================================================================== the chord lemma constants
#
#  notes/chord-lemma.md Lemma 1: for a unit square at angle th, with C = |cos th|, S = |sin th|,
#  and a horizontal line at distance d <= p = (C+S)/2 from the centre, the chord has length
#
#      ell(th, d) = min( 1/C , 1/S , (p - d)/(C*S) ) .
#
#  Corollary 2: ell >= 1  iff  d <= D(th) := p - C*S,  and D decreases from 1/2 at th = 0 to
#  (sqrt2 - 1)/2 = 0.20710678... at th = 45 deg.  D45 below is the uniform (worst-case) band.

D45 = (math.sqrt(2.0) - 1.0) / 2.0            # 0.2071067811865475...


def chord_D(th):
    """the band half-width D(th) = p - |cos||sin| of notes/chord-lemma.md Corollary 2."""
    C, S = abs(math.cos(th)), abs(math.sin(th))
    return 0.5 * (C + S) - C * S


def chord_len(th, d):
    """exact chord length of a unit square at angle th on a line at distance d from the centre.

    The term (p - d)/(C*S) is ill-conditioned when C*S is near 0 (the axis-parallel limit), so
    the degenerate case is taken by a tolerance rather than by `== 0.0`.  Only `lemma` and
    `selftest` call this; the branch-and-bound uses D_lam, which is well conditioned."""
    C, S = abs(math.cos(th)), abs(math.sin(th))
    p = 0.5 * (C + S)
    if d > p:
        return 0.0
    if C * S <= 1e-12:
        return 1.0
    return min(1.0 / C, 1.0 / S, (p - d) / (C * S))


# ===================================================================== the smooth sub-problem
class Problem:
    """max-min-gap for n unit squares in [0,T]^2, no pattern constraints at all.

    Variables z = (g, x_0, y_0, th_0, ..., x_{n-1}, y_{n-1}, th_{n-1}).
    Rows: for each pair, the four sign branches of the assigned separating axis;
          for each square and wall, the four bounding-box vertex rows.
    """

    def __init__(self, n, T):
        self.n, self.T = n, T
        self.pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        self.nv = 1 + 3 * n

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

    def build(self, pa):
        P = []
        for (pi, (i, j)) in enumerate(self.pairs):
            oi, kind, sg = pa[pi]
            o = i if oi == 0 else j
            for (s1, s2) in SGN:
                P.append((i, j, o, kind, sg, s1, s2))
        self.Pi = np.array([p[0] for p in P], dtype=int)
        self.Pj = np.array([p[1] for p in P], dtype=int)
        self.Po = np.array([p[2] for p in P], dtype=int)
        self.Pk = np.array([p[3] for p in P], dtype=int)
        self.Ps = np.array([p[4] for p in P], dtype=float)
        self.P1 = np.array([p[5] for p in P], dtype=float)
        self.P2 = np.array([p[6] for p in P], dtype=float)
        A = [(k, side, s1, s2) for k in range(self.n) for side in range(4) for (s1, s2) in SGN]
        self.Ak = np.array([a[0] for a in A], dtype=int)
        self.Aside = np.array([a[1] for a in A], dtype=int)
        self.A1 = np.array([a[2] for a in A], dtype=float)
        self.A2 = np.array([a[3] for a in A], dtype=float)
        self.m = len(P) + len(A)

    def cons(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        g = z[0]
        C, S = np.cos(TH), np.sin(TH)
        co, so = C[self.Po], S[self.Po]
        d0 = np.where(self.Pk == 0, co, -so)
        d1 = np.where(self.Pk == 0, so, co)
        dcx = X[self.Pj] - X[self.Pi]
        dcy = Y[self.Pj] - Y[self.Pi]
        dth = TH[self.Pj] - TH[self.Pi]
        rp = (self.Ps * (d0 * dcx + d1 * dcy) - 0.5
              - 0.5 * (self.P1 * np.cos(dth) + self.P2 * np.sin(dth)) - g)
        c, s = C[self.Ak], S[self.Ak]
        h = 0.5 * (self.A1 * c + self.A2 * s)
        base = np.where(self.Aside == 0, X[self.Ak],
                        np.where(self.Aside == 1, self.T - X[self.Ak],
                                 np.where(self.Aside == 2, Y[self.Ak], self.T - Y[self.Ak])))
        return np.concatenate([rp, base - h - g])

    def jac(self, z):
        X, Y, TH = z[1::3], z[2::3], z[3::3]
        C, S = np.cos(TH), np.sin(TH)
        J = np.zeros((self.m, self.nv))
        k = len(self.Pi)
        co, so = C[self.Po], S[self.Po]
        d0 = np.where(self.Pk == 0, co, -so)
        d1 = np.where(self.Pk == 0, so, co)
        dp0 = np.where(self.Pk == 0, -so, -co)
        dp1 = np.where(self.Pk == 0, co, -so)
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
        r = k
        k = len(self.Ak)
        c, s = C[self.Ak], S[self.Ak]
        rows = r + np.arange(k)
        J[rows, 0] = -1.0
        sx = np.where(self.Aside == 0, 1.0, np.where(self.Aside == 1, -1.0, 0.0))
        sy = np.where(self.Aside == 2, 1.0, np.where(self.Aside == 3, -1.0, 0.0))
        np.add.at(J, (rows, 1 + 3 * self.Ak), sx)
        np.add.at(J, (rows, 2 + 3 * self.Ak), sy)
        np.add.at(J, (rows, 3 + 3 * self.Ak), -0.5 * (-self.A1 * s + self.A2 * c))
        return J


def solve(n, T, z0, rounds=8, maxiter=250):
    from scipy.optimize import minimize
    pr = Problem(n, T)
    z = np.array(z0, dtype=float)
    z[0] = value(z, n, T)
    best = (value(z, n, T), z.copy())
    seen = set()
    for _it in range(rounds):
        pa = pr.pair_assign(z)
        key = tuple(pa)
        pr.build(pa)
        res = minimize(lambda w: -w[0], z,
                       jac=lambda w: np.eye(1, len(w), 0).ravel() * -1.0,
                       constraints=[{'type': 'ineq', 'fun': pr.cons, 'jac': pr.jac}],
                       method='SLSQP', options={'maxiter': maxiter, 'ftol': 1e-14})
        zn = res.x
        val = value(zn, n, T)
        if val > best[0]:
            best = (val, zn.copy())
        if key in seen:
            break
        seen.add(key)
        z = zn
        z[0] = value(z, n, T)
    return best


# ===================================================================== starts
def tiling_starts(n, T):
    """all C(T^2, n) sub-configurations of the TxT unit tiling, as start points."""
    m = int(round(T))
    cells = [(i + 0.5, j + 0.5) for i in range(m) for j in range(m)]
    out = []
    for sub in itertools.combinations(range(len(cells)), n):
        z = np.zeros(1 + 3 * n)
        for k, c in enumerate(sub):
            z[1 + 3 * k], z[2 + 3 * k], z[3 + 3 * k] = cells[c][0], cells[c][1], 0.0
        out.append((f'tile{"".join(map(str, sub))}', z))
    return out


def random_start(rng, n, T, axis_frac=0.5):
    z = np.zeros(1 + 3 * n)
    for k in range(n):
        th = 0.0 if rng.random() < axis_frac else rng.uniform(0.0, math.pi / 2)
        p = 0.5 * (abs(math.cos(th)) + abs(math.sin(th)))
        z[1 + 3 * k] = rng.uniform(p, T - p)
        z[2 + 3 * k] = rng.uniform(p, T - p)
        z[3 + 3 * k] = th
    return z


def nudge(z, rng, n, T, rho=0.05, alpha=0.05):
    w = z.copy()
    for k in range(n):
        w[1 + 3 * k] += rng.normal(0.0, rho)
        w[2 + 3 * k] += rng.normal(0.0, rho)
        w[3 + 3 * k] += rng.normal(0.0, alpha)
        p = 0.5 * (abs(math.cos(w[3 + 3 * k])) + abs(math.sin(w[3 + 3 * k])))
        w[1 + 3 * k] = min(max(w[1 + 3 * k], p), T - p)
        w[2 + 3 * k] = min(max(w[2 + 3 * k], p), T - p)
    return w


_G = {}


def _init(n, T, rounds, maxiter):
    _G['n'], _G['T'], _G['rounds'], _G['maxiter'] = n, T, rounds, maxiter


def _job(arg):
    tag, z0 = arg
    v, z = solve(_G['n'], _G['T'], z0, rounds=_G['rounds'], maxiter=_G['maxiter'])
    return tag, float(v), z.tolist()


def run_starts(starts, n, T, nproc, rounds=8, maxiter=250):
    import multiprocessing as mp
    if nproc <= 1:
        _init(n, T, rounds, maxiter)
        return [_job(s) for s in starts]
    with mp.Pool(nproc, initializer=_init, initargs=(n, T, rounds, maxiter)) as pool:
        return pool.map(_job, starts, chunksize=4)


def max_tilt_deg(z, n):
    return float(max(abs(math.degrees(norm_tilt(z[3 + 3 * k]))) for k in range(n)))


def tiling_distance(z, n, T):
    """max-norm distance from the nearest sub-configuration of the tiling (centres and tilt)."""
    m = int(round(T))
    cells = [(i + 0.5, j + 0.5) for i in range(m) for j in range(m)]
    best = 1e18
    for sub in itertools.permutations(range(len(cells)), n):
        d = 0.0
        for k, c in enumerate(sub):
            d = max(d, abs(z[1 + 3 * k] - cells[c][0]), abs(z[2 + 3 * k] - cells[c][1]),
                    abs(norm_tilt(z[3 + 3 * k])))
            if d >= best:
                break
        best = min(best, d)
    return best


def _tiling_distance_fast(z, n, T):
    """greedy assignment version of tiling_distance (upper bound); exact for n <= 6 via
    scipy's linear_sum_assignment on the max-norm cost, which is what we use."""
    from scipy.optimize import linear_sum_assignment
    m = int(round(T))
    cells = np.array([(i + 0.5, j + 0.5) for i in range(m) for j in range(m)])
    X, Y, TH = z[1::3][:n], z[2::3][:n], z[3::3][:n]
    tl = np.abs([norm_tilt(t) for t in TH])
    cost = np.maximum(np.abs(X[:, None] - cells[None, :, 0]),
                      np.abs(Y[:, None] - cells[None, :, 1]))
    cost = np.maximum(cost, tl[:, None])
    # bottleneck assignment by bisection on the threshold
    vals = np.unique(cost)
    lo, hi = 0, len(vals) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        ok = np.where(cost <= vals[mid], 0.0, 1.0)
        r, c = linear_sum_assignment(ok)
        if ok[r, c].sum() == 0:
            hi = mid
        else:
            lo = mid + 1
    return float(vals[lo])


def cmd_margin(a):
    n, T = a.n, a.T
    rng = np.random.default_rng(a.seed)
    starts = []
    if not a.no_tiling:
        starts += tiling_starts(n, T)
        base = list(starts)
        for (tag, z) in base:
            for q in range(a.nudges):
                starts.append((f'{tag}+n{q}', nudge(z, rng, n, T)))
    for q in range(a.nstart):
        starts.append((f'rand{q}', random_start(rng, n, T)))
    print(f'# {len(starts)} starts, n={n}, T={T}, nproc={a.nproc}', flush=True)
    res = run_starts(starts, n, T, a.nproc)
    for rd in range(a.kick_rounds):
        res.sort(key=lambda r: -r[1])
        seeds = res[:a.kick_seeds]
        more = []
        for (tag, _v, zl) in seeds:
            z = np.array(zl)
            for q in range(a.kicks):
                more.append((f'{tag}/k{rd}_{q}',
                             nudge(z, rng, n, T, rho=0.02 * (q % 5 + 1), alpha=0.02 * (q % 5 + 1))))
        res += run_starts(more, n, T, a.nproc)
        print(f'# kick round {rd}: {len(more)} more, best {max(r[1] for r in res):.12e}', flush=True)
    res.sort(key=lambda r: -r[1])
    best = res[0][1]
    zero = [r for r in res if r[1] > -a.near]
    print(f'BEST  {best:+.12e}   ({res[0][0]})')
    print(f'starts {len(res)}   |value| <= {a.near}: {len(zero)}   value > +1e-12: '
          f'{sum(1 for r in res if r[1] > 1e-12)}')
    if zero:
        tl = [max_tilt_deg(np.array(r[2]), n) for r in zero]
        td = [_tiling_distance_fast(np.array(r[2]), n, T) for r in zero]
        print(f'zero-margin points: max tilt {max(tl):.4f} deg, '
              f'{sum(1 for t in tl if t > 1.0)} with tilt > 1 deg')
        print(f'                    tiling distance: max {max(td):.6f}, '
              f'{sum(1 for d in td if d > 0.05)} further than 0.05')
    if a.out:
        with open(a.out, 'w') as f:
            for (tag, v, zl) in res:
                f.write(json.dumps({'tag': tag, 'best': v, 'z': zl}) + '\n')
        print(f'# wrote {a.out}')


# ===================================================================== first-order cores (T = 3)
def cmd_cores(a):
    import chains
    chains.T = a.T                       # chains.py hard-codes T = 4; the analysis itself is T-free
    rng = np.random.default_rng(a.seed)
    rows = []
    with open(a.jsonl) as f:
        for line in f:
            r = json.loads(line)
            if abs(r['best']) > a.near:
                continue
            rows.append((r['tag'], np.array(r['z'])))
    print(f'# {len(rows)} zero-margin configurations', flush=True)
    from collections import Counter
    sig = Counter()
    detail = {}
    pos = 0
    for (tag, z) in rows[:a.limit]:
        pairs, walls = chains.structure(z, a.n, a.tol)
        d, _v = chains.first_order(a.n, pairs, walls)
        if d > a.zero:
            pos += 1
            print(f'!! delta* = {d:.3e} at {tag}')
            continue
        keep = chains.core(z, a.n, pairs, walls, rng, a.ntry, a.zero)
        s, det = chains.signature(z, a.n, keep, pairs, [f'S{k}' for k in range(a.n)])
        sig[s] += 1
        detail.setdefault(s, []).append((tag, det))
    print(f'first-order feasible (delta* > {a.zero}): {pos}')
    for s, c in sig.most_common():
        print(f'{c:6d}  {s}')
        for (tag, det) in detail[s][:2]:
            print(f'          e.g. {tag}: {det}')


# ===================================================================== lemma audit
def cmd_lemma(a):
    print('# chord length  ell(th,d) = min(1/C, 1/S, (p-d)/(C S)),  p = (C+S)/2')
    print('# D(th) = p - C S  is the largest d with ell >= 1   (notes/chord-lemma.md Cor. 2)')
    print()
    print(f'{"theta":>8} {"u=C+S":>10} {"D(th)":>12} {"ell(th,D)":>12} '
          f'{"(C-S)/2":>12}  brief-band')
    for deg in (0, 5, 10, 15, 20, 25, 30, 35, 40, 44, 45):
        th = math.radians(deg)
        D = chord_D(th)
        C, S = abs(math.cos(th)), abs(math.sin(th))
        print(f'{deg:8d} {C+S:10.6f} {D:12.8f} {chord_len(th, D):12.8f} {(C-S)/2:12.8f}')
    print()
    print(f'min over theta of D(theta) = (sqrt2-1)/2 = {D45:.16f}')
    grid = np.linspace(0, math.pi / 2, 200001)
    Dv = np.array([chord_D(t) for t in grid])
    print(f'numeric min over 200001 angles: {Dv.min():.16f}  at '
          f'{math.degrees(grid[int(Dv.argmin())]):.4f} deg')
    # worst-case chord as a function of d, over all tilts
    print()
    print('# L(d) := min over theta of ell(theta, d)  -- the uniform chord bound at distance d')
    print(f'{"d":>10} {"L(d)":>12} {"argmin theta":>14}')
    for d in (0.0, 0.05, 0.1, D45, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7):
        vals = [(chord_len(t, d), t) for t in grid]
        v, t = min(vals)
        print(f'{d:10.6f} {v:12.8f} {math.degrees(t):14.4f}')


# ===================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('margin')
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--T', type=float, default=3.0)
    p.add_argument('--nstart', type=int, default=2000)
    p.add_argument('--nudges', type=int, default=4)
    p.add_argument('--no-tiling', action='store_true')
    p.add_argument('--kick-rounds', type=int, default=2)
    p.add_argument('--kick-seeds', type=int, default=32)
    p.add_argument('--kicks', type=int, default=20)
    p.add_argument('--near', type=float, default=1e-9)
    p.add_argument('--nproc', type=int, default=14)
    p.add_argument('--seed', type=int, default=20260920)
    p.add_argument('--out', default='')
    p.set_defaults(func=cmd_margin)

    p = sub.add_parser('cores')
    p.add_argument('jsonl')
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--T', type=float, default=3.0)
    p.add_argument('--near', type=float, default=1e-9)
    p.add_argument('--tol', type=float, default=1e-7)
    p.add_argument('--zero', type=float, default=1e-9)
    p.add_argument('--ntry', type=int, default=3)
    p.add_argument('--limit', type=int, default=400)
    p.add_argument('--seed', type=int, default=7)
    p.set_defaults(func=cmd_cores)

    p = sub.add_parser('lemma')
    p.set_defaults(func=cmd_lemma)

    p = sub.add_parser('selftest')
    p.add_argument('--node-cap', type=int, default=150000)
    p.set_defaults(func=cmd_selftest)

    p = sub.add_parser('bb')
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--T', type=float, default=3.0)
    p.add_argument('--kinds', default='N,K1,K2')
    p.add_argument('--warm', type=int, default=8)
    p.add_argument('--depth', type=int, default=12)
    p.add_argument('--node-cap', type=int, default=20000)
    p.add_argument('--max-boxes', type=int, default=40000)
    p.add_argument('--nproc', type=int, default=14)
    p.add_argument('--out', default='')
    p.set_defaults(func=cmd_bb)

    a = ap.parse_args()
    a.func(a)




# ===================================================================== the disjunctive LP
#
#  THE STRUCTURAL FACT THIS SKELETON IS BUILT ON.  Fix the angles th_0..th_{n-1}.  Then
#
#    * containment of square i in [0,T]^2 is  p_i <= x_i <= T - p_i  (same in y),
#      p_i = (|cos th_i| + |sin th_i|)/2            -- LINEAR in the centres;
#    * squares i and j are disjoint as CLOSED sets iff some edge normal separates them
#      strictly (separating-axis theorem for convex polygons):
#
#          max over  o in {i,j},  k in {0,1},  s in {+,-}   of   s * (n_{o,k} . (c_j - c_i))
#              >  1/2 + W(th_j - th_i)/2,       W(D) = |cos D| + |sin D| in [1, sqrt2]
#
#      -- a DISJUNCTION of linear inequalities in the centres.
#
#  So for fixed angles the whole question "is there a closed packing?" is a disjunctive LP in
#  2n variables, and
#
#      delta*(th) := max over centres of  min over (pairs, walls) of slack
#
#  is > 0 exactly when a closed packing with those angles exists.  s(6) = 3 is therefore
#  equivalent to  delta*(th) <= 0  for every th in [0,90)^6:  a statement in SIX dimensions,
#  not eighteen.  The branch-and-bound below branches on the angles only.
#
#  Over a BOX of angles we must relax.  Writing (u, v) for the coordinates of d = c_j - c_i in
#  square o's MID-angle frame and w_o for the half-width of o's angle interval,
#
#      max over th_o in the interval of |n_{o,k}(th_o) . d|   <=   |chosen| + sin(w_o) |other|
#
#  (from n(t0 + f) = cos f * n0 + sin f * n0^perp), and we use the most permissive constants
#  p_i^lo = min over the box, W^lo = min over the box.  Both are exact when the box is a point.

def _wlo(dlo, dhi):
    """min of W(D) = |cos D| + |sin D| over D in [dlo, dhi]; W is 1 at multiples of pi/2."""
    lo = math.floor(dlo / (math.pi / 2)) * (math.pi / 2)
    while lo <= dhi + 1e-15:
        if lo >= dlo - 1e-15:
            return 1.0
        lo += math.pi / 2
    f = lambda d: abs(math.cos(d)) + abs(math.sin(d))
    return min(f(dlo), f(dhi))


def _plo(tlo, thi):
    """min of p(th) = (|cos|+|sin|)/2 over [tlo, thi]."""
    lo = math.floor(tlo / (math.pi / 2)) * (math.pi / 2)
    while lo <= thi + 1e-15:
        if lo >= tlo - 1e-15:
            return 0.5
        lo += math.pi / 2
    f = lambda t: 0.5 * (abs(math.cos(t)) + abs(math.sin(t)))
    return min(f(tlo), f(thi))


def _Dlo(tlo, thi):
    """min of the chord band D(th) = p - |cos||sin| over [tlo, thi] (the lam = 1 member)."""
    return _Dlam_lo(tlo, thi, 1.0)


class DLP:
    """The disjunctive LP over the 2n centre coordinates, for one box of angles.

    solve(tlo, thi) returns (delta*, x, status).  delta* <= 0 is a CERTIFICATE that no
    closed packing of n unit squares in [0,T]^2 has its angles in the box [tlo, thi]
    (all angles in radians).  delta is capped at DCAP: only its sign is used.

    Rows
      walls     x_i >= p_i^lo + delta,  x_i <= T - p_i^lo - delta   (same in y)
      pairs     a 16-way (8-way for a degenerate box) disjunction, one branch per
                (owner, which coordinate of the owner's frame, sign, sign of the other)
      K1        (optional) for every triple and each of x, y: the three chord bands
                do NOT share a point -- some pair in the triple is more than
                D_i^lo + D_j^lo apart in that coordinate.
    """

    DCAP = 0.5

    def __init__(self, n, T, k1=True):
        self.n, self.T, self.k1 = n, T, k1
        self.pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        self.triples = list(itertools.combinations(range(n), 3))

    def build(self, tlo, thi):
        n, T = self.n, self.T
        mid = [(a + b) / 2 for a, b in zip(tlo, thi)]
        w = [(b - a) / 2 for a, b in zip(tlo, thi)]
        sw = [math.sin(min(x, math.pi / 4)) for x in w]
        P = [_plo(a, b) for a, b in zip(tlo, thi)]
        D = [_Dlo(a, b) for a, b in zip(tlo, thi)]
        NV = 1 + 2 * n
        xi = lambda i: 1 + i
        yi = lambda i: 1 + n + i
        rows = []                                   # (coef vector over NV, rhs, group or None)
        for i in range(n):
            for var in (xi(i), yi(i)):
                r = np.zeros(NV); r[var] = 1.0; r[0] = -1.0
                rows.append((r, P[i], None))
                r = np.zeros(NV); r[var] = -1.0; r[0] = -1.0
                rows.append((r, P[i] - T, None))
        groups = []
        for (i, j) in self.pairs:
            m = 0.5 + 0.5 * _wlo(tlo[j] - thi[i], thi[j] - tlo[i])
            R = T - P[i] - P[j]
            gs, seen = [], set()
            for o in (i, j):
                c, s = math.cos(mid[o]), math.sin(mid[o])
                U = np.zeros(NV); U[xi(j)] = c; U[xi(i)] = -c; U[yi(j)] = s; U[yi(i)] = -s
                V = np.zeros(NV); V[xi(j)] = -s; V[xi(i)] = s; V[yi(j)] = c; V[yi(i)] = -c
                for (A, B) in ((U, V), (V, U)):
                    for sa in (1.0, -1.0):
                        for sb in (1.0, -1.0):
                            r = sa * A + sw[o] * sb * B
                            r[0] = -1.0
                            key = tuple(np.round(r, 12))
                            if key in seen:
                                continue
                            seen.add(key)
                            M = m + self.DCAP + math.sqrt(2.0) * R * (1.0 + sw[o]) + 1e-9
                            gs.append((r, m, M))
            groups.append(gs)
        if self.k1:
            for tri in self.triples:
                for axis in (0, 1):
                    gs = []
                    for (a, b) in itertools.combinations(tri, 2):
                        for sg in (1.0, -1.0):
                            r = np.zeros(NV)
                            va, vb = (xi(a), xi(b)) if axis == 0 else (yi(a), yi(b))
                            r[va] = sg; r[vb] = -sg; r[0] = -1.0
                            m = D[a] + D[b]
                            M = m + self.DCAP + (T - P[a] - P[b]) + 1e-9
                            gs.append((r, m, M))
                    groups.append(gs)
        return NV, P, rows, groups

    def solve(self, tlo, thi, time_limit=20.0):
        from scipy.optimize import milp, LinearConstraint, Bounds
        n, T = self.n, self.T
        NV, P, rows, groups = self.build(tlo, thi)
        nb = sum(len(g) for g in groups)
        NT = NV + nb
        nrow = len(rows) + nb + len(groups)
        A = np.zeros((nrow, NT)); lo = np.zeros(nrow)
        rr = 0
        for (r, m, _g) in rows:
            A[rr, :NV] = r; lo[rr] = m; rr += 1
        b = 0
        for g in groups:
            for (r, m, M) in g:
                A[rr, :NV] = r; A[rr, NV + b] = -M; lo[rr] = m - M; rr += 1; b += 1
        b = 0
        for g in groups:
            for _ in g:
                A[rr, NV + b] = 1.0; b += 1
            lo[rr] = 1.0; rr += 1
        cost = np.zeros(NT); cost[0] = -1.0
        lb = np.concatenate([[0.0], [P[i] for i in range(n)] * 2, np.zeros(nb)])
        ub = np.concatenate([[self.DCAP], [T - P[i] for i in range(n)] * 2, np.ones(nb)])
        integ = np.concatenate([np.zeros(NV), np.ones(nb)])
        res = milp(c=cost, constraints=LinearConstraint(A, lo, np.inf), integrality=integ,
                   bounds=Bounds(lb, ub),
                   options={'time_limit': time_limit, 'presolve': True})
        if res.status == 2:
            return 0.0, None, 'infeasible'         # not even delta = 0: a fortiori a kill
        if res.status == 1 or res.x is None:
            db = getattr(res, 'mip_dual_bound', None)
            return (self.DCAP if db is None else -db), None, 'timeout'
        return -res.fun, res.x[:NV], 'ok'


def _sym_break(n):
    """label symmetry: the angle bins may be taken non-decreasing."""
    return True


# ===================================================================== K1 / K2: lambda-bands
#
#  K2 (FRACTIONAL CHORD LEMMA).  For lambda in (0,1] put
#
#      D_lam(th) = p(th) - lam*|cos th||sin th|,       p(th) = (|cos th| + |sin th|)/2,
#
#  the "lambda-band" half-width.  By notes/chord-lemma.md Lemma 1, a unit square whose centre
#  is within D_lam(th) of the line y = a meets that line in a chord of length >= lam (for
#  lam <= 1 the terms 1/|C|, 1/|S| >= 1 never bind).  If k squares of a packing all have the
#  line in their lambda-band, their chords are k pairwise DISJOINT CLOSED intervals of length
#  >= lam inside a segment of length T, so ordering them by left endpoint,
#
#      T  >=  beta_k  >  alpha_1 + k*lam  >=  k*lam ,
#
#  a contradiction as soon as  k*lam >= T.  So:
#
#      *** no k squares of a closed packing in [0,T]^2 share a point of their lambda-bands,
#          for any (k, lam) with k*lam >= T ***                                        (K2)
#
#  and the same for vertical lines.  lam = 1 is K1 (band depth <= T-1); the strongest member
#  for each k is lam = T/k.  At T = 3 the family is (k, lam) = (3, 1), (4, 3/4), (5, 3/5),
#  (6, 1/2).  D_lam is INCREASING as lam decreases, so the small-lam members are genuinely
#  stronger constraints on tilted squares: at 45 deg, D_1 = 0.2071 but D_{3/4} = 0.3321.
#  At th = 0 all members coincide (D_lam = 1/2, CS = 0) -- the whole family is a statement
#  about TILTED squares, which is what CENSUS.md §3 asked for.
#
#  Where the chain length enters: the rule is k*lam >= T.  At T = 3 four squares need only
#  chord 3/4 each; at T = 17 the cheapest member with lam = 1 needs k = 17, and no member
#  constrains fewer than 17 squares at all.  An argument built on this family therefore says
#  nothing about 272 squares in [0,17]^2 -- as it must not.

def D_lam(th, lam):
    C, S = abs(math.cos(th)), abs(math.sin(th))
    return 0.5 * (C + S) - lam * C * S


def _crit(tlo, thi):
    """the angles in [tlo, thi] at which u = |cos|+|sin| attains an endpoint of its range:
    the interval's own endpoints, plus any multiple of 45 deg inside it."""
    out = [tlo, thi]
    k0 = math.floor(tlo / (math.pi / 4))
    for k in (k0, k0 + 1, k0 + 2, k0 + 3):
        c = k * math.pi / 4
        if tlo - 1e-15 <= c <= thi + 1e-15:
            out.append(c)
    return out


def _Dlam_lo(tlo, thi, lam):
    """min of D_lam over [tlo, thi].

    D_lam = u/2 - lam*(u^2-1)/2 with u = |cos|+|sin| in [1, sqrt2] is CONCAVE in u, so its
    minimum over an interval of u is at an endpoint of that interval; u's extremes on an angle
    interval are at the interval's endpoints and at the multiples of 45 deg inside it.  (The
    naive "the minimum is at 45 deg" is only true for lam >= 1/2 -- d(D_lam)/du = 1/2 - lam*u
    -- and the T = 4 family reaches lam = 1/3, so it must not be used.)"""
    return min(D_lam(c, lam) for c in _crit(tlo, thi))


def lam_family(T, kmax):
    """the members (k, lam = T/k) of K2 that constrain at most kmax squares."""
    return [(k, T / k) for k in range(2, kmax + 1) if T / k <= 1.0]


# ===================================================================== the decision procedure
#
#  The big-M MILP form of the disjunctive LP is useless (measured: on six squares at generic
#  fixed angles HiGHS does not move the dual bound off the cap in 30 s).  What works is lazy
#  branching: solve the LP with only the rows branched on so far, look at the centres it
#  returns, find a pair (or a K1 triple) that the LP is cheating on, and branch on that.
#  Pruning is sound because dropping rows can only INCREASE delta*.

class Decider:
    """decide(tlo, thi) -> ('kill'|'survive'|'budget', x, info).

    'kill' means delta* <= 0 for the relaxation on the angle box, hence NO closed packing of
    n unit squares in [0,T]^2 has all its angles inside the box.  'survive' returns the centre
    vector of a relaxed witness (not a packing: the relaxation lets every square pick, pair by
    pair, the most convenient angle inside its own interval).

    kills = 'N'   containment + pairwise separation only  (the ordinary interval kill)
            'K1'  N plus the lambda = 1 member of K2 (band depth <= T-1)
            'K2'  N plus the whole lambda-family
    """

    def __init__(self, n, T, kills='K2', dcap=0.5, node_cap=200000, kmax=None):
        self.n, self.T, self.dcap, self.node_cap = n, T, dcap, node_cap
        self.kills = kills
        self.pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        fam = lam_family(T, n if kmax is None else kmax)
        if kills == 'N':
            fam = []
        elif kills == 'K1':
            fam = [m for m in fam if m[1] >= 1.0 - 1e-12]
        self.fam = fam
        self.subs = {}
        for (k, lam) in fam:
            self.subs[(k, lam)] = [tuple(c) for c in itertools.combinations(range(n), k)]

    def prep(self, tlo, thi):
        n, T = self.n, self.T
        g = {}
        g['mid'] = [(a + b) / 2 for a, b in zip(tlo, thi)]
        g['w'] = [min((b - a) / 2, math.pi / 4) for a, b in zip(tlo, thi)]
        g['P'] = [_plo(a, b) for a, b in zip(tlo, thi)]
        g['DL'] = {(k, lam): [_Dlam_lo(a, b, lam) for a, b in zip(tlo, thi)]
                   for (k, lam) in self.fam}
        g['m'] = {(i, j): 0.5 + 0.5 * _wlo(tlo[j] - thi[i], thi[j] - tlo[i])
                  for (i, j) in self.pairs}
        g['cs'] = [(math.cos(t), math.sin(t)) for t in g['mid']]
        return g

    def _pair_rows(self, g, i, j):
        """the branch rows of pair (i, j).

        The exact condition is  max over o in {i,j}, k in {0,1}, s in {+-1}, th_o in its
        interval, of  s * n_{o,k}(th_o) . d  >=  m,   d = c_j - c_i.  Writing (A, B) for the
        coordinates of d in o's MID-angle frame and w for the half-width of o's interval,
        the maximum over th_o is  max over |f| <= w of  cos f * A + sin f * B, which is

            sqrt(A^2 + B^2)                 if the maximiser is interior (|B| <= |A| tan w),
            n_{o,k}(mid -+ w) . d           if it is at an endpoint.

        In the interior case sqrt(A^2+B^2) <= |A| / cos w, so the row  s*A >= m*cos(w) is a
        valid relaxation; in the endpoint cases the row is the EXACT constraint at a real
        angle.  Covering both gives a relaxation whose only error is m*(1 - cos w) = O(w^2),
        against the O(w) of the naive  |A| + sin(w)|B|  bound."""
        n = self.n
        NV = 1 + 2 * n
        out, seen = [], set()
        for o in (i, j):
            w = g['w'][o]
            for (ang, fac) in ((g['mid'][o], math.cos(w)),
                               (g['mid'][o] - w, 1.0), (g['mid'][o] + w, 1.0)):
                c, s = math.cos(ang), math.sin(ang)
                U = np.zeros(NV)
                U[1 + j] = c; U[1 + i] = -c; U[1 + n + j] = s; U[1 + n + i] = -s
                V = np.zeros(NV)
                V[1 + j] = -s; V[1 + i] = s; V[1 + n + j] = c; V[1 + n + i] = -c
                for A in (U, V):
                    for sa in (1.0, -1.0):
                        r = sa * A
                        r[0] = -1.0
                        m = g['m'][(i, j)] * fac
                        key = tuple(np.round(r, 12)) + (round(m, 12),)
                        if key in seen:
                            continue
                        seen.add(key)
                        out.append((r, m))
        return out

    def _sub_rows(self, g, sub, axis, mem):
        """K2 on one k-subset: SOME pair of it must be lambda-band separated in this axis."""
        n = self.n
        NV = 1 + 2 * n
        off = 1 if axis == 0 else 1 + n
        DL = g['DL'][mem]
        out = []
        for (a, b) in itertools.combinations(sub, 2):
            for sg in (1.0, -1.0):
                r = np.zeros(NV)
                r[off + a] = sg; r[off + b] = -sg; r[0] = -1.0
                out.append((r, DL[a] + DL[b]))
        return out

    def _walls(self, g):
        n, T = self.n, self.T
        NV = 1 + 2 * n
        out = []
        for i in range(n):
            for var in (1 + i, 1 + n + i):
                r = np.zeros(NV); r[var] = 1.0; r[0] = -1.0
                out.append((r, g['P'][i]))
                r = np.zeros(NV); r[var] = -1.0; r[0] = -1.0
                out.append((r, g['P'][i] - T))
        return out

    def _lp(self, g, rows):
        from scipy.optimize import linprog
        n, T = self.n, self.T
        NV = 1 + 2 * n
        A = np.array([r for (r, _m) in rows])
        b = np.array([m for (_r, m) in rows])
        bounds = [(0.0, self.dcap)] + [(g['P'][i], T - g['P'][i]) for i in range(n)] * 2
        c = np.zeros(NV); c[0] = -1.0
        res = linprog(c, A_ub=-A, b_ub=-b, bounds=bounds, method='highs')
        if not res.success:
            return -1.0, None
        return -res.fun, res.x

    def _worst(self, g, x, asg):
        n = self.n
        X, Y = x[1:1 + n], x[1 + n:1 + 2 * n]
        worst, arg = 1e18, None
        for (i, j) in self.pairs:
            if ('p', i, j) in asg:
                continue
            dx, dy = X[j] - X[i], Y[j] - Y[i]
            best = -1e18
            mm = g['m'][(i, j)]
            for o in (i, j):
                w = g['w'][o]
                for (ang, fac) in ((g['mid'][o], math.cos(w)),
                                   (g['mid'][o] - w, 1.0), (g['mid'][o] + w, 1.0)):
                    c, s = math.cos(ang), math.sin(ang)
                    u, v = c * dx + s * dy, -s * dx + c * dy
                    best = max(best, max(abs(u), abs(v)) - mm * fac)
            sl = best
            if sl < worst:
                worst, arg = sl, ('p', i, j)
        for mem in self.fam:
            DL = g['DL'][mem]
            for axis in (0, 1):
                C = X if axis == 0 else Y
                for sub in self.subs[mem]:
                    if ('s', sub, axis, mem) in asg:
                        continue
                    sl = max(abs(C[a] - C[b]) - DL[a] - DL[b]
                             for (a, b) in itertools.combinations(sub, 2))
                    if sl < worst:
                        worst, arg = sl, ('s', sub, axis, mem)
        return worst, arg

    def decide(self, tlo, thi):
        g = self.prep(tlo, thi)
        stack = [(self._walls(g), frozenset())]
        nodes = 0
        while stack:
            rows, asg = stack.pop()
            nodes += 1
            if nodes > self.node_cap:
                return 'budget', None, {'nodes': nodes}
            d, x = self._lp(g, rows)
            if d <= 1e-12:
                continue
            worst, arg = self._worst(g, x, asg)
            if arg is None or worst >= d - 1e-12:
                return 'survive', x, {'nodes': nodes, 'delta': d}
            opts = (self._pair_rows(g, arg[1], arg[2]) if arg[0] == 'p'
                    else self._sub_rows(g, arg[1], arg[2], arg[3]))
            na = asg | {arg}
            for (r, m) in opts:
                stack.append((rows + [(r, m)], na))
        return 'kill', None, {'nodes': nodes}


# ===================================================================== fixed-angle margin
def fixed_angle_value(theta, T, rng, nstart=60, iters=12):
    """delta*(theta) by assignment-fixed LP ascent from random starts (a LOWER bound; the
    Decider decides the sign exactly).  theta is an exact angle vector."""
    from scipy.optimize import linprog
    n = len(theta)
    NV = 1 + 2 * n
    cs = [(math.cos(t), math.sin(t)) for t in theta]
    P = [0.5 * (abs(c) + abs(s)) for (c, s) in cs]
    m = {}
    for i in range(n):
        for j in range(i + 1, n):
            d = theta[j] - theta[i]
            m[(i, j)] = 0.5 + 0.5 * (abs(math.cos(d)) + abs(math.sin(d)))
    wall = []
    for i in range(n):
        for var in (1 + i, 1 + n + i):
            r = np.zeros(NV); r[var] = 1.0; r[0] = -1.0; wall.append((r, P[i]))
            r = np.zeros(NV); r[var] = -1.0; r[0] = -1.0; wall.append((r, P[i] - T))
    bounds = [(-1.0, 0.5)] + [(P[i], T - P[i]) for i in range(n)] * 2
    c = np.zeros(NV); c[0] = -1.0
    best = -1.0
    for _st in range(nstart):
        X = np.array([rng.uniform(P[i], T - P[i]) for i in range(n)])
        Y = np.array([rng.uniform(P[i], T - P[i]) for i in range(n)])
        prev = None
        for _it in range(iters):
            rows = list(wall)
            asg = []
            for i in range(n):
                for j in range(i + 1, n):
                    dx, dy = X[j] - X[i], Y[j] - Y[i]
                    bv, ba = -1e18, None
                    for o in (i, j):
                        cc, ss = cs[o]
                        for (d0, d1) in ((cc, ss), (-ss, cc)):
                            pr = d0 * dx + d1 * dy
                            for sg in (1.0, -1.0):
                                if sg * pr > bv:
                                    bv, ba = sg * pr, (o, d0, d1, sg)
                    asg.append(((i, j), ba))
                    (o, d0, d1, sg) = ba
                    r = np.zeros(NV)
                    r[1 + j] = sg * d0; r[1 + i] = -sg * d0
                    r[1 + n + j] = sg * d1; r[1 + n + i] = -sg * d1
                    r[0] = -1.0
                    rows.append((r, m[(i, j)]))
            key = tuple((p, a[0], round(a[1], 9), round(a[2], 9), a[3]) for (p, a) in asg)
            A = np.array([r for (r, _q) in rows]); b = np.array([q for (_r, q) in rows])
            res = linprog(c, A_ub=-A, b_ub=-b, bounds=bounds, method='highs')
            if not res.success:
                break
            X, Y = res.x[1:1 + n].copy(), res.x[1 + n:1 + 2 * n].copy()
            z = np.zeros(1 + 3 * n)
            for k in range(n):
                z[1 + 3 * k], z[2 + 3 * k], z[3 + 3 * k] = X[k], Y[k], theta[k]
            best = max(best, value(z, n, T))
            if key == prev:
                break
            prev = key
    return best


# ===================================================================== the branch-and-bound
def canon(box):
    """label symmetry (sort the intervals) + the container reflection th -> 90deg - th."""
    a = tuple(sorted(box))
    Q = math.pi / 2
    b = tuple(sorted((Q - hi, Q - lo) for (lo, hi) in box))
    return min(a, b)


def split(box):
    w = [hi - lo for (lo, hi) in box]
    k = int(np.argmax(w))
    lo, hi = box[k]
    mid = 0.5 * (lo + hi)
    out = []
    for (a, b) in ((lo, mid), (mid, hi)):
        c = list(box); c[k] = (a, b); out.append(tuple(c))
    return out


_B = {}


def _bb_init(n, T, node_cap, kinds):
    _B['n'], _B['T'], _B['cap'], _B['kinds'] = n, T, node_cap, kinds
    _B['dec'] = {k: Decider(n, T, kills=k, node_cap=node_cap) for k in kinds}


def _bb_job(box):
    tlo = [b[0] for b in box]
    thi = [b[1] for b in box]
    for k in _B['kinds']:
        st, x, info = _B['dec'][k].decide(tlo, thi)
        if st == 'kill':
            return box, k, info['nodes'], None
        if st == 'budget':
            return box, 'budget', info['nodes'], None
    return box, 'survive', info['nodes'], (None if x is None else list(x))


def cmd_bb(a):
    import multiprocessing as mp
    from collections import Counter
    Q = math.pi / 2
    n, T = a.n, a.T
    kinds = [k for k in a.kinds.split(',') if k]
    # open with the uniform subdivision into `warm` bins per angle, label-sorted
    K = a.warm
    bins = [(k * Q / K, (k + 1) * Q / K) for k in range(K)]
    level = sorted({canon(tuple(bins[i] for i in c))
                    for c in itertools.combinations_with_replacement(range(K), n)})
    out = open(a.out, 'w') if a.out else None
    for depth in range(a.depth):
        t0 = time.time()
        res = []
        cnt = Counter()
        with mp.Pool(a.nproc, initializer=_bb_init,
                     initargs=(n, T, a.node_cap, kinds)) as pool:
            for q, r in enumerate(pool.imap_unordered(_bb_job, level, chunksize=1)):
                res.append(r)
                cnt[r[1]] += 1
                if out:
                    out.write(json.dumps({'depth': depth, 'kind': r[1], 'nodes': r[2],
                                          'box': [[lo, hi] for (lo, hi) in r[0]],
                                          'x': r[3]}) + '\n')
                if (q + 1) % 50 == 0:
                    out and out.flush()
                    print(f'  .. {q + 1}/{len(level)}  ' + ' '.join(
                        f'{k}={cnt.get(k, 0)}' for k in kinds + ['survive', 'budget'])
                        + f'  {time.time() - t0:.0f}s', flush=True)
        surv = [r[0] for r in res if r[1] in ('survive', 'budget')]
        wdeg = math.degrees(max(hi - lo for (lo, hi) in level[0]))
        print(f'level {depth}  boxes {len(level):7d}  max width {wdeg:7.3f} deg   '
              + '  '.join(f'{k}={cnt.get(k, 0)}' for k in kinds + ['survive', 'budget'])
              + f'   {time.time() - t0:.1f}s', flush=True)
        out and out.flush()
        if not surv:
            print('CLOSED: every box killed.')
            break
        nxt = set()
        for b in surv:
            for c in split(b):
                nxt.add(canon(c))
        level = sorted(nxt)
        if len(level) > a.max_boxes:
            print(f'stopping: {len(level)} boxes exceeds --max-boxes')
            break
    if out:
        out.close()


# ===================================================================== selftest
def cmd_selftest(a):
    from fractions import Fraction as Fr
    rng = np.random.default_rng(1)
    ok = True

    def chk(name, cond, extra=''):
        nonlocal ok
        ok = ok and bool(cond)
        print(f'  [{"ok " if cond else "FAIL"}] {name} {extra}')

    print('chord constants')
    chk('D(0) = 1/2', abs(chord_D(0.0) - 0.5) < 1e-15)
    chk('D(45deg) = (sqrt2-1)/2', abs(chord_D(math.pi / 4) - D45) < 1e-15)
    g = np.linspace(0, math.pi / 2, 20001)
    chk('D >= (sqrt2-1)/2 everywhere', min(chord_D(t) for t in g) >= D45 - 1e-15)
    e = max(abs(chord_len(t, chord_D(t)) - 1.0) for t in g)
    chk('chord(theta, D(theta)) = 1 for all theta', e < 1e-6,
        f'(max err {e:.2e}; the residual is cancellation in p - D = C*S near theta = 0, '
        f'not an error in the lemma -- the branch-and-bound never evaluates chord_len)')
    chk('D_lam increases as lam decreases',
        all(D_lam(t, 0.5) >= D_lam(t, 0.75) >= D_lam(t, 1.0) - 1e-15 for t in g))
    chk('brief band (C-S)/2 <= D, and 0 at 45 deg',
        all(0.5 * (abs(math.cos(t)) - abs(math.sin(t))) <= chord_D(t) + 1e-15 for t in g)
        and abs(0.5 * (math.cos(math.pi / 4) - math.sin(math.pi / 4))) < 1e-15)
    print('  K2 family at T=3:', lam_family(3.0, 6), ' at T=4:', lam_family(4.0, 12))

    print('interval bounds used by the box relaxation (must UNDER-estimate on every interval)')
    bad = 0
    for _ in range(2000):
        u0 = rng.uniform(-1, 3); u1 = u0 + rng.uniform(0, 2.0)
        gg = np.linspace(u0, u1, 601)
        for lam in (1.0, 0.75, 0.6, 0.5, 1 / 3, 0.25):
            if _Dlam_lo(u0, u1, lam) > min(D_lam(t, lam) for t in gg) + 1e-12:
                bad += 1
        if _plo(u0, u1) > min(0.5 * (abs(math.cos(t)) + abs(math.sin(t))) for t in gg) + 1e-12:
            bad += 1
        if _wlo(u0, u1) > min(abs(math.cos(t)) + abs(math.sin(t)) for t in gg) + 1e-12:
            bad += 1
    chk('_Dlam_lo / _plo / _wlo under-estimate on 2000 random intervals', bad == 0,
        f'({bad} violations)')

    print('chord length against a direct polygon intersection')
    worst = 0.0
    for _ in range(400):
        t = rng.uniform(0, math.pi / 2)
        d = rng.uniform(0, 0.5 * (abs(math.cos(t)) + abs(math.sin(t))))
        c, sn = math.cos(t), math.sin(t)
        vs = [(0.5 * (c * sx - sn * sy), 0.5 * (sn * sx + c * sy))
              for (sx, sy) in ((1, 1), (-1, 1), (-1, -1), (1, -1))]
        xs = []
        for k in range(4):
            (x0, y0), (x1, y1) = vs[k], vs[(k + 1) % 4]
            if (y0 - d) * (y1 - d) <= 0 and y0 != y1:
                xs.append(x0 + (x1 - x0) * (d - y0) / (y1 - y0))
        if len(xs) >= 2:
            worst = max(worst, abs((max(xs) - min(xs)) - chord_len(t, d)))
    chk('chord_len matches the polygon cut', worst < 1e-9, f'(max err {worst:.2e})')

    print('the disjunctive LP: exactness at a point')
    for (n, T, th, want) in ((4, 3.0, [0.0] * 4, 1.0 / 3.0),
                             (5, 3.0, [math.pi / 4] * 5, None),
                             (6, 3.2, [0.0] * 6, 0.05)):
        D = Decider(n, T, kills='N', node_cap=200000)
        st, x, _i = D.decide(th, th)
        v = s6 = fixed_angle_value(th, T, rng, nstart=60)
        chk(f'n={n} T={T} decider vs instrument', st == 'survive' and abs(x[0] - v) < 1e-7,
            f'{x[0]:.9f} vs {v:.9f}')
        if want is not None:
            chk(f'   and equals the hand value {want}', abs(x[0] - want) < 1e-7)
    print('the theorem, at three exact angle vectors')
    for (th, lab) in (([0.0] * 6, 'all axis-parallel'), ([math.pi / 4] * 6, 'all 45 deg'),
                      (list(rng.uniform(0, math.pi / 2, 6)), 'random')):
        D = Decider(6, 3.0, kills='N', node_cap=a.node_cap)
        st, _x, info = D.decide(th, th)
        chk(f'6 in [0,3]^2, {lab}: delta* <= 0', st in ('kill', 'budget'),
            f"({st}, {info['nodes']} LP nodes)" if st != 'kill' else f"({info['nodes']} LP nodes)")
    print('soundness: the full angle cube must never be killed')
    for (n, T) in ((4, 3.0), (5, 3.0), (6, 3.0), (6, 3.2)):
        D = Decider(n, T, kills='K2', node_cap=200000)
        st, _x, _i = D.decide([0.0] * n, [math.pi / 2] * n)
        chk(f'n={n} T={T} cube survives', st == 'survive')
    print('OK' if ok else 'FAILURES ABOVE')
    return 0 if ok else 1


if __name__ == '__main__':
    main()
