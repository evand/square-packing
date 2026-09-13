#!/usr/bin/env python3
"""bentz.py -- Bentz's case analysis as incidence-pattern branching on point regions.

Task: `tasks/bentz-incidence/README.md`.  Note: `search/BENTZ.md`.

The branching object is the **pattern** of a pose `S` with respect to a finite point set
`P0 subset [0,4]^2`:

    pi(S) = { p in P0 : p in S }        (CLOSED containment; a point on dS counts)

Everything here that decides a pattern is done in exact integer arithmetic: a point is
`(X, Y, D)` with integer `X, Y, D`, a pose is `(p, q, cx, cy)` with integer `p, q` and rational
`cx, cy`, and membership is `leaf_ceiling.sq_contains`, four integer half-plane tests in the
square's own frame.  No floats decide a pattern anywhere.

Facts used (`tasks/bentz-incidence/README.md` "Semantics"):

  * a packing is a family of closed unit squares, pairwise disjoint as CLOSED sets, so two
    packing squares never contain a common point;
  * hence for `pi != empty`, `R_pi = {S : pi(S) = pi}` is contained in the point region
    `R_p = {S : p in S}` for any `p in pi`, and `R_p` is a clique: `mu(R_pi) <= 1`
    (`card_filter_clique_le_one` in Lean, `notes/level2-design.md` §2.2);
  * two patterns used by the same packing are disjoint subsets of `P0`.

Sub-commands
------------
  points   : print `P0` (exact) and its D4 orbit structure
  cover    : is `P0` a closed cover of `[0,4]^2`?  (emit a unit-weight cover file for
             `search/zeromargin.py` / `verify2/zmcheck`, and run a lattice pre-screen here)
  patterns : the pattern census of a pose set / measure file (exact), with D4 classes
  leaves   : enumerate the corner-pattern level up to D4 and, for each leaf, the set-packings
             of the remaining points that Bentz's counting `k + u <= 4` allows
  spec     : write the JSON leaf spec that `cliquelever.py --pose-filter` consumes
  seed     : write a pose file (cliquelever `--load-poses` format) of poses realising each
             admitted pattern of a leaf, including the nudged grid squares of
             `search/family_rows.py`'s item-3 family
  filter   : (library entry point, also a CLI self-check) build the admit predicate

The LP itself is `search/cliquelever.py` in the `PGCORN` configuration of
`runs/launch_2026-09-12_pgonly.sh`, with `--pose-filter SPEC.json` (a ~20-line hook added by this
task) restricting BOTH the loaded pose set and every priced column to the leaf's patterns, and
`--count-rows` adding the leaf's `mu(R_pi) = k_pi` equalities.
"""
import argparse
import json
import math
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_ceiling as lc                                             # noqa: E402

T = Fr(4)

# ===================================================================== P0: Bentz's 16 points
# B10 §3 / notes/proof-anatomy.md §2.3, at m = 4:
#   A(1, 0.914), B(0.914, 1), C(0.914, 2), D(1.65, 1.65) and their D4 images.
# 0.914 = 457/500 (Bentz's own decimal; the exact constant behind it is sqrt(2) - 1/2 = 0.914213...,
# and he rounds DOWN so that the lemma hypotheses hold -- we keep his decimal, exactly).
BENTZ_A = (Fr(1), Fr(457, 500))
BENTZ_C = (Fr(457, 500), Fr(2))
BENTZ_D = (Fr(33, 20), Fr(33, 20))


def d4_point(x, y, t=T):
    """the D4 orbit of a point of [0,t]^2 (as a sorted tuple of distinct points)"""
    o = {(x, y), (t - x, y), (x, t - y), (t - x, t - y),
         (y, x), (t - y, x), (y, t - x), (t - y, t - x)}
    return tuple(sorted(o))


def corner_of(x, y, t=T):
    """which corner box [0,1]^2 / [3,4]x[0,1] / [0,1]x[3,4] / [3,4]^2 a corner point belongs to.
    Matches `leaf_ceiling.region_boxes`' C0..C3 numbering."""
    i = 0
    if x > t / 2:
        i |= 1
    if y > t / 2:
        i |= 2
    return i


def bentz_points():
    """-> [(name, x, y)] in a fixed, documented order.

    a0,b0,a1,b1,a2,b2,a3,b3 : the eight corner points, a_i/b_i the pair of corner i
                              (a_i is the one on a vertical grid line x in {1,3}, b_i on y in {1,3})
    c0..c3                  : the four mid-wall points (D4 orbit of C(0.914, 2))
    d0..d3                  : the four interior points (D4 orbit of D(1.65, 1.65))
    """
    out = []
    corner = {}
    for (x, y) in d4_point(*BENTZ_A):
        i = corner_of(x, y)
        corner.setdefault(i, []).append((x, y))
    for i in range(4):
        pts = sorted(corner[i])
        # a_i: the point whose x is an integer (on a vertical grid line); b_i: the other
        ai = [p for p in pts if p[0].denominator == 1]
        bi = [p for p in pts if p[0].denominator != 1]
        assert len(ai) == 1 and len(bi) == 1, (i, pts)
        out.append((f'a{i}', ai[0][0], ai[0][1]))
        out.append((f'b{i}', bi[0][0], bi[0][1]))
    for j, (x, y) in enumerate(d4_point(*BENTZ_C)):
        out.append((f'c{j}', x, y))
    for j, (x, y) in enumerate(d4_point(*BENTZ_D)):
        out.append((f'd{j}', x, y))
    return out


def as_int_points(pts):
    """[(name, x, y)] -> [(name, X, Y, D)] with a COMMON integer denominator D"""
    D = 1
    for (_, x, y) in pts:
        D = D * x.denominator // math.gcd(D, x.denominator)
        D = D * y.denominator // math.gcd(D, y.denominator)
    return [(n, int(x * D), int(y * D), D) for (n, x, y) in pts]


# ===================================================================== exact pattern of a pose
def pattern_of(sq, ipts):
    """exact pattern of the square record `sq` (leaf_ceiling.make_square) as a frozenset of the
    indices of `ipts = [(name, X, Y, D)]` it contains.  Integer arithmetic only."""
    return frozenset(k for k, (_, X, Y, D) in enumerate(ipts) if lc.sq_contains(sq, X, Y, D))


def pattern_of_pose(p, q, cx, cy, ipts, t=T):
    return pattern_of(lc.make_square(Fr(cx), Fr(cy), p, q, t), ipts)


def pname(pat, ipts):
    return '{' + ','.join(ipts[k][0] for k in sorted(pat)) + '}' if pat else '{}'


# ===================================================================== D4 action on P0
def d4_maps(t=T):
    """the eight elements of D4 on [0,t]^2 as functions (x, y) -> (x', y'), with names"""
    return [('e', lambda x, y: (x, y)),
            ('mx', lambda x, y: (t - x, y)),
            ('my', lambda x, y: (x, t - y)),
            ('r2', lambda x, y: (t - x, t - y)),
            ('d', lambda x, y: (y, x)),
            ('r1', lambda x, y: (t - y, x)),
            ('r3', lambda x, y: (y, t - x)),
            ('ad', lambda x, y: (t - y, t - x))]


def d4_perm(pts, t=T):
    """-> [perm] : for each g in D4, the permutation of the indices of `pts` it induces"""
    idx = {(x, y): k for k, (_, x, y) in enumerate(pts)}
    out = []
    for (_, g) in d4_maps(t):
        out.append(tuple(idx[g(x, y)] for (_, x, y) in pts))
    return out


def canon(obj, perms):
    """canonical form of a frozenset-of-frozensets (or a tuple) under the D4 permutations"""
    best = None
    for pm in perms:
        img = frozenset(frozenset(pm[k] for k in s) for s in obj)
        key = tuple(sorted(tuple(sorted(s)) for s in img))
        if best is None or key < best:
            best = key
    return best


# ===================================================================== pose filters (the hook)
class Filter:
    """the predicate `cliquelever.Poses.add` consults.  Built from a JSON spec so that a run is
    reproducible from its command line alone."""

    def __init__(self, spec):
        self.spec = spec
        self.kind = spec['kind']
        self.t = Fr(spec.get('t', '4'))
        self.n_admit = 0
        self.n_reject = 0
        if self.kind == 'axis':
            # min(|cos|, |sin|) <= num/den, exactly: min(|a|,|b|) * den <= num * r,
            # with (a, b, r) = (q^2 - p^2, 2pq, q^2 + p^2).  num = 0 <=> exactly axis-parallel.
            self.num = int(spec['sin_num'])
            self.den = int(spec['sin_den'])
        elif self.kind == 'pattern':
            self.ipts = [(n, int(X), int(Y), int(D)) for (n, X, Y, D) in spec['points']]
            self.allow = set(frozenset(s) for s in spec['patterns'])
        else:
            raise SystemExit(f'bentz.Filter: unknown kind {self.kind!r}')

    def __call__(self, p, q, cx, cy):
        ok = self._test(p, q, cx, cy)
        if ok:
            self.n_admit += 1
        else:
            self.n_reject += 1
        return ok

    def _test(self, p, q, cx, cy):
        if self.kind == 'axis':
            a, b, r = lc.cos_sin(p, q)
            return min(abs(a), abs(b)) * self.den <= self.num * r
        sq = lc.make_square(Fr(cx), Fr(cy), p, q, self.t)
        pat = pattern_of(sq, self.ipts)
        self._last = ((p, q, cx, cy), pat)
        return pat in self.allow

    def pattern_index(self, p, q, cx, cy):
        """for the count rows: the index of this pose's pattern in `spec['patterns']`, or -1.
        Reuses the pattern just computed by `_test` for the same pose (the caller is
        `cliquelever.Poses.add`, which tests admission immediately before)."""
        if self.kind != 'pattern':
            return -1
        key = (p, q, cx, cy)
        if getattr(self, '_last', (None, None))[0] == key:
            pat = self._last[1]
        else:
            pat = pattern_of(lc.make_square(Fr(cx), Fr(cy), p, q, self.t), self.ipts)
        return self.order.get(pat, -1)

    @property
    def order(self):
        if not hasattr(self, '_order'):
            self._order = {frozenset(s): i for i, s in enumerate(self.spec.get('patterns', []))}
        return self._order

    @property
    def counts(self):
        """[(pattern index, rhs)] for the equality rows, or [] if the spec pins none"""
        return list(enumerate(self.spec.get('counts', []))) if self.spec.get('counts') else []

    def describe(self):
        if self.kind == 'axis':
            if self.num == 0:
                return 'axis-parallel exactly (theta = 0 mod 90 deg)'
            th = math.degrees(math.asin(self.num / self.den))
            return (f'within {th:.7f} deg of axis-parallel '
                    f'(min(|cos|,|sin|) <= {self.num}/{self.den}, exact)')
        return (f'{len(self.allow)} admitted patterns over {len(self.ipts)} points'
                + (f'; {len(self.spec["counts"])} count rows' if self.spec.get('counts') else ''))


def load_filter(path_or_json):
    spec = json.loads(open(path_or_json).read()) if os.path.exists(path_or_json) \
        else json.loads(path_or_json)
    return Filter(spec)


# ===================================================================== CLI
def cmd_points(a):
    pts = bentz_points()
    ip = as_int_points(pts)
    print(f'P0 = Bentz\'s 16 points at m = 4 (notes/proof-anatomy.md 2.3), common denominator '
          f'{ip[0][3]}')
    for (n, x, y), (_, X, Y, D) in zip(pts, ip):
        print(f'  {n:>3}  ({x}, {y}) = ({float(x):.6f}, {float(y):.6f})   = ({X}/{D}, {Y}/{D})')
    perms = d4_perm(pts)
    print(f'D4 orbits: a/b = {{{",".join(n for n,_,_ in pts[:8])}}} (one orbit of 8), '
          f'c = 4, d = 4')
    for (nm, _), pm in zip(d4_maps(), perms):
        print(f'  {nm:>3}: ' + ' '.join(f'{pts[k][0]}' for k in pm))


def cmd_patterns(a):
    """exact pattern census of a cliquelever measure / pose file (`pose p q cx cy mass`)"""
    pts = bentz_points()
    ipts = as_int_points(pts)
    t, sym, poses = lc.read_measure(a.FILE)
    assert sym == 1, 'D4-symmetrised files are not handled here'
    cen = {}
    tot = Fr(0)
    for (p, q, cx, cy, m) in poses:
        if a.min_mass and float(m) < a.min_mass:
            continue
        pat = pattern_of(lc.make_square(cx, cy, p, q, t), ipts)
        e = cen.setdefault(pat, [Fr(0), 0, []])
        e[0] += m
        e[1] += 1
        if len(e[2]) < 4:
            e[2].append((p, q, cx, cy, m))
        tot += m
    print(f'# {a.FILE}: {len(poses)} poses, mass {float(tot):.9f}, {len(cen)} distinct patterns')
    print(f'# pattern              mass        poses   example pose (p q cx cy) mass')
    for pat, (m, n, ex) in sorted(cen.items(), key=lambda kv: (-float(kv[1][0]), sorted(kv[0]))):
        p, q, cx, cy, mm = ex[0]
        th = math.degrees(2 * math.atan2(p, q))
        print(f'{pname(pat, ipts):>22} {float(m):11.6f} {n:7d}   '
              f'({float(cx):.6f}, {float(cy):.6f}, {th:+.3f} deg) m={float(mm):.6f}')
    return cen, ipts


def cmd_seed(a):
    """write a cliquelever `--load-poses` file (`pose cx cy theta_deg mu`).

    Two families, both of which the `0.04 / 2.5 deg` pricing lattice steps over
    (`search/HONEST.md` §0 item 3, `search/family_rows.py`):

      * **item 3** -- `c_x` (or `c_y`) on a half-integer line, so BOTH edges of that axis sit on
        grid lines, the other coordinate on a fine pitch, at a geometric ladder of small angles;
      * **nudged grid** -- the 16 squares of the `4 x 4` tiling with each edge pushed `eps` off
        its grid line, for a ladder of `eps` down to the `1e-6` the pose denominator can hold.
        These are the poses whose PATTERN differs from the on-grid one: `[1,2]x[0,1]` contains
        `a1 = (3, 0.914)`? no -- but `[1+eps, 2+eps] x [0,1]` contains `c1 = (2, 0.914)` and the
        un-nudged one contains it too, while `[2,3]x[0,1]` contains BOTH `c1` and `a1` and
        `[2+eps,3+eps]x[0,1]` contains only `a1`.  Pattern branching lives or dies on these.
    """
    import numpy as np
    s = 4.0
    angles = [0.0] + [10 ** e for e in (-6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5)]
    angles = angles + [-t for t in angles[1:]]
    out = []
    lines = [k + 0.5 for k in range(int(s))]
    for th in angles:
        w = (abs(math.cos(th)) + abs(math.sin(th))) / 2
        lo, hi = w + 1e-7, s - w - 1e-7
        for c in lines:
            cc = min(max(c, lo), hi)
            for v in np.arange(lo, hi + 1e-12, a.pitch):
                out.append((cc, float(v), th))
                out.append((float(v), cc, th))
    for e in (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        for sx in (1, -1):
            for sy in (1, -1):
                for i in range(4):
                    for j in range(4):
                        cx = min(max(i + 0.5 + sx * e, 0.5), s - 0.5)
                        cy = min(max(j + 0.5 + sy * e, 0.5), s - 0.5)
                        out.append((cx, cy, 0.0))
    filt = load_filter(a.spec) if a.spec else None
    kept = 0
    with open(a.out, 'w') as f:
        f.write(f'# bentz.py seed; t=4; item-3 family (pitch {a.pitch}) + nudged grid; '
                f'spec={a.spec}\n')
        for (cx, cy, th) in out:
            if filt is not None:
                p, q, x, y = lc.snap_pose(cx, cy, th, T, a.Q, a.Dc)
                if not filt(p, q, x, y):
                    continue
            f.write(f'pose {cx:.17g} {cy:.17g} {math.degrees(th):.17g} 0\n')
            kept += 1
    print(f'{kept} of {len(out)} poses written to {a.out}'
          + (f' (filter: {filt.describe()})' if filt else ''))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('points')
    c = sub.add_parser('patterns')
    c.add_argument('FILE')
    c.add_argument('--min-mass', type=float, default=0.0)
    c = sub.add_parser('seed')
    c.add_argument('out')
    c.add_argument('--pitch', type=float, default=0.02)
    c.add_argument('--spec', default=None)
    c.add_argument('--Q', type=int, default=10 ** 5)
    c.add_argument('--Dc', type=int, default=10 ** 6)
    a = ap.parse_args()
    if a.cmd == 'points':
        cmd_points(a)
    elif a.cmd == 'patterns':
        cmd_patterns(a)
    elif a.cmd == 'seed':
        cmd_seed(a)


if __name__ == '__main__':
    main()
