#!/usr/bin/env python3
"""rung2_bound.py -- the m^2 obstruction of search/RUNG2.md, and its check on a certificate.

A *monotone witness certificate* for a pose box B is a set S of cover points with total weight
>= 1 such that EVERY p in S lies in the closed unit square Q(c, theta) at EVERY admissible pose
(c, theta) of B.  CORE, P1 and the ADM primitive of search/zeromargin.py all produce exactly this
kind of certificate, and so does any union of them; the triangle lemma TRI does NOT (it is
disjunctive: it says *some* vertex is captured, not that a fixed one is).

THEOREM (search/RUNG2.md sec 2).  Fix a tile (i, j) of [0,m]^2 and signs sx, sy with sx = +1 if
i < m-1 (the centre may move right) and sx = -1 if i = m-1 (it may only move left), likewise sy.
The pose (i+1/2, j+1/2, 0) is admissible and, the leaves being closed, finitely many and space-
filling, one leaf box B contains it together with poses (i+1/2 + sx e, j+1/2 + sy e, theta) for
arbitrarily small e, theta > 0.  B's monotone witness set S therefore lies in

    PIN(i, j, sx, sy) = intersection of the closed unit squares over that whole family,

so w(PIN(i, j, sx, sy)) >= 1 for every tile and every admissible sign pair.  Minimising the total
weight subject to those constraints is an LP; on a 1/20 or 1/40 lattice of candidate points its
value is exactly m^2 = 16, and the LP value on a lattice is an upper bound for the continuum one,
so the computed 16 is the bound for lattice covers (every cover this project builds is one: the
columns of search/closed4.py live on a 1/1000 lattice).

So no cover of [0,4]^2 with total weight below 16 -- in particular none with W < 13 -- can be
certified by CORE / P1 / ADM and their unions alone, at any subdivision depth.  Rung 2 needs a
DISJUNCTIVE primitive (the weighted analogue of TRI / the two-region argument).

Modes
  bound  [--m 4] [--pitch 20] [--octants]   the LP behind the theorem, on a lattice: with the
         one-sided (octant) constraints the value is exactly m^2; without them (i.e. if only the
         *two-sided* information at each tile pose is used) it drops to m^2/2, which is why the
         one-sided reading matters.
  check  CERT                               the direct test on a certificate-format cover: the
         exact (Fraction) weight of PIN(i,j,sx,sy) for every tile and sign pair, and the list of
         those below 1 -- each of which is a pose the box checker can never certify.

usage: python3 search/rung2_bound.py bound [--m 4] [--pitch 20] [--octants]
       python3 search/rung2_bound.py check runs/closed4_best_x103.txt
"""
import argparse, math, sys
from fractions import Fraction as F
import numpy as np


def pin_mask(X, Y, m, i, j, thetas, sx=0, sy=0, tol=1e-12):
    """Points captured at every pose of the admissible curve through the tile-(i,j) pose
    (i+1/2, j+1/2, 0).  sx/sy != 0 additionally move the centre in that direction, i.e. restrict
    to the poses of a leaf box that touches the tile pose from one side."""
    ok = np.ones(len(X), dtype=bool)
    for th in thetas:
        w = math.cos(th) + math.sin(th)
        for dcx in ((0.0,) if sx == 0 else (0.0, sx * 1e-3)):
            for dcy in ((0.0,) if sy == 0 else (0.0, sy * 1e-3)):
                cx = min(max(i + 0.5 + dcx, w / 2), m - w / 2)
                cy = min(max(j + 0.5 + dcy, w / 2), m - w / 2)
                c, s = math.cos(th), math.sin(th)
                dx, dy = X - cx, Y - cy
                ok &= (np.abs(dx * c + dy * s) <= 0.5 + tol) & (np.abs(-dx * s + dy * c) <= 0.5 + tol)
    return ok


def bound(a):
    from scipy.optimize import linprog
    m, K = a.m, a.pitch
    g = np.arange(0, m * K + 1) / K
    GX, GY = np.meshgrid(g, g, indexing='ij')
    X, Y = GX.ravel(), GY.ravel()
    thetas = [0.0, 1e-4, 3e-4, 1e-3, 3e-3]
    rows, labels = [], []
    for i in range(m):
        for j in range(m):
            sxs = [0] if not a.octants else ([1] if i == 0 else [-1] if i == m - 1 else [1, -1])
            sys_ = [0] if not a.octants else ([1] if j == 0 else [-1] if j == m - 1 else [1, -1])
            for sx in sxs:
                for sy in sys_:
                    msk = pin_mask(X, Y, m, i, j, thetas, sx, sy)
                    rows.append(msk.astype(float)); labels.append((i, j, sx, sy, int(msk.sum())))
    A = np.array(rows)
    res = linprog(c=np.ones(len(X)), A_ub=-A, b_ub=-np.ones(len(rows)), bounds=(0, None), method='highs')
    print(f"container [0,{m}]^2, lattice pitch 1/{K} ({len(X)} candidate points), "
          f"{len(rows)} PIN constraints, octants={a.octants}")
    print(f"LP value (lower bound on any monotone-witness-certifiable cover on this lattice) = {res.fun:.6f}")
    w = res.x; sup = np.nonzero(w > 1e-9)[0]
    print(f"support: {len(sup)} points, total {w.sum():.6f}")
    for k in sup[np.argsort(-w[sup])][:20]:
        print(f"   ({X[k]:.3f}, {Y[k]:.3f}) w={w[k]:.6f}")


def check(path):
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    ms = F(sn, sd)
    assert ms.denominator == 1, "container side must be an integer for the tile test"
    m = int(ms)
    Xi = np.empty(n, dtype=np.int64); Yi = np.empty(n, dtype=np.int64); Wi = [0] * n
    for k in range(n):
        x, y, w = map(int, tok[5 + 3 * k: 8 + 3 * k])
        Xi[k] = x; Yi[k] = y; Wi[k] = w
    tot = F(sum(Wi), W)
    Xf = Xi / D; Yf = Yi / D
    print(f"container [0,{m}]^2, {n} points, total weight {tot} = {float(tot):.9f}")
    print("w(PIN(i,j,sx,sy)) for every tile and every admissible sign pair -- each must be >= 1")
    print("for a monotone witness certificate to exist at the pose (i+1/2, j+1/2, 0):")
    thetas = [0.0, 1e-5, 1e-4, 1e-3]
    bad = []
    for j in range(m - 1, -1, -1):
        cells = []
        for i in range(m):
            sxs = [1] if i == 0 else [-1] if i == m - 1 else [1, -1]
            sys_ = [1] if j == 0 else [-1] if j == m - 1 else [1, -1]
            worst = None
            for sx in sxs:
                for sy in sys_:
                    msk = pin_mask(Xf, Yf, m, i, j, thetas, sx, sy)
                    v = F(sum(Wi[k] for k in np.nonzero(msk)[0]), W)
                    if worst is None or v < worst[0]: worst = (v, sx, sy)
                    if v < 1: bad.append((i, j, sx, sy, v))
            cells.append(worst[0])
        print("  j=%d: " % j + "  ".join(f"{float(v):8.5f}" for v in cells))
    print(f"\nPIN sets below 1: {len(bad)}")
    for i, j, sx, sy, v in bad:
        print(f"   tile ({i},{j}) sx={sx:+d} sy={sy:+d}: weight {v} = {float(v):.6f}  -> the pose "
              f"({i + 0.5}, {j + 0.5}, 0) can never be certified by a monotone witness set")
    print("VERDICT:", "no monotone-witness obstruction found at the tile poses" if not bad else
          "OBSTRUCTED -- this cover cannot be certified by CORE/P1/ADM at any depth")


def dual(a):
    """The explicit dual certificate of Theorem 1, verified on a complete set of cell
    representatives.  PIN-membership depends only on the position of a point relative to the
    half-integer grid, so it is constant on every cell of that arrangement, and the 1/4-lattice
    contains a representative of every cell (2-cells at (k+1/2)/2, 1-cells at their midpoints,
    0-cells themselves).  Checking `multiplicity <= 1` on the 1/4-lattice therefore checks it on
    ALL of [0,m]^2 -- the bound that follows holds for every cover, not only lattice ones."""
    m = a.m
    sx = [1, 1] + [-1] * (m - 2)       # the sign pattern of Theorem 1: the +/- switch must sit
                                       # between two INTERIOR columns (Lemma 2, case (iv))
    sy = list(sx)
    if a.switch is not None:
        sx = [1] * (a.switch + 1) + [-1] * (m - a.switch - 1); sy = list(sx)
    g = np.arange(0, m * 4 + 1) / 4
    GX, GY = np.meshgrid(g, g, indexing='ij')
    X, Y = GX.ravel(), GY.ravel()
    thetas = [0.0, 1e-5, 1e-4, 1e-3]
    mult = np.zeros(len(X), dtype=int)
    for i in range(m):
        for j in range(m):
            mult += pin_mask(X, Y, m, i, j, thetas, sx[i], sy[j]).astype(int)
    print(f"container [0,{m}]^2; dual certificate lambda = 1 on the {m*m} constraints")
    print(f"  sigma_x per column: {sx}")
    print(f"  sigma_y per row:    {sy}")
    print(f"sum of lambda = {m*m}")
    print(f"max multiplicity over the 1/4-lattice ({len(X)} cell representatives) = {mult.max()}")
    bad = np.nonzero(mult > 1)[0]
    for k in bad[:10]: print(f"   OVERLAP at ({X[k]}, {Y[k]}): multiplicity {mult[k]}")
    print("VERDICT:", f"the {m*m} PIN sets are pairwise disjoint  =>  every monotone-witness-"
          f"certifiable cover of [0,{m}]^2 has total weight >= {m*m}" if mult.max() <= 1 else
          "NOT disjoint for this sign pattern")



# ------------------------------------------------------------------ T2: credited certificates
# A *monotone-credited* certificate (certificates/FORMAT.md "Clique certificates", "Anchor
# cliques", "Branch certificates"; notes/clique-family.md sec 5) carries points p with weights
# w_p, cliques K with weights w_K, and optionally corner-region multipliers lambda_j.  It asserts
# that every admissible pose captures  sum_{p in S} w_p + sum_{K ni S} w_K >= 1  (>= 1 + lambda_j
# if the centre is in corner box j), and its accounting total is
#     Theta = sum_p w_p + sum_K w_K - sum_j lambda_j K_j.
# *Monotone* crediting is the rule the sweep verifier and xcheck.py implement: a cell is credited
# a clique (or the region threshold) only when EVERY pose of the cell is a member -- "a cell that
# only partially satisfies a predicate gets no credit" (clique-family.md sec 5).
#
# Lemma 0 of sec 2 applies verbatim: the leaf that carries the germ of the tile pose can only be
# credited a clique K with germ ⊆ K.  A clique is a pairwise-intersecting family of poses, so a
# clique may serve two tile germs only if every square of one germ meets every square of the
# other.  Lemma T (search/RUNG2.md sec 9) settles which pairs those are:
#   * same tile, any signs               -- always compatible;
#   * tiles (i,j), (i+1,j)               -- compatible iff sigma_x = +1 on the left tile and
#                                           -1 on the right one (and likewise in y);
#   * diagonal or farther apart          -- never.
# Since the compatibility graph on tiles has no triangle, a clique serves at most two tiles, and
# a clique column of the LP below is a maximal pairwise-compatible set of (tile, sign) constraints.

def _compatible(t1, s1, t2, s2):
    """Lemma T: may one clique contain the germs of (tile t1, signs s1) and (tile t2, signs s2)?"""
    (i, j), (i2, j2) = t1, t2
    if (i, j) == (i2, j2): return True
    di, dj = i2 - i, j2 - j
    if abs(di) > 1 or abs(dj) > 1: return False
    if di != 0 and dj != 0: return False                  # diagonal: separated by (1, ±1)
    if di == 1: return s1[0] == 1 and s2[0] == -1
    if di == -1: return s2[0] == 1 and s1[0] == -1
    if dj == 1: return s1[1] == 1 and s2[1] == -1
    return s2[1] == 1 and s1[1] == -1


def _constraints(m):
    """the (tile, signs) constraints: sigma_x is forced at the wall columns, free inside"""
    out = []
    for i in range(m):
        for j in range(m):
            sxs = [1] if i == 0 else [-1] if i == m - 1 else [1, -1]
            sys_ = [1] if j == 0 else [-1] if j == m - 1 else [1, -1]
            for sx in sxs:
                for sy in sys_:
                    out.append(((i, j), (sx, sy)))
    return out


def _clique_columns(cons):
    """every maximal pairwise-compatible set of constraints (at most two tiles, by Lemma T)"""
    idx = {c: k for k, c in enumerate(cons)}
    cands = set()
    tiles = sorted({c[0] for c in cons})
    for t in tiles:                                        # single-tile cliques
        cands.add(frozenset(idx[c] for c in cons if c[0] == t))
    for t in tiles:                                        # two-tile cliques
        for t2 in tiles:
            if t2 <= t: continue
            S = [c for c in cons if c[0] in (t, t2)]
            grp = frozenset(idx[c] for c in S
                            if all(_compatible(c[0], c[1], d[0], d[1]) for d in S
                                   if idx[d] in [idx[e] for e in S]))
            # build properly: greedily keep the constraints that are pairwise compatible
            keep = []
            for c in S:
                if all(_compatible(c[0], c[1], d[0], d[1]) for d in keep):
                    keep.append(c)
            # only two-tile sets that really use both tiles are interesting
            if len({c[0] for c in keep}) == 2:
                cands.add(frozenset(idx[c] for c in keep))
            # the sign-symmetric variant (start from the other tile)
            keep = []
            for c in reversed(S):
                if all(_compatible(c[0], c[1], d[0], d[1]) for d in keep):
                    keep.append(c)
            if len({c[0] for c in keep}) == 2:
                cands.add(frozenset(idx[c] for c in keep))
            # and the canonical one: all (t, sx=+1) plus all (t2, sx=-1) when they are adjacent
            for A, Bt in ((t, t2), (t2, t)):
                di, dj = Bt[0] - A[0], Bt[1] - A[1]
                if (abs(di), abs(dj)) not in ((1, 0), (0, 1)): continue
                sel = []
                for c in cons:
                    if c[0] == A and ((di == 1 and c[1][0] == 1) or (dj == 1 and c[1][1] == 1)):
                        sel.append(c)
                    if c[0] == Bt and ((di == 1 and c[1][0] == -1) or (dj == 1 and c[1][1] == -1)):
                        sel.append(c)
                if len({c[0] for c in sel}) == 2:
                    cands.add(frozenset(idx[c] for c in sel))
    # drop non-maximal ones
    out = [c for c in cands if not any(c < d for d in cands)]
    return sorted(out, key=lambda s: (-len(s), sorted(s)))


def credit(a):
    from scipy.optimize import linprog
    m = a.m
    g = np.arange(0, m * 4 + 1) / 4          # a complete set of cell representatives (sec 2.3)
    GX, GY = np.meshgrid(g, g, indexing='ij')
    X, Y = GX.ravel(), GY.ravel()
    thetas = [0.0, 1e-5, 1e-4, 1e-3]
    cons = _constraints(m)
    A_pt = np.array([pin_mask(X, Y, m, t[0], t[1], thetas, s[0], s[1]).astype(float)
                     for t, s in cons])                     # (ncon, npoints)
    cols = _clique_columns(cons)
    A_cl = np.zeros((len(cons), len(cols)))
    for c, S in enumerate(cols):
        for r in S: A_cl[r, c] = 1.0
    corners = [(0, 0), (m - 1, 0), (0, m - 1), (m - 1, m - 1)]
    Ks = [int(x) for x in a.k.split(',')] if a.k else [0, 0, 0, 0]
    # variables: point weights (npoints), clique weights (ncols), lambda_1..4 (free)
    npt, ncl = A_pt.shape[1], len(cols)
    nl = 4
    Aub = np.hstack([-A_pt.T.T, -A_cl, np.zeros((len(cons), nl))])
    for r, (t, s) in enumerate(cons):
        if t in corners: Aub[r, npt + ncl + corners.index(t)] = 1.0   # ... - lambda_j >= 1
    bub = -np.ones(len(cons))
    c = np.r_[np.ones(npt), np.ones(ncl), -np.array(Ks, dtype=float)]
    # with no region trailer there are no multipliers at all: lambda = 0.  With one, lambda may
    # be any sign (FORMAT.md: "a multiplier ... may be negative").
    bounds = [(0, None)] * (npt + ncl) + ([(None, None)] * nl if a.k else [(0, 0)] * nl)
    res = linprog(c=c, A_ub=Aub, b_ub=bub, bounds=bounds, method='highs')
    print(f"container [0,{m}]^2, monotone-credited certificates")
    print(f"  {len(cons)} (tile, sign) constraints, {npt} point columns, {ncl} clique columns "
          f"(max size {max(len(s) for s in cols)}), "
          + (f"branch occupancy k = {Ks}" if a.k else "no region trailer (lambda = 0)"))
    print(f"  MINIMUM accounting total Theta = sum w_p + sum w_K - sum lambda_j K_j "
          f"= {res.fun:.6f}")
    lam = res.x[npt + ncl:]
    print(f"  lambda = {np.round(lam, 6).tolist()}")
    pw = res.x[:npt].sum(); cw = res.x[npt:npt + ncl].sum()
    print(f"  points {pw:.6f}, cliques {cw:.6f}")
    used = [(len(cols[k]), float(res.x[npt + k]), sorted(cons[r] for r in cols[k]))
            for k in range(ncl) if res.x[npt + k] > 1e-9]
    print(f"  clique columns used: {len(used)}")
    for sz, v, S in used[:12]:
        tt = sorted({t for t, _ in S})
        print(f"    w={v:.4f} over tiles {tt} ({sz} constraints)")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['bound', 'check', 'dual', 'credit'])
    ap.add_argument('cert', nargs='?')
    ap.add_argument('--m', type=int, default=4)
    ap.add_argument('--pitch', type=int, default=20)
    ap.add_argument('--octants', action='store_true')
    ap.add_argument('--switch', type=int, default=None,
                    help='dual: the last column with sigma_x = +1 (default 1: the switch sits '
                         'between columns 1 and 2, both interior, which Lemma 2 requires)')
    ap.add_argument('--k', type=str, default=None,
                    help="credit: the branch occupancy pattern K_1,..,K_4 (e.g. '1,1,1,1' for the "
                         "corner leaf k = 4; default '0,0,0,0', the no-branch case)")
    a = ap.parse_args()
    if a.mode == 'credit': credit(a)
    elif a.mode == 'bound': bound(a)
    elif a.mode == 'dual': dual(a)
    else: check(a.cert)


if __name__ == '__main__':
    main()
