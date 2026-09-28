#!/usr/bin/env python3
"""corner_deficit.py -- can a corner-local modification of Lebesgue measure save anything?  (No.)

Note: search/CORNER_DEFICIT.md.  2026-09-27.

The question (task corner-deficit): D*(R) = sup R^2 - mu([0,R]^2) over nonnegative measures mu on C = [0,R]^2 such
that every closed unit square S in the quarter plane Q = [0,inf)^2 captures mu(S cap C) + area(S minus C) >= 1,
i.e. mu(S cap C) >= area(S cap C).

Modes
  dual  R [--eps E]     PROVED part, exact arithmetic: the dilated grid family (pairwise disjoint closed unit squares
                        [i(1+e), i(1+e)+1] x [j(1+e), j(1+e)+1], 0 <= i, j < ceil R) and the uncovered area of C it
                        leaves; any feasible mu has mu(C) >= R^2 - uncovered, so D*(R) <= uncovered -> 0 as e -> 0.
  lp    R --q Q         HEURISTIC LP: columns = cell densities (pitch 1/Q), lattice points (pitch 1/Q), uniform
                        densities on the pieces of the grid lines x, y in {1, 2, ..} inside C; rows = ALL axis-parallel
                        closed unit squares with lattice corners (pitch 1/Q) meeting C.  Rows are a subset of the true
                        constraints, so the LP value is an UPPER bound on D over these columns; Lebesgue is feasible, so
                        it is >= 0.
  profile FILE [--rs ...]  EXACT deficit profile of a box certificate (plain or mixed format): for each R,
                        R^2 - mu([0,R]^2) (closed corner square, all four corners), and ring / centre deficits.
"""
import sys, os, math, argparse, time
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr


# ---------------------------------------------------------------- dual (proved)
def dual_bound(R, eps):
    """Exact: uncovered area of C=[0,R]^2 by the dilated family; also asserts disjointness / containment in Q."""
    R, eps = Fr(R), Fr(eps)
    n = math.ceil(R)
    starts = [i * (1 + eps) for i in range(n)]
    # 1-D: covered part of [0,R] by the intervals [s, s+1]; intervals pairwise disjoint (gap eps > 0)
    for a, b in zip(starts, starts[1:]):
        assert a + 1 < b
    cov = sum(max(Fr(0), min(s + 1, R) - s) for s in starts if s <= R)
    assert all(s <= R for s in starts), 'some square misses C (harmless, but ceil R should avoid it)'
    assert starts[-1] + 1 >= R, 'family does not reach x = R'
    unc1 = R - cov                      # uncovered length in [0,R]
    covered_area = cov * cov            # product family
    return dict(n=n, squares=n * n, uncovered_len=unc1, uncovered_area=R * R - covered_area,
                bound=R * R - covered_area)


# ---------------------------------------------------------------- LP (heuristic)
def build_lp(R, q):
    import numpy as np
    import scipy.sparse as sp
    r = round(R * q)
    assert abs(r - R * q) < 1e-9, 'R must be a multiple of 1/q'
    h = 1.0 / q
    cols, cost = [], []                 # each column: function (a, b) -> coefficient for the square [a h, a h + 1] x [b h, ..]
    # rows: lower-left corner (a, b), a, b in 0..r (square meets C iff a h <= R, b h <= R)
    rows = [(a, b) for a in range(r + 1) for b in range(r + 1)]
    ridx = {rb: k for k, rb in enumerate(rows)}
    I, J, V = [], [], []
    nc = 0

    def add_col(entries, c):
        nonlocal nc
        for (k, v) in entries:
            I.append(k); J.append(nc); V.append(v)
        cost.append(c)
        nc += 1

    # cells [c h, (c+1) h] x [d h, ..], density variable (mass = density * h^2); coefficient = overlap area
    # (aligned: the cell is inside the closed square iff a <= c <= a + q - 1)
    kinds = []
    for c in range(r):
        for d in range(r):
            ent = [(ridx[(a, b)], h * h) for a in range(max(0, c - q + 1), c + 1)
                   for b in range(max(0, d - q + 1), d + 1)]
            add_col(ent, h * h); kinds.append(('cell', c, d))
    # points (i h, j h), 0 <= i, j <= r: in the closed square iff a <= i <= a + q
    for i in range(r + 1):
        for j in range(r + 1):
            ent = [(ridx[(a, b)], 1.0) for a in range(max(0, i - q), min(i, r) + 1)
                   for b in range(max(0, j - q), min(j, r) + 1)]
            add_col(ent, 1.0); kinds.append(('pt', i, j))
    # line pieces on x = L (L = 1..floor R, L < R or = R), piece [t h, (t+1) h] in y, unit total mass
    for L in range(1, int(math.floor(R + 1e-9)) + 1):
        iL = L * q
        for t in range(r):
            ent = [(ridx[(a, b)], 1.0) for a in range(max(0, iL - q), min(iL, r) + 1)
                   for b in range(max(0, t - q + 1), t + 1)]
            add_col(ent, 1.0); kinds.append(('vx', L, t))
            ent = [(ridx[(b, a)], 1.0) for a in range(max(0, iL - q), min(iL, r) + 1)
                   for b in range(max(0, t - q + 1), t + 1)]
            add_col(ent, 1.0); kinds.append(('hy', L, t))
    A = sp.csr_matrix((V, (I, J)), shape=(len(rows), nc))
    # rhs: area(S cap C) for S = [a h, a h + 1] x [b h, b h + 1]
    rhs = np.array([min(1.0, R - a * h) * min(1.0, R - b * h) for (a, b) in rows])
    return A, np.array(cost), rhs, kinds, rows


def solve_lp(R, q):
    import numpy as np
    from scipy.optimize import linprog
    t0 = time.time()
    A, c, rhs, kinds, rows = build_lp(R, q)
    res = linprog(c, A_ub=-A, b_ub=-rhs, bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    x = res.x
    mass = float(c @ x)
    cap = A @ x - rhs
    by = {}
    for (k, *_), xi, ci in zip(kinds, x, c):
        by[k] = by.get(k, 0.0) + xi * ci
    ratio, worst = scan_axis(R, q, kinds, x, c, sub=4)
    return dict(R=R, q=q, rows=len(rows), cols=len(c), mass=mass, D=R * R - mass, minslack=float(cap.min()),
                by_kind=by, secs=time.time() - t0, ratio=ratio, worst=worst,
                D_rep=R * R - mass - max(0.0, 1 - 1 / ratio) * R * R)   # mu + delta Leb(C), delta = 1 - min(mu/area)


def _ov(a0, a1, b0, b1):
    return max(0.0, min(a1, b1) - max(a0, b0))


def scan_axis(R, q, kinds, x, c, sub=4):
    """Float scan of the LP solution over axis-parallel closed unit squares with lower-left corner on the pitch
    1/(sub q) lattice (off the LP rows).  Returns (max over squares of area(S cap C) / mu(S cap C), worst pose): the
    factor by which the modification must be scaled up to be valid on these poses (a necessary repair only)."""
    import numpy as np
    h = 1.0 / q
    r = round(R * q)
    dens = np.zeros((r, r)); pts = np.zeros((r + 1, r + 1)); vx = {}; hy = {}
    for (k, *ij), xi in zip(kinds, x):
        if k == 'cell': dens[ij[0], ij[1]] = xi
        elif k == 'pt': pts[ij[0], ij[1]] = xi
        elif k == 'vx': vx.setdefault(ij[0], np.zeros(r))[ij[1]] = xi
        else: hy.setdefault(ij[0], np.zeros(r))[ij[1]] = xi
    m = round(R * q * sub)
    g = np.arange(m + 1) / (q * sub)
    eps = 1e-12
    OX = np.array([[_ov(x0, x0 + 1, i * h, (i + 1) * h) for i in range(r)] for x0 in g])      # cell overlap lengths
    IX = np.array([[1.0 if x0 - eps <= i * h <= x0 + 1 + eps else 0.0 for i in range(r + 1)] for x0 in g])
    FX = OX / h                                                                                 # piece fractions
    M = OX @ dens @ OX.T + IX @ pts @ IX.T
    for L, v in vx.items():                       # vertical line x = L: indicator in x, fraction in y
        ind = np.array([1.0 if x0 - eps <= L <= x0 + 1 + eps else 0.0 for x0 in g])
        M += np.outer(ind, FX @ v)
    for L, v in hy.items():
        ind = np.array([1.0 if x0 - eps <= L <= x0 + 1 + eps else 0.0 for x0 in g])
        M += np.outer(FX @ v, ind)
    A = np.outer(np.minimum(1.0, R - g), np.minimum(1.0, R - g))
    rat = A / np.maximum(M, 1e-300)
    k = np.unravel_index(np.argmax(rat), rat.shape)
    return float(rat[k]), (float(g[k[0]]), float(g[k[1]]), float(M[k]), float(A[k]))


# ---------------------------------------------------------------- certificate profile (exact)
def _seg_frac_in_rect(seg, D, x0, y0, x1, y1):
    X0, Y0, X1, Y1, w = seg
    ax, ay, bx, by = Fr(X0, D), Fr(Y0, D), Fr(X1, D), Fr(Y1, D)
    lo, hi = Fr(0), Fr(1)
    for p, dp, lo_b, hi_b in ((ax, bx - ax, x0, x1), (ay, by - ay, y0, y1)):
        if dp == 0:
            if p < lo_b or p > hi_b:
                return Fr(0)
        else:
            t0, t1 = (lo_b - p) / dp, (hi_b - p) / dp
            if t0 > t1:
                t0, t1 = t1, t0
            lo, hi = max(lo, t0), min(hi, t1)
    return max(Fr(0), hi - lo)


def mass_in_rect(cover, x0, y0, x1, y1):
    """Exact mu of the CLOSED rectangle [x0,x1] x [y0,y1]."""
    D, W = cover['D'], cover['W']
    assert not cover['polygons']
    X0, Y0, X1, Y1 = x0 * D, y0 * D, x1 * D, y1 * D
    m = sum(w for (X, Y, w) in cover['points'] if X0 <= X <= X1 and Y0 <= Y <= Y1)
    m = Fr(m, W)
    for sgm in cover['segments']:
        f = _seg_frac_in_rect(sgm, D, x0, y0, x1, y1)
        if f:
            m += f * Fr(sgm[4], W)
    return m


def profile(path, rs):
    import mixed_cover as MC
    cv = MC.load(path)
    s = cv['s']
    tot = MC.total(cv)
    print('# %s  s=%s  total=%.6f  savings s^2-total=%.6f' % (os.path.basename(path), s, float(tot), float(s * s - tot)))
    print('# R    corner deficits R^2-mu([0,R]^2) at the 4 corners (closed)      [proved <= 0 for INTEGER R <= s-1]')
    out = []
    for R in rs:
        R = Fr(R)
        if R > s:
            continue
        ds = []
        for (x0, y0) in ((0, 0), (s - R, 0), (0, s - R), (s - R, s - R)):
            ds.append(R * R - mass_in_rect(cv, x0, y0, x0 + R, y0 + R))
        out.append((R, ds))
        print('%5.2f  ' % float(R) + '  '.join('%+.6f' % float(d) for d in ds))
    print('# centre [w, s-w]^2 deficit and ring (complement) deficit = savings - centre deficit')
    for w in (Fr(1, 2), Fr(3, 4), Fr(1), Fr(5, 4), Fr(3, 2), Fr(2)):
        if 2 * w >= s:
            continue
        L = s - 2 * w
        # closed centre; ring = box minus closed centre  (mu(box) - mu(centre))
        dc = L * L - mass_in_rect(cv, w, w, s - w, s - w)
        dr = (s * s - tot) - dc
        print('w=%.2f  centre %+.6f   ring %+.6f   [centre proved <= 0 for integer w >= 1]' % (float(w), float(dc), float(dr)))
    print('# strips: bottom strip [0,s] x [0,w] deficit')
    for w in (Fr(1, 2), Fr(1), Fr(3, 2), Fr(2)):
        if w > s:
            continue
        print('w=%.2f  strip %+.6f' % (float(w), float(s * w - mass_in_rect(cv, 0, 0, s, w))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode')
    ap.add_argument('arg')
    ap.add_argument('--eps', default='1/1000')
    ap.add_argument('--q', type=int, default=10)
    ap.add_argument('--rs', default='0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7')
    a = ap.parse_args()
    if a.mode == 'dual':
        for e in a.eps.split(','):
            d = dual_bound(Fr(a.arg), Fr(e))
            print('R=%s eps=%s: %d disjoint closed unit squares in Q, uncovered area of C = %s = %.3g  => D*(R) <= %.3g'
                  % (a.arg, e, d['squares'], d['bound'], float(d['bound']), float(d['bound'])))
    elif a.mode == 'lp':
        r = solve_lp(float(a.arg), a.q)
        print('R=%.2f q=%d rows=%d cols=%d  LP mass=%.6f  area=%.4f  D_LP=%.3e  (min slack %.1e, %.1fs)  mass by kind %s'
              % (r['R'], r['q'], r['rows'], r['cols'], r['mass'], r['R'] ** 2, r['D'], r['minslack'], r['secs'],
                 {k: round(float(v), 4) for k, v in r['by_kind'].items()}))
        w = r['worst']
        print('   off-lattice axis scan (pitch 1/%d): min capture ratio mu/area = %.4f at corner (%.4f, %.4f) '
              '[mu=%.4f area=%.4f];  repair mu + (1-min ratio) Leb(C): D = %.4f (on these poses only)' % (4 * a.q, 1 / r['ratio'], w[0], w[1], w[2], w[3], r['D_rep']))
    elif a.mode == 'profile':
        profile(a.arg, [Fr(x) for x in a.rs.split(',')])


if __name__ == '__main__':
    main()
