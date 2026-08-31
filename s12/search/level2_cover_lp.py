#!/usr/bin/env python3
"""Cover-side level-2 leaf value on the `k = 4` leaf's own saved instance.

`branch.py`'s LP is  min W - sum_j lam_j K_j  subject to "every pose row captures
`1 + sum_j lam_j [row in region j]`", over a sampled set of pose rows and a point column set.  The
`k = 4` run `t398hk4` saved both: `runs/branch_t398hk4_cols.txt` (the point columns, integers over
`D`, one representative per D4 orbit) and `runs/branch_t398hk4_dual_it16.txt` (the pose rows that
carry dual mass, with the verifier's shrunk half-side `h`).  This script rebuilds that LP with
**per-region multipliers** (4 corner boxes + 8 wall slots) and re-solves it for a given occupancy
pattern, so the value is directly comparable with the published `12.000024` of the `k = 4` leaf.

Semantics is the verifier's, exactly as in `branch.py`: a row is the concentric shrunk square of
half-side `h` at the bin's base angle, and a point is captured iff it lies in that square.  The
certificate is asymmetric (unequal multipliers), so every point is its own variable and every
dihedral image of a pose is its own row.

Direction of the error: the pose ROWS are a sample (only the rows that carried dual mass in the
`k = 4` run), so the minimum is **too small** -- an optimistic leaf value.  A value `>= 12` is
therefore conclusive ("this leaf does not close"), a value `< 12` is not.  The point columns are
the full saved set, which is the same set the `12.000024` was computed with, so the comparison
between patterns on this fixed instance is meaningful even though the level is not.

Usage
    python3 search/level2_cover_lp.py --cols RUNS/branch_t398hk4_cols.txt \\
        --rows RUNS/branch_t398hk4_dual_it16.txt --patterns ........,11111111,01010101,...
"""
import argparse, math, os, sys, time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
for _alt in ('/home/evand/math/square-packing/s12/search',):
    if os.path.exists(_alt):
        sys.path.append(_alt)
from level2_regions import classify, images                      # noqa: E402


def read_cols(path):
    with open(path) as f:
        K, D, sym = f.readline().split()
        K, D, sym = int(K), int(D), int(sym)
        pts = np.array([[int(a) for a in line.split()] for line in f if line.strip()], dtype=float)
    return K / D, pts / D, bool(sym)


def read_rows(path):
    t = None
    rows = []
    for line in open(path):
        if line.startswith('#'):
            for tok in line.split():
                if tok.startswith('s=') and '/' in tok:
                    a, b = tok[2:].split('/')
                    t = float(a) / float(b)
            continue
        q = line.split()
        if len(q) == 6:
            rows.append((float(q[0]), float(q[1]), float(q[2]), float(q[3]), float(q[5])))
    return t, rows


def expand(t, items, kind):
    """all 8 dihedral images"""
    out = []
    for it in items:
        if kind == 'pt':
            x, y = it
            out += [(x, y), (t - x, y), (x, t - y), (t - x, t - y),
                    (y, x), (t - y, x), (y, t - x), (t - y, t - x)]
        else:
            cx, cy, th, h = it[:4]
            out += [(x, y, u, h) for (x, y, u) in images(cx, cy, th, t)]
    return out


def incidence(rows, pts, chunk=200):
    """rows x pts 0/1: point in the shrunk square of half-side h"""
    R = np.asarray(rows, dtype=float)
    P = np.asarray(pts, dtype=float)
    blocks = []
    for i in range(0, len(R), chunk):
        r = R[i:i + chunk]
        ct = np.cos(r[:, 2])[:, None]
        st = np.sin(r[:, 2])[:, None]
        dx = P[None, :, 0] - r[:, 0][:, None]
        dy = P[None, :, 1] - r[:, 1][:, None]
        u = dx * ct + dy * st
        v = -dx * st + dy * ct
        hh = r[:, 3][:, None]
        blocks.append(sp.csr_matrix((np.abs(u) <= hh + 1e-12) & (np.abs(v) <= hh + 1e-12)))
    return sp.vstack(blocks, format='csr').astype(float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cols', required=True)
    ap.add_argument('--rows', required=True)
    ap.add_argument('--r', type=float, default=1.0)
    ap.add_argument('--corners', default='1111')
    ap.add_argument('--patterns', default='........,11111111,01111111,01010101,00110011,00001111,00000011,00000000')
    ap.add_argument('--method', default='highs')
    ap.add_argument('--col-stride', type=int, default=1,
                    help='keep every k-th point orbit as a column (a coarser cover: values err UP)')
    ap.add_argument('--row-top', type=int, default=0,
                    help='keep only the heaviest N pose orbits of the dual as rows (0 = all)')
    ap.add_argument('--sym-corner', action='store_true',
                    help='one shared multiplier for the four corner boxes (as branch.py does)')
    a = ap.parse_args()

    t1, pts0, sym = read_cols(a.cols)
    t2, rows0 = read_rows(a.rows)
    t = t2 or t1
    if a.col_stride > 1:
        pts0 = pts0[::a.col_stride]
    pts = np.array(sorted(set(map(tuple, np.round(expand(t, pts0, 'pt'), 9)))))
    if a.row_top:
        rows0 = sorted(rows0, key=lambda z: -z[4])[:a.row_top]
    rows = list(dict.fromkeys(map(tuple, np.round(expand(t, rows0, 'row'), 9))))
    print(f't = {t:.9f}   points {len(pts0)} orbits -> {len(pts)} columns   '
          f'poses {len(rows0)} -> {len(rows)} rows')

    reg = []
    for (cx, cy, th, h) in rows:
        k, j = classify(cx, cy, t, a.r)
        reg.append(j if k == 'C' else (4 + j if k == 'W' else 12))
    reg = np.array(reg)
    print('   rows by region: corners %d, slots %s, interior %d'
          % ((reg < 4).sum(), [int((reg == 4 + j).sum()) for j in range(8)], (reg == 12).sum()))

    t0 = time.time()
    A = incidence(rows, pts)
    empty = int((np.asarray(A.sum(axis=1)).ravel() == 0).sum())
    print(f'   incidence {A.shape} nnz={A.nnz} empty rows {empty} ({time.time()-t0:.0f}s)')

    kc = [None if c == '.' else float(c) for c in a.corners]
    for pat in a.patterns.split(','):
        kw = [None if c == '.' else float(c) for c in pat]
        # variables: w (points) then the free multipliers lam (split into lam+ - lam-)
        cols = []
        counts = []
        if a.sym_corner and all(v is not None for v in kc):
            cols.append(((reg < 4)).astype(float)); counts.append(sum(kc))
        else:
            for i in range(4):
                if kc[i] is not None:
                    cols.append((reg == i).astype(float)); counts.append(kc[i])
        for j in range(8):
            if kw[j] is not None:
                cols.append((reg == 4 + j).astype(float)); counts.append(kw[j])
        nl = len(cols)
        L = np.array(cols).T if nl else np.zeros((len(rows), 0))
        # min sum w - sum_j lam_j K_j   s.t.  A w - L lam >= 1
        Aub = sp.hstack([-A, sp.csr_matrix(L), sp.csr_matrix(-L)], format='csr')
        c = np.concatenate([np.ones(len(pts)), -np.array(counts), np.array(counts)]) if nl else np.ones(len(pts))
        res = linprog(c=c, A_ub=Aub, b_ub=-np.ones(len(rows)), bounds=(0, None), method=a.method)
        if not res.success:
            print(f'   pattern {pat}: LP failed ({res.message[:60]})')
            continue
        W = float(res.x[:len(pts)].sum())
        lam = res.x[len(pts):len(pts) + nl] - res.x[len(pts) + nl:] if nl else np.zeros(0)
        print(f'   pattern {pat:>8}  value = {res.fun:.6f}   W = {W:.6f}   '
              f'lam = [{" ".join(f"{v:+.3f}" for v in lam)}]')


if __name__ == '__main__':
    main()
