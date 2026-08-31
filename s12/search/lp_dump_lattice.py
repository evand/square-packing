#!/usr/bin/env python3
"""lp_dump_lattice.py -- build a large branch.py LP instance without running the loop, for lp_bench.py.

Columns = exactly the points of a column checkpoint (runs/branch_TAG_cols.txt; one column per line,
NOT symmetrised -- for an asymmetric leaf that file *is* the column set of the run), weights from a
probe/certificate file for the lattice filter.  Rows = lattice poses capturing < thr under those
weights (branch.lattice_rows: h = 1/2, closed unit squares, [0,90deg) for asymmetric leaves),
`perang` lowest-capture poses per angle.

    python3 search/lp_dump_lattice.py PROBE --cols COLS --k 1110 --eta 0.005 --dt 0.01 --perang 1000 --out runs/lp1110_pool

`--sub OUT ROWS COLS` afterwards (lp_bench.py --subset) cuts sub-instances with the first ROWS rows and
the first COLS point columns (row keys and the lambda columns are kept), so that a pair of dumps
(it0 ⊂ it1) mimics one round of the cutting-plane loop: rows appended, columns appended.
"""
import sys, os, argparse, time
import numpy as np, scipy.sparse as sp
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import branch as B
import tighten as T


def build_from_cols(probe, cols, r, kreg):
    sn, sd, D, WD, rows = T.read_cert(probe)
    K = Fraction(sn, sd) * D; assert K.denominator == 1; K = int(K)
    with open(cols) as f:
        K2, D2, sym = (int(v) for v in f.readline().split()); assert (K2, D2) == (K, D)
        pts = [tuple(int(v) for v in l.split()) for l in f if l.strip()]
    okey = {}; orbs = []
    for p in pts:
        if p in okey: continue
        okey[p] = len(orbs); orbs.append(np.array([p], dtype=np.int64))
    w0 = np.zeros(len(orbs))
    for X, Y, W in rows:
        p = (int(X), int(Y))
        if p not in okey: okey[p] = len(orbs); orbs.append(np.array([p], dtype=np.int64)); w0 = np.append(w0, 0.0)
        w0[okey[p]] += W / WD
    m = B.BModel(K, D, orbs, r, kreg, sym=False); m.okey = okey
    return m, w0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('probe'); ap.add_argument('--cols', required=True); ap.add_argument('--k', required=True); ap.add_argument('--r', default='1')
    ap.add_argument('--eta', type=float, default=0.005); ap.add_argument('--dt', type=float, default=0.01); ap.add_argument('--perang', type=int, default=1000)
    ap.add_argument('--thr', type=float, default=3.0); ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--margin', type=float, default=2e-6)
    ap.add_argument('--row-grid', type=float, nargs=2, default=(0.002, 0.005)); ap.add_argument('--shuffle', type=int, default=0, help='shuffle rows with this seed before dumping (so that a prefix is a uniform sample)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    t = time.time()
    m, w0 = build_from_cols(a.probe, a.cols, Fraction(a.r), B.parse_k(a.k)); m.row_grid = tuple(a.row_grid)
    print(f"model: {len(m.orbits)} columns, container {m.K}/{m.D}, k={m.kvec}, weight {float(m.sizes @ w0):.6f}", flush=True)
    lr = B.lattice_rows(m, w0[m.own], eta=a.eta, dt=a.dt, thr=a.thr, perang=a.perang, nproc=a.nproc)
    print(f"lattice: {len(lr)} candidate rows ({time.time()-t:.0f}s)", flush=True)
    if a.shuffle:
        rng = np.random.default_rng(a.shuffle); lr = [lr[i] for i in rng.permutation(len(lr))]
    n = m.add_rows(lr); print(f"added {n} -> {len(m.rows)} rows ({time.time()-t:.0f}s)", flush=True)
    B.write_lp_dump(m, a.margin, a.out)
    M = m.matrix(); print(f"wrote {a.out}: rows={M.shape[0]} cols={M.shape[1]+m.nl} nnz={M.nnz} ({time.time()-t:.0f}s)")


if __name__ == '__main__':
    main()
