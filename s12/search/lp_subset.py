#!/usr/bin/env python3
"""lp_subset.py -- cut benchmark instances out of a pool dump (lp_dump_lattice.py / branch.py --dump-lp).

    python3 search/lp_subset.py POOL --out OUT --rows R [--drop-cols C]
        the first R rows, all point columns except the last C (lambda columns kept)
    python3 search/lp_subset.py POOL --round OUT --rows R --add A --drop-cols C [--threads 8]
        one cutting-plane round: OUT_it0 = first R rows without the last C point columns; solve it
        (highspy dual simplex); OUT_it1 = OUT_it0's rows + the A rows of the pool most violated by that
        solution + all columns.  it0 ⊂ it1 exactly as branch.py appends rows and columns.
"""
import sys, os, argparse, time
import numpy as np, scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lp_bench as LB


def sub(d, rows, cols):
    """rows: index array; cols: index array over ALL columns (points + lambda)"""
    e = {}
    e['M'] = d['M'][rows][:, cols].tocsr()
    for k in ('c', 'lo', 'hi'): e[k] = d[k][cols]
    e['b'] = d['b'][rows]; e['rows'] = d['rows'][rows]; e['rflag'] = d['rflag'][rows]; e['rowkeys'] = d['rowkeys'][rows]
    e['n'] = len(cols) - int(d['nl']); e['nl'] = int(d['nl']); e['nr'] = len(rows)
    for k in ('K', 'D', 'sym', 'kvec', 'r'): e[k] = d[k]
    e['val'] = np.nan; e['x'] = np.zeros(0); e['lam'] = np.zeros(0); e['y'] = np.zeros(0)
    return e


def write(e, path):
    sp.save_npz(path + "_A.npz", e['M'], compressed=False)
    np.savez(path + ".npz", **{k: v for k, v in e.items() if k != 'M'})
    print(f"wrote {path}: rows={e['M'].shape[0]} cols={e['M'].shape[1]} nnz={e['M'].nnz}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pool'); ap.add_argument('--out', default=None); ap.add_argument('--round', default=None)
    ap.add_argument('--rows', type=int, required=True); ap.add_argument('--drop-cols', type=int, default=0); ap.add_argument('--add', type=int, default=3000)
    ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    d = LB.load(a.pool); nr, nc = d['M'].shape; nl = int(d['nl']); n = nc - nl
    rows0 = np.arange(min(a.rows, nr))
    cols0 = np.concatenate([np.arange(n - a.drop_cols), np.arange(n, nc)])
    if a.out:
        write(sub(d, rows0, cols0), a.out); return
    e0 = sub(d, rows0, cols0); write(e0, a.round + "_it0")
    t = time.time(); r = LB.run_hs(e0, 'simplex', a.threads, 1e-7)
    print(f"it0 solved: obj={r['obj']:.9f} t={r['t']:.1f}s it={r['it']} {r['status']}")
    z = np.zeros(nc); z[cols0] = r['z']
    slack = d['M'] @ z - d['b']; slack[rows0] = np.inf
    viol = np.nonzero(slack < -1e-9)[0]; print(f"violated among the remaining {nr - len(rows0)} pool rows: {len(viol)} (worst {slack[viol].min() if len(viol) else 0:.3e})")
    pick = viol[np.argsort(slack[viol])[:a.add]]
    rows1 = np.concatenate([rows0, np.sort(pick)])
    e1 = sub(d, rows1, np.arange(nc)); e1['val'] = np.nan; write(e1, a.round + "_it1")
    # record it0's solution in it0's dump (for sift --warm-from: tight rows)
    e0['x'] = r['z'][:e0['n']]; e0['lam'] = r['z'][e0['n']:]; e0['y'] = r['y']; e0['val'] = r['obj']; write(e0, a.round + "_it0")


if __name__ == '__main__':
    main()
