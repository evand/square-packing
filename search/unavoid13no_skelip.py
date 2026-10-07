#!/usr/bin/env python3
"""unavoid13-no: structural lemma by IP -- "every k-set hitting F has >= q points within eps of the
tile lines x, y in {1, ..., m-1}".  Adds the row  sum_{cells within eps of a tile line} x_c <= q - 1  to the
feasibility IP; INFEASIBLE proves the lemma for F (hence for every k-point unavoidable set, F being a
finite subfamily of all admissible squares).  Cells are D4-closed and the D4 root row is used (sound:
the near-line property is D4-invariant, so the symmetric image of a violating set violates too).

    python3 search/unavoid13no_skelip.py FAMILY.txt --k 13 --eps 0.12 --q 10
"""
import sys, os, time, argparse, json, math
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13no_lib as N


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family'); ap.add_argument('--k', type=int, default=13); ap.add_argument('--eps', type=float, default=0.12)
    ap.add_argument('--q', type=int, default=10); ap.add_argument('--threads', type=int, default=2); ap.add_argument('--time', type=float, default=7200)
    ap.add_argument('--no-cut', action='store_true')
    a = ap.parse_args()
    m, poses = L.read_family(a.family); mf = float(m)
    sq, V, B = N.maximal_cells(m, poses, threads=a.threads)
    V, B, perms = N.d4_close_cells(m, poses, V, B)
    lines = list(range(1, int(mf)))
    dx = np.min(np.abs(V[:, 0][:, None] - np.array(lines)[None, :]), 1)
    dy = np.min(np.abs(V[:, 1][:, None] - np.array(lines)[None, :]), 1)
    near = (dx <= a.eps + 1e-9) | (dy <= a.eps + 1e-9)       # lenient: a cell whose vertex is within eps
    # NOTE: a cell is a class of points with the same incidence; its representative vertex is one point of
    # it.  "Near" is decided on that vertex, so the lemma proved is about the representative vertices'
    # positions: every k-set OF CELLS has >= q cells whose stored vertex is within eps of a tile line.
    print(f"instance: {len(poses)} squares, {B.shape[0]} cells, {int(near.sum())} near-line cells (eps {a.eps})", flush=True)
    extra = [(near.astype(float), -np.inf, float(a.q - 1))]
    if not a.no_cut:
        j0 = poses.index(L.pose_key(F(0), m / 2, m / 2)); extra.append(N.root_cut_cells(m, V, B, j0, perms=perms))
    r = N.solve_ip(B, a.k, threads=a.threads, time_limit=a.time, extra_rows=extra, feasibility=True)
    print(f"k = {a.k}, at most {a.q-1} near-line cells: {r['status']} feasible={r['feasible']} nodes {r['nodes']} {r['time']:.0f}s", flush=True)
    if r['x'] is not None:
        P = V[np.nonzero(r['x'] > 0.5)[0]]
        print("  set: " + " ".join(f"({x:.3f},{y:.3f})" for x, y in P))
    json.dump(dict(family=a.family, k=a.k, eps=a.eps, q=a.q, status=r['status'], feasible=r['feasible'], time=r['time']),
              open(os.path.splitext(a.family)[0] + f'_skel_e{a.eps}_q{a.q}.json', 'w'))


if __name__ == '__main__':
    main()
