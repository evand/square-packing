#!/usr/bin/env python3
"""unavoid13-no step 4: a purpose-built branch-and-bound for "is there a k-point hitting set of F?"

    python3 search/unavoid13no_bb.py FAMILY.txt [--k 13] [--threads 2] [--no-cut] [--node-limit N]

Instance: maximal cells of the family (unavoid13no_lib.maximal_cells, tol 1e-11), D4-closed.
Search: node = set of chosen cells; rows hit so far.  Branch on an UNHIT ROW with the fewest cells
(every hitting set contains one of its cells); children = those cells, LP-value-first.  Bound at
each node: the LP relaxation of the residual cover (all cells, unhit rows only; highspy simplex,
warm-started by changing row bounds), pruned when  depth + ceil(LP - 1e-6) > k;  plus a free greedy
packing bound (pairwise cell-disjoint unhit rows) tried first.

D4 root rule [proved in unavoid13no_lib docstring]: at the root the branching row is the
D4-invariant square Q0 (axis-parallel at the centre); its cells form a D4-invariant set and only one
representative per D4-orbit is branched on (if a hitting set X exists, some gX contains a
representative, and gX is a hitting set because F is D4-closed and the cells are D4-closed).
Nothing else is assumed about symmetry deeper in the tree.
"""
import sys, os, time, math, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13no_lib as N


class BB:
    def __init__(self, B, k, threads=2, perms=None, j0=None, verbose=True, node_limit=10 ** 7,
                 time_limit=None, lp_every=1):
        import highspy
        self.B = B; self.K, self.n = B.shape; self.k = k; self.verbose = verbose
        self.node_limit = node_limit; self.time_limit = time_limit; self.lp_every = lp_every
        self.rows_of = [np.nonzero(B[c])[0] for c in range(self.K)]
        self.cells_of = [np.nonzero(B[:, r])[0] for r in range(self.n)]
        self.rowcnt = B.sum(0)
        self.perms = perms; self.j0 = j0
        starts = np.zeros(self.K + 1, dtype=np.int64); starts[1:] = np.cumsum(B.sum(1))
        index = np.concatenate(self.rows_of).astype(np.int32)
        h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', int(threads))
        h.setOptionValue('presolve', 'off')
        lp = highspy.HighsLp(); lp.num_col_ = self.K; lp.num_row_ = self.n
        lp.col_cost_ = np.ones(self.K); lp.col_lower_ = np.zeros(self.K); lp.col_upper_ = np.ones(self.K)
        lp.row_lower_ = np.ones(self.n); lp.row_upper_ = np.full(self.n, highspy.kHighsInf)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = np.ones(len(index))
        h.passModel(lp); self.h = h; self.hs = highspy
        self.nodes = 0; self.lps = 0; self.pruned_lp = 0; self.pruned_pack = 0; self.pruned_dom = 0; self.last_rc = None; self.t0 = time.time()
        self.solution = None; self.maxdepth = 0

    def lp_bound(self, cnt):
        """LP over unhit rows (cnt == 0): set row lower bounds 1 for unhit, 0 for hit."""
        lo = (cnt == 0).astype(float)
        self.h.changeRowsBounds(self.n, np.arange(self.n, dtype=np.int32), lo, np.full(self.n, self.hs.kHighsInf))
        self.h.run(); self.lps += 1
        if self.h.getModelStatus() != self.hs.HighsModelStatus.kOptimal:
            return float('inf'), None
        sol = self.h.getSolution()
        self.last_rc = np.array(sol.col_dual)
        return self.h.getInfo().objective_function_value, np.array(sol.col_value)

    def pack_bound(self, cnt):
        """Greedy set of pairwise cell-disjoint unhit rows (each needs its own point)."""
        unhit = np.nonzero(cnt == 0)[0]
        if len(unhit) == 0: return 0
        order = unhit[np.argsort(self.rowcnt[unhit])]
        used = np.zeros(self.K, dtype=bool); m = 0
        for r in order:
            cs = self.cells_of[r]
            if not used[cs].any():
                used[cs] = True; m += 1
        return m

    def run(self):
        cnt = np.zeros(self.n, dtype=np.int32)
        try:
            self._node([], cnt, 0)
        except KeyboardInterrupt:
            pass
        return self.solution

    def _node(self, chosen, cnt, depth):
        self.nodes += 1
        if self.nodes > self.node_limit: raise KeyboardInterrupt
        if self.time_limit and time.time() - self.t0 > self.time_limit: raise KeyboardInterrupt
        self.maxdepth = max(self.maxdepth, depth)
        unhit = np.nonzero(cnt == 0)[0]
        if len(unhit) == 0:
            self.solution = list(chosen); raise KeyboardInterrupt
        if depth >= self.k: return
        rem = self.k - depth
        if self.pack_bound(cnt) > rem:
            self.pruned_pack += 1; return
        x = None; val = None
        if depth % self.lp_every == 0:
            val, x = self.lp_bound(cnt)
            if math.ceil(val - 1e-6) > rem:
                self.pruned_lp += 1; return
            rc = self.last_rc
        if self.verbose and self.nodes % 200 == 0:
            print(f"    nodes {self.nodes}, lps {self.lps}, depth {depth}, pruned lp {self.pruned_lp} pack {self.pruned_pack}, {time.time()-self.t0:.0f}s", flush=True)
        # branching row: unhit row with fewest cells (root: Q0 with orbit representatives)
        if depth == 0 and self.j0 is not None and self.perms is not None and cnt[self.j0] == 0:
            r = self.j0
            cells = [c for c in self.cells_of[r] if c == self.perms[:, c].min()]
        else:
            r = unhit[np.argmin(self.rowcnt[unhit])]
            cells = list(self.cells_of[r])
        if x is not None:
            # reduced-cost pruning: fixing x_c = 1 raises the LP bound by at least rc(c) >= 0
            cells = [c for c in cells if math.ceil(val + max(rc[c], 0.0) - 1e-6) <= rem]
            cells.sort(key=lambda c: -x[c])
        # residual dominance among the children: drop c if (S(c) ∩ unhit) ⊆ (S(c') ∩ unhit) for an
        # earlier-kept c' (a hitting set using c can use c' instead)
        if len(cells) > 1:
            unhitmask = (cnt == 0)
            res = [self.B[c] & unhitmask for c in cells]
            resp = [np.packbits(r) for r in res]
            keep = []
            for i, c in enumerate(cells):
                dom = False
                for j in keep:
                    if not (resp[i] & ~resp[j]).any():
                        dom = True; break
                if not dom: keep.append(i)
            self.pruned_dom += len(cells) - len(keep)
            cells = [cells[i] for i in keep]
        for c in cells:
            cnt[self.rows_of[c]] += 1
            self._node(chosen + [int(c)], cnt, depth + 1)
            cnt[self.rows_of[c]] -= 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family'); ap.add_argument('--k', type=int, default=13); ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--no-cut', action='store_true'); ap.add_argument('--node-limit', type=int, default=10 ** 7)
    ap.add_argument('--time-limit', type=float, default=None); ap.add_argument('--highs', action='store_true', help='also run highspy on the same instance')
    ap.add_argument('--lp-every', type=int, default=1)
    a = ap.parse_args()
    m, poses = L.read_family(a.family)
    sq, V, B = N.maximal_cells(m, poses, threads=a.threads)
    V, B, perms = N.d4_close_cells(m, poses, V, B)
    j0 = poses.index(L.pose_key(F(0), m / 2, m / 2)) if not a.no_cut else None
    print(f"instance: {len(poses)} squares, {B.shape[0]} cells; root row Q0 has {int(B[:, j0].sum()) if j0 is not None else '-'} cells", flush=True)
    bb = BB(B, a.k, threads=a.threads, perms=perms if not a.no_cut else None, j0=j0, node_limit=a.node_limit, time_limit=a.time_limit, lp_every=a.lp_every)
    t0 = time.time(); sol = bb.run(); dt = time.time() - t0
    verdict = ('FEASIBLE' if sol is not None else ('INFEASIBLE' if bb.nodes <= a.node_limit and (a.time_limit is None or dt < a.time_limit) else 'UNDECIDED'))
    print(f"B&B: {verdict}, nodes {bb.nodes}, LPs {bb.lps}, pruned lp {bb.pruned_lp} pack {bb.pruned_pack} dom {bb.pruned_dom}, max depth {bb.maxdepth}, {dt:.0f}s", flush=True)
    if sol is not None:
        assert B[sol].any(0).all()
        print("  cells: " + " ".join(f"({V[c][0]:.3f},{V[c][1]:.3f})" for c in sol))
    if a.highs:
        extra = [N.root_cut_cells(m, V, B, j0, perms=perms)] if j0 is not None else []
        r = N.solve_ip(B, a.k, threads=a.threads, extra_rows=extra, feasibility=False)
        print(f"highspy: {r['status']} feasible={r['feasible']} nodes {r['nodes']} {r['time']:.0f}s")
    json.dump(dict(family=a.family, k=a.k, verdict=verdict, nodes=bb.nodes, lps=bb.lps, time=dt),
              open(os.path.splitext(a.family)[0] + f'_bb{a.k}.json', 'w'))


if __name__ == '__main__':
    main()
