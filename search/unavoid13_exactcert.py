#!/usr/bin/env python3
"""Generate an exact branch-and-bound infeasibility certificate for an unavoid13 hitting-set claim.

    python3 search/unavoid13_exactcert.py FAMILY.txt --claim plain --k 6  --out CERT.txt
    python3 search/unavoid13_exactcert.py FAMILY.txt --claim C2    --k 13 --out CERT.txt

The instance (exact vertices, exact incidence, candidate reduction) is built with the existing
exploration code (unavoid13_recheck / unavoid13_lib); the certificate is checked by
search/unavoid13_exactcheck.py, which recomputes all of the geometry itself.

Candidates.  plain: distinct incidence sets of the exact arrangement vertices, dominated ones
dropped.  C2: the half-turn closure of the family is required; columns are orbit incidences
(p in Q or s(p) in Q), cost 2, plus the centre (m/2, m/2) with cost 1; orbit columns dominated by
another cost-2 column or by the centre's column are dropped.

Search.  Depth-first B&B on the 0/1 candidate variables.  At each node the LP
min c.x, A x >= 1, bounds from the fixings, is solved by highspy (floats); its row duals are
clamped to >= 0 and rounded down to multiples of 1/DEN, and the leaf is closed only if the EXACT
bound computed from those rationals (same formula as the checker) excludes cost <= k, using the
gcd rounding of the checker.  Otherwise branch (strong branching over the most promising
fractional variables; for C2 the centre first).  Uncoverable rows give `U` leaves.
"""
import sys, os, time, math, argparse, json, hashlib
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L
import unavoid13_recheck as R
import highspy


def build_plain(m, poses):
    V = R.exact_vertices(m, poses)
    B = R.exact_incidence(m, poses, V)
    uniq, idx = R.reduce_candidates(B)
    pts = [V[i] for i in idx]
    return uniq.astype(np.int64), np.ones(len(pts), dtype=np.int64), pts


def build_c2(m, poses):
    m = Fr(m)
    pindex = {p: i for i, p in enumerate(poses)}
    img = [L.pose_key(u, m - cx, m - cy) for (u, cx, cy) in poses]
    assert all(q in pindex for q in img), "family is not closed under the half-turn"
    perm = np.array([pindex[q] for q in img])
    centre = (m / 2, m / 2)
    V = sorted(set(R.exact_vertices(m, poses)) | {centre})
    B = R.exact_incidence(m, poses, V)
    OB = B | B[:, perm]
    ic = V.index(centre)
    cen_col = OB[ic].copy()
    rest = np.array([i for i in range(len(V)) if i != ic])
    uniq, idx = R.reduce_candidates(OB[rest])
    idx = rest[idx]
    # drop cost-2 columns dominated by the centre's column
    keep = [i for i in idx if not ((OB[i] & ~cen_col) == 0).all()]
    cols = [ic] + keep
    A = OB[cols].astype(np.int64)
    cost = np.array([1] + [2] * len(keep), dtype=np.int64)
    pts = [V[i] for i in cols]
    return A, cost, pts


class BB:
    def __init__(self, A, cost, k, den, threads, strong, first=None, verbose=True):
        self.A = A; self.c = cost; self.k = k; self.den = den; self.strong = strong
        self.first = first
        N, Rr = A.shape; self.N = N; self.R = Rr
        self.colsets = [np.nonzero(A[:, i])[0] for i in range(Rr)]
        h = highspy.Highs(); h.setOptionValue('output_flag', False)
        h.setOptionValue('threads', threads)
        lp = highspy.HighsLp(); lp.num_col_ = N; lp.num_row_ = Rr
        lp.col_cost_ = cost.astype(float); lp.col_lower_ = np.zeros(N); lp.col_upper_ = np.ones(N)
        lp.row_lower_ = np.ones(Rr); lp.row_upper_ = np.full(Rr, highspy.kHighsInf)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        cnt = A.sum(1)
        lp.a_matrix_.start_ = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int32)
        lp.a_matrix_.index_ = np.concatenate([np.nonzero(A[j])[0] for j in range(N)]).astype(np.int32)
        lp.a_matrix_.value_ = np.ones(int(cnt.sum()))
        h.passModel(lp); self.h = h
        self.state = np.zeros(N, dtype=np.int8)
        self.out = []; self.nodes = 0; self.lps = 0; self.t0 = time.time(); self.verbose = verbose
        self.minslack = None

    def solve(self, state):
        lo = (state == 1).astype(float); hi = (state != -1).astype(float)
        self.h.changeColsBounds(self.N, np.arange(self.N, dtype=np.int32), lo, hi)
        self.h.run(); self.lps += 1
        st = self.h.getModelStatus()
        if st == highspy.HighsModelStatus.kInfeasible:
            return None, None, None
        assert st == highspy.HighsModelStatus.kOptimal, self.h.modelStatusToString(st)
        sol = self.h.getSolution()
        return self.h.getInfo().objective_function_value, np.array(sol.col_value), np.array(sol.row_dual)

    def exact_close(self, state, ydual):
        """rationalise y; return (closed, L as Fraction, ynum)."""
        for den in (self.den, self.den * 1000):
            ynum = np.floor(np.maximum(ydual, 0) * den).astype(np.int64)
            rD = self.c * den - self.A @ ynum
            f1 = state == 1; fr = state == 0
            LD = int(ynum.sum()) + int(rD[f1].sum()) + int(np.minimum(rD[fr], 0).sum())
            Lb = Fr(LD, den)
            C1 = int(self.c[f1].sum())
            if fr.any():
                g = 0
                for cj in set(self.c[fr].tolist()): g = math.gcd(g, int(cj))
                least = C1 + g * max(0, math.ceil((Lb - C1) / g))
            else:
                least = C1 if Lb <= C1 else math.inf
            if least > self.k:
                return True, Lb, ynum, den
        return False, Lb, ynum, den

    def uncoverable(self, state):
        for i in range(self.R):
            if (state[self.colsets[i]] == -1).all():
                return i
        return None

    def node(self, depth=0):
        self.nodes += 1
        if self.verbose and self.nodes % 500 == 0:
            print(f"    {self.nodes} nodes, {self.lps} LPs, {len(self.out)} records, depth {depth}, "
                  f"{time.time()-self.t0:.0f}s", flush=True)
        st = self.state
        i = self.uncoverable(st)
        if i is not None:
            self.out.append(f"U {i}"); return
        val, x, y = self.solve(st)
        if val is None:
            raise RuntimeError("LP infeasible without an uncoverable row")
        closed, Lb, ynum, den = self.exact_close(st, y)
        if closed:
            self.out.append("D " + str(den) + " " + " ".join(f"{i}:{int(v)}" for i, v in enumerate(ynum) if v > 0))
            sl = float(Lb) - self.k
            self.minslack = sl if self.minslack is None else min(self.minslack, sl)
            return
        j = self.choose(st, x)
        if j is None:
            raise RuntimeError(f"integral LP solution of cost {val} <= k: the claim is FALSE")
        self.out.append(f"B {j}")
        st[j] = 1; self.node(depth + 1)
        st[j] = -1; self.node(depth + 1)
        st[j] = 0

    def choose(self, st, x):
        if self.first is not None and st[self.first] == 0:
            return self.first
        free = np.nonzero(st == 0)[0]
        xf = x[free]
        frac = np.minimum(xf, 1 - xf)
        cand = free[frac > 1e-7]
        if len(cand) == 0:
            # integral LP solution; if its cost > k we would have closed; so it is a counterexample
            return None
        if self.strong <= 1:
            return int(cand[np.argmax(x[cand] * 0 + np.minimum(x[cand], 1 - x[cand]))])
        # strong branching over the top candidates by LP value
        order = cand[np.argsort(-x[cand])][:self.strong]
        best = None; bestsc = -1
        for j in order:
            sc = []
            for v in (1, -1):
                st[j] = v
                if self.uncoverable(st) is not None:
                    sc.append(1e9)
                else:
                    val, _, _ = self.solve(st)
                    sc.append(1e9 if val is None else val)
                st[j] = 0
            s = min(sc) * 1000 + max(sc)
            if s > bestsc: bestsc = s; best = int(j)
        return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family'); ap.add_argument('--claim', required=True, choices=['plain', 'C2'])
    ap.add_argument('--k', type=int, required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--den', type=int, default=10 ** 6)
    ap.add_argument('--strong', type=int, default=8)
    ap.add_argument('--threads', type=int, default=1)
    a = ap.parse_args()
    sys.setrecursionlimit(100000)
    t0 = time.time()
    m, poses = L.read_family(a.family)
    if a.claim == 'plain':
        A, cost, pts = build_plain(m, poses); first = None
    else:
        A, cost, pts = build_c2(m, poses); first = 0
    print(f"  instance: {A.shape[0]} candidates x {A.shape[1]} rows, costs "
          f"{dict(zip(*np.unique(cost, return_counts=True)))}, {time.time()-t0:.1f}s", flush=True)
    bb = BB(A, cost, a.k, a.den, a.threads, a.strong, first)
    val, _, _ = bb.solve(bb.state)
    print(f"  root LP {val:.6f}", flush=True)
    bb.node()
    nb = sum(1 for r in bb.out if r[0] == 'B')
    print(f"  tree: {len(bb.out)} records, {nb} branchings, {len(bb.out) - nb} leaves, "
          f"{bb.lps} LPs, min slack {bb.minslack}, {time.time()-t0:.1f}s", flush=True)
    sha = hashlib.sha256(open(a.family, 'rb').read()).hexdigest()
    with open(a.out, 'w') as f:
        f.write("# unavoid13 exact hitting-set infeasibility certificate "
                "(format and meaning: search/unavoid13_exactcheck.py docstring)\n")
        f.write(f"# generated by search/unavoid13_exactcert.py; root LP {val:.6f}; "
                f"{nb} branchings, {len(bb.out) - nb} leaves\n")
        f.write(f"family {os.path.basename(a.family)}\nfamily_sha256 {sha}\n")
        f.write(f"claim {a.claim}\nm {m}\nk {a.k}\n")
        for (x, y) in pts:
            f.write(f"cand {x} {y}\n")
        f.write("tree\n")
        for r in bb.out: f.write(r + "\n")
        f.write("end\n")
    print(f"  wrote {a.out}")


if __name__ == '__main__':
    main()
