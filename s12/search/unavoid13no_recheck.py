#!/usr/bin/env python3
"""unavoid13-no: theorem-grade exact re-solve of "h(F) >= k+1" for a family too large for
search/unavoid13_recheck.py's dominance pass.

    python3 search/unavoid13no_recheck.py FAMILY.txt --k 13 [--threads 4] [--out PREFIX]

Everything load-bearing is exact (fractions.Fraction):
  1. exact arrangement vertices (unavoid13_recheck.exact_vertices: corners + intersections of every
     pair of non-parallel closed edges, parallel = exact zero cross product);
  2. exact incidence of every vertex (unavoid13_recheck.exact_incidence: float pre-filter, every pair
     within 1e-6 of a boundary decided in Fractions);
  3. exact candidate set K: start from K0 = the exact rows equal to the incidence rows of the float
     maximal cells (unavoid13no_lib.maximal_cells -- a HINT only); then for EVERY exact vertex v
     check, on exact incidence rows, that some member of K contains v's row (posting list of v's
     rarest square, packed-bit subset test); if none does, v is added to K.  So K dominates every
     arrangement vertex, hence every point of the container, by the vertex-maximality argument
     (notes/unavoid13.md s1) -- independently of what the float filter did;
  4. the hitting-set IP on K (exact rows), NO symmetry cut, solved by highspy MIP (gap 0), by
     scipy.optimize.milp, and by the purpose-built B&B (unavoid13no_bb.BB, no D4 root rule);
     all three must report "no k-set" / "optimum >= k+1".
The instance (family, exact candidate coordinates, incidence) is dumped next to the family.
"""
import sys, os, time, json, argparse
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13_recheck as R, unavoid13no_lib as N


def exact_candidates(m, poses, threads=4, verbose=True):
    t0 = time.time()
    V = R.exact_vertices(m, poses, verbose=verbose)
    B = R.exact_incidence(m, poses, V, verbose=verbose)
    n = len(poses); K = len(V)
    packed = np.packbits(B, axis=1)
    # dedupe exact rows
    _, idx = np.unique(packed, axis=0, return_index=True); idx = np.sort(idx)
    nz = B[idx].any(1); idx = idx[nz]
    Bu = B[idx]; Pu = packed[idx]; Vf = np.array([(float(V[i][0]), float(V[i][1])) for i in idx])
    if verbose: print(f"  {len(idx)} distinct nonempty exact incidence rows", flush=True)
    # hint: float maximal cells
    sq, Vc, Bc = N.maximal_cells(m, poses, threads=threads, verbose=verbose)
    inK = np.zeros(len(idx), dtype=bool)
    # match float cells to exact rows by incidence set (the float cell's vertex may be another
    # vertex of the same incidence class than the one the dedupe kept); unmatched cells are ignored
    rowindex = {Pu[t].tobytes(): t for t in range(len(idx))}
    Pc = np.packbits(Bc, axis=1)
    for k in range(len(Bc)):
        t = rowindex.get(Pc[k].tobytes())
        if t is not None: inK[t] = True
    if verbose: print(f"  hint: {int(inK.sum())} exact rows matched to {len(Vc)} float cells", flush=True)
    # witness check for every row not in K; add undominated ones
    cnt = Bu.sum(1)
    colcnt = Bu.sum(0)
    added = 0
    changed = True
    while changed:
        changed = False
        Kidx = np.nonzero(inK)[0]
        post = [Kidx[Bu[Kidx, j]] for j in range(n)]
        PK = Pu[Kidx]; cntK = cnt[Kidx]
        pos_in_K = {int(k): i for i, k in enumerate(Kidx)}
        for t in np.nonzero(~inK)[0]:
            row = np.nonzero(Bu[t])[0]
            j = row[np.argmin(colcnt[row])]
            cand = post[j]
            cand = cand[cnt[cand] > cnt[t]]
            if len(cand) == 0:
                inK[t] = True; added += 1; changed = True; continue
            ci = np.array([pos_in_K[int(c)] for c in cand])
            sub = ((Pu[t][None, :] & ~PK[ci]) == 0).all(1)
            if not sub.any():
                inK[t] = True; added += 1; changed = True
        if verbose: print(f"  witness pass: |K| = {int(inK.sum())} (+{added} undominated rows added), {time.time()-t0:.0f}s", flush=True)
        added = 0
    Kidx = np.nonzero(inK)[0]
    # final exact dominance among K itself (posting lists on a small set)
    keep = L.drop_dominated_fast(Bu[Kidx])
    Kidx = Kidx[keep]
    Vex = [V[idx[t]] for t in Kidx]
    if verbose: print(f"  exact candidates: {len(Kidx)} ({time.time()-t0:.0f}s)", flush=True)
    return Vex, Bu[Kidx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family'); ap.add_argument('--k', type=int, default=13); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--out', default=None); ap.add_argument('--no-bb', action='store_true'); ap.add_argument('--no-scipy', action='store_true')
    ap.add_argument('--time', type=float, default=None)
    a = ap.parse_args()
    m, poses = L.read_family(a.family)
    out = a.out or os.path.splitext(a.family)[0] + '_exact'
    print(f"family {a.family}: m = {m}, {len(poses)} squares", flush=True)
    bad = [p for p in poses if not L.admissible_exact(m, *p)]
    assert not bad, f"inadmissible poses: {bad[:3]}"
    assert len(set(poses)) == len(poses), "duplicate poses"
    Vex, B = exact_candidates(m, poses, threads=a.threads)
    np.save(out + '_B.npy', np.packbits(B, axis=1))
    with open(out + '_cand.txt', 'w') as f:
        f.write("# candidate index; exact vertex (x, y) as Fractions; squares containing it (indices into the family)\n")
        for k in range(len(B)):
            f.write(f"{k} {Vex[k][0]} {Vex[k][1]} : " + " ".join(map(str, np.nonzero(B[k])[0])) + "\n")
    res = dict(family=a.family, m=str(m), n_squares=len(poses), n_cand=int(len(B)), k=a.k)
    t0 = time.time()
    r1 = N.solve_ip(B, a.k, threads=a.threads, feasibility=False, time_limit=a.time)
    print(f"  highspy MIP (min sum x, cutoff {a.k}+1/2): {r1['status']}, feasible(k)={r1['feasible']}, nodes {r1['nodes']}, {time.time()-t0:.0f}s", flush=True)
    res['highspy_feasible_k'] = r1['feasible']
    if not a.no_scipy:
        t0 = time.time()
        r2 = L.solve_hitting_set_scipy(B, time_limit=a.time)
        print(f"  scipy milp: status {r2.status} ({r2.message}), h = {r2.fun}, {time.time()-t0:.0f}s", flush=True)
        res['scipy_h'] = None if r2.fun is None else int(round(r2.fun))
    if not a.no_bb:
        import unavoid13no_bb as BBM
        t0 = time.time()
        bb = BBM.BB(B, a.k, threads=a.threads, perms=None, j0=None, verbose=False, time_limit=a.time)
        sol = bb.run()
        verdict = 'FEASIBLE' if sol is not None else ('INFEASIBLE' if (a.time is None or time.time() - t0 < a.time) else 'UNDECIDED')
        print(f"  B&B (no symmetry rule): {verdict}, nodes {bb.nodes}, LPs {bb.lps}, {time.time()-t0:.0f}s", flush=True)
        res['bb'] = verdict
    json.dump(res, open(out + '_result.json', 'w'), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
