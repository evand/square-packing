#!/usr/bin/env python3
"""unavoid13: independent, theorem-grade re-solve of a hitting-set instance from a dumped family.

    python3 search/unavoid13_recheck.py FAMILY.txt [--out PREFIX] [--threads 4] [--shrink]

Steps (everything load-bearing is exact):
  1. read the family F (rational poses), check admissibility exactly;
  2. enumerate the arrangement vertices EXACTLY (Fractions): all square corners and the
     intersection point of every pair of non-parallel closed edges of two different squares
     (no near-parallel cut-off: parallel edges are detected by an exact zero cross product);
  3. exact incidence of every vertex in every square (Fraction test |x'|, |y'| <= 1/2, with a
     float pre-filter that only decides pairs at distance > 1e-6 from the boundary);
  4. dedupe incidence sets, drop dominated ones (subsets), dump the incidence matrix;
  5. solve the hitting-set IP twice: highspy MIP (gap 0) and scipy.optimize.milp; and, when the
     instance is small enough, a third time with a hand-written branch-and-bound that uses only
     LP bounds from highspy's simplex (no MIP code path).
  6. --shrink: greedily drop squares while h stays >= target (float-free: exact incidence is
     already computed; the IP is re-run on the reduced row set), writing the reduced family.
"""
import sys, os, time, json, argparse, math, itertools
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L

HALF = F(1, 2)


def exact_vertices(m, poses, verbose=True):
    t0 = time.time()
    n = len(poses)
    corners = [L.square_corners_exact(*p) for p in poses]
    ctr = np.array([(float(cx), float(cy)) for (_, cx, cy) in poses])
    V = set()
    for i in range(n):
        for c in corners[i]: V.add(c)
    npairs = 0
    for i in range(n):
        Ci = corners[i]
        d2 = ((ctr[i + 1:] - ctr[i]) ** 2).sum(1)
        for jj in np.nonzero(d2 <= 2.0 + 1e-9)[0]:
            j = i + 1 + int(jj); Cj = corners[j]; npairs += 1
            for a in range(4):
                P = Ci[a]; R = (Ci[(a + 1) % 4][0] - P[0], Ci[(a + 1) % 4][1] - P[1])
                for b in range(4):
                    Q = Cj[b]; S = (Cj[(b + 1) % 4][0] - Q[0], Cj[(b + 1) % 4][1] - Q[1])
                    den = R[0] * S[1] - R[1] * S[0]
                    if den == 0: continue
                    PQ = (Q[0] - P[0], Q[1] - P[1])
                    t = (PQ[0] * S[1] - PQ[1] * S[0]) / den
                    if t < 0 or t > 1: continue
                    s = (PQ[0] * R[1] - PQ[1] * R[0]) / den
                    if s < 0 or s > 1: continue
                    V.add((P[0] + t * R[0], P[1] + t * R[1]))
    V = sorted(V)
    if verbose:
        print(f"  exact vertices: {len(V)} from {n} squares ({npairs} close pairs), {time.time()-t0:.1f}s", flush=True)
    return V


def exact_incidence(m, poses, V, verbose=True):
    """(K, n) bool, exact.  Float pre-filter decides pairs whose float depth is > 1e-6 from 0."""
    t0 = time.time()
    sq = L.Squares(poses)
    Vf = np.array([(float(x), float(y)) for x, y in V])
    n = len(poses); K = len(V)
    B = np.zeros((K, n), dtype=bool)
    nexact = 0
    trig = [L.trig(p[0]) for p in poses]
    for i0 in range(0, K, 2000):
        D = sq.depth(Vf[i0:i0 + 2000])
        B[i0:i0 + 2000] = D >= 0
        amb = np.nonzero(np.abs(D) <= 1e-6)
        for r, j in zip(*amb):
            k = i0 + int(r); j = int(j)
            u, cx, cy = poses[j]; c, s = trig[j]
            px, py = V[k]
            dx, dy = px - cx, py - cy
            x = dx * c + dy * s; y = -dx * s + dy * c
            B[k, j] = (-HALF <= x <= HALF) and (-HALF <= y <= HALF)
            nexact += 1
    if verbose:
        print(f"  exact incidence: {K} x {n}, {nexact} borderline pairs decided in Fractions, "
              f"{B.sum()} incidences, {time.time()-t0:.1f}s", flush=True)
    return B


def reduce_candidates(B, dominance_max=20000000):
    # dedupe on the packed rows (8x less data to sort than the bool rows)
    packed = np.packbits(B, axis=1)
    _, idx = np.unique(packed, axis=0, return_index=True)
    idx = np.sort(idx)
    uniq = B[idx]
    nz = uniq.any(1); uniq = uniq[nz]; idx = idx[nz]
    if len(uniq) <= dominance_max:
        keep = L.drop_dominated_fast(uniq)
        uniq = uniq[keep]; idx = idx[keep]
    return uniq, idx


def bb_lp(B, threads=2, verbose=True, node_limit=10 ** 6):
    """Hand branch-and-bound for min sum x, B^T x >= 1, x in {0,1}: LP bound by highspy simplex
    (LP only), branching on the most fractional variable; returns the optimum."""
    import highspy
    K, n = B.shape
    starts = np.zeros(K + 1, dtype=np.int64); starts[1:] = np.cumsum(B.sum(1))
    index = np.concatenate([np.nonzero(B[k])[0] for k in range(K)]).astype(np.int32)
    value = np.ones(len(index))
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', threads)
    lp = highspy.HighsLp(); lp.num_col_ = K; lp.num_row_ = n
    lp.col_cost_ = np.ones(K); lp.col_lower_ = np.zeros(K); lp.col_upper_ = np.ones(K)
    lp.row_lower_ = np.ones(n); lp.row_upper_ = np.full(n, highspy.kHighsInf)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = value
    h.passModel(lp)
    best = [K + 1]; nodes = [0]
    # greedy incumbent
    cov = np.zeros(n, dtype=bool); x0 = np.zeros(K)
    while not cov.all():
        gain = (B & ~cov[None, :]).sum(1); k = int(gain.argmax()); x0[k] = 1; cov |= B[k]
    best[0] = int(x0.sum())

    def node(fix0, fix1):
        nodes[0] += 1
        if nodes[0] > node_limit: raise RuntimeError('node limit')
        lo = np.zeros(K); hi = np.ones(K)
        for k in fix0: hi[k] = 0
        for k in fix1: lo[k] = 1
        h.changeColsBounds(K, np.arange(K, dtype=np.int32), lo, hi)
        h.run()
        st = h.getModelStatus()
        if st != highspy.HighsModelStatus.kOptimal:
            return
        val = h.getInfo().objective_function_value
        if math.ceil(val - 1e-9) >= best[0]:
            return
        x = np.array(h.getSolution().col_value)
        frac = np.abs(x - np.round(x))
        k = int(frac.argmax())
        if frac[k] < 1e-9:
            best[0] = min(best[0], int(round(x.sum())))
            return
        node(fix0, fix1 + [k])
        node(fix0 + [k], fix1)

    node([], [])
    if verbose: print(f"  hand B&B (LP bounds only): optimum {best[0]}, {nodes[0]} nodes", flush=True)
    return best[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family')
    ap.add_argument('--out', default=None)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--shrink', type=int, default=None, help='target h; greedily drop squares keeping h >= target')
    ap.add_argument('--no-bb', action='store_true')
    ap.add_argument('--bb-max', type=int, default=20000, help='hand B&B only if candidates <= this')
    a = ap.parse_args()
    m, poses = L.read_family(a.family)
    out = a.out or os.path.splitext(a.family)[0] + '_recheck'
    print(f"family {a.family}: m = {m}, {len(poses)} squares")
    bad = [p for p in poses if not L.admissible_exact(m, *p)]
    assert not bad, f"inadmissible poses: {bad[:3]}"
    assert len(set(poses)) == len(poses), "duplicate poses"
    V = exact_vertices(m, poses)
    B = exact_incidence(m, poses, V)
    Bc, idx = reduce_candidates(B)
    print(f"  candidates after dedupe/dominance: {len(Bc)}")
    np.save(out + '_B.npy', np.packbits(Bc, axis=1))
    with open(out + '_cand.txt', 'w') as f:
        f.write("# candidate index; exact vertex (x, y) as Fractions; squares containing it (indices into the family)\n")
        for k in range(len(Bc)):
            x, y = V[idx[k]]
            f.write(f"{k} {x} {y} : " + " ".join(map(str, np.nonzero(Bc[k])[0])) + "\n")
    t0 = time.time()
    r1 = L.solve_hitting_set(Bc, threads=a.threads)
    print(f"  highspy MIP: {r1['status']}, h = {r1['obj']}, LP = {r1['lp']:.6f}, {time.time()-t0:.1f}s")
    t0 = time.time()
    r2 = L.solve_hitting_set_scipy(Bc)
    print(f"  scipy milp: status {r2.status} ({r2.message}), h = {r2.fun}, {time.time()-t0:.1f}s")
    r3 = None
    if not a.no_bb and len(Bc) <= a.bb_max:
        t0 = time.time()
        try:
            r3 = bb_lp(Bc, threads=a.threads)
        except RuntimeError as e:
            print("  hand B&B:", e)
        print(f"  ({time.time()-t0:.1f}s)")
    res = dict(family=a.family, m=str(m), n_squares=len(poses), n_vertices=len(V), n_cand=int(len(Bc)),
               highspy=r1['obj'], highspy_lp=r1['lp'], scipy=(None if r2.fun is None else int(round(r2.fun))),
               bb=r3)
    if a.shrink is not None and r1['obj'] is not None and r1['obj'] >= a.shrink:
        # greedy shrink: try dropping squares one at a time (in a fixed order), keep if h >= target
        # pass 1: drop whole D4 orbits (keeps F symmetric, 8x fewer IPs); pass 2: single squares
        pose_index = {p: k for k, p in enumerate(poses)}
        seen = set(); orbits = []
        for k, p in enumerate(poses):
            if k in seen: continue
            orb = sorted(pose_index[q] for q in L.d4_images(m, *p) if q in pose_index)
            seen.update(orb); orbits.append(orb)
        keep = set(range(len(poses)))
        for pas, groups in ((1, orbits[::-1]), (2, [[k] for k in sorted(keep, reverse=True)])):
            for grp in groups:
                if not set(grp) <= keep: continue
                trial = sorted(keep - set(grp))
                # dropping squares only makes more candidates dominated, so the once-reduced Bc
                # restricted to the remaining squares is still an exact candidate set
                Bt2 = np.unique(Bc[:, trial], axis=0); Bt2 = Bt2[Bt2.any(1)]
                r = L.solve_hitting_set(Bt2, threads=a.threads)
                if r['obj'] is not None and r['obj'] >= a.shrink:
                    keep = set(trial)
                    print(f"    pass {pas}: dropped {grp}: |F| = {len(keep)}, h = {r['obj']}", flush=True)
            if pas == 1:
                print(f"  after orbit pass: |F| = {len(keep)}", flush=True)
        keep = sorted(keep)
        fam = [poses[k] for k in keep]
        L.write_family(out + '_shrunk_family.txt', m, fam, header=f"shrunk: h >= {a.shrink} with {len(fam)} squares")
        res['shrunk_n'] = len(fam)
        print(f"  shrunk family: {len(fam)} squares -> {out}_shrunk_family.txt")
    with open(out + '_result.json', 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
