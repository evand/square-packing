#!/usr/bin/env python3
"""unavoid13: exact (Fraction) re-check of a symmetric-IP infeasibility claim.

    python3 search/unavoid13_symcheck.py FAMILY.txt --group D4 --k 13

Claim checked: no point set invariant under the group G hits every square of the family F with at
most k points.  Candidates: the exact arrangement vertices of F (unavoid13_recheck.exact_vertices)
plus the exact intersections of every square edge with every symmetry axis of G and the centre
(a point on an axis is dominated, within the axis, by such an intersection); orbit candidates with
exact orbit sizes; one row per G-orbit of F; feasibility IP  sum |O| x_O <= k  solved by highspy and
by scipy.optimize.milp.  Infeasible in both => the claim is a theorem for this F (and hence for
every G-symmetric k-point set, since F is a finite subfamily of all admissible squares).
"""
import sys, os, time, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L
import unavoid13_recheck as R
import unavoid13_sym as S


def exact_axis_points(m, poses, group):
    m = F(m); c = m / 2
    lines = []
    if group in ('D4', 'V', 'C4', 'C2'):
        lines += [((c, F(0)), (F(0), F(1))), ((F(0), c), (F(1), F(0)))]
    if group == 'Rx':
        lines += [((c, F(0)), (F(0), F(1)))]
    if group in ('D4', 'Rd'):
        lines += [((F(0), F(0)), (F(1), F(1)))]
    if group == 'D4':
        lines += [((F(0), m), (F(1), F(-1)))]
    pts = {(c, c)}
    for p in poses:
        C = L.square_corners_exact(*p)
        for a in range(4):
            A = C[a]; Rr = (C[(a + 1) % 4][0] - A[0], C[(a + 1) % 4][1] - A[1])
            for (P0, d) in lines:
                den = Rr[0] * d[1] - Rr[1] * d[0]
                if den == 0: continue
                t = ((P0[0] - A[0]) * d[1] - (P0[1] - A[1]) * d[0]) / den
                if 0 <= t <= 1:
                    pts.add((A[0] + t * Rr[0], A[1] + t * Rr[1]))
    return sorted(pts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family'); ap.add_argument('--group', required=True); ap.add_argument('--k', type=int, default=13)
    ap.add_argument('--threads', type=int, default=2)
    a = ap.parse_args()
    m, poses = L.read_family(a.family)
    G = S.group_maps(m, a.group)
    pindex = {p: i for i, p in enumerate(poses)}
    perms = []
    for (pt, ps) in G:
        img = [ps(*p) for p in poses]
        assert all(q in pindex for q in img), "family is not closed under the group"
        perms.append(np.array([pindex[q] for q in img]))
    V = R.exact_vertices(m, poses)
    Vax = exact_axis_points(m, poses, a.group)
    V = sorted(set(V) | set(Vax))
    print(f"  {len(Vax)} exact axis points; {len(V)} candidate points in total")
    B = R.exact_incidence(m, poses, V)
    # orbit incidence and exact orbit sizes
    OB = np.zeros_like(B)
    for perm in perms:
        OB[:, perm] |= B
    sizes = np.array([len({pt(x, y) for (pt, _) in G}) for (x, y) in V], dtype=np.int64)
    # dedupe orbit sets keeping the smallest orbit
    order = np.argsort(sizes, kind='stable'); keyed = {}
    P = np.packbits(OB, axis=1)
    for k in order:
        key = P[k].tobytes()
        if key not in keyed: keyed[key] = k
    cand = np.array(sorted(keyed.values())); OB = OB[cand]; sizes = sizes[cand]
    seen = set(); rows = []
    for i in range(len(poses)):
        if i in seen: continue
        orb = {int(perm[i]) for perm in perms}; seen |= orb; rows.append(i)
    Bs = OB[:, np.array(rows)]
    nz = Bs.any(1); Bs = Bs[nz]; sizes = sizes[nz]
    print(f"  {len(Bs)} orbit candidates (sizes: " + ", ".join(f"{s}:{int((sizes==s).sum())}" for s in sorted(set(sizes.tolist()))) + f"), {Bs.shape[1]} square orbits")
    # highspy feasibility
    import highspy
    Kc, n = Bs.shape
    starts = np.zeros(Kc + 1, dtype=np.int64); starts[1:] = np.cumsum(Bs.sum(1) + 1)
    index = np.zeros(starts[-1], dtype=np.int32); value = np.zeros(starts[-1])
    for k in range(Kc):
        nzk = np.nonzero(Bs[k])[0]
        index[starts[k]:starts[k + 1] - 1] = nzk; value[starts[k]:starts[k + 1] - 1] = 1.0
        index[starts[k + 1] - 1] = n; value[starts[k + 1] - 1] = float(sizes[k])
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', a.threads)
    h.setOptionValue('mip_feasibility_tolerance', 1e-9)
    lp = highspy.HighsLp(); lp.num_col_ = Kc; lp.num_row_ = n + 1
    lp.col_cost_ = sizes.astype(float); lp.col_lower_ = np.zeros(Kc); lp.col_upper_ = np.ones(Kc)
    lp.row_lower_ = np.concatenate([np.ones(n), [-highspy.kHighsInf]])
    lp.row_upper_ = np.concatenate([np.full(n, highspy.kHighsInf), [float(a.k)]])
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = value
    lp.integrality_ = [highspy.HighsVarType.kInteger] * Kc
    h.passModel(lp); h.run(); st1 = h.modelStatusToString(h.getModelStatus())
    # also the minimum symmetric hitting number (optimisation form)
    lp.row_upper_[n] = highspy.kHighsInf; h2 = highspy.Highs(); h2.setOptionValue('output_flag', False)
    h2.setOptionValue('threads', a.threads); h2.setOptionValue('mip_rel_gap', 0.0); h2.passModel(lp); h2.run()
    hmin = h2.getInfo().objective_function_value
    # scipy
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import csc_matrix
    A = csc_matrix(np.vstack([Bs.T.astype(float), sizes[None, :].astype(float)]))
    lb = np.concatenate([np.ones(n), [-np.inf]]); ub = np.concatenate([np.full(n, np.inf), [float(a.k)]])
    r2 = milp(np.zeros(Kc), constraints=LinearConstraint(A, lb, ub), integrality=np.ones(Kc), bounds=Bounds(0, 1))
    print(f"  highspy feasibility (k = {a.k}): {st1};  minimum {a.group}-symmetric hitting number of F: {hmin:.0f};  "
          f"scipy milp: status {r2.status} ({r2.message})")
    verdict = (st1 == 'Infeasible') and (r2.status == 2)
    print(f"VERDICT: {'no ' + a.group + '-symmetric ' + str(a.k) + '-point set hits F -- PROVED (exact candidates, two solvers)' if verdict else 'NOT confirmed'}")
    json.dump(dict(family=a.family, group=a.group, k=a.k, n_squares=len(poses), n_cand=int(Kc), highspy=st1,
                   sym_min=hmin, scipy_status=int(r2.status), proved=bool(verdict)),
              open(os.path.splitext(a.family)[0] + f'_symcheck_{a.group}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
