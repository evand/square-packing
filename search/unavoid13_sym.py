#!/usr/bin/env python3
"""unavoid13: symmetry-reduced loop (the brief's structured family (ii)) -- an UPPER-BOUND search.

    python3 search/unavoid13_sym.py --m 4 --k 13 --group D4|C4|V --tag sym_d4 [--init ...] [--resume F]

A point set invariant under a group G <= D4 (of the container) hits a G-symmetric family F iff it
hits one square per G-orbit of F.  Variables: G-orbits of candidate vertices (cost = orbit size);
rows: one per G-orbit of squares; feasibility IP  sum_O |O| x_O <= k.  Infeasible means only that
no G-SYMMETRIC k-point set hits F (not a lower bound on the general problem).  Feasible: polish,
symmetrise, search violations, add them (D4-closed), repeat; when no float violation is found the
set is written as a certificate and checked exactly with unavoid13_check.py --tri --seg.
Dominance pruning is NOT applied (a vertex on a symmetry axis has a smaller, cheaper orbit than a
dominating vertex off it), only deduplication of incidence sets.
"""
import sys, os, math, time, json, argparse, subprocess
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L
import unavoid13_loop as LOOP


def group_maps(m, name):
    """List of (point_map, pose_map) for the group elements."""
    m = F(m)
    def rot(k):
        def pt(x, y):
            for _ in range(k): x, y = m - y, x
            return x, y
        def ps(u, cx, cy):
            cx, cy = pt(cx, cy); return L.pose_key(u, cx, cy)
        return pt, ps
    def refl_x():
        return (lambda x, y: (m - x, y)), (lambda u, cx, cy: L.pose_key(L.norm_u(-u), m - cx, cy))
    def refl_y():
        return (lambda x, y: (x, m - y)), (lambda u, cx, cy: L.pose_key(L.norm_u(-u), cx, m - cy))
    def compose(a, b):   # a after b
        return (lambda x, y: a[0](*b[0](x, y))), (lambda u, cx, cy: a[1](*b[1](u, cx, cy)))
    if name == 'C4': return [rot(k) for k in range(4)]
    if name == 'D4': return [rot(k) for k in range(4)] + [compose(rot(k), refl_x()) for k in range(4)]
    if name == 'V':  return [rot(0), refl_x(), refl_y(), compose(refl_x(), refl_y())]
    if name == 'C2': return [rot(0), rot(2)]
    if name == 'Rx': return [rot(0), refl_x()]                       # one axis reflection
    if name == 'Rd':                                                 # one diagonal reflection (x,y) -> (y,x)
        return [rot(0), ((lambda x, y: (y, x)), (lambda u, cx, cy: L.pose_key(L.norm_u(-u), cy, cx)))]
    raise ValueError(name)


def axis_points(sq, mf, group):
    """Intersections of all square edges with the symmetry axes of the group, plus the centre."""
    c = mf / 2
    lines = []   # (point, direction)
    if group in ('D4', 'V', 'C4', 'C2'):
        lines += [((c, 0.0), (0.0, 1.0)), ((0.0, c), (1.0, 0.0))]
    if group == 'Rx':
        lines += [((c, 0.0), (0.0, 1.0))]
    if group in ('D4', 'Rd'):
        lines += [((0.0, 0.0), (1.0, 1.0))]
    if group == 'D4':
        lines += [((0.0, mf), (1.0, -1.0))]
    pts = [np.array([[c, c]])]
    C = sq.corners
    for (P0, d) in lines:
        P0 = np.array(P0); d = np.array(d)
        for a in range(4):
            A = C[:, a]; R = C[:, (a + 1) % 4] - A
            den = R[:, 0] * d[1] - R[:, 1] * d[0]
            ok = np.abs(den) > 1e-12
            t = np.where(ok, ((P0 - A)[:, 0] * d[1] - (P0 - A)[:, 1] * d[0]) / np.where(ok, den, 1.0), -1.0)
            ok &= (t >= -1e-9) & (t <= 1 + 1e-9)
            if ok.any():
                pts.append(A[ok] + t[ok, None] * R[ok])
    out = np.concatenate(pts, 0)
    return np.clip(out, 0.0, mf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=str, default='4')
    ap.add_argument('--k', type=int, default=13)
    ap.add_argument('--group', type=str, default='D4')
    ap.add_argument('--tag', type=str, required=True)
    ap.add_argument('--init', type=str, default='tiles,grid2,support:40')
    ap.add_argument('--support', type=str, default='search/cover4_exact_support.txt')
    ap.add_argument('--resume', type=str, default=None)
    ap.add_argument('--rounds', type=int, default=80)
    ap.add_argument('--threads', type=int, default=3)
    ap.add_argument('--add', type=int, default=10)
    ap.add_argument('--pitch', type=float, default=0.01)
    ap.add_argument('--dth', type=float, default=0.5)
    ap.add_argument('--pden', type=int, default=1000)
    ap.add_argument('--zm-depth', type=int, default=14)
    ap.add_argument('--ip-time', type=float, default=None)
    a = ap.parse_args()
    m = F(a.m); mf = float(m)
    G = group_maps(m, a.group)
    out = f"runs/unavoid13_{a.tag}"; os.makedirs(out, exist_ok=True)
    Fam = L.read_family(a.resume)[1] if a.resume else LOOP.initial_family(m, a.init, a.support)
    Fset = set(Fam)
    log = open(f"{out}/round_log.txt", 'a')
    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()
    say(f"# unavoid13 symmetric loop: m = {m}, k = {a.k}, group {a.group}, |F0| = {len(Fam)}, {time.ctime()}")
    for rnd in range(a.rounds):
        t0 = time.time()
        say(f"--- round {rnd}: |F| = {len(Fam)}")
        # candidates (deduped, no dominance)
        sq = L.Squares(Fam)
        V = np.clip(L.arrangement_vertices(sq), 0.0, mf)
        # symmetric cells: a symmetric set may have points ON the axes of G (stabiliser > 1), whose
        # dominating cell must itself lie on the axis: add every intersection of a square edge with
        # each axis line, and the centre.  (For a point p on axis l, P(p) ∩ l is a segment whose
        # endpoints are such intersections and dominate p.)
        V = np.concatenate([V, axis_points(sq, mf, a.group)], 0)
        packed = L.incidence_packed(sq, V, tol=1e-7)
        uniq, idx = L.dedupe_rows(packed)
        B = L.unpack_rows(uniq, len(sq)); nz = B.any(1); B = B[nz]; idx = idx[nz]; rep = V[idx]
        # permutations of F under G
        pindex = {p: i for i, p in enumerate(Fam)}
        perms = []
        for (pt, ps) in G:
            perm = np.array([pindex[ps(*p)] for p in Fam])
            perms.append(perm)
        # orbit incidence: union over g of g.S(v)  (g.S(v) = {g.Q : Q in S(v)})
        K = len(B)
        OB = np.zeros_like(B)
        for perm in perms:
            OB[:, perm] |= B          # column Q of B maps to column perm[Q]
        # orbit sizes from the representative coordinates
        sizes = np.zeros(K, dtype=np.int64)
        for k in range(K):
            x, y = rep[k]
            pts = {(round(float(pt(F(x), F(y))[0]), 6), round(float(pt(F(x), F(y))[1]), 6)) for (pt, _) in G}
            sizes[k] = len(pts)     # 6-digit rounding: float vertices carry ~1e-8 noise
        # dedupe orbit incidence sets keeping the cheapest orbit
        order = np.lexsort((sizes, ))
        keyed = {}
        for k in order:
            key = np.packbits(OB[k]).tobytes()
            if key not in keyed: keyed[key] = k
        cand = np.array(sorted(keyed.values()))
        OB = OB[cand]; sizes = sizes[cand]; rep_c = rep[cand]
        # one row per square orbit
        seen = set(); rows = []
        for i in range(len(Fam)):
            if i in seen: continue
            orb = {int(perm[i]) for perm in perms}; seen |= orb; rows.append(i)
        rows = np.array(rows)
        Bs = OB[:, rows]
        say(f"  candidates: {len(V)} vertices, {K} distinct, {len(cand)} orbit candidates; "
            f"{len(rows)} square orbits; {time.time()-t0:.1f}s")
        # feasibility IP with orbit costs
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
        if a.ip_time: h.setOptionValue('time_limit', float(a.ip_time))
        lp = highspy.HighsLp(); lp.num_col_ = Kc; lp.num_row_ = n + 1
        lp.col_cost_ = sizes.astype(float); lp.col_lower_ = np.zeros(Kc); lp.col_upper_ = np.ones(Kc)
        lp.row_lower_ = np.concatenate([np.ones(n), [-highspy.kHighsInf]])
        lp.row_upper_ = np.concatenate([np.full(n, highspy.kHighsInf), [float(a.k)]])
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = value
        lp.integrality_ = [highspy.HighsVarType.kInteger] * Kc
        h.passModel(lp); h.run()
        st = h.modelStatusToString(h.getModelStatus())
        L.write_family(f"{out}/F_round{rnd:02d}.txt", m, Fam, header=f"round {rnd}; symmetric IP {st}")
        if st != 'Optimal':
            say(f"  symmetric IP: {st} -> no {a.group}-symmetric {a.k}-point set hits F ({time.time()-t0:.1f}s)")
            break
        x = np.array(h.getSolution().col_value)
        chosen = np.nonzero(x > 0.5)[0]
        obj = int(round(h.getInfo().objective_function_value))
        say(f"  symmetric IP: feasible with {obj} points ({len(chosen)} orbits), {time.time()-t0:.1f}s")
        # the full point set
        P = []
        for k in chosen:
            x0, y0 = F(rep_c[k][0]), F(rep_c[k][1])
            pts = {(round(float(pt(x0, y0)[0]), 6), round(float(pt(x0, y0)[1]), 6)) for (pt, _) in G}
            P += sorted(pts)
        P = np.array(P, float)
        say(f"  points (raw): " + " ".join(f"({x:.4f},{y:.4f})" for x, y in P))
        Pp, tstar = L.polish_positions(P, m, [sq], rounds=5, verbose=False)
        # symmetrise: r0 = (1/|G|) sum_g g^{-1}(polished image g(x0,y0)); averaging over all group
        # elements (not one per distinct image) keeps an on-axis orbit on its axis.
        Ps = []
        pos = 0
        for k in chosen:
            x0, y0 = F(rep_c[k][0]), F(rep_c[k][1])
            imgs = sorted({(round(float(pt(x0, y0)[0]), 6), round(float(pt(x0, y0)[1]), 6)) for (pt, _) in G})
            nimg = len(imgs)
            pol = {imgs[i]: Pp[pos + i] for i in range(nimg)}
            acc = np.zeros(2)
            for (pt, _) in G:
                gx, gy = pt(x0, y0)
                key = (round(float(gx), 6), round(float(gy), 6))
                q = pol[key]
                # inverse of g as a group element: g'(g(z)) = z for a generic test point z
                tx, ty = pt(F(3, 10), F(7, 10))
                for (pt2, _) in G:
                    bx, by = pt2(tx, ty)
                    if bx == F(3, 10) and by == F(7, 10):
                        qx, qy = pt2(F(float(q[0])), F(float(q[1])))
                        acc += (float(qx), float(qy)); break
            r0 = acc / len(G)
            imgs2 = sorted({(round(float(pt(F(r0[0]), F(r0[1]))[0]), 6), round(float(pt(F(r0[0]), F(r0[1]))[1]), 6)) for (pt, _) in G})
            Ps += imgs2; pos += nimg
        P = np.array(Ps, float)
        say(f"  polished (t* = {tstar:+.2e}) & symmetrised: {len(P)} points: " + " ".join(f"({x:.4f},{y:.4f})" for x, y in P))
        np.save(f"{out}/P_round{rnd:02d}.npy", P)
        viol, gmin = L.find_violations(P, m, dth_deg=a.dth, pitch=a.pitch, top=400, want=a.add)
        if not viol:
            Pex = L.snap_points(P, a.pden)
            cert = f"{out}/cand_round{rnd:02d}.txt"; L.write_cert(cert, m, Pex)
            say(f"  no float violation (grid min {gmin:+.2e}); certifying {cert}")
            ok, txt = LOOP.run_zeromargin(cert, m, f"{out}/zm_round{rnd:02d}.log", full=True, tri=True,
                                          depth=a.zm_depth, nproc=a.threads, oracle=f"{out}/oracle_round{rnd:02d}.txt")
            tail = [l for l in txt.strip().split('\n') if l.strip()][-4:]
            say("  checker: " + " | ".join(tail))
            if ok:
                say(f"  CERTIFIED: {len(Pex)}-point unavoidable set, {cert}"); break
            seeds = LOOP.read_oracle(f"{out}/oracle_round{rnd:02d}.txt")
            Pf = np.array([(float(x), float(y)) for x, y in Pex])
            viol, _ = L.find_violations(Pf, m, dth_deg=a.dth, pitch=a.pitch, top=50, want=a.add,
                                        seeds=seeds[:5000], cutoff=-1e-13)
            if not viol:
                say("  no violation found from the oracle seeds either: STALL"); break
        added = 0
        Pex = L.snap_points(P, 10 ** 6)
        for (f, cx, cy, th) in viol:
            key = L.rationalise_pose(m, cx, cy, th, den=1000)
            adm, hits, worst = L.exact_violation(m, Pex, *key)
            fresh = [kk for kk in L.d4_images(m, *key) if kk not in Fset]
            say(f"  violation f = {f:+.3e} at ({cx:.5f},{cy:.5f},{math.degrees(th):.3f}deg); exact adm {adm}, hits {hits}; +{len(fresh)}")
            for kk in fresh:
                Fset.add(kk); Fam.append(kk); added += 1
        say(f"  added {added}; round time {time.time()-t0:.1f}s")
        if added == 0:
            say("  nothing new: STALL"); break
    say(f"# done {time.ctime()}")


if __name__ == '__main__':
    main()
