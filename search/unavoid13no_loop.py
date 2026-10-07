#!/usr/bin/env python3
"""unavoid13-no: the family loop on exact maximal cells (C helper), aiming at h(F) >= 14.

    python3 search/unavoid13no_loop.py --tag A --init support,t4L [--k 13] [--add 40] [--ip-time 3600]

Round: cells of F (unavoid13no_lib.maximal_cells, lenient float) -> is there a k-column hitting set?
  1. repair heuristic from the previous k-set (1- and 2-swaps restricted to cells covering the
     missed rows) -- cheap, finds a k-set when one is near;
  2. else HiGHS feasibility IP (root D4 cut, optional), time limit; INFEASIBLE -> dump, stop.
Then polish the k-set (max-margin LP over F), float violation search (many distinct local minima),
rationalise, add poses + D4 images.  Logs: runs/unavoid13no_<tag>/round_log.txt.
"""
import sys, os, math, time, json, argparse
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L
import unavoid13_loop as LOOP
import unavoid13no_lib as N


def repair(B, P_idx, k, max_pairs=4000000, verbose=False):
    """B (K, n) bool.  P_idx: indices of a previous k-set (may miss rows).  Try to find k columns
    covering all rows by replacing 1 or 2 members with columns that cover the missed rows.
    Returns index array or None."""
    K, n = B.shape
    P = [int(i) for i in P_idx if 0 <= i < K]
    if len(P) < k:
        P = P + [0] * (k - len(P))
    cov = B[P].sum(0)
    missed = np.nonzero(cov == 0)[0]
    if len(missed) == 0: return np.array(P)
    Pb = np.packbits(B, axis=1)
    def covered(idx):
        return B[list(idx)].any(0).all()
    # 1-swap: remove i, add c where c covers all missed rows and the rows only i covered
    cand_all = np.nonzero(B[:, missed].all(1))[0]
    for pos in range(k):
        rest = P[:pos] + P[pos + 1:]
        cov_r = B[rest].any(0)
        need = np.nonzero(~cov_r)[0]
        cs = np.nonzero(B[:, need].all(1))[0]
        if len(cs):
            return np.array(rest + [int(cs[0])])
    # 2-swap: remove i, j; add c1, c2
    for pi in range(k):
        for pj in range(pi + 1, k):
            rest = [P[t] for t in range(k) if t != pi and t != pj]
            cov_r = B[rest].any(0)
            need = np.nonzero(~cov_r)[0]
            if len(need) == 0: return np.array(rest + [0, 0])
            # c1 must cover some of need, c2 the rest: pick c1 among columns covering need[0]
            c1s = np.nonzero(B[:, need[0]])[0]
            Bn = B[:, need]
            for c1 in c1s[:400]:
                left = need[~Bn[c1]]
                if len(left) == 0: return np.array(rest + [int(c1), int(c1)])
                c2 = np.nonzero(B[:, left].all(1))[0]
                if len(c2): return np.array(rest + [int(c1), int(c2[0])])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=str, default='4')
    ap.add_argument('--k', type=int, default=13)
    ap.add_argument('--tag', type=str, required=True)
    ap.add_argument('--init', type=str, default='support,t4L')
    ap.add_argument('--resume', type=str, default=None)
    ap.add_argument('--rounds', type=int, default=500)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--add', type=int, default=40, help='violated poses added per round (before D4)')
    ap.add_argument('--pitch', type=float, default=0.01)
    ap.add_argument('--dth', type=float, default=0.5)
    ap.add_argument('--ip-time', type=float, default=3600)
    ap.add_argument('--no-cut', action='store_true')
    ap.add_argument('--no-repair', action='store_true')
    ap.add_argument('--feas', action='store_true', help='pure feasibility IP (else min sum x with cutoff)')
    ap.add_argument('--tol', type=float, default=1e-11)
    ap.add_argument('--min-sep', type=float, default=0.05)
    ap.add_argument('--no-bb', action='store_true')
    ap.add_argument('--race', type=int, default=1, help='portfolio of this many 1-thread HiGHS runs; first verdict wins')
    ap.add_argument('--min-viol', type=float, default=0.003, help='add only poses violated by at least this depth (always the worst 3)')
    ap.add_argument('--bb-nodes', type=int, default=3000)
    ap.add_argument('--bb-time', type=float, default=600)
    ap.add_argument('--snap', type=int, default=0, help='snap violated poses to the family angles and the 1/snap centre grid when still violated')
    ap.add_argument('--angles', type=str, default='', help='comma list of extra rational u values allowed for snapping')
    a = ap.parse_args()
    m = F(a.m); mf = float(m)
    out = f"runs/unavoid13no_{a.tag}"; os.makedirs(out, exist_ok=True)
    log = open(f"{out}/round_log.txt", 'a')
    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()
    if a.resume:
        Fam = L.read_family(a.resume)[1]
    else:
        parts = []
        for item in a.init.split(','):
            if item == 'support': parts.append(N.support_poses(m))
            elif item == 't4L': parts.append(L.read_family('runs/unavoid13_t4L/F_round01.txt')[1])
            elif item.startswith('file:'): parts.append(L.read_family(item[5:])[1])
            else: parts.append(LOOP.initial_family(m, item, 'search/cover4_exact_support.txt'))
        Fam = N.family_union(m, [LOOP.initial_family(m, 'tiles,grid2')] + parts)
        Fam = L.d4_closure(m, Fam)
    Fset = set(Fam)
    say(f"# unavoid13no loop: m = {m}, k = {a.k}, |F0| = {len(Fam)}, {time.ctime()}")
    prevP = None
    for rnd in range(a.rounds):
        t0 = time.time()
        say(f"--- round {rnd}: |F| = {len(Fam)}")
        L.write_family(f"{out}/F_round{rnd:02d}.txt", m, Fam, header=f"round {rnd}")
        sq, V, B = N.maximal_cells(m, Fam, tol=a.tol, threads=a.threads, verbose=False)
        say(f"  cells: {B.shape[0]} candidates x {B.shape[1]} rows, {time.time()-t0:.0f}s")
        sol = None; how = ''
        if prevP is not None and not a.no_repair:
            # map previous points to cells: nearest cell containing... use cells whose vertex is closest
            t1 = time.time()
            idx = []
            for (x, y) in prevP:
                d = (V[:, 0] - x) ** 2 + (V[:, 1] - y) ** 2
                idx.append(int(d.argmin()))
            r = repair(B, idx, a.k)
            if r is not None:
                sol = r; how = f'repair ({time.time()-t1:.0f}s)'
        if sol is None and not a.no_bb:
            # B&B dive first (cheap when a k-set exists)
            import unavoid13no_bb as BBM
            t1 = time.time()
            try:
                Vc, Bc, perms = N.d4_close_cells(m, Fam, V, B, verbose=False)
                j0 = Fam.index(L.pose_key(F(0), m / 2, m / 2))
            except AssertionError as e:
                Vc, Bc, perms, j0 = V, B, None, None
            bb = BBM.BB(Bc, a.k, threads=a.threads, perms=perms, j0=j0, verbose=False, node_limit=a.bb_nodes, time_limit=a.bb_time)
            r = bb.run()
            verdict = 'FEASIBLE' if r is not None else ('INFEASIBLE' if bb.nodes <= a.bb_nodes and time.time() - t1 < a.bb_time else 'UNDECIDED')
            say(f"  B&B: {verdict}, nodes {bb.nodes}, LPs {bb.lps}, pruned lp {bb.pruned_lp} pack {bb.pruned_pack}, {time.time()-t1:.0f}s")
            if r is not None:
                V, B = Vc, Bc; sol = np.array(r); how = 'B&B dive'
            elif verdict == 'INFEASIBLE':
                lpv, _ = N.lp_value(Bc, threads=a.threads)
                say(f"  INFEASIBLE (B&B): no {a.k}-column hitting set -> h(F) >= {a.k + 1}.  Dumping; confirming with highspy.")
                L.write_family(f"{out}/final_family.txt", m, Fam, header=f"h(F) >= {a.k + 1} (round {rnd})")
                L.dump_instance(f"{out}/final", m, Fam, Bc, Vc, extra=dict(k=a.k, round=rnd, lp=lpv))
                rr = N.solve_ip(Bc, a.k, threads=a.threads, feasibility=False, verbose=False)
                say(f"  highspy on the same instance: {rr['status']} feasible={rr['feasible']} {rr['time']:.0f}s")
                if rr['feasible'] is False:
                    say(f"# done {time.ctime()}"); return
                sol = np.nonzero(rr['x'] > 0.5)[0]; how = 'highspy (B&B disagreed!)'
                V, B = Vc, Bc
        if sol is None:
            extra = []
            if not a.no_cut:
                j0 = Fam.index(L.pose_key(F(0), m / 2, m / 2))
                # D4 on cells: match by mapped coordinates
                try:
                    extra = [N.root_cut_cells(m, V, B, j0, poses=Fam)]
                except AssertionError as e:
                    say(f"  root cut unavailable: {e}")
            t1 = time.time()
            if a.race > 1:
                rr = N.solve_ip_race(B, a.k, extra_rows=extra, nproc=a.race, time_limit=a.ip_time)
                if rr is None:
                    r = dict(status='race timeout', feasible=None, x=None, time=time.time() - t1, nodes=0, bound=None)
                else:
                    x = np.zeros(B.shape[0]); 
                    if rr['sol'] is not None: x[rr['sol']] = 1.0
                    r = dict(status=f"{rr['status']} ({'feas' if rr['cfg']['feasibility'] else 'opt'} seed {rr['cfg']['seed']})",
                             feasible=rr['feasible'], x=(x if rr['sol'] is not None else None), time=rr['time'], nodes=rr['nodes'], bound=None)
            else:
                r = N.solve_ip(B, a.k, threads=a.threads, time_limit=a.ip_time, extra_rows=extra,
                               feasibility=a.feas, verbose=False)
            lpv, _ = N.lp_value(B, threads=a.threads)
            say(f"  IP: {r['status']}, LP {lpv:.4f}, {r['time']:.0f}s, nodes {r['nodes']}, bound {r['bound']}")
            if r['feasible'] is False:
                say(f"  INFEASIBLE: no {a.k}-column hitting set -> h(F) >= {a.k + 1}.  Dumping.")
                L.write_family(f"{out}/final_family.txt", m, Fam, header=f"h(F) >= {a.k + 1} (round {rnd})")
                L.dump_instance(f"{out}/final", m, Fam, B, V, extra=dict(k=a.k, round=rnd, lp=lpv))
                say(f"# done {time.ctime()}"); return
            if r['feasible'] is None:
                say("  IP time limit without a verdict; continuing with a bigger limit")
                a.ip_time *= 2
                continue
            sol = np.nonzero(r['x'] > 0.5)[0]; how = 'IP'
        P = V[sol].copy()
        assert B[sol].any(0).all(), "solution does not cover"
        say(f"  {a.k}-set by {how}: " + " ".join(f"({x:.3f},{y:.3f})" for x, y in P))
        P2, tstar = L.polish_positions(P, m, [sq], rounds=6, verbose=False)
        say(f"  polished: t* = {tstar:+.2e}; " + " ".join(f"({x:.4f},{y:.4f})" for x, y in P2))
        if tstar is not None and tstar >= -1e-9:
            P = P2
        np.save(f"{out}/P_round{rnd:02d}.npy", P)
        prevP = P
        viol, gmin = L.find_violations(P, m, dth_deg=a.dth, pitch=a.pitch, top=800, want=a.add,
                                       min_sep=a.min_sep, verbose=False)
        say(f"  violation search: grid min {gmin:+.3e}, {len(viol)} distinct local minima" + (f", worst {viol[0][0]:+.3e}" if viol else ''))
        if not viol:
            Pex = L.snap_points(P, 1000)
            cert = f"{out}/cand_round{rnd:02d}.txt"; L.write_cert(cert, m, Pex)
            say(f"  no float violation: certifying {cert}")
            ok, txt = LOOP.run_zeromargin(cert, m, f"{out}/zm_round{rnd:02d}.log", full=True, tri=True,
                                          depth=14, nproc=a.threads, oracle=f"{out}/oracle_round{rnd:02d}.txt")
            tail = [l for l in txt.strip().split('\n') if l.strip()][-4:]
            say("  checker: " + " | ".join(tail))
            if ok:
                say(f"  CERTIFIED: {len(Pex)}-point unavoidable set, {cert}"); break
            seeds = LOOP.read_oracle(f"{out}/oracle_round{rnd:02d}.txt")
            Pf = np.array([(float(x), float(y)) for x, y in Pex])
            viol, _ = L.find_violations(Pf, m, dth_deg=a.dth, pitch=a.pitch, top=50, want=a.add,
                                        seeds=seeds[:5000], cutoff=-1e-13, verbose=False)
            if not viol:
                say("  no violation from the oracle seeds either: STALL"); break
        viol = [v for i, v in enumerate(viol) if i < 3 or v[0] <= -a.min_viol]
        added = 0
        Pex = L.snap_points(P, 10 ** 6)
        us = sorted(set(p[0] for p in Fam) | set(F(t) for t in a.angles.split(',') if t))
        nsnap = 0
        for (f, cx, cy, th) in viol:
            key = L.rationalise_pose(m, cx, cy, th, den=1000)
            if a.snap:
                # nearest family angle (normalised), centre on the grid, exact clamp; keep if still violated
                u0 = L.norm_u(F(math.tan(th / 2)).limit_denominator(10 ** 6))
                ub = min(us, key=lambda u: abs(float(u) - float(u0)))
                c, s_ = L.trig(ub); w = abs(c) + abs(s_); lo, hi = w / 2, m - w / 2
                cxs = min(max(F(round(cx * a.snap), a.snap), lo), hi); cys = min(max(F(round(cy * a.snap), a.snap), lo), hi)
                ks = L.pose_key(ub, cxs, cys)
                adm, hits, worst = L.exact_violation(m, Pex, *ks)
                if adm and hits == 0:
                    key = ks; nsnap += 1
            for kk in L.d4_images(m, *key):
                if kk not in Fset:
                    Fset.add(kk); Fam.append(kk); added += 1
        say(f"  added {added} squares ({nsnap} snapped) (worst violations: " + ", ".join(f"{v[0]:+.3f}" for v in viol[:6]) + f"); round time {time.time()-t0:.0f}s")
        if added == 0:
            say("  nothing new: STALL"); break
    say(f"# done {time.ctime()}")


if __name__ == '__main__':
    main()
