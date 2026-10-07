#!/usr/bin/env python3
"""unavoid13: the cutting-plane loop, second form (lazy rows + feasibility IP).

    python3 search/unavoid13_loop2.py --m 4 --k 13 --tag t4L --resume runs/unavoid13_t4b/F_round08.txt

Compared with unavoid13_loop.py:
  * the IP is the FEASIBILITY question "is there a k-point hitting set of the active family?"
    (sum x <= k, no objective): any incumbent is enough to generate violations; only a proof of
    infeasibility stops the loop -- and infeasibility on the active family F' ⊆ F is
    infeasibility on F (fewer rows to hit), so it is the lower bound h(F') >= k + 1;
  * lazy rows: the IP runs on an active subset F' of the whole family F; the k-set it returns is
    then checked against all of F (float, lenient) and the squares it misses are added to F'
    (inner loop) until it hits all of F; only then are new geometric violations searched and
    added to F and F';
  * a time limit on the IP: an incumbent within the limit is used; a time-out without incumbent
    doubles the limit.
Round log: runs/unavoid13_<tag>/round_log.txt.  Families: F_round*.txt (whole) and
Fact_round*.txt (active); on infeasibility the active family and the instance are dumped as
final_* for unavoid13_recheck.py.
"""
import sys, os, math, time, json, argparse
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L
import unavoid13_loop as LOOP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=str, default='4')
    ap.add_argument('--k', type=int, default=13)
    ap.add_argument('--tag', type=str, required=True)
    ap.add_argument('--init', type=str, default='tiles,grid2,support:40')
    ap.add_argument('--support', type=str, default='search/cover4_exact_support.txt')
    ap.add_argument('--resume', type=str, default=None)
    ap.add_argument('--resume-active', type=str, default=None)
    ap.add_argument('--rounds', type=int, default=200)
    ap.add_argument('--threads', type=int, default=3)
    ap.add_argument('--add', type=int, default=20)
    ap.add_argument('--lazy-add', type=int, default=40, help='missed squares of F added to the active set per inner step')
    ap.add_argument('--pitch', type=float, default=0.01)
    ap.add_argument('--dth', type=float, default=0.5)
    ap.add_argument('--pden', type=int, default=1000)
    ap.add_argument('--ip-time', type=float, default=600)
    ap.add_argument('--zm-depth', type=int, default=14)
    a = ap.parse_args()
    m = F(a.m); mf = float(m)
    out = f"runs/unavoid13_{a.tag}"; os.makedirs(out, exist_ok=True)
    Fam = L.read_family(a.resume)[1] if a.resume else LOOP.initial_family(m, a.init, a.support)
    Fset = set(Fam)
    if a.resume_active:
        Fact = [p for p in L.read_family(a.resume_active)[1] if p in Fset]
    else:
        # start the active set with the tiles and the axis-parallel grid squares that are in F
        Fact = [p for p in Fam if p[0] == 0]
    Aset = set(Fact)
    log = open(f"{out}/round_log.txt", 'a')
    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()
    say(f"# unavoid13 loop2 (lazy rows, feasibility IP): m = {m}, k = {a.k}, |F| = {len(Fam)}, |F'| = {len(Fact)}, {time.ctime()}")
    ip_time = a.ip_time
    for rnd in range(a.rounds):
        t0 = time.time()
        say(f"--- round {rnd}: |F| = {len(Fam)}, |F'| = {len(Fact)}")
        L.write_family(f"{out}/F_round{rnd:02d}.txt", m, Fam, header=f"round {rnd}")
        # ---- inner loop: IP on the active set until its k-set hits all of F
        sqF = L.Squares(Fam)
        inner = 0
        while True:
            inner += 1
            t1 = time.time()
            sq, V, B, rep = L.build_candidates(m, Fact, verbose=False)
            r = L.solve_feasible(B, a.k, threads=a.threads, time_limit=ip_time)
            lpv = L.lp_relaxation(B, a.threads)
            say(f"  [{inner}] IP on |F'| = {len(Fact)} ({B.shape[0]} cand): {r['status']}, LP {lpv:.4f}, {time.time()-t1:.0f}s")
            if r['feasible'] is False:
                say(f"  INFEASIBLE: no {a.k}-point set hits the active family -> h >= {a.k + 1}.  Dumping.")
                L.write_family(f"{out}/final_family.txt", m, Fact, header=f"active family with h >= {a.k + 1} (round {rnd})")
                L.dump_instance(f"{out}/final", m, Fact, B, rep, extra=dict(k=a.k, round=rnd, lp=lpv))
                say(f"# done {time.ctime()}"); return
            if r['feasible'] is None:
                ip_time *= 2
                say(f"  time limit without incumbent; doubling the limit to {ip_time:.0f}s")
                continue
            P = rep[np.nonzero(r['x'] > 0.5)[0]].copy()
            D = sqF.depth(P).max(0)              # (|F|,) best depth per square
            missed = np.nonzero(D < -1e-7)[0]
            if len(missed) == 0:
                break
            order = missed[np.argsort(D[missed])][:a.lazy_add]
            for j in order:
                if Fam[j] not in Aset:
                    Aset.add(Fam[j]); Fact.append(Fam[j])
            say(f"      {a.k}-set misses {len(missed)} squares of F (worst depth {D[missed].min():+.3e}); +{len(order)} to F'")
        L.write_family(f"{out}/Fact_round{rnd:02d}.txt", m, Fact, header=f"active set, round {rnd}")
        say(f"  {a.k}-set hits all of F after {inner} inner steps; raw: " + " ".join(f"({x:.3f},{y:.3f})" for x, y in P))
        P, tstar = L.polish_positions(P, m, [sqF], rounds=5, verbose=False)
        say(f"  polished: t* = {tstar:+.2e}; " + " ".join(f"({x:.4f},{y:.4f})" for x, y in P))
        np.save(f"{out}/P_round{rnd:02d}.npy", P)
        viol, gmin = L.find_violations(P, m, dth_deg=a.dth, pitch=a.pitch, top=600, want=a.add, verbose=False)
        say(f"  violation search: grid min {gmin:+.3e}, {len(viol)} distinct local minima" + (f", worst {viol[0][0]:+.3e}" if viol else ''))
        if not viol:
            Pex = L.snap_points(P, a.pden)
            cert = f"{out}/cand_round{rnd:02d}.txt"; L.write_cert(cert, m, Pex)
            say(f"  no float violation: certifying {cert}")
            ok, txt = LOOP.run_zeromargin(cert, m, f"{out}/zm_round{rnd:02d}.log", full=True, tri=True,
                                          depth=a.zm_depth, nproc=a.threads, oracle=f"{out}/oracle_round{rnd:02d}.txt")
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
        added = 0
        for (f, cx, cy, th) in viol:
            key = L.rationalise_pose(m, cx, cy, th, den=1000)
            for kk in L.d4_images(m, *key):
                if kk not in Fset:
                    Fset.add(kk); Fam.append(kk); added += 1
                if kk not in Aset:
                    Aset.add(kk); Fact.append(kk)
        say(f"  added {added} squares to F (worst violations: " + ", ".join(f"{v[0]:+.3f}" for v in viol[:5]) + f"); round time {time.time()-t0:.0f}s")
        if added == 0:
            say("  nothing new: STALL"); break
    say(f"# done {time.ctime()}")


if __name__ == '__main__':
    main()
