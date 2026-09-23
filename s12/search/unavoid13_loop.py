#!/usr/bin/env python3
"""unavoid13: the hitting-set cutting-plane loop (tasks/unavoid13/README.md).

    python3 search/unavoid13_loop.py --m 3 --target 8 --tag t3 [--init tiles,grid2,rot] [--rounds 40]
    python3 search/unavoid13_loop.py --m 4 --target 14 --tag t4 --init tiles,grid2,support:40

Each round: candidates = arrangement vertices of F (float, lenient incidence) -> exact-gap IP
(highspy) -> h(F).  If h(F) >= target: stop (lower-bound certificate; re-solve with
unavoid13_recheck.py).  Else polish the h(F)-point solution's positions (max-margin LP), search
for violated poses (float grid + Nelder-Mead), rationalise, add them and their D4 images to F.
If the float search finds nothing, write the set as a certificate and run zeromargin.py cert.
Everything is logged to runs/unavoid13_<tag>/.
"""
import sys, os, math, time, json, argparse, subprocess
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L


def initial_family(m, spec, support_path=None):
    m = F(m); poses = []
    for item in spec.split(','):
        item = item.strip()
        if item == 'tiles':
            for i in range(int(m)):
                for j in range(int(m)):
                    poses.append((F(0), i + F(1, 2), j + F(1, 2)))
        elif item.startswith('grid'):          # gridK: axis-parallel squares, centres on a 1/K grid
            K = int(item[4:])
            lo, hi = F(1, 2), m - F(1, 2)
            g = [lo + F(k, K) for k in range(int((hi - lo) * K) + 1)]
            for x in g:
                for y in g:
                    poses.append((F(0), x, y))
        elif item.startswith('rot'):           # rotK: a few rational angles on a coarse 1/K grid
            K = int(item[3:]) if len(item) > 3 else 2
            for u in (F(1, 8), F(1, 4), F(1, 3), F(2, 5), F(12, 29), F(1, 2), F(3, 5), F(3, 4), F(7, 8)):
                c, s = L.trig(u); w = abs(c) + abs(s)
                lo, hi = w / 2, m - w / 2
                g = [lo] + [F(k, K) for k in range(int(lo * K) + 1, int(hi * K) + 1)] + [hi]
                g = sorted(set(x for x in g if lo <= x <= hi))
                for x in g:
                    for y in g:
                        poses.append((u, x, y))
        elif item.startswith('support'):       # support:N  -> top-N poses by mass of the exact support
            N = int(item.split(':')[1]) if ':' in item else 10 ** 9
            sup = L.read_support(support_path, m)
            sup.sort(key=lambda t: -t[1])
            poses += [p for p, _ in sup[:N]]
        elif item.startswith('file:'):
            _, fam = L.read_family(item[5:]); poses += fam
        else:
            raise ValueError(item)
    poses = [p for p in poses if L.admissible_exact(m, *p)]
    return L.d4_closure(m, poses)


def run_zeromargin(cert_path, m, log_path, full=True, tri=True, depth=14, nproc=4, oracle=None):
    cmd = [sys.executable, 'search/unavoid13_check.py', 'cert', cert_path, '--depth', str(depth),
           '--nproc', str(nproc), '--seg']
    if tri: cmd.append('--tri')
    if full: cmd.append('--full')
    if oracle: cmd += ['--oracle', oracle]
    with open(log_path, 'w') as f:
        f.write('$ ' + ' '.join(cmd) + '\n'); f.flush()
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, text=True)
    txt = open(log_path).read()
    return ('VERIFIED' in txt.split('\n')[-2:][0] or txt.rstrip().endswith('VERIFIED')) and 'NOT VERIFIED' not in txt, txt


def read_oracle(path):
    out = []
    for line in open(path):
        if line.startswith('#'): continue
        t = line.split()
        if len(t) == 3: out.append((float(t[0]), float(t[1]), float(t[2])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=str, default='4')
    ap.add_argument('--target', type=int, default=14, help='stop when h(F) >= target')
    ap.add_argument('--tag', type=str, required=True)
    ap.add_argument('--init', type=str, default='tiles,grid2')
    ap.add_argument('--support', type=str, default='search/cover4_exact_support.txt')
    ap.add_argument('--rounds', type=int, default=60)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--add', type=int, default=10, help='violated poses added per round (before D4)')
    ap.add_argument('--pitch', type=float, default=0.01)
    ap.add_argument('--dth', type=float, default=0.5)
    ap.add_argument('--den', type=int, default=1000, help='denominator for snapping violated poses')
    ap.add_argument('--pden', type=int, default=1000, help='denominator for snapping points into a certificate')
    ap.add_argument('--tol', type=float, default=1e-7)
    ap.add_argument('--ip-time', type=float, default=None)
    ap.add_argument('--no-polish', action='store_true')
    ap.add_argument('--resume', type=str, default=None, help='family file to start from')
    ap.add_argument('--zm-depth', type=int, default=14)
    ap.add_argument("--dominance-max", type=int, default=2000000)
    a = ap.parse_args()

    m = F(a.m)
    out = f"runs/unavoid13_{a.tag}"; os.makedirs(out, exist_ok=True)
    if a.resume:
        _, Fam = L.read_family(a.resume)
    else:
        Fam = initial_family(m, a.init, a.support)
    Fset = set(Fam)
    log = open(f"{out}/round_log.txt", 'a')
    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()
    say(f"# unavoid13 loop: m = {m}, target {a.target}, init '{a.init}', |F0| = {len(Fam)}, {time.ctime()}")
    hist = []
    for rnd in range(a.rounds):
        t0 = time.time()
        say(f"--- round {rnd}: |F| = {len(Fam)}")
        sq, V, B, rep = L.build_candidates(m, Fam, tol=a.tol, dominance_max=a.dominance_max)
        r = L.solve_hitting_set(B, threads=a.threads, time_limit=a.ip_time)
        say(f"  IP: status {r['status']}, h(F) = {r['obj']}, LP relaxation = {r['lp']:.6f}, "
            f"bound {r['bound']}, {time.time()-t0:.1f}s")
        hist.append(dict(round=rnd, nF=len(Fam), ncand=int(B.shape[0]), nV=int(len(V)), h=r['obj'],
                         lp=r['lp'], status=r['status']))
        L.write_family(f"{out}/F_round{rnd:02d}.txt", m, Fam,
                       header=f"round {rnd}; h(F) = {r['obj']}; LP = {r['lp']:.6f}; candidates {B.shape[0]}")
        if r['obj'] is None:
            say("  IP not solved to optimality; stopping"); break
        if r['obj'] >= a.target:
            say(f"  h(F) = {r['obj']} >= {a.target}: LOWER BOUND REACHED.  Dumping instance.")
            L.dump_instance(f"{out}/final", m, Fam, B, rep,
                            extra=dict(h=r['obj'], lp=r['lp'], round=rnd))
            break
        chosen = np.nonzero(r['x'] > 0.5)[0]
        P = rep[chosen].copy()
        say(f"  points (raw vertices): " + " ".join(f"({x:.4f},{y:.4f})" for x, y in P))
        if not a.no_polish:
            P, tstar = L.polish_positions(P, m, [sq], rounds=5, verbose=False)
            say(f"  polished: t* = {tstar:+.2e}; " + " ".join(f"({x:.4f},{y:.4f})" for x, y in P))
        np.save(f"{out}/P_round{rnd:02d}.npy", P)
        # violation search
        viol, gmin = L.find_violations(P, m, dth_deg=a.dth, pitch=a.pitch, top=400, want=a.add)
        if not viol:
            # nothing found by the float search: certify exactly
            Pex = L.snap_points(P, a.pden)
            cert = f"{out}/cand_round{rnd:02d}.txt"
            L.write_cert(cert, m, Pex)
            say(f"  no float violation (grid min {gmin:+.2e}); running zeromargin cert on {cert}")
            ok, txt = run_zeromargin(cert, m, f"{out}/zm_round{rnd:02d}.log", full=True, tri=True,
                                     depth=a.zm_depth, nproc=a.threads, oracle=f"{out}/oracle_round{rnd:02d}.txt")
            tail = [l for l in txt.strip().split('\n') if l.strip()][-6:]
            say("  zeromargin: " + " | ".join(tail))
            if ok:
                say(f"  CERTIFIED: {len(Pex)}-point unavoidable set, {cert}")
                break
            seeds = read_oracle(f"{out}/oracle_round{rnd:02d}.txt")
            Pf = np.array([(float(x), float(y)) for x, y in Pex])
            viol, _ = L.find_violations(Pf, m, dth_deg=a.dth, pitch=a.pitch, top=50, want=a.add,
                                        seeds=seeds[:5000], cutoff=-1e-13)
            if not viol:
                say("  no violation found even from the oracle seeds: STALL (uncertified boxes, no witness)")
                break
        # rationalise and add
        added = 0; new = []
        Pex = L.snap_points(P, 10 ** 6)
        for (f, cx, cy, th) in viol:
            key = L.rationalise_pose(m, cx, cy, th, den=a.den)
            adm, hits, worst = L.exact_violation(m, Pex, *key)
            if not adm:
                key = L.rationalise_pose(m, cx, cy, th, den=a.den * 1000)
                adm, hits, worst = L.exact_violation(m, Pex, *key)
            imgs = L.d4_images(m, *key)
            fresh = [k for k in imgs if k not in Fset]
            say(f"  violation f = {f:+.3e} at ({cx:.5f},{cy:.5f},{math.degrees(th):.3f}deg) -> "
                f"u={key[0]} cx={key[1]} cy={key[2]}; exact: adm {adm}, hits {hits}, "
                f"maxdepth {float(worst) if worst is not None else None:+.3e}; +{len(fresh)} squares")
            for k in fresh:
                Fset.add(k); Fam.append(k); added += 1
        say(f"  added {added} squares; round time {time.time()-t0:.1f}s")
        if added == 0:
            say("  nothing new to add: STALL"); break
    with open(f"{out}/history.json", 'w') as f:
        json.dump(hist, f, indent=1)
    say(f"# done {time.ctime()}")


if __name__ == '__main__':
    main()
