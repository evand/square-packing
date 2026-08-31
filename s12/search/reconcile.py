#!/usr/bin/env python3
"""Reconcile the packing-side and cover-side anchor-clique gains  (task J).

On ONE finite instance -- a finite pose set P (cover rows = packing variables) and a finite
family C of anchor cliques (cover columns = packing constraints, with the point cliques
K(p,{p}) among them) -- the cover LP and the packing LP are exact duals:

    packing:  max  1.mu    s.t.  mu(K) <= 1  for every K in C,  mu >= 0     (mu on P)
    cover:    min  1.y     s.t.  sum_{K in C : S in K} y_K >= 1  for every S in P,  y >= 0

so they have the SAME value.  Task G measured 0.10-0.12 on the packing side at t = 3.99 and
task I measured 0.0004 on the cover side at t = 3.98; this module puts both sides on one
instance and finds out which pose set / separation depth explains the gap.

Subcommands
    build     rebuild task G's converged anchor-clique instance from its support + cut dump,
              converge rows and cuts, solve BOTH LPs and print the duality table
    price     from that instance, add poses by the TRUE reduced cost of the clique LP
              (1 - capture - sum_K y_K [S in K]) rather than by capture alone, re-converge
              rows and cuts, and report the ladder
    escape    diagnostic: for the converged instance, how much reduced cost the sweep pricer's
              candidates really have once the cut duals are charged

Everything here is float and heuristic; a number that is meant as a bound says so.
"""
import sys, os, math, time, json, argparse
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
MAIN_RUNS = '/home/evand/math/square-packing/s12/runs'
sys.path.insert(0, HERE)
import packing_dual as PD
import clique_family as CF
import clique_continuum as CC

TOL = 1e-9


# ------------------------------------------------------------------ instance construction
def build_model(t, threads, log, supports, cutfile, row_pitch=0.05, extra_poses=None):
    M = CC.CModel(t, threads, log)
    seeds = []
    for path in supports:
        ps = CC.read_support_any(path)
        log(f"# warm {path}: {len(ps)} poses")
        seeds += ps
    if extra_poses:
        seeds += list(extra_poses)
    M.add_poses(seeds)
    M.add_points(PD.seed_points(t, row_pitch))
    if cutfile:
        CC.load_cuts(M, cutfile, log)
    log(f"# built: {len(M.poses)} orbits, {M.m.A.shape[0]} rows, {len(M.cuts)} cuts")
    return M


class SepArgs:
    """the separator knobs of clique_continuum.separate, as a plain object"""
    mode = 'anchor'
    sep_pts = 40
    nseed = 4
    maxstep = 12
    neps = 6
    nrho = 5
    ndir = 8
    cut_per_round = 12
    clique_tl = 20.0
    kcut_nth = 120


def converge(M, log, sargs, maxit=60, tag='', tlimit=1e9, pure_too=True):
    """inner loop on a fixed pose set: add violated coverage rows and violated anchor cuts
    until neither exists.  Returns (mu, y, ycut, obj, converged).

    `pure_too` also certifies the CUTS-OFF solution on the same poses and adds the rows it
    violates.  Without that the matched pure value is computed on a row set converged only for
    the clique solution, so it is an UPPER bound on the pure LP's value on those poses and the
    reported clique gain is an over-estimate -- which is what `clique_continuum.py` does
    (`pure_obj` is one `M.solve()` on the clique loop's rows)."""
    t0 = time.time()
    last = None
    for it in range(maxit):
        r = M.solve()
        if r is None:
            log("  LP FAILED")
            return None
        mu, y, ycut, obj = r
        last = r
        Mx, bad, nv, ovf = M.m.certify(mu, thresh=1.0 + 1e-7, cap=200000)
        nadd = M.add_points(bad[:, :2]) if len(bad) else 0
        npure = 0
        if pure_too and M.cuts:
            rp = pure_solution(M)
            if rp is not None:
                Mp, badp, _, _ = M.m.certify(rp[0], thresh=1.0 + 1e-7, cap=200000)
                npure = M.add_points(badp[:, :2]) if len(badp) else 0
        ncut, mcl = CC.separate(M, mu, sargs, log)
        if it % 10 == 0 or (nadd == 0 and ncut == 0 and npure == 0):
            log(f"  {tag}it{it}: obj {obj:.6f} cov {Mx:.6f} rows {M.m.A.shape[0]} "
                f"cols {len(M.poses)} cuts {len(M.cuts)} +rows {nadd} +pure-rows {npure} "
                f"+cuts {ncut} maxanchor {mcl:.4f}")
        if nadd == 0 and ncut == 0 and npure == 0:
            return last + (True,)
        if time.time() - t0 > tlimit:
            log(f"  {tag}time limit in inner loop")
            break
    return last + (False,)


def pure_solution(M):
    """the same LP on the same poses and rows with the cuts removed (the matched pair)"""
    saved = M.cuts, M.Ccols
    M.cuts, M.Ccols = [], []
    r = M.solve()
    M.cuts, M.Ccols = saved
    return r


def pure_value(M):
    r = pure_solution(M)
    return None if r is None else r[3]


# ------------------------------------------------------------------ the cover LP (the dual)
def cover_lp(M):
    """solve the COVER LP of the same instance from scratch:

        min sum_i y_i   s.t.  A^T y >= 1 (one constraint per pose orbit),  y >= 0

    where A is the packing LP's constraint matrix (point rows then cut rows).  Column i of A
    is the pose orbit i; row j is the clique K_j (a point clique for a point row, an anchor
    clique for a cut row).  Returns (value, y).  This is a genuinely separate LP solve, not a
    reading of the packing LP's duals."""
    A = M.m.A
    if M.cuts:
        A = sp.vstack([A, sp.csr_matrix(np.array(M.Ccols))], format='csr')
    m, n = A.shape
    res = linprog(c=np.ones(m), A_ub=(-A.T).tocsr(), b_ub=-np.ones(n),
                  bounds=(0, None), method='highs')
    if not res.success:
        return None, None
    return float(res.fun), np.maximum(res.x, 0.0)


# ------------------------------------------------------------------ true-reduced-cost pricing
def cut_penalty(M, ycut, poses, active=None):
    """sum over cuts of y_K * [pose in K], evaluated on single poses (n,3)"""
    poses = np.asarray(poses, dtype=float).reshape(-1, 3)
    pen = np.zeros(len(poses))
    if not M.cuts:
        return pen
    idx = range(len(M.cuts)) if active is None else active
    for i in idx:
        d = ycut[i]
        if d <= 1e-12:
            continue
        pen += d * M.cuts[i]['obj'].members(poses)
    return pen


def orbit_reduced_cost(M, lib, ax, ay, aw, ycut, poses, threads):
    """the reduced cost of the ORBIT of each pose: the packing column is the orbit, whose
    coefficient in a cut row is (#images in the clique)/8 and in a point row is the summed
    incidence/8.  Equivalently: rc = 1 - mean over the 8 images of (capture + cut penalty)."""
    poses = np.asarray(poses, dtype=float).reshape(-1, 3)
    if len(poses) == 0:
        return np.zeros(0)
    imgs = []
    for (cx, cy, th) in poses:
        imgs += PD.pose_images(cx, cy, th, M.t)
    imgs = np.array(imgs, dtype=float).reshape(-1, 3)
    cap = PD.capture(lib, ax, ay, aw, imgs, threads)
    pen = cut_penalty(M, ycut, imgs)
    return 1.0 - (cap + pen).reshape(len(poses), 8).mean(axis=1)


def price_true(M, lib, ax, ay, aw, ycut, angs, h, threads, want, log, chunk_angles=1):
    """price the centre lattice at pitch h against the TRUE reduced cost of the clique LP.

    The point of the exercise: clique_continuum's h = 0 path ranks candidates by capture only
    (`got.sort(key=capture)`), which is the pure LP's reduced cost, so with 900 cut duals in
    play the columns it adds need not improve the clique LP at all."""
    best = []
    npos = 0
    for th in angs:
        P = CC.lattice_slice(M.t, h, th)
        if len(P) == 0:
            continue
        npos += len(P)
        cap = PD.capture(lib, ax, ay, aw, P, threads)
        rc = 1.0 - cap
        keep = np.nonzero(rc > -1e-9)[0]           # cut duals only lower rc, so this is a superset
        if len(keep) == 0:
            continue
        Pk = P[keep]
        pen = cut_penalty(M, ycut, Pk)
        rck = rc[keep] - pen
        good = np.nonzero(rck > 1e-9)[0]
        if len(good):
            sel = good[np.argsort(-rck[good])[:want]]
            best.append(np.c_[Pk[sel], rck[sel], rc[keep][sel]])
    if not best:
        log(f"  price_true: no improving pose on {npos} lattice poses (h={h})")
        return np.zeros((0, 5))
    B = np.vstack(best)
    B = B[np.argsort(-B[:, 3])][:want]
    return B


def price_neighbours(M, lib, ax, ay, aw, ycut, mu, threads, want, log,
                     steps=(0.04, 0.02, 0.01, 0.004, 0.001), dths=(2.0, 1.0, 0.5, 0.2, 0.05)):
    """perturbations of the CURRENT SUPPORT, scored by the true reduced cost.  These are the
    poses a cut has to be escaped by: a support pose that sits inside {S : A subset S} has a
    neighbour just outside it with almost the same capture and no cut penalty."""
    sup = [M.poses[i] for i in np.nonzero(mu > 1e-11)[0]]
    if not sup:
        return np.zeros((0, 5))
    cand = []
    for (cx, cy, th) in sup:
        for d, dd in zip(steps, dths):
            dt = math.radians(dd)
            for (a, b, c) in ((d, 0, 0), (-d, 0, 0), (0, d, 0), (0, -d, 0), (0, 0, dt), (0, 0, -dt),
                              (d, d, 0), (-d, -d, 0), (d, -d, 0), (-d, d, 0),
                              (d, 0, dt), (-d, 0, -dt), (0, d, dt), (0, -d, -dt)):
                cand.append((cx + a, cy + b, th + c))
    cand = [PD.clamp_pose(*p, M.t) for p in cand]
    cand = np.array([p for p in cand if p is not None], dtype=float)
    if len(cand) == 0:
        return np.zeros((0, 5))
    rc = orbit_reduced_cost(M, lib, ax, ay, aw, ycut, cand, threads)
    good = np.nonzero(rc > 1e-9)[0]
    if len(good) == 0:
        log(f"  neighbours: none of {len(cand)} support perturbations improves")
        return np.zeros((0, 5))
    sel = good[np.argsort(-rc[good])[:want]]
    return np.c_[cand[sel], rc[sel], rc[sel]]


# ------------------------------------------------------------------ commands
def cmd_build(a, log):
    t = a.t
    sargs = SepArgs()
    sargs.cut_per_round = a.cut_per_round
    M = build_model(t, a.threads, log, a.support, a.cuts, row_pitch=a.row_pitch)
    r = converge(M, log, sargs, maxit=a.maxit, tag='b')
    mu, y, ycut, obj, conv = r
    pure = pure_value(M)
    log(f"* PACKING: clique {obj:.6f}  pure(same poses) {pure:.6f}  gain {pure-obj:.6f} "
        f"rows {M.m.A.shape[0]} cols {len(M.poses)} cuts {len(M.cuts)} converged={conv}")
    t0 = time.time()
    cv, ycov = cover_lp(M)
    log(f"* COVER  : value {cv:.6f}  ({time.time()-t0:.0f}s)   |packing - cover| = {abs(cv-obj):.3e}")
    save_state(M, mu, a.tag, t, log)
    return M, mu, y, ycut, obj


def save_state(M, mu, tag, t, log):
    sp_path = os.path.join(RUNS, f'rc_{tag}_support.txt')
    with open(sp_path, 'w') as f:
        f.write(f"# t={t}\n")
        for i, (cx, cy, th) in enumerate(M.poses):
            if mu[i] > 1e-12:
                f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} {mu[i]:.12f}\n")
    pose_path = os.path.join(RUNS, f'rc_{tag}_poses.txt')
    with open(pose_path, 'w') as f:
        f.write(f"# t={t} full pose set ({len(M.poses)} orbits)\n")
        for (cx, cy, th) in M.poses:
            f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} 0\n")
    CC.dump_cuts(M, tag)
    log(f"# wrote {sp_path}, {pose_path}, runs/cq_{tag}_cuts.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('build', 'price', 'escape'))
    ap.add_argument('tag')
    ap.add_argument('--t', type=float, default=3.99)
    ap.add_argument('--support', action='append', default=[])
    ap.add_argument('--cuts', default=None)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--row-pitch', type=float, default=0.05)
    ap.add_argument('--maxit', type=int, default=60)
    ap.add_argument('--cut-per-round', type=int, default=12)
    ap.add_argument('--rounds', type=int, default=20)
    ap.add_argument('--h', type=float, default=0.02, help='pricing lattice pitch')
    ap.add_argument('--dth', type=float, default=1.0)
    ap.add_argument('--want', type=int, default=400)
    ap.add_argument('--maxatoms', type=int, default=8000)
    ap.add_argument('--time', type=float, default=1e9)
    ap.add_argument('--no-neighbours', action='store_true')
    ap.add_argument('--sweep-too', action='store_true',
                    help="also add task G's sweep candidates (capture-ranked), for comparison")
    a = ap.parse_args()

    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'rc_{a.tag}.log'), 'w')

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + '\n')
        logf.flush()

    log(f"# reconcile {a.cmd} {a.tag}: t={a.t} threads={a.threads} h={a.h}")
    if a.cmd == 'build':
        cmd_build(a, log)
        return

    # ---- price / escape both need the built instance
    sargs = SepArgs()
    sargs.cut_per_round = a.cut_per_round
    M = build_model(a.t, a.threads, log, a.support, a.cuts, row_pitch=a.row_pitch)
    lib = M.m.lib
    r = converge(M, log, sargs, maxit=a.maxit, tag='b')
    mu, y, ycut, obj, conv = r
    pure = pure_value(M)
    cv, _ = cover_lp(M)
    log(f"* r-1: STAGE clique {obj:.6f}  pure {pure:.6f}  gain {pure-obj:.6f}  "
        f"cover {cv:.6f}  |dual gap| {abs(cv-obj):.3e}  rows {M.m.A.shape[0]} "
        f"cols {len(M.poses)} cuts {len(M.cuts)} conv={conv}")
    angs = CC.lattice_angles(a.dth)

    if a.cmd == 'escape':
        ax, ay, aw = PD.atoms_from_dual(M.m, y, maxrows=a.maxatoms)
        got = PD.price_sweep(M.m, lib, ax, ay, aw, angs, 6, 0.02, a.threads)
        got = [g for g in got if g[3] < 1.0 - 1e-9]
        got.sort(key=lambda g: g[3])
        cand = np.array([[g[0], g[1], g[2]] for g in got[:540]], dtype=float)
        capt = np.array([g[3] for g in got[:540]])
        rc_true = orbit_reduced_cost(M, lib, ax, ay, aw, ycut, cand, a.threads)
        log(f"# task G's sweep candidates: {len(cand)}")
        log(f"#   capture-only reduced cost 1-cap: min {1-capt.max():.6f} "
            f"med {np.median(1-capt):.6f} max {1-capt.min():.6f}")
        log(f"#   TRUE orbit reduced cost:         min {rc_true.min():.6f} "
            f"med {np.median(rc_true):.6f} max {rc_true.max():.6f}")
        log(f"#   improving (rc > 1e-9): {int((rc_true > 1e-9).sum())} of {len(cand)}")
        nb = price_neighbours(M, lib, ax, ay, aw, ycut, mu, a.threads, 4000, log)
        if len(nb):
            log(f"#   support perturbations with rc > 0: {len(nb)}, best {nb[0,3]:.6f}")
        # how much dual mass sits on the cuts
        tot_pt = float(y.sum())
        tot_cut = float(np.sum(ycut))
        log(f"#   dual mass: points {tot_pt:.6f}  cuts {tot_cut:.6f}  "
            f"(sum {tot_pt+tot_cut:.6f} = LP value)")
        nact = int((ycut > 1e-9).sum())
        log(f"#   cuts with positive dual: {nact} of {len(M.cuts)}")
        return

    # ---- price: the ladder with a correct pricer
    t0 = time.time()
    hist = []
    for rnd in range(a.rounds):
        ax, ay, aw = PD.atoms_from_dual(M.m, y, maxrows=a.maxatoms)
        cands = []
        nb = np.zeros((0, 5))
        if not a.no_neighbours:
            nb = price_neighbours(M, lib, ax, ay, aw, ycut, mu, a.threads, a.want, log)
            if len(nb):
                log(f"  neighbours: {len(nb)} improving, best true rc {nb[0,3]:.6f}")
                cands += [(c[0], c[1], c[2]) for c in nb]
        pl = price_true(M, lib, ax, ay, aw, ycut, angs, a.h, a.threads, a.want, log)
        if len(pl):
            log(f"  lattice(h={a.h}): {len(pl)} improving, best true rc {pl[0,3]:.6f} "
                f"(capture-only rc there {pl[0,4]:.6f})")
            cands += [(c[0], c[1], c[2]) for c in pl]
        if a.sweep_too:
            got = PD.price_sweep(M.m, lib, ax, ay, aw, angs, 6, 0.02, a.threads)
            got = [g for g in got if g[3] < 1.0 - 1e-9]
            got.sort(key=lambda g: g[3])
            cands += [(g[0], g[1], g[2]) for g in got[:a.want]]
        if not cands:
            log("* PRICING CLEAN on the true reduced cost: no improving pose found")
            break
        nnew = M.add_poses(cands)
        log(f"  added {nnew} new orbits -> {len(M.poses)}")
        if nnew == 0:
            log("* no new columns")
            break
        r = converge(M, log, sargs, maxit=a.maxit, tag=f'r{rnd} ')
        if r is None:
            break
        mu, y, ycut, obj, conv = r
        pure = pure_value(M)
        mcl = CC.residual_clique(M, mu, sargs)
        log(f"* r{rnd}: STAGE clique {obj:.6f}  pure {pure:.6f}  gain {pure-obj:.6f}  "
            f"rows {M.m.A.shape[0]} cols {len(M.poses)} cuts {len(M.cuts)} "
            f"maxclique {mcl:.4f} conv={conv}  ({time.time()-t0:.0f}s)")
        hist.append(dict(rnd=rnd, obj=obj, pure=pure, gain=pure - obj,
                         rows=int(M.m.A.shape[0]), cols=len(M.poses), cuts=len(M.cuts),
                         conv=bool(conv)))
        json.dump(hist, open(os.path.join(RUNS, f'rc_{a.tag}.json'), 'w'), indent=1)
        save_state(M, mu, a.tag, a.t, log)
        if time.time() - t0 > a.time:
            log("# time limit")
            break
    cv, _ = cover_lp(M)
    log(f"# FINAL {a.tag}: clique {obj:.6f} cover {cv:.6f} |gap| {abs(cv-obj):.3e} "
        f"cols {len(M.poses)} cuts {len(M.cuts)}")


if __name__ == '__main__':
    main()
