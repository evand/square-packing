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
def segments_only(log):
    """restrict the separator to SEGMENT anchors.

    `clique_family.anchor_separate` builds `A` as the intersection of support squares, i.e. a
    polygon.  A polygon anchor is a valid clique but it is (a) not representable in a certificate
    -- `certificates/FORMAT.md` has `anchorP` (a point) and `anchorS` (a segment), nothing else --
    and (b) typically of tiny or zero area, so `{S : A subseteq S}` is a measure-zero set of poses:
    it can cut a finite pose list and can never pay for itself in a cover of the continuum.
    `anchor_local` is the Lemma-2 segment family, which is representable."""
    orig = CF.anchor_separate_all

    def seg_all(IM, p, t, nseed=8, maxstep=40, tol=1e-12, neps=10, nrho=8, ndir=16):
        return CF.anchor_local(IM, p, t, neps=neps, nrho=nrho, tol=tol, ndir=ndir)

    CF.anchor_separate_all = seg_all
    CC.CF.anchor_separate_all = seg_all
    log("# separator restricted to SEGMENT anchors (certificate-representable)")
    return orig


def build_model(t, threads, log, supports, cutfile, row_pitch=0.05, extra_poses=None,
                seg_only=False):
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
        if seg_only:
            keep = [len(c['obj'].A) == 2 for c in M.cuts]
            log(f"# dropping {len(keep) - sum(keep)} polygon-anchored cuts of {len(keep)}")
            M.drop_cuts(keep)
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
def refine_true(M, lib, ax, ay, aw, ycut, cand, threads, rounds=3,
                steps=(0.02, 0.005, 0.0012), dths=(1.0, 0.25, 0.06)):
    """coordinate descent on the TRUE orbit reduced cost, so the priced columns are not stuck on
    the pricing lattice (the analogue of packing_dual.refine, which maximises 1 - capture)."""
    cur = np.asarray(cand, dtype=float).reshape(-1, 3)
    if len(cur) == 0:
        return cur, np.zeros(0)
    val = orbit_reduced_cost(M, lib, ax, ay, aw, ycut, cur, threads)
    for h, dd in zip(steps, dths):
        dt = math.radians(dd)
        moves = np.array([[h, 0, 0], [-h, 0, 0], [0, h, 0], [0, -h, 0], [0, 0, dt], [0, 0, -dt],
                          [h, h, 0], [-h, -h, 0], [h, -h, 0], [-h, h, 0]])
        for _ in range(rounds):
            trial = (cur[:, None, :] + moves[None, :, :]).reshape(-1, 3)
            cl = [PD.clamp_pose(*p, M.t) for p in trial]
            ok = np.array([p is not None for p in cl])
            T = np.array([p if p is not None else (0.0, 0.0, 0.0) for p in cl], dtype=float)
            v = orbit_reduced_cost(M, lib, ax, ay, aw, ycut, T, threads)
            v[~ok] = -9e9
            v = v.reshape(len(cur), len(moves))
            T = T.reshape(len(cur), len(moves), 3)
            j = np.argmax(v, axis=1)
            best = v[np.arange(len(cur)), j]
            better = best > val + 1e-12
            if not better.any():
                break
            cur[better] = T[np.arange(len(cur)), j][better]
            val = np.maximum(val, best)
    return cur, val


def cmd_build(a, log):
    t = a.t
    sargs = SepArgs()
    sargs.cut_per_round = a.cut_per_round
    if a.seg_only:
        segments_only(log)
    M = build_model(t, a.threads, log, a.support, a.cuts, row_pitch=a.row_pitch,
                    seg_only=a.seg_only)
    r = converge(M, log, sargs, maxit=a.maxit, tag='b', pure_too=not a.no_pure_rows)
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


def cmd_xmember(a, log):
    """compare mu(K) for IDENTICAL cliques computed by the two code paths.

    G's path (`clique_family.contains_np` / `meets_np` via `clique_continuum.ACut`) tests the
    exact closed unit square of a single pose.  I's path (`anchorclique.member`) tests a *cell*
    of the verifier's sweep: `contains` by the concentric sigma_k-square (smaller than the unit
    square) and `meets` by the hexagon `A + [-h,h]^2` (larger).  With h = 1/2 the two should
    agree; with h = sigma_k(N)/2 the cover side is conservative and credits less.
    """
    import anchorclique as AC
    from fractions import Fraction as Fr
    import tighten as T
    rng = np.random.default_rng(a.seed)
    t = a.t
    D = 100000
    Ns = [int(v) for v in a.Ns.split(',')]
    hs = [('exact 1/2', 0.5)] + [(f'sigma_k/2 N={N}', None) for N in Ns]
    log("# clique                                  poses  G-in  I-in  G\\I   I\\G   mass_G   mass_I")
    tot = {k: [0, 0, 0] for k, _ in hs}
    ncl = 0
    for trial in range(a.n):
        # a wall point p and the Lemma-2 anchor, exactly as anchorsep/branch.py builds them
        d = float(rng.uniform(0.35, 0.995))
        py = float(rng.uniform(1.0, t / 2))
        X, Y = int(round(d * D)), int(round(py * D))
        params = AC.cand_params(Fr(t).limit_denominator(10000), D, X, Y,
                                frac=float(rng.choice([0.95, 0.7, 0.4])))
        params = [q for q in params if q[0] == 0]
        if not params:
            continue
        wall, en, rn, dd = params[0]
        cl = AC.kpa(Fr(t).limit_denominator(10000), D, X, Y, wall, en, rn)
        if cl is None:
            continue
        ncl += 1
        p = (X / D, Y / D)
        A = [(float(cl[0][1][1]), float(cl[0][1][2])), (float(cl[0][1][3]), float(cl[0][1][4]))]
        acut = CC.ACut(p, A)
        # poses: random admissible poses concentrated near p (where the clique lives)
        n = a.poses
        th = rng.uniform(0, math.pi / 2, n)
        w2 = np.array([PD.wid_of(v) / 2 for v in th])
        cx = np.clip(p[0] + rng.uniform(-1.2, 1.2, n), w2, t - w2)
        cy = np.clip(p[1] + rng.uniform(-1.2, 1.2, n), w2, t - w2)
        P = np.c_[cx, cy, th]
        mG = acut.members(P)
        for name, hv in hs:
            if hv is None:
                N = int(name.split('=')[1])
                # the verifier's net is uniform in u = tan(theta/2), k = floor(u N), u in [0,1)
                k = np.minimum(np.floor(np.tan(th / 2) * N).astype(int), N - 1)
                h = np.array([T.sigma_k(N, int(kk)) / 2 for kk in k])
            else:
                h = np.full(n, hv)
            rows = np.c_[P, h]
            mI = AC.member(cl, rows)
            tot[name][0] += int(mG.sum())
            tot[name][1] += int(mI.sum())
            tot[name][2] += int((mG & ~mI).sum())
        if trial < a.show:
            for name, hv in hs:
                pass
    log(f"# {ncl} cliques x {a.poses} random poses each")
    for name, _ in hs:
        g, i, gi = tot[name]
        log(f"  {name:20s}  G credits {g:8d}   I credits {i:8d}   G\\I {gi:7d}  "
            f"({100.0*gi/max(g,1):.4f} % of G's)")


def probe_cover(M, ycov, poses, threads):
    """covered weight of each pose under the cover solution ycov (point rows then cut rows):
    sum over point cliques p with p in S of y_p, plus sum over anchor cliques K with S in K.

    The packing LP's value is the cover's cost ONLY on the poses that are columns of the LP;
    this is what the same cover pays on an arbitrary pose, i.e. what the verifier would see."""
    poses = np.asarray(poses, dtype=float).reshape(-1, 3)
    npt = M.m.A.shape[0]
    y = ycov[:npt]
    ycut = ycov[npt:]
    IM = np.c_[poses, np.zeros(len(poses))]
    tot = np.zeros(len(poses))
    sel = np.nonzero(y > 1e-12)[0]
    t = M.t
    for i in sel:
        px, py = M.m.pts[i]
        # WITH multiplicity: packing_dual.images_array emits all 8 images of a pose whether or
        # not they coincide, and the incidence matrix sums the duplicates, so the point side
        # must count coincident images too (a point on a symmetry axis has fewer distinct
        # images but the same coefficient)
        for (qx, qy) in [(px, py), (t - px, py), (px, t - py), (t - px, t - py),
                         (py, px), (t - py, px), (py, t - px), (t - py, t - px)]:
            tot += (y[i] / 8.0) * CF.contains_np(IM, [(qx, qy)])
    act = np.nonzero(ycut > 1e-12)[0]
    if len(act):
        # the cut row's coefficient for a pose orbit is (#images of the orbit in K)/8
        imgs = []
        for (cx, cy, th) in poses:
            imgs += PD.pose_images(cx, cy, th, t)
        imgs = np.ascontiguousarray(np.array(imgs, dtype=float).reshape(-1, 3))
        for j in act:
            m = M.cuts[j]['obj'].members(imgs).reshape(len(poses), 8)
            tot += ycut[j] * m.mean(axis=1)
    return tot


def cmd_probe(a, log):
    """the cover-side probe: what does the converged cover of the finite instance pay on poses
    that are not columns of it?"""
    sargs = SepArgs()
    sargs.cut_per_round = a.cut_per_round
    if a.seg_only:
        segments_only(log)
    M = build_model(a.t, a.threads, log, a.support, a.cuts, row_pitch=a.row_pitch,
                    seg_only=a.seg_only)
    r = converge(M, log, sargs, maxit=a.maxit, tag='b')
    mu, y, ycut, obj, conv = r
    cv, ycov = cover_lp(M)
    log(f"* instance: packing {obj:.6f}  cover {cv:.6f}  cols {len(M.poses)} cuts {len(M.cuts)}")
    npt = M.m.A.shape[0]
    log(f"  cover uses {int((ycov[:npt] > 1e-12).sum())} point cliques (weight "
        f"{ycov[:npt].sum():.4f}) and {int((ycov[npt:] > 1e-12).sum())} anchor cliques (weight "
        f"{ycov[npt:].sum():.4f})")
    # what SHAPE are the anchors the cover actually pays for?  Only a segment anchor is
    # representable in a certificate (`anchorS`); an anchor that is a polygon -- the intersection
    # of support squares, which is what clique_family.anchor_separate returns -- has
    # {S : A subseteq S} of tiny or zero measure, so it can cut a finite pose list and can never
    # help a cover of the continuum.
    seg = poly = 0
    wseg = wpoly = 0.0
    areas = []
    for j in np.nonzero(ycov[npt:] > 1e-12)[0]:
        A = M.cuts[j]['obj'].A
        if len(A) == 2:
            seg += 1
            wseg += ycov[npt + j]
        else:
            poly += 1
            wpoly += ycov[npt + j]
            s = 0.0
            for i in range(len(A)):
                k = (i + 1) % len(A)
                s += A[i][0] * A[k][1] - A[k][0] * A[i][1]
            areas.append(abs(s) / 2)
    log(f"  of those anchor cliques: {seg} segment anchors (weight {wseg:.4f}), "
        f"{poly} polygon anchors (weight {wpoly:.4f})"
        + (f", polygon areas min {min(areas):.3e} med {sorted(areas)[len(areas)//2]:.3e} "
           f"max {max(areas):.3e}" if areas else ""))
    # (a) the LP's own columns: must all be covered to weight >= 1
    own = np.array(M.poses, dtype=float)
    cown = probe_cover(M, ycov, own, a.threads)
    log(f"  covered weight on the LP's OWN {len(own)} pose orbits: min {cown.min():.6f} "
        f"med {np.median(cown):.6f}")
    # (b) random admissible poses
    rng = np.random.default_rng(a.seed)
    n = a.probe_n
    th = rng.uniform(0, math.pi / 2, n)
    w2 = np.array([PD.wid_of(v) / 2 for v in th])
    cx = w2 + rng.uniform(0, 1, n) * (a.t - 2 * w2)
    cy = w2 + rng.uniform(0, 1, n) * (a.t - 2 * w2)
    P = np.c_[cx, cy, th]
    cr = probe_cover(M, ycov, P, a.threads)
    log(f"  covered weight on {n} RANDOM admissible poses: min {cr.min():.6f} "
        f"1% {np.quantile(cr, 0.01):.6f} med {np.median(cr):.6f} max {cr.max():.6f}")
    log(f"  fraction of random poses with covered weight < 1: {float((cr < 1).mean()):.4f}")
    # (c) perturbations of the LP's own support
    nb = []
    sup = [M.poses[i] for i in np.nonzero(mu > 1e-11)[0]]
    for (px, py, pt) in sup:
        for d in (0.001, 0.005, 0.02, 0.05):
            for (ax_, ay_, at_) in ((d, 0, 0), (-d, 0, 0), (0, d, 0), (0, -d, 0),
                                    (0, 0, math.radians(d * 50)), (0, 0, -math.radians(d * 50))):
                q = PD.clamp_pose(px + ax_, py + ay_, pt + at_, a.t)
                if q is not None:
                    nb.append(q)
    nb = np.array(nb, dtype=float)
    cn = probe_cover(M, ycov, nb, a.threads)
    log(f"  covered weight on {len(nb)} SUPPORT PERTURBATIONS: min {cn.min():.6f} "
        f"1% {np.quantile(cn, 0.01):.6f} med {np.median(cn):.6f}")
    log(f"  fraction < 1: {float((cn < 1).mean()):.4f}")


def cmd_dualsep(a, log):
    """Separate anchor cliques on a cover-LP dual (`runs/branch_TAG_dual.txt`) over the FULL
    family: any anchor point `p` (not only wall points) and any direction (not only the wall
    normal), with the cover side's own credit rule (`anchorclique.member`).

    `search/anchorsep.py` -- the separator task I's cover runs use -- scans only wall points
    `0.05 < dist < 1` with the wall-PERPENDICULAR Lemma-2 segment.  That is the family for which
    `K(p,A) ⊇ P_p` (Lemma 2), i.e. the family whose columns dominate a point column.  It is not
    the family that violates these packings hardest: the best violated cliques of the certified
    measures at `t = 3.99` sit in the INTERIOR, where `K(p,A)` is not a superset of `P_p` at all
    but "the squares through `p` that reach `A`, plus the squares that contain `A`".  Lemma 1
    does not forbid those -- it says a *superset* of an interior point clique is not a clique,
    and this is not a superset."""
    import anchorclique as AC
    from fractions import Fraction as Fr
    import anchorsep as AS
    R, y, s = AS.load_dual(a.dual)
    log(f"# {a.dual}: {len(R)} rows, dual mass {y.sum():.6f}, s = {s}")
    D = a.D
    sfr = Fr(int(round(s * D)), D)
    rot = AC._rot(R)
    ct, st, u0, u1 = rot
    h = R[:, 3]

    def ybar(cl):
        imgs = AC.images(sfr, cl)
        return float(AC.coeff(imgs, R, rot) @ y) / len(imgs), len(imgs)

    def point_clique(X, Y):
        return ((('P', Fr(X, D), Fr(Y, D)),), ((0, ()),))

    def anchor_clique(X, Y, X0, Y0, X1, Y1):
        return ((('P', Fr(X, D), Fr(Y, D)), ('S', Fr(X0, D), Fr(Y0, D), Fr(X1, D), Fr(Y1, D))),
                ((0, (1,)), (1, ())))

    # candidate anchor points: the grid points with the heaviest dual coverage (a point clique's
    # own ybar), over the WHOLE container, not just the wall band
    n = int(round(s / a.pitch)) + 1
    ax = np.linspace(0, s, n)
    G = np.array([(x, yv) for x in ax for yv in ax if x <= yv + 1e-12 and yv <= s / 2 + 1e-12])
    gs = [lambda p: p, lambda p: np.c_[s - p[:, 0], p[:, 1]], lambda p: np.c_[p[:, 0], s - p[:, 1]],
          lambda p: np.c_[s - p[:, 0], s - p[:, 1]], lambda p: np.c_[p[:, 1], p[:, 0]],
          lambda p: np.c_[s - p[:, 1], p[:, 0]], lambda p: np.c_[p[:, 1], s - p[:, 0]],
          lambda p: np.c_[s - p[:, 1], s - p[:, 0]]]
    cov = sum(AS.cover_grid(R, y, g(G)) for g in gs) / 8.0
    order = np.argsort(-cov)[:a.top]
    log(f"# {len(G)} grid points at pitch {a.pitch}; best point-clique ybar {cov.max():.6f}")
    dirs = [(math.cos(t), math.sin(t)) for t in
            np.linspace(0, 2 * math.pi, a.ndir, endpoint=False)]
    best = []
    t0 = time.time()
    for kk, i in enumerate(order):
        px, py = float(G[i][0]), float(G[i][1])
        X, Y = int(round(px * D)), int(round(py * D))
        yp = float(cov[i])
        for (dx, dy) in dirs:
            for ie in range(1, a.neps + 1):
                eps = a.epsmax * ie / a.neps
                mx, my = px + eps * dx, py + eps * dy
                for ir in range(1, a.nrho + 1):
                    rho = 0.49 * ir / a.nrho
                    x0, y0 = mx + rho * dy, my - rho * dx
                    x1, y1 = mx - rho * dy, my + rho * dx
                    if min(x0, x1, y0, y1) < 0 or max(x0, x1, y0, y1) > s:
                        continue
                    cl = anchor_clique(X, Y, int(round(x0 * D)), int(round(y0 * D)),
                                       int(round(x1 * D)), int(round(y1 * D)))
                    v, ni = ybar(cl)
                    if v > 1.0 + 1e-9:
                        best.append((v, yp, (px, py), (dx, dy), eps, rho, ni))
        if kk % 20 == 0:
            log(f"  {kk}/{len(order)}  best ybar(K) so far "
                f"{max((b[0] for b in best), default=0):.6f}  ({time.time()-t0:.0f}s)")
    best.sort(key=lambda z: -z[0])
    log(f"{'ybar(K)':>9} {'ybar(P_p)':>10} {'extra':>9}  {'p':>20} {'walldist':>9} "
        f"{'eps':>7} {'rho':>7} {'|A|':>7}")
    seen = set()
    shown = 0
    for (v, yp, p, d, eps, rho, ni) in best:
        k = (round(p[0], 6), round(p[1], 6))
        if k in seen:
            continue
        seen.add(k)
        wd = min(p[0], p[1], s - p[0], s - p[1])
        log(f"{v:9.5f} {yp:10.5f} {v-yp:9.5f}  ({p[0]:8.4f},{p[1]:8.4f}) {wd:9.4f} "
            f"{eps:7.4f} {rho:7.4f} {2*rho:7.4f}")
        shown += 1
        if shown >= a.show:
            break
    if best:
        v, yp, p, d, eps, rho, ni = best[0]
        wd = min(p[0], p[1], s - p[0], s - p[1])
        log(f"\n# BEST ybar(K) = {v:.6f} at p = ({p[0]:.4f}, {p[1]:.4f}), wall distance {wd:.4f}, "
            f"eps {eps:.4f}, |A| {2*rho:.4f}; the point clique there is {yp:.6f}")
        nw = sum(1 for b in best if min(b[2][0], b[2][1], s - b[2][0], s - b[2][1]) >= 1.0)
        log(f"# violated cliques found: {len(best)}, of which {nw} sit at wall distance >= 1 "
            f"(unreachable by search/anchorsep.py)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('build', 'price', 'escape', 'xmember', 'probe', 'dualsep'))
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
    ap.add_argument('--no-pure-rows', action='store_true',
                    help="converge rows for the clique solution only (task G's protocol), so the "
                         "matched pure value is the upper bound clique_continuum.py reports")
    ap.add_argument('--dual', default=None, help='dualsep: runs/branch_TAG_dual.txt')
    ap.add_argument('--D', type=int, default=784000, help='dualsep: coordinate denominator')
    ap.add_argument('--pitch', type=float, default=0.02, help='dualsep: anchor-point grid')
    ap.add_argument('--top', type=int, default=200, help='dualsep: anchor points to try')
    ap.add_argument('--ndir', type=int, default=16, help='dualsep: anchor directions')
    ap.add_argument('--neps', type=int, default=8, help='dualsep: anchor offsets')
    ap.add_argument('--nrho', type=int, default=8, help='dualsep: anchor half-lengths')
    ap.add_argument('--epsmax', type=float, default=0.9, help='dualsep: largest anchor offset')
    ap.add_argument('--seg-only', action='store_true',
                    help='keep and separate only SEGMENT anchors (the representable family)')
    ap.add_argument('--refine', action='store_true',
                    help='coordinate-descent the priced candidates on the true reduced cost')
    ap.add_argument('--n', type=int, default=200, help='xmember: number of cliques')
    ap.add_argument('--poses', type=int, default=20000, help='xmember: poses per clique')
    ap.add_argument('--Ns', default='2000,6000', help='xmember: verifier nets to compare')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--probe-n', type=int, default=200000, help='probe: random poses')
    ap.add_argument('--show', type=int, default=0)
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
    if a.cmd == 'xmember':
        cmd_xmember(a, log)
        return
    if a.cmd == 'probe':
        cmd_probe(a, log)
        return
    if a.cmd == 'dualsep':
        cmd_dualsep(a, log)
        return

    # ---- price / escape both need the built instance
    sargs = SepArgs()
    sargs.cut_per_round = a.cut_per_round
    if a.seg_only:
        segments_only(log)
    M = build_model(a.t, a.threads, log, a.support, a.cuts, row_pitch=a.row_pitch,
                    seg_only=a.seg_only)
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
        if a.refine:
            C0 = np.array(cands, dtype=float)
            C1, v1 = refine_true(M, lib, ax, ay, aw, ycut, C0, a.threads)
            log(f"  refine: best true rc {v1.max():.6f} on {len(C1)} candidates")
            cands = [(c[0], c[1], c[2]) for c in C1] + cands
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
