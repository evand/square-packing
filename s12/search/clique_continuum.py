#!/usr/bin/env python3
"""The clique-strengthened packing LP on a FIXED pose lattice  (task G, part 2).

    max  sum mu_S   s.t.   coverage(p) <= 1 at every point p of [0,t]^2,
                           sum_{S in K} mu_S <= 1 for every clique K,     mu >= 0

restricted to poses on a lattice: centres at pitch h (or exact, h = 0, via the arrangement
sweep pricer) and angles at pitch dtheta.  Two facts make the numbers readable:

  * a measure supported on the lattice that satisfies every clique of the lattice's own
    closed-intersection graph satisfies every clique constraint of the continuum, so its mass
    is a LOWER bound on the continuum clique-LP value V;
  * the same LP without clique cuts must reproduce the pure value L(3.99) = 12.008 -- the
    calibration.  A lattice that fails calibration says nothing.

Modes
    pure     no clique cuts (calibration)
    sat      cuts from the exact max-mass clique of the support's closed-intersection graph
             (the family task A used: any pairwise closed-intersecting set of support poses)
    kcut     cuts  mu(K(p)) <= 1  where K(p) = {S : S meets every admissible S' with p in S'}
             (the maximal point-anchored family; membership is geometric, so column generation
             cannot dodge it -- see search/clique_family.py and notes/clique-family.md)

Everything here is floating point and heuristic; the certification of a final measure is
exact and lives in dual_exact.py / clique_family.py.
"""
import sys, os, math, time, json, argparse
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
sys.path.insert(0, HERE)
import packing_dual as PD
import clique_family as CF

TOL = 1e-9


# ============================================================================ cliques (float)
# the closed-intersection graph and the max-weight-clique branch and bound live in
# clique_family.py (one implementation, used by both scripts)
closed_adj = CF.closed_adj
max_weight_clique = CF.max_weight_clique


# ============================================================================ K(p) cuts
class KCut:
    """the maximal point-anchored family K(p) as a geometric membership test.

    K_th'(p) = intersection over S in U(p) of (S (+) R_th' Q), with U(p) approximated from
    outside by its extreme members E (vertices of the centre polygons over an angle grid),
    so the family is a SUPERSET of the true K(p) -- the constraint is at least as strong.
    Membership of (c', th'):  <c', n> - w_n(th')/2 <= G(n) for every candidate normal n,
    G(n) = min over S in E of [ <c_S, n> + w_n(th_S)/2 ]."""

    def __init__(self, p, t, nth=180):
        self.p = p
        self.t = t
        E = CF.extreme_squares(p, t, nth=nth)
        self.E = E
        ths = np.unique(np.round(E[:, 2], 12))
        ct, st = np.cos(ths), np.sin(ths)
        nx = np.concatenate([ct, -ct, -st, st])
        ny = np.concatenate([st, -st, ct, -ct])
        self.nx, self.ny = nx, ny
        self.G = self._G(nx, ny)

    def _G(self, nx, ny):
        cS = self.E[:, :2]
        ctS, stS = np.cos(self.E[:, 2]), np.sin(self.E[:, 2])
        d1 = np.abs(np.outer(nx, ctS) + np.outer(ny, stS))
        d2 = np.abs(np.outer(-nx, stS) + np.outer(ny, ctS))
        proj = np.outer(nx, cS[:, 0]) + np.outer(ny, cS[:, 1])
        return np.min(proj + 0.5 * (d1 + d2), axis=1)

    def members(self, poses, tol=1e-9, chunk=200000):
        """boolean mask: which of the poses (n, 3) lie in K(p)"""
        poses = np.asarray(poses, dtype=float).reshape(-1, 3)
        out = np.zeros(len(poses), dtype=bool)
        for a in range(0, len(poses), chunk):
            b = min(a + chunk, len(poses))
            cx, cy, th = poses[a:b, 0], poses[a:b, 1], poses[a:b, 2]
            ct, st = np.cos(th), np.sin(th)
            # cheap reject: K(p) is inside a ball of radius ~sqrt(2) around p
            ok = (np.abs(cx - self.p[0]) < 1.5) & (np.abs(cy - self.p[1]) < 1.5)
            idx = np.nonzero(ok)[0]
            if len(idx) == 0:
                continue
            cxs, cys, cts, sts = cx[idx], cy[idx], ct[idx], st[idx]
            m = np.ones(len(idx), dtype=bool)
            # fixed normals
            wn = 0.5 * (np.abs(np.outer(self.nx, cts) + np.outer(self.ny, sts))
                        + np.abs(np.outer(-self.nx, sts) + np.outer(self.ny, cts)))
            lhs = np.outer(self.nx, cxs) + np.outer(self.ny, cys) - wn
            m &= np.all(lhs <= self.G[:, None] + tol, axis=0)
            # the pose's own normals
            for (ex, ey) in ((cts, sts), (-cts, sts), (-sts, cts), (sts, -cts)):
                g = self._G(ex, ey)
                m &= (cxs * ex + cys * ey - 0.5) <= g + tol
            out[a + idx] = m
        return out


class ACut:
    """the anchor clique  K(p, A) = {S : p in S, S meets A} u {S : A subset S}.

    A clique for EVERY point p and every nonempty compact convex A, and membership is a
    geometric test, so newly generated columns are charged correctly."""

    def __init__(self, p, A):
        self.p = p
        self.A = [tuple(v) for v in A]

    def members(self, poses):
        poses = np.asarray(poses, dtype=float).reshape(-1, 3)
        IM = np.c_[poses, np.zeros(len(poses))]
        thru = CF.contains_np(IM, [self.p])
        return (thru & CF.meets_np(IM, self.A)) | ((~thru) & CF.contains_np(IM, self.A))


# ============================================================================ LP model
class CModel:
    """packing LP: orbit columns, point rows in the fundamental domain, plus cut rows"""

    def __init__(self, t, threads, log, kmass=None, r=1.0):
        self.m = PD.Model(t, PD.build_lib(), threads, log)
        self.t = t
        self.kmass = kmass
        self.r = r
        self.lam = 0.0
        self.cuts = []                     # list of dicts: kind, data, matrix column values
        self.Ccols = []                    # list of 1-d arrays (len = ncols) per cut
        self.log = log

    @property
    def poses(self):
        return self.m.poses

    def add_poses(self, cand):
        n0 = len(self.m.poses)
        k = self.m.add_poses(cand)
        if k and self.cuts:
            new = np.array(self.m.poses[n0:], dtype=float)
            for i, cut in enumerate(self.cuts):
                self.Ccols[i] = np.concatenate([self.Ccols[i], self.cut_row(cut, new)])
        return k

    def add_points(self, pts):
        return self.m.add_points(pts)

    def cut_row(self, cut, poses):
        """coefficient of each orbit in the cut: (#images in the clique)/8"""
        poses = np.asarray(poses, dtype=float).reshape(-1, 3)
        if len(poses) == 0:
            return np.zeros(0)
        imgs = []
        for (cx, cy, th) in poses:
            imgs += PD.pose_images(cx, cy, th, self.t)
        imgs = np.array(imgs, dtype=float).reshape(len(poses), 8, 3)
        flat = imgs.reshape(-1, 3)
        if cut['kind'] in ('kcut', 'anchor'):
            mask = cut['obj'].members(flat)
        else:
            mask = self.sat_members(cut, flat)
        return mask.reshape(len(poses), 8).sum(axis=1) / 8.0

    def sat_members(self, cut, flat):
        """membership in an explicit clique: the pose must equal one of the listed images"""
        key = cut['keys']
        k = np.round(flat[:, 0] * 1e7).astype(np.int64) * 1000000000000 \
            + np.round(flat[:, 1] * 1e7).astype(np.int64) * 1000000 \
            + np.round(np.degrees(flat[:, 2]) * 1e4).astype(np.int64)
        return np.isin(k, key)

    def add_cut(self, cut):
        poses = np.array(self.m.poses, dtype=float)
        row = self.cut_row(cut, poses)
        self.cuts.append(cut)
        self.Ccols.append(row)

    def drop_cuts(self, keep):
        self.cuts = [c for c, k in zip(self.cuts, keep) if k]
        self.Ccols = [c for c, k in zip(self.Ccols, keep) if k]

    def flags(self):
        """orbits whose poses are centred in one of the four corner boxes [0,r]^2 (leaf mode)"""
        t = self.t
        P = np.array(self.m.poses, dtype=float).reshape(-1, 3)
        if len(P) == 0:
            return np.zeros(0, dtype=bool)
        return (np.minimum(P[:, 0], t - P[:, 0]) <= self.r + 1e-12) & \
               (np.minimum(P[:, 1], t - P[:, 1]) <= self.r + 1e-12)

    def solve(self):
        n = len(self.m.poses)
        A = self.m.A
        b = np.ones(A.shape[0])
        if self.cuts:
            C = sp.csr_matrix(np.array(self.Ccols))
            A = sp.vstack([A, C], format='csr')
            b = np.ones(A.shape[0])
        Aeq = beq = None
        if self.kmass is not None:
            Aeq = sp.csr_matrix(self.flags().astype(float).reshape(1, -1))
            beq = [float(self.kmass)]
        res = linprog(c=-np.ones(n), A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq, bounds=(0, None),
                      method='highs')
        if not res.success:
            res = linprog(c=-np.ones(n), A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq, bounds=(0, None),
                          method='highs-ds')
        if not res.success:
            return None
        mu = np.maximum(res.x, 0.0)
        y = np.maximum(-res.ineqlin.marginals, 0.0)
        self.lam = float(-res.eqlin.marginals[0]) if Aeq is not None else 0.0
        npt = self.m.A.shape[0]
        return mu, y[:npt], y[npt:], -res.fun


# ============================================================================ lattice pricing
def lattice_angles(dth_deg, thmax=90.0):
    return np.radians(np.arange(0.0, thmax - 1e-9, dth_deg))


def lattice_slice(t, h, th):
    """the admissible centre lattice at one angle (generated on demand: memory)"""
    w2 = PD.wid_of(th) / 2
    lo, hi = w2 + 1e-9, t - w2 - 1e-9
    if hi <= lo:
        return np.zeros((0, 3))
    n = int(math.floor((hi - lo) / h))
    g = lo + h * np.arange(n + 1)
    g = np.unique(np.concatenate([g, [hi]]))
    G0, G1 = np.meshgrid(g, g, indexing='ij')
    return np.ascontiguousarray(np.c_[G0.ravel(), G1.ravel(), np.full(G0.size, th)])


def lattice_size(t, h, angs):
    return sum(len(lattice_slice(t, h, th)) for th in angs)


def price_lattice(lib, ax, ay, aw, t, h, angs, cutduals, cuts, threads, want,
                  lam=0.0, rr=1.0):
    tt = t
    """reduced cost of every lattice pose = 1 - capture - sum of duals of the cuts it is in.
    Returns the `want` best (cx, cy, th, rc)."""
    best = []
    for th in angs:
        P = lattice_slice(t, h, th)
        if len(P) == 0:
            continue
        cap = PD.capture(lib, ax, ay, aw, P, threads)
        rc = 1.0 - cap
        if lam:
            corner = (np.minimum(P[:, 0], tt - P[:, 0]) <= rr + 1e-12) & \
                     (np.minimum(P[:, 1], tt - P[:, 1]) <= rr + 1e-12)
            rc = rc + lam * corner
        for (cut, dual) in zip(cuts, cutduals):
            if dual <= 1e-12:
                continue
            if cut['kind'] in ('kcut', 'anchor'):
                rc = rc - dual * cut['obj'].members(P)
        idx = np.nonzero(rc > 1e-9)[0]
        if len(idx):
            if len(idx) > want:
                idx = idx[np.argsort(-rc[idx])[:want]]
            best.append(np.c_[P[idx], rc[idx]])
    if not best:
        return np.zeros((0, 4))
    B = np.vstack(best)
    B = B[np.argsort(-B[:, 3])]
    return B[:want]


# ============================================================================ main loop
def read_support_any(path, t_from=None, t=None):
    poses = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        q = line.split()
        if q[0] != 'pose':
            continue
        if len(q) == 6:                       # exact: p q cx cy mass
            from fractions import Fraction as Fr
            th = 2 * math.atan2(int(q[1]), int(q[2]))
            poses.append((float(Fr(q[3])), float(Fr(q[4])), th))
        else:
            poses.append((float(q[1]), float(q[2]), math.radians(float(q[3]))))
    return poses


def snap_to_lattice(poses, t, h, dth_deg):
    out = []
    for (cx, cy, th) in poses:
        thd = round(math.degrees(th) / dth_deg) * dth_deg
        thd = thd % 90.0
        th2 = math.radians(thd)
        w2 = PD.wid_of(th2) / 2
        lo, hi = w2 + 1e-9, t - w2 - 1e-9
        if hi <= lo:
            continue
        f = lambda z: min(max(lo + round((z - lo) / h) * h, lo), hi)
        out.append((f(cx), f(cy), th2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('T', type=str)
    ap.add_argument('TAG')
    ap.add_argument('--mode', default='pure', choices=('pure', 'sat', 'kcut', 'anchor'))
    ap.add_argument('--sep-pts', type=int, default=40)
    ap.add_argument('--nseed', type=int, default=4)
    ap.add_argument('--maxstep', type=int, default=12)
    ap.add_argument('--h', type=float, default=0.01, help='centre pitch (0 = exact sweep pricer)')
    ap.add_argument('--dth', type=float, default=1.0, help='angle pitch in degrees')
    ap.add_argument('--rounds', type=int, default=60)
    ap.add_argument('--time', type=float, default=7200)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--warm', action='append', default=[])
    ap.add_argument('--row-pitch', type=float, default=0.05)
    ap.add_argument('--cg-want', type=int, default=400)
    ap.add_argument('--cut-per-round', type=int, default=12)
    ap.add_argument('--cut-inner', type=int, default=40, help='cut-only iterations per column round')
    ap.add_argument('--clique-tl', type=float, default=20.0)
    ap.add_argument('--maxatoms', type=int, default=8000)
    ap.add_argument('--kcut-nth', type=int, default=120)
    ap.add_argument('--kmass', type=float, default=None,
                    help='leaf mode: total mass of poses centred in the four corner boxes [0,r]^2')
    ap.add_argument('--r', type=float, default=1.0)
    args = ap.parse_args()

    t = float(eval(args.T)) if '/' in args.T else float(args.T)
    logf = open(os.path.join(RUNS, f'cq_{args.TAG}.log'), 'w')

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + '\n')
        logf.flush()

    log(f"# clique_continuum {args.TAG}: t={t} mode={args.mode} h={args.h} dth={args.dth} "
        f"threads={args.threads}")
    M = CModel(t, args.threads, log, kmass=args.kmass, r=args.r)
    lib = M.m.lib

    # ---- lattice
    angs = lattice_angles(args.dth)
    if args.h > 0:
        log(f"# lattice: {lattice_size(t, args.h, angs)} poses "
            f"(h={args.h}, dth={args.dth} deg, {len(angs)} angles)")
    else:
        log(f"# exact-centre sweep pricer over {len(angs)} angles (dth={args.dth} deg)")

    # ---- seeds
    seeds = []
    for wpath in args.warm:
        ps = read_support_any(wpath)
        seeds += ps
        log(f"# warm start {wpath}: {len(ps)} poses")
    if args.h > 0 and seeds:
        seeds = snap_to_lattice(seeds, t, args.h, args.dth)
    seeds += PD.seed_poses(t, 0.25, max(args.dth, 5.0))
    if args.kmass is not None:
        for thd in np.arange(0.0, 90.0 - 1e-9, max(args.dth, 2.0)):
            th = math.radians(thd)
            w2 = PD.wid_of(th) / 2
            for a in np.arange(w2, min(args.r, t / 2) + 1e-9, 0.02):
                for b2 in np.arange(w2, min(args.r, t / 2) + 1e-9, 0.02):
                    seeds.append((float(a), float(b2), th))
    if args.h > 0:
        seeds = snap_to_lattice(seeds, t, args.h, args.dth)
    M.add_poses(seeds)
    M.add_points(PD.seed_points(t, args.row_pitch))
    log(f"# start: {len(M.poses)} orbits, {M.m.A.shape[0]} rows")

    t0 = time.time()
    hist = []
    for rnd in range(args.rounds):
        if time.time() - t0 > args.time:
            log("# time limit")
            break
        # ---------- inner loop: rows + cuts on a fixed pose set
        for it in range(args.cut_inner):
            r = M.solve()
            if r is None:
                log("  LP FAILED")
                break
            mu, y, ycut, obj = r
            # row separation: arrangement vertices with coverage > 1
            Mx, bad, nv, ovf = M.m.certify(mu, thresh=1.0 + 1e-7, cap=200000)
            nadd = 0
            if len(bad):
                nadd = M.add_points(bad[:, :2])
            # clique separation
            ncut = 0
            mcl = 0.0
            if args.mode != 'pure':
                ncut, mcl = separate(M, mu, args, log)
            if it % 10 == 0 or (nadd == 0 and ncut == 0):
                log(f"  r{rnd} it{it}: obj {obj:.6f} cov {Mx:.6f} rows {M.m.A.shape[0]} "
                    f"cols {len(M.poses)} cuts {len(M.cuts)} +rows {nadd} +cuts {ncut} "
                    f"maxclique {mcl:.4f}")
            if nadd == 0 and ncut == 0:
                break
            if time.time() - t0 > args.time:
                break
        # ---------- pricing over the lattice
        r = M.solve()
        if r is None:
            break
        mu, y, ycut, obj = r
        Mx, bad, nv, ovf = M.m.certify(mu, thresh=1.0 + 1e-7, cap=200000)
        mcl = residual_clique(M, mu, args) if args.mode != 'pure' else 0.0
        hist.append(dict(round=rnd, obj=obj, cov=Mx, rows=int(M.m.A.shape[0]),
                         cols=len(M.poses), cuts=len(M.cuts), maxclique=mcl))
        log(f"* r{rnd}: STAGE obj {obj:.6f}  maxcov {Mx:.8f}  rows {M.m.A.shape[0]} "
            f"cols {len(M.poses)} cuts {len(M.cuts)} maxclique {mcl:.5f}  "
            f"({time.time()-t0:.0f}s)")
        write_support(M, mu, obj, Mx, args.TAG, t)
        json.dump(dict(tag=args.TAG, t=t, mode=args.mode, h=args.h, dth=args.dth, obj=obj,
                       maxcov=Mx, hist=hist),
                  open(os.path.join(RUNS, f'cq_{args.TAG}.json'), 'w'), indent=1)
        ax, ay, aw = PD.atoms_from_dual(M.m, y, maxrows=args.maxatoms)
        if len(ax) == 0:
            log("  no dual atoms")
            break
        if args.h > 0:
            cand = price_lattice(lib, ax, ay, aw, t, args.h, angs, ycut, M.cuts,
                                 args.threads, args.cg_want, lam=M.lam, rr=args.r)
            if len(cand) == 0:
                log("* PRICING CLEAN: LP solved to optimality over the lattice")
                break
            log(f"  pricing: best reduced cost {cand[0, 3]:.6f}, adding {len(cand)}")
            nnew = M.add_poses([(c[0], c[1], c[2]) for c in cand])
        else:
            got = PD.price_sweep(M.m, lib, ax, ay, aw, angs, 6, 0.02, args.threads)
            got = [g for g in got if g[3] < 1.0 - 1e-9]
            if not got:
                log("* PRICING CLEAN (sweep): optimal over all centres at these angles")
                break
            got.sort(key=lambda g: g[3])
            log(f"  sweep pricing: best capture {got[0][3]:.6f}, adding {min(len(got), args.cg_want)}")
            nnew = M.add_poses([(g[0], g[1], g[2]) for g in got[:args.cg_want]])
        if nnew == 0:
            log("* no new columns (all duplicates) -- stop")
            break

    r = M.solve()
    mu, y, ycut, obj = r
    Mx, bad, nv, ovf = M.m.certify(mu, thresh=1.0 + 1e-7, cap=200000)
    mcl = residual_clique(M, mu, args) if args.mode != 'pure' else 0.0
    log(f"# FINAL {args.TAG}: obj {obj:.6f}  maxcov {Mx:.9f}  L = obj/max(1,cov) "
        f"{obj/max(1.0, Mx):.6f}  rows {M.m.A.shape[0]} cols {len(M.poses)} cuts {len(M.cuts)} "
        f"residual max clique {mcl:.5f}")
    write_support(M, mu, obj, Mx, args.TAG, t)
    json.dump(dict(tag=args.TAG, t=t, mode=args.mode, h=args.h, dth=args.dth, obj=obj,
                   maxcov=Mx, L=obj / max(1.0, Mx), rows=int(M.m.A.shape[0]),
                   cols=len(M.poses), cuts=len(M.cuts), maxclique=mcl, hist=hist),
              open(os.path.join(RUNS, f'cq_{args.TAG}.json'), 'w'), indent=1)


def support_images(M, mu, thr=1e-11):
    sup = np.nonzero(mu > thr)[0]
    rows = []
    for i in sup:
        cx, cy, th = M.poses[i]
        for (x, y, u) in PD.pose_images(cx, cy, th, M.t):
            rows.append((x, y, u, mu[i] / 8.0))
    return np.array(rows) if rows else np.zeros((0, 4))


def residual_clique(M, mu, args):
    im = support_images(M, mu)
    if len(im) == 0 or len(im) > 4000:
        im = im[np.argsort(-im[:, 3])[:4000]] if len(im) else im
    if len(im) == 0:
        return 0.0
    adj = closed_adj(im[:, :3])
    w, cl = max_weight_clique(adj, im[:, 3], tlimit=args.clique_tl)
    return w


def separate(M, mu, args, log):
    """add violated clique cuts; returns (#added, max clique mass found)"""
    im = support_images(M, mu)
    if len(im) == 0:
        return 0, 0.0
    if len(im) > 3000:
        im = im[np.argsort(-im[:, 3])[:3000]]
    if args.mode == 'sat':
        adj = closed_adj(im[:, :3])
        w0, cl0 = max_weight_clique(adj, im[:, 3], tlimit=args.clique_tl)
        cands = [(w0, cl0)]
        # plus a greedy clique through each of the heaviest images (task A's separation)
        order = np.argsort(-im[:, 3])[:args.cut_per_round * 3]
        for v in order:
            ww, cc = CF.greedy_clique(adj, im[:, 3], [int(v)])
            if ww > 1.0 + 1e-6:
                cands.append((ww, cc))
        cands.sort(key=lambda z: -z[0])
        seen = set()
        added = 0
        for (ww, cc) in cands:
            if ww <= 1.0 + 1e-6 or added >= args.cut_per_round:
                break
            key = tuple(sorted(cc))
            if key in seen:
                continue
            seen.add(key)
            keys = np.round(im[cc, 0] * 1e7).astype(np.int64) * 1000000000000 \
                + np.round(im[cc, 1] * 1e7).astype(np.int64) * 1000000 \
                + np.round(np.degrees(im[cc, 2]) * 1e4).astype(np.int64)
            M.add_cut(dict(kind='sat', keys=keys))
            added += 1
        return added, w0
    if args.mode == 'anchor':
        pts = tight_row_points(M, args, npts=args.sep_pts)
        recs = []
        for p in pts:
            r = CF.anchor_separate(im, p, nseed=args.nseed, maxstep=args.maxstep)
            if r is not None:
                recs.append(r)
        recs.sort(key=lambda r: -r['mass'])
        added = 0
        worst = recs[0]['mass'] if recs else 0.0
        for r in recs:
            if r['mass'] <= 1.0 + 1e-6 or added >= args.cut_per_round:
                break
            if any(c['kind'] == 'anchor' and abs(c['obj'].p[0] - r['p'][0]) < 1e-9
                   and abs(c['obj'].p[1] - r['p'][1]) < 1e-9
                   and len(c['obj'].A) == len(r['A'])
                   and max(abs(u[0] - v[0]) + abs(u[1] - v[1])
                           for u, v in zip(c['obj'].A, r['A'])) < 1e-9 for c in M.cuts):
                continue
            M.add_cut(dict(kind='anchor', obj=ACut(r['p'], r['A'])))
            added += 1
        return added, worst
    # ---- kcut: anchor points where mu(K(p)) > 1
    cand = kcut_candidates(M, im, args)
    added = 0
    worst = 0.0
    for (val, p) in cand:
        worst = max(worst, val)
        if val <= 1.0 + 1e-6 or added >= args.cut_per_round:
            continue
        if any(c['kind'] == 'kcut' and abs(c['p'][0] - p[0]) < 1e-6 and abs(c['p'][1] - p[1]) < 1e-6
               for c in M.cuts):
            continue
        M.add_cut(dict(kind='kcut', p=p, obj=KCut(p, M.t, nth=args.kcut_nth)))
        added += 1
    return added, worst


def tight_row_points(M, args, npts=60):
    """the row points with the highest coverage (candidate anchor points)"""
    r = M.solve()
    if r is None:
        return []
    mu = r[0]
    cov = M.m.A.dot(mu)
    pts = M.m.pts
    t = M.t
    sel = np.nonzero(cov > 0.95)[0]
    if len(sel) == 0:
        return []
    sel = sel[np.argsort(-cov[sel])]
    if len(sel) > npts:
        sel = sel[::max(1, len(sel) // npts)][:npts]
    out = []
    for i in sel:
        x, y = float(pts[i, 0]), float(pts[i, 1])
        for (a, b) in ((x, y), (y, x), (t - x, y), (x, t - y)):
            out.append((a, b))
    return out[:4 * npts]


def kcut_candidates(M, im, args, npts=60):
    """points p with the largest mu(K(p)):  candidates are the tight rows near a wall"""
    t = M.t
    pts = M.m.pts
    # coverage at each row point
    A = M.m.A
    r = M.solve()
    if r is None:
        return []
    mu = r[0]
    cov = A.dot(mu)
    near = np.minimum(np.minimum(pts[:, 0], pts[:, 1]), np.minimum(t - pts[:, 0], t - pts[:, 1]))
    sel = np.nonzero((cov > 0.97) & (near < 1.0))[0]
    if len(sel) == 0:
        return []
    sel = sel[np.argsort(-cov[sel])[:npts]]
    out = []
    for i in sel:
        p = (float(pts[i, 0]), float(pts[i, 1]))
        try:
            kc = KCut(p, t, nth=args.kcut_nth)
        except Exception:
            continue
        mask = kc.members(im[:, :3])
        out.append((float(im[mask, 3].sum()), p))
    out.sort(reverse=True)
    return out


def write_support(M, mu, obj, Mx, tag, t):
    path = os.path.join(RUNS, f'cq_{tag}_support.txt')
    with open(path, 'w') as f:
        f.write(f"# t={t} mass={obj:.9f} M={Mx:.9f}\n")
        f.write("# pose cx cy theta_deg mu (mass mu/8 on each of the 8 dihedral images)\n")
        for i, (cx, cy, th) in enumerate(M.poses):
            if mu[i] > 1e-12:
                f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} {mu[i]:.12f}\n")
    return path


if __name__ == '__main__':
    main()
