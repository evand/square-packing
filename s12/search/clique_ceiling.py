#!/usr/bin/env python3
"""clique_ceiling.py -- the ceiling of the clique-strengthened method, computed on the packing side.

    max  sum_S mu_S
    s.t. coverage(p)      = sum_S mu_S [p in S]        <= 1   for every point p of the container,
         clique mass(K)   = sum_{S in K} mu_S          <= 1   for every clique K of the CLOSED-intersection
                                                              graph of poses (pairwise closed squares meet),
         [leaf mode]      sum_{S centred in a corner box} mu_S = k,
         mu >= 0,

over closed unit squares S inside [0,t]^2 (D4-symmetrised: mass mu/8 on each dihedral image).  Column
generation on the poses (packing_dual.py's pricers), row generation on the points (violated
arrangement vertices) and on the cliques (max-mass clique of the current measure, plus the best clique
through each heavy pose; existing cliques are LIFTED by new poses that meet all their members).

Why closed intersection: a certificate at t refutes packings at side t' < t; rescaled to [0,t]^2 those
are pairwise DISJOINT closed unit squares, so any set of poses that pairwise closed-intersect (touching
counts) holds at most one square of a packing.  A feasible measure must therefore respect every such
clique; strict-overlap cliques are a subfamily and a measure feasible only for them is not a ceiling.

By weak duality, a measure feasible for all these constraints has mass <= the value of ANY clique
certificate (points, box cliques, or any pairwise-closed-intersecting pose family) at that t, so a
certified mass >= 12 at t says: no clique certificate below 12 exists at any t' >= t.

Certification (--exact): poses snapped to rational rotations 2 arctan(p/q) and rational centres,
exactly admissible; coverage at every arrangement vertex of the support exactly (dual_exact.py);
closed intersection of every pair of images decided exactly in integers (separating-axis test on
integer corner coordinates, <=); maximum clique mass by branch and bound on INTEGER masses; masses
rounded down and, if needed, scaled so that both maxima are <= 1 exactly.  Floats choose the measure.

    python3 search/clique_ceiling.py 3.99 C99 --warm /path/runs/dual_PA2_support.txt --rounds 40 --exact
    python3 search/clique_ceiling.py 3.98 L98 --warm ... --branch-dual /path/runs/branch_t398hk4_dual_it16.txt --kmass 4 --r 1 --exact
"""
import sys, os, math, time, json, argparse
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction as Fr
from math import gcd, lcm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import packing_dual as PD
import clique_check as CC
import dual_exact as DE

RUNS = PD.RUNS
EPS_CLOSED = 1e-9          # float SAT: closed intersection with this tolerance (cut CHOICE only)


# ============================================================================== float geometry
def sat_closed(P, Q, eps=EPS_CLOSED):
    """bool |P| x |Q|: closed unit squares (cx, cy, th) pairwise intersect (separating-axis test, <=).
    Over-inclusive by eps (a pair that barely misses may be called intersecting); every clique chosen
    with it is re-checked exactly at certification, so this can only cost mass, never soundness."""
    P = np.asarray(P, float).reshape(-1, 3); Q = np.asarray(Q, float).reshape(-1, 3)
    if len(P) == 0 or len(Q) == 0: return np.zeros((len(P), len(Q)), bool)
    dx = Q[None, :, 0] - P[:, None, 0]; dy = Q[None, :, 1] - P[:, None, 1]
    ok = (dx * dx + dy * dy) <= 2.0 + 1e-6            # centre distance <= sqrt2 is necessary
    for src in (0, 1):
        th = P[:, None, 2] if src == 0 else Q[None, :, 2]
        for k in range(2):
            a = th + k * math.pi / 2; ux = np.cos(a); uy = np.sin(a)
            d = np.abs(dx * ux + dy * uy)
            eP = 0.5 * (np.abs(np.cos(P[:, None, 2] - a)) + np.abs(np.sin(P[:, None, 2] - a)))
            eQ = 0.5 * (np.abs(np.cos(Q[None, :, 2] - a)) + np.abs(np.sin(Q[None, :, 2] - a)))
            ok &= d <= eP + eQ + eps
    return ok


def expand_images(poses, mu, t):
    """all 8 images (cx, cy, th) of each pose, parent index, image index g, mass mu/8"""
    X = []; par = []; gs = []
    for k, (cx, cy, th) in enumerate(poses):
        for g, (x, y, u) in enumerate(PD.pose_images(cx, cy, th, t)):
            X.append((x, y, u)); par.append(k); gs.append(g)
    return np.array(X, float).reshape(-1, 3), np.array(par, int), np.array(gs, int)


def in_box(cx, cy, th, t, r):
    """centre of some image in a corner box [0,r]^2 (equivalently: both coordinates within r of a wall)"""
    return min(cx, t - cx) <= r + 1e-12 and min(cy, t - cy) <= r + 1e-12


# ============================================================================== box cliques
SQ2 = math.sqrt(2.0)


def sat_margin(X):
    """m x m: min over the 4 edge normals of (e_i + e_j - |d.u|) for unit squares X (cx, cy, th);
    >= 0 iff the closed squares intersect; the amount by which both may shrink (in half-extent) and still meet"""
    X = np.asarray(X, float).reshape(-1, 3); n = len(X)
    dx = X[None, :, 0] - X[:, None, 0]; dy = X[None, :, 1] - X[:, None, 1]
    mg = np.full((n, n), np.inf)
    for src in (0, 1):
        th = X[:, None, 2] if src == 0 else X[None, :, 2]
        for k in range(2):
            a = th + k * math.pi / 2; ux = np.cos(a); uy = np.sin(a)
            d = np.abs(dx * ux + dy * uy)
            eP = 0.5 * (np.abs(np.cos(X[:, None, 2] - a)) + np.abs(np.sin(X[:, None, 2] - a)))
            eQ = 0.5 * (np.abs(np.cos(X[None, :, 2] - a)) + np.abs(np.sin(X[None, :, 2] - a)))
            mg = np.minimum(mg, eP + eQ - d)
    return mg


def boxes_from_clique(X, w, dmin=2e-4, safety=0.9):
    """member images X (m x 3) with masses w -> pose boxes (cx, cy, th, h, phi) such that the cores of
    any two boxes intersect (so every pose in the union of boxes pairwise meets every other: a clique).
    Core of a box: the square at the box centre pose shrunk to side 1 - delta_i, where
    delta_i = phi_i + 2*sqrt2*h_i (sigma(phi) >= 1 - phi and the centre-box erosion is <= sqrt2*h per side);
    two cores meet if margin_ij >= (delta_i + delta_j)/sqrt2.  Members whose delta would be < dmin are
    dropped (a sub-clique is a clique) while the remaining mass stays > 1."""
    X = np.asarray(X, float).reshape(-1, 3); w = np.asarray(w, float); idx = np.arange(len(X))
    while len(idx) > 1:
        mg = sat_margin(X[idx]); np.fill_diagonal(mg, np.inf)
        delta = safety * SQ2 * mg.min(axis=1) / 2.0        # equal split: delta_i + delta_j <= sqrt2 * margin_ij
        delta = np.minimum(delta, 0.05)
        if delta.min() >= dmin or w[idx].sum() - w[idx[np.argmin(delta)]] <= 1.0: break
        idx = np.delete(idx, np.argmin(delta))
    if len(idx) <= 1: return None
    keep = delta >= 2e-9
    if keep.sum() <= 1: return None
    idx = idx[keep]; delta = delta[keep]
    phi = np.maximum(delta / 2.0, 1e-9); h = np.maximum(delta / (4.0 * SQ2), 1e-9)
    return np.c_[X[idx], h, phi], idx


def in_boxes(X, B, chunk=4000):
    """bool n: image (cx, cy, th) lies in some box (cx, cy, th, h, phi) (angles mod pi/2)"""
    X = np.asarray(X, float).reshape(-1, 3); out = np.zeros(len(X), bool)
    if len(B) == 0: return out
    for s0 in range(0, len(X), chunk):
        Xc = X[s0:s0 + chunk]
        ok = np.abs(Xc[:, None, 0] - B[None, :, 0]) <= B[None, :, 3]
        ok &= np.abs(Xc[:, None, 1] - B[None, :, 1]) <= B[None, :, 3]
        da = (Xc[:, None, 2] - B[None, :, 2] + math.pi / 4) % (math.pi / 2) - math.pi / 4
        ok &= np.abs(da) <= B[None, :, 4]
        out[s0:s0 + chunk] = ok.any(axis=1)
    return out


def greedy_clique(adj, w, v):
    """clique through v: add neighbours by decreasing weight while pairwise adjacent"""
    cand = np.nonzero(adj[v])[0]; cand = cand[np.argsort(-w[cand])]
    K = [v]; mask = np.ones(len(cand), bool)
    for i, u in enumerate(cand):
        if not mask[i]: continue
        K.append(u); mask &= adj[u][cand] | (np.arange(len(cand)) <= i)
    return K


def find_cliques(adj, w, time_limit, per_round, per_time, greedy_top, thresh=1 + 1e-7):
    """(max clique mass, its size, seconds, list of violated cliques as index lists)"""
    bw, bc, nodes, dt = CC.max_weight_clique(adj, w, time_limit=time_limit)
    found = []; seen = set()
    def push(val, mem):
        key = tuple(sorted(int(i) for i in mem))
        if val > thresh and key not in seen: seen.add(key); found.append(list(key))
    push(bw, bc)
    order = np.argsort(-w)
    for v in order[:greedy_top]:
        K = greedy_clique(adj, w, v); push(w[K].sum(), K)
    t1 = time.time()
    for v in order[:per_round]:
        if time.time() - t1 > per_time * per_round: break
        nb = np.nonzero(adj[v])[0]
        if nb.size == 0: continue
        bw2, bc2, _, _ = CC.max_weight_clique(adj[np.ix_(nb, nb)], w[nb], time_limit=per_time)
        push(bw2 + w[v], [v] + [nb[i] for i in bc2])
    return bw, len(bc), dt, found


# ============================================================================== model with box-clique cuts
class CModel(PD.Model):
    def __init__(self, t, lib, threads, log, kmass=None, r=1.0, margin=0.0):
        super().__init__(t, lib, threads, log)
        self.cuts = []            # list of box arrays (m x 5): union of pose boxes whose cores pairwise meet
        self.K = sp.csr_matrix((0, 0))   # cut rows x columns, coefficient = (#images of the column in the union)/8
        self.kmass = kmass; self.r = r; self.margin = margin
        self.z = np.zeros(0); self.lam = 0.0; self.greedy_top = 300

    def flags(self):
        return np.array([in_box(cx, cy, th, self.t, self.r) for (cx, cy, th) in self.poses], bool)

    def cut_rows(self, cuts, poses):
        """csr len(cuts) x len(poses)"""
        if not cuts or not poses: return sp.csr_matrix((len(cuts), len(poses)))
        X, par, gs = expand_images(poses, None, self.t)
        rows = []; cols = []; vals = []
        for i, B in enumerate(cuts):
            hit = in_boxes(X, B)
            if hit.any():
                c = np.bincount(par[hit], minlength=len(poses))
                nz = np.nonzero(c)[0]; rows += [i] * len(nz); cols += nz.tolist(); vals += (c[nz] / 8.0).tolist()
        return sp.csr_matrix((vals, (rows, cols)), shape=(len(cuts), len(poses)))

    def add_poses(self, cand):
        n0 = len(self.poses); k = super().add_poses(cand)
        if k and self.cuts:
            Knew = self.cut_rows(self.cuts, self.poses[n0:])
            self.K = sp.hstack([self.K, Knew], format='csr') if self.K.shape[1] else Knew
        elif self.cuts and self.K.shape[1] != len(self.poses):
            self.K = self.cut_rows(self.cuts, self.poses)
        return k

    def add_cuts(self, cuts):
        if not cuts: return 0
        Knew = self.cut_rows(cuts, self.poses)
        self.K = sp.vstack([self.K, Knew], format='csr') if self.K.shape[0] else Knew
        self.cuts += cuts
        return len(cuts)

    def solve(self):
        n = len(self.poses); m = len(self.pts)
        A = self.A if not self.cuts else sp.vstack([self.A, self.K], format='csr')
        b = np.full(A.shape[0], 1.0 - self.margin)
        Aeq = beq = None
        if self.kmass is not None:
            Aeq = sp.csr_matrix(self.flags().astype(float).reshape(1, -1)); beq = [self.kmass]
        res = linprog(c=-np.ones(n), A_ub=A, b_ub=b, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs',
                      options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
        if not res.success: return None
        mu = np.maximum(res.x, 0.0); yz = np.maximum(-res.ineqlin.marginals, 0.0)
        y = yz[:m]; self.z = yz[m:]
        self.lam = float(-res.eqlin.marginals[0]) if Aeq is not None else 0.0
        return mu, y, -res.fun

    # ---- clique separation on the current measure -> box cuts
    def separate(self, mu, time_limit, per_round, per_time, log, dmin=2e-4):
        sup = np.nonzero(mu > 1e-9)[0]
        poses = [self.poses[i] for i in sup]
        X, par, gs = expand_images(poses, None, self.t)
        w = mu[sup][par] / 8.0
        adj = sat_closed(X, X); np.fill_diagonal(adj, False)
        bw, bsz, dt, found = find_cliques(adj, w, time_limit, per_round, per_time, self.greedy_top)
        cuts = []; dropped = 0
        for K in found:
            K = np.array(K); r = boxes_from_clique(X[K], w[K], dmin=dmin)
            if r is None: dropped += 1; continue
            B, idx = r; cuts.append(B)
        new = self.add_cuts(cuts)
        return bw, bsz, dt, new

    # ---- clique part of the reduced cost of candidate poses (binding cuts only)
    def clique_cost(self, cand, zmin=1e-9):
        if not self.cuts or not len(self.z): return np.zeros(len(cand))
        X, par, gs = expand_images(cand, None, self.t); out = np.zeros(len(cand))
        for ci in np.nonzero(self.z > zmin)[0]:
            hit = in_boxes(X, self.cuts[ci])
            if hit.any(): out += self.z[ci] * np.bincount(par[hit], minlength=len(cand)) / 8.0
        return out


# ============================================================================== pricing helpers
def read_warm(path, t):
    lam = None; poses = []
    for line in open(path):
        if line.startswith('# t='): lam = t / float(line.split('=')[1].split()[0])
        q = line.split()
        if q and q[0] == 'pose': poses.append((float(q[1]), float(q[2]), math.radians(float(q[3]))))
    lam = lam or 1.0
    return [(cx * lam, cy * lam, th) for cx, cy, th in poses]


def read_branch_dual(path, t):
    poses = []; s = None
    for line in open(path):
        if line.startswith('#'):
            for tok in line.split():
                if tok.startswith('s='):
                    a, b = tok[2:].split('/'); s = int(a) / int(b)
            continue
        q = line.split(); poses.append((float(q[0]), float(q[1]), float(q[2])))
    lam = t / s if s else 1.0
    return [(cx * lam, cy * lam, th) for cx, cy, th in poses]


def corner_grid(t, r, pitch, dth_deg):
    out = []
    for th in np.arange(0.0, math.pi / 4 + 1e-9, math.radians(dth_deg)):
        th = PD.snap_angle(min(th, math.pi / 4)); w2 = PD.wid_of(th) / 2
        g = np.arange(w2 + PD.ADM, r + 1e-12, pitch)
        for cx in g:
            for cy in g: out.append((float(cx), float(cy), th))
    return out


# ============================================================================== main loop (refinement ladder)
def inner_converge(m, a, log, T0, stage):
    """fixed pose set: alternate LP / violated vertices / violated cliques until both are (nearly) exhausted.
    returns (mu, y, mass, M, kmax, feasible-scaled value, iterations)"""
    last = None
    for it in range(a.inner):
        t0 = time.time(); out = m.solve()
        if out is None: log("LP failed"); return last
        mu, y, mass = out; tlp = time.time() - t0
        t0 = time.time(); M, bad, nv, ovf = m.certify(mu, cap=2_000_000); tc = time.time() - t0
        t0 = time.time(); kmax, ksize, kdt, ncut = m.separate(mu, a.clique_time, a.per_round, a.per_time, log, dmin=a.dmin); ts = time.time() - t0
        Lf = mass / max(M, kmax, 1.0)
        nsup = int((mu > 1e-12).sum()); zpos = int((m.z > 1e-9).sum()); zsum = float(m.z.sum())
        cm = float(mu[m.flags()].sum()) if a.kmass is not None else None
        log(f"  stage {stage} it {it}: LP mass={mass:.6f} M={M:.6f} maxclique={kmax:.6f}(size {ksize}) -> feasible-scaled {Lf:.6f} | support={nsup} "
            f"cols={len(m.poses)} rows={len(m.pts)} cuts={len(m.cuts)} (binding {zpos}, dual mass {zsum:.4f}) lam={m.lam:.4f}{'' if cm is None else f' corner={cm:.6f}'} "
            f"| lp {tlp:.0f}s cert {tc:.0f}s sep {ts:.0f}s (+{ncut}) total {time.time()-T0:.0f}s")
        last = (mu, y, mass, M, kmax, Lf, it + 1)
        nrow = 0
        if len(bad):
            o = np.argsort(bad[:, 2])[::-1][:a.row_cap]; nrow = m.add_points(bad[o, :2])
        # age rows
        n_old_rows = len(y); slack = 1.0 - np.asarray(m.A[:n_old_rows] @ mu).ravel()
        stale = (y <= 1e-12) & (slack > 0.02)
        m.row_zero[:n_old_rows] = np.where(stale, m.row_zero[:n_old_rows] + 1, 0)
        keepr = m.row_zero < a.row_age
        if (~keepr).sum() > 0: m.drop_rows(keepr)
        if nrow == 0 and kmax <= 1 + a.ktol: log(f"  stage {stage}: inner converged (no violated vertex, max clique {kmax:.7f})"); break
        if time.time() - T0 > a.time: log("  time limit (inner)"); break
    return last


def price_columns(m, a, lib, y, mu, angs, t, log):
    ax, ay, aw = PD.atoms_from_dual(m, y, ymin=1e-12)
    if not len(ax): return 0, None, None
    cand = []; vmin = None
    if a.N > 0:
        try:
            pruned = (y > 1e-12).sum() > a.maxatoms // 8
            px, py, pw = PD.atoms_from_dual(m, y, ymin=1e-12, maxrows=a.maxatoms // 8) if pruned else (ax, ay, aw)
            vmin, tot, wit = PD.price_verifier(m, px, py, pw, a.N, a.topk, 'cc' + a.TAG, a.maxatoms, a.threads, log)
            cand += [(cx, cy, th) for cx, cy, th, v in wit]
        except Exception as e:
            log(f"   [verifier failed: {e}]")
    sx, sy, sww = PD.atoms_from_dual(m, y, ymin=1e-12, maxrows=a.sweep_rows)
    sw = PD.price_sweep(m, lib, sx, sy, sww, angs, a.sweep_k, 0.05, a.threads)
    cand += [(cx, cy, th) for cx, cy, th, v in sw]
    sup_poses = [m.poses[i] for i in np.nonzero(mu > 1e-9)[0]]
    nb = PD.neighbours(sup_poses, t)
    if a.kmass is not None: cand += corner_grid(t, a.r, a.corner_pitch, a.corner_dth)
    raw = [(cx, cy, th, v) for (cx, cy, th), v in zip(cand, PD.capture(lib, ax, ay, aw, cand, a.threads))] if cand else []
    ref = PD.refine(lib, ax, ay, aw, cand[:a.refine_cap], t, a.threads) if cand else []
    nbv = [(cx, cy, th, v) for (cx, cy, th), v in zip(nb, PD.capture(lib, ax, ay, aw, nb, a.threads))] if nb else []
    allc = raw + ref + nbv
    if not allc: return 0, None, vmin
    cc = m.clique_cost([(cx, cy, th) for cx, cy, th, v in allc])
    allc = [(cx, cy, th, v + c + ((m.lam if in_box(cx, cy, th, t, a.r) else 0.0) if a.kmass is not None else 0.0))
            for (cx, cy, th, v), c in zip(allc, cc)]
    allc.sort(key=lambda r: r[3]); rmin = allc[0][3]
    newc = [(cx, cy, th) for cx, cy, th, v in allc if v < 1.0 - a.price_tol][:a.cg_want]
    ncol = m.add_poses(newc)
    log(f"   pricing: candidates {len(allc)} improving {sum(1 for r in allc if r[3] < 1 - a.price_tol)} min reduced cost {rmin:.6f} verifier min {vmin} | +cols {ncol}")
    return ncol, rmin, vmin


def main_loop(a):
    t = float(Fr(a.T)); tag = a.TAG
    os.makedirs(RUNS, exist_ok=True)
    lf = open(os.path.join(RUNS, f'cc_{tag}.log'), 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"clique_ceiling t={t} tag={tag} args={vars(a)}")
    lib = PD.build_lib(); T0 = time.time()
    m = CModel(t, lib, a.threads, log, kmass=a.kmass, r=a.r, margin=a.margin); m.greedy_top = a.greedy_top
    seeds = PD.seed_poses(t, a.seed_pitch, a.seed_dth) if a.seed_pitch > 0 else []
    for w in a.warm: seeds += read_warm(w, t); log(f"warm: {w}")
    if a.branch_dual: seeds += read_branch_dual(a.branch_dual, t); log(f"branch dual: {a.branch_dual}")
    if a.kmass is not None: seeds += corner_grid(t, a.r, 0.05, 5.0)
    if a.seed_neighbours: seeds += PD.neighbours(seeds[:3000], t)
    m.add_poses(seeds)
    m.add_points(PD.seed_points(t, a.row_pitch)); m.add_points(PD.corners(m.poses, t))
    log(f"seed: {len(m.poses)} poses, {len(m.pts)} rows, nnz={m.A.nnz} ({time.time()-T0:.0f}s)")
    if a.kmass is not None: log(f"leaf mode: corner box r={a.r}, corner mass = {a.kmass}, flagged poses {int(m.flags().sum())}")
    hist = []; best = dict(L=0.0, mass=0.0); angs = np.arange(0.0, math.pi / 4 + 1e-9, math.radians(a.raster_dth))
    for stage in range(a.stages):
        res = inner_converge(m, a, log, T0, stage)
        if res is None: break
        mu, y, mass, M, kmax, Lf, its = res
        out = m.solve()                                   # duals consistent with the current rows (ageing may have dropped some)
        if out is not None: mu, y, _ = out
        nsup = int((mu > 1e-12).sum())
        write_support(m, mu, y, t, tag + '_last', mass, M, kmax, stage)
        if Lf > best['L'] and kmax <= 1 + a.ktol and M <= 1 + 1e-6:
            best = dict(L=Lf, mass=mass, M=M, kmax=kmax, stage=stage, support=nsup); write_support(m, mu, y, t, tag, mass, M, kmax, stage)
        ncol, rmin, vmin = price_columns(m, a, lib, y, mu, angs, t, log)
        rec = dict(stage=stage, t=time.time() - T0, mass=mass, M=M, kmax=kmax, Lf=Lf, inner_its=its, support=nsup, cols=len(m.poses), rows=len(m.pts),
                   cuts=len(m.cuts), new_cols=ncol, rmin=rmin, vmin=vmin, lam=m.lam, corner=(float(mu[m.flags()].sum()) if a.kmass is not None else None))
        hist.append(rec)
        log(f"STAGE {stage}: value {Lf:.6f} (LP {mass:.6f}, M {M:.6f}, maxclique {kmax:.6f}) cols {len(m.poses)} cuts {len(m.cuts)} min reduced cost {rmin} +cols {ncol}")
        json.dump(dict(t=t, tag=tag, args=vars(a), best=best, hist=hist), open(os.path.join(RUNS, f'cc_{tag}.json'), 'w'), indent=1)
        if ncol == 0: log("converged: no improving column"); break
        if time.time() - T0 > a.time: log("time limit"); break
    log(f"BEST t={t} feasible-scaled mass={best['L']:.6f} (stage {best.get('stage')}, LP mass {best.get('mass')}, M {best.get('M')}, maxclique {best.get('kmax')}) wall {time.time()-T0:.0f}s")
    return m, best, log


def write_support(m, mu, y, t, tag, mass, M, kmax, rd):
    with open(os.path.join(RUNS, f'cc_{tag}_support.txt'), 'w') as f:
        f.write(f"# t={t} round={rd} mass={mass:.9f} M={M:.9f} maxclique={kmax:.9f} L={mass/max(M,kmax,1.0):.9f}\n")
        f.write("# D4-symmetrised measure: pose cx cy theta_deg mu  (mass mu/8 on each of the 8 dihedral images); clique cuts below\n")
        for i in np.nonzero(mu > 1e-9)[0]:
            cx, cy, th = m.poses[i]; f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} {mu[i]:.12f}\n")
        for ci, B in enumerate(m.cuts):
            z = m.z[ci] if ci < len(m.z) else 0.0
            f.write(f"cut {z:.9f} " + " ".join(f"{b[0]:.6f},{b[1]:.6f},{b[2]:.6f},{b[3]:.2e},{b[4]:.2e}" for b in B) + "\n")
        for i in range(len(m.pts)):
            f.write(f"point {m.pts[i,0]:.10f} {m.pts[i,1]:.10f} {y[i]:.10f}\n")


# ============================================================================== exact certification
def sq_proj(sq, ux, uy):
    """min and max of the integer dot products of the square's corners with (ux, uy), over Dq"""
    d = [ux * X + uy * Y for X, Y in sq[6]]
    return min(d), max(d), sq[7]


def closed_intersect_exact(si, sj):
    """closed squares meet <=> not separated on any of the 4 edge normals (integers only)"""
    for s in (si, sj):
        a, b = s[3], s[4]
        for ux, uy in ((a, b), (-b, a)):
            lo1, hi1, D1 = sq_proj(si, ux, uy); lo2, hi2, D2 = sq_proj(sj, ux, uy)
            if hi1 * D2 < lo2 * D1 or hi2 * D1 < lo1 * D2: return False
    return True


def max_clique_int(nbrset, w, time_limit):
    """maximum-weight clique with INTEGER weights (exact bound arithmetic); returns (weight, members, complete)"""
    n = len(w); best = [0, []]; t0 = time.time(); timed_out = [False]
    def colour_bound(cands):
        classes = []; cw = 0
        for v in sorted(cands, key=lambda v: -w[v]):
            for cl in classes:
                if not (nbrset[v] & cl): cl.add(v); break
            else: classes.append({v}); cw += w[v]
        return cw
    def expand(cands, cur, curw):
        if time.time() - t0 > time_limit: timed_out[0] = True; return
        if not cands:
            if curw > best[0]: best[0] = curw; best[1] = list(cur)
            return
        if curw + colour_bound(cands) <= best[0]: return
        cl = sorted(cands, key=lambda v: -w[v])
        while cl:
            v = cl.pop(0)
            expand([u for u in cl if u in nbrset[v]], cur + [v], curw + w[v])
            if timed_out[0]: return
            if curw + colour_bound(cl) <= best[0]: return
    expand([v for v in range(n) if w[v] > 0], [], 0)
    return best[0], best[1], not timed_out[0]


def exact_certify(a, m, log):
    """snap the best measure, then exact rows + exact clique graph, LP polish with clique cuts, integerize, check"""
    t = Fr(a.T); DE.TN, DE.TD = t.numerator, t.denominator; DE.T = t
    DE.in_F = lambda v: v[0] >= 0 and v[0] <= v[1] and 2 * DE.TD * v[1] <= DE.TN * v[2]
    src = os.path.join(RUNS, f'cc_{a.TAG}_support.txt')
    poses_f = []; mus_f = []
    for line in open(src):
        q = line.split()
        if q and q[0] == 'pose': poses_f.append((float(q[1]), float(q[2]), math.radians(float(q[3])))); mus_f.append(float(q[4]))
    log(f"exact: {len(poses_f)} float poses from {src}; t = {DE.TN}/{DE.TD}; snap Q={a.Q} Dc={a.Dc}")
    poses = [DE.snap_pose(cx, cy, th, a.Q, a.Dc) for cx, cy, th in poses_f]
    squares, pose_of = DE.build_squares(poses)
    log(f"  {len(squares)} images, all exactly admissible")
    verts = DE.enumerate_vertices(squares, full=False, procs=a.procs)
    inc = DE.incidences(a.procs); pat = DE.patterns(inc, pose_of); keys = list(pat.keys())
    log(f"  {len(keys)} distinct incidence patterns (point rows)")
    # exact closed-intersection graph on the images (prefilter by float centre distance)
    t0 = time.time(); n = len(squares)
    C = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2]))] for s in squares])
    nbr = [set() for _ in range(n)]; npairs = 0
    for i in range(n):
        d2 = (C[:, 0] - C[i, 0]) ** 2 + (C[:, 1] - C[i, 1]) ** 2
        for j in np.nonzero((d2 <= 2.0 + 1e-6) & (np.arange(n) > i))[0]:
            npairs += 1
            if closed_intersect_exact(squares[i], squares[int(j)]): nbr[i].add(int(j)); nbr[int(j)].add(i)
    log(f"  exact closed-intersection graph: {n} images, {sum(map(len, nbr))//2} edges ({npairs} candidate pairs, {time.time()-t0:.0f}s)")
    adj = np.zeros((n, n), bool)
    for i in range(n):
        for j in nbr[i]: adj[i, j] = True
    nP = len(poses); DM = a.DM
    flags = None
    if a.kmass is not None:
        flags = np.array([in_box(float(cx), float(cy), 0.0, float(t), a.r) for (p, q, cx, cy) in poses], bool)
    # LP polish with exact rows and clique cuts found on the exact graph
    rows = []; cols = []; vals = []
    for i, key in enumerate(keys):
        for k, c in key: rows.append(i); cols.append(k); vals.append(c / 8.0)
    Apt = sp.csr_matrix((vals, (rows, cols)), shape=(len(keys), nP))
    cuts = []                  # lists of image indices
    mu = None
    for rd in range(a.exact_rounds):
        A = Apt if not cuts else sp.vstack([Apt, sp.csr_matrix((np.concatenate([np.full(len(c), 1 / 8) for c in cuts]),
                                            (np.concatenate([np.full(len(c), k) for k, c in enumerate(cuts)]), np.concatenate([[pose_of[i] for i in c] for c in cuts]))), shape=(len(cuts), nP))], format='csr')
        A.sum_duplicates()
        Aeq = beq = None
        if flags is not None: Aeq = sp.csr_matrix(flags.astype(float).reshape(1, -1)); beq = [a.kmass]
        res = linprog(-np.ones(nP), A_ub=A, b_ub=np.full(A.shape[0], 1.0 - a.margin), A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs',
                      options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
        assert res.status == 0, res.message
        mu = np.maximum(res.x, 0.0)
        w = mu[pose_of] / 8.0
        keep = np.nonzero(w > 1e-12)[0]
        bw, bsz, dt, found = find_cliques(adj[np.ix_(keep, keep)], w[keep], a.clique_time, a.per_round, a.per_time, a.greedy_top, thresh=1 + 1e-9)
        log(f"  polish {rd}: LP mass {-res.fun:.9f}, cuts {len(cuts)}, max clique (float) {bw:.9f} size {bsz} ({dt:.0f}s), +{len(found)} cuts")
        if bw <= 1 + 1e-9: break
        for K in found: cuts.append([int(keep[i]) for i in K])
    # integerize
    MU = [int(math.floor(float(x) * DM)) for x in mu]
    if flags is not None:
        S = sum(MU[k] for k in range(nP) if flags[k]); d = a.kmass * DM - S
        h = max((k for k in range(nP) if flags[k]), key=lambda k: MU[k]); MU[h] += d
        log(f"  leaf: corner mass adjusted by {d}/{DM} on pose {h} to make it exactly {a.kmass}")
    def exact_maxima(MU):
        M, key = DE.exact_max(pat, MU, DM)
        wi = [MU[pose_of[i]] for i in range(n)]
        Kw, Kc, complete = max_clique_int(nbr, wi, a.clique_time * 4)
        return M, Fr(Kw, 8 * DM), Kc, complete
    M, K, Kc, complete = exact_maxima(MU)
    log(f"  exact: mass {sum(MU)/DM:.9f}, M = {float(M):.12f}, max clique = {float(K):.12f} (size {len(Kc)}, B&B complete: {complete})")
    if not complete: log("  WARNING: clique branch-and-bound timed out; NOT certified"); return None
    s = max(M, K)
    if s > 1:
        if flags is not None: log("  leaf: maxima exceed 1 after rounding; scaling breaks the corner equality -- reported as scaled")
        MU = [(x * s.denominator) // s.numerator for x in MU]
        M, K, Kc, complete = exact_maxima(MU)
        log(f"  scaled by 1/{float(s):.12f}: mass {sum(MU)/DM:.9f}, M = {float(M):.12f}, max clique = {float(K):.12f} (complete: {complete})")
        if not complete: return None
    assert M <= 1 and K <= 1
    mass = Fr(sum(MU), DM); cm = Fr(sum(MU[k] for k in range(nP) if flags[k]), DM) if flags is not None else None
    log(f"  CERTIFIED (exact): t = {DE.TN}/{DE.TD}, mass = {mass} = {float(mass):.9f}, coverage <= {M} <= 1, closed-clique mass <= {K} <= 1"
        + (f", corner mass = {cm} = {float(cm):.9f}" if cm is not None else "") + f"  ({'>=' if mass >= 12 else '<'} 12)")
    out = os.path.join(RUNS, f'cc_{a.TAG}_exact.txt')
    with open(out, 'w') as f:
        f.write(f"# t = {DE.TN}/{DE.TD} exact; D4-symmetrised measure; closed-intersection clique-feasible\n")
        f.write(f"# mass = {mass} = {float(mass):.12f}; M = {M}; maxclique = {K}" + (f"; corner mass = {cm}" if cm is not None else "") + "\n")
        f.write("# pose p q cx cy mass : theta = 2 arctan(p/q); mass/8 on each of the 8 dihedral images (dual_exact.py convention)\n")
        for (p, q, cx, cy), x in zip(poses, MU):
            if x > 0: f.write(f"pose {p} {q} {cx} {cy} {Fr(x, DM)}\n")
    log(f"  wrote {out}")
    return dict(mass=str(mass), mass_float=float(mass), M=str(M), K=str(K), K_size=len(Kc), corner=str(cm) if cm is not None else None,
                images=n, vertices_F=len(verts), patterns=len(keys), cuts=len(cuts))


def cmd_check(a):
    """independent exact re-check of an exact support file: coverage and closed-clique mass"""
    t = Fr(a.T); DE.TN, DE.TD = t.numerator, t.denominator; DE.T = t
    DE.in_F = lambda v: v[0] >= 0 and v[0] <= v[1] and 2 * DE.TD * v[1] <= DE.TN * v[2]
    sup = DE.read_exact_support(a.check)
    poses = [(p, q, cx, cy) for (p, q, cx, cy, m_) in sup]; masses = [m_ for (_, _, _, _, m_) in sup]
    DM = 1
    for x in masses: DM = lcm(DM, x.denominator)
    MU = [x.numerator * (DM // x.denominator) for x in masses]
    squares, pose_of = DE.build_squares(poses)
    DE.enumerate_vertices(squares, full=a.full, procs=a.procs); inc = DE.incidences(a.procs)
    best = max(sum(MU[pose_of[s]] for s in lst) for lst in inc); M = Fr(best, 8 * DM)
    n = len(squares); C = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2]))] for s in squares]); nbr = [set() for _ in range(n)]
    for i in range(n):
        d2 = (C[:, 0] - C[i, 0]) ** 2 + (C[:, 1] - C[i, 1]) ** 2
        for j in np.nonzero((d2 <= 2.0 + 1e-6) & (np.arange(n) > i))[0]:
            if closed_intersect_exact(squares[i], squares[int(j)]): nbr[i].add(int(j)); nbr[int(j)].add(i)
    Kw, Kc, complete = max_clique_int(nbr, [MU[pose_of[i]] for i in range(n)], 3600)
    K = Fr(Kw, 8 * DM); mass = Fr(sum(MU), DM)
    print(f"CHECK t={t}: {len(poses)} poses, {n} images; mass = {mass} = {float(mass):.9f}; M = {M} = {float(M):.12f}; "
          f"closed max clique = {K} = {float(K):.12f} (size {len(Kc)}, complete {complete}); "
          f"{'FEASIBLE' if M <= 1 and K <= 1 and complete else 'NOT FEASIBLE / INCOMPLETE'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('T'); ap.add_argument('TAG')
    ap.add_argument('--warm', action='append', default=[], help='support file(s) to seed from (centres rescaled)')
    ap.add_argument('--branch-dual', default=None); ap.add_argument('--kmass', type=float, default=None); ap.add_argument('--r', type=float, default=1.0)
    ap.add_argument('--margin', type=float, default=0.0, help='rows use 1 - margin (leaf mode: room for the exact corner adjustment)')
    ap.add_argument('--stages', type=int, default=12); ap.add_argument('--inner', type=int, default=30); ap.add_argument('--ktol', type=float, default=1e-5)
    ap.add_argument('--price-tol', type=float, default=1e-4)
    ap.add_argument('--time', type=float, default=3600)
    ap.add_argument('--seed-pitch', type=float, default=0.12); ap.add_argument('--seed-dth', type=float, default=15.0); ap.add_argument('--row-pitch', type=float, default=0.04)
    ap.add_argument('--seed-neighbours', action='store_true')
    ap.add_argument('--threads', type=int, default=4); ap.add_argument('--N', type=int, default=1000); ap.add_argument('--topk', type=int, default=8)
    ap.add_argument('--maxatoms', type=int, default=6000); ap.add_argument('--cg-want', type=int, default=600); ap.add_argument('--refine-cap', type=int, default=3000)
    ap.add_argument('--row-cap', type=int, default=20000); ap.add_argument('--row-age', type=int, default=8)
    ap.add_argument('--raster-dth', type=float, default=0.25); ap.add_argument('--sweep-rows', type=int, default=350); ap.add_argument('--sweep-k', type=int, default=6)
    ap.add_argument('--corner-pitch', type=float, default=0.02); ap.add_argument('--corner-dth', type=float, default=2.0)
    ap.add_argument('--clique-time', type=float, default=120); ap.add_argument('--per-round', type=int, default=40); ap.add_argument('--per-time', type=float, default=5)
    ap.add_argument('--greedy-top', type=int, default=300); ap.add_argument('--dmin', type=float, default=2e-4)
    ap.add_argument('--exact', action='store_true'); ap.add_argument('--exact-rounds', type=int, default=40)
    ap.add_argument('--Q', type=int, default=100000); ap.add_argument('--Dc', type=int, default=1000000); ap.add_argument('--DM', type=int, default=10 ** 9); ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--check', default=None, help='exact support file to re-check (no search)'); ap.add_argument('--full', action='store_true')
    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    if a.check: return cmd_check(a)
    m, best, log = main_loop(a)
    if a.exact and best['L'] > 0:
        res = exact_certify(a, m, log)
        js = json.load(open(os.path.join(RUNS, f'cc_{a.TAG}.json'))); js['exact'] = res
        json.dump(js, open(os.path.join(RUNS, f'cc_{a.TAG}.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
