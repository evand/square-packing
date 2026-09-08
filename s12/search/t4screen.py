#!/usr/bin/env python3
"""t4screen.py -- packing-side value of one level-2 leaf, strengthened by anchor cliques.

    max  sum_S mu_S
    s.t. coverage(p) = sum_S mu_S [p in S]  <= 1        for every point p of [0,t]^2
         mu({S : centre(S) in C_i}) = kc_i              i = 0..3   (corner boxes [0,r]^2 and images)
         mu({S : centre(S) in W_j}) = kw_j              j = 0..7   (wall slots)
         mu(K)            = sum_{S in K} mu_S   <= 1    for every anchor clique K = K(p, A)
         mu >= 0

over CLOSED unit squares S contained in the CLOSED container [0,t]^2 (`admissible`: the centre may
sit anywhere with the square still inside, so a square may touch a wall; a point on the boundary of
a square counts as being in it).  The measure is NOT D4-symmetrised: an asymmetric slot pattern
needs the whole container, and every clique image is a separate constraint (stronger than the
orbit-averaged clique column a symmetric model can carry).

This is `search/level2_lp.py` (regions, row generation, reduced-cost pricing) glued to the anchor
clique family of `search/anchorclique.py` / `search/anchorsep.py` (`kpa` = the wall-band
Lemma-2 family, `kseg` = the general interior segment anchor of `search/RECONCILE.md` 2).
`level2_lp.py` has no cliques; `clique_ceiling.py` has cliques but no region equalities and is
D4-symmetric, so it cannot express an asymmetric slot pattern.  Nothing in the certificate
pipeline is touched: this is a packing-side screening instrument only.

Direction of every error, which is the whole point:

  * POSES are a finite set (seed grid + warm starts + pricing).  Restricting the columns can only
    LOWER the value: the reported number is an optimistic (low) estimate of the leaf's true value.
  * ROWS (points) are a finite set, which INFLATES the value.  The loop adds every violated vertex
    of the arrangement of the support squares (`packing_dual.py`'s exact vertex kernel) until the
    certified maximum coverage `M` is `<= 1 + 1e-9`.  Read `LP` together with `M`: only when
    `M <= 1` is `LP` a value the current poses really attain.
  * CLIQUES are separated heuristically (a grid over anchor point, direction, offset and
    half-length), which also INFLATES the value.  `kmax` is the largest `mu(K)` the separator
    could find at the end of the inner loop; `kmax <= 1` means the measure is feasible for
    everything the separator can see, not for the whole family.
  * Anchor membership (`anchorclique.member`) drops poses within 1e-9 of a piece boundary.  On the
    packing side that under-counts members, i.e. imposes `mu(K') <= 1` for a slightly smaller
    `K' subset K` -- a valid but marginally weaker constraint.

So: `L = LP / max(M, kmax, 1)` is the value of a measure feasible for everything checked, and a
leaf value close to or above 12 after aggressive refinement is evidence (not proof) that the leaf
does not close; a value well below 12 is the weaker, optimistic direction.

  python3 search/t4screen.py 4.0 TAG --corners 1111 --patterns 01010101 --cliques --cq-interior
  python3 search/t4screen.py 3.98 CAL --corners 1111 --cliques --cq-interior \\
      --warm /path/runs/branch_J16i_dual.txt --cq-load /path/runs/branch_J16i_cliques.txt
"""
import argparse, json, math, os, sys, time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import packing_dual as pd                                            # noqa: E402
import level2_lp as L2                                               # noqa: E402
import anchorclique as AC                                            # noqa: E402
import anchorsep as ASEP                                             # noqa: E402

RUNS = os.path.join(REPO, 'runs')
pd.RUNS = RUNS
L2.RUNS = RUNS
H_UNIT = 0.5            # half-side of a CLOSED UNIT square: the packing-side credit rule


# ====================================================================== the model
class CLeaf(L2.Leaf):
    """`level2_lp.Leaf` (full-container poses x points LP with region equalities) plus anchor
    cliques as extra `<= 1` rows.  Every clique is a single clique of the container, not an
    orbit."""

    def __init__(self, t, r, lib, threads, log, D=1000000):
        super().__init__(t, r, lib, threads, log)
        self.D = D
        self.sfr = Fr(int(round(t * D)), D)
        self.cliques = []                # list of clique objects (anchors, pieces)
        self.cparams = []                # generating integer parameters, for the checkpoint
        self.ckeys = set()
        self.K = sp.csr_matrix((0, 0))   # cliques x poses, entry 1 if the pose is in the clique
        self.z = np.zeros(0)             # clique duals of the last solve
        self.cage = np.zeros(0, dtype=int)   # consecutive solves with zero dual, for sifting
        self.pool = []                   # loaded candidate cliques, added only when violated
        self.rage = np.zeros(0, dtype=int)   # consecutive solves a point row has been slack and dual-free
        self.rfix = 0                        # the first `rfix` rows (the seed grid) are never dropped
        self.lp_secs = 0.0

    # ---------------------------------------------------------------- row ageing
    def add_points(self, pts):
        n = super().add_points(pts)
        if n:
            self.rage = np.concatenate([self.rage, np.zeros(n, dtype=int)])
        return n

    def drop_rows(self, keep):
        idx = np.nonzero(keep)[0]
        self.rage = self.rage[idx]
        super().drop_rows(keep)

    def sift_rows(self, mu, y, age, floor=4000):
        """drop point rows that have been slack and dual-free for `age` consecutive solves
        (`clique_ceiling.py`'s rule).  Coverage is still certified over the whole continuum by
        `certify`, so a wrongly dropped row simply comes back as a violated arrangement vertex."""
        if age <= 0 or len(self.pts) <= floor:
            return 0
        slack = 1.0 - np.asarray(self.A @ mu).ravel()
        stale = (y <= 1e-12) & (slack > 0.02)
        self.rage = np.where(stale, self.rage + 1, 0)
        keep = (self.rage < age) | (np.arange(len(self.rage)) < self.rfix)
        if keep.all() or keep.sum() < floor:
            return 0
        n = int((~keep).sum())
        self.drop_rows(keep)
        return n

    # ---------------------------------------------------------------- geometry helpers
    def rowsarr(self, poses=None):
        """the (cx, cy, theta, h) rows `anchorclique.member` wants, with h = 1/2 (unit squares)"""
        P = np.asarray(poses if poses is not None else self.poses, dtype=float).reshape(-1, 3)
        return np.column_stack([P[:, 0], P[:, 1], P[:, 2], np.full(len(P), H_UNIT)])

    def clique_matrix(self, cliques, poses):
        if not len(cliques) or not len(poses):
            return sp.csr_matrix((len(cliques), len(poses)))
        Q = self.rowsarr(poses)
        rot = AC._rot(Q)
        rows, cols = [], []
        for i, cl in enumerate(cliques):
            nz = np.nonzero(AC.member(cl, Q, rot))[0]
            rows.append(np.full(len(nz), i)); cols.append(nz)
        if not rows:
            return sp.csr_matrix((len(cliques), len(poses)))
        R = np.concatenate(rows); C = np.concatenate(cols)
        return sp.csr_matrix((np.ones(len(R)), (R, C)), shape=(len(cliques), len(poses)))

    # ---------------------------------------------------------------- columns
    def add_poses(self, cand):
        n0 = len(self.poses)
        k = super().add_poses(cand)
        if self.cliques:
            if k:
                Kn = self.clique_matrix(self.cliques, self.poses[n0:])
                self.K = sp.hstack([self.K, Kn], format='csr') if self.K.shape[1] else Kn
            elif self.K.shape[1] != len(self.poses):
                self.K = self.clique_matrix(self.cliques, self.poses)
        return k

    def add_cliques(self, cls, params=None):
        new, npar = [], []
        for i, cl in enumerate(cls):
            k = AC.key(cl)
            if k in self.ckeys:
                continue
            self.ckeys.add(k); new.append(cl); npar.append(None if params is None else params[i])
        if not new:
            return 0
        Kn = self.clique_matrix(new, self.poses)
        self.K = sp.vstack([self.K, Kn], format='csr') if self.K.shape[0] else Kn
        self.cliques += new; self.cparams += npar
        self.cage = np.concatenate([self.cage, np.zeros(len(new), dtype=int)])
        return len(new)

    def sift_poses(self, keep):
        """drop pose columns (a restriction: it can only LOWER the value).  Used to keep the LP
        solvable across restarts -- the pose set otherwise grows by ~3000 columns per pricing
        stage.  Poses carrying mass are always kept."""
        idx = np.nonzero(keep)[0]
        if len(idx) == len(self.poses):
            return 0
        n = len(self.poses) - len(idx)
        self.poses = [self.poses[i] for i in idx]
        self.reg = [self.reg[i] for i in idx]
        self.key = {}
        for j, q in enumerate(self.poses):
            self.key[(round(q[0] * 1e7), round(q[1] * 1e7), round(q[2] / pd.ANG_UNIT))] = j
        self.A = self.A[:, idx]
        if self.cliques:
            self.K = self.K[:, idx]
        return n

    def sift_cliques(self, age):
        """drop clique rows that have carried no dual for `age` consecutive solves.  A dropped
        clique that is still violated is separated again on the next call, so this only trims
        rows the LP is not using -- clique rows are dense and a few thousand of them make the
        LP unusable (search/RECONCILE.md, branch.py's restricted master)."""
        if not self.cliques:
            return 0
        keep = self.cage < age
        if keep.all():
            return 0
        idx = np.nonzero(keep)[0]
        self.cliques = [self.cliques[i] for i in idx]
        self.cparams = [self.cparams[i] for i in idx]
        self.ckeys = set(AC.key(cl) for cl in self.cliques)
        self.K = self.K[idx]
        self.cage = self.cage[idx]
        self.z = self.z[idx] if len(self.z) == len(keep) else np.zeros(len(idx))
        return int((~keep).sum())

    def clique_cost(self, cand, zmin=1e-9, topz=500):
        """the clique part of the reduced cost of candidate poses: sum_K z_K [pose in K]"""
        out = np.zeros(len(cand))
        if not self.cliques or not len(self.z):
            return out
        idx = np.nonzero(self.z > zmin)[0]
        if len(idx) > topz:
            idx = idx[np.argsort(-self.z[idx])[:topz]]
        if not len(idx):
            return out
        Q = self.rowsarr(cand); rot = AC._rot(Q)
        for ci in idx:
            m = AC.member(self.cliques[ci], Q, rot)
            if m.any():
                out[m] += self.z[ci]
        return out

    # ---------------------------------------------------------------- LP
    def solve(self, kc, kw, method='highs', no_cliques=False):
        n = len(self.poses)
        reg = np.array(self.reg)
        erows, evals = [], []
        for i in range(4):
            if kc[i] is not None:
                erows.append((reg == i).astype(float)); evals.append(kc[i])
        for j in range(8):
            if kw[j] is not None:
                erows.append((reg == 4 + j).astype(float)); evals.append(kw[j])
        Aeq = np.array(erows) if erows else None
        beq = np.array(evals) if erows else None
        m = len(self.pts)
        use_cq = bool(self.cliques) and not no_cliques
        A = sp.vstack([self.A, self.K], format='csr') if use_cq else self.A
        t0 = time.time()
        res = None
        for meth in (method, 'highs-ds', 'highs-ipm'):
            res = linprog(c=-np.ones(n), A_ub=A, b_ub=np.ones(A.shape[0]),
                          A_eq=Aeq, b_eq=beq, bounds=(0, None), method=meth)
            if res.success:
                break
            self.log(f'   [LP {meth} failed: status {res.status} {res.message}]')
        self.lp_secs = time.time() - t0
        if not res.success:
            return None
        mu = np.maximum(res.x, 0.0)
        yz = np.maximum(-res.ineqlin.marginals, 0.0)
        y = yz[:m]
        if not no_cliques:
            self.z = yz[m:] if use_cq else np.zeros(0)
            if use_cq and len(self.z) == len(self.cage):
                self.cage = np.where(self.z > 1e-9, 0, self.cage + 1)
        lam = -res.eqlin.marginals if erows else np.zeros(0)
        return mu, y, -res.fun, lam

    # ---------------------------------------------------------------- clique separation
    def mass_of(self, cl, Q, rot, w):
        return float(w[AC.member(cl, Q, rot)].sum())

    def separate(self, mu, a, log, want=None):
        """separate anchor cliques against the current measure; returns (kmax, n_added, best)"""
        sup = np.nonzero(mu > 1e-9)[0]
        if not len(sup):
            return 0.0, 0, None
        Q = self.rowsarr([self.poses[i] for i in sup])
        w = mu[sup]
        rot = AC._rot(Q)
        want = a.cq_want if want is None else want
        found = []                       # (mass, mass of the point clique, params, clique)
        # (a) the loaded pool: cheap, and it is the J16i warm start's own cliques
        for cl, par in self.pool:
            v = self.mass_of(cl, Q, rot, w)
            if v > 1.0 + a.ktol:
                found.append((v, 0.0, par, cl))
        # (b) the wall band family (anchorclique.kpa, Lemma 2)
        if a.cq_wall:
            found += ASEP.separate(Q, w, self.t, self.D, pitch=a.cq_pitch, top=a.cq_top,
                                   band=1.0, sym=False)
        # (c) the general interior segment family (anchorclique.kseg)
        if a.cq_interior:
            found += ASEP.separate_interior(Q, w, self.t, self.D, pitch=a.cq_int_pitch,
                                            top=a.cq_int_top, ndir=a.cq_ndir, neps=a.cq_neps,
                                            nrho=a.cq_nrho, sym=False)
        found.sort(key=lambda r: -r[0])
        kmax = found[0][0] if found else 0.0
        add, addp = [], []
        for (mk, mp, par, cl) in found:
            if mk <= 1.0 + a.ktol or len(add) >= want:
                break
            add.append(cl); addp.append(par)
        n = self.add_cliques(add, addp)
        return kmax, n, (found[0] if found else None)

    def load_pool(self, path, log):
        """read a branch.py anchor-clique checkpoint (`runs/branch_TAG_cliques.txt`) as CANDIDATES.
        The file's container may differ from ours: anchors are rescaled by t/s (Lemma 0 holds for
        every p and A, so a rescaled clique is still a clique -- only its usefulness changes)."""
        with open(path) as f:
            head = f.readline().split()
            Ks, Ds, sym = int(head[0]), int(head[1]), int(head[2])
            s0 = Ks / Ds
            lam = self.t / s0
            n = 0
            for line in f:
                q = [int(v) for v in line.split()]
                if len(q) not in (5, 6):
                    continue
                # rescale to our container and our denominator
                qq = [int(round(v / Ds * lam * self.D)) for v in q]
                if len(q) == 5:
                    qq[2] = q[2]                    # the wall index is not a coordinate
                    cl0 = AC.kpa(self.sfr, self.D, qq[0], qq[1], qq[2], qq[3], qq[4])
                else:
                    cl0 = AC.kseg(self.sfr, self.D, *qq)
                if cl0 is None:
                    continue
                for cl in (AC.images(self.sfr, cl0) if sym else [cl0]):
                    k = AC.key(cl)
                    if k in self.ckeys:
                        continue
                    self.ckeys.add(k)               # reserve the key: pool members are candidates
                    self.pool.append((cl, None)); n += 1
        # pool members are candidates, not constraints: release their keys so `add_cliques` works
        for cl, _ in self.pool:
            self.ckeys.discard(AC.key(cl))
        log(f'   clique pool from {os.path.basename(path)}: {n} candidate cliques '
            f'(container {s0:.6f} -> {self.t}, sym={sym})')
        return n


# ====================================================================== pricing
def price(leaf, y, lam_of_region, a, log, tag=''):
    """best reduced-cost poses on a grid: rc = 1 - capture(point duals) - lambda(region) - clique duals"""
    sel = np.nonzero(y > 1e-9)[0]
    if not len(sel):
        return [], None, None
    ax = np.ascontiguousarray(leaf.pts[sel, 0])
    ay = np.ascontiguousarray(leaf.pts[sel, 1])
    aw = np.ascontiguousarray(y[sel])
    o = np.argsort(ax, kind='stable')
    ax, ay, aw = np.ascontiguousarray(ax[o]), np.ascontiguousarray(ay[o]), np.ascontiguousarray(aw[o])

    def rc_of(poses):
        cap = pd.capture(leaf.lib, ax, ay, aw, poses, leaf.threads)
        reg = np.array([leaf.region_of(p) for p in poses])
        lam = np.array([lam_of_region.get(int(rr), 0.0) for rr in reg])
        return 1.0 - cap - lam - leaf.clique_cost(poses)

    cand = L2.seed_poses_full(leaf.t, a.price_pitch, a.price_dth)
    rc = rc_of(cand)
    # STRATIFY BY REGION.  With an equality `mu(R) = k` the region's multiplier `lam_R` is free and
    # can be strongly negative (the branch forces mass into a region the LP does not want it in),
    # so the reduced cost of a candidate centred there is `1 - capture + |lam_R|` and can be in the
    # hundreds -- while adding such a column only redistributes the region's fixed mass `k` and
    # moves the objective by nothing.  A plain top-`cg_want` by reduced cost is then entirely
    # frame poses, and the only region that can actually raise the objective in an `m = 4` leaf --
    # the interior, whose multiplier is 0 because it carries no equality -- never gets priced.
    creg = np.array([leaf.region_of(p) for p in cand])
    idx = []
    regs = sorted(set(int(v) for v in creg))
    quota = max(a.cg_want // max(len(regs), 1), 1)
    for rr in regs:
        w = np.nonzero(creg == rr)[0]
        if not len(w):
            continue
        o = w[np.argsort(-rc[w])[:quota]]
        idx += [int(i) for i in o if rc[i] > a.price_tol]
    extra = [int(i) for i in np.argsort(-rc) if rc[i] > a.price_tol and int(i) not in set(idx)]
    idx = (idx + extra)[:a.cg_want]
    best = [cand[i] for i in idx]
    # coordinate descent on the point-capture part only (the clique part is not differentiable
    # on a grid; the refined poses are re-priced honestly below)
    ref = pd.refine(leaf.lib, ax, ay, aw, best, leaf.t, leaf.threads) if best else []
    out = list(best) + list(ref)
    gap = float(rc.max())
    gi = np.nonzero(creg == 12)[0]
    gint = float(rc[gi].max()) if len(gi) else float('nan')
    log(f'   pricing{tag}: {len(cand)} candidates, best rc {gap:+.6f} (interior {gint:+.6f}), '
        f'kept {len(best)} over {len(regs)} regions (+{len(ref)} refined)')
    return out, gap, (ax, ay, aw)


def pose_rc(leaf, y, lam_of_region):
    """reduced cost of every column of the model at the current duals"""
    sel = np.nonzero(y > 1e-12)[0]
    if not len(sel):
        return np.zeros(len(leaf.poses))
    ax = np.ascontiguousarray(leaf.pts[sel, 0]); ay = np.ascontiguousarray(leaf.pts[sel, 1])
    aw = np.ascontiguousarray(y[sel])
    o = np.argsort(ax, kind='stable')
    ax, ay, aw = np.ascontiguousarray(ax[o]), np.ascontiguousarray(ay[o]), np.ascontiguousarray(aw[o])
    cap = pd.capture(leaf.lib, ax, ay, aw, leaf.poses, leaf.threads)
    lam = np.array([lam_of_region.get(int(r), 0.0) for r in leaf.reg])
    return 1.0 - cap - lam - leaf.clique_cost(leaf.poses)


def rc_check(leaf, mu, y, lam_of_region, log):
    """diagnostic: at an optimum the reduced cost of every column carrying mass must be 0.
    This catches a sign error in the region multipliers or the clique duals."""
    sup = np.nonzero(mu > 1e-9)[0]
    if not len(sup):
        return None
    sel = np.nonzero(y > 1e-12)[0]
    if not len(sel):
        return None
    ax = np.ascontiguousarray(leaf.pts[sel, 0]); ay = np.ascontiguousarray(leaf.pts[sel, 1])
    aw = np.ascontiguousarray(y[sel])
    o = np.argsort(ax, kind='stable')
    ax, ay, aw = np.ascontiguousarray(ax[o]), np.ascontiguousarray(ay[o]), np.ascontiguousarray(aw[o])
    P = [leaf.poses[i] for i in sup]
    cap = pd.capture(leaf.lib, ax, ay, aw, P, leaf.threads)
    lam = np.array([lam_of_region.get(int(leaf.region_of(p)), 0.0) for p in P])
    rc = 1.0 - cap - lam - leaf.clique_cost(P)
    return float(np.abs(rc).max())


# ====================================================================== checkpoint
def save_poses(leaf, mu, path, t):
    """`pose cx cy theta_deg mu`, i.e. the support-file format `level2_regions.py` and
    `packing_dual.read_support` already read -- but with EVERY pose, not just the support, so the
    file is a pose-set checkpoint as well as a measure."""
    with open(path + '.tmp', 'w') as f:
        f.write(f'# t4screen poses t={t} n={len(leaf.poses)} mass={float(mu.sum()):.9f}\n')
        for i, (cx, cy, th) in enumerate(leaf.poses):
            # 17 significant digits: the LP parks mass on region boundaries (cy = t/2 to within
            # 1e-14), so a checkpoint printed to fewer digits changes a pose's region on read-back
            f.write(f'pose {cx:.17g} {cy:.17g} {math.degrees(th):.17g} {mu[i]:.12g}\n')
    os.replace(path + '.tmp', path)


def load_poses(path):
    out = []
    for line in open(path):
        q = line.split()
        if q and q[0] == 'pose':
            out.append((float(q[1]), float(q[2]), math.radians(float(q[3]))))
    return out


def save_cliques(leaf, path):
    """the cliques as exact integer anchors over `leaf.D` (a superset of branch.py's parameter
    format: `A` lines carry the segment endpoints directly, so a separated clique of unknown
    provenance round-trips)."""
    with open(path + '.tmp', 'w') as f:
        f.write(f'{int(round(leaf.t * leaf.D))} {leaf.D} 0\n')
        for cl in leaf.cliques:
            (p, A) = cl[0]
            v = [int(Fr(x) * leaf.D) for x in (p[1], p[2], A[1], A[2], A[3], A[4])]
            f.write(' '.join(str(x) for x in v) + '\n')
    os.replace(path + '.tmp', path)


# ====================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('T'); ap.add_argument('TAG')
    ap.add_argument('--r', type=float, default=1.0, help='corner box / wall strip width')
    ap.add_argument('--corners', default='1111', help='mass in each corner box ("." = free, "...." = pure)')
    ap.add_argument('--patterns', default='........',
                    help='comma-separated slot patterns valued on the shared model ("........" = control)')
    ap.add_argument('--stages', type=int, default=4, help='pricing stages')
    ap.add_argument('--warm-stages', type=int, default=1, help='control-only stages before the patterns')
    ap.add_argument('--rowloops', type=int, default=14, help='row/clique generation loops per value')
    ap.add_argument('--seed-pitch', type=float, default=0.25)
    ap.add_argument('--seed-dth', type=float, default=15.0)
    ap.add_argument('--row-pitch', type=float, default=0.08)
    ap.add_argument('--price-pitch', type=float, default=0.04)
    ap.add_argument('--price-dth', type=float, default=2.5)
    ap.add_argument('--price-tol', type=float, default=1e-7)
    ap.add_argument('--cg-want', type=int, default=400)
    ap.add_argument('--pose-max', type=int, default=6000, help='cap on pose columns after pricing (0 = no cap)')
    ap.add_argument('--corner-rows', type=int, default=1200)
    ap.add_argument('--row-cap', type=int, default=20000)
    ap.add_argument('--row-age', type=int, default=0, help='drop a slack, dual-free point row after this many solves (0 = never)')
    ap.add_argument('--warm', action='append', default=[], help='support / branch-dual file(s) to seed poses from')
    ap.add_argument('--load-poses', default=None, help='t4screen pose checkpoint to restart from')
    ap.add_argument('--price-pattern', default='........',
                    help='the slot pattern whose dual drives column pricing (default: the control)')
    ap.add_argument('--method', default='highs-ipm')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--time', type=float, default=780, help='wall-clock budget (s); checkpoints then stops')
    # --- cliques
    ap.add_argument('--cliques', action='store_true', help='separate anchor cliques at all')
    ap.add_argument('--cq-wall', action='store_true', help='the wall-band family (anchorclique.kpa)')
    ap.add_argument('--cq-interior', action='store_true', help='the general segment family (anchorclique.kseg)')
    ap.add_argument('--cq-load', action='append', default=[], help='branch.py clique checkpoint as candidates')
    ap.add_argument('--cq-want', type=int, default=150, help='cliques added per separation call')
    ap.add_argument('--cq-max', type=int, default=1200, help='cap on live clique rows')
    ap.add_argument('--cq-age', type=int, default=5, help='drop a clique row after this many solves with zero dual')
    ap.add_argument('--cq-pitch', type=float, default=0.02)
    ap.add_argument('--cq-top', type=int, default=250)
    ap.add_argument('--cq-int-pitch', type=float, default=0.04)
    ap.add_argument('--cq-int-top', type=int, default=60)
    ap.add_argument('--cq-ndir', type=int, default=12)
    ap.add_argument('--cq-neps', type=int, default=6)
    ap.add_argument('--cq-nrho', type=int, default=6)
    ap.add_argument('--ktol', type=float, default=1e-6)
    ap.add_argument('--D', type=int, default=1000000)
    a = ap.parse_args()
    if a.cliques and not (a.cq_wall or a.cq_interior):
        a.cq_interior = True                                   # interior is the family that bites
    if not a.cliques and not a.cq_load:
        a.cq_wall = a.cq_interior = False

    t = float(Fr(a.T))
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f't4_{a.TAG}.log'), 'a')

    def log(m):
        print(m, flush=True); logf.write(m + '\n'); logf.flush()

    T0 = time.time()
    lib = pd.build_lib()
    leaf = CLeaf(t, a.r, lib, a.threads, log, D=a.D)
    kc = [None if c == '.' else float(c) for c in a.corners]
    patterns = a.patterns.split(',')
    log(f'# t4screen t={t} r={a.r} corners={a.corners} patterns={a.patterns} '
        f'cliques={a.cliques} wall={a.cq_wall} interior={a.cq_interior} stages={a.stages} args={vars(a)}')

    # ---- rows and columns
    g = np.arange(a.row_pitch / 2, t, a.row_pitch)
    GX, GY = np.meshgrid(g, g, indexing='ij')
    leaf.add_points(np.column_stack([GX.ravel(), GY.ravel()]))
    cand = L2.seed_poses_full(t, a.seed_pitch, a.seed_dth)
    for w in a.warm:
        cand += L2.read_measure_poses(w, t)
    if a.load_poses:
        cand += load_poses(a.load_poses)
    n = leaf.add_poses(cand)
    leaf.add_points(L2.pose_corners(leaf.poses[:a.corner_rows]))
    leaf.rfix = len(leaf.pts)              # the seed grid + square corners: never aged out
    log(f'   seed: {n} poses, {len(leaf.pts)} rows ({time.time()-T0:.0f}s)')
    for cf in a.cq_load:
        leaf.load_pool(cf, log)

    hist = []
    _hj = os.path.join(RUNS, f't4_{a.TAG}.json')
    if os.path.exists(_hj):                       # keep the trajectory across chunked restarts
        try:
            hist = json.load(open(_hj)).get('hist', [])
        except Exception:
            hist = []

    def value(kc, kw, tag):
        """row + clique generation to convergence; returns the trajectory's last record.

        Ordering matters: the sifts (which REINDEX rows and cliques) run at the TOP of the loop,
        against the previous solve's mu / y, so that after the final solve only APPENDS have
        happened.  Otherwise the returned point duals `y` no longer line up with `leaf.pts` and
        the pricer sees garbage -- which is exactly what the `rc-check` diagnostic caught
        (`rc = +1.0`, pricing gap `+8.0`, on the first row-aged runs)."""
        rec = None
        pmu = py = None
        for it in range(a.rowloops):
            ndrop = ndrop_r = 0
            if a.cliques or leaf.pool:
                ndrop = leaf.sift_cliques(a.cq_age)
            if pmu is not None and py is not None and len(py) <= len(leaf.pts):
                # rows added since that solve have no dual yet and are never stale
                ypad = np.concatenate([py, np.zeros(len(leaf.pts) - len(py))])
                ndrop_r = leaf.sift_rows(pmu, ypad, a.row_age)
            sol = leaf.solve(kc, kw, a.method)
            if sol is None:
                log(f'   [{tag}.{it}] LP infeasible/failed')
                return None
            mu, y, val, lam = sol
            pmu, py = mu, y
            M, bad, nv, ovf = leaf.certify(mu)
            kmax, ncq, best = (0.0, 0, None)
            if a.cliques or leaf.pool:
                room = max(a.cq_max - len(leaf.cliques), 0)
                kmax, ncq, best = leaf.separate(mu, a, log, want=min(a.cq_want, room))
            L = val / max(M, kmax, 1.0)
            log(f'   [{tag}.{it}] cols={len(leaf.poses)} rows={len(leaf.pts)} cq={len(leaf.cliques)} '
                f'LP={val:.6f} M={M:.9f} kmax={kmax:.6f} L={L:.6f} bad={len(bad)} +cq={ncq} -cq={ndrop} -row={ndrop_r} '
                f'(lp {leaf.lp_secs:.0f}s, {time.time()-T0:.0f}s)')
            rec = dict(tag=tag, it=it, LP=val, M=M, kmax=kmax, L=L, cols=len(leaf.poses),
                       rows=len(leaf.pts), cq=len(leaf.cliques), bad=len(bad))
            nrow = 0
            if M > 1 + 1e-9 and len(bad):
                o = np.argsort(-bad[:, 2])[:a.row_cap]
                nrow = leaf.add_points(bad[o, :2])
            if nrow == 0 and ncq == 0 and M <= 1 + 1e-9 and kmax <= 1 + a.ktol:
                break
            if nrow == 0 and ncq == 0:
                break
            if time.time() - T0 > a.time:
                log('   [time limit inside the inner loop]')
                break
        # matched pair: the same poses and the same rows with the clique constraints switched off.
        # This is the honest measurement of what the clique family is worth (RECONCILE.md 3): one
        # row set, one column set, cliques on/off.  `M_pure` is certified for that measure too.
        if leaf.cliques:
            s2 = leaf.solve(kc, kw, a.method, no_cliques=True)
            if s2 is not None:
                mu2, _, val2, _ = s2
                M2, _, _, _ = leaf.certify(mu2)
                rec['LP_pure'] = val2
                rec['M_pure'] = M2
                rec['gain'] = val2 - rec['LP']
        return rec, mu, y, lam

    def lam_map(kw, lam):
        """region index -> multiplier, in the order `solve` builds the equality rows"""
        out = {}
        idx = 0
        for i in range(4):
            if kc[i] is not None:
                out[i] = float(lam[idx]); idx += 1
        for j in range(8):
            if kw[j] is not None:
                out[4 + j] = float(lam[idx]); idx += 1
        return out

    results = {}
    ppat = a.price_pattern
    for stage in range(a.stages):
        pats = [ppat] if stage < a.warm_stages else ([ppat] + [p for p in patterns if p != ppat])
        pmu = py = plam = None
        for pat in pats:
            kw = [None if c == '.' else float(c) for c in pat]
            out = value(kc, kw, f's{stage}.{pat}')
            if out is None:
                results[pat] = None
                continue
            rec, mu, y, lam = out
            rec['stage'] = stage; rec['pattern'] = pat
            hist.append(rec)
            results[pat] = rec
            if pat == ppat:
                pmu, py, plam = mu, y, lam_map(kw, lam)
            log(f'   STAGE {stage} {pat}: LP={rec["LP"]:.6f} M={rec["M"]:.9f} kmax={rec["kmax"]:.6f} '
                f'L={rec["L"]:.6f} cols={rec["cols"]} rows={rec["rows"]} cq={rec["cq"]}'
                + (f' | matched pure LP={rec["LP_pure"]:.6f} (M={rec["M_pure"]:.6f}) gain={rec["gain"]:+.6f}'
                   if 'LP_pure' in rec else ''))
            if time.time() - T0 > a.time:
                break
        # ---- checkpoint
        if pmu is not None:
            save_poses(leaf, pmu, os.path.join(RUNS, f't4_{a.TAG}_poses.txt'), t)
        if leaf.cliques:
            save_cliques(leaf, os.path.join(RUNS, f't4_{a.TAG}_cliques.txt'))
        json.dump(dict(t=t, tag=a.TAG, args=vars(a), hist=hist), open(os.path.join(RUNS, f't4_{a.TAG}.json'), 'w'), indent=1)
        if time.time() - T0 > a.time or stage == a.stages - 1 or pmu is None:
            log('   [stopping: budget or last stage]')
            break
        # ---- pricing, against the dual of `--price-pattern`
        chk = rc_check(leaf, pmu, py, plam, log)
        if chk is not None:
            log(f'   rc-check (max |reduced cost| over columns carrying mass, must be ~0): {chk:.2e}')
        new, gap, _ = price(leaf, py, plam, a, log)
        heavy = [leaf.poses[i] for i in np.argsort(-pmu)[:100]]
        nn = leaf.add_poses(new + pd.neighbours(heavy, t))
        ndp = 0
        if a.pose_max and len(leaf.poses) > a.pose_max:
            rc_all = pose_rc(leaf, py, plam)
            keep = np.zeros(len(leaf.poses), bool)
            keep[:len(pmu)] = pmu > 1e-9                  # everything carrying mass survives
            room = a.pose_max - int(keep.sum())
            if room > 0:                                  # stratified, for the reason in `price`
                preg = np.array(leaf.reg)
                regs = sorted(set(int(v) for v in preg))
                q = max(room // max(len(regs), 1), 1)
                for rr in regs:
                    w = np.nonzero((preg == rr) & ~keep)[0]
                    if len(w):
                        keep[w[np.argsort(-rc_all[w])[:q]]] = True
                room2 = a.pose_max - int(keep.sum())
                if room2 > 0:
                    w = np.nonzero(~keep)[0]
                    keep[w[np.argsort(-rc_all[w])[:room2]]] = True
            ndp = leaf.sift_poses(keep)
        log(f'   stage {stage}: +{nn} columns, -{ndp} sifted -> {len(leaf.poses)} '
            f'(pricing gap {gap:+.6f})')
        if nn == 0:
            log('   [no improving column]')
            break

    log('RESULT tag=%s t=%s corners=%s cliques=%d ' % (a.TAG, t, a.corners, int(a.cliques))
        + ' '.join((f'{k}=LP{v["LP"]:.6f}/M{v["M"]:.6f}/K{v["kmax"]:.4f}/L{v["L"]:.6f}'
                    + (f'/pure{v["LP_pure"]:.6f}' if 'LP_pure' in v else '')) if v else f'{k}=FAIL'
                   for k, v in results.items())
        + f' poses={len(leaf.poses)} rows={len(leaf.pts)} cq={len(leaf.cliques)} secs={time.time()-T0:.0f}')


if __name__ == '__main__':
    main()
