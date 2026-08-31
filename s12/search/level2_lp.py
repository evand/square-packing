#!/usr/bin/env python3
"""Packing-side value of a level-2 (corner + wall-slot occupancy) branch leaf.

    max  sum mu_S
    s.t. coverage(p) <= 1                       for every point p of [0,t]^2
         mu({S : centre(S) in C_i}) = kc_i      i = 1..4   (corner boxes)
         mu({S : centre(S) in W_j}) = kw_j      j = 1..8   (wall slots)
         mu >= 0

over closed unit squares in [0,t]^2.  By LP duality this is exactly the cover-side leaf value
`min W - sum_i lam_i kc_i - sum_j lam_j kw_j` of a branch certificate with those regions and
per-region multipliers (`packing_le_weight_regions`), so a value `< 12` means the leaf closes and
a value `>= 12` means it does not.  The leaf is `open` in the LP sense at exactly the value 12.

The measure is NOT D4-symmetrised: an asymmetric occupancy pattern needs the full container.
(For a symmetric pattern the optimum can be symmetrised, so `--slots 11111111` etc. could also be
run in `packing_dual.py`'s reduced model; this script always uses the full one.)

Direction of the errors, which is the whole point of using it as an oracle:

  * the pose set is finite (a seed grid + warm starts + pricing), so the LP maximum over it is a
    **lower** bound on the leaf's true value -- an *optimistic* leaf value.  A leaf whose value
    here is already `>= 12` certainly does not close.
  * the rows are a finite set of points, so before row generation converges the value is inflated.
    The loop adds every violated vertex of the arrangement of the support squares (the exact
    certification of `packing_dual.py`), and stops when the certified max coverage is `<= 1 + 1e-9`;
    the printed `L` is then `sum(mu)/M` with `M` the certified maximum, which is again a genuine
    lower bound on the leaf value over the current poses.

Usage
    python3 search/level2_lp.py 3.98 TAG --slots 11110000 [--corners 1111] [--rounds 8]
    python3 search/level2_lp.py 4.0  TAG --wall 4          # total wall mass only (shared lambda)
"""
import argparse, ctypes, math, os, sys, time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
for _alt in ('/home/evand/math/square-packing/s12/search',):        # worktree fallback
    if not os.path.exists(os.path.join(HERE, 'packing_dual.py')) and os.path.exists(_alt):
        sys.path.insert(0, _alt)
import packing_dual as pd                                            # noqa: E402
from level2_regions import classify                                  # noqa: E402

RUNS = os.path.join(REPO, 'runs')
pd.RUNS = RUNS                                                       # compile the kernels into OUR runs/
TOL = pd.TOL


def _dp(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def sq_array(poses, weights=None):
    """n x 5 (cx, cy, ct, st, col) or n x 6 (+w) -- one row per pose, NO symmetry images"""
    if len(poses) == 0:
        return np.zeros((0, 6 if weights is not None else 5))
    P = np.asarray(poses, dtype=float)
    out = [P[:, 0], P[:, 1], np.cos(P[:, 2]), np.sin(P[:, 2]), np.arange(len(P), dtype=float)]
    if weights is not None:
        out.append(np.asarray(weights, dtype=float))
    return np.ascontiguousarray(np.column_stack(out))


class Leaf:
    """the full-container LP for one leaf"""

    def __init__(self, t, r, lib, threads, log):
        self.t, self.r, self.lib, self.threads, self.log = t, r, lib, threads, log
        self.poses = []
        self.key = {}
        self.reg = []                      # 0 = corner box i (0..3), 4+j = wall slot j, 12 = interior
        self.pts = np.zeros((0, 2))
        self.pkeys = set()
        self.A = None

    # ---------------------------------------------------------------- columns
    def region_of(self, p):
        k, j = classify(p[0], p[1], self.t, self.r)
        return j if k == 'C' else (4 + j if k == 'W' else 12)

    def add_poses(self, cand):
        new = []
        for c in cand:
            q = pd.clamp_pose(c[0], c[1], pd.snap_angle(c[2]), self.t)
            if q is None or not pd.admissible(*q, self.t):
                continue
            k = (round(q[0] * 1e7), round(q[1] * 1e7), round(q[2] / pd.ANG_UNIT))
            if k in self.key:
                continue
            self.key[k] = len(self.poses)
            self.poses.append(q)
            self.reg.append(self.region_of(q))
            new.append(q)
        if new and len(self.pts):
            An = self.incidence(self.pts, new)
            self.A = An if self.A is None else sp.hstack([self.A, An], format='csr')
        return len(new)

    def add_points(self, pts):
        pts = np.asarray(pts, dtype=float).reshape(-1, 2)
        keep = []
        for (x, y) in pts:
            if not (-1e-9 <= x <= self.t + 1e-9 and -1e-9 <= y <= self.t + 1e-9):
                continue
            k = (round(x * 1e9), round(y * 1e9))
            if k in self.pkeys:
                continue
            self.pkeys.add(k)
            keep.append((x, y))
        if not keep:
            return 0
        P = np.array(keep)
        if self.poses:
            An = self.incidence(P, self.poses)
            self.A = An if self.A is None else sp.vstack([self.A, An], format='csr')
        self.pts = np.vstack([self.pts, P])
        return len(P)

    def incidence(self, P, poses):
        P = np.ascontiguousarray(np.asarray(P, dtype=float))
        S = sq_array(poses)
        n = len(P)
        cnt = np.zeros(n, dtype=np.int64)
        off = np.zeros(n, dtype=np.int64)
        self.lib.incidence(n, _dp(P), len(S), _dp(S), self.t, TOL, 0,
                           cnt.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
                           off.ctypes.data_as(ctypes.POINTER(ctypes.c_long)), None, None, self.threads)
        off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
        nnz = int(off[-1])
        ri = np.zeros(max(nnz, 1), dtype=np.int32)
        ci = np.zeros(max(nnz, 1), dtype=np.int32)
        self.lib.incidence(n, _dp(P), len(S), _dp(S), self.t, TOL, 1,
                           cnt.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
                           off.ctypes.data_as(ctypes.POINTER(ctypes.c_long)),
                           ri.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                           ci.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), self.threads)
        return sp.coo_matrix((np.ones(nnz), (ri[:nnz], ci[:nnz])), shape=(n, len(poses))).tocsr()

    def drop_rows(self, keep):
        idx = np.nonzero(keep)[0]
        self.pts = self.pts[idx]
        self.A = self.A[idx]
        self.pkeys = set((round(x * 1e9), round(y * 1e9)) for x, y in self.pts)

    # ---------------------------------------------------------------- LP
    def solve(self, kc, kw, method='highs'):
        n = len(self.poses)
        reg = np.array(self.reg)
        rows, vals = [], []
        for i in range(4):
            if kc[i] is not None:
                rows.append((reg == i).astype(float)); vals.append(kc[i])
        for j in range(8):
            if kw[j] is not None:
                rows.append((reg == 4 + j).astype(float)); vals.append(kw[j])
        Aeq = np.array(rows) if rows else None
        beq = np.array(vals) if rows else None
        t0 = time.time()
        res = linprog(c=-np.ones(n), A_ub=self.A, b_ub=np.ones(len(self.pts)),
                      A_eq=Aeq, b_eq=beq, bounds=(0, None), method=method)
        self.lp_secs = time.time() - t0
        if not res.success:
            return None
        mu = np.maximum(res.x, 0.0)
        y = np.maximum(-res.ineqlin.marginals, 0.0)
        lam = -res.eqlin.marginals if rows else np.zeros(0)
        return mu, y, -res.fun, lam

    # ---------------------------------------------------------------- certification
    def certify(self, mu, thresh=1.0 + 1e-9, cap=400_000):
        sup = np.nonzero(mu > 1e-12)[0]
        S = sq_array([self.poses[i] for i in sup], weights=mu[sup])
        out = np.zeros((cap, 3))
        nout = ctypes.c_long(0)
        M = ctypes.c_double(0.0)
        ovf = ctypes.c_int(0)
        nv = self.lib.vertex_cov(len(S), _dp(S), self.t, TOL, 0, thresh, _dp(out), cap,
                                 ctypes.byref(nout), ctypes.byref(M), ctypes.byref(ovf), self.threads)
        return M.value, out[:nout.value], nv, bool(ovf.value)


# ====================================================================== seeds and pricing
def pose_corners(poses):
    """the four corners of each square (no symmetry images)"""
    S = sq_array(poses)
    out = []
    for e1 in (-0.5, 0.5):
        for e2 in (-0.5, 0.5):
            u0 = S[:, 0] * S[:, 2] + S[:, 1] * S[:, 3] + e1
            u1 = -S[:, 0] * S[:, 3] + S[:, 1] * S[:, 2] + e2
            out.append(np.column_stack([S[:, 2] * u0 - S[:, 3] * u1, S[:, 3] * u0 + S[:, 2] * u1]))
    return np.concatenate(out)


def seed_poses_full(t, pitch, dth_deg):
    out = []
    for th in np.arange(0.0, math.pi / 2 - 1e-9, math.radians(dth_deg)):
        th = pd.snap_angle(th)
        w2 = pd.wid_of(th) / 2
        lo, hi = w2 + pd.ADM, t - w2 - pd.ADM
        if hi < lo:
            continue
        k = max(int(math.floor((hi - lo) / pitch)), 1)
        g = np.linspace(lo, hi, k + 1)
        for cx in g:
            for cy in g:
                out.append((float(cx), float(cy), float(th)))
    return out


def read_measure_poses(path, t):
    """poses of a support / branch-dual file, rescaled to t, with all 8 dihedral images"""
    from level2_regions import read_measure, images
    t0, meas = read_measure(path)
    out = []
    for (cx, cy, th, m) in meas:
        s = t / t0
        for im in images(cx * s, cy * s, th, t):
            out.append(im)
    return out


def price(leaf, y, lam_of_region, pitch, dth_deg, want, log):
    """best reduced-cost poses on a grid: rc = 1 - capture(dual points) - lambda(region)"""
    sel = np.nonzero(y > 1e-9)[0]
    if len(sel) == 0:
        return []
    ax = np.ascontiguousarray(leaf.pts[sel, 0])
    ay = np.ascontiguousarray(leaf.pts[sel, 1])
    aw = np.ascontiguousarray(y[sel])
    o = np.argsort(ax, kind='stable')
    ax, ay, aw = np.ascontiguousarray(ax[o]), np.ascontiguousarray(ay[o]), np.ascontiguousarray(aw[o])
    cand = seed_poses_full(leaf.t, pitch, dth_deg)
    cap = pd.capture(leaf.lib, ax, ay, aw, cand, leaf.threads)
    reg = np.array([leaf.region_of(p) for p in cand])
    lam = np.array([lam_of_region.get(int(rr), 0.0) for rr in reg])
    rc = 1.0 - cap - lam
    idx = np.argsort(-rc)[:want]
    idx = idx[rc[idx] > 1e-7]
    best = [cand[i] for i in idx]
    ref = pd.refine(leaf.lib, ax, ay, aw, best, leaf.t, leaf.threads)
    log(f'   pricing: {len(cand)} candidates, best rc {rc.max():+.4f}, kept {len(best)}')
    return list(best) + list(ref)


# ====================================================================== main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('T')
    ap.add_argument('TAG')
    ap.add_argument('--r', type=float, default=1.0, help='corner box / wall strip width')
    ap.add_argument('--corners', default='1111', help='mass in each of the four corner boxes ("." = free)')
    ap.add_argument('--slots', default=None, help='mass in each of the eight wall slots ("." = free)')
    ap.add_argument('--wall', type=float, default=None, help='total wall mass instead of per-slot')
    ap.add_argument('--rounds', type=int, default=6, help='pricing rounds')
    ap.add_argument('--rowloops', type=int, default=12, help='row-generation loops per round')
    ap.add_argument('--seed-pitch', type=float, default=0.10)
    ap.add_argument('--seed-dth', type=float, default=7.5)
    ap.add_argument('--row-pitch', type=float, default=0.05)
    ap.add_argument('--price-pitch', type=float, default=0.04)
    ap.add_argument('--price-dth', type=float, default=2.5)
    ap.add_argument('--cg-want', type=int, default=400)
    ap.add_argument('--warm', action='append', default=[])
    ap.add_argument('--warm-neigh', type=int, default=0, help='perturb this many warm poses as extra columns')
    ap.add_argument('--corner-rows', type=int, default=1200, help='seed rows: the corners of this many seed squares')
    ap.add_argument('--method', default='highs')
    ap.add_argument('--patterns', default=None,
                    help='comma-separated slot patterns to value on the same model, e.g. 11111111,10101010')
    ap.add_argument('--row-cap', type=int, default=20000)
    ap.add_argument('--prune-rows', type=int, default=0, help='prune slack rows once the row count exceeds this')
    ap.add_argument('--warm-rounds', type=int, default=3, help='control-only pricing rounds before the leaves')
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--time', type=float, default=3600)
    a = ap.parse_args()

    t = float(eval(a.T)) if '/' in a.T else float(a.T)
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'l2_{a.TAG}.log'), 'w')

    def log(m):
        print(m, flush=True)
        logf.write(m + '\n')
        logf.flush()

    t0 = time.time()
    lib = pd.build_lib()
    leaf = Leaf(t, a.r, lib, a.threads, log)

    kc = [None if c == '.' else float(c) for c in a.corners]
    if a.slots:
        kw = [None if c == '.' else float(c) for c in a.slots]
    else:
        kw = [None] * 8
    log(f'# level-2 leaf  t={t}  r={a.r}  corners={a.corners}  slots={a.slots or "-"}  wall={a.wall}')

    # ---- seeds
    g = np.arange(a.row_pitch / 2, t, a.row_pitch)
    GX, GY = np.meshgrid(g, g, indexing='ij')
    n = leaf.add_points(np.column_stack([GX.ravel(), GY.ravel()]))
    log(f'   seed rows {n}')
    cand = seed_poses_full(t, a.seed_pitch, a.seed_dth)
    for w in a.warm:
        wp = read_measure_poses(w, t)
        cand += wp
        if a.warm_neigh:
            cand += pd.neighbours(wp[:a.warm_neigh], t)
    n = leaf.add_poses(cand)
    log(f'   seed poses {n}')
    leaf.add_points(pose_corners(leaf.poses[:a.corner_rows]))
    log(f'   rows after square corners {len(leaf.pts)}')

    def leaf_value(kc, kw, rowloops, tag):
        """row-generate to a certified measure and return (L, LP value, M, duals)"""
        val = L = M = None
        for it in range(rowloops):
            sol = leaf.solve(kc, kw, a.method)
            if sol is None:
                log(f'   [{tag}.{it}] LP infeasible/failed')
                return None
            mu, y, val, lam = sol
            M, bad, nv, ovf = leaf.certify(mu)
            L = val / max(M, 1.0)
            log(f'   [{tag}.{it}] cols={len(leaf.poses)} rows={len(leaf.pts)} LP={val:.6f} '
                f'M={M:.9f} L={L:.6f} bad={len(bad)} (lp {leaf.lp_secs:.0f}s, {time.time()-t0:.0f}s)')
            if a.prune_rows and len(leaf.pts) > a.prune_rows:
                slack = 1.0 - leaf.A @ mu
                keep = (y > 1e-12) | (slack < 0.02)
                if 200 < keep.sum() < len(leaf.pts):
                    leaf.drop_rows(keep)
                    log(f'        pruned rows -> {len(leaf.pts)}')
            if M <= 1 + 1e-9:
                break
            o = np.argsort(-bad[:, 2])[:a.row_cap]
            leaf.add_points(bad[o, :2])
        return L, val, M, mu, y, lam

    patterns = [p for p in (a.patterns.split(',') if a.patterns else [a.slots or '........'])]
    results = {}
    # ---- warm-up: control only, row generation + pricing, to build a decent pose set
    for wr in range(a.warm_rounds):
        r = leaf_value(kc, [None] * 8, a.rowloops, f'w{wr}.ctl')
        if r is None:
            return
        Lc, valc, Mc, mu, y, lam = r
        log(f'   warm {wr}: control L = {Lc:.6f} ({len(leaf.poses)} poses, {len(leaf.pts)} rows)')
        lam_of = {}
        idx = 0
        for i in range(4):
            if kc[i] is not None:
                lam_of[i] = lam[idx]; idx += 1
        new = price(leaf, y, lam_of, a.price_pitch, a.price_dth, a.cg_want, log)
        nn = leaf.add_poses(new + pd.neighbours([leaf.poses[i] for i in np.argsort(-mu)[:100]], t))
        log(f'   warm {wr}: +{nn} columns')
        if time.time() - t0 > a.time:
            break
    for rd in range(a.rounds):
        # the control (corner constraints only) fixes the row set for this round
        r = leaf_value(kc, [None] * 8, a.rowloops, f'{rd}.ctl')
        if r is None:
            return
        Lc, valc, Mc, mu, y, lam = r
        results['control'] = Lc
        log(f'   round {rd} CONTROL L = {Lc:.6f}  ({len(leaf.poses)} poses, {len(leaf.pts)} rows)')
        for pat in patterns:
            kwp = [None if c == '.' else float(c) for c in pat]
            rr = leaf_value(kc, kwp, a.rowloops, f'{rd}.{pat}')
            if rr is None:
                results[pat] = None
                continue
            results[pat] = rr[0]
            log(f'   round {rd} LEAF {pat}: L = {rr[0]:.6f}   gap to control = {Lc - rr[0]:+.6f}')
        log('   round %d summary  ' % rd + '  '.join(
            f'{k}={v:.4f}' if v is not None else f'{k}=INF' for k, v in results.items()))
        if time.time() - t0 > a.time or rd == a.rounds - 1:
            break
        # ---- pricing against the control dual
        lam_of = {}
        idx = 0
        for i in range(4):
            if kc[i] is not None:
                lam_of[i] = lam[idx]; idx += 1
        new = price(leaf, y, lam_of, a.price_pitch, a.price_dth, a.cg_want, log)
        nn = leaf.add_poses(new + pd.neighbours([leaf.poses[i] for i in np.argsort(-mu)[:100]], t))
        log(f'   round {rd}: +{nn} columns')
        if nn == 0:
            break

    log('RESULT tag=%s t=%s corners=%s ' % (a.TAG, t, a.corners) + ' '.join(
        f'{k}={v:.6f}' if v is not None else f'{k}=INFEASIBLE' for k, v in results.items())
        + f' poses={len(leaf.poses)} rows={len(leaf.pts)} secs={time.time()-t0:.0f}')


if __name__ == '__main__':
    main()
