#!/usr/bin/env python3
"""A3 (far field): the fractional-packing LP at container side t = T with a CAP on near-axis mass.

    V(T, eps, K) = sup  mu(Pi)
                   s.t. mu({Q : p in Q}) <= 1           for every p in [0,T]^2      (coverage)
                        mu({Q : tilt(Q) < eps}) <= K                                 (the cap)
                        mu >= 0

`tilt` of a pose is the distance of its angle to 0 mod 90 deg.  Lemma A3 of
`notes/proof-architecture.md` asks whether V(T, eps, T-1) < n = T^2 - T, i.e. whether capping the
near-axis mass at T-1 makes the relaxation drop below the number of squares.  Here K = T-1 and
n = T^2 - T by default (T = 4: K = 3, n = 12;  T = 3: K = 2, n = 6).

WHICH SIDE EVERY NUMBER IS ON (`notes/review-2026-09-13b.md` "Clarifications"):

  * a measure supported on a FINITE pose set that satisfies the coverage rows and the cap is
    feasible for the continuum problem, so its mass is a LOWER bound on V -- the packing side.
    Restricting the pose set can only LOWER the LP value, so every number produced by `cg`/`exact`
    is a lower bound on V.  A value >= n therefore REFUTES A3 at that eps; a value < n proves
    nothing at all (the pose set may simply be too poor).
  * only a weighted point set y >= 0 and a scalar lam >= 0 with  y(S) >= 1 for every admissible
    pose of tilt >= eps  and  y(S) >= 1 - lam for every near-axis pose  is an upper bound
    (= |y| + K lam), and that needs every pose in the continuum, not a sample.  `dualcost` reports
    such a cost from the LP dual, but its verification here is an angle-SAMPLED sweep, so it is a
    HEURISTIC upper bound, never a certificate.

Modes
    hist FILE...                     tilt histogram of a support file (float or exact), and the
                                     "restrict-and-refill" bound mass(tilt >= eps) + min(near, K)
    cg T TAG --eps E [--cap K]       capped column generation (reuses search/packing_dual.py:
                                     its C kernels, Model, sweep pricer, refiner and the
                                     arrangement-vertex certification), -> runs/arch_farfield_TAG*
    exact --src FILE --t T --eps E   exact snap of a float support + capped LP over the EXACT
                                     arrangement rows + exact certificate (reuses
                                     search/dual_exact.py), -> runs/arch_farfield_TAG_exact.txt

Nothing in this file edits or imports anything that writes to the repo's existing outputs.
"""
import sys, os, math, time, json, argparse, ctypes
from fractions import Fraction as Fr
from math import gcd, lcm

import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
sys.path.insert(0, HERE)
import packing_dual as pdl                     # float engine (C kernels, Model, pricers)
import dual_exact as dex                       # exact arrangement / certification


# ============================================================================== tilt
def tilt_deg(theta_deg):
    """distance of an angle (degrees) to 0 mod 90"""
    a = abs(theta_deg) % 90.0
    return min(a, 90.0 - a)


def sin_eps_rational(eps_deg, up=True):
    """exact rational bound on sin(eps): rounded UP, so {tilt < eps} is over-estimated, which makes
    the cap constraint stronger and the resulting value a valid lower bound on V(T, eps, K)."""
    D = 10 ** 12
    s = math.sin(math.radians(eps_deg))
    n = math.ceil(s * D) if up else math.floor(s * D)
    return Fr(n, D)


def near_axis_exact(p, q, se):
    """exact: is the pose with theta = 2 arctan(p/q) of tilt < eps?   sin(tilt) = min(|cos|,|sin|)"""
    a, b, r = dex.cos_sin(p, q)
    return min(abs(a), abs(b)) * se.denominator <= se.numerator * r


# ============================================================================== support files
def read_support_any(path):
    """-> (t, [(cx, cy, theta_deg, mass)]).  Accepts packing_dual float supports and dual_exact
    exact supports (p q cx cy mass)."""
    t = None; out = []
    for line in open(path):
        if line.startswith('#'):
            if '# t=' in line or '# t =' in line:
                try: t = float(Fr(line.split('=')[1].split()[0]))
                except Exception: pass
            continue
        w = line.split('#')[0].split()
        if not w or w[0] != 'pose': continue
        if len(w) == 5:                                     # float: cx cy theta_deg mu
            out.append((float(w[1]), float(w[2]), float(w[3]), float(w[4])))
        elif len(w) == 6:                                   # exact: p q cx cy mass
            p, q = int(w[1]), int(w[2])
            th = math.degrees(2 * math.atan2(p, q))
            out.append((float(Fr(w[3])), float(Fr(w[4])), th, float(Fr(w[5]))))
    return t, out


def cmd_hist(a):
    for path in a.files:
        t, poses = read_support_any(path)
        if not poses: print(f"{path}: no poses"); continue
        tl = np.array([tilt_deg(th) for _, _, th, _ in poses])
        mu = np.array([m for _, _, _, m in poses])
        print(f"\n{path}:  t = {t}, {len(poses)} poses, mass {mu.sum():.6f}")
        edges = [0, 1e-9, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 45.001]
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = (tl >= lo) & (tl < hi)
            if s.any(): print(f"   tilt [{lo:7.3f},{hi:7.3f}) deg: mass {mu[s].sum():9.6f}  ({100*mu[s].sum()/mu.sum():5.1f}%)  poses {s.sum()}")
        K = a.cap
        print(f"   restrict-and-refill lower bounds (drop the near-axis part, put back at most K = {K}):")
        for eps in a.eps:
            near = mu[tl < eps].sum(); far = mu[tl >= eps].sum()
            print(f"     eps = {eps:5.2f} deg: near {near:9.6f}  far {far:9.6f}  ->  far + min(near,K) = {far + min(near, K):9.6f}")


def cmd_table(a):
    """markdown table of the exact results in runs/arch_farfield_<TAG>_exact.json"""
    import json as J
    rows = {}; caps = []
    for tag in a.tags:
        d = J.load(open(os.path.join(RUNS, f'arch_farfield_{tag}_exact.json')))
        n = d['n']
        for r in d['results']:
            c = 'none' if r['cap'] is None else str(int(r['cap']))
            if c not in caps: caps.append(c)
            key = r['eps']
            prev = rows.get((key, c))
            if prev is None or r['mass_float'] > prev[0]:
                rows[(key, c)] = (r['mass_float'], tag, r['near_float'], r['M_float'])
    caps = sorted(caps, key=lambda c: (c == 'none', float(c) if c != 'none' else 0))
    print('| `eps` | ' + ' | '.join(f'cap `{c}`' for c in caps) + ' |')
    print('|---|' + '---|' * len(caps))
    for eps in sorted(set(k[0] for k in rows)):
        cells = []
        for c in caps:
            v = rows.get((eps, c))
            cells.append('—' if v is None else
                         (f"**{v[0]:.6f}**" if v[0] >= n else f"{v[0]:.6f}") + f" ({v[1]})")
        print(f'| `{eps}°` | ' + ' | '.join(cells) + ' |')


def cmd_bands(a):
    """K1 (the transversal chord lemma, `search/S6_SKELETON.md` §2.1) says at most T-1 squares of a
    packing share a point of their D-bands,  D(t) = (|cos t|+|sin t|)/2 - |cos t sin t|.  So
    mu(band_a) <= T-1 is a SOUND row for the relaxation.  This reports max_a mu(band_a) for a
    measure: if it is <= T-1 the measure survives the strengthened relaxation as well."""
    t, poses = read_support_any(a.support)
    t = float(Fr(a.t)) if a.t else t
    ev = {'y': [], 'x': []}
    for cx, cy, th, m in poses:
        for (x, y, u) in pdl.pose_images(cx, cy, math.radians(th), t):
            C = abs(math.cos(u)); S = abs(math.sin(u)); D = (C + S) / 2 - C * S
            ev['y'].append((y - D, m / 8)); ev['y'].append((y + D, -m / 8))
            ev['x'].append((x - D, m / 8)); ev['x'].append((x + D, -m / 8))
    for k in ('y', 'x'):
        e = sorted(ev[k]); cur = 0.0; best = 0.0; at = None
        for pos, d in e:
            cur += d
            if cur > best: best = cur; at = pos
        print(f"  max over {k}-band heights a of mu(band_a) = {best:.6f} at a = {at:.6f} "
              f"(K1 allows {round(t)-1} for a packing)")
    # region rows that the corner/wall branches of leaf_ceiling.py use, all sound for a packing:
    # at most one square has its centre in a corner box [0,1]^2 (BRANCH.md), and at most T-1 have
    # their centre in a wall strip [0,1] x [0,t] (K1 at a = 1/2).
    Tn = round(t); corner = [0.0] * 4; strip = [0.0] * 4; interior = 0.0
    for cx, cy, th, m in poses:
        for (x, y, u) in pdl.pose_images(cx, cy, math.radians(th), t):
            w = m / 8
            for i, (a_, b_) in enumerate(((0, 0), (1, 0), (0, 1), (1, 1))):
                if (x <= 1 if a_ == 0 else x >= t - 1) and (y <= 1 if b_ == 0 else y >= t - 1):
                    corner[i] += w
            for j, c in enumerate((x <= 1, x >= t - 1, y <= 1, y >= t - 1)):
                if c: strip[j] += w
            if 1 < x < t - 1 and 1 < y < t - 1: interior += w
    print(f"  corner-box masses (each <= 1 for a packing): " + ", ".join(f"{c:.6f}" for c in corner))
    print(f"  wall-strip masses (each <= {Tn-1} for a packing): " + ", ".join(f"{c:.6f}" for c in strip))
    print(f"  mass with centre in the open interior: {interior:.6f}")


def cmd_rotseed(a):
    """Rigidly rotate a support about the container centre by alpha.  A rigid motion preserves
    coverage exactly, so the sub-family that still fits in [0,t]^2 is a feasible measure whose
    angles have all been shifted by alpha -- i.e. the near-axis mass of the original sits at tilt
    exactly alpha.  This is the natural far-field candidate and is used as seed columns."""
    t = float(Fr(a.t))
    tt, poses = read_support_any(a.src)
    if tt and abs(tt - t) > 1e-9:
        lam = t / tt; poses = [(cx * lam, cy * lam, th, m) for cx, cy, th, m in poses]
    lines = []; kept = 0.0; tot = sum(m for *_, m in poses)
    for alpha in a.alpha:
        for sgn in (1, -1):
            A = math.radians(alpha * sgn); ca, sa = math.cos(A), math.sin(A)
            k = 0.0
            for cx, cy, th, m in poses:
                x = cx - t / 2; y = cy - t / 2
                nx = x * ca - y * sa + t / 2; ny = x * sa + y * ca + t / 2
                nth = math.radians(th) + A
                w2 = pdl.wid_of(nth) / 2
                if w2 <= nx <= t - w2 and w2 <= ny <= t - w2:
                    lines.append(f"pose {nx:.13f} {ny:.13f} {math.degrees(nth):.11f} {m:.12f}")
                    k += m
            print(f"  alpha = {alpha*sgn:+7.3f} deg: {k:.6f} of {tot:.6f} mass still fits "
                  f"({100*k/tot:.1f}%)")
            kept = max(kept, k)
    with open(a.out, 'w') as f:
        f.write(f"# t={t} rotated seeds from {a.src}, alphas {a.alpha} deg (both signs)\n")
        f.write('\n'.join(lines) + '\n')
    print(f"wrote {a.out}: {len(lines)} poses; best single-rotation surviving mass {kept:.6f}")


# ============================================================================== capped float LP
class CapModel(pdl.Model):
    """packing_dual.Model + the row  sum_{tilt(pose) < eps} mu <= K."""
    def __init__(self, t, lib, threads, log, eps_deg, cap):
        super().__init__(t, lib, threads, log)
        self.eps = eps_deg; self.cap = float(cap)

    def near_mask(self):
        # canon_pose() folds the angle into [0, pi/4], so the stored angle IS the tilt
        return np.array([math.degrees(p[2]) < self.eps - 1e-12 for p in self.poses], dtype=bool)

    def solve(self):
        n = len(self.poses); m = len(self.pts)
        near = self.near_mask().astype(float)
        cap_row = sp.csr_matrix(near.reshape(1, n))
        A = sp.vstack([self.A, cap_row], format='csr')
        b = np.concatenate([np.ones(m), [self.cap]])
        res = linprog(c=-np.ones(n), A_ub=A, b_ub=b, bounds=(0, None), method='highs')
        if not res.success: return None
        mu = np.maximum(res.x, 0.0)
        marg = np.maximum(-res.ineqlin.marginals, 0.0)
        return mu, marg[:m], float(marg[m]), -res.fun


def capped_value(mass, near, M, cap):
    """the certified value of a float measure: scale by 1/M (coverage <= 1 exactly) and then shrink
    the near-axis part back down to the cap.  Both operations keep coverage <= 1."""
    if M <= 0: return 0.0
    return (mass - near) / M + min(near / M, cap)


def seed_poses_at(t, pitch, angles_rad):
    out = []
    for th in angles_rad:
        th = pdl.snap_angle(th); w2 = pdl.wid_of(th) / 2
        k = int(math.floor((t / 2 - w2 - pdl.ADM) / pitch))
        g = t / 2 + pitch * np.arange(-k, k + 1)
        g = np.unique(np.concatenate([g, [w2 + pdl.ADM, t - w2 - pdl.ADM]]))
        for cx in g:
            for cy in g: out.append((float(cx), float(cy), th))
    return out


def write_support(path, t, poses, mu, pts, y, header):
    with open(path, 'w') as f:
        f.write(header + '\n')
        f.write("# D4-symmetrised measure: pose cx cy theta_deg mu  (mass mu/8 on each of the 8 dihedral images)\n")
        for i in np.nonzero(mu > 1e-9)[0]:
            cx, cy, th = poses[i]
            f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} {mu[i]:.12f}\n")
        f.write("# constraint points (fundamental domain 0<=x<=y<=t/2) with dual weight y\n")
        for i in range(len(pts)):
            f.write(f"point {pts[i,0]:.10f} {pts[i,1]:.10f} {y[i]:.10f}\n")


def cmd_cg(a):
    t = float(Fr(a.T)); tag = a.TAG
    os.makedirs(RUNS, exist_ok=True)
    lf = open(os.path.join(RUNS, f'arch_farfield_{tag}.log'), 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    cap = a.cap if a.cap is not None else round(t) - 1
    n_target = round(t) ** 2 - round(t)
    log(f"\narch_farfield cg t={t} tag={tag} eps={a.eps} deg cap={cap} n={n_target} args={vars(a)}")
    lib = pdl.build_lib(); T0 = time.time()
    m = CapModel(t, lib, a.threads, log, a.eps, cap)

    seeds = pdl.seed_poses(t, a.seed_pitch, a.seed_dth)
    # extra angles hugging the band edge (the LP wants tilt just above eps) and the near-axis side
    extra = [math.radians(x) for x in
             (a.eps, a.eps * 1.02 + 0.02, a.eps * 1.2, a.eps * 1.5, a.eps * 2, a.eps * 3,
              max(a.eps * 0.5, 0.0), 0.0)] if a.eps > 0 else []
    extra = sorted(set(e for e in extra if 0 <= e <= math.pi / 4))
    seeds += seed_poses_at(t, a.seed_pitch, extra)
    if a.seed_file:
        for path in a.seed_file:
            tt, ps = read_support_any(path)
            lam = t / tt if tt else 1.0
            seeds += [(cx * lam, cy * lam, math.radians(th)) for cx, cy, th, _ in ps]
            log(f"seed file {path}: {len(ps)} poses (centres scaled by {lam:.6f})")
    m.add_poses(seeds)
    m.add_points(pdl.seed_points(t, a.row_pitch))
    m.add_points(pdl.corners(m.poses, t))
    log(f"seed: {len(m.poses)} poses ({int(m.near_mask().sum())} near-axis), {len(m.pts)} rows, nnz={m.A.nnz} ({time.time()-T0:.0f}s)")

    angs = np.unique(np.concatenate([
        np.arange(0.0, math.pi / 4 + 1e-9, math.radians(a.raster_dth)),
        np.array(extra) if extra else np.zeros(0),
        np.arange(math.radians(a.eps), math.radians(min(a.eps * 4, 45.0)) + 1e-9,
                  math.radians(max(a.eps / 4, 0.05))) if a.eps > 0 else np.zeros(0)]))
    hist = []; best = dict(L=0.0)
    for rd in range(a.rounds):
        t0 = time.time(); out = m.solve()
        if out is None: log("LP failed"); break
        mu, y, lam, mass = out; tlp = time.time() - t0
        near_mask = m.near_mask(); near = float(mu[near_mask].sum())
        t0 = time.time(); M, bad, nv, ovf = m.certify(mu, cap=2_000_000); tc = time.time() - t0
        L = capped_value(mass, near, M, cap)
        nsup = int((mu > 1e-12).sum())
        log(f"round {rd}: LP mass={mass:.6f} near={near:.6f} lam={lam:.6f} M={M:.6f} -> L={L:.6f} "
            f"| support={nsup} cols={len(m.poses)} rows={len(m.pts)} vertices={nv} bad={len(bad)} "
            f"| lp {tlp:.0f}s cert {tc:.0f}s total {time.time()-T0:.0f}s")
        hdr = (f"# t={t} round={rd} eps={a.eps} cap={cap} mass={mass:.9f} near={near:.9f} "
               f"M={M:.9f} L={L:.9f}")
        write_support(os.path.join(RUNS, f'arch_farfield_{tag}_last_support.txt'), t, m.poses, mu, m.pts, y, hdr)
        if L > best['L']:
            best = dict(L=L, M=M, mass=mass, near=near, lam=lam, round=rd, support=nsup)
            write_support(os.path.join(RUNS, f'arch_farfield_{tag}_support.txt'), t, m.poses, mu, m.pts, y, hdr)

        n_old_cols = len(mu); n_old_rows = len(y)
        slack = 1.0 - np.asarray(m.A @ mu).ravel()
        nrow = m.add_points(bad[np.argsort(bad[:, 2])[::-1][:a.row_cap], :2]) if len(bad) else 0

        ncol = 0; rmin = None
        polish = (rd >= a.polish_after) or (time.time() - T0 > a.time)
        if not polish:
            ax, ay, aw = pdl.atoms_from_dual(m, y, ymin=1e-12)
            sx, sy, sww = pdl.atoms_from_dual(m, y, ymin=1e-12, maxrows=a.sweep_rows)
            cand = []
            if len(sx):
                cand += [(cx, cy, th) for cx, cy, th, v in
                         pdl.price_sweep(m, lib, sx, sy, sww, angs, a.sweep_k, 0.05, a.threads)]
            sup_poses = [m.poses[i] for i in np.nonzero(mu > 1e-9)[0]]
            cand += pdl.neighbours(sup_poses, t)
            if len(ax) == 0: cand = []
            allc = []
            if cand:
                v = pdl.capture(lib, ax, ay, aw, cand, a.threads)
                allc += [(cx, cy, th, vv) for (cx, cy, th), vv in zip(cand, v)]
                allc += pdl.refine(lib, ax, ay, aw, cand, t, a.threads)
            # reduced cost: 1 - capture - lam*[near-axis] > 0
            scored = []
            for cx, cy, th, v in allc:
                p = pdl.clamp_pose(cx, cy, th, t)
                if p is None: continue
                p = pdl.canon_pose(*p, t)
                rc = 1.0 - v - (lam if math.degrees(p[2]) < a.eps - 1e-12 else 0.0)
                if rc > 1e-7: scored.append((rc, p))
            scored.sort(key=lambda r: -r[0])
            rmin = scored[0][0] if scored else None
            ncol = m.add_poses([p for _, p in scored[:a.cg_want]])
            log(f"   candidates {len(allc)} improving {len(scored)} best rc {rmin}")

        for i in range(n_old_cols): m.col_zero[i] = m.col_zero[i] + 1 if mu[i] <= 1e-12 else 0
        keepc = np.array([cz < a.col_age for cz in m.col_zero], dtype=bool)
        if (~keepc).sum() > 0: m.drop_columns(keepc)
        stale = (y <= 1e-12) & (slack > 0.02)
        m.row_zero[:n_old_rows] = np.where(stale, m.row_zero[:n_old_rows] + 1, 0)
        keepr = m.row_zero < a.row_age
        if (~keepr).sum() > 0: m.drop_rows(keepr)
        hist.append(dict(round=rd, mass=mass, near=near, lam=lam, M=M, L=L, cols=len(m.poses),
                         rows=len(m.pts), new_rows=nrow, new_cols=ncol, best_rc=rmin, t=time.time() - T0))
        json.dump(dict(t=t, tag=tag, eps=a.eps, cap=cap, n=n_target, args=vars(a), best=best, hist=hist),
                  open(os.path.join(RUNS, f'arch_farfield_{tag}.json'), 'w'), indent=1)
        log(f"   +rows {nrow} +cols {ncol} -> cols={len(m.poses)} rows={len(m.pts)} ({time.time()-t0:.0f}s)")
        if time.time() - T0 > a.time + a.polish_time: log("time limit"); break
        if nrow == 0 and ncol == 0:
            log("converged: no violated vertex" + ("" if polish else ", no improving column")); break
    log(f"BEST t={t} eps={a.eps} cap={cap}: L={best['L']:.6f} (n = {n_target}; "
        f"{'>= n -- A3 REFUTED at this eps' if best['L'] >= n_target - 1e-9 else '< n so far -- inconclusive, it is only a lower bound'}) "
        f"round {best.get('round')} mass {best.get('mass')} near {best.get('near')} M {best.get('M')} wall {time.time()-T0:.0f}s")


# ============================================================================== exact
def cmd_exact(a):
    t = Fr(a.t); dex.set_t(t, a.tag)
    dex.LOG = open(os.path.join(RUNS, f'arch_farfield_{a.tag}_exact.log'), 'w')
    log = dex.log
    T0 = time.time()
    Tint = int(round(float(t)))
    n_target = Tint ** 2 - Tint
    src = dex.read_float_support(a.src)
    log(f"exact: t = {t}, caps {a.cap}, n = {n_target}; eps list {a.eps} deg")
    log(f"  {len(src)} float poses from {a.src}; snap Q={a.Q}, Dc={a.Dc}")
    poses = [dex.snap_pose(cx, cy, th, a.Q, a.Dc) for cx, cy, th, mu in src]
    seen = {}; keep = []
    for i, p in enumerate(poses):                       # dedupe after snapping
        if p not in seen: seen[p] = len(keep); keep.append(p)
    poses = keep
    log(f"  {len(poses)} distinct snapped poses")
    squares, pose_of = dex.build_squares(poses)
    log(f"  {len(squares)} images, all admissible (exact)")
    dex.enumerate_vertices(squares, full=False, procs=a.procs)
    inc = dex.incidences(a.procs)
    pat = dex.patterns(inc, pose_of)
    log(f"  {len(pat)} distinct incidence patterns (LP rows)")
    keys = list(pat.keys()); nP = len(poses)
    rows = []; cols = []; vals = []
    for i, key in enumerate(keys):
        for k, c in key: rows.append(i); cols.append(k); vals.append(c / 8.0)
    A0 = sp.csr_matrix((vals, (rows, cols)), shape=(len(keys), nP))
    DM = a.DM
    out_all = []
    caps = a.cap if a.cap else [Tint - 1]
    for eps in a.eps:
        se = sin_eps_rational(eps, up=True)
        near = [near_axis_exact(p, q, se) for (p, q, _, _) in poses]
        log(f"\n  --- eps = {eps} deg (sin(eps) <= {se}; the near-axis set is a SUPERSET, so the cap "
            f"is if anything too strong): {sum(near)} of {nP} poses near-axis")
        nearrow = sp.csr_matrix(np.array([1.0 if z else 0.0 for z in near]).reshape(1, nP))
        for cap in caps:
            uncapped = (eps <= 0) or (cap >= nP * 100)
            if uncapped:
                A = A0; b = np.ones(len(keys))
            else:
                A = sp.vstack([A0, nearrow], format='csr')
                b = np.concatenate([np.ones(len(keys)), [float(cap)]])
            t0 = time.time()
            res = linprog(-np.ones(nP), A_ub=A, b_ub=b, bounds=(0, None), method='highs',
                          options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
            assert res.status == 0, res.message
            MU = [int(math.floor(float(x) * DM)) for x in np.maximum(res.x, 0.0)]
            M, key = dex.exact_max(pat, MU, DM)
            if M > 1:
                MU = [(x * M.denominator) // M.numerator for x in MU]
                M, key = dex.exact_max(pat, MU, DM)
            assert M <= 1
            NEAR = sum(x for x, z in zip(MU, near) if z)
            if not uncapped and Fr(NEAR, DM) > cap:
                num, den = Fr(int(cap * DM), NEAR).as_integer_ratio() if NEAR else (1, 1)
                MU = [((x * num) // den if z else x) for x, z in zip(MU, near)]
                M, key = dex.exact_max(pat, MU, DM)
                NEAR = sum(x for x, z in zip(MU, near) if z)
            mass = Fr(sum(MU), DM); nearm = Fr(NEAR, DM)
            assert M <= 1 and (uncapped or nearm <= cap)
            log(f"    cap = {'none' if uncapped else cap}: LP {-res.fun:.9f} ({time.time()-t0:.0f} s) "
                f"-> EXACT mass = {float(mass):.9f}, near-axis {float(nearm):.6f}, M = {float(M):.12f} <= 1 "
                f"=>  V(t={t}, eps={eps}, cap={'none' if uncapped else cap}) >= {float(mass):.9f} "
                f"({'>=' if mass >= n_target else '<'} n = {n_target})")
            tag2 = f"{a.tag}_e{str(eps).replace('.','p')}_c{cap}"
            out = os.path.join(RUNS, f'arch_farfield_{tag2}_exact_support.txt')
            with open(out, 'w') as f:
                f.write(f"# t = {dex.TN}/{dex.TD} exact; D4-symmetrised measure of closed unit squares; "
                        f"eps = {eps} deg (sin(eps) <= {se}), cap = {'none' if uncapped else cap}\n")
                f.write(f"# mass = {mass} = {float(mass):.12f};  near-axis mass = {nearm} = {float(nearm):.9f};  "
                        f"M = max coverage = {M} = {float(M):.15f}\n")
                f.write("# pose p q cx cy mass : theta = 2*arctan(p/q); mass/8 on each of the 8 dihedral images\n")
                for (pp, qq, cx, cy), x, z in zip(poses, MU, near):
                    if x > 0: f.write(f"pose {pp} {qq} {cx} {cy} {Fr(x, DM)}   {'# near-axis' if z else ''}\n")
            out_all.append(dict(eps=eps, cap=None if uncapped else cap, mass=str(mass), mass_float=float(mass),
                                near=str(nearm), near_float=float(nearm), M=str(M), M_float=float(M),
                                lp=float(-res.fun), poses=sum(1 for x in MU if x > 0), sin_eps=str(se),
                                support=out))
            json.dump(dict(t=str(t), n=n_target, src=a.src, poses=nP, rows=len(keys), results=out_all,
                           seconds=time.time() - T0), open(os.path.join(RUNS, f'arch_farfield_{a.tag}_exact.json'), 'w'), indent=1)
    log(f"\n  total {time.time()-T0:.1f} s")


# ============================================================================== dual cost (cover side)
def cmd_dualcost(a):
    """heuristic (angle-SAMPLED) upper bound from a run's dual: |y| + K*lam, valid only if the
    sampled minimum capture is >= 1 over tilt >= eps and >= 1 - lam over tilt < eps."""
    t, _ = read_support_any(a.support)
    pts = []; ys = []
    for line in open(a.support):
        w = line.split()
        if w and w[0] == 'point':
            pts.append((float(w[1]), float(w[2]))); ys.append(float(w[3]))
    pts = np.array(pts); ys = np.array(ys)
    lib = pdl.build_lib()
    m = pdl.Model(t, lib, a.threads, print); m.pts = pts
    ax, ay, aw = pdl.atoms_from_dual(m, ys, ymin=1e-15)
    tot = float(aw.sum())
    angs_far = np.arange(math.radians(a.eps), math.pi / 4 + 1e-9, math.radians(a.dth))
    angs_near = np.arange(0.0, math.radians(a.eps), math.radians(max(a.eps / 40, 1e-3)))
    def minsweep(angs):
        if len(angs) == 0: return None
        sw = pdl.price_sweep(m, lib, ax, ay, aw, np.ascontiguousarray(angs), 1, 0.02, a.threads)
        return min([v for _, _, _, v in sw], default=None)
    mf = minsweep(angs_far); mn = minsweep(angs_near)
    cap = a.cap if a.cap is not None else round(t) - 1
    print(f"{a.support}: t={t} total dual weight |y| = {tot:.6f}")
    print(f"  sampled min capture over tilt >= {a.eps} deg ({len(angs_far)} angles): {mf}")
    print(f"  sampled min capture over tilt <  {a.eps} deg ({len(angs_near)} angles): {mn}")
    if mf and mf > 0:
        lam = max(0.0, 1.0 - (mn if mn is not None else 1.0) / mf)
        print(f"  scale y by 1/{mf:.6f}; lam = {lam:.6f}; HEURISTIC upper bound |y|/minfar + cap*lam "
              f"= {tot/mf + cap*lam:.6f}  (NOT a certificate: angles are sampled)")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    h = sub.add_parser('hist'); h.add_argument('files', nargs='+')
    h.add_argument('--eps', type=float, nargs='+', default=[0.5, 1.0, 2.0, 5.0, 10.0])
    h.add_argument('--cap', type=float, default=3.0)

    c = sub.add_parser('cg')
    c.add_argument('T'); c.add_argument('TAG')
    c.add_argument('--eps', type=float, required=True, help='near-axis threshold, degrees')
    c.add_argument('--cap', type=float, default=None, help='cap on near-axis mass (default T-1)')
    c.add_argument('--rounds', type=int, default=60); c.add_argument('--time', type=float, default=3000)
    c.add_argument('--polish-time', type=float, default=600)
    c.add_argument('--polish-after', type=int, default=10 ** 9)
    c.add_argument('--seed-pitch', type=float, default=0.05); c.add_argument('--seed-dth', type=float, default=5.0)
    c.add_argument('--row-pitch', type=float, default=0.02)
    c.add_argument('--seed-file', nargs='*', default=None)
    c.add_argument('--threads', type=int, default=6)
    c.add_argument('--cg-want', type=int, default=1500)
    c.add_argument('--row-cap', type=int, default=30000)
    c.add_argument('--raster-dth', type=float, default=0.25)
    c.add_argument('--sweep-rows', type=int, default=350); c.add_argument('--sweep-k', type=int, default=6)
    c.add_argument('--col-age', type=int, default=4); c.add_argument('--row-age', type=int, default=6)

    e = sub.add_parser('exact')
    e.add_argument('--src', required=True); e.add_argument('--t', required=True)
    e.add_argument('--eps', type=float, nargs='+', required=True, help='0 or negative = uncapped control')
    e.add_argument('--cap', type=float, nargs='+', default=None, help='list of caps; default T-1')
    e.add_argument('--tag', required=True)
    e.add_argument('--Q', type=int, default=10 ** 7); e.add_argument('--Dc', type=int, default=10 ** 8)
    e.add_argument('--DM', type=int, default=10 ** 9); e.add_argument('--procs', type=int, default=6)

    r = sub.add_parser('rotseed')
    r.add_argument('--src', required=True); r.add_argument('--t', required=True)
    r.add_argument('--alpha', type=float, nargs='+', required=True, help='rotation angles, degrees')
    r.add_argument('--out', required=True)

    tb = sub.add_parser('table'); tb.add_argument('tags', nargs='+')

    bd = sub.add_parser('bands'); bd.add_argument('support'); bd.add_argument('--t', default=None)

    d = sub.add_parser('dualcost')
    d.add_argument('support'); d.add_argument('--eps', type=float, required=True)
    d.add_argument('--cap', type=float, default=None); d.add_argument('--dth', type=float, default=0.05)
    d.add_argument('--threads', type=int, default=6)

    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    dict(hist=cmd_hist, cg=cmd_cg, exact=cmd_exact, dualcost=cmd_dualcost, rotseed=cmd_rotseed, bands=cmd_bands, table=cmd_table)[a.cmd](a)


if __name__ == '__main__':
    main()
