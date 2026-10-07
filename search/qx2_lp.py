#!/usr/bin/env python3
"""qx2_lp.py -- exact-friendly re-solve of the (R = 2, w = 2) quadrant LP (task quadrant-exact-w2, 2026-09-28).
HEURISTIC (float LP, float oracle).  Write-up: search/QUADRANT_EXACT.md.  Model and oracle: search/quadrant_lp.py
(imported, not edited).

Differences from quadrant_lp.py `lp` mode:
  * Lebesgue layer: the measure is exactly Lebesgue on U = [a, inf)^2 (a = w - delta), i.e. the profile is Lebesgue on
    y in [a, w] and the corner module is Lebesgue on [a, R]^2; no atom lies in the interior of U.  A square inside U
    then has mass exactly 1 by a trivial lemma (QUADRANT_EXACT.md, Lemma U).
  * optional: no point atoms (--no-points), finer vertical pieces on the seam (--seam-hp).
  * exact theta = 0 face: rows at every one-sided limit corner of the axis-parallel arrangement (at theta = 0 the mass is
    multilinear on each open cell of the breakpoint grid, so these finitely many rows imply validity at theta = 0).
  * tilt margin: every row at angle th asks for 1 + kappa * min(th', th0) * (1 - area(Q cap U)), th' = distance of th
    to 0 (mod 90 deg); squares inside U have mass exactly 1, so the margin must vanish there.
  * vertex solutions (HiGHS simplex) by default.
Outputs runs/qx2_<tag>/ : log.txt, hist.json, sol.npz (x, D, rows, args).
"""
import sys, os, math, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import numpy as np
import multiprocessing as mp
import quadrant_lp as Q

BIG = 1e6
RUNS = Q.RUNS

# STRICT containment.  quadrant_lp's TOL = 1e-9 counts chord ends up to 1e-9 outside the square; divided by sin(theta)
# that over-counts the mass by ~1e-9/theta at tiny tilts (2e-4 at theta = 1e-5), which hid a genuine germ violation of
# 6.4e-5 at the corner germ (1.5, 1.5) from the LP's oracle (found by qx2_zm.py, confirmed exactly, QUADRANT_EXACT.md
# sec 6).  With tol = 0 the float error is ~1e-16/theta.
Q.TOL = 0.0
_chord_lenient = Q.chord
def _chord_strict(orient, a, cx, cy, th, tol=0.0):
    return _chord_lenient(orient, a, cx, cy, th, tol)
Q.chord = _chord_strict


class XModel(Q.Model):
    def __init__(self, R, w, a, kappa=0.0, th0=0.0):
        super().__init__(R, w)
        self.a = float(a); self.kappa = float(kappa); self.th0 = float(th0); self.use_margin = True

    def margin(self, th):
        t = np.mod(th, math.pi / 2); t = np.minimum(t, math.pi / 2 - t)
        return self.kappa * np.minimum(t, self.th0)

    def lebesgue(self, P):
        V = Q.corners(P); n = len(P); a = self.a
        if self.band:
            ar = Q.area_rect(V, np.full(n, -BIG), np.full(n, BIG), np.full(n, a), np.full(n, BIG))
        else:
            ar = Q.area_rect(V, np.full(n, a), np.full(n, BIG), np.full(n, a), np.full(n, BIG))
        if self.use_margin and self.kappa > 0: ar = ar - self.margin(P[:, 2]) * np.clip(1.0 - ar, 0.0, 1.0)
        return ar


def add_corner(m, hc, hca, points=True):
    R = m.R; a = m.a; g = Q.lattice(0, R, hc)
    inU = lambda x, y: x > a + 1e-12 and y > a + 1e-12
    if points:
        for x in g:
            for y in g:
                if x <= 0 or y <= 0 or y > x or inU(x, y): continue
                m.corner_point(x, y, var=-2)
    for y in g:
        if y <= 0: continue
        for i in range(len(g) - 1):
            if y > a + 1e-12 and g[i] >= a - 1e-12: continue
            m.corner_hseg(y, g[i], g[i + 1], var=-2)
    if hca:
        ga = Q.lattice(0, R, hca)
        for i in range(len(ga) - 1):
            for j in range(len(ga) - 1):
                if ga[j] > ga[i]: continue
                if ga[i + 1] > a + 1e-12 and ga[j + 1] > a + 1e-12: continue
                m.corner_cell(ga[i], ga[i + 1], ga[j], ga[j + 1], var=-2)
    m.const_obj += (R - a) ** 2           # the Lebesgue piece [a,R]^2 of the corner module


def add_profile(m, hp, hpa, points=True, seam_hp=0.0):
    w = m.w; a = m.a
    phs = Q.lattice(0, 1, hp)[:-1]; ys = Q.lattice(0, w, hp)
    if points:
        for phi in phs:
            if (-phi) % 1.0 < phi - 1e-12: continue
            for y in ys:
                if y <= 0 or y > a + 1e-12: continue
                m.prof_point(phi, y, var=-2)
    for y in ys:
        if y <= 0 or y > a + 1e-12: continue
        for phi in phs:
            p0, p1 = phi, round(phi + hp, 10)
            if 1 - p1 < p0 - 1e-12: continue
            m.prof_hseg(y, p0, p1, var=-2)
    for phi in phs:
        if (-phi) % 1.0 < phi - 1e-12: continue
        yy = Q.lattice(0, w, seam_hp) if (seam_hp and phi == 0) else ys
        for i in range(len(yy) - 1):
            if yy[i + 1] > a + 1e-12: continue
            m.prof_vseg(phi, yy[i], yy[i + 1], var=-2)
    if hpa:
        pa = Q.lattice(0, 1, hpa); ya = Q.lattice(0, w, hpa)
        for i in range(len(pa) - 1):
            if 1 - pa[i + 1] < pa[i] - 1e-12: continue
            for j in range(len(ya) - 1):
                if ya[j + 1] > a + 1e-12: continue
                m.prof_cell(pa[i], pa[i + 1], ya[j], ya[j + 1], var=-2)
    m.const_sig += (w - a)               # the Lebesgue layer [a,w] of the profile, per period


def build(a):
    m = XModel(a.R, a.w, a.w - a.delta, a.kappa, math.radians(a.th0))
    if a.mode == 'band':
        m.band = True; m.G = ('band', float(a.w)); m.R = 0
        add_profile(m, a.hp, a.hpa, points=not a.no_points, seam_hp=a.seam_hp)
        m.finalize(); m.cobj = -m.cobj
    else:
        add_corner(m, a.hc, a.hca, points=not a.no_points)
        add_profile(m, a.hp, a.hpa, points=not a.no_points, seam_hp=a.seam_hp)
        m.finalize()
    m.hlat = min(a.hc, a.hp) if a.mode != 'band' else a.hp
    return m


def breakpoints(m, lo, hi):
    """all x- and y-coordinates at which some element (or the Lebesgue region) has a breakpoint, +- 1/2, in [lo,hi]."""
    xs = set([m.a, 0.0])
    for k, arr in m.A.items():
        if not len(arr): continue
        if k == 'pts': cols = (0, 1)
        elif k in ('hs', 'vs'): cols = (0, 1, 2)
        else: cols = (0, 1, 2, 3)
        for c in cols: xs.update(np.round(arr[:, c], 9).tolist())
    out = set()
    for x in xs:
        for s in (-0.5, 0.5):
            v = round(x + s, 9)
            if lo - 1e-9 <= v <= hi + 1e-9: out.add(v)
    return np.array(sorted(out))


def axis_rows(m, eps=1e-7):
    """theta = 0 one-sided limit poses at every corner of every open cell of the breakpoint grid."""
    if m.band:
        ex, ey = 1.0, m.w + 0.75
        gx = np.r_[breakpoints(m, -1.0, 2.0)]; gy = breakpoints(m, 0.5, ey)
    else:
        ex, ey = m.R + 1.75, m.R + 0.75
        gx = breakpoints(m, 0.5, ex); gy = breakpoints(m, 0.5, ey)
    out = []
    for sx in (-1, 1):
        for sy in (-1, 1):
            X, Y = np.meshgrid(gx + sx * eps, gy + sy * eps, indexing='ij')
            X = X.ravel(); Y = Y.ravel()
            k = (Y >= 0.5) & (X >= 0.5 if not m.band else True)
            if not m.band: k &= (Y <= X + 1e-6)
            out.append(np.c_[X[k], Y[k], np.zeros(k.sum())])
    P = np.concatenate(out)
    if m.band: P[:, 0] = np.mod(P[:, 0], 1.0)
    return P


def add_axis_rows(lp, m):
    """the theta = 0 one-sided-limit rows, with EXACT coefficients (qx2_exact.ExactModel.coeffs), so that the LP
    satisfies the theta = 0 face to LP tolerance (Lemma Z)."""
    import qx2_exact as E
    from fractions import Fraction as F
    EM = E.ExactModel(m)
    R = m.R; ex = R + F(7, 4); ey = R + F(3, 4); H = F(1, 2)
    g = E.breakgrid(EM, H, ex + 1)
    starts, idx, vals, lo = [], [], [], []
    nnz = 0
    for cx in g:
        for cy in g:
            if cy > cx or cy > ey + 1: continue
            for sx in (-1, 1):
                for sy in (-1, 1):
                    if cx == H and sx < 0: continue
                    if cy == H and sy < 0: continue
                    co, ar = EM.coeffs(cx, cy, 0, sx, sy)
                    rhs = 1 - ar
                    if rhs <= 0: continue
                    starts.append(nnz)
                    for j, c in sorted(co.items()):
                        idx.append(j); vals.append(float(c)); nnz += 1
                    lo.append(float(rhs))
    n = len(lo)
    lp.h.addRows(n, np.array(lo), np.full(n, lp.inf), nnz, np.array(starts, dtype=np.int32),
                 np.array(idx, dtype=np.int32), np.array(vals, dtype=float))
    lp.naxis = n
    return n


class GermRows:
    """the exact tile-germ limit rows (qx2_exact.germ_rows_exact) as a float sparse matrix; rows are handed to the LP
    only when violated by the current solution (QUADRANT_EXACT.md 1.5)."""

    def __init__(self, m, nproc):
        import qx2_exact as E
        import scipy.sparse as sp
        EM = E.ExactModel(m)
        rows = E.germ_rows_exact(EM, m.R, nproc=nproc)
        r, c, v, rhs = [], [], [], []
        seen = set()
        k = 0
        for pose, co, ar in rows:
            key = (tuple(sorted((j, c_) for j, c_ in co.items())), ar)
            if key in seen or 1 - ar <= 0: continue
            seen.add(key)
            for j, c_ in co.items(): r.append(k); c.append(j); v.append(float(c_))
            rhs.append(float(1 - ar)); k += 1
        self.A = sp.csr_matrix((v, (r, c)), shape=(k, m.nvar)); self.rhs = np.array(rhs)
        self.added = np.zeros(k, bool)

    def violated(self, x, tol=1e-10):
        val = self.A @ x - self.rhs
        return np.nonzero((val < -tol) & ~self.added)[0], float(val.min()) if len(val) else 0.0

    def add(self, lp, idx):
        sub = self.A[idx].tocsr(); sub.sort_indices()
        n = sub.shape[0]
        lp.h.addRows(n, self.rhs[idx], np.full(n, lp.inf), int(sub.nnz), sub.indptr[:-1].astype(np.int32),
                     sub.indices.astype(np.int32), sub.data.astype(float))
        self.added[idx] = True


def run(a):
    tag = a.tag
    od = os.path.join(RUNS, 'qx2_' + tag); os.makedirs(od, exist_ok=True)
    logf = open(os.path.join(od, 'log.txt'), 'a')

    def log(*s):
        msg = ' '.join(str(t) for t in s); print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    m = build(a)
    log(f"== qx2_{tag} mode={a.mode} R={m.R} w={m.w} a={m.a} kappa={m.kappa} th0={a.th0}deg nvar={m.nvar} "
        f"elements={ {k: len(v) for k, v in m.A.items()} } args={vars(a)}")
    rng = np.random.default_rng(a.seed)
    lp = Q.LP(m); lp.solver = a.solver; lp.crossover = a.crossover
    lp.h.setOptionValue('primal_feasibility_tolerance', 1e-10); lp.h.setOptionValue('dual_feasibility_tolerance', 1e-10)
    t0 = time.time()
    P0 = np.r_[Q.lattice_poses(m.G, a.pitch, a.degs), Q.near_lattice_poses(m.G, m.hlat, full=False)]
    n0 = lp.add(P0)
    na = add_axis_rows(lp, m)
    log(f"initial rows: {n0} lattice + {na} axis-limit rows [{time.time()-t0:.1f}s]")
    if a.warm:
        for wf in a.warm.split(','):
            W = np.load(wf)['rows']; nw = lp.add(W); log(f"warm rows from {wf}: {nw} of {len(W)}")
    hist = []; pool = None; x = None
    GR = None
    if a.germs:
        t1 = time.time(); GR = GermRows(m, max(a.nproc, 1))
        log(f"germ-limit rows: {GR.A.shape[0]} distinct [{time.time()-t1:.0f}s]")
    for rnd in range(a.rounds + 1):
        t1 = time.time()
        x, obj, st = lp.solve()
        if x is None: log(f"round {rnd}: LP status {st}"); break
        if GR is not None:
            for it in range(50):
                idx, gmin = GR.violated(x)
                if len(idx) == 0: break
                GR.add(lp, idx); x, obj, st = lp.solve()
                log(f"   germ rows: added {len(idx)} violated (min {gmin:.3e}) -> LP {st}")
                if x is None: break
            if x is None: log(f"round {rnd}: LP status {st}"); break
        if a.mode == 'band':
            D = -obj; mv = D
        else:
            D = m.R ** 2 - obj - m.const_obj
            mv = float(sum(m.cobj[i] * x[i] for i in range(m.nvar) if m.label[i][0][0] == 'P'))
        tlp = time.time() - t1
        if a.nproc > 1:
            if pool is not None: pool.terminate()
            pool = mp.Pool(a.nproc, initializer=Q._init, initargs=(m, x))
        else:
            Q._init(m, x); pool = Q.Serial()
        t2 = time.time()
        mn, fp, fv, wp = Q.oracle(pool, m, x, rng, nrand=a.nrand, lat_pitch=a.oracle_pitch, npol=a.npol, nproc=a.nproc)
        rows = Q.pick_rows(fp, fv, maxn=a.maxadd)
        hist.append(dict(round=rnd, D=D, mv=mv, rows=lp.nrows, oracle_min=mn, nviol=int(len(fp)),
                         worst=[float(t) for t in wp], tlp=tlp, tor=time.time() - t2))
        log(f"round {rnd}: D={D:.7f} m_v={mv:.6f} rows={lp.nrows} | oracle min(-margin)={mn:.9f} at "
            f"({wp[0]:.6f},{wp[1]:.6f},{math.degrees(wp[2]):.5f}deg) viol={len(fp)} add={len(rows)} "
            f"[lp {tlp:.1f}s, oracle {time.time()-t2:.1f}s]")
        json.dump(hist, open(os.path.join(od, 'hist.json'), 'w'), indent=1)
        np.savez(os.path.join(od, 'sol.npz'), x=x, D=D, round=rnd, rows=lp.P, args=json.dumps(vars(a)))
        if mn >= 1 - a.stop_tol:
            clean = True
            for extra in range(2):
                mn2, fp2, fv2, wp2 = Q.oracle(pool, m, x, rng, nrand=2 * a.nrand, lat_pitch=a.oracle_pitch * 0.77,
                                              npol=a.npol, nproc=a.nproc)
                log(f"   confirm pass {extra}: oracle min={mn2:.9f} viol={len(fp2)}")
                if mn2 < 1 - a.stop_tol:
                    rows = Q.pick_rows(fp2, fv2, maxn=a.maxadd); clean = False; break
            if clean: log(f"stable: oracle min >= 1 - {a.stop_tol} in 3 passes"); break
        if rnd == a.rounds or len(rows) == 0: break
        lp.add(rows)
    if pool is not None: pool.terminate()
    if x is not None:
        sup = [(m.label[i], float(x[i])) for i in np.nonzero(x > 1e-12)[0]]
        with open(os.path.join(od, 'support.txt'), 'w') as f:
            for lab, val in sorted(sup, key=lambda t: str(t[0])): f.write(f"{lab}\t{val:.12f}\n")
        log(f"support size {len(sup)}")
    return hist


def parser():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['lp', 'band'])
    ap.add_argument('--R', type=int, default=2); ap.add_argument('--w', type=float, default=2.0)
    ap.add_argument('--delta', type=float, default=0.2)
    ap.add_argument('--hc', type=float, default=0.1); ap.add_argument('--hca', type=float, default=0.25)
    ap.add_argument('--hp', type=float, default=0.1); ap.add_argument('--hpa', type=float, default=0.25)
    ap.add_argument('--seam-hp', type=float, default=0.0)
    ap.add_argument('--no-points', action='store_true')
    ap.add_argument('--kappa', type=float, default=0.0); ap.add_argument('--th0', type=float, default=0.0)
    ap.add_argument('--pitch', type=float, default=0.1)
    ap.add_argument('--degs', type=lambda s: [float(t) for t in s.split(',')],
                    default=[0, 1e-4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85])
    ap.add_argument('--rounds', type=int, default=40)
    ap.add_argument('--nrand', type=int, default=1000000); ap.add_argument('--oracle-pitch', type=float, default=0.02)
    ap.add_argument('--npol', type=int, default=2000); ap.add_argument('--maxadd', type=int, default=4000)
    ap.add_argument('--stop-tol', type=float, default=1e-9)
    ap.add_argument('--nproc', type=int, default=1); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tag', required=True); ap.add_argument('--warm', default='')
    ap.add_argument('--solver', default='simplex'); ap.add_argument('--crossover', default='on')
    ap.add_argument('--germs', action='store_true', help='exact tile-germ limit rows (added when violated)')
    return ap


if __name__ == '__main__':
    run(parser().parse_args())
