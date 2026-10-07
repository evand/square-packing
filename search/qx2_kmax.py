#!/usr/bin/env python3
"""qx2_kmax.py -- the tilt-margin ceiling of a qx2 fixed-profile family (2026-10-01).  HEURISTIC (float LP, float oracle).

qx2_lp.py asks every row at angle th for mass >= 1 + kappa * min(th', th0) * (1 - area(Q cap U)).  At R = w = 3 the
k2m3 setting kappa = 0.2 is infeasible (runs/qx2_k4_R3w3.out, round 4).  This script makes kappa a variable and
MAXIMISES it (optionally subject to D >= --dmin), with the same row generation (lattice + near-lattice + exact axis rows
+ exact germ rows + float oracle at the current kappa).  Output: kappa*, and the binding rows weighted by their duals
(which squares cap the slope).  kappa*(dmin) traced over dmin is the D(kappa) frontier; kappa*(w) at dmin = -inf is
the slope ceiling of width w.  qx2_lp.py / quadrant_lp.py are imported, not edited.

Rows:  A x - kappa g(P) >= 1 - area_U(P),   g(P) = min(th', th0) * clip(1 - area_U(P), 0, 1)   (th' = dist(th, 0 mod 90deg))
Axis rows (th = 0) and germ rows (th -> 0+) carry no margin, as in qx2_lp.
Outputs runs/qxk_<tag>/ : log.txt, sol.npz, binding.txt.
"""
import sys, os, math, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import numpy as np
import multiprocessing as mp
import qx2_lp as X          # sets Q.TOL = 0 (strict containment)
Q = X.Q


class KLP:
    """HiGHS LP over (x, kappa); minimise -kappa."""

    def __init__(self, m, kcap, dmin):
        import highspy
        self.hp = highspy; self.m = m
        h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('random_seed', 0)
        h.setOptionValue('primal_feasibility_tolerance', 1e-10); h.setOptionValue('dual_feasibility_tolerance', 1e-10)
        inf = highspy.kHighsInf; self.inf = inf
        n = m.nvar; self.n = n
        h.addVars(n + 1, np.zeros(n + 1), np.r_[np.full(n, inf), kcap])
        h.changeColsCost(n + 1, np.arange(n + 1, dtype=np.int32), np.r_[np.zeros(n), -1.0])
        self.info = []
        idx = np.nonzero(m.csig)[0].astype(np.int32)
        if len(idx):
            rhs = m.w - m.const_sig
            h.addRow(rhs, rhs, len(idx), idx, m.csig[idx].astype(float)); self.info.append(('sigma', None))
        if dmin is not None:      # D = R^2 - cobj.x - const_obj >= dmin
            ub = m.R ** 2 - m.const_obj - dmin
            idx = np.nonzero(m.cobj)[0].astype(np.int32)
            h.addRow(-inf, ub, len(idx), idx, m.cobj[idx].astype(float)); self.info.append(('dmin', None))
        self.h = h; self.nrows = len(self.info); self.P = np.zeros((0, 3))

    def rows_of(self, P):
        m = self.m
        m.use_margin = False
        A, ar = m.contrib(P)
        m.use_margin = True
        g = m.margin(P[:, 2]) / max(m.kappa, 1e-300) * np.clip(1.0 - ar, 0.0, 1.0)   # unit-kappa margin
        return A.tocsr(), ar, g

    def add(self, P):
        A, ar, g = self.rows_of(P)
        lo = 1.0 - ar
        keep = (lo > 1e-12) | (g > 0)
        A = A[keep]; lo = lo[keep]; g = g[keep]; P = P[keep]
        if A.shape[0] == 0: return 0
        import scipy.sparse as sp
        A = sp.hstack([A, sp.csr_matrix(-g.reshape(-1, 1))]).tocsr(); A.sort_indices(); A.eliminate_zeros()
        self.h.addRows(A.shape[0], lo, np.full(A.shape[0], self.inf), int(A.nnz), A.indptr[:-1].astype(np.int32),
                       A.indices.astype(np.int32), A.data.astype(float))
        self.nrows += A.shape[0]; self.P = np.r_[self.P, P]
        self.info += [('pose', tuple(p)) for p in P]
        return A.shape[0]

    def solve(self):
        h = self.h; h.run(); st = h.getModelStatus()
        if st != self.hp.HighsModelStatus.kOptimal:      # retry from scratch: IPM + crossover, then simplex again
            h.clearSolver(); h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', 'on')
            h.run(); st = h.getModelStatus(); h.setOptionValue('solver', 'simplex')
        if st != self.hp.HighsModelStatus.kOptimal: return None, None, str(st)
        sol = h.getSolution()
        v = np.array(sol.col_value); x = np.maximum(v[:self.n], 0.0)
        return x, float(v[self.n]), 'ok'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--R', type=int, default=3); ap.add_argument('--w', type=float, default=3.0)
    ap.add_argument('--delta', type=float, default=0.2)
    ap.add_argument('--hc', type=float, default=0.2); ap.add_argument('--hp', type=float, default=0.2)
    ap.add_argument('--th0', type=float, default=3.0)
    ap.add_argument('--dmin', type=float, default=None)
    ap.add_argument('--kcap', type=float, default=1.0)
    ap.add_argument('--pitch', type=float, default=0.1)
    ap.add_argument('--degs', type=lambda s: [float(t) for t in s.split(',')],
                    default=[0, 0.0001, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85])
    ap.add_argument('--rounds', type=int, default=40)
    ap.add_argument('--nrand', type=int, default=1000000); ap.add_argument('--oracle-pitch', type=float, default=0.02)
    ap.add_argument('--npol', type=int, default=2000); ap.add_argument('--maxadd', type=int, default=4000)
    ap.add_argument('--stop-tol', type=float, default=1e-9)
    ap.add_argument('--nproc', type=int, default=2); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tag', required=True); ap.add_argument('--warm', default='')
    ap.add_argument('--band', action='store_true',
                    help='band-only model (far band, no corner); --dmin then bounds m_v; float axis rows, no exact germ rows')
    a = ap.parse_args()
    # build the model exactly as qx2_lp does (lines only, Lebesgue layer, kappa placeholder 1)
    b = argparse.Namespace(mode='band' if a.band else 'lp', R=a.R, w=a.w, delta=a.delta, kappa=1.0, th0=a.th0, hc=a.hc, hca=0.0,
                           hp=a.hp, hpa=0.0, no_points=True, seam_hp=0.0)
    m = X.build(b)
    od = os.path.join(Q.RUNS, 'qxk_' + a.tag); os.makedirs(od, exist_ok=True)
    logf = open(os.path.join(od, 'log.txt'), 'a')

    def log(*s):
        msg = ' '.join(str(t) for t in s); print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    log(f"== qxk_{a.tag} R={m.R} w={m.w} a={m.a} th0={a.th0}deg dmin={a.dmin} nvar={m.nvar} args={vars(a)}")
    rng = np.random.default_rng(a.seed)
    lp = KLP(m, a.kcap, a.dmin)
    t0 = time.time()
    P0 = np.r_[Q.lattice_poses(m.G, a.pitch, a.degs), Q.near_lattice_poses(m.G, m.hlat, full=False)]
    n0 = lp.add(P0)
    if a.band:      # qx2_exact has no band model: float theta = 0 rows at the same one-sided grid corners
        na = lp.add(X.axis_rows(m))
    else:
        na = X.add_axis_rows(lp, m); lp.info += [('axis', None)] * na; lp.nrows += na
    log(f"initial rows: {n0} lattice + {na} axis-limit rows [{time.time()-t0:.1f}s]")
    if a.warm:
        for wf in a.warm.split(','):
            W = np.load(wf)['rows']; nw = lp.add(W); log(f"warm rows from {wf}: {nw} of {len(W)}")
    GR = None
    if not a.band:
        t1 = time.time(); GR = X.GermRows(m, max(a.nproc, 1))
        log(f"germ-limit rows: {GR.A.shape[0]} distinct [{time.time()-t1:.0f}s]")
    pool = None; x = None; kap = None
    for rnd in range(a.rounds + 1):
        t1 = time.time()
        x, kap, st = lp.solve()
        if x is None: log(f"round {rnd}: LP status {st}"); break
        for it in range(50 if GR is not None else 0):
            idx, gmin = GR.violated(x)
            if len(idx) == 0: break
            GR.add(lp, idx); lp.info += [('germ', int(i)) for i in idx]; lp.nrows += len(idx)
            x, kap, st = lp.solve()
            log(f"   germ rows: added {len(idx)} violated (min {gmin:.3e}) -> LP {st}")
            if x is None: break
        if x is None: log(f"round {rnd}: LP status {st}"); break
        D = m.R ** 2 - float(m.cobj @ x) - m.const_obj
        m.kappa = max(kap, 1e-12)
        if pool is not None: pool.terminate()
        pool = mp.Pool(a.nproc, initializer=Q._init, initargs=(m, x)) if a.nproc > 1 else Q.Serial()
        if a.nproc <= 1: Q._init(m, x)
        t2 = time.time()
        mn, fp, fv, wp = Q.oracle(pool, m, x, rng, nrand=a.nrand, lat_pitch=a.oracle_pitch, npol=a.npol, nproc=a.nproc)
        rows = Q.pick_rows(fp, fv, maxn=a.maxadd)
        log(f"round {rnd}: kappa={kap:.6f} D={D:.6f} rows={lp.nrows} | oracle min={mn:.9f} at "
            f"({wp[0]:.5f},{wp[1]:.5f},{math.degrees(wp[2]):.5f}deg) viol={len(fp)} add={len(rows)} "
            f"[lp {t2-t1:.1f}s, oracle {time.time()-t2:.1f}s]")
        np.savez(os.path.join(od, 'sol.npz'), x=x, kappa=kap, D=D, round=rnd, rows=lp.P, args=json.dumps(vars(a)))
        if mn >= 1 - a.stop_tol or len(rows) == 0:
            log("stable: oracle clean"); break
        if rnd == a.rounds: break
        lp.add(rows)
    if pool is not None: pool.terminate()
    if x is None: return
    # binding rows by dual weight
    dual = np.array(lp.h.getSolution().row_dual)
    order = np.argsort(-np.abs(dual))
    with open(os.path.join(od, 'binding.txt'), 'w') as f:
        f.write(f"# kappa*={kap:.8f} D={D:.6f}; rows by |dual| (pose rows: cx cy theta_deg; germ rows: index into GermRows)\n")
        for i in order:
            if abs(dual[i]) < 1e-9: break
            kind, val = lp.info[i] if i < len(lp.info) else ('?', None)
            if kind == 'pose':
                f.write(f"{dual[i]:+.6e} pose {val[0]:.6f} {val[1]:.6f} {math.degrees(val[2]):.6f}\n")
            else:
                f.write(f"{dual[i]:+.6e} {kind} {val}\n")
    log(f"final: kappa*={kap:.8f} D={D:.6f}; binding rows -> {od}/binding.txt")


if __name__ == '__main__':
    main()
