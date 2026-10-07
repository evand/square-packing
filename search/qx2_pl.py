#!/usr/bin/env python3
"""qx2_pl.py -- continuous piecewise-linear ("hat") line densities for the qx2 fixed-profile LP (2026-10-01).
HEURISTIC (float LP, float oracle).  Write-up: search/QX2_PL.md.  Model: search/qx2_lp.py / search/quadrant_lp.py
(imported, not edited; Q.chord is the strict-containment version installed by qx2_lp).

Element: a RAMP = a segment [p,q] on an axis-parallel line with linear density, normalised to the element's mass x;
CDF F(t) = ((t-p)/(q-p))^2 (up, dir=+1) or 1 - ((q-t)/(q-p))^2 (down, dir=-1).  Captured fraction for the chord
[lo,hi] of the line inside the square: F(min(hi,q)) - F(max(lo,p)).
A HAT variable at node t_i = up-ramp on [t_{i-1}, t_i] + down-ramp on [t_i, t_{i+1}], each piece of mass x (peak
density 2x/h, so equal-pitch hats give a continuous density); at the ends of a line's allowed range a half-hat.
  * corner module: lines y = c (c on the hc-grid, 0 < c <= R), nodes on the hc-grid of the same range qx2_lp's
    uniform pieces cover ([0,R], or [0, ~a] for c > a), plus the diagonal mirror x = c;
  * profile horizontal lines (y on the hp-grid, 0 < y <= a): periodic hats in phase (no ends), orbit {phi, 1-phi}
    (the mirror maps up-ramps to down-ramps and the hat at phi to the hat at 1-phi); copies j = R..R+J on both walls,
    every ramp living in its period cell [j, j+1] (so the band is pi restricted to x >= R, as in the uniform model);
  * profile vertical lines (phase on the hp-grid): hats in y on [0, a] (half-hats at y = 0 and y = a).
Bookkeeping: newvar(obj, sig): obj = number of pieces counted in M_C (corner) or m_v (phase-0 vertical), sig = number
of pieces per period, exactly as the uniform code counts its pieces.  Lebesgue layer and tilt margin: XModel.

Comparison driver (`run`): uniform (qx2_lp elements) or hat elements, FLOAT rows only: no germ rows, and the exact
axis rows replaced by float theta = 0 rows at the one-sided corners (+-1e-7) of the breakpoint grid plus cell and edge
midpoints.  Oracle: quadrant_lp.oracle with qx2_lp's settings.  Stop: oracle min >= 1 - 1e-6 in two consecutive
rounds, or D flat to 1e-4 over 5 rounds, or 40 rounds.  `kmax` mode: maximise the tilt-margin slope kappa
(qx2_kmax.KLP) on the same float rows.
Outputs runs/qpl_<tag>/ : log.txt, hist.json, sol.npz, support.txt, binding.txt.
"""
import sys, os, math, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import numpy as np
import scipy.sparse as sp
import multiprocessing as mp
import qx2_lp as X          # strict containment (Q.TOL = 0, Q.chord strict)
Q = X.Q

RAMPS = ('rh', 'rv')        # rh: line y = c, param x; rv: line x = c, param y.  Row: (c, p, q, dir, var, fm)


def ramp_frac(lo, hi, p, q, d):
    """captured fraction of a ramp [p,q] (dir d = +1 up / -1 down) by the chord [lo,hi]."""
    L = q - p
    s = np.clip((np.maximum(lo, p) - p) / L, 0.0, 1.0)
    t = np.clip((np.minimum(hi, q) - p) / L, 0.0, 1.0)
    dt = np.maximum(t - s, 0.0)
    return np.where(d > 0, dt * (t + s), dt * (2.0 - t - s))


class HatModel(X.XModel):
    def __init__(self, R, w, a, kappa=0.0, th0=0.0):
        super().__init__(R, w, a, kappa, th0)
        for k in RAMPS: self.E[k] = []; self.K[k] = []

    def ramp(self, k, c, p, q, d, var, fm=0.0):
        self._add(k, (c, p, q, d, var, fm))

    def finalize(self):
        RE = {k: self.E.pop(k) for k in RAMPS}; RK = {k: self.K.pop(k) for k in RAMPS}
        super().finalize()
        for k in RAMPS:
            self.E[k] = RE[k]; self.K[k] = RK[k]
            self.A[k] = np.array(RE[k], float).reshape(-1, 6); self.KA[k] = np.array(RK[k], int)

    def contrib(self, P, xval=None, keep=None):
        out = super().contrib(P, xval=xval)
        if xval is not None: tot = out
        else: A0, const = out; rows, cols, vals = [], [], []
        n = len(P); CH = 128
        for orient, k in ((1, 'rh'), (0, 'rv')):
            arr = self.A.get(k)
            if arr is None or not len(arr): continue
            if xval is not None:
                var = arr[:, 4].astype(int); mm = np.where(var >= 0, xval[np.maximum(var, 0)], arr[:, 5])
                arr = arr[mm > 0]
                if not len(arr): continue
            for i in range(0, n, CH):
                B = P[i:i + CH]; ib = np.arange(i, min(i + CH, n))
                cc = B[:, 1:2] if orient == 1 else B[:, 0:1]
                pc = B[:, 0:1] if orient == 1 else B[:, 1:2]
                near = ((np.abs(arr[None, :, 0] - cc) <= 0.7072) & (arr[None, :, 1] <= pc + 0.7072) &
                        (arr[None, :, 2] >= pc - 0.7072))
                a, b = np.nonzero(near)
                if not len(a): continue
                lo, hi = Q.chord(orient, arr[b, 0], B[a, 0], B[a, 1], B[a, 2])
                fr = ramp_frac(lo, hi, arr[b, 1], arr[b, 2], arr[b, 3])
                ok = fr > 0
                pi = ib[a[ok]]; ei = b[ok]; fr = fr[ok]
                var = arr[ei, 4].astype(int); fm = arr[ei, 5]
                if xval is not None:
                    np.add.at(tot, pi, fr * np.where(var >= 0, xval[np.maximum(var, 0)], fm))
                else:
                    fx = var < 0
                    if fx.any(): np.add.at(const, pi[fx], fr[fx] * fm[fx])
                    rows.append(pi[~fx]); cols.append(var[~fx]); vals.append(fr[~fx])
        if xval is not None: return tot
        if rows:
            A1 = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                               shape=(n, self.nvar))
            A0 = (A0 + A1).tocsr()
        return A0, const


# ---------------------------------------------------------------------------------------------- builders
def hat_pieces(nodes, i):
    out = []
    if i > 0: out.append((nodes[i - 1], nodes[i], +1))
    if i < len(nodes) - 1: out.append((nodes[i], nodes[i + 1], -1))
    return out


def add_corner_hats(m, hc):
    R = m.R; a = m.a; g = Q.lattice(0, R, hc)
    m.kind = 0
    for y in g:
        if y <= 0: continue
        # the uniform model's pieces on this line: [g_i, g_i+1] unless y > a and g_i >= a
        last = max(i for i in range(len(g) - 1) if not (y > a + 1e-12 and g[i] >= a - 1e-12))
        nodes = g[:last + 2]
        for i in range(len(nodes)):
            pcs = hat_pieces(nodes, i)
            var = m.newvar(2 * len(pcs), 0, ('Ch', y, nodes[i]))
            for (p, q, d) in pcs:
                m.ramp('rh', y, p, q, d, var); m.ramp('rv', y, p, q, d, var)
    m.const_obj += (R - a) ** 2


def _orbit(phi):
    return sorted(set([round(phi % 1.0, 10), round((-phi) % 1.0, 10)]))


def add_profile_hats(m, hp):
    w = m.w; a = m.a
    phs = Q.lattice(0, 1, hp); nodes_p = phs[:-1]; ys = Q.lattice(0, w, hp)
    m.kind = 1
    # horizontal lines: periodic hats in phase
    for y in ys:
        if y <= 0 or y > a + 1e-12: continue
        for phi in nodes_p:
            orb = _orbit(phi)
            if orb[0] != round(phi, 10): continue
            pcs = []
            for p_ in orb:
                k = min(range(len(nodes_p)), key=lambda i: abs(nodes_p[i] - p_))
                pcs.append((phs[k], phs[k + 1], -1))
                pcs.append((phs[k - 1], phs[k], +1) if k > 0 else (phs[-2], 1.0, +1))
            var = m.newvar(0, len(pcs), ('Hh', y, phi))
            for (p, q, d) in pcs:
                for j in m.js():
                    m.ramp('rh', y, j + p, j + q, d, var)
                    if not m.band: m.ramp('rv', y, j + p, j + q, d, var)
    # vertical lines: hats in y on [0, a]
    for phi in nodes_p:
        orb = _orbit(phi)
        if orb[0] != round(phi, 10): continue
        last = max(i for i in range(len(ys) - 1) if ys[i + 1] <= a + 1e-12)
        nodes = ys[:last + 2]
        mv = 1 if orb == [0.0] else 0
        for i in range(len(nodes)):
            pcs = hat_pieces(nodes, i)
            var = m.newvar(mv * len(pcs), len(orb) * len(pcs), ('Hv', phi, nodes[i]))
            for p_ in orb:
                for (y0, y1, d) in pcs:
                    for j in m.js():
                        m.ramp('rv', j + p_, y0, y1, d, var)
                        if not m.band: m.ramp('rh', j + p_, y0, y1, d, var)
    m.const_sig += (w - a)


def build(a):
    """a: namespace with R, w, delta, kappa, th0 (deg), hc, hp, elem in {'uniform','hat'}."""
    if a.elem == 'uniform':
        b = argparse.Namespace(mode='lp', R=a.R, w=a.w, delta=a.delta, kappa=a.kappa, th0=a.th0, hc=a.hc, hca=0.0,
                               hp=a.hp, hpa=0.0, no_points=True, seam_hp=0.0)
        return X.build(b)
    m = HatModel(a.R, a.w, a.w - a.delta, a.kappa, math.radians(a.th0))
    if a.elem == 'both':          # uniform pieces AND hats (superset of both bases); Lebesgue constants counted once
        X.add_corner(m, a.hc, 0.0, points=False); X.add_profile(m, a.hp, 0.0, points=False)
        m.const_obj -= (m.R - m.a) ** 2; m.const_sig -= (m.w - m.a)
    add_corner_hats(m, a.hc); add_profile_hats(m, a.hp)
    m.finalize(); m.hlat = min(a.hc, a.hp)
    return m


# ---------------------------------------------------------------------------------------------- float axis rows
def breakgrid(m, lo, hi):
    xs = set([m.a, 0.0])
    for k, arr in m.A.items():
        if not len(arr): continue
        if k in ('hs', 'vs') or k in RAMPS: cols = (0, 1, 2)
        elif k == 'pts': cols = (0, 1)
        else: cols = (0, 1, 2, 3)
        for c in cols: xs.update(np.round(arr[:, c], 9).tolist())
    out = set()
    for x in xs:
        for s in (-0.5, 0.5):
            v = round(x + s, 9)
            if lo - 1e-9 <= v <= hi + 1e-9: out.add(v)
    return np.array(sorted(out))


def axis_poses(m, eps=1e-7, mids=True):
    """theta = 0 poses: one-sided corners (+-eps) of every cell of the breakpoint grid, plus (mids) the cell midpoints
    and edge midpoints (all combinations of {g - eps, g + eps, midpoints} in x and y)."""
    ex, ey = m.R + 1.75, m.R + 0.75
    gx = breakgrid(m, 0.5, ex); gy = breakgrid(m, 0.5, ey)
    cx = [gx - eps, gx + eps]; cy = [gy - eps, gy + eps]
    if mids: cx.append(0.5 * (gx[1:] + gx[:-1])); cy.append(0.5 * (gy[1:] + gy[:-1]))
    Xs = np.concatenate(cx); Ys = np.concatenate(cy)
    Xg, Yg = np.meshgrid(Xs, Ys, indexing='ij'); Xg = Xg.ravel(); Yg = Yg.ravel()
    k = (Yg >= 0.5) & (Xg >= 0.5) & (Yg <= Xg + 1e-6)
    return np.c_[Xg[k], Yg[k], np.zeros(k.sum())]


# ---------------------------------------------------------------------------------------------- driver
def binding_summary(lp, m, x, od, log):
    dual = np.array(lp.h.getSolution().row_dual)
    nfix = len(dual) - lp.nrows            # sigma row(s) first
    d = dual[nfix:]; P = lp.P
    order = np.argsort(-np.abs(d)); act = order[np.abs(d[order]) > 1e-9]
    th = np.degrees(np.minimum(P[act, 2], math.pi / 2 - P[act, 2]))
    with open(os.path.join(od, 'binding.txt'), 'w') as f:
        f.write("# dual cx cy theta_deg  (pose rows, by |dual|)\n")
        for i in act: f.write(f"{d[i]:+.6e} {P[i,0]:.7f} {P[i,1]:.7f} {math.degrees(P[i,2]):.7f}\n")
    W = np.abs(d[act]); Wt = W.sum() if len(W) else 1.0
    bins = [(0, 0), (0, 1e-3), (1e-3, 0.6), (0.6, 5), (5, 46)]
    parts = []
    for lo, hi in bins:
        k = (th == 0) if hi == 0 else ((th > lo) & (th <= hi))
        parts.append(f"th'{'=0' if hi == 0 else f' in ({lo},{hi}]'}: {k.sum()} rows, dual share {W[k].sum()/Wt:.3f}")
    # are theta = 0 binding poses at grid corners (within 1e-6 of a breakpoint in both coords) or interior?
    gx = breakgrid(m, 0.0, m.R + 3); z = act[th == 0]
    near = lambda v: np.min(np.abs(v[:, None] - gx[None, :]), axis=1) < 1e-6
    if len(z):
        c2 = near(P[z, 0]) & near(P[z, 1]); c1 = near(P[z, 0]) ^ near(P[z, 1])
        parts.append(f"theta=0 binding: {c2.sum()} at grid corners, {c1.sum()} on grid edges, "
                     f"{(~c2 & ~c1).sum()} interior")
    log(f"binding rows: {len(act)};  " + ';  '.join(parts))


def run(a):
    od = os.path.join(Q.RUNS, 'qpl_' + a.tag); os.makedirs(od, exist_ok=True)
    logf = open(os.path.join(od, 'log.txt'), 'a')

    def log(*s):
        msg = ' '.join(str(t) for t in s); print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    T0 = time.time()
    m = build(a)
    log(f"== qpl_{a.tag} elem={a.elem} R={m.R} w={m.w} a={m.a} kappa={m.kappa} th0={a.th0}deg nvar={m.nvar} "
        f"elements={ {k: len(v) for k, v in m.A.items()} } args={vars(a)}")
    rng = np.random.default_rng(a.seed)
    lp = Q.LP(m); lp.solver = 'simplex'; lp.crossover = 'on'
    lp.h.setOptionValue('primal_feasibility_tolerance', 1e-10); lp.h.setOptionValue('dual_feasibility_tolerance', 1e-10)
    P0 = np.r_[Q.lattice_poses(m.G, a.pitch, a.degs), Q.near_lattice_poses(m.G, m.hlat, full=False)]
    n0 = lp.add(P0); PA = axis_poses(m); na = lp.add(PA)
    log(f"initial rows: {n0} lattice + {na} float axis rows (of {len(PA)} poses) [{time.time()-T0:.1f}s]")
    hist = []; pool = None; x = None; status = 'rounds'
    for rnd in range(a.rounds + 1):
        t1 = time.time()
        x, obj, st = lp.solve()
        if x is None:
            log(f"round {rnd}: LP status {st}"); status = 'LP ' + st; break
        D = m.R ** 2 - obj - m.const_obj
        tlp = time.time() - t1
        if pool is not None: pool.terminate()
        if a.nproc > 1: pool = mp.Pool(a.nproc, initializer=Q._init, initargs=(m, x))
        else: Q._init(m, x); pool = Q.Serial()
        t2 = time.time()
        mn, fp, fv, wp = Q.oracle(pool, m, x, rng, nrand=a.nrand, lat_pitch=a.oracle_pitch, npol=a.npol, nproc=a.nproc)
        rows = Q.pick_rows(fp, fv, maxn=a.maxadd)
        hist.append(dict(round=rnd, D=D, rows=lp.nrows, oracle_min=mn, nviol=int(len(fp)), worst=[float(t) for t in wp],
                         tlp=tlp, tor=time.time() - t2))
        log(f"round {rnd}: D={D:.7f} rows={lp.nrows} | oracle min(-margin)={mn:.9f} at ({wp[0]:.6f},{wp[1]:.6f},"
            f"{math.degrees(wp[2]):.5f}deg) viol={len(fp)} add={len(rows)} [lp {tlp:.1f}s, oracle {time.time()-t2:.1f}s]")
        json.dump(hist, open(os.path.join(od, 'hist.json'), 'w'), indent=1)
        np.savez(os.path.join(od, 'sol.npz'), x=x, D=D, round=rnd, rows=lp.P, args=json.dumps(vars(a)))
        if len(hist) >= 2 and all(h['oracle_min'] >= 1 - a.stop_tol for h in hist[-2:]):
            status = 'clean2'; break
        if len(hist) >= 6 and abs(hist[-1]['D'] - hist[-6]['D']) < a.flat_tol:
            status = 'flat5'; break
        if rnd == a.rounds or len(rows) == 0:
            status = 'rounds' if len(rows) else 'norows'; break
        lp.add(rows)
    if pool is not None: pool.terminate()
    if x is not None and hist:
        sup = [(m.label[i], float(x[i])) for i in np.nonzero(x > 1e-12)[0]]
        with open(os.path.join(od, 'support.txt'), 'w') as f:
            for lab, val in sorted(sup, key=lambda t: str(t[0])): f.write(f"{lab}\t{val:.12f}\n")
        try: binding_summary(lp, m, x, od, log)
        except Exception as e: log(f"binding summary failed: {e}")
    h = hist[-1] if hist else {}
    log(f"FINAL elem={a.elem} kappa={a.kappa} th0={a.th0} status={status} D={h.get('D')} "
        f"oracle_min={h.get('oracle_min')} rounds={len(hist)} time={time.time()-T0:.0f}s")
    json.dump(dict(status=status, D=h.get('D'), oracle_min=h.get('oracle_min'), rounds=len(hist),
                   time=time.time() - T0, args=vars(a)), open(os.path.join(od, 'final.json'), 'w'))
    return hist


def run_kmax(a):
    """maximise kappa (tilt-margin slope at th0) on float rows; qx2_kmax.KLP, oracle at the current kappa."""
    import qx2_kmax as K
    od = os.path.join(Q.RUNS, 'qplk_' + a.tag); os.makedirs(od, exist_ok=True)
    logf = open(os.path.join(od, 'log.txt'), 'a')

    def log(*s):
        msg = ' '.join(str(t) for t in s); print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    T0 = time.time()
    a.kappa = 1.0
    m = build(a)
    log(f"== qplk_{a.tag} elem={a.elem} R={m.R} w={m.w} th0={a.th0}deg nvar={m.nvar} args={vars(a)}")
    rng = np.random.default_rng(a.seed)
    lp = K.KLP(m, a.kcap, a.dmin)
    P0 = np.r_[Q.lattice_poses(m.G, a.pitch, a.degs), Q.near_lattice_poses(m.G, m.hlat, full=False), axis_poses(m)]
    log(f"initial rows: {lp.add(P0)}")
    pool = None; hist = []; status = 'rounds'; kap = None; D = None
    for rnd in range(a.rounds + 1):
        x, kap, st = lp.solve()
        if x is None: log(f"round {rnd}: LP status {st}"); status = 'LP ' + st; break
        D = m.R ** 2 - float(m.cobj @ x) - m.const_obj
        m.kappa = max(kap, 1e-12)
        if pool is not None: pool.terminate()
        if a.nproc > 1: pool = mp.Pool(a.nproc, initializer=Q._init, initargs=(m, x))
        else: Q._init(m, x); pool = Q.Serial()
        t2 = time.time()
        mn, fp, fv, wp = Q.oracle(pool, m, x, rng, nrand=a.nrand, lat_pitch=a.oracle_pitch, npol=a.npol, nproc=a.nproc)
        rows = Q.pick_rows(fp, fv, maxn=a.maxadd)
        hist.append(dict(round=rnd, kappa=kap, D=D, oracle_min=mn))
        log(f"round {rnd}: kappa={kap:.6f} D={D:.6f} rows={lp.nrows} | oracle min={mn:.9f} at ({wp[0]:.5f},{wp[1]:.5f},"
            f"{math.degrees(wp[2]):.5f}deg) viol={len(fp)} add={len(rows)} [oracle {time.time()-t2:.1f}s]")
        np.savez(os.path.join(od, 'sol.npz'), x=x, kappa=kap, D=D, round=rnd, rows=lp.P)
        if len(hist) >= 2 and all(h['oracle_min'] >= 1 - a.stop_tol for h in hist[-2:]): status = 'clean2'; break
        if len(hist) >= 6 and abs(hist[-1]['kappa'] - hist[-6]['kappa']) < a.flat_tol: status = 'flat5'; break
        if rnd == a.rounds or len(rows) == 0: break
        lp.add(rows)
    if pool is not None: pool.terminate()
    log(f"FINAL kmax elem={a.elem} th0={a.th0} dmin={a.dmin} status={status} kappa*={kap} D={D} "
        f"oracle_min={hist[-1]['oracle_min'] if hist else None} rounds={len(hist)} time={time.time()-T0:.0f}s")


def parser():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['lp', 'kmax'])
    ap.add_argument('--elem', choices=['uniform', 'hat', 'both'], default='hat')
    ap.add_argument('--R', type=int, default=3); ap.add_argument('--w', type=float, default=3.0)
    ap.add_argument('--delta', type=float, default=0.2)
    ap.add_argument('--hc', type=float, default=0.2); ap.add_argument('--hp', type=float, default=0.2)
    ap.add_argument('--kappa', type=float, default=0.0); ap.add_argument('--th0', type=float, default=3.0)
    ap.add_argument('--dmin', type=float, default=None); ap.add_argument('--kcap', type=float, default=1.0)
    ap.add_argument('--pitch', type=float, default=0.1)
    ap.add_argument('--degs', type=lambda s: [float(t) for t in s.split(',')],
                    default=[0, 1e-4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85])
    ap.add_argument('--rounds', type=int, default=40)
    ap.add_argument('--nrand', type=int, default=1000000); ap.add_argument('--oracle-pitch', type=float, default=0.02)
    ap.add_argument('--npol', type=int, default=2000); ap.add_argument('--maxadd', type=int, default=4000)
    ap.add_argument('--stop-tol', type=float, default=1e-6); ap.add_argument('--flat-tol', type=float, default=1e-4)
    ap.add_argument('--nproc', type=int, default=2); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tag', required=True)
    return ap


if __name__ == '__main__':
    args = parser().parse_args()
    if args.mode == 'kmax': run_kmax(args)
    else: run(args)
