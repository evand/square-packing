#!/usr/bin/env python3
"""Bracket the fractional packing number nu_f(s)  (TODO item B1: the ceiling of the method).

Primal / dual pair (container C = [0,s]^2, Pi = placements of a CLOSED unit square inside C):

    COVER(s)   = min  sum_p w_p          s.t.  w(Q) >= 1  for every Q in Pi,   w >= 0
    nu_f(s)    = sup  mu(Pi)             s.t.  mu({Q : p in Q}) <= 1  for every point p in C,  mu >= 0

Weak duality:  for any feasible w, mu:   mu(Pi) <= int w(Q) dmu = sum_p w_p mu(Q ∋ p) <= sum_p w_p.
So   nu_f(s) <= COVER(s),   and a certificate (a feasible w) of weight < 12 exists at s  only if
nu_f(s) < 12.  Both functions are monotone non-decreasing in s.

This script produces, for a given s,

  L(s)  -- a RIGOROUS LOWER BOUND on nu_f(s): an explicit finitely supported measure y on exact
           placements, whose closed-square coverage is bounded above by M at EVERY point of C
           (M is computed by adaptive subdivision of C with a dilation that dominates the closed
           coverage on each cell), so y/M is feasible and  L = sum(y)/M <= nu_f(s).
           Every LP iterate gives such a bound; they are monotone in the LP's cut/column set only
           heuristically, so we keep the best.

  U(s') -- a RIGOROUS UPPER BOUND on COVER(s') (hence on nu_f(s')) for a range of s' >= s: the LP
           weights w are exported as a certificate, scaled to container s', and the exact i128
           verifier (verify/) returns the exact minimum m(s') of the captured weight over ALL
           closed unit squares in [0,s']^2 (its angle net checks a slightly SHRUNK square, so m is
           a lower bound on the true minimum, which is the safe direction).  Then w/m(s') is a
           cover of weight  U(s') = sum(w)/m(s').

The LP in between is only heuristic (exact points, exact placements, no erosion/dilation): it is
the finite restriction  min sum w  s.t.  w(Q) >= 1 for Q in S  over a finite point set P, grown by
cutting planes (placements on an (eta, dt) lattice with w(Q) < 1, plus the exact worst placements
reported by the verifier) and column generation (points where the dual coverage exceeds 1, found
by the same subdivision that computes M).  Rigor never depends on the LP: only on M and on the
verifier.

Usage
    python3 search/nu_f.py run  S  TLIMIT_SEC  TAG  [eta dt]
        -> runs/nuf_TAG.log, runs/nuf_TAG_best.txt (certificate), runs/nuf_TAG.json
    python3 search/nu_f.py check CERT  s1 s2 ...  [--N 6000] [--threads 8]
        -> prints U(s') = total/m(s') for each s'
"""
import sys, os, math, time, json, subprocess
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lp_search as LS                       # reuse scan_angle / admissible (unchanged)
VERIFY = os.path.join(HERE, '..', 'verify', 'target', 'release', 'verify')


# ----------------------------------------------------------------------------- geometry helpers
def rot(px, py, tm):
    ct, st = math.cos(tm), math.sin(tm)
    return px * ct + py * st, -px * st + py * ct


def d4(x, y, s):
    return [(x, y), (s - x, y), (x, s - y), (s - x, s - y), (y, x), (s - y, x), (y, s - x), (s - y, s - x)]


class ExactLP:
    """min sum_k |orbit_k| x_k  s.t.  for every placement r in S:  sum_k A_rk x_k >= 1,
    A_rk = # atoms of orbit k inside the CLOSED unit square at placement r (float tol 1e-9)."""

    def __init__(self, s, D0):
        self.s = s; self.D0 = D0                     # atoms live on the lattice Z^2 / D0
        self.orbits = []; self.okey = {}; self._ca = None
        self.rows = []; self.rkey = set()
        self.R = []; self.C = []; self.V = []

    # -- columns
    def add_point(self, x, y):
        X = int(round(x * self.D0)); Y = int(round(y * self.D0)); SD = int(round(self.s * self.D0))
        X = min(max(X, 0), SD); Y = min(max(Y, 0), SD)
        imgs = sorted(set([(X, Y), (SD - X, Y), (X, SD - Y), (SD - X, SD - Y), (Y, X), (SD - Y, X), (Y, SD - X), (SD - Y, SD - X)]))
        key = imgs[0]
        if key in self.okey: return False
        k = len(self.orbits); self.okey[key] = k
        Pk = np.array(imgs, dtype=float) / self.D0
        self.orbits.append(Pk); self._ca = None
        if self.rows:
            Q = np.array(self.rows); tot = np.zeros(len(Q))
            for (px, py) in Pk:
                q0 = px * Q[:, 3] + py * Q[:, 4]; q1 = -px * Q[:, 4] + py * Q[:, 3]
                tot += (np.maximum(np.abs(q0 - Q[:, 5]), np.abs(q1 - Q[:, 6])) <= 0.5 + 1e-9)
            hit = np.nonzero(tot)[0]
            if hit.size:
                self.R.append(hit.astype(np.int32)); self.C.append(np.full(hit.size, k, dtype=np.int32)); self.V.append(tot[hit])
        return True

    def atoms(self):
        if self._ca is None or self._ca[2] != len(self.orbits):
            P = np.concatenate(self.orbits); own = np.concatenate([np.full(len(Q), k) for k, Q in enumerate(self.orbits)])
            self._ca = (P, own, len(self.orbits))
        return self._ca[0], self._ca[1]

    # -- rows: exact admissible placements (cx, cy, tm)
    def add_rows(self, pl):
        P, own = self.atoms(); added = 0; s = self.s
        for cx, cy, tm in pl:
            ct, st = math.cos(tm), math.sin(tm); wid = abs(ct) + abs(st)
            if not (wid / 2 - 1e-12 <= cx <= s - wid / 2 + 1e-12 and wid / 2 - 1e-12 <= cy <= s - wid / 2 + 1e-12): continue
            key = (round(cx * 2e5), round(cy * 2e5), round(tm * 2e5))
            if key in self.rkey: continue
            self.rkey.add(key)
            u0, u1 = rot(cx, cy, tm)
            r = len(self.rows); self.rows.append((cx, cy, tm, ct, st, u0, u1))
            q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
            ok = np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= 0.5 + 1e-9
            cnt = np.bincount(own[ok], minlength=len(self.orbits)); idx = np.nonzero(cnt)[0]
            if idx.size:
                self.R.append(np.full(idx.size, r, dtype=np.int32)); self.C.append(idx.astype(np.int32)); self.V.append(cnt[idx].astype(float))
            added += 1
        return added

    def prune(self, keep):
        idx = np.full(len(self.rows), -1, dtype=np.int64); idx[np.nonzero(keep)[0]] = np.arange(int(keep.sum()))
        R = np.concatenate(self.R); C = np.concatenate(self.C); V = np.concatenate(self.V)
        nr = idx[R]; sel = nr >= 0
        self.R = [nr[sel].astype(np.int32)]; self.C = [C[sel]]; self.V = [V[sel]]
        self.rows = [rw for rw, k in zip(self.rows, keep) if k]
        self.rkey = set((round(rw[0] * 2e5), round(rw[1] * 2e5), round(rw[2] * 2e5)) for rw in self.rows)

    def solve(self):
        R = np.concatenate(self.R); C = np.concatenate(self.C); V = np.concatenate(self.V)
        A = sp.coo_matrix((V, (R, C)), shape=(len(self.rows), len(self.orbits))).tocsr()
        sizes = np.array([len(P) for P in self.orbits], dtype=float)
        res = linprog(c=sizes, A_ub=-A, b_ub=-np.ones(len(self.rows)), bounds=(0, None), method='highs')
        if not res.success: return None
        return res.fun, res.x, np.maximum(-res.ineqlin.marginals, 0.0)


# ----------------------------------------------------------------------------- rigorous M (max closed coverage)
def max_coverage(rows, y, s, h0=0.02, depth=8, want=400, budget=600_000):
    """Upper bound M on  max_{p in C} sum_r y_r [p in closed Q_r]  by adaptive subdivision.
    A cell of half-side h centred at g is dominated by the square dilated by h(|cos|+|sin|) in the
    square's own frame (any p in the cell has |R(p)-R(g)|_inf <= h(|cos t|+|sin t|)).
    Returns (M, worst_points) where worst_points are leaf centres with the highest dilated coverage."""
    rows = np.asarray(rows); y = np.asarray(y)
    sup = y > 1e-12; rows = rows[sup]; y = y[sup]
    if len(y) == 0: return 0.0, np.zeros((0, 2))
    # the LP's dual constraint bounds only the ORBIT-AVERAGED coverage, i.e. the coverage of the
    # D4-symmetrised measure (mass y_r/8 on each of the 8 images of placement r); bound that one.
    cx0, cy0, t0 = rows[:, 0], rows[:, 1], rows[:, 2]
    imgs = [(cx0, cy0, t0), (s - cx0, cy0, -t0), (cx0, s - cy0, -t0), (s - cx0, s - cy0, t0),
            (cy0, cx0, -t0), (s - cy0, cx0, t0), (cy0, s - cx0, t0), (s - cy0, s - cx0, -t0)]
    R = []
    for a, b, th in imgs:
        ct_, st_ = np.cos(th), np.sin(th)
        R.append(np.c_[a, b, th, ct_, st_, a * ct_ + b * st_, -a * st_ + b * ct_])
    rows = np.concatenate(R); y = np.tile(y / 8.0, 8)
    n0 = int(math.ceil(s / h0)); h = s / n0 / 2
    gx = (np.arange(n0) + 0.5) * (s / n0)
    G0, G1 = np.meshgrid(gx, gx, indexing='ij'); cells = np.c_[G0.ravel(), G1.ravel()]
    cx, cy, tm, ct, st, u0, u1 = [rows[:, i] for i in range(7)]
    wid = (np.abs(ct) + np.abs(st))
    certified = 0.0; worst = np.zeros((0, 2)); M = 0.0
    for lev in range(depth + 1):
        o = np.argsort(cells[:, 0]); cells = cells[o]; xs = cells[:, 0]
        cov = np.zeros(len(cells))
        dil = h * wid + 1e-9
        reach = 0.5 * math.sqrt(2) + dil
        a = np.searchsorted(xs, cx - reach, 'left'); b = np.searchsorted(xs, cx + reach, 'right')
        for r in range(len(y)):
            if b[r] <= a[r]: continue
            px = cells[a[r]:b[r], 0]; py = cells[a[r]:b[r], 1]
            q0 = px * ct[r] + py * st[r]; q1 = -px * st[r] + py * ct[r]
            ins = np.maximum(np.abs(q0 - u0[r]), np.abs(q1 - u1[r])) <= 0.5 + dil[r]
            cov[a[r]:b[r]] += y[r] * ins
        bad = cov > 1.0 + 1e-9
        if lev == 0: M = cov.max()
        good = ~bad
        if good.any(): certified = max(certified, cov[good].max())
        if not bad.any():
            M = certified; break
        if lev == depth or 4 * int(bad.sum()) > budget:
            # stop here: dilated coverage at this level is already a valid upper bound everywhere
            M = max(certified, cov[bad].max())
            o2 = np.argsort(cov[bad])[::-1][:want]; worst = cells[bad][o2]
            break
        # refine bad cells
        bc = cells[bad]; h = h / 2
        off = np.array([[-h, -h], [h, -h], [-h, h], [h, h]])
        cells = (bc[:, None, :] + off[None, :, :]).reshape(-1, 2)
        cells = cells[(cells[:, 0] >= 0) & (cells[:, 0] <= s) & (cells[:, 1] >= 0) & (cells[:, 1] <= s)]
    return M, worst


# ----------------------------------------------------------------------------- certificates & verifier
def write_cert_scaled(P, w, s, D0, path, s_target=None, WD=10 ** 7):
    """Write atoms (on lattice Z^2/D0, container s = SD/D0) as an exact certificate, scaled to
    container s_target (points scale with the container).  Returns (exported total, actual s')."""
    SD = int(round(s * D0)); X = np.round(P[:, 0] * D0).astype(np.int64); Y = np.round(P[:, 1] * D0).astype(np.int64)
    assert np.max(np.abs(P[:, 0] * D0 - X)) < 1e-6
    keep = w > 1e-10; X = X[keep]; Y = Y[keep]; W = np.ceil(w[keep] * WD).astype(np.int64)
    if s_target is None: s_target = s
    MUL = 1000
    Dp = int(round(SD * MUL / s_target))                 # D' with  s' = SD*MUL / D'
    sp_ = Fraction(SD * MUL, Dp)
    with open(path, 'w') as f:
        f.write(f"{sp_.numerator} {sp_.denominator}\n{Dp}\n{WD}\n{len(X)}\n")
        for a, b, c in zip(X * MUL, Y * MUL, W): f.write(f"{a} {b} {c}\n")
    return W.sum() / WD, float(sp_)


def read_cert(path):
    t = open(path).read().split(); it = iter(t)
    sn, sd, D, WD, m = [int(next(it)) for _ in range(5)]
    rows = np.array([[int(next(it)), int(next(it)), int(next(it))] for _ in range(m)], dtype=np.int64)
    return sn, sd, D, WD, rows


def rescale_cert(path, s_target, out):
    """Rewrite an existing certificate at container s_target (integer coords times 1000)."""
    sn, sd, D, WD, rows = read_cert(path)
    K = Fraction(sn, sd) * D; assert K.denominator == 1; K = int(K)   # s*D = container in grid units
    MUL = 1000; Dp = int(round(K * MUL / s_target)); sp_ = Fraction(K * MUL, Dp)
    with open(out, 'w') as f:
        f.write(f"{sp_.numerator} {sp_.denominator}\n{Dp}\n{WD}\n{len(rows)}\n")
        for a, b, c in rows: f.write(f"{a * MUL} {b * MUL} {c}\n")
    return rows[:, 2].sum() / WD, float(sp_)


def run_verifier(path, N=2000, threads=4, topk=0, sep=None, n=12):
    args = [VERIFY, path, str(n), str(N), str(threads), str(topk)] + ([sep] if sep else [])
    r = subprocess.run(args, capture_output=True, text=True, timeout=3600)
    ln = [l for l in r.stdout.split('\n') if l.startswith('min covered')]
    if not ln: raise RuntimeError(r.stdout[-500:] + r.stderr[-500:])
    frac = ln[0].split('=')[1].strip().split()[0]; a, b = frac.split('/')
    return Fraction(int(a), int(b))


# ----------------------------------------------------------------------------- separation
def separate(P, w, s, thetas, eta, perang=120, B=10, pool=None, thr=1.0):
    """placements on the (eta, dt) lattice where the CLOSED unit square captures < 1; worst per BxB block."""
    if pool is not None and len(thetas) > 8:
        k = pool._processes; chunks = [list(thetas[i::k]) for i in range(k)]
        outs = pool.starmap(separate, [(P, w, s, ch, eta, perang, B, None, thr) for ch in chunks])
        return min(o[0] for o in outs), sum(o[1] for o in outs), sum([o[2] for o in outs], [])
    pend = []; gmin = 9e9; nviol = 0
    for tm in thetas:
        r = LS.scan_angle(P, w, s, tm, 0.5 + 1e-9, eta)
        if r is None: continue
        vals, g0, g1, box = r
        G0, G1 = np.meshgrid(g0, g1, indexing='ij')
        ok = LS.admissible(G0.ravel(), G1.ravel(), box, -1e-9)
        v = np.where(ok, vals.ravel(), 9e9)
        gmin = min(gmin, v.min()); bad = v < thr - 1e-7; nviol += int(bad.sum())
        if not bad.any(): continue
        VV = v.reshape(vals.shape); n0, n1 = VV.shape
        p0 = (-n0) % B; p1 = (-n1) % B
        VP = np.pad(VV, ((0, p0), (0, p1)), constant_values=9e9)
        bl = VP.reshape(VP.shape[0] // B, B, VP.shape[1] // B, B).transpose(0, 2, 1, 3).reshape(-1, B * B)
        mn = bl.min(axis=1); am = bl.argmin(axis=1)
        sb = np.nonzero(mn < thr - 1e-7)[0]
        if len(sb) > perang: sb = sb[np.argsort(mn[sb])[:perang]]
        nb1 = VP.shape[1] // B
        bi, bj = np.divmod(sb, nb1); ii, jj = np.divmod(am[sb], B)
        r0 = bi * B + ii; r1 = bj * B + jj
        okk = (r0 < n0) & (r1 < n1); sel = r0[okk] * n1 + r1[okk]
        U0 = G0.ravel()[sel]; U1 = G1.ravel()[sel]
        lo, hi, ct, st = box
        cx = U0 * ct - U1 * st; cy = U0 * st + U1 * ct
        pend += [(float(a), float(b), float(tm)) for a, b in zip(cx, cy)]
    return gmin, nviol, pend


# ----------------------------------------------------------------------------- main loop
def run(s_frac, tlimit, tag, eta=0.01, dt=0.01, log=print, nproc=3):
    import multiprocessing as mp
    pool = mp.get_context('fork').Pool(nproc) if nproc > 1 else None
    s = float(s_frac)
    D0 = s_frac.denominator * int(math.ceil(10000 / s_frac.denominator))          # lattice ~1e-4
    t0 = time.time(); m = ExactLP(s, D0)
    thetas = np.linspace(0.0, math.pi / 4, int(round(math.pi / 4 / dt)) + 1)
    # initial columns: coarse grid (as lp_search), plus lattice points at 1/4 from the sides
    g = np.arange(0.0625, s, 0.125)
    for x in g:
        for yv in g: m.add_point(x, yv)
    P, own = m.atoms()
    # initial rows: a sparse sample of placements at 6 angles
    for tm in thetas[::max(1, len(thetas) // 6)]:
        gmin, nv, pend = separate(P, np.zeros(len(P)), s, [tm], eta, perang=100000, B=40)
        m.add_rows(pend)
    best = dict(L=0.0, U=9e9, Lit=None, Uit=None)
    hist = []
    sep_path = f"runs/nuf_{tag}_sep.txt"
    for it in range(10000):
        out = m.solve()
        if out is None: log("LP failed"); break
        val, x, y = out
        P, own = m.atoms(); w = x[own]
        # --- rigorous lower bound
        M, worst = max_coverage(np.array(m.rows), y, s)
        Lb = y.sum() / M if M > 0 else 0.0
        # --- cuts from the lattice scan
        gmin, nviol, pend = separate(P, w, s, thetas, eta, pool=pool)
        # --- exact worst placements from the verifier + exact min (rigorous U at this s)
        exact = None; Ub = None
        try:
            tw, sp_ = write_cert_scaled(P, w, s, D0, f"runs/nuf_{tag}_cur.txt")
            mv = run_verifier(f"runs/nuf_{tag}_cur.txt", N=2000, threads=4, topk=6, sep=sep_path)
            exact = float(mv)
            if mv > 0: Ub = tw / float(mv)
            pts = []
            for line in open(sep_path):
                q = line.split()
                if len(q) == 4 and float(q[0]) < 1.0: pts.append((float(q[2]), float(q[3]), float(q[1])))
            os.remove(sep_path)
        except Exception as e:
            log(f"   [verifier: {e}]"); pts = []
        if Lb > best['L']: best['L'] = Lb; best['Lit'] = it
        if Ub is not None and Ub < best['U']:
            best['U'] = Ub; best['Uit'] = it
            os.replace(f"runs/nuf_{tag}_cur.txt", f"runs/nuf_{tag}_best.txt")
        hist.append(dict(it=it, t=time.time() - t0, LP=val, M=M, L=Lb, exact=exact, U=Ub, rows=len(m.rows), orbits=len(m.orbits), gmin=gmin, nviol=nviol))
        log(f"it{it} LP={val:.5f} M={M:.5f} L={Lb:.5f} exact_min={exact} U={Ub if Ub is None else round(Ub,5)} "
            f"rows={len(m.rows)} orb={len(m.orbits)} atoms={len(P)} latmin={gmin:.4f} viol={nviol} cuts={len(pend)}+{len(pts)} t={time.time()-t0:.0f}s")
        json.dump(dict(s=str(s_frac), eta=eta, dt=dt, D0=D0, best=best, hist=hist), open(f"runs/nuf_{tag}.json", 'w'), indent=1)
        if time.time() - t0 > tlimit: break
        # --- column generation at the worst coverage leaves
        nnew = 0
        for (px, py) in worst[:300]: nnew += m.add_point(px, py)
        # --- prune non-binding rows, add cuts
        if len(m.rows) > 40000:
            keep = y > 1e-12; keep[-12000:] = True; m.prune(keep)
        m.add_rows(pend); m.add_rows(pts[:4000])
        if nviol == 0 and not pts and nnew == 0 and M < 1 + 1e-7:
            log("converged"); break
    return best


# ----------------------------------------------------------------------------- post-hoc rigorous lower bound
def lower_from_cert(cert, tag, h=0.01, rounds=7, eta=0.01, dt=0.01, log=print, nproc=4, K=3000):
    """Rigorous L(s) from a certificate's POINT SET: (1) LP dual over that point set and a broad row
    set (lattice placements with w(Q) < 1.05, plus the verifier's exact worst placements);
    (2) re-optimise the packing on the dual's support against every grid point (pitch h) of the
    D4 fundamental domain, closed squares; (3) certify with max_coverage, adding the worst leaves
    as further constraints, a few rounds.  Every round yields a valid L = sum(mu)/M."""
    import multiprocessing as mp
    t0 = time.time()
    sn, sd, D, WD, rows = read_cert(cert); s_frac = Fraction(sn, sd); s = float(s_frac)
    m = ExactLP(s, D)
    for X, Y, W in rows: m.add_point(X / D, Y / D)
    P, own = m.atoms()
    wat = np.zeros(len(P)); keyw = {}
    for X, Y, W in rows: keyw[(int(X), int(Y))] = W / WD
    for i, (px, py) in enumerate(P): wat[i] = keyw.get((int(round(px * D)), int(round(py * D))), 0.0)
    thetas = np.linspace(0.0, math.pi / 4, int(round(math.pi / 4 / dt)) + 1)
    pool = mp.get_context('fork').Pool(nproc)
    gmin, nviol, pend = separate(P, wat, s, thetas, eta, perang=400, B=6, pool=pool, thr=1.05)
    pool.close()
    m.add_rows(pend)
    sep_path = f"runs/lower_{tag}_sep.txt"
    try:
        mv = run_verifier(cert, N=2000, threads=4, topk=8, sep=sep_path)
        pts = [(float(q[2]), float(q[3]), float(q[1])) for q in (l.split() for l in open(sep_path)) if len(q) == 4]
        os.remove(sep_path); m.add_rows(pts)
    except Exception as e: log(f"[verifier: {e}]")
    out = m.solve(); val, x, y = out
    log(f"[{tag}] s={s:.6f} atoms={len(P)} rows={len(m.rows)} cover-LP over cert points = {val:.5f} (cert total {rows[:,2].sum()/WD:.5f}) t={time.time()-t0:.0f}s")
    # support for the packing: the dual's support plus the K tightest placements of the cover LP
    Rm = np.concatenate(m.R); Cm = np.concatenate(m.C); Vm = np.concatenate(m.V)
    Am = sp.coo_matrix((Vm, (Rm, Cm)), shape=(len(m.rows), len(m.orbits))).tocsr()
    slack = Am @ x - 1.0
    sup = np.nonzero(y > 1e-9)[0]
    extra = np.argsort(slack)[:K]
    sup = np.unique(np.concatenate([sup, extra])); R = np.array(m.rows)[sup]
    log(f"[{tag}] packing support: {len(sup)} placements (dual support {(y > 1e-9).sum()}, K={K})")
    cx0, cy0, t0_ = R[:, 0], R[:, 1], R[:, 2]
    imgs = [(cx0, cy0, t0_), (s - cx0, cy0, -t0_), (cx0, s - cy0, -t0_), (s - cx0, s - cy0, t0_),
            (cy0, cx0, -t0_), (s - cy0, cx0, t0_), (cy0, s - cx0, t0_), (s - cy0, s - cx0, -t0_)]
    n0 = int(math.ceil(s / h)); gx = (np.arange(n0) + 0.5) * (s / n0)
    G0, G1 = np.meshgrid(gx, gx, indexing='ij'); G = np.c_[G0.ravel(), G1.ravel()]
    G = G[(G[:, 0] <= G[:, 1] + 1e-12) & (G[:, 1] <= s / 2 + 1e-12)]          # D4 fundamental domain
    best = dict(L=0.0, M=None, mass=None); hist = []
    for rd in range(rounds):
        Rr, Cc, Vv = [], [], []
        for (a, b, th) in imgs:
            ct_, st_ = np.cos(th), np.sin(th); u0 = a * ct_ + b * st_; u1 = -a * st_ + b * ct_
            for r in range(len(sup)):
                q0 = G[:, 0] * ct_[r] + G[:, 1] * st_[r]; q1 = -G[:, 0] * st_[r] + G[:, 1] * ct_[r]
                ins = np.nonzero(np.maximum(np.abs(q0 - u0[r]), np.abs(q1 - u1[r])) <= 0.5 + 1e-9)[0]
                if ins.size: Rr.append(ins); Cc.append(np.full(ins.size, r)); Vv.append(np.full(ins.size, 0.125))
        A = sp.coo_matrix((np.concatenate(Vv), (np.concatenate(Rr), np.concatenate(Cc))), shape=(len(G), len(sup))).tocsr()
        res = linprog(c=-np.ones(len(sup)), A_ub=A, b_ub=np.ones(len(G)), bounds=(0, None), method='highs')
        if not res.success: log("packing LP failed"); break
        mu = res.x; mass = mu.sum()
        M, worst = max_coverage(R, mu, s, h0=0.02, depth=9, want=2000, budget=1_500_000)
        L = mass / M
        hist.append(dict(round=rd, grid=len(G), mass=mass, M=M, L=L, t=time.time() - t0))
        log(f"[{tag}] round{rd} grid={len(G)} packing mass={mass:.5f} M={M:.5f} -> L={L:.5f} t={time.time()-t0:.0f}s")
        if L > best['L']: best = dict(L=L, M=M, mass=mass, round=rd)
        if M < 1 + 1e-6 or len(worst) == 0: break
        G = np.concatenate([G, worst])
    json.dump(dict(cert=cert, s=str(s_frac), h=h, coverLP=val, best=best, hist=hist), open(f"runs/lower_{tag}.json", 'w'), indent=1)
    return best


if __name__ == "__main__":
    if sys.argv[1] == 'run':
        sf = Fraction(sys.argv[2]); tl = float(sys.argv[3]); tag = sys.argv[4]
        eta = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
        dt = float(sys.argv[6]) if len(sys.argv) > 6 else 0.01
        os.makedirs('runs', exist_ok=True)
        lf = open(f"runs/nuf_{tag}.log", 'a')
        def log(msg):
            print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
        log(f"nu_f run s={sf} ({float(sf):.6f}) eta={eta} dt={dt} tlimit={tl}")
        best = run(sf, tl, tag, eta, dt, log)
        log(f"RESULT s={float(sf):.6f}  L={best['L']:.5f} (it{best['Lit']})  U(s)={best['U']:.5f} (it{best['Uit']})")
    elif sys.argv[1] == 'lower':
        cert = sys.argv[2]; tag = sys.argv[3]; h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
        K = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
        os.makedirs('runs', exist_ok=True); lf = open(f"runs/lower_{tag}.log", 'a')
        def log(msg):
            print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
        best = lower_from_cert(cert, tag, h=h, log=log, K=K)
        log(f"LOWER {cert} L={best['L']:.5f} M={best.get('M')} mass={best.get('mass')}")
    elif sys.argv[1] == 'check':
        cert = sys.argv[2]; N = 6000; th = 8; targets = []
        a = sys.argv[3:]
        i = 0
        while i < len(a):
            if a[i] == '--N': N = int(a[i + 1]); i += 2
            elif a[i] == '--threads': th = int(a[i + 1]); i += 2
            else: targets.append(float(a[i])); i += 1
        tmp = cert + '.scaled.tmp'
        for st in targets:
            tot, sp_ = rescale_cert(cert, st, tmp)
            mv = run_verifier(tmp, N=N, threads=th)
            U = tot / float(mv) if mv > 0 else float('inf')
            print(f"EVAL {cert} s'={sp_:.6f} total={tot:.7f} min={float(mv):.7f} U={U:.5f}", flush=True)
        os.remove(tmp)
