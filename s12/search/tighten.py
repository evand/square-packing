#!/usr/bin/env python3
"""Re-optimise / sparsify the weights of an existing certificate at fixed point locations.

The point set of a certificate is kept (as D4 orbits, so the verifier's [0,45deg] reduction
still applies) and only the weights are recomputed by a cutting-plane loop whose separation
oracle is the EXACT verifier itself:

    repeat
        x  <- LP:  min sum_k |orbit_k| x_k   s.t.  A x >= 1 + m_lp,  x >= 0
              (A_rk = number of atoms of orbit k inside the closed square of row r)
        export x/(1+m_probe) at weight denominator 1e12 (rounded DOWN), run
              verify <cert> 12 N <threads> <topk> <witness-file>
        every witness (a placement of the verifier's shrunk sigma_k-square, at bin angle
              theta_k, whose captured weight is < 1) becomes a new LP row with exactly that
              square (half-side sigma_k/2, computed with the verifier's own integer formula)
    until the probe reports no violated placement

Because the probe weights are x/(1+m_probe) < x and rounded down, "no violation" means the
true weights x cover >= 1+m_probe (up to 1e-9) on the whole N-net, and the final export
(rounded UP at W = 1e7) can only cover more.  The final file is then handed to the verifier
again at N and at 2N; the second net is not implied by the first (its squares are shrunk by
less but sit at other angles), so its witnesses are fed back too until both are clean.

Container.  The integer coordinates never change; the container is K/D where K = s*D of the
input (7840 for the shipped file).  `--Dp` picks another D (smaller D = bigger container),
optionally with `--mul` (coordinates multiplied by mul, so D can be chosen to 1/mul).

Modes
    python3 search/tighten.py reopt    CERT TAG [--Dp D --mul M --N 6000 --topk 6 --margin 2e-6]
        -> runs/tight_TAG.txt   (weights re-optimised, W = 1e7)
    python3 search/tighten.py sparsify CERT TAG --budget 11.99 [--rounds 8 ...]
        -> runs/sparse_TAG.txt  (reweighted-L1 heuristic: min sum c_k |orbit_k| x_k with
           c_k = 1/(x_k+eps) under total weight <= budget, then min total weight on the support)
    python3 search/tighten.py scan     CERT TAG --Dps 1994,1992,...   (reopt at each D, report)
    python3 search/tighten.py reopt    CERT TAG --colgen 30 [--cg-want 200]
        -> as reopt, but for the first 30 iterations new points are also generated where the
           dual (packing) measure's D4-symmetrised closed-square coverage exceeds 1 (exact
           pricing on a --cg-pitch grid, refined on the integer grid K/D); the point set grows
           and the LP value can drop below the fixed-point optimum.

Rows from the lattice separation of nu_f.py (closed unit squares on an (eta, dt) lattice with
captured weight < thr) are used as a warm start; they are valid constraints (h = 1/2).
"""
import sys, os, math, time, subprocess, argparse, json
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lp_search as LS
import nu_f as NF

REPO = os.path.dirname(HERE)
VERIFY = os.path.join(REPO, 'verify', 'target', 'release', 'verify')
NPROC = os.cpu_count() or 4
HIGHS = {'random_seed': 0}
TOL = 1e-9          # atoms within TOL of the square boundary are NOT counted (safe direction)


def _lp(c, A, b, A_eq=None, b_eq=None, A_ub2=None, b_ub2=None):
    """min c.x  s.t. A x >= b, [A_ub2 x <= b_ub2], x >= 0"""
    Aub = -A; bub = -b
    if A_ub2 is not None:
        Aub = sp.vstack([Aub, A_ub2]).tocsr(); bub = np.concatenate([bub, b_ub2])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return linprog(c=c, A_ub=Aub, b_ub=bub, bounds=(0, None), method='highs', options=dict(HIGHS))


def read_cert(path):
    t = open(path).read().split(); it = iter(t)
    sn, sd, D, WD, m = [int(next(it)) for _ in range(5)]
    rows = np.array([[int(next(it)), int(next(it)), int(next(it))] for _ in range(m)], dtype=np.int64)
    return sn, sd, D, WD, rows


def sigma_k(N, k):
    """the verifier's shrunk side for bin k of the N-net (integer formula, rounded down to 1e-6)"""
    c0, s0, g0 = N * N - k * k, 2 * k * N, N * N + k * k
    k2 = k + 1
    c1, s1, g1 = N * N - k2 * k2, 2 * k2 * N, N * N + k2 * k2
    cd = c0 * c1 + s0 * s1; sd = c0 * s1 - s0 * c1
    return (g0 * g1 * 10 ** 6) // (cd + sd) / 10 ** 6


class Model:
    """orbits of integer points (grid K), container s = K/D; rows = (cx, cy, theta, h)."""

    def __init__(self, K, D, orbits_int):
        self.K = K; self.D = D; self.s = K / D
        self.orb_int = list(orbits_int)                                   # list of int arrays (n_k, 2)
        self.orbits = [o / D for o in orbits_int]                   # float positions
        self.sizes = np.array([len(o) for o in self.orbits], dtype=float)
        self.P = np.concatenate(self.orbits)
        self.own = np.concatenate([np.full(len(o), k) for k, o in enumerate(self.orbits)])
        self.rows = []; self.rkey = set(); self.R = []; self.C = []; self.V = []
        self.okey = {tuple(o[0]): k for k, o in enumerate(orbits_int)}

    def add_rows(self, pl):
        """pl: iterable of (cx, cy, theta, h).  Returns number added."""
        P, own = self.P, self.own; added = 0
        for cx, cy, tm, h in pl:
            key = (round(cx * 1e7), round(cy * 1e7), round(tm * 1e8), round(h * 1e7))
            if key in self.rkey: continue
            self.rkey.add(key)
            ct, st = math.cos(tm), math.sin(tm)
            u0 = cx * ct + cy * st; u1 = -cx * st + cy * ct
            q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
            ok = np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= h - TOL
            cnt = np.bincount(own[ok], minlength=len(self.orbits)); idx = np.nonzero(cnt)[0]
            r = len(self.rows); self.rows.append((cx, cy, tm, h))
            if idx.size:
                self.R.append(np.full(idx.size, r, dtype=np.int32)); self.C.append(idx.astype(np.int32)); self.V.append(cnt[idx].astype(float))
            added += 1
        return added

    def add_orbit(self, X, Y):
        """add the D4 orbit of the integer point (X, Y); entries against all existing rows."""
        K = self.K; X = min(max(int(X), 0), K); Y = min(max(int(Y), 0), K)
        imgs = sorted(set([(X, Y), (K - X, Y), (X, K - Y), (K - X, K - Y), (Y, X), (K - Y, X), (Y, K - X), (K - Y, K - X)]))
        key = imgs[0]
        if key in self.okey: return False
        k = len(self.orbits); self.okey[key] = k
        oi = np.array(imgs, dtype=np.int64); self.orb_int.append(oi); of = oi / self.D; self.orbits.append(of)
        self.sizes = np.append(self.sizes, float(len(of)))
        self.P = np.concatenate([self.P, of]); self.own = np.concatenate([self.own, np.full(len(of), k)])
        if self.rows:
            Q = np.array(self.rows); ct = np.cos(Q[:, 2]); st = np.sin(Q[:, 2])
            u0 = Q[:, 0] * ct + Q[:, 1] * st; u1 = -Q[:, 0] * st + Q[:, 1] * ct; h = Q[:, 3]
            tot = np.zeros(len(Q))
            for (px, py) in of:
                q0 = px * ct + py * st; q1 = -px * st + py * ct
                tot += (np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= h - TOL)
            hit = np.nonzero(tot)[0]
            if hit.size:
                self.R.append(hit.astype(np.int32)); self.C.append(np.full(hit.size, k, dtype=np.int32)); self.V.append(tot[hit])
        return True

    def matrix(self):
        if not self.R: return sp.csr_matrix((len(self.rows), len(self.orbits)))
        R = np.concatenate(self.R); C = np.concatenate(self.C); V = np.concatenate(self.V)
        return sp.coo_matrix((V, (R, C)), shape=(len(self.rows), len(self.orbits))).tocsr()

    def prune(self, keep):
        idx = np.full(len(self.rows), -1, dtype=np.int64); idx[np.nonzero(keep)[0]] = np.arange(int(keep.sum()))
        R = np.concatenate(self.R); C = np.concatenate(self.C); V = np.concatenate(self.V)
        nr = idx[R]; sel = nr >= 0
        self.R = [nr[sel].astype(np.int32)]; self.C = [C[sel]]; self.V = [V[sel]]
        self.rows = [rw for rw, k in zip(self.rows, keep) if k]
        self.rkey = set((round(a * 1e7), round(b * 1e7), round(t * 1e8), round(h * 1e7)) for a, b, t, h in self.rows)

    def solve(self, margin, cost=None, budget=None, fixed_zero=None):
        """min cost.x  s.t.  A x >= 1+margin  [, sizes.x <= budget] [, x_k = 0 for k in fixed_zero]"""
        A = self.matrix(); n = len(self.orbits)
        c = self.sizes.copy() if cost is None else np.asarray(cost, dtype=float)
        A2 = b2 = None
        if budget is not None:
            A2 = sp.csr_matrix(self.sizes.reshape(1, -1)); b2 = np.array([budget])
        if fixed_zero is not None and len(fixed_zero):
            Z = sp.csr_matrix((np.ones(len(fixed_zero)), (np.arange(len(fixed_zero)), np.asarray(fixed_zero))), shape=(len(fixed_zero), n))
            A2 = Z if A2 is None else sp.vstack([A2, Z]).tocsr()
            b2 = np.zeros(len(fixed_zero)) if b2 is None else np.concatenate([b2, np.zeros(len(fixed_zero))])
        res = _lp(c, A, np.full(A.shape[0], 1.0 + margin), A_ub2=A2, b_ub2=b2)
        if not res.success: return None
        y = -res.ineqlin.marginals[:A.shape[0]]
        return res.fun, res.x, np.maximum(y, 0.0)


def build_model(cert, Dp=None, mul=1):
    sn, sd, D, WD, rows = read_cert(cert)
    K = Fraction(sn, sd) * D; assert K.denominator == 1; K = int(K)
    if Dp is None: Dp = D * mul
    K2 = K * mul
    okey = {}; orbs = []; w0 = []
    for X, Y, W in rows:
        X = int(X) * mul; Y = int(Y) * mul
        imgs = sorted(set([(X, Y), (K2 - X, Y), (X, K2 - Y), (K2 - X, K2 - Y), (Y, X), (K2 - Y, X), (Y, K2 - X), (K2 - Y, K2 - X)]))
        key = imgs[0]
        if key not in okey:
            okey[key] = len(orbs); orbs.append(np.array(imgs, dtype=np.int64)); w0.append(W / WD)
    m = Model(K2, Dp, orbs)
    return m, np.array(w0), Fraction(K2, Dp)


def export(m, x, path, WD=10 ** 7, factor=1.0, up=True):
    """write the orbits with x_k > 0 as a certificate at container K/D; weights x_k*factor over WD"""
    s = Fraction(m.K, m.D)
    lines = []; tot = 0
    for k, o in enumerate(m.orb_int):
        if x[k] <= 0: continue
        v = x[k] * factor * WD
        w = int(math.ceil(v - 1e-9)) if up else int(math.floor(v))
        if w <= 0: continue
        for X, Y in o: lines.append(f"{X} {Y} {w}"); tot += w
    with open(path, 'w') as f:
        f.write(f"{s.numerator} {s.denominator}\n{m.D}\n{WD}\n{len(lines)}\n" + "\n".join(lines) + "\n")
    return tot / WD, len(lines)


def run_verifier(path, N, topk=0, sep=None, threads=NPROC, n=12):
    args = [VERIFY, path, str(n), str(N), str(threads), str(topk)] + ([sep] if sep else [])
    r = subprocess.run(args, capture_output=True, text=True)
    ln = [l for l in r.stdout.split('\n') if l.startswith('min covered')]
    if not ln: raise RuntimeError(r.stdout[-500:] + r.stderr[-500:])
    frac = ln[0].split('=')[1].strip().split()[0]; a, b = frac.split('/')
    ok = ('VERIFIED' in r.stdout) and ('NOT VERIFIED' not in r.stdout)
    return Fraction(int(a), int(b)), ok


def read_witnesses(path, N):
    """verifier witness file -> rows (cx, cy, theta_k, sigma_k/2), canonical order"""
    out = []
    for line in open(path):
        q = line.split()
        if len(q) != 4: continue
        v, th, cx, cy = map(float, q)
        k = int(round(N * math.tan(th / 2)))
        out.append((v, th, cx, cy, k))
    out.sort()
    return [(cx, cy, th, sigma_k(N, k) / 2) for v, th, cx, cy, k in out]


def lattice_rows(m, w, eta=0.01, dt=0.01, thr=1.05, perang=400, B=6, nproc=8):
    import multiprocessing as mp
    thetas = np.linspace(0.0, math.pi / 4, int(round(math.pi / 4 / dt)) + 1)
    pool = mp.get_context('fork').Pool(nproc)
    gmin, nviol, pend = NF.separate(m.P, w, m.s, thetas, eta, perang=perang, B=B, pool=pool, thr=thr)
    pool.close(); pool.join()
    s = m.s; out = []
    for cx, cy, tm in pend:
        wid = abs(math.cos(tm)) + abs(math.sin(tm))
        if wid / 2 - 1e-12 <= cx <= s - wid / 2 + 1e-12 and wid / 2 - 1e-12 <= cy <= s - wid / 2 + 1e-12:
            out.append((cx, cy, tm, 0.5))
    return out


def price(m, y, pitch=0.01, want=200, ysup=1e-9):
    """exact reduced-cost pricing on a grid: the reduced cost of the orbit of p is
    |orbit| (1 - covsym(p)) with covsym(p) = (1/8) sum_g cov(g p), cov(q) = sum_r y_r [q in Q_r]
    (closed squares of the rows, half-side h_r).  Evaluates covsym on a pitch-grid of the D4
    fundamental domain, refines the best candidates on the integer grid 1/D around them, and
    returns integer points (X, Y) with covsym > 1, best first."""
    Q = np.array(m.rows); sup = y > ysup; Q = Q[sup]; yy = y[sup]
    if len(yy) == 0: return []
    s = m.s; cx0, cy0, t0, h0 = Q[:, 0], Q[:, 1], Q[:, 2], Q[:, 3]
    imgs = [(cx0, cy0, t0), (s - cx0, cy0, -t0), (cx0, s - cy0, -t0), (s - cx0, s - cy0, t0),
            (cy0, cx0, -t0), (s - cy0, cx0, t0), (cy0, s - cx0, t0), (s - cy0, s - cx0, -t0)]
    R = []
    for a, b, th in imgs:
        ct_, st_ = np.cos(th), np.sin(th)
        R.append(np.c_[a * ct_ + b * st_, -a * st_ + b * ct_, ct_, st_, h0, yy / 8.0])
    R = np.concatenate(R)                                   # u0 u1 ct st h y
    def cov_at(G):
        out = np.zeros(len(G))
        for i in range(0, len(R), 2000):
            B = R[i:i + 2000]
            q0 = G[:, 0:1] * B[:, 2] + G[:, 1:2] * B[:, 3]; q1 = -G[:, 0:1] * B[:, 3] + G[:, 1:2] * B[:, 2]
            ins = np.maximum(np.abs(q0 - B[:, 0]), np.abs(q1 - B[:, 1])) <= B[:, 4] + 1e-9
            out += ins @ B[:, 5]
        return out
    n0 = int(math.ceil(s / pitch)); gx = np.arange(n0 + 1) * (s / n0)
    G0, G1 = np.meshgrid(gx, gx, indexing='ij'); G = np.c_[G0.ravel(), G1.ravel()]
    G = G[(G[:, 0] <= G[:, 1] + 1e-12) & (G[:, 1] <= s / 2 + 1e-12)]
    c = cov_at(G)
    o = np.argsort(c)[::-1][:want]
    cand = []
    D = m.D; step = s / n0
    for i in o:
        if c[i] <= 1.0 + 1e-7: break
        px, py = G[i]
        # refine on the integer grid within +-step
        X0, Y0 = int(round(px * D)), int(round(py * D)); r_ = int(math.ceil(step * D))
        xs = np.arange(max(0, X0 - r_), min(m.K, X0 + r_) + 1); ys_ = np.arange(max(0, Y0 - r_), min(m.K, Y0 + r_) + 1)
        if len(xs) * len(ys_) > 4000:
            xs = xs[::max(1, len(xs) // 60)]; ys_ = ys_[::max(1, len(ys_) // 60)]
        GX, GY = np.meshgrid(xs, ys_, indexing='ij'); Gf = np.c_[GX.ravel(), GY.ravel()] / D
        cf = cov_at(Gf); j = int(np.argmax(cf))
        if cf[j] > 1.0 + 1e-7: cand.append((float(cf[j]), int(GX.ravel()[j]), int(GY.ravel()[j])))
    cand.sort(reverse=True)
    return cand


def tighten(m, tag, margin=2e-6, N=6000, topk=6, max_iters=200, log=print, cost=None, budget=None,
            fixed_zero=None, probe_margin=None, final_check=True, x0=None, prune_at=120000,
            colgen=0, cg_want=200, cg_pitch=0.01, n=12):
    """cutting-plane loop.  Returns (x, val, info).
    colgen=k: for the first k iterations also generate columns: the dual y is a (D4-symmetrised)
    packing measure on the rows; points where its closed-square coverage exceeds 1 have negative
    reduced cost |orbit|(1 - coverage) < 0.  `price` evaluates the symmetrised coverage on a
    pitch-grid, refines on the integer grid, and the best `cg_want` points are added."""
    if probe_margin is None: probe_margin = margin / 2
    tmp = f"runs/tight_{tag}_probe.txt"; sep = f"runs/tight_{tag}_sep.txt"
    t0 = time.time(); x = x0; val = None; it = 0
    while it < max_iters:
        out = m.solve(margin, cost=cost, budget=budget, fixed_zero=fixed_zero)
        if out is None:
            log(f"  it{it} LP infeasible/failed (rows={len(m.rows)})"); return None, None, dict(status='infeasible')
        val, x, y = out; tw = float(m.sizes @ x)
        export(m, x, tmp, WD=10 ** 12, factor=1.0 / (1.0 + probe_margin), up=False)
        mv, ok = run_verifier(tmp, N, topk=topk, sep=sep, n=n)
        wit = read_witnesses(sep, N) if os.path.exists(sep) else []
        if os.path.exists(sep): os.remove(sep)
        nz = int((x > 1e-12).sum()); natoms = int(m.sizes[x > 1e-12].sum())
        ncol = 0; Mcov = None
        if it < colgen:                                  # y belongs to the rows BEFORE this iteration's cuts
            cand = price(m, y, pitch=cg_pitch, want=cg_want)
            Mcov = cand[0][0] if cand else 1.0
            for (cv, X, Y) in cand:
                ncol += m.add_orbit(X, Y)
        added = m.add_rows(wit)
        log(f"  it{it} total={tw:.6f} obj={val:.6f} probe_min={float(mv):.7f} viol={len(wit)} new_rows={added} rows={len(m.rows)} support={nz} orbits/{natoms} atoms"
            + (f" dualcov={Mcov:.4f} new_orbits={ncol} orbits={len(m.orbits)}" if Mcov is not None else "") + f" t={time.time()-t0:.0f}s")
        if added == 0 and mv >= 1 and ncol == 0:
            break
        if len(m.rows) > prune_at:
            keep = y > 1e-12; keep[-30000:] = True; m.prune(keep); log(f"   pruned to {len(m.rows)} rows")
        it += 1
    return x, val, dict(iters=it, rows=len(m.rows), t=time.time() - t0)


def finalize(m, x, tag, out_path, N=6000, log=print, WD=10 ** 7, margin=2e-6, cost=None, budget=None, fixed_zero=None, x0=None, n=12):
    """export at WD (rounded up); check at N and 2N, feeding back 2N witnesses if needed."""
    for rnd in range(20):
        tw, npts = export(m, x, out_path, WD=WD, up=True)
        mv1, ok1 = run_verifier(out_path, N, n=n)
        sep = f"runs/tight_{tag}_sep2.txt"
        mv2, ok2 = run_verifier(out_path, 2 * N, topk=6, sep=sep, n=n)
        wit = read_witnesses(sep, 2 * N) if os.path.exists(sep) else []
        if os.path.exists(sep): os.remove(sep)
        log(f"  final[{rnd}] {out_path}: points={npts} total={tw:.7f} min@{N}={float(mv1):.7f} {'OK' if ok1 else 'FAIL'}  min@{2*N}={float(mv2):.7f} {'OK' if ok2 else 'FAIL'} (viol {len(wit)})")
        if ok1 and ok2: return tw, npts, mv1, mv2
        if tw >= n and mv1 >= 1 and mv2 >= 1:
            log(f"  final: covering holds at both nets but total weight {tw:.6f} >= {n}: not a certificate")
            return tw, npts, mv1, mv2
        m.add_rows(wit)
        x, val, info = tighten(m, tag, margin=margin, N=N, log=log, cost=cost, budget=budget, fixed_zero=fixed_zero, x0=x, n=n)
        if x is None: return None
    return None


def reopt(cert, tag, Dp=None, mul=1, N=6000, topk=6, margin=2e-6, log=print, warm=True, out=None, colgen=0, cg_want=200, cg_pitch=0.01, n=12):
    m, w0, s = build_model(cert, Dp, mul)
    log(f"[{tag}] {cert}: {len(m.orbits)} orbits / {len(m.P)} atoms, container {s} = {float(s):.7f}, input total {float(m.sizes @ w0):.6f}")
    t0 = time.time()
    if warm:
        lr = lattice_rows(m, w0[m.own], thr=1.05)
        m.add_rows(lr); log(f"[{tag}] warm start: {len(m.rows)} lattice rows ({time.time()-t0:.0f}s)")
    x, val, info = tighten(m, tag, margin=margin, N=N, topk=topk, log=log, colgen=colgen, cg_want=cg_want, cg_pitch=cg_pitch, n=n)
    if x is None: log(f"[{tag}] FAILED"); return None
    out = out or f"runs/tight_{tag}.txt"
    res = finalize(m, x, tag, out, N=N, log=log, margin=margin, n=n)
    if res is None: log(f"[{tag}] finalize FAILED"); return None
    tw, npts, mv1, mv2 = res
    log(f"[{tag}] RESULT s={s} ({float(s):.7f}) points={npts} total={tw:.7f} min@{N}={mv1} min@{2*N}={mv2} t={time.time()-t0:.0f}s")
    return dict(cert=out, s=str(s), sf=float(s), points=npts, total=tw, min1=str(mv1), min2=str(mv2), t=time.time() - t0, model=m, x=x)


def sparsify(cert, tag, budget, Dp=None, mul=1, N=6000, topk=6, margin=2e-6, rounds=8, eps=1e-4, log=print, out=None, seed_from=None, n=12):
    """reweighted-L1 under a total-weight budget; then min-weight on the support."""
    m, w0, s = build_model(cert, Dp, mul)
    log(f"[{tag}] sparsify {cert}: {len(m.orbits)} orbits / {len(m.P)} atoms, s={float(s):.7f}, budget {budget}")
    t0 = time.time()
    lr = lattice_rows(m, w0[m.own], thr=1.05); m.add_rows(lr)
    log(f"[{tag}] warm start: {len(m.rows)} lattice rows ({time.time()-t0:.0f}s)")
    # first: plain min-weight, to have a clean row set and a starting x
    x, val, info = tighten(m, tag, margin=margin, N=N, topk=topk, log=log, n=n)
    if x is None: return None
    best = None
    fixed = np.zeros(0, dtype=int)
    for rd in range(rounds):
        cost = m.sizes / (x + eps)
        x2, val2, info = tighten(m, tag, margin=margin, N=N, topk=topk, log=log, cost=cost, budget=budget, fixed_zero=fixed, n=n)
        if x2 is None:
            log(f"[{tag}] round {rd}: infeasible under budget with {len(fixed)} fixed zeros"); break
        x = x2
        sup = np.nonzero(x > 1e-9)[0]; zero = np.nonzero(x <= 1e-9)[0]
        fixed = zero                                    # once zero, stays zero (keeps rows consistent)
        log(f"[{tag}] round {rd}: support {len(sup)} orbits / {int(m.sizes[sup].sum())} atoms, total {float(m.sizes @ x):.6f}")
        if best is not None and len(sup) >= best[0] and rd >= 2 and len(sup) == best[0]: break
        best = (len(sup), x.copy(), fixed.copy())
    # polish: min total weight on the final support
    nsup, x, fixed = best
    x, val, info = tighten(m, tag, margin=margin, N=N, topk=topk, log=log, fixed_zero=fixed, n=n)
    if x is None: return None
    out = out or f"runs/sparse_{tag}.txt"
    res = finalize(m, x, tag, out, N=N, log=log, margin=margin, fixed_zero=fixed, n=n)
    if res is None: log(f"[{tag}] finalize FAILED"); return None
    tw, npts, mv1, mv2 = res
    log(f"[{tag}] RESULT s={s} ({float(s):.7f}) points={npts} total={tw:.7f} min@{N}={mv1} min@{2*N}={mv2} t={time.time()-t0:.0f}s")
    return dict(cert=out, s=str(s), sf=float(s), points=npts, total=tw, min1=str(mv1), min2=str(mv2), t=time.time() - t0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['reopt', 'sparsify', 'scan'])
    ap.add_argument('cert'); ap.add_argument('tag')
    ap.add_argument('--Dp', type=int, default=None); ap.add_argument('--mul', type=int, default=1)
    ap.add_argument('--Dps', default=None, help='comma-separated D values for scan')
    ap.add_argument('--N', type=int, default=6000); ap.add_argument('--topk', type=int, default=6)
    ap.add_argument('--margin', type=float, default=2e-6)
    ap.add_argument('--budget', type=float, default=11.99); ap.add_argument('--rounds', type=int, default=8)
    ap.add_argument('--eps', type=float, default=1e-4)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--out', default=None)
    ap.add_argument('--colgen', type=int, default=0, help='number of column-generation iterations (reopt)')
    ap.add_argument('--cg-want', type=int, default=200); ap.add_argument('--cg-pitch', type=float, default=0.01)
    ap.add_argument('--n', type=int, default=12, help='number of squares to rule out (only the weight threshold: verifier n and the "not a certificate" test)')
    a = ap.parse_args()
    HIGHS['random_seed'] = a.seed
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/tight_{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"tighten.py {' '.join(sys.argv[1:])}")
    if a.mode == 'reopt':
        r = reopt(a.cert, a.tag, Dp=a.Dp, mul=a.mul, N=a.N, topk=a.topk, margin=a.margin, log=log, out=a.out, colgen=a.colgen, cg_want=a.cg_want, cg_pitch=a.cg_pitch, n=a.n)
        if r: r.pop('model'); r.pop('x'); json.dump(r, open(f"runs/tight_{a.tag}.json", 'w'), indent=1)
    elif a.mode == 'sparsify':
        r = sparsify(a.cert, a.tag, a.budget, Dp=a.Dp, mul=a.mul, N=a.N, topk=a.topk, margin=a.margin, rounds=a.rounds, eps=a.eps, log=log, out=a.out, n=a.n)
        if r: json.dump(r, open(f"runs/tight_{a.tag}.json", 'w'), indent=1)
    elif a.mode == 'scan':
        res = []
        for Dp in [int(v) for v in a.Dps.split(',')]:
            r = reopt(a.cert, f"{a.tag}_D{Dp}", Dp=Dp, mul=a.mul, N=a.N, topk=a.topk, margin=a.margin, log=log, n=a.n)
            if r: r.pop('model'); r.pop('x'); res.append(r)
            json.dump(res, open(f"runs/tight_{a.tag}_scan.json", 'w'), indent=1)
        for r in res: log(f"SCAN s={r['sf']:.7f} points={r['points']} total={r['total']:.7f} min={r['min1']}")


if __name__ == '__main__':
    main()
