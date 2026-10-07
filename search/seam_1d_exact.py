"""Exact (rational) certificate checker for the 1D seam covering problem (FRIEDMAN.md section 8, SEAM_1D.md).

Measure on the circle [0,1): atoms w_k at x = k/N (k = 0 is the seam, g = w_0) plus a uniform density eta.
Claim: F(th, cx) = eta + sum_k w_k G_th(k/N - cx) >= 1 for all th in [0, 45 deg], cx in [0, 1), where
G_th(u) = sum_j V_th(u + j) and V_th(d) is the vertical chord of the closed unit square (centre 0, angle th) at offset d.

Exactness: th is parametrised by t = tan(th/2) (c, s rational); every quantity in the checker is a Python int or
Fraction.  Floats appear only in `prepare` (choosing the measure) and in `test`.  The B&B's split rule and search order
are heuristics; they cannot affect soundness, since every leaf is certified by an exact bound.

Subcommands
  prepare  --npy FILE | --N N --e E   [--cuts R] [--margin M] --out FILE.json
           float LP (search/seam_1d.py) + cutting planes, symmetrise, round down to denominator Q, add top-up eta.
  check    FILE.json [--maxdepth D]   exact th=0 check + exact branch and bound; writes FILE.check.txt
  test                                 sanity tests (formula vs seam_1d.V, closed form vs definition, adversarial)
"""
import sys, os, json, time, math, hashlib, argparse, random
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
TMAX = Fr(107, 256)          # tan(th_max / 2) with th_max = 45.37 deg >= 45 deg; checked: TMAX^2 + 2 TMAX - 1 >= 0


# ----------------------------------------------------------------------------------------------- exact primitives
def trig(t):
    """c = cos th, s = sin th for t = tan(th/2) (rational in, rational out)."""
    d = 1 + t * t
    return (1 - t * t) / d, 2 * t / d


def V_def(t, d):
    """Vertical chord of the closed unit square at angle th = 2 atan t in [0, 45 deg], at horizontal offset d.
    Straight from the statement (SEAM_1D.md section 1)."""
    d = abs(d)
    if t == 0:
        return Fr(1) if d <= Fr(1, 2) else Fr(0)
    c, s = trig(t); a = (c + s) / 2
    return max(Fr(0), min(1 / c, (a - d) / (s * c)))


def G_def(t, u):
    """Periodised chord: sum over images |j| <= 3 (|j| <= 1 would suffice: the support of V is |d| <= a <= 0.708)."""
    return sum(V_def(t, u + j) for j in range(-3, 4))


def rdist(u):
    """distance from u to the nearest integer, in [0, 1/2]."""
    f = u - math.floor(u)
    return min(f, 1 - f)


def G_closed(t, r):
    """Lemma 2 (SEAM_1D.md): for 0 < th <= 45 deg and r = dist(u, Z) in [0, 1/2]:
    G_th(u) = max(P, min(1/c, (a - r)/(s c))),  P = (c + s - 1)/(s c) = 1/c - s/(c (1 + c))."""
    c, s = trig(t); a = (c + s) / 2
    P = (c + s - 1) / (s * c)
    return max(P, min(1 / c, (a - r) / (s * c)))


class Measure:
    """atoms: dict k -> integer W (weight W/Q at x = k/N); eta: Fraction."""
    def __init__(self, N, Q, atoms, eta):
        self.N, self.Q, self.eta = N, Q, Fr(eta)
        self.atoms = {k % N: W for k, W in atoms.items() if W != 0}
        assert all(isinstance(W, int) and W > 0 for W in self.atoms.values()), 'weights must be positive ints'
        assert self.eta >= 0

    def weight(self, k): return Fr(self.atoms.get(k % self.N, 0), self.Q)
    def total(self): return Fr(sum(self.atoms.values()), self.Q) + self.eta
    def g(self): return self.weight(0)
    def e(self): return self.total() - 1
    def symmetric(self): return all(self.atoms.get((-k) % self.N, 0) == W for k, W in self.atoms.items())

    def F(self, t, cx, closed=True):
        """exact mass of the closed unit square at angle 2 atan t, centre x-coordinate cx."""
        tot = self.eta
        for k, W in self.atoms.items():
            u = Fr(k, self.N) - cx
            if closed and t > 0:
                tot += Fr(W, self.Q) * G_closed(t, rdist(u))
            else:
                tot += Fr(W, self.Q) * G_def(t, u)
        return tot

    # ---- certificate I/O
    def to_json(self, extra=None):
        d = {'format': 'seam1d-v1', 'N': self.N, 'Q': self.Q,
             'atoms': [[k, self.atoms[k]] for k in sorted(self.atoms)],
             'atoms_note': 'pair [k, W] means an atom of weight W/Q at x = k/N',
             'eta': str(self.eta), 'g': str(self.g()), 'e': str(self.e()),
             'g_float': float(self.g()), 'e_float': float(self.e()), 'kappa_float': float(self.e() / self.g()),
             'eta_float': float(self.eta)}
        if extra: d.update(extra)
        return d

    @staticmethod
    def from_json(d):
        m = Measure(d['N'], d['Q'], {k: W for k, W in d['atoms']}, Fr(d['eta']))
        assert Fr(d['g']) == m.g() and Fr(d['e']) == m.e(), 'stored g/e inconsistent with atoms'
        return m


def write_cert(m, path, extra=None):
    s = json.dumps(m.to_json(extra), indent=1)
    with open(path, 'w') as f: f.write(s + '\n')
    h = hashlib.sha256((s + '\n').encode()).hexdigest()
    with open(path + '.sha256', 'w') as f: f.write(f'{h}  {os.path.basename(path)}\n')
    return h


def read_cert(path):
    raw = open(path, 'rb').read()
    return Measure.from_json(json.loads(raw)), hashlib.sha256(raw).hexdigest()


# ------------------------------------------------------------------------------------------------- th = 0, exact
def check_theta0(m):
    """F(0, .) is piecewise constant, upper semicontinuous (closed squares).  Its breakpoints are cx = x_k +- 1/2.
    Evaluate F exactly at every breakpoint and at the midpoint of every gap between consecutive breakpoints."""
    bps = sorted({(Fr(k, m.N) + h) % 1 for k in m.atoms for h in (Fr(1, 2), Fr(-1, 2))} | {Fr(0)})
    pts = list(bps) + [(bps[i] + (bps[i + 1] if i + 1 < len(bps) else bps[0] + 1)) / 2 for i in range(len(bps))]
    vals = [(m.F(Fr(0), p % 1, closed=False), p % 1) for p in pts]
    return min(vals), len(pts)


# --------------------------------------------------------------------------------------------- branch and bound
class BnB:
    """Boxes [tlo, thi] x [Clo/L, Chi/L] with L = N 2^K (cx on a dyadic refinement of the atom grid, ints)."""
    def __init__(self, m, K=64):
        self.m = m; self.K = K; self.L = m.N << K; self.half = self.L // 2
        self.at = [(k << K, W) for k, W in sorted(m.atoms.items())]
        self.cache = {}
        self.Q = m.Q

    def tq(self, tlo, thi):
        """Lemma 3 quantities for the th-interval T = [2 atan tlo, 2 atan thi], valid for th in T cap [0, 45 deg] (the
        part above 45 deg, t > sqrt2 - 1, is redundant by th -> 90 - th and is not claimed):
        A = 1/c(tlo) <= 1/c,  aL <= a,  scU >= s c,  PL <= P;  R1 = floor(L (aL - scU A)), R2 = floor(L (aL - PL scU))."""
        key = (tlo, thi)
        q = self.cache.get(key)
        if q is None:
            clo, slo = trig(tlo); chi, shi = trig(thi)
            A = 1 / clo
            aL = (clo + slo) / 2                     # a increasing on [0, 45 deg]: a >= a(tlo) on T cap [0, 45]
            scU = shi * chi if thi * thi + 2 * thi - 1 <= 0 else Fr(1, 2)
            PL = A - shi / (chi * (1 + chi))
            R1 = math.floor((aL - scU * A) * self.L); R2 = math.floor((aL - PL * scU) * self.L)
            q = (A, aL, scU, PL, R1, R2)
            if len(self.cache) > 200000: self.cache.clear()
            self.cache[key] = q
        return q

    def lb(self, tlo, thi, Clo, Chi):
        """exact lower bound of F over the box (Lemma 3 + Lemma 2 per atom, summed)."""
        A, aL, scU, PL, R1, R2 = self.tq(tlo, thi)
        L, h, w = self.L, self.half, Chi - Clo
        Sp = Sr = SrR = SP = 0
        for p, W in self.at:
            lo = (p - Chi) % L; hi = lo + w             # u ranges over [lo, hi] / L  (hi < lo + L)
            if lo <= h <= hi or lo <= h + L <= hi:
                R = h
            else:
                R = max(min(lo, L - lo), min(hi % L, L - hi % L))
            if R <= R1: Sp += W
            elif R <= R2: Sr += W; SrR += W * R
            else: SP += W
        return self.m.eta + (A * Sp + (aL * Sr - Fr(SrR, L)) / scU + PL * SP) / self.Q

    def run(self, cx_max=None, maxdepth=90, log=None, progress=0, sample=0):
        m = self.m; self.samples = []
        if cx_max is None: cx_max = Fr(1, 2) if m.symmetric() else Fr(1)
        Cmax = cx_max * self.L; assert Cmax.denominator == 1; Cmax = int(Cmax)
        stack = [(Fr(0), TMAX, 0, Cmax, 0, 0, self.lb(Fr(0), TMAX, 0, Cmax))]
        leaves = 0; redundant = 0; maxd = 0; nodes = 0; minleaf = None; t0 = time.time(); dt_hist = {}
        while stack:
            tlo, thi, Clo, Chi, dt, dc, b = stack.pop(); nodes += 1
            if tlo * tlo + 2 * tlo - 1 >= 0:         # tlo >= sqrt2 - 1: the whole box has th >= 45 deg, redundant
                redundant += 1; continue
            if b >= 1:
                leaves += 1; maxd = max(maxd, dt + dc); dt_hist[dt] = dt_hist.get(dt, 0) + 1
                if minleaf is None or b < minleaf: minleaf = b
                if sample and leaves % sample == 0: self.samples.append((tlo, thi, Clo, Chi, b))
                if progress and leaves % progress == 0:
                    print(f'  leaves {leaves} nodes {nodes} stack {len(stack)} {time.time()-t0:.0f}s', flush=True)
                continue
            # failing box: is the centre an actual counterexample?
            tm, Cm = (tlo + thi) / 2, (Clo + Chi) // 2
            if dt + dc >= 24 or dt + dc >= maxdepth:
                tw = tm if tm * tm + 2 * tm - 1 <= 0 else tlo      # keep witnesses at th <= 45 deg (formula range)
                Fc = m.F(tw, Fr(Cm, self.L))
                if Fc < 1:
                    return dict(ok=False, reason='counterexample', t=tw, cx=Fr(Cm, self.L), F=Fc,
                                box=(tlo, thi, Fr(Clo, self.L), Fr(Chi, self.L)), nodes=nodes, leaves=leaves)
            if dt + dc >= maxdepth or Chi - Clo < 2:
                pts = [(m.F(tt, Fr(CC, self.L)) if tt > 0 else m.F(tt, Fr(CC, self.L), closed=False), tt, Fr(CC, self.L))
                       for tt in (tlo, tm, thi) for CC in (Clo, Cm, Chi) if tt * tt + 2 * tt - 1 <= 0]
                Fp, tp, cp = min(pts)
                return dict(ok=False, reason='counterexample (depth)' if Fp < 1 else 'depth', box=(tlo, thi, Fr(Clo, self.L), Fr(Chi, self.L)), lb=b,
                            t=tp, cx=cp, F=Fp, nodes=nodes, leaves=leaves)
            # split the coordinate whose width costs more: compare the bound with th (resp. cx) collapsed to the
            # midpoint.  Heuristic only; soundness does not depend on it.
            loss_t = self.lb(tm, tm, Clo, Chi) - b
            loss_c = self.lb(tlo, thi, Cm, Cm) - b
            if loss_t >= loss_c:
                kids = [(tlo, tm, Clo, Chi, dt + 1, dc), (tm, thi, Clo, Chi, dt + 1, dc)]
            else:
                kids = [(tlo, thi, Clo, Cm, dt, dc + 1), (tlo, thi, Cm, Chi, dt, dc + 1)]
            kids = [k + (self.lb(*k[:4]),) for k in kids]
            kids.sort(key=lambda k: -k[6])          # the worse child is popped first (finds counterexamples fast)
            stack += kids
        return dict(ok=True, leaves=leaves, redundant=redundant, nodes=nodes, maxdepth=maxd, minleaf=minleaf, cx_max=cx_max,
                    secs=time.time() - t0, dt_hist=dict(sorted(dt_hist.items())))


def certify(m, maxdepth=90, progress=0):
    out = {}
    v0, npts = check_theta0(m)
    out['theta0'] = dict(min=v0[0], at=v0[1], points=npts, ok=v0[0] >= 1)
    assert TMAX * TMAX + 2 * TMAX - 1 >= 0                      # TMAX >= tan(22.5 deg) = sqrt2 - 1
    out['sym'] = m.symmetric()
    out['bnb'] = BnB(m).run(maxdepth=maxdepth, progress=progress)
    out['ok'] = out['theta0']['ok'] and out['bnb']['ok'] and m.e() >= 0
    return out


# ------------------------------------------------------------------------------------------- float tools (prepare)
def fmin_theta(X, W, eta, th):
    """float: min over cx of F at angle th > 0.  F(th, .) is continuous piecewise linear in cx with kinks only where an
    atom's r crosses b or 1 - a (Lemma 2), so the min is at one of those cx.  Returns (min, cx)."""
    import numpy as np
    c, s = np.cos(th), np.sin(th); a = (c + s) / 2; b = (c - s) / 2
    P = (c + s - 1) / (s * c)
    crit = (X[:, None] + np.array([b, -b, 1 - a, a - 1])[None, :]).ravel() % 1
    u = X[None, :] - crit[:, None]; f = u - np.floor(u); r = np.minimum(f, 1 - f)
    G = np.maximum(P, np.minimum(1 / c, (a - r) / (s * c)))
    tot = eta + G @ W
    i = tot.argmin()
    return tot[i], crit[i]


def fscan(X, W, eta, dth=0.005, refine=True):
    """float: min of F over th in (0, 45] (grid + local refinement) and cx (exact breakpoints).  Returns list of
    (F, th_deg, cx) local minima sorted."""
    import numpy as np
    ths = np.radians(np.r_[np.arange(dth, 45.0 + 1e-9, dth), [1e-5, 1e-4, 1e-3, 0.01, 0.05]])
    ths.sort()
    vals = np.array([fmin_theta(X, W, eta, t)[0] for t in ths])
    loc = [i for i in range(len(ths)) if (i == 0 or vals[i] <= vals[i - 1]) and (i == len(ths) - 1 or vals[i] <= vals[i + 1])]
    res = []
    for i in loc:
        lo, hi = ths[max(i - 1, 0)], ths[min(i + 1, len(ths) - 1)]
        best = (vals[i], ths[i])
        if refine:
            for _ in range(3):
                g = np.linspace(lo, hi, 41)
                v = [fmin_theta(X, W, eta, t)[0] for t in g]
                j = int(np.argmin(v)); best = min(best, (v[j], g[j]))
                lo, hi = g[max(j - 1, 0)], g[min(j + 1, 40)]
        f, cx = fmin_theta(X, W, eta, best[1])
        res.append((f, np.degrees(best[1]), cx))
    res.sort()
    return res


def row_for(N, th, cx):
    import numpy as np
    xs = np.arange(N) / N
    u = xs - cx; f = u - np.floor(u); r = np.minimum(f, 1 - f)
    c, s = np.cos(th), np.sin(th); a = (c + s) / 2
    if s < 1e-12: return (r <= 0.5).astype(float)
    return np.maximum((c + s - 1) / (s * c), np.minimum(1 / c, (a - r) / (s * c)))


def prepare(args):
    import numpy as np
    sys.path.insert(0, HERE); import seam_1d
    if args.npy:
        x = np.load(args.npy); N = len(x); e_lp = float(x.sum() - 1)
    else:
        N = args.N; e_lp = args.e
    xs = np.arange(N) / N
    extra = []
    A = None
    margin = args.margin if args.margin is not None else 0.02 * e_lp
    best = None
    for rnd in range(args.cuts + 1):
        if rnd > 0 or not args.npy:
            if A is None:
                t0 = time.time(); A = seam_1d.build_rows(N); print(f'  rows {A.shape} {time.time()-t0:.0f}s', flush=True)
            t0 = time.time()
            try:
                x = seam_1d.solve(e_lp, A, np.array(extra) if extra else None)
            except RuntimeError as err:
                print('  ', err, '-- stopping cut rounds'); break
            print(f'  LP round {rnd}: g = {x[0]:.6f}  extra rows {len(extra)}  {time.time()-t0:.0f}s', flush=True)
        xsym = np.maximum((x + x[(-np.arange(N)) % N]) / 2, 0)    # mirror average: min F can only go up
        sup = np.nonzero(xsym > 1e-13)[0]
        mins = fscan(xs[sup], xsym[sup], 0.0)
        kp = (e_lp + max(0, 1 - mins[0][0]) + margin) / xsym[0]
        print(f'  float min F = {mins[0][0]:.7f} at th = {mins[0][1]:.4f} deg, cx = {mins[0][2]:.6f}; '
              f'{sum(1 for f in mins if f[0] < 1 - 1e-7)} local minima below 1; predicted kappa {kp:.5f}', flush=True)
        if best is None or kp < best[0]: best = (kp, xsym.copy(), rnd)
        if rnd == args.cuts or mins[0][0] >= 1 - 1e-7: break
        poses = [(np.radians(thd), cx) for f, thd, cx in mins if f < 1 - 1e-9]
        for th in np.radians(np.arange(0.05, 45.0001, 0.05)):            # dense cuts: argmin cx at each th
            f, cx = fmin_theta(xs[sup], xsym[sup], 0.0, th)
            if f < 1 - 1e-9: poses.append((th, cx))
        for th, cx in poses:
            for cc in (cx, -cx):                                       # mirror pair keeps the LP symmetric
                extra.append(row_for(N, th, cc % 1))
    kp, xsym, rnd = best
    print(f'  using round {rnd} (predicted kappa {kp:.5f})')
    if args.save_npy: np.save(args.save_npy, xsym)
    # exact rounding: W_k = floor(Q x_k) for k <= N/2, mirrored
    Q = args.Q
    Wd = {}
    for k in range(N // 2 + 1):
        W = int(math.floor(xsym[k] * Q))
        if W > 0: Wd[k] = W; Wd[(-k) % N] = W
    m0 = Measure(N, Q, Wd, 0)
    Xr = np.array([k / N for k in sorted(Wd)]); Wr = np.array([Wd[k] / Q for k in sorted(Wd)])
    mins = fscan(Xr, Wr, 0.0)
    deficit = max(0.0, 1 - mins[0][0])
    eta = Fr(math.ceil((deficit + margin) * 10 ** 9), 10 ** 9)
    m = Measure(N, Q, Wd, eta)
    print(f'  rounded: float min F (eta=0) = {mins[0][0]:.8f} at th {mins[0][1]:.4f}; eta = {float(eta):.3e} '
          f'(deficit {deficit:.2e} + margin {margin:.2e});  g = {float(m.g()):.6f}  e = {float(m.e()):.6f}  '
          f'kappa = {float(m.e()/m.g()):.5f}  eta/e = {float(eta/m.e()):.3f}', flush=True)
    h = write_cert(m, args.out, dict(source=args.npy or f'seam_1d LP N={N} e={e_lp}', e_lp=e_lp, cut_rounds=args.cuts, round_used=rnd,
                                     extra_rows=len(extra), float_min_F_before_topup=float(mins[0][0]),
                                     margin=margin))
    print('  wrote', args.out, 'sha256', h)


# ---------------------------------------------------------------------------------------------------------- tests
def run_tests(args):
    import numpy as np
    sys.path.insert(0, HERE); import seam_1d
    rng = random.Random(1)
    nfail = 0
    # 1. TMAX covers 45 deg; V depends on th only via |th| mod 90 and th -> 90 - th (checked numerically on the
    #    definition from the square's geometry: chord of a polygon, independent of the formula)
    assert TMAX * TMAX + 2 * TMAX - 1 >= 0
    def chord_geom(th, d):     # vertical chord of the rotated closed unit square by polygon clipping (float)
        c, s = math.cos(th), math.sin(th)
        P = [(c * x - s * y, s * x + c * y) for x, y in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]
        ys = []
        for i in range(4):
            (x1, y1), (x2, y2) = P[i], P[(i + 1) % 4]
            if (x1 - d) * (x2 - d) <= 0 and x1 != x2:
                ys.append(y1 + (y2 - y1) * (d - x1) / (x2 - x1))
            elif x1 == x2 == d: ys += [y1, y2]
        return max(ys) - min(ys) if ys else 0.0
    worst = 0
    for _ in range(20000):
        th = rng.uniform(-math.pi, math.pi); d = rng.uniform(-0.8, 0.8)
        g = chord_geom(th, d)
        th0 = abs(th) % (math.pi / 2); th0 = min(th0, math.pi / 2 - th0)          # reduce to [0, 45 deg]
        f = float(V_def(Fr(math.tan(th0 / 2)), Fr(d)))
        worst = max(worst, abs(f - g))
    print(f'[1] V formula (reduced to [0,45]) vs polygon chord at random th in (-180,180): max |diff| = {worst:.2e}')
    nfail += worst > 1e-9
    # 2. closed form (Lemma 2) == definition, exactly, at random rational poses incl. near-singular ones
    bad = 0
    for i in range(3000):
        t = Fr(rng.randint(1, 10 ** 6), 10 ** 6) * TMAX if i % 3 else Fr(1, rng.randint(10, 10 ** 7))
        u = Fr(rng.randint(-10 ** 6, 10 ** 6), 10 ** 6)
        if i % 5 == 0:       # atoms near u = +-1/2 and near the kinks
            c, s = trig(t); u = rng.choice([Fr(1, 2), (c - s) / 2, 1 - (c + s) / 2, (c + s) / 2]) + Fr(rng.randint(-3, 3), 10 ** 9)
        if G_closed(t, rdist(u)) != G_def(t, u): bad += 1
    print(f'[2] G closed form vs periodised definition, exact, 3000 rational poses: {bad} mismatches')
    nfail += bad > 0
    # 3. exact F vs seam_1d.V (float) at random poses, for a certificate measure
    if args.cert:
        m, _ = read_cert(args.cert)
        X = np.array([k / m.N for k in m.atoms]); Wf = np.array([float(Fr(W, m.Q)) for W in m.atoms.values()])
        worst = 0
        for i in range(300):
            t = Fr(rng.randint(0, 10 ** 6), 10 ** 6) * TMAX; cx = Fr(rng.randint(0, 10 ** 6), 10 ** 6)
            ex = m.F(t, cx)
            th = 2 * math.atan(float(t)); d = (X - float(cx) + 0.5) % 1 - 0.5
            fl = float(m.eta) + (Wf * (seam_1d.V(th, d) + seam_1d.V(th, d + 1) + seam_1d.V(th, d - 1))).sum()
            worst = max(worst, abs(float(ex) - fl))
        print(f'[3] exact F vs seam_1d.V at 300 random poses ({args.cert}): max |diff| = {worst:.2e}')
        nfail += worst > 1e-7
        # 4. B&B soundness spot check: random boxes, lb <= exact F at random interior rational points
        bb = BnB(m); viol = 0
        for i in range(300):
            depth = rng.randint(2, 30)
            tw = TMAX / 2 ** rng.randint(0, depth); tlo = tw * rng.randint(0, int(TMAX / tw) - 1) if tw < TMAX else Fr(0)
            thi = min(tlo + tw, TMAX)
            cw = bb.L >> (rng.randint(1, depth + 1)); Clo = cw * rng.randint(0, bb.L // cw - 1); Chi = Clo + cw
            b = bb.lb(tlo, thi, Clo, Chi)
            for _ in range(5):
                t = tlo + (thi - tlo) * Fr(rng.randint(0, 1000), 1000); cx = Fr(Clo, bb.L) + Fr(cw, bb.L) * Fr(rng.randint(0, 1000), 1000)
                if 0 < t and t * t + 2 * t - 1 <= 0 and m.F(t, cx, closed=False) < b: viol += 1   # claim: th <= 45
        print(f'[4] box lower bound <= F (by definition) at 1500 random points of 300 random boxes: {viol} violations')
        nfail += viol > 0
        # 4b. the actual certificate leaves: every 50th leaf, F (definition) at 4 random points >= its bound >= 1
        bb = BnB(m); r = bb.run(sample=50) if not args.no4b else {'ok': None}; viol = 0
        for tlo, thi, Clo, Chi, b in getattr(bb, 'samples', []):
            for _ in range(4):
                t = tlo + (thi - tlo) * Fr(rng.randint(1, 1000), 1000); cx = Fr(Clo, bb.L) + Fr(Chi - Clo, bb.L) * Fr(rng.randint(0, 1000), 1000)
                if t * t + 2 * t - 1 <= 0 and not (m.F(t, cx, closed=False) >= b >= 1): viol += 1
        print(f'[4b] {len(getattr(bb, "samples", []))} sampled certificate leaves x 4 random points: {viol} violations (B&B ok={r["ok"]})')
        nfail += viol > 0
        # 5. adversarial: remove the top-up; perturb one atom.  Both must be rejected near the float minimum.
        Xs = np.array([k / m.N for k in sorted(m.atoms)]); Ws = np.array([m.atoms[k] / m.Q for k in sorted(m.atoms)])
        mins = fscan(Xs, Ws, 0.0)
        print(f'    float min of F without eta: {mins[0][0]:.8f} at th = {mins[0][1]:.4f} deg, cx = {mins[0][2]:.6f}')
        # 5a: top-up removed and every atom scaled by (1 - 1e-4), so min F ~ 1 - 1e-4.  (Removing eta alone leaves
        #     F within ~1e-6 of 1 on whole 2D patches of LP-tight poses; rejecting that needs boxes of size ~1e-7,
        #     as hard as certifying at margin 1e-7: measured > 1 h, not run.)
        bad1 = Measure(m.N, m.Q, {kk: W - W // 10000 for kk, W in m.atoms.items()}, 0)
        t0 = time.time(); r = BnB(bad1).run(maxdepth=70)
        print(f'[5a] eta = 0, atoms x (1 - 1e-4): ok={r["ok"]} reason={r.get("reason")} at th = '
              f'{math.degrees(2*math.atan(float(r.get("t", 0)))):.4f} deg, cx = {float(r.get("cx", 0)):.6f}, '
              f'F - 1 = {float(r.get("F", 1) - 1):.3e}  (nodes {r["nodes"]}, {time.time()-t0:.0f}s)')
        nfail += r['ok']
        # the atom contributing most at the float-minimum pose, lowered so that F there becomes ~1 - 1e-4
        th_s, cx_s = math.radians(mins[0][1]), mins[0][2]
        Fs = float(m.eta) + mins[0][0]
        contrib = {kk: (W / m.Q) * float(G_closed(Fr(math.tan(th_s / 2)), Fr(rdist(kk / m.N - cx_s)))) for kk, W in m.atoms.items()}
        k = max((kk for kk in contrib if kk), key=lambda kk: contrib[kk])
        Gk = contrib[k] / (m.atoms[k] / m.Q)
        dec = int(((Fs - 1 + 1e-4) / Gk) * m.Q) + 1           # F at that pose drops to ~1 - 1e-4
        at2 = dict(m.atoms); at2[k] -= dec; at2[(-k) % m.N] = at2[k]
        bad2 = Measure(m.N, m.Q, at2, m.eta)
        tq_, cq_ = Fr(math.tan(th_s / 2)).limit_denominator(10 ** 12), Fr(cx_s).limit_denominator(10 ** 12)
        print(f'    perturbed measure at the float-min pose (th = {mins[0][1]:.4f}, cx = {cx_s:.6f}): exact F - 1 = '
              f'{float(bad2.F(tq_, cq_) - 1):.3e}')
        t0 = time.time(); r = BnB(bad2).run(maxdepth=70)
        print(f'[5b] atom k={k} (x={k/m.N}) lowered by {dec/m.Q:.2e} (mirror too; float min pose th={mins[0][1]:.4f}, cx={cx_s:.6f}): '
              f'ok={r["ok"]} reason={r.get("reason")} th = {math.degrees(2*math.atan(float(r.get("t", 0)))):.4f} deg, '
              f'cx = {float(r.get("cx", 0)):.6f}, F - 1 = {float(r.get("F", 1) - 1):.3e}  ({time.time()-t0:.0f}s)')
        nfail += r['ok']
        # 5c: one-sided perturbation of the same atom: symmetry check must fail, cx in [0,1] is used
        at3 = dict(m.atoms); at3[k] -= dec
        bad3 = Measure(m.N, m.Q, at3, m.eta)
        print(f'[5c] one-sided perturbation: symmetric() = {bad3.symmetric()} (must be False)')
        nfail += bad3.symmetric()
        r = BnB(bad3).run(maxdepth=70)
        print(f'     ok={r["ok"]} reason={r.get("reason")} th = {math.degrees(2*math.atan(float(r.get("t", 0)))):.4f} deg, '
              f'cx = {float(r.get("cx", 0)):.6f}, F - 1 = {float(r.get("F", 1) - 1):.3e}')
        nfail += r['ok']
        # 5d: theta = 0 check must reject a measure with total < 1 (e.g. scale everything down)
        bad4 = Measure(m.N, m.Q, {k: W // 2 for k, W in m.atoms.items()}, 0)
        v0, _ = check_theta0(bad4)
        print(f'[5d] half measure: theta=0 min F = {float(v0[0]):.4f} (must be < 1)')
        nfail += v0[0] >= 1
    print('TESTS', 'FAILED' if nfail else 'PASSED', nfail)
    return nfail


# ----------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('prepare'); p.add_argument('--npy'); p.add_argument('--N', type=int, default=200)
    p.add_argument('--e', type=float); p.add_argument('--cuts', type=int, default=0)
    p.add_argument('--margin', type=float); p.add_argument('--Q', type=int, default=10 ** 9); p.add_argument('--out', required=True)
    p.add_argument('--save-npy', dest='save_npy')
    p = sp.add_parser('check'); p.add_argument('cert'); p.add_argument('--maxdepth', type=int, default=90)
    p.add_argument('--progress', type=int, default=0)
    p = sp.add_parser('test'); p.add_argument('--cert'); p.add_argument('--no4b', action='store_true')
    a = ap.parse_args()
    if a.cmd == 'prepare':
        prepare(a)
    elif a.cmd == 'check':
        m, h = read_cert(a.cert)
        t0 = time.process_time()
        r = certify(m, a.maxdepth, a.progress)
        cpu = time.process_time() - t0
        b = r['bnb']
        lines = [f'certificate {a.cert}  sha256 {h}',
                 f'N = {m.N}  atoms = {len(m.atoms)}  Q = {m.Q}  symmetric = {r["sym"]}',
                 f'g = {m.g()}  ({float(m.g()):.8f})', f'e = {m.e()}  ({float(m.e()):.8f})',
                 f'eta = {m.eta}  ({float(m.eta):.3e}, eta/e = {float(m.eta / m.e()):.4f})',
                 f'kappa = e/g = {float(m.e() / m.g()):.6f}',
                 f'theta = 0: min F = {r["theta0"]["min"]} over {r["theta0"]["points"]} breakpoints+midpoints  ok={r["theta0"]["ok"]}',
                 f'B&B over t = tan(th/2) in [0, {TMAX}] (th <= {math.degrees(2*math.atan(float(TMAX))):.3f} deg), '
                 f'cx in [0, {b.get("cx_max")}]: ok={b["ok"]}']
        if b['ok']:
            lines += [f'  leaves {b["leaves"]} (+ {b["redundant"]} boxes with th >= 45 deg skipped)  nodes {b["nodes"]}  max depth {b["maxdepth"]}  '
                      f'min leaf bound {float(b["minleaf"]):.9f}', f'  leaves by th-depth {b["dt_hist"]}']
        else:
            lines += [f'  FAILED: {b}']
        lines += [f'CPU {cpu:.1f}s', f'RESULT {"CERTIFIED" if r["ok"] else "NOT CERTIFIED"}']
        txt = '\n'.join(lines)
        print(txt)
        with open(a.cert.replace('.json', '') + '.check.txt', 'w') as f: f.write(txt + '\n')
        sys.exit(0 if r['ok'] else 1)
    elif a.cmd == 'test':
        sys.exit(1 if run_tests(a) else 0)


if __name__ == '__main__':
    main()
