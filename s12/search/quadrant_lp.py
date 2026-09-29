#!/usr/bin/env python3
"""quadrant_lp.py -- HEURISTIC quadrant LP for fixed-profile families (task quadrant-lp, 2026-09-28).
search/QUADRANT.md is the write-up.  Nothing here is a proof: float LP over generated rows, float oracle.

Family (quadrant Q = [0,inf)^2, corner at the origin, diagonal symmetry (x,y) <-> (y,x)):
  * corner module: free measure on C = [0,R]^2 (points, axis-parallel segment pieces, uniform cells), R integer;
  * wall band [0,w]: a 1-periodic profile (phase cell [0,1) x [0,w], phase 0 at integer x) repeated along the
    bottom wall at x = j + phi for integers j >= R, and along the left wall by the diagonal reflection.
    The profile is mirror symmetric, phi -> -phi (mod 1), so that the D4 box is consistent (right corner = mirror);
  * Lebesgue measure on L = {x > w, y > w} \\ [w,R]^2   (lebwall mode: w_leb = 0, i.e. L = Q \\ C).
Objective (QUADRANT.md sec 2):   D = R^2 - M_C - m_v,   M_C = total corner-module mass,
  m_v = profile mass on the phase-0 cross-section {0} x [0,w].   Box saving k^2 - mu([0,k]^2) = 4 D.
sigma = 0:  profile mass per period cell = w  (equality row).
Rows: closed unit squares in Q, pose (cx, cy, th), th in [0, pi/2), restricted to cx >= cy (diagonal symmetry),
  cy <= R + 0.75, cx <= R + 1.75 (the far band is periodic: one period beyond the corner's reach suffices).
  Poses are kept WALL = 1e-7 inside Q (closed-semantics: wall-touching squares are limits of interior ones).

Modes
  lp      --R R --w W [--hc .1 --hca .25 --hp .1 --hpa .25] [--rounds N] [--tag T]   row generation
  lebwall --R R  [...]          band = Lebesgue, corner free (sanity: D -> 0)
  naga                          Nagamochi's structure (R=2, w=1) fixed: oracle min + D
  nagaprof --R R                Nagamochi's profile fixed, corner free
  band --w W                    band only (no corner): max m_v s.t. sigma = 0 and validity in the half-plane
  (checks: quadrant_tools.py check / box)
Outputs runs/quad_<tag>/ : log, per-round json, solution npz.
"""
import sys, os, math, time, json, argparse
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
import numpy as np
import scipy.sparse as sp
import multiprocessing as mp

TOL = 1e-9          # closed containment (safe direction)
WALL = 1e-7
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, '..', 'runs')
BIG = 1e6


# ============================================================================ geometry
def corners(P):
    """P (n,3) poses -> (n,4,2) square vertices (counter-clockwise)."""
    c = np.cos(P[:, 2]); s = np.sin(P[:, 2])
    ex = np.stack([c, s], 1) * 0.5; ey = np.stack([-s, c], 1) * 0.5
    ctr = P[:, :2]
    return np.stack([ctr - ex - ey, ctr + ex - ey, ctr + ex + ey, ctr - ex + ey], 1)


def _clip(V, nx, ny, c):
    """Sutherland-Hodgman, keep nx*x + ny*y <= c.  V (n,M,2) closed polygons (duplicates allowed); nx, ny, c (n,)."""
    n, M, _ = V.shape
    s = V[:, :, 0] * nx[:, None] + V[:, :, 1] * ny[:, None] - c[:, None]
    Vn = np.roll(V, -1, axis=1); sn = np.roll(s, -1, axis=1)
    ins = s <= 0; insn = sn <= 0
    cross = ins != insn
    den = np.where(cross, s - sn, 1.0)
    t = np.where(cross, s / den, 0.0)
    I = V + t[:, :, None] * (Vn - V)
    cand = np.stack([V, I], 2).reshape(n, 2 * M, 2)
    mask = np.stack([ins, cross], 2).reshape(n, 2 * M)
    order = np.argsort(~mask, axis=1, kind='stable')
    cand = np.take_along_axis(cand, order[:, :, None], 1)
    cnt = mask.sum(1)
    K = M + 1
    cand = cand[:, :K]
    idx = np.minimum(np.arange(K)[None, :], np.maximum(cnt - 1, 0)[:, None])
    out = np.take_along_axis(cand, idx[:, :, None], 1)
    out[cnt == 0] = 0.0
    return out


def area_rect(V, x0, x1, y0, y1):
    """area of polygon V (n,M,2) intersected with [x0,x1]x[y0,y1] (arrays (n,))."""
    n = len(V)
    one = np.ones(n); zero = np.zeros(n)
    W = _clip(V, one, zero, x1)
    W = _clip(W, -one, zero, -x0)
    W = _clip(W, zero, one, y1)
    W = _clip(W, zero, -one, -y0)
    x = W[:, :, 0]; y = W[:, :, 1]
    return 0.5 * np.abs((x * np.roll(y, -1, 1) - np.roll(x, -1, 1) * y).sum(1))


def _interval(k, A, B):
    pos = k > 1e-14; neg = k < -1e-14; zer = ~(pos | neg)
    ks = np.where(zer, 1.0, k)
    lo = np.where(pos, A / ks, np.where(neg, B / ks, np.where((A <= 0) & (B >= 0), -np.inf, np.inf)))
    hi = np.where(pos, B / ks, np.where(neg, A / ks, np.where((A <= 0) & (B >= 0), np.inf, -np.inf)))
    return lo, hi


def chord(orient, a, cx, cy, th, tol=TOL):
    """(line_cover.chord) parameter interval of the line x=a (orient 0, param y) / y=a (orient 1, param x)
    inside the closed unit square (cx,cy,th).  Empty: lo > hi."""
    h = 0.5 + tol
    c = np.cos(th); s = np.sin(th)
    if orient == 0:
        d = a - cx; base = cy
        l1, h1 = _interval(s, -h - c * d, h - c * d)
        l2, h2 = _interval(c, -h + s * d, h + s * d)
    else:
        d = a - cy; base = cx
        l1, h1 = _interval(c, -h - s * d, h - s * d)
        l2, h2 = _interval(s, -h + c * d, h + c * d)
    return base + np.maximum(l1, l2), base + np.minimum(h1, h2)


def halfwidth(th):
    return (np.abs(np.cos(th)) + np.abs(np.sin(th))) / 2


# ============================================================================ model
class Model:
    """elements: pts (x,y,var,fixedmass), hs (y,x0,x1,var,fm), vs (x,y0,y1,var,fm), cells (x0,x1,y0,y1,var,fm).
    var = -1 for fixed elements.  Each element carries the full mass of its variable (mass per element)."""

    def __init__(self, R, w, wleb=None):
        self.R = int(R); self.w = float(w); self.wleb = self.w if wleb is None else float(wleb)
        self.E = {'pts': [], 'hs': [], 'vs': [], 'cells': []}
        self.K = {'pts': [], 'hs': [], 'vs': [], 'cells': []}; self.kind = 0; self.base = []     # 0 corner, 1 profile copy
        self.nvar = 0; self.cobj = []; self.csig = []; self.label = []; self.const_obj = 0.0; self.const_sig = 0.0
        self.J = 4; self.G = self.R; self.band = False; self.jlo = self.R

    def _add(self, k, tup):
        self.E[k].append(tup); self.K[k].append(self.kind)

    def newvar(self, obj, sig, label):
        self.cobj.append(obj); self.csig.append(sig); self.label.append(label); self.nvar += 1
        return self.nvar - 1

    # -- corner atoms (orbit under diagonal reflection) ---------------------------------------------------------
    def corner_point(self, x, y, var=-1, fm=0.0):
        self.kind = 0
        orb = sorted(set([(x, y), (y, x)]))
        if var == -2: var = self.newvar(len(orb), 0, ('Cp', x, y))
        for (a, b) in orb: self._add('pts', (a, b, var, fm))
        return var

    def corner_hseg(self, y, x0, x1, var=-1, fm=0.0):
        """horizontal piece y, [x0,x1] and its mirror (vertical x=y, [x0,x1]); fm = mass per piece."""
        self.kind = 0
        if var == -2: var = self.newvar(2, 0, ('Cs', y, x0, x1))
        self._add('hs', (y, x0, x1, var, fm)); self._add('vs', (y, x0, x1, var, fm))
        return var

    def corner_cell(self, x0, x1, y0, y1, var=-1, fm=0.0):
        self.kind = 0
        orb = sorted(set([(x0, x1, y0, y1), (y0, y1, x0, x1)]))
        if var == -2: var = self.newvar(len(orb), 0, ('Cc', x0, x1, y0, y1))
        for o in orb: self._add('cells', o + (var, fm))
        return var

    # -- profile atoms (phase cell; mirror phi -> -phi mod 1; copies j = R..R+J on both walls) --------------------
    def js(self):
        return range(-3, 5) if self.band else range(self.R, self.R + self.J + 1)

    def _copies_pt(self, phi, y, var, fm):
        for j in self.js():
            self._add('pts', (j + phi, y, var, fm))
            if not self.band: self._add('pts', (y, j + phi, var, fm))

    def prof_point(self, phi, y, var=-1, fm=0.0):
        self.kind = 1
        orb = sorted(set([round(phi % 1.0, 12), round((-phi) % 1.0, 12)]))
        mv = 1 if orb == [0.0] else 0
        if var == -2: var = self.newvar(mv, len(orb), ('Pp', phi, y))
        else: self.const_obj += mv * fm; self.const_sig += len(orb) * fm
        for p in orb: self._copies_pt(p, y, var, fm); self.base.append(('p', p, p, y, y, var, fm))
        return var

    def prof_hseg(self, y, p0, p1, var=-1, fm=0.0):
        """horizontal piece at height y, phases [p0,p1] subset [0,1]; mirror [1-p1, 1-p0]."""
        self.kind = 1
        orb = sorted(set([(round(p0, 12), round(p1, 12)), (round(1 - p1, 12), round(1 - p0, 12))]))
        if var == -2: var = self.newvar(0, len(orb), ('Ph', y, p0, p1))
        else: self.const_sig += len(orb) * fm
        for (a, b) in orb:
            self.base.append(('h', a, b, y, y, var, fm))
            for j in self.js():
                self._add('hs', (y, j + a, j + b, var, fm))
                if not self.band: self._add('vs', (y, j + a, j + b, var, fm))
        return var

    def prof_vseg(self, phi, y0, y1, var=-1, fm=0.0):
        self.kind = 1
        orb = sorted(set([round(phi % 1.0, 12), round((-phi) % 1.0, 12)]))
        mv = 1 if orb == [0.0] else 0
        if var == -2: var = self.newvar(mv, len(orb), ('Pv', phi, y0, y1))
        else: self.const_obj += mv * fm; self.const_sig += len(orb) * fm
        for p in orb:
            self.base.append(('v', p, p, y0, y1, var, fm))
            for j in self.js():
                self._add('vs', (j + p, y0, y1, var, fm))
                if not self.band: self._add('hs', (j + p, y0, y1, var, fm))
        return var

    def prof_cell(self, p0, p1, y0, y1, var=-1, fm=0.0):
        self.kind = 1
        orb = sorted(set([(round(p0, 12), round(p1, 12)), (round(1 - p1, 12), round(1 - p0, 12))]))
        if var == -2: var = self.newvar(0, len(orb), ('Pc', p0, p1, y0, y1))
        else: self.const_sig += len(orb) * fm
        for (a, b) in orb:
            self.base.append(('c', a, b, y0, y1, var, fm))
            for j in self.js():
                self._add('cells', (j + a, j + b, y0, y1, var, fm))
                if not self.band: self._add('cells', (y0, y1, j + a, j + b, var, fm))
        return var

    def finalize(self):
        self.A = {}
        for k, v in self.E.items():
            n = {'pts': 4, 'hs': 5, 'vs': 5, 'cells': 6}[k]
            self.A[k] = np.array(v, float).reshape(-1, n)
        self.KA = {k: np.array(v, int) for k, v in self.K.items()}
        self.cobj = np.array(self.cobj, float); self.csig = np.array(self.csig, float)

    # -- coefficients --------------------------------------------------------------------------------------------
    def contrib(self, P, xval=None, keep=None):
        """P (n,3) poses.  Returns (Acoef csr (n, nvar), const (n,)) where const = fixed elements + Lebesgue.
        If xval is given, returns just the total captured mass (n,) with element mass x[var] (fixed: fm)."""
        n = len(P); rows, cols, vals = [], [], []
        const = self.lebesgue(P)
        tot = const.copy() if xval is not None else None

        def emit(pi, ei, frac, arr):
            var = arr[ei, -2].astype(int); fm = arr[ei, -1]
            if xval is not None:
                m = np.where(var >= 0, xval[np.maximum(var, 0)], fm)
                np.add.at(tot, pi, frac * m)
                return
            fx = var < 0
            if fx.any(): np.add.at(const, pi[fx], frac[fx] * fm[fx])
            nf = ~fx
            rows.append(pi[nf]); cols.append(var[nf]); vals.append(frac[nf])

        def sel(arr):
            if xval is None or len(arr) == 0: return arr
            var = arr[:, -2].astype(int)
            m = np.where(var >= 0, xval[np.maximum(var, 0)], arr[:, -1])
            return arr[m > 0]

        pts = sel(self.A['pts']); hs = sel(self.A['hs']); vs = sel(self.A['vs']); cl = sel(self.A['cells'])
        CH = 128
        ct = np.cos(P[:, 2]); st = np.sin(P[:, 2])
        for i in range(0, n, CH):
            B = P[i:i + CH]; c_ = ct[i:i + CH]; s_ = st[i:i + CH]; ib = np.arange(i, min(i + CH, n))
            if len(pts):
                dx = pts[None, :, 0] - B[:, 0:1]; dy = pts[None, :, 1] - B[:, 1:2]
                u = dx * c_[:, None] + dy * s_[:, None]; v = -dx * s_[:, None] + dy * c_[:, None]
                ins = np.maximum(np.abs(u), np.abs(v)) <= 0.5 + TOL
                a, b = np.nonzero(ins)
                emit(ib[a], b, np.ones(len(a)), pts)
            for orient, arr in ((1, hs), (0, vs)):
                if not len(arr): continue
                # near filter: line coordinate within 0.71 of centre
                cc = B[:, 1:2] if orient == 1 else B[:, 0:1]
                pc = B[:, 0:1] if orient == 1 else B[:, 1:2]
                near = (np.abs(arr[None, :, 0] - cc) <= 0.7072) & (arr[None, :, 1] <= pc + 0.7072) & (arr[None, :, 2] >= pc - 0.7072)
                a, b = np.nonzero(near)
                if not len(a): continue
                lo, hi = chord(orient, arr[b, 0], B[a, 0], B[a, 1], B[a, 2])
                L = np.minimum(hi, arr[b, 2]) - np.maximum(lo, arr[b, 1])
                ln = arr[b, 2] - arr[b, 1]
                fr = np.clip(L / ln, 0.0, 1.0)
                ok = fr > 1e-15
                emit(ib[a[ok]], b[ok], fr[ok], arr)
            if len(cl):
                near = ((cl[None, :, 0] <= B[:, 0:1] + 0.7072) & (cl[None, :, 1] >= B[:, 0:1] - 0.7072) &
                        (cl[None, :, 2] <= B[:, 1:2] + 0.7072) & (cl[None, :, 3] >= B[:, 1:2] - 0.7072))
                a, b = np.nonzero(near)
                if len(a):
                    V = corners(B[a])
                    ar = area_rect(V, cl[b, 0], cl[b, 1], cl[b, 2], cl[b, 3])
                    fr = ar / ((cl[b, 1] - cl[b, 0]) * (cl[b, 3] - cl[b, 2]))
                    ok = fr > 1e-15
                    emit(ib[a[ok]], b[ok], fr[ok], cl)
        if xval is not None: return tot
        if rows:
            r = np.concatenate(rows); c = np.concatenate(cols); v = np.concatenate(vals)
        else:
            r = c = np.zeros(0, int); v = np.zeros(0)
        A = sp.csr_matrix((v, (r, c)), shape=(n, self.nvar))
        return A, const

    def lebesgue(self, P):
        V = corners(P); n = len(P); wl = self.wleb; R = self.R
        if self.band:
            return area_rect(V, np.full(n, -BIG), np.full(n, BIG), np.full(n, wl), np.full(n, BIG))
        a1 = area_rect(V, np.full(n, wl), np.full(n, BIG), np.full(n, wl), np.full(n, BIG))
        a2 = area_rect(V, np.full(n, wl), np.full(n, float(R)), np.full(n, wl), np.full(n, float(R)))
        return a1 - a2

    def D(self, x):
        """D = R^2 - M_C - m_v (cobj carries both; band mode: -m_v)."""
        return self.R ** 2 - (self.cobj @ x + self.const_obj)


# ============================================================================ builders
def lattice(a, b, h):
    n = int(round((b - a) / h))
    return [round(a + i * h, 10) for i in range(n + 1)]


def add_free_corner(m, hc, hca):
    R = m.R
    g = lattice(0, R, hc)
    for x in g:
        for y in g:
            if x <= 0 or y <= 0 or y > x: continue            # wall points useless (closed semantics); orbit rep y<=x
            m.corner_point(x, y, var=-2)
    for y in g:
        if y <= 0: continue
        for i in range(len(g) - 1):
            m.corner_hseg(y, g[i], g[i + 1], var=-2)
    if hca:
        ga = lattice(0, R, hca)
        for i in range(len(ga) - 1):
            for j in range(len(ga) - 1):
                if ga[j] > ga[i]: continue
                m.corner_cell(ga[i], ga[i + 1], ga[j], ga[j + 1], var=-2)


def add_free_profile(m, hp, hpa):
    w = m.w
    phs = lattice(0, 1, hp)[:-1]
    ys = lattice(0, w, hp)
    for phi in phs:
        if (-phi) % 1.0 < phi - 1e-12: continue
        for y in ys:
            if y <= 0: continue
            m.prof_point(phi, y, var=-2)
    for y in ys:
        if y <= 0: continue
        for i, phi in enumerate(phs):
            p0, p1 = phi, round(phi + hp, 10)
            if 1 - p1 < p0 - 1e-12: continue
            m.prof_hseg(y, p0, p1, var=-2)
    for phi in phs:
        if (-phi) % 1.0 < phi - 1e-12: continue
        for i in range(len(ys) - 1):
            m.prof_vseg(phi, ys[i], ys[i + 1], var=-2)
    if hpa:
        pa = lattice(0, 1, hpa); ya = lattice(0, w, hpa)
        for i in range(len(pa) - 1):
            if 1 - pa[i + 1] < pa[i] - 1e-12: continue
            for j in range(len(ya) - 1):
                m.prof_cell(pa[i], pa[i + 1], ya[j], ya[j + 1], var=-2)


def add_naga_corner(m, fixed=True):
    """Nagamochi's corner in quadrant coordinates, R = 2: Lebesgue on [1,2]^2, segments y=1 / x=1 on [0.9,2]
    (0.5 per length), Q points (0.9,1), (1,0.9) (0.45 each)."""
    m.corner_cell(1.0, 2.0, 1.0, 2.0, fm=1.0)
    m.corner_hseg(1.0, 0.9, 2.0, fm=0.55)
    m.corner_point(1.0, 0.9, fm=0.45)


def add_naga_profile(m):
    """line y = 1 at 0.5 per length (per period: 0.5), point (phase 0, 0.9) mass 0.5.  Band [0,1], w = 1."""
    m.prof_hseg(1.0, 0.0, 0.5, fm=0.25)      # the orbit {[0,.5],[.5,1]} -> 2 pieces x 0.25 = 0.5 per period
    m.prof_point(0.0, 0.9, fm=0.5)


# ============================================================================ poses
def clipQ(P, R):
    if isinstance(R, tuple):             # band-only geometry ('band', w): periodic in x, half-plane y >= 0
        w = R[1]; hw = halfwidth(P[:, 2])
        P[:, 0] = np.mod(P[:, 0], 1.0)
        P[:, 1] = np.minimum(np.maximum(P[:, 1], hw + WALL), w + 0.75)
        return P
    hw = halfwidth(P[:, 2])
    P[:, 0] = np.maximum(P[:, 0], hw + WALL); P[:, 1] = np.maximum(P[:, 1], hw + WALL)
    P[:, 1] = np.minimum(P[:, 1], R + 0.75)
    P[:, 0] = np.minimum(P[:, 0], R + 1.75)
    sw = P[:, 1] > P[:, 0]
    if sw.any():   # diagonal reflection: (cx,cy,th) -> (cy,cx,pi/2-th)
        a = P[sw, 0].copy(); P[sw, 0] = P[sw, 1]; P[sw, 1] = a
        P[sw, 2] = (math.pi / 2 - P[sw, 2]) % (math.pi / 2)
    return P


def _ext(R):
    """(x extent, y extent) of the pose window"""
    if isinstance(R, tuple): return 1.0, R[1] + 0.75
    return R + 1.75, R + 0.75


def lattice_poses(R, pitch, degs):
    if isinstance(R, tuple):
        xs = np.arange(0, 1.0, pitch); ys = np.arange(0, R[1] + 0.75 + 1e-9, pitch); out = []
        for d in degs:
            th = math.radians(d); hw = halfwidth(np.array([th]))[0]
            X, Y = np.meshgrid(xs, ys[ys >= hw - 1e-12], indexing='ij')
            out.append(np.c_[X.ravel(), Y.ravel(), np.full(X.size, th)])
        return clipQ(np.concatenate(out), R)
    xs = np.arange(0, R + 1.75 + 1e-9, pitch); ys = np.arange(0, R + 0.75 + 1e-9, pitch)
    out = []
    for d in degs:
        th = math.radians(d); hw = halfwidth(np.array([th]))[0]
        X, Y = np.meshgrid(xs[xs >= hw - 1e-12], ys[ys >= hw - 1e-12], indexing='ij')
        X = X.ravel(); Y = Y.ravel(); k = Y <= X + 1e-12
        out.append(np.c_[X[k], Y[k], np.full(k.sum(), th)])
    P = np.concatenate(out)
    return clipQ(P, R)


def random_poses(R, n, rng):
    ex, ey = _ext(R)
    P = np.c_[rng.uniform(0, ex, n), rng.uniform(0, ey, n), rng.uniform(0, math.pi / 2, n)]
    r = rng.random(n)
    P[r < 0.25, 2] = rng.choice([0.0, 1e-6, 1e-4, 1e-2, math.pi / 4, math.pi / 2 - 1e-4], size=int((r < 0.25).sum()))
    P = clipQ(P, R)
    r = rng.random(n); hw = halfwidth(P[:, 2])
    P[r < 0.3, 1] = hw[r < 0.3] + WALL                    # resting on the bottom wall
    if not isinstance(R, tuple): k = r < 0.08; P[k, 0] = hw[k] + WALL                   # in the corner
    return clipQ(P, R)


def near_lattice_poses(R, h, rng=None, nrand=0, full=True):
    """axis-parallel (and infinitesimally tilted) squares whose edges sit just off the atom lattice h*Z (the dilated-grid
    / epsilon-gap poses of CORNER_DEFICIT.md sec 0-1): centres at lattice + 1/2 + (+-d, +-d), d in {1e-7, 1e-5}, angles
    0, 1e-6, pi/2 - 1e-6; plus (nrand) random centres within 2e-3 of lattice + 1/2 at angles within 2e-3 of 0."""
    ex, ey = _ext(R); band = isinstance(R, tuple)
    g = np.arange(0, ex + 1e-9, h)
    out = []
    for d in ((1e-7, 1e-5) if full else (1e-7,)):
        for sx in (-1, 1):
            for sy in (-1, 1):
                X, Y = np.meshgrid(g + 0.5 + sx * d, np.arange(0, ey + 1e-9, h) + 0.5 + sy * d, indexing='ij')
                X = X.ravel(); Y = Y.ravel(); k = (Y <= ey) & ((Y <= X + 1e-6) | band)
                for th in ((0.0, 1e-6, math.pi / 2 - 1e-6) if full else (0.0,)):
                    out.append(np.c_[X[k], Y[k], np.full(k.sum(), th)])
    if nrand and rng is not None:
        X = rng.choice(g, nrand) + 0.5 + rng.uniform(-2e-3, 2e-3, nrand)
        gy = np.arange(0, ey + 1e-9, h)
        Y = rng.choice(gy[gy <= ey - 0.45], nrand) + 0.5 + rng.uniform(-2e-3, 2e-3, nrand)
        th = np.mod(rng.uniform(-2e-3, 2e-3, nrand), math.pi / 2)
        out.append(np.c_[X, Y, th])
    return clipQ(np.concatenate(out), R)


def germ_poses(R, h, thetas=(1e-5, -1e-5, 1e-3, -1e-3)):
    """tile-germ blow-up family (line_cover.germ_poses, LINE_COVER.md sec 2): centres (g_i + 1/2 + th t_x,
    g_j + 1/2 + th t_y) at the tiny angle th, t_x, t_y on the h-grid of [-1/2, 1/2].  At th -> 0 the capture of
    piecewise-uniform line densities is piecewise linear in (t_x, t_y) with vertices on this grid."""
    ex, ey = _ext(R); band = isinstance(R, tuple)
    g = np.arange(0, ex + 1e-9, h); gy = np.arange(0, ey - 0.5 + 1e-9, h)
    ts = np.arange(-0.5, 0.5 + 1e-9, h)
    TX, TY = np.meshgrid(ts, ts, indexing='ij'); TX = TX.ravel(); TY = TY.ravel()
    X, Y = np.meshgrid(g + 0.5, gy + 0.5, indexing='ij'); X = X.ravel(); Y = Y.ravel()
    k = (Y <= X + 1e-6) | band; X = X[k]; Y = Y[k]
    out = []
    for th in thetas:
        cx = (X[:, None] + abs(th) * TX[None, :]).ravel(); cy = (Y[:, None] + abs(th) * TY[None, :]).ravel()
        out.append(np.c_[cx, cy, np.full(len(cx), th % (math.pi / 2))])
    return clipQ(np.concatenate(out), R)


# ============================================================================ oracle (worker functions)
_G = {}


def _init(m, x):
    _G['m'] = m; _G['x'] = x


def _eval(P):
    return _G['m'].contrib(P, xval=_G['x'])


def _polish(args):
    P0, seed, rounds = args
    m = _G['m']; x = _G['x']; rng = np.random.default_rng(seed)
    cur = P0.copy(); curv = m.contrib(cur, xval=x); found = []
    nper = 24
    for r in range(rounds):
        rad = 0.05 * 0.4 ** r; arad = math.radians(2.0) * 0.4 ** r
        Q = np.repeat(cur, nper, 0)
        d = rng.normal(size=Q.shape); d[:, :2] *= rad; d[:, 2] *= arad
        k = rng.random(len(Q)) < 0.3; d[k, 2] = 0
        Q = Q + d
        k = rng.random(len(Q)) < 0.2; Q[k, 1] = halfwidth(Q[k, 2]) + WALL
        if not isinstance(m.G, tuple): k = rng.random(len(Q)) < 0.05; Q[k, 0] = halfwidth(Q[k, 2]) + WALL
        Q[:, 2] = np.mod(Q[:, 2], math.pi / 2)
        Q = clipQ(Q, m.G)
        v = m.contrib(Q, xval=x)
        V = v.reshape(len(cur), nper); j = V.argmin(1)
        better = V[np.arange(len(cur)), j] < curv
        cur[better] = Q.reshape(len(cur), nper, 3)[np.arange(len(cur)), j][better]
        curv[better] = V[np.arange(len(cur)), j][better]
        bad = v < 1 - 1e-7
        if bad.any(): found.append((Q[bad], v[bad]))
    fp = np.concatenate([f[0] for f in found]) if found else np.zeros((0, 3))
    fv = np.concatenate([f[1] for f in found]) if found else np.zeros(0)
    return cur, curv, fp, fv


def oracle(pool, m, x, rng, nrand=200000, lat_pitch=0.02, lat_degs=None, npol=600, prounds=8, nproc=5):
    """float separation: lattice + random poses, then local polish of the worst.  Returns (min, violated poses,
    their values, worst pose)."""
    if lat_degs is None: lat_degs = [0, 1e-5, 1e-3, 0.5, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 88, 89.5]
    P = np.r_[lattice_poses(m.G, lat_pitch, lat_degs), random_poses(m.G, nrand, rng)]
    # shifted copy of the lattice (off-lattice), axis-parallel
    P = np.r_[P, clipQ(lattice_poses(m.G, lat_pitch, [0, 1e-5]) + np.array([lat_pitch * 0.37, lat_pitch * 0.61, 0]), m.G)]
    P = np.r_[P, near_lattice_poses(m.G, m.hlat, rng, nrand // 5), germ_poses(m.G, m.hlat)]
    ch = np.array_split(P, max(nproc * 4, 8))
    v = np.concatenate(pool.map(_eval, ch))
    order = np.argsort(v)
    # diverse seeds: worst per 0.05-bucket
    key = np.floor(P[order, :2] / 0.05).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
    _, first = np.unique(kk, return_index=True)
    seeds = P[order[np.sort(first)][:npol]]
    res = pool.map(_polish, [(s, i + int(rng.integers(1 << 30)), prounds) for i, s in enumerate(np.array_split(seeds, max(nproc * 2, 4)))])
    cur = np.concatenate([r[0] for r in res]); curv = np.concatenate([r[1] for r in res])
    fp = np.concatenate([P[v < 1 - 1e-7]] + [r[2] for r in res]); fv = np.concatenate([v[v < 1 - 1e-7]] + [r[3] for r in res])
    allv = np.r_[v, curv]; allp = np.r_[P, cur]
    j = int(allv.argmin())
    return float(allv[j]), fp, fv, allp[j]


class Serial:
    def map(self, f, it): return [f(t) for t in it]
    def terminate(self): pass


def pick_rows(fp, fv, maxn=4000, block=0.01):
    if len(fp) == 0: return fp
    o = np.argsort(fv); fp = fp[o]
    key = np.c_[np.floor(fp[:, :2] / block), np.floor(fp[:, 2] / math.radians(0.5))].astype(np.int64)
    kk = (key[:, 0] * 1000 + key[:, 1]) * 1000 + key[:, 2]
    _, first = np.unique(kk, return_index=True)
    return fp[np.sort(first)][:maxn]


# ============================================================================ LP
class LP:
    def __init__(self, m):
        import highspy
        self.hp = highspy; self.m = m
        h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('random_seed', 0)
        inf = highspy.kHighsInf; self.inf = inf
        n = m.nvar
        h.addVars(n, np.zeros(n), np.full(n, inf))
        h.changeColsCost(n, np.arange(n, dtype=np.int32), m.cobj.astype(float))
        # sigma = 0: profile mass per period = w
        idx = np.nonzero(m.csig)[0].astype(np.int32)
        if len(idx):
            rhs = m.w - m.const_sig
            h.addRow(rhs, rhs, len(idx), idx, m.csig[idx].astype(float))
        self.h = h; self.nrows = 0; self.P = np.zeros((0, 3)); self.solver = 'simplex'; self.crossover = 'off'

    def add(self, P):
        A, const = self.m.contrib(P)
        A = A.tocsr(); A.sort_indices()
        lo = 1.0 - const
        need = lo > 1e-12
        # rows with no variable and const < 1 are infeasible outright (fixed structure violated)
        A = A[need]; lo = lo[need]; P = P[need]
        if A.shape[0] == 0: return 0
        self.h.addRows(A.shape[0], lo, np.full(A.shape[0], self.inf), int(A.nnz), A.indptr[:-1].astype(np.int32),
                       A.indices.astype(np.int32), A.data.astype(float))
        self.nrows += A.shape[0]; self.P = np.r_[self.P, P]
        return A.shape[0]

    def solve(self):
        h = self.h
        if self.solver == 'ipm':
            h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', self.crossover)
        elif self.solver == 'pdlp':
            h.setOptionValue('solver', 'pdlp')
        h.run()
        st = h.getModelStatus()
        if st != self.hp.HighsModelStatus.kOptimal and self.solver == 'ipm':
            h.setOptionValue('run_crossover', 'on'); h.run(); st = h.getModelStatus()
            h.setOptionValue('run_crossover', self.crossover)
        if st != self.hp.HighsModelStatus.kOptimal:
            h.setOptionValue('solver', 'simplex'); h.run(); st = h.getModelStatus()
        if st != self.hp.HighsModelStatus.kOptimal:
            return None, None, str(st)
        x = np.maximum(np.array(h.getSolution().col_value), 0.0)
        if self.solver != 'simplex': x[x < 1e-10] = 0.0
        return x, float(h.getInfo().objective_function_value), 'ok'


# ============================================================================ drivers
def build(a):
    mode = a.mode
    if mode == 'lebwall':
        m = Model(a.R, 0.0, wleb=0.0)
        add_free_corner(m, a.hc, a.hca)
    elif mode == 'naga':
        m = Model(2, 1.0)
        add_naga_corner(m); add_naga_profile(m)
    elif mode == 'nagaprof':
        m = Model(a.R, 1.0)
        add_naga_profile(m); add_free_corner(m, a.hc, a.hca)
    elif mode == 'band':
        m = Model(0, a.w); m.band = True; m.G = ('band', float(a.w))
        add_free_profile(m, a.hp, a.hpa)
        m.cobj = [-c for c in m.cobj]            # maximise m_v
    else:
        m = Model(a.R, a.w)
        add_free_corner(m, a.hc, a.hca); add_free_profile(m, a.hp, a.hpa)
    m.finalize()
    m.hlat = min(a.hc, a.hp) if mode != 'naga' else 0.1
    if mode == 'nagaprof': m.hlat = min(a.hc, 0.1)
    return m


def describe(m, x, thr=1e-7):
    out = []
    for i in np.nonzero(x > thr)[0]:
        out.append((m.label[i], float(x[i])))
    return out


def run_lp(a):
    tag = a.tag or f"{a.mode}_R{a.R}_w{a.w:g}_hc{a.hc:g}_hp{a.hp:g}"
    od = os.path.join(RUNS, 'quad_' + tag); os.makedirs(od, exist_ok=True)
    logf = open(os.path.join(od, 'log.txt'), 'a')

    def log(*s):
        msg = ' '.join(str(t) for t in s); print(msg, flush=True); logf.write(msg + '\n'); logf.flush()
    m = build(a)
    log(f"== {tag}  mode={a.mode} R={m.R} w={m.w} wleb={m.wleb} nvar={m.nvar} elements="
        f"{ {k: len(v) for k, v in m.A.items()} }  args={vars(a)}")
    rng = np.random.default_rng(a.seed)
    lp = LP(m); lp.solver = a.solver; lp.crossover = a.crossover
    t0 = time.time()
    P0 = np.r_[lattice_poses(m.G, a.pitch, a.degs), near_lattice_poses(m.G, m.hlat, full=False)]
    n0 = lp.add(P0)
    log(f"initial rows: {n0} (pitch {a.pitch}, degs {a.degs})  [{time.time()-t0:.1f}s]")
    if a.warm:
        for wf in a.warm.split(','):
            W = np.load(wf)['rows']; nw = lp.add(W)
            log(f"warm rows from {wf}: {nw}")
    hist = []
    pool = None
    x = None
    for rnd in range(a.rounds + 1):
        t1 = time.time()
        x, obj, st = lp.solve()
        if x is None:
            log(f"round {rnd}: LP status {st}"); break
        D = m.R ** 2 - obj - m.const_obj
        MC = float(sum(m.cobj[i] * x[i] for i in range(m.nvar) if m.label[i][0][0] == 'C'))
        mv = float(sum(m.cobj[i] * x[i] for i in range(m.nvar) if m.label[i][0][0] == 'P')) + m.const_obj
        tlp = time.time() - t1
        if rnd == a.rounds and not a.final_oracle:
            hist.append(dict(round=rnd, D=D, MC=MC, mv=mv, rows=lp.nrows)); log(f"round {rnd}: D={D:.6f} (no oracle)"); break
        if a.nproc > 1:
            if pool is not None: pool.terminate()
            pool = mp.Pool(a.nproc, initializer=_init, initargs=(m, x))
        else:
            _init(m, x); pool = Serial()
        t2 = time.time()
        mn, fp, fv, wp = oracle(pool, m, x, rng, nrand=a.nrand, lat_pitch=a.oracle_pitch, npol=a.npol, nproc=a.nproc)
        rows = pick_rows(fp, fv, maxn=a.maxadd)
        hist.append(dict(round=rnd, D=D, MC=MC, mv=mv, rows=lp.nrows, oracle_min=mn, nviol=int(len(fp)),
                         worst=[float(t) for t in wp], tlp=tlp, tor=time.time() - t2))
        log(f"round {rnd}: D={D:.6f} M_C={MC:.5f} m_v={mv:.5f} rows={lp.nrows} | oracle min={mn:.6f} at "
            f"({wp[0]:.5f},{wp[1]:.5f},{math.degrees(wp[2]):.4f}deg) viol={len(fp)} add={len(rows)} "
            f"[lp {tlp:.1f}s, oracle {time.time()-t2:.1f}s]")
        json.dump(hist, open(os.path.join(od, 'hist.json'), 'w'), indent=1)
        np.savez(os.path.join(od, 'sol.npz'), x=x, D=D, round=rnd, rows=lp.P)
        if mn >= 1 - a.stop_tol:
            clean = True
            for extra in range(2):
                mn2, fp2, fv2, wp2 = oracle(pool, m, x, rng, nrand=2 * a.nrand, lat_pitch=a.oracle_pitch * 0.77, npol=a.npol, nproc=a.nproc)
                log(f"   confirm pass {extra}: oracle min={mn2:.7f} viol={len(fp2)}")
                if mn2 < 1 - a.stop_tol:
                    rows = pick_rows(fp2, fv2, maxn=a.maxadd); clean = False; break
            if clean:
                log(f"stable: oracle min >= 1 - {a.stop_tol} in 3 passes"); break
        if rnd == a.rounds: break
        if len(rows) == 0: break
        lp.add(rows)
    if pool is not None: pool.terminate()
    if x is not None:
        sup = describe(m, x)
        with open(os.path.join(od, 'support.txt'), 'w') as f:
            for lab, val in sorted(sup, key=lambda t: str(t[0])):
                f.write(f"{lab}\t{val:.9f}\n")
        log(f"support size {len(sup)}; top: " + '; '.join(f"{l}={v:.4f}" for l, v in sorted(sup, key=lambda t: -t[1])[:25]))
    return hist


def run_fixed(a):
    """evaluate a fully fixed structure (naga): oracle min and D."""
    m = build(a)
    x = np.zeros(max(m.nvar, 1))
    rng = np.random.default_rng(a.seed)
    if a.nproc > 1: pool = mp.Pool(a.nproc, initializer=_init, initargs=(m, x))
    else: _init(m, x); pool = Serial()
    mn, fp, fv, wp = oracle(pool, m, x, rng, nrand=a.nrand, lat_pitch=a.oracle_pitch, npol=a.npol, nproc=a.nproc)
    pool.terminate()
    print(f"fixed structure {a.mode}: oracle min = {mn:.6f} at ({wp[0]:.6f},{wp[1]:.6f},{math.degrees(wp[2]):.5f} deg); "
          f"violations found {len(fp)};  const_sig={m.const_sig} const_obj(m_v)={m.const_obj}")
    return m, mn, wp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['lp', 'lebwall', 'naga', 'nagaprof', 'band'])
    ap.add_argument('--R', type=int, default=2); ap.add_argument('--w', type=float, default=1.0)
    ap.add_argument('--hc', type=float, default=0.1); ap.add_argument('--hca', type=float, default=0.25)
    ap.add_argument('--hp', type=float, default=0.1); ap.add_argument('--hpa', type=float, default=0.25)
    ap.add_argument('--pitch', type=float, default=0.1)
    ap.add_argument('--degs', type=lambda s: [float(t) for t in s.split(',')], default=[0, 1e-4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85])
    ap.add_argument('--rounds', type=int, default=30); ap.add_argument('--min-rounds', type=int, default=3)
    ap.add_argument('--final-oracle', action='store_true', default=True)
    ap.add_argument('--nrand', type=int, default=1000000); ap.add_argument('--oracle-pitch', type=float, default=0.02)
    ap.add_argument('--npol', type=int, default=2000); ap.add_argument('--maxadd', type=int, default=4000)
    ap.add_argument('--stop-tol', type=float, default=1e-6)
    ap.add_argument('--nproc', type=int, default=5); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tag', default='')
    ap.add_argument('--warm', default='')
    ap.add_argument('--solver', default='ipm'); ap.add_argument('--crossover', default='off')
    a = ap.parse_args()
    if a.mode == 'naga': run_fixed(a)
    else: run_lp(a)


if __name__ == '__main__':
    main()
