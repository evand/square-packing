#!/usr/bin/env python3
"""closed4.py -- HEURISTIC LP for the minimum-weight point cover of [0,s]^2 in CLOSED semantics.

Question (search/CLOSED4.md).  At the container [0,4]^2 exactly: what is the least total weight of
a weighted point set such that every CLOSED unit square contained in [0,4]^2, at every angle,
captures weight >= 1, a point on the boundary of the square counting as captured?  Friedman's
14-point set and Nagamochi's weighted set give 14; a value < 12 would prove s(12) = 4.

Semantics.  The captured weight  w(Q) = sum_p w_p [p in Q]  (Q closed) is upper semicontinuous in
the pose (c, theta) of Q, so "w(Q) >= 1 for every closed unit square with centre in the OPEN
admissible box" already implies it for the wall-touching squares (limits).  Rows here are therefore
placed in the open admissible box (offset 1e-7 from its boundary, so the near-wall limit is
sampled); a consequence is that points ON the container walls are captured by no row and get
weight 0 -- they are useless in closed semantics for exactly this reason (a wall-touching square
is dominated by its interior neighbours, which do not contain the wall).  They are still offered as
columns (--wall-points, default on) so that the LP confirms this.

This is a SAMPLED LP: rows are finitely many poses (a lattice of centres of pitch --pitch in the
rotated frame of each angle, angles = 0, a geometric sequence of small angles, every degree up to
45, and 45 - small; plus near-wall band rows, plus cutting planes found by a local search around
the worst poses), columns are a lattice of pitch --col-pitch that contains the integer and
half-integer lines exactly, the points of the three literature sets, and points added by pricing
on the dual (as tighten.py).  Nothing here is a certificate: the LP value is neither an upper nor
a lower bound on COVER(s).  The honest number for a finished point set is
        cost = (total weight) / (min captured weight over a dense scan)
which `stress` computes (still a scan, not a proof; the exact verifier in verify/ would be the
proof, and at s = 4 it is not needed unless the total is < 12).

Usage
    python3 search/closed4.py run   [--s 4] [--tag closed4] [--time 5400] [--nproc 8] ...
    python3 search/closed4.py sanity          # checks (i)-(iii) and the axis-aligned reference
    python3 search/closed4.py stress runs/closed4_best.txt [--nproc 8]
    python3 search/closed4.py axis   [--s 4]  # LP with axis-aligned rows only (expect 9)

Outputs (runs/): closed4_<tag>.log, closed4_<tag>.json (history), closed4_<tag>_best.txt (certificate
format, D = 1000, W = 1e7, weights rounded UP), closed4_<tag>_best.json, closed4_<tag>_best.png.
"""
import sys, os, math, time, json, argparse, warnings
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')   # shared machine: one thread per process
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction

TOL = 1e-9           # closed containment: |q - u|_inf <= 1/2 + TOL counts (safe direction for the cover)
# --margin MU sets TOL = -MU: a point then counts only if it lies in the CONCENTRIC square of side
# 1 - 2*MU, i.e. at distance >= MU inside every edge of the unit square.  This is the "stable"
# (margin-MU) cover problem of search/RUNG2.md: a solution of it is a closed cover with a uniform
# geometric margin, which is exactly what the box checker search/zeromargin.py needs in order to
# terminate (a witness set fixed over a whole pose box).  Equivalently, by scaling, it is the
# ordinary closed cover problem for the container [0, s/(1-2MU)]^2.
WALL = 1e-7          # rows are placed >= WALL inside the open admissible box
HIGHS = {'random_seed': 0}

# ----------------------------------------------------------------------------- literature sets
FRIEDMAN14 = [(1, 1), (1.6, 1), (2.4, 1), (3, 1), (1, 1.8), (2, 1.8), (3, 1.8), (1, 2.2), (2, 2.2), (3, 2.2),
              (1, 3), (1.6, 3), (2.4, 3), (3, 3)]                                   # DS7 Thm 4 (n = 15), [0,4]^2


def d4(x, y, s):
    return [(x, y), (s - x, y), (x, s - y), (s - x, s - y), (y, x), (s - y, x), (y, s - x), (s - y, s - x)]


def bentz16(s=4.0):
    """B10 Fig. 2: A(1,0.914), C(0.914,2), D(1.65,1.65) and their images under x->4-x, y->4-y, x<->y."""
    pts = set()
    for (x, y) in [(1, 0.914), (0.914, 2), (1.65, 1.65)]:
        for (a, b) in d4(x, y, s): pts.add((round(a, 6), round(b, 6)))
    return sorted(pts)


def nagamochi_points(s=4.0):
    """N05 Sec. 3 for a = b = 4: Q = 8 segment endpoints (weight 0.45), P = 4 points (weight 0.5).
    (The set also has the four segments L_i on the lines x,y in {1,3}, |x|,|y| in [0.9,3.1] at 0.5 per
    unit length, and the area [1,3]^2 at 1 per unit area -- those are sampled by the column lattice.)"""
    Q = [(0.9, 1), (3.1, 1), (0.9, 3), (3.1, 3), (1, 0.9), (1, 3.1), (3, 0.9), (3, 3.1)]
    P = [(2, 0.9), (2, 3.1), (0.9, 2), (3.1, 2)]
    return [(x, y, 0.45) for x, y in Q] + [(x, y, 0.5) for x, y in P]


def angle_list(deg_step=1.0):
    """degrees: 0, geometric small angles, every degree to 45, and 45 - small (x<->y maps 45+d to 45-d)."""
    small = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3]
    degs = [0.0] + small + [0.5] + list(np.arange(1.0, 45.0 + 1e-9, deg_step)) + [45 - d for d in (0.5, 0.1, 0.03, 0.01, 0.003, 0.001)]
    return sorted(set(round(d, 6) for d in degs))


# ----------------------------------------------------------------------------- geometry
def captured(P, w, poses, tol=None, chunk=256):
    """captured weight of the CLOSED unit square at each pose (cx, cy, theta); brute force."""
    if tol is None: tol = TOL
    poses = np.asarray(poses, dtype=float).reshape(-1, 3); out = np.empty(len(poses))
    for i in range(0, len(poses), chunk):
        B = poses[i:i + chunk]; ct = np.cos(B[:, 2]); st = np.sin(B[:, 2])
        u0 = B[:, 0] * ct + B[:, 1] * st; u1 = -B[:, 0] * st + B[:, 1] * ct
        q0 = P[:, 0:1] * ct[None, :] + P[:, 1:2] * st[None, :]        # (n, m)
        q1 = -P[:, 0:1] * st[None, :] + P[:, 1:2] * ct[None, :]
        ins = np.maximum(np.abs(q0 - u0[None, :]), np.abs(q1 - u1[None, :])) <= 0.5 + tol
        out[i:i + chunk] = w @ ins
    return out


def admissible_box(s, th):
    wid = abs(math.cos(th)) + abs(math.sin(th))
    return wid / 2, s - wid / 2


def scan_angle(P, w, s, th, pitch, tol=None, band=True):
    """captured weight on a lattice of centres of pitch `pitch` in the rotated frame (u = R(-th) c),
    restricted to the OPEN admissible box (offset WALL), plus (band=True) explicit rows along the
    four walls at offset WALL (c-space lattice), i.e. the near-wall limits.
    Returns (values, poses (m,3))."""
    if tol is None: tol = TOL
    ct, st = math.cos(th), math.sin(th)
    q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
    lo, hi = admissible_box(s, th); lo += WALL; hi -= WALL
    if hi <= lo: return np.zeros(0), np.zeros((0, 3))
    corners = np.array([[lo, lo], [hi, lo], [hi, hi], [lo, hi]])
    cu0 = corners[:, 0] * ct + corners[:, 1] * st; cu1 = -corners[:, 0] * st + corners[:, 1] * ct
    g0 = np.arange(cu0.min(), cu0.max() + 1e-12, pitch); g1 = np.arange(cu1.min(), cu1.max() + 1e-12, pitch)
    G0, G1 = np.meshgrid(g0, g1, indexing='ij')
    cx = G0 * ct - G1 * st; cy = G0 * st + G1 * ct
    ok = (cx >= lo) & (cx <= hi) & (cy >= lo) & (cy <= hi)
    # range sums: for each u0 row, atoms with |q0-u0| <= h, cumsum in q1 order
    h = 0.5 + tol
    oy = np.argsort(q1); q1s = q1[oy]; wy = w[oy]; q0y = q0[oy]
    lo1 = np.searchsorted(q1s, g1 - h, 'left'); hi1 = np.searchsorted(q1s, g1 + h, 'right')
    vals = np.empty((len(g0), len(g1)))
    for i in range(len(g0)):
        if not ok[i].any(): vals[i] = 9e9; continue
        mask = np.abs(q0y - g0[i]) <= h
        cw = np.r_[0.0, np.cumsum(np.where(mask, wy, 0.0))]
        vals[i] = cw[hi1] - cw[lo1]
    vals = vals[ok]; poses = np.c_[cx[ok], cy[ok], np.full(int(ok.sum()), th)]
    if band:
        t = np.arange(lo, hi + 1e-12, pitch)
        if t[-1] < hi - 1e-9: t = np.r_[t, hi]
        bp = np.concatenate([np.c_[np.full(len(t), lo), t], np.c_[np.full(len(t), hi), t],
                             np.c_[t, np.full(len(t), lo)], np.c_[t, np.full(len(t), hi)]])
        bp = np.c_[bp, np.full(len(bp), th)]
        bv = captured(P, w, bp, tol)
        vals = np.r_[vals, bv]; poses = np.r_[poses, bp]
    return vals, poses


def _scan_worker(args):
    P, w, s, th, pitch, thr, perang, block = args
    vals, poses = scan_angle(P, w, s, th, pitch)
    if len(vals) == 0: return 9e9, 0, np.zeros((0, 3)), np.zeros(0)
    gmin = float(vals.min()); bad = vals < thr; nviol = int(bad.sum())
    if nviol == 0: return gmin, 0, np.zeros((0, 3)), np.zeros(0)
    bv = vals[bad]; bp = poses[bad]
    # thin: worst per block of side `block` in c-space
    if block > 0:
        key = np.floor(bp[:, :2] / block).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
        o = np.argsort(bv, kind='stable'); seen = {}; keep = []
        for j in o:
            k = kk[j]
            if k in seen: continue
            seen[k] = 1; keep.append(j)
            if len(keep) >= perang: break
        bp = bp[keep]; bv = bv[keep]
    else:
        o = np.argsort(bv)[:perang]; bp = bp[o]; bv = bv[o]
    return gmin, nviol, bp, bv


def separate(P, w, s, thetas, pitch, thr=1.0 - 1e-9, perang=150, block=0.1, pool=None):
    """lattice separation over all angles: (global min, #violated, poses (m,3), values)."""
    jobs = [(P, w, s, th, pitch, thr, perang, block) for th in thetas]
    outs = pool.map(_scan_worker, jobs) if pool is not None else [_scan_worker(j) for j in jobs]
    gmin = min(o[0] for o in outs); nviol = sum(o[1] for o in outs)
    poses = np.concatenate([o[2] for o in outs]) if outs else np.zeros((0, 3))
    vals = np.concatenate([o[3] for o in outs]) if outs else np.zeros(0)
    return gmin, nviol, poses, vals


def _polish_worker(args):
    P, w, s, seeds, rounds, nper, seed, thmax = args
    return polish(P, w, s, seeds, rounds=rounds, nper=nper, rng=np.random.default_rng(seed), thmax=thmax,
                  arad0=(0.0 if thmax == 0 else math.radians(1.0)))


def polish(P, w, s, seeds, rounds=7, nper=48, rng=None, radius0=0.03, arad0=math.radians(1.0), thmax=math.pi / 4):
    """local minimisation of the captured weight around seed poses by random perturbation of
    (cx, cy, theta) with shrinking radii; also tries theta -> 0 / 45deg limits and exact grid-aligned
    angles.  Returns (poses, values) of every distinct pose visited with value < 1 - 1e-9, plus the
    best pose per seed (for diagnostics)."""
    if rng is None: rng = np.random.default_rng(0)
    seeds = np.asarray(seeds, dtype=float).reshape(-1, 3)
    if len(seeds) == 0: return np.zeros((0, 3)), np.zeros(0), np.zeros((0, 3)), np.zeros(0)
    cur = seeds.copy(); curv = captured(P, w, cur)
    found = {}
    def clip(Q):
        Q[:, 2] = np.clip(Q[:, 2], 0.0, thmax)
        wid = np.abs(np.cos(Q[:, 2])) + np.abs(np.sin(Q[:, 2])); lo = wid / 2 + WALL; hi = s - wid / 2 - WALL
        Q[:, 0] = np.clip(Q[:, 0], lo, hi); Q[:, 1] = np.clip(Q[:, 1], lo, hi)
        return Q
    for r in range(rounds):
        rad = radius0 * (0.3 ** r); arad = arad0 * (0.3 ** r)
        m = len(cur)
        Q = np.repeat(cur, nper, axis=0)
        d = rng.normal(size=(m * nper, 3)); d[:, :2] *= rad; d[:, 2] *= arad
        # a share of the trials keeps the angle exactly, or jumps to the special angles
        keep = rng.random(m * nper) < 0.3; d[keep, 2] = 0.0
        Q = Q + d
        spec = rng.random(m * nper) < (0.08 if thmax > 0 else 0.0)
        specials = [0.0, math.pi / 4, 1e-6, 1e-5, 1e-4, 1e-3, math.pi / 4 - 1e-6, math.pi / 4 - 1e-4]
        if thmax > math.pi / 4: specials += [math.pi / 2 - v for v in specials]
        Q[spec, 2] = rng.choice(specials, size=int(spec.sum()))
        Q = clip(Q); v = captured(P, w, Q)
        V = v.reshape(m, nper); j = V.argmin(axis=1); better = V[np.arange(m), j] < curv
        cur[better] = Q.reshape(m, nper, 3)[np.arange(m), j][better]; curv[better] = V[np.arange(m), j][better]
        bad = v < 1 - 1e-9
        for q, vv in zip(Q[bad], v[bad]):
            key = (round(q[0], 6), round(q[1], 6), round(q[2], 8))
            if key not in found or vv < found[key][1]: found[key] = (q, vv)
    if found:
        fp = np.array([f[0] for f in found.values()]); fv = np.array([f[1] for f in found.values()])
    else:
        fp = np.zeros((0, 3)); fv = np.zeros(0)
    return fp, fv, cur, curv


# ----------------------------------------------------------------------------- LP model
class Model:
    """columns = D4 orbits (sym=True) or single points (sym=False); rows = poses (cx, cy, theta).
    A_rk = number of atoms of column k inside the closed square of row r."""

    def __init__(self, s, D=1000, sym=True):
        self.s = s; self.D = D; self.sym = sym; self.K = int(round(s * D))
        assert abs(self.K - s * D) < 1e-9, "s must be a multiple of 1/D"
        self.orbits = []; self.okey = {}; self.sizes = np.zeros(0)
        self.P = np.zeros((0, 2)); self.own = np.zeros(0, dtype=np.int64)
        self.rows = []; self.rkey = set(); self.R = []; self.C = []; self.V = []
        self.tag = {}

    def add_point(self, x, y, tag=''):
        X = int(round(x * self.D)); Y = int(round(y * self.D)); K = self.K
        X = min(max(X, 0), K); Y = min(max(Y, 0), K)
        if self.sym:
            imgs = sorted(set([(X, Y), (K - X, Y), (X, K - Y), (K - X, K - Y), (Y, X), (K - Y, X), (Y, K - X), (K - Y, K - X)]))
        else:
            imgs = [(X, Y)]
        key = imgs[0]
        if key in self.okey: return False
        k = len(self.orbits); self.okey[key] = k; self.tag[k] = tag
        of = np.array(imgs, dtype=float) / self.D; self.orbits.append(of)
        self.sizes = np.append(self.sizes, float(len(of)))
        self.P = np.concatenate([self.P, of]); self.own = np.concatenate([self.own, np.full(len(of), k)])
        if self.rows:
            Q = np.array(self.rows); ct = np.cos(Q[:, 2]); st = np.sin(Q[:, 2])
            u0 = Q[:, 0] * ct + Q[:, 1] * st; u1 = -Q[:, 0] * st + Q[:, 1] * ct
            tot = np.zeros(len(Q))
            for (px, py) in of:
                q0 = px * ct + py * st; q1 = -px * st + py * ct
                tot += (np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= 0.5 + TOL)
            hit = np.nonzero(tot)[0]
            if hit.size:
                self.R.append(hit.astype(np.int32)); self.C.append(np.full(hit.size, k, dtype=np.int32)); self.V.append(tot[hit])
        return True

    def canon(self, cx, cy, th):
        """row key modulo the C4 rotations (which fix theta); reflections were used for theta in [0,45]."""
        s = self.s
        if self.sym:
            imgs = [(cx, cy), (s - cy, cx), (s - cx, s - cy), (cy, s - cx)]
            a, b = min((round(x, 7), round(y, 7)) for x, y in imgs)
        else:
            a, b = round(cx, 7), round(cy, 7)
        return (a, b, round(th, 9))

    def add_rows(self, poses):
        P, own = self.P, self.own; added = 0
        for cx, cy, th in np.asarray(poses, dtype=float).reshape(-1, 3):
            key = self.canon(cx, cy, th)
            if key in self.rkey: continue
            self.rkey.add(key)
            ct, st = math.cos(th), math.sin(th)
            u0 = cx * ct + cy * st; u1 = -cx * st + cy * ct
            q0 = P[:, 0] * ct + P[:, 1] * st; q1 = -P[:, 0] * st + P[:, 1] * ct
            ok = np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= 0.5 + TOL
            cnt = np.bincount(own[ok], minlength=len(self.orbits)); idx = np.nonzero(cnt)[0]
            r = len(self.rows); self.rows.append((float(cx), float(cy), float(th)))
            if idx.size:
                self.R.append(np.full(idx.size, r, dtype=np.int32)); self.C.append(idx.astype(np.int32)); self.V.append(cnt[idx].astype(float))
            added += 1
        return added

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
        self.rkey = set(self.canon(*rw) for rw in self.rows)

    def solve(self, fixed=None, method='highs'):
        """min sizes.x  s.t.  A x >= 1, x >= 0  [x = fixed if given: just evaluate]."""
        A = self.matrix()
        if fixed is not None:
            return float(self.sizes @ fixed), fixed, None, A @ fixed
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c=self.sizes, A_ub=-A, b_ub=-np.ones(A.shape[0]), bounds=(0, None), method=method, options=dict(HIGHS))
        if not res.success: return None
        y = np.maximum(-res.ineqlin.marginals, 0.0)
        return res.fun, res.x, y, A @ res.x


def price(m, y, pitch=0.01, want=200, ysup=1e-9):
    """reduced-cost pricing on the D4-symmetrised dual (as tighten.price): symmetrised closed-square
    coverage of the dual measure on a pitch-grid (includes the grid lines exactly), refined on the
    1/D grid; returns (coverage, X, Y) with coverage > 1, best first."""
    Q = np.array(m.rows); sup = y > ysup; Q = Q[sup]; yy = y[sup]
    if len(yy) == 0: return []
    s = m.s; cx0, cy0, t0 = Q[:, 0], Q[:, 1], Q[:, 2]
    if m.sym:
        imgs = [(cx0, cy0, t0), (s - cx0, cy0, -t0), (cx0, s - cy0, -t0), (s - cx0, s - cy0, t0),
                (cy0, cx0, -t0), (s - cy0, cx0, t0), (cy0, s - cx0, t0), (s - cy0, s - cx0, -t0)]
        yy8 = yy / 8.0
    else:
        imgs = [(cx0, cy0, t0)]; yy8 = yy
    R = []
    for a, b, th in imgs:
        ct_, st_ = np.cos(th), np.sin(th)
        R.append(np.c_[a * ct_ + b * st_, -a * st_ + b * ct_, ct_, st_, yy8])
    R = np.concatenate(R)
    def cov_at(G):
        out = np.zeros(len(G))
        for i in range(0, len(R), 2000):
            B = R[i:i + 2000]
            q0 = G[:, 0:1] * B[:, 2] + G[:, 1:2] * B[:, 3]; q1 = -G[:, 0:1] * B[:, 3] + G[:, 1:2] * B[:, 2]
            ins = np.maximum(np.abs(q0 - B[:, 0]), np.abs(q1 - B[:, 1])) <= 0.5 + TOL
            out += ins @ B[:, 4]
        return out
    n0 = int(round(s / pitch)); gx = np.arange(n0 + 1) * pitch
    G0, G1 = np.meshgrid(gx, gx, indexing='ij'); G = np.c_[G0.ravel(), G1.ravel()]
    if m.sym: G = G[(G[:, 0] <= G[:, 1] + 1e-12) & (G[:, 1] <= s / 2 + 1e-12)]
    c = cov_at(G); o = np.argsort(c)[::-1][:want]
    cand = []; D = m.D; r_ = int(round(pitch * D))
    for i in o:
        if c[i] <= 1.0 + 1e-7: break
        X0, Y0 = int(round(G[i, 0] * D)), int(round(G[i, 1] * D))
        xs = np.arange(max(0, X0 - r_), min(m.K, X0 + r_) + 1); ys_ = np.arange(max(0, Y0 - r_), min(m.K, Y0 + r_) + 1)
        GX, GY = np.meshgrid(xs, ys_, indexing='ij'); Gf = np.c_[GX.ravel(), GY.ravel()] / D
        cf = cov_at(Gf); j = int(np.argmax(cf))
        if cf[j] > 1.0 + 1e-7: cand.append((float(cf[j]), int(GX.ravel()[j]), int(GY.ravel()[j])))
    cand.sort(reverse=True)
    return cand


# ----------------------------------------------------------------------------- columns
def column_lattice(s, pitch):
    """pitch-lattice from both walls (so the integer/half-integer lines from either wall are exact)."""
    g = np.arange(0.0, s + 1e-9, pitch)
    g = np.unique(np.round(np.concatenate([g, s - g]), 6))
    return [(x, y) for x in g for y in g]


def build_columns(m, s, col_pitch, wall_points=True, literature=True, lattice=True):
    n = 0
    if lattice:
        for x, y in column_lattice(s, col_pitch):
            if not wall_points and (min(x, y) < 1e-9 or max(x, y) > s - 1e-9): continue
            n += m.add_point(x, y, 'lattice')
    if literature:
        for x, y in FRIEDMAN14: n += m.add_point(x, y, 'friedman')
        for x, y in bentz16(4.0):
            if max(x, y) <= s: n += m.add_point(x, y, 'bentz')
        for x, y, wt in nagamochi_points(4.0):
            if max(x, y) <= s: n += m.add_point(x, y, 'nagamochi')
    return n


# ----------------------------------------------------------------------------- export
def export(m, x, path, WD=10 ** 7, up=True):
    s = Fraction(m.K, m.D); lines = []; tot = 0
    for k, o in enumerate(m.orbits):
        if x[k] <= 1e-12: continue
        v = x[k] * WD; wq = int(math.ceil(v - 1e-9)) if up else int(math.floor(v))
        if wq <= 0: continue
        for px, py in o:
            lines.append((int(round(px * m.D)), int(round(py * m.D)), wq)); tot += wq
    with open(path, 'w') as f:
        f.write(f"{s.numerator} {s.denominator}\n{m.D}\n{WD}\n{len(lines)}\n")
        for X, Y, W in lines: f.write(f"{X} {Y} {W}\n")
    return tot / WD, len(lines)


def read_cert(path):
    t = open(path).read().split(); it = iter(t)
    sn, sd, D, WD, n = [int(next(it)) for _ in range(5)]
    rows = np.array([[int(next(it)), int(next(it)), int(next(it))] for _ in range(n)], dtype=np.int64)
    return Fraction(sn, sd), D, WD, rows


def weight_breakdown(P, w, s, D=1000):
    """where the weight sits: walls, integer lines (x or y integer, 1..s-1), half-integer lines, other."""
    X = np.round(P[:, 0] * D).astype(np.int64); Y = np.round(P[:, 1] * D).astype(np.int64); K = int(round(s * D))
    wall = (X == 0) | (Y == 0) | (X == K) | (Y == K)
    onint = lambda Z: (Z % D == 0)
    onhalf = lambda Z: (Z % (D // 2) == 0) & ~(Z % D == 0)
    integ = ~wall & (onint(X) | onint(Y))
    both = ~wall & onint(X) & onint(Y)
    half = ~wall & ~integ & (onhalf(X) | onhalf(Y))
    other = ~wall & ~integ & ~half
    tot = w.sum()
    return dict(total=float(tot), wall=float(w[wall].sum()), integer_lines=float(w[integ].sum()),
                integer_points=float(w[both].sum()), half_integer_lines=float(w[half].sum()), elsewhere=float(w[other].sum()),
                n_points=int((w > 0).sum()))


# ----------------------------------------------------------------------------- the LP loop
def run(a, log=print):
    import multiprocessing as mp
    s = float(a.s); t0 = time.time()
    m = Model(s, D=a.D, sym=not a.nosym)
    thetas = [math.radians(d) for d in (angle_list(a.deg_step) if not a.axis else [0.0])]
    thmax = 0.0 if a.axis else math.pi / 4
    pool = mp.get_context('fork').Pool(a.nproc) if a.nproc > 1 else None
    if a.from_cert:
        # refinement: columns are exactly the points of an existing certificate-format file (no lattice,
        # no literature, pricing off unless asked); initial rows = every lattice pose whose captured
        # weight under the file's own weights is < 1 + warm_thr (thinned per block)
        s0, P0, w0c = load_points(a.from_cert); assert abs(s0 - s) < 1e-12, "container mismatch"
        ncol = sum(m.add_point(x, y, 'cert') for x, y in P0)
        wmap = {(int(round(x * m.D)), int(round(y * m.D))): wt for (x, y), wt in zip(P0, w0c)}
        w0 = np.array([wmap.get((int(round(px * m.D)), int(round(py * m.D))), 0.0) for px, py in m.P])
        log(f"[{a.tag}] s={s} columns from {a.from_cert}: {ncol} orbits / {len(m.P)} atoms, file total {w0.sum():.6f}")
        gmin, nviol, poses, vals = separate(m.P, w0, s, thetas, a.pitch, thr=1.0 + a.warm_thr, perang=a.perang * 5, block=a.block / 2, pool=pool)
        log(f"[{a.tag}] file weights on the {a.pitch} lattice: min {gmin:.6f}, {nviol} poses below {1 + a.warm_thr}")
    else:
        ncol = build_columns(m, s, a.col_pitch, wall_points=not a.no_wall_points, literature=not a.no_literature)
        log(f"[{a.tag}] s={s} columns: {ncol} orbits / {len(m.P)} atoms (col pitch {a.col_pitch}, sym={m.sym}); "
            f"{len(thetas)} angles, row pitch {a.pitch}, nproc {a.nproc}")
        # initial rows: coarse lattice at every angle (nothing is covered yet, so take a thinned sample)
        w0 = np.zeros(len(m.P))
        gmin, nviol, poses, vals = separate(m.P, w0, s, thetas, a.pitch * 5, perang=a.perang, block=a.pitch * 5, pool=pool)
    m.add_rows(poses)
    log(f"[{a.tag}] initial rows: {len(m.rows)}")
    hist = []; best = None
    for it in range(a.max_iters):
        out = m.solve()
        if out is None: log("LP failed"); break
        val, x, y, Ax = out; w = x[m.own]
        # --- cutting planes: lattice scan
        gmin, nviol, poses, vals = separate(m.P, w, s, thetas, a.pitch, perang=a.perang, block=a.block, pool=pool)
        # --- local search around the worst lattice poses and the tightest rows
        seeds = []
        if len(vals): seeds.append(poses[np.argsort(vals)[:a.nseeds]])
        tight = np.argsort(Ax)[:a.nseeds]; seeds.append(np.array(m.rows)[tight])
        if it % 3 == 0:                                            # random restarts in the whole pose space
            rng = np.random.default_rng(it); th = rng.choice(thetas, size=a.nseeds)
            lo = (np.cos(th) + np.sin(th)) / 2 + WALL
            seeds.append(np.c_[lo + rng.random(a.nseeds) * (s - 2 * lo), lo + rng.random(a.nseeds) * (s - 2 * lo), th])
        seeds = np.concatenate(seeds)
        if pool is not None and len(seeds) > 8:
            chunks = np.array_split(seeds, a.nproc)
            outs = pool.map(_polish_worker, [(m.P, w, s, ch, a.polish_rounds, a.polish_n, it * 100 + i, thmax) for i, ch in enumerate(chunks)])
            fp = np.concatenate([o[0] for o in outs]); fv = np.concatenate([o[1] for o in outs])
            pmin = min(o[3].min() for o in outs if len(o[3]))
        else:
            fp, fv, cur, curv = polish(m.P, w, s, seeds, rounds=a.polish_rounds, nper=a.polish_n, thmax=thmax, arad0=(0.0 if thmax == 0 else math.radians(1.0)))
            pmin = curv.min() if len(curv) else 9e9
        # thin the polished poses: worst per (block/4 cell, angle bucket)
        if len(fv):
            o = np.argsort(fv); keep = []; seen = set()
            for j in o:
                k = (int(fp[j, 0] / (a.block / 4)), int(fp[j, 1] / (a.block / 4)), int(round(math.degrees(fp[j, 2]) * 20)))
                if k in seen: continue
                seen.add(k); keep.append(j)
                if len(keep) >= a.perang * 4: break
            fp = fp[keep]; fv = fv[keep]
        # --- column generation
        ncols = 0; Mcov = None
        if it >= a.cg_start and it % a.cg_every == 0 and not a.no_colgen:
            cand = price(m, y, pitch=a.cg_pitch, want=a.cg_want)
            Mcov = cand[0][0] if cand else 1.0
            for cv, X, Y in cand: ncols += m.add_point(X / m.D, Y / m.D, 'priced')
        n1 = m.add_rows(poses); n2 = m.add_rows(fp)
        rec = dict(it=it, t=round(time.time() - t0, 1), LP=val, rows=len(m.rows), orbits=len(m.orbits), atoms=len(m.P),
                   lattice_min=gmin, lattice_viol=nviol, polish_min=float(pmin), new_rows=n1 + n2, new_cols=ncols, dualcov=Mcov,
                   support=int((x > 1e-9).sum()))
        hist.append(rec)
        log(f"[{a.tag}] it{it} LP={val:.6f} rows={len(m.rows)} orb={len(m.orbits)} atoms={len(m.P)} support={rec['support']} "
            f"latmin={gmin:.6f} viol={nviol} polmin={pmin:.6f} +rows={n1}+{n2} +cols={ncols}"
            + (f" dualcov={Mcov:.4f}" if Mcov is not None else "") + f" t={time.time()-t0:.0f}s")
        if n1 + n2 == 0: best = (val, it)                       # LP value with no violated pose found
        json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/closed4_{a.tag}.json", 'w'), indent=1)
        # checkpoint every round (2026-08-30 coordinator note): so a kill/timeout mid-run never
        # loses the current LP weights -- NOT itself a certified cover (rows are a finite sample).
        # (column generation may have widened the model since `x` was solved: pad with zeros)
        xe = x if len(x) == len(m.orbits) else np.r_[x, np.zeros(len(m.orbits) - len(x))]
        export(m, xe, f"runs/closed4_{a.tag}_last.txt", WD=10 ** 7, up=True)
        if n1 + n2 == 0 and ncols == 0:
            log(f"[{a.tag}] converged (no violated pose found, no priced column)"); break
        if time.time() - t0 > a.time: log(f"[{a.tag}] time limit"); break
        if len(m.rows) > a.prune_at:
            keep = y > 1e-12; keep[-a.prune_keep:] = True; m.prune(keep); log(f"   pruned to {len(m.rows)} rows")
    # final solve on all rows (the last iteration added rows after solving)
    out = m.solve(); val, x, y, Ax = out
    w = x[m.own]
    log(f"[{a.tag}] final LP over all {len(m.rows)} rows: {val:.6f}")
    if pool is not None: pool.close(); pool.join()
    tw, npts = export(m, x, f"runs/closed4_{a.tag}_best.txt", WD=10 ** 7, up=True)
    bd = weight_breakdown(m.P, w, s, m.D)
    tags = {}
    for k in range(len(m.orbits)):
        if x[k] > 1e-9: tags[m.tag[k]] = tags.get(m.tag[k], 0.0) + x[k] * m.sizes[k]
    js = dict(s=s, LP=val, exported_total=tw, points=npts, breakdown=bd, weight_by_origin=tags, rows=len(m.rows), orbits=len(m.orbits),
              hist=hist, points_list=[(float(px), float(py), float(wp)) for (px, py), wp in zip(m.P, w) if wp > 1e-9])
    json.dump(js, open(f"runs/closed4_{a.tag}_best.json", 'w'), indent=1)
    log(f"[{a.tag}] exported runs/closed4_{a.tag}_best.txt: {npts} points, total {tw:.6f}; breakdown {bd}; by origin {tags}")
    try: plot_points(m.P, w, s, f"runs/closed4_{a.tag}_best.png", title=f"s={s} LP={val:.4f}")
    except Exception as e: log(f"   [plot failed: {e}]")
    return val, m, x


def plot_points(P, w, s, path, title=''):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    sel = w > 1e-9
    fig, ax = plt.subplots(figsize=(7, 7))
    for k in range(1, int(s)): ax.axhline(k, color='0.85', lw=0.6); ax.axvline(k, color='0.85', lw=0.6)
    ax.scatter(P[sel, 0], P[sel, 1], s=400 * w[sel], c='tab:blue', alpha=0.6, edgecolors='k', linewidths=0.3)
    ax.set_xlim(0, s); ax.set_ylim(0, s); ax.set_aspect('equal'); ax.set_title(title + f"  total={w.sum():.4f}  ({sel.sum()} points)")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


# ----------------------------------------------------------------------------- stress test
def _stress_worker(args):
    P, w, s, th, pitch = args
    vals, poses = scan_angle(P, w, s, th, pitch)
    if len(vals) == 0: return 9e9, None
    j = int(vals.argmin()); return float(vals[j]), poses[j]


def stress(P, w, s, pitch=0.005, nproc=8, nrand=200, seed=1, log=print, full_range=True):
    """dense scan: centres at pitch (rotated-frame lattice + wall bands), angles 0, 1e-4, 1e-3, ..., 45
    plus nrand random angles in [0, 90) (full range: does not rely on the D4 symmetry).  Returns
    (min, worst poses list)."""
    import multiprocessing as mp
    degs = [0.0, 1e-4, 1e-3, 1e-2, 1e-1] + [float(d) for d in range(1, 46)] + [45 - 1e-4, 45 - 1e-3, 45 - 1e-2, 45 - 1e-1]
    rng = np.random.default_rng(seed)
    rand = rng.random(nrand) * (90.0 if full_range else 45.0)
    thetas = [math.radians(d) for d in degs] + [math.radians(d) for d in rand]
    pool = mp.get_context('fork').Pool(nproc)
    outs = pool.map(_stress_worker, [(P, w, s, th, pitch) for th in thetas]); pool.close(); pool.join()
    res = sorted([(v, p) for v, p in outs if p is not None], key=lambda t: t[0])
    vmin = res[0][0]
    # polish the worst few to see how low a local search can push them
    seeds = np.array([p for v, p in res[:16]])
    fp, fv, cur, curv = polish(P, w, s, seeds, rounds=8, nper=64)
    pmin = float(curv.min())
    worst = [(float(v), [float(p[0]), float(p[1]), math.degrees(float(p[2]))]) for v, p in res[:10]]
    log(f"stress: {len(thetas)} angles, pitch {pitch}: lattice min = {vmin:.7f}; after local polish = {pmin:.7f}")
    for v, p in worst[:6]: log(f"   {v:.7f} at centre ({p[0]:.4f}, {p[1]:.4f}) angle {p[2]:.5f} deg")
    j = int(curv.argmin()); log(f"   polished worst: {pmin:.7f} at ({cur[j,0]:.5f}, {cur[j,1]:.5f}) angle {math.degrees(cur[j,2]):.6f} deg")
    return min(vmin, pmin), worst, (float(cur[j, 0]), float(cur[j, 1]), math.degrees(float(cur[j, 2])))


def load_points(path):
    sfr, D, WD, rows = read_cert(path)
    P = rows[:, :2] / D; w = rows[:, 2] / WD
    return float(sfr), P, w


# ----------------------------------------------------------------------------- sanity checks
def sanity(a, log=print):
    import multiprocessing as mp
    s = 4.0; pool = mp.get_context('fork').Pool(a.nproc)
    thetas = [math.radians(d) for d in angle_list(a.deg_step)]
    # (i) Friedman's 14 points, unit weights, NOT symmetrised: check all rows over the full angle range
    P = np.array(FRIEDMAN14, dtype=float); w = np.ones(len(P))
    th_full = sorted(set([t for t in thetas] + [math.pi / 2 - t for t in thetas]))     # x<->y images
    gmin, nviol, poses, vals = separate(P, w, s, th_full, a.pitch, perang=50, block=0.1, pool=pool)
    log(f"(i) Friedman 14, unit weights, lattice pitch {a.pitch}, {len(th_full)} angles in [0,90]: min captured = {gmin:.4f}, violated poses = {nviol}")
    # also the mirror image set (x<->y) and a polish from the worst places
    rng = np.random.default_rng(0); th = rng.choice(th_full, size=64)
    lo = (np.abs(np.cos(th)) + np.abs(np.sin(th))) / 2 + WALL
    seeds = np.c_[lo + rng.random(64) * (s - 2 * lo), lo + rng.random(64) * (s - 2 * lo), th]
    fp, fv, cur, curv = polish(P, w, s, seeds, rounds=8, nper=64, thmax=math.pi / 2)
    log(f"    local search from 64 random poses: min captured = {curv.min():.4f} (poses < 1 found: {len(fv)})")
    ok1 = gmin >= 1 - 1e-9 and curv.min() >= 1 - 1e-9
    # (ii) LP over the 14 points alone (unsymmetrised columns, rows over [0,90])
    m = Model(s, D=a.D, sym=False)
    for x, y in FRIEDMAN14: m.add_point(x, y, 'friedman')
    g0, nv0, poses0, v0 = separate(m.P, np.zeros(len(m.P)), s, th_full, a.pitch * 5, perang=100000, block=0.0, pool=pool)
    m.add_rows(poses0)
    val = None
    for it in range(12):
        out = m.solve(); val, x, y, Ax = out
        gmin, nviol, poses, vals = separate(m.P, x[m.own], s, th_full, a.pitch, perang=200, block=0.05, pool=pool)
        fp, fv, cur, curv = polish(m.P, x[m.own], s, poses[np.argsort(vals)[:32]] if len(vals) else np.zeros((0, 3)), rounds=6, nper=48, thmax=math.pi / 2)
        n = m.add_rows(poses) + m.add_rows(fp)
        log(f"(ii) LP over Friedman's 14 points: it{it} value {val:.6f} rows {len(m.rows)} latmin {gmin:.5f} viol {nviol} +rows {n}")
        if n == 0: break
    ok2 = val <= 14 + 1e-6
    log(f"(ii) LP over Friedman's 14 points alone = {val:.6f} (weights: {np.round(x, 4).tolist()})  {'OK (<= 14)' if ok2 else 'FAIL'}")
    # Bentz 16 and Nagamochi points (with the segments/area sampled by a 0.01 lattice on the lines) as a check of the transcription
    P = np.array(bentz16(4.0)); w = np.ones(len(P))
    gmin, nviol, poses, vals = separate(P, w, s, th_full, a.pitch, perang=50, block=0.1, pool=pool)
    fp, fv, cur, curv = polish(P, w, s, poses[np.argsort(vals)[:32]] if len(vals) else seeds, rounds=8, nper=64, thmax=math.pi / 2)
    log(f"(i') Bentz 16, unit weights: lattice min = {gmin:.4f}, violated = {nviol}, polished min = {curv.min():.4f}")
    pool.close(); pool.join()
    return ok1 and ok2


def nagamochi_set(s=4.0, seg_pitch=0.01):
    """Nagamochi's weighted set for [0,4]^2 discretised: points as given; each segment L_i as atoms of
    weight 0.5*seg_pitch at pitch seg_pitch; the area [1,3]^2 as atoms of weight seg_pitch^2 on a
    seg_pitch grid (cell centres).  Total 14 up to discretisation."""
    pts = [(x, y, wt) for x, y, wt in nagamochi_points(s)]
    n = int(round(2.2 / seg_pitch))
    for k in range(n + 1):
        t = 0.9 + k * seg_pitch; wt = 0.5 * seg_pitch * (0.5 if k in (0, n) else 1.0)
        pts += [(t, 1, wt), (t, 3, wt), (1, t, wt), (3, t, wt)]
    na = int(round(2 / seg_pitch))
    for i in range(na):
        for j in range(na):
            pts.append((1 + (i + 0.5) * seg_pitch, 1 + (j + 0.5) * seg_pitch, seg_pitch ** 2))
    return pts


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['run', 'sanity', 'stress', 'axis', 'nagamochi', 'lower'])
    ap.add_argument('cert', nargs='?', default=None, help='certificate file (stress)')
    ap.add_argument('--s', type=float, default=4.0, help='container side (multiple of 1/D)')
    ap.add_argument('--D', type=int, default=1000, help='coordinate denominator of the columns')
    ap.add_argument('--tag', default='closed4')
    ap.add_argument('--time', type=float, default=5400, help='wall-clock limit (s) for run')
    ap.add_argument('--max-iters', type=int, default=200)
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--pitch', type=float, default=0.02, help='row centre pitch')
    ap.add_argument('--col-pitch', type=float, default=0.05, help='column lattice pitch')
    ap.add_argument('--deg-step', type=float, default=1.0, help='angle step (degrees) for 1..45')
    ap.add_argument('--perang', type=int, default=150, help='max lattice cuts per angle per round')
    ap.add_argument('--block', type=float, default=0.1, help='thinning block for lattice cuts')
    ap.add_argument('--nseeds', type=int, default=48, help='local-search seeds per source per round')
    ap.add_argument('--polish-rounds', type=int, default=7); ap.add_argument('--polish-n', type=int, default=48)
    ap.add_argument('--cg-start', type=int, default=2); ap.add_argument('--cg-every', type=int, default=1)
    ap.add_argument('--cg-want', type=int, default=150); ap.add_argument('--cg-pitch', type=float, default=0.01)
    ap.add_argument('--no-colgen', action='store_true')
    ap.add_argument('--from-cert', default=None, help='run: columns = the points of this certificate file (refinement of a finished set)')
    ap.add_argument('--warm-thr', type=float, default=0.02, help='--from-cert: initial rows are lattice poses with captured weight < 1 + warm_thr')
    ap.add_argument('--no-wall-points', action='store_true'); ap.add_argument('--no-literature', action='store_true')
    ap.add_argument('--nosym', action='store_true', help='columns are single points, no D4 orbits (rows then need the full angle range: not implemented for run)')
    ap.add_argument('--axis', action='store_true', help='axis-aligned rows only (angle 0)')
    ap.add_argument('--prune-at', type=int, default=90000); ap.add_argument('--prune-keep', type=int, default=30000)
    ap.add_argument('--stress-pitch', type=float, default=0.005); ap.add_argument('--nrand', type=int, default=200)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--margin', type=float, default=None,
                    help='solve the STABLE (margin-mu) cover problem: a point counts only if it is at '
                         'distance >= mu inside every edge of the unit square (sets TOL = -mu).  A cover '
                         'found this way is a closed cover with uniform geometric margin mu, which the '
                         'exact box checker search/zeromargin.py can certify; see search/RUNG2.md.')
    ap.add_argument('--K', type=int, default=3000, help='lower: tightest placements added to the packing support')
    ap.add_argument('--seg-pitch', type=float, default=0.01, help='discretisation pitch of Nagamochi segments/area (nagamochi mode)')
    a = ap.parse_args()
    if a.margin is not None:
        global TOL
        TOL = -float(a.margin)
        print(f"MARGIN MODE: a point counts only if inside by >= {a.margin} (TOL = {TOL})")
    HIGHS['random_seed'] = a.seed
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/closed4_{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"closed4.py {' '.join(sys.argv[1:])}")
    if a.mode == 'sanity':
        ok = sanity(a, log); log(f"SANITY {'OK' if ok else 'FAIL'}")
    elif a.mode == 'run':
        val, m, x = run(a, log)
        s, P, w = load_points(f"runs/closed4_{a.tag}_best.txt")
        vmin, worst, wp = stress(P, w, s, pitch=a.stress_pitch, nproc=a.nproc, nrand=a.nrand, log=log)
        log(f"[{a.tag}] RESULT s={s} LP={val:.6f} exported total={w.sum():.6f} stress min={vmin:.7f} cost=total/min={w.sum()/vmin:.6f}")
        js = json.load(open(f"runs/closed4_{a.tag}_best.json")); js['stress'] = dict(min=vmin, worst=worst, polished_worst=wp, cost=w.sum() / vmin)
        json.dump(js, open(f"runs/closed4_{a.tag}_best.json", 'w'), indent=1)
    elif a.mode == 'axis':
        a.axis = True; a.no_colgen = True
        val, m, x = run(a, log)
        log(f"[{a.tag}] AXIS-ALIGNED LP = {val:.6f}")
    elif a.mode == 'stress':
        s, P, w = load_points(a.cert)
        vmin, worst, wp = stress(P, w, s, pitch=a.stress_pitch, nproc=a.nproc, nrand=a.nrand, log=log)
        log(f"STRESS {a.cert}: total={w.sum():.6f} min={vmin:.7f} cost={w.sum()/vmin:.6f}")
    elif a.mode == 'lower':
        # rigorous (float-checked) lower bound on nu_f(s) from the point set of CERT, via nu_f.lower_from_cert
        # (cover LP over the file's points -> dual support -> packing LP on the 0.01 grid -> max_coverage
        # certification by subdivision; see search/CEILING.md).  Output: runs/closed4_<tag>_lower.json
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import nu_f as NF
        best = NF.lower_from_cert(a.cert, f"closed4tmp_{a.tag}", h=0.01, log=log, nproc=min(a.nproc, 4), K=a.K)
        if os.path.exists(f"runs/lower_closed4tmp_{a.tag}.json"): os.replace(f"runs/lower_closed4tmp_{a.tag}.json", f"runs/closed4_{a.tag}_lower.json")
        log(f"LOWER {a.cert}: L={best['L']:.5f} (packing mass {best.get('mass')}, M={best.get('M')})")
    elif a.mode == 'nagamochi':
        pts = nagamochi_set(4.0, a.seg_pitch); P = np.array([(x, y) for x, y, _ in pts]); w = np.array([wt for _, _, wt in pts])
        log(f"Nagamochi set discretised at {a.seg_pitch}: {len(pts)} atoms, total {w.sum():.4f}")
        vmin, worst, wp = stress(P, w, 4.0, pitch=a.stress_pitch, nproc=a.nproc, nrand=a.nrand, log=log)
        log(f"NAGAMOCHI total={w.sum():.6f} min={vmin:.7f} cost={w.sum()/vmin:.6f}")


if __name__ == '__main__':
    main()
