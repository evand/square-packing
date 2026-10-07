#!/usr/bin/env python3
"""line_cover.py -- HEURISTIC mixed covers (points + uniform densities on axis-parallel segments) of [0,s]^2.

Task line-cover (A side), 2026-09-27; search/LINE_COVER.md.  Format: tasks/line-cover/FORMAT.md, reader
search/mixed_cover.py.  Nothing here is a proof: float LP over sampled rows, float scans (the exact word is B's
search/zm_mixed.py).

Columns (all D4 orbits):
  * points (as closed4.Model), from a lattice (--col-pitch), from --cols-from files, and priced in on the dual;
  * segments: the pieces [k/q, (k+1)/q] of the axis-parallel lines x = a / y = a for a in --lines (default: the
    interior grid lines 1..s-1), total mass spread uniformly by length.  A row's coefficient is the length fraction of
    the piece inside the CLOSED square (a piece on an edge of Q counts fully).
Float evaluation of a mixed cover on a lattice: points by closed4's cumulative-sum scan (close_shifted.scan_shift),
each line by its cumulative mass function F (piecewise linear) at the two ends of the chord of the line in Q.

Modes
  loop  TAG --s S [--q 50] [--cols-from C,..] [--warm FILE] [--lp-rounds N] [--hs-pitch 0.002] ...
        phase A (cutting planes on the 0.01 lattice + polish + point pricing) for --lp-rounds rounds, then
        phase B (the close_shifted recipe: every round MEASURE by the strict protocol -- hardscan at --hs-pitch +
        item-3 family --, plus a near-tile blow-up scan; permanent dip rows; shifted-lattice separation; re-solve)
        until the strict min is stable (< --stable-tol over --stable rounds).
  eval  FILE [--pitch 0.004] [--tiles]    strict protocol on a mixed (or plain certificate) file
  confirm FILE [--dips D] [--pitch 0.001]  pitch-0.001 confirmation (full container at the hardscan angles and the
        interleaved ones, near-tile boxes, boxes around the dips), then polish
  scale FILE OUT --factor F               multiply every mass by F (rounded up)
Outputs runs/lc_TAG.log / .json / _r{R}.txt (mixed format, masses rounded UP) / _dips.txt.
"""
import sys, os, math, time, json, argparse
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
import multiprocessing as mp
from fractions import Fraction
import closed4 as C
import family_rows as F
import close_shifted as CS
import mixed_cover as MC

TOL = C.TOL
NEAR = CS.NEAR
TILT = CS.TILT
EPS_K = 1e-14


# ============================================================================== float mixed cover
def _interval(k, A, B):
    """u with k*u in [A, B] (arrays): returns (lo, hi); empty -> lo > hi."""
    pos = k > EPS_K; neg = k < -EPS_K; zer = ~(pos | neg)
    ks = np.where(zer, 1.0, k)
    lo = np.where(pos, A / ks, np.where(neg, B / ks, np.where((A <= 0) & (B >= 0), -np.inf, np.inf)))
    hi = np.where(pos, B / ks, np.where(neg, A / ks, np.where((A <= 0) & (B >= 0), np.inf, -np.inf)))
    return lo, hi


def chord(orient, a, cx, cy, th, tol=None):
    """parameter interval (along the line) of the line x = a (orient 0; parameter y) or y = a (orient 1; parameter x)
    inside the closed unit square with centre (cx, cy), angle th (arrays broadcast).  Empty: lo > hi."""
    if tol is None: tol = TOL
    h = 0.5 + tol
    c = np.cos(th); s = np.sin(th)
    if orient == 0:
        d = a - cx; base = cy
        l1, h1 = _interval(s, -h - c * d, h - c * d)          # |c d + s u| <= h
        l2, h2 = _interval(c, -h + s * d, h + s * d)          # |-s d + c u| <= h
    else:
        d = a - cy; base = cx
        l1, h1 = _interval(c, -h - s * d, h - s * d)          # |c u + s d| <= h
        l2, h2 = _interval(s, -h + c * d, h + c * d)          # |-s u + c d| <= h
    return base + np.maximum(l1, l2), base + np.minimum(h1, h2)


class Cover:
    """float mixed cover: points P (n,2) with masses w, and axis-parallel segment pieces
    segs (m,4) = [orient, a, t0, t1] with masses sw (total mass of the piece, uniform by length)."""

    def __init__(self, s, P, w, segs=None, sw=None):
        self.s = float(s)
        self.P = np.asarray(P, float).reshape(-1, 2); self.w = np.asarray(w, float)
        segs = np.zeros((0, 4)) if segs is None else np.asarray(segs, float).reshape(-1, 4)
        sw = np.zeros(0) if sw is None else np.asarray(sw, float)
        nz = self.w > 0; self.P = self.P[nz]; self.w = self.w[nz]
        nz = sw > 0; self.segs = segs[nz]; self.sw = sw[nz]
        self.lines = []
        keys = sorted(set((int(o), round(float(a), 9)) for o, a in self.segs[:, :2])) if len(self.segs) else []
        for (o, a) in keys:
            sel = (self.segs[:, 0] == o) & (np.abs(self.segs[:, 1] - a) < 1e-9)
            t0 = self.segs[sel, 2]; t1 = self.segs[sel, 3]; ww = self.sw[sel]
            bk = np.unique(np.r_[t0, t1])
            cum = np.zeros(len(bk))
            for x0, x1, m_ in zip(t0, t1, ww):
                cum += m_ * np.clip((bk - x0) / (x1 - x0), 0.0, 1.0)
            self.lines.append((o, a, bk, cum))

    def total(self):
        return float(self.w.sum() + self.sw.sum())

    def lines_mass(self, cx, cy, th):
        cx = np.asarray(cx, float); out = np.zeros(np.broadcast(cx, np.asarray(cy), np.asarray(th)).shape)
        for (o, a, bk, cum) in self.lines:
            lo, hi = chord(o, a, cx, cy, th)
            ok = hi > lo
            v = np.interp(hi, bk, cum) - np.interp(lo, bk, cum)
            out += np.where(ok, v, 0.0)
        return out

    def mass(self, poses, chunk=200000):
        poses = np.asarray(poses, float).reshape(-1, 3); out = np.zeros(len(poses))
        for i in range(0, len(poses), chunk):
            B = poses[i:i + chunk]
            if len(self.w): out[i:i + chunk] += C.captured(self.P, self.w, B)
            if self.lines: out[i:i + chunk] += self.lines_mass(B[:, 0], B[:, 1], B[:, 2])
        return out

    def scan(self, th, pitch, sh=(0.0, 0.0), box=None):
        """close_shifted.scan_shift (points; rotated-frame lattice + wall bands) plus the lines at the same poses."""
        P, w = (self.P, self.w) if len(self.w) else (np.zeros((1, 2)) - 10.0, np.zeros(1))
        v, p = CS.scan_shift(P, w, self.s, th, pitch, sh, box)
        if len(v) and self.lines:
            for i in range(0, len(v), 500000):
                v[i:i + 500000] += self.lines_mass(p[i:i + 500000, 0], p[i:i + 500000, 1], th)
        return v, p


def cover_from_file(path):
    cv = MC.load(path); MC.validate(cv)
    D, W = cv['D'], cv['W']
    P = np.array([[x / D, y / D] for x, y, _ in cv['points']]).reshape(-1, 2)
    w = np.array([ww / W for *_, ww in cv['points']])
    segs, sw = [], []
    for (X0, Y0, X1, Y1, ww) in cv['segments']:
        if X0 == X1: segs.append((0, X0 / D, min(Y0, Y1) / D, max(Y0, Y1) / D))
        elif Y0 == Y1: segs.append((1, Y0 / D, min(X0, X1) / D, max(X0, X1) / D))
        else: raise ValueError('line_cover.Cover handles axis-parallel segments only')
        sw.append(ww / W)
    if cv['polygons']: raise ValueError('line_cover.Cover does not handle polygons')
    return Cover(float(cv['s']), P, w, segs, sw), cv


# ============================================================================== polish (closed4.polish on a Cover)
def polish(cov, seeds, rounds=9, nper=64, rng=None, radius0=0.03, arad0=math.radians(1.0), thmax=math.pi / 4):
    if rng is None: rng = np.random.default_rng(0)
    s = cov.s; WALL = C.WALL
    seeds = np.asarray(seeds, dtype=float).reshape(-1, 3)
    if len(seeds) == 0: return np.zeros((0, 3)), np.zeros(0), np.zeros((0, 3)), np.zeros(0)
    cur = seeds.copy(); curv = cov.mass(cur); found = {}

    def clip(Q):
        Q[:, 2] = np.clip(Q[:, 2], 0.0, thmax)
        wid = np.abs(np.cos(Q[:, 2])) + np.abs(np.sin(Q[:, 2])); lo = wid / 2 + WALL; hi = s - wid / 2 - WALL
        Q[:, 0] = np.clip(Q[:, 0], lo, hi); Q[:, 1] = np.clip(Q[:, 1], lo, hi)
        return Q
    specials = [0.0, math.pi / 4, 1e-6, 1e-5, 1e-4, 1e-3, math.pi / 4 - 1e-6, math.pi / 4 - 1e-4]
    for r in range(rounds):
        rad = radius0 * (0.3 ** r); arad = arad0 * (0.3 ** r); m = len(cur)
        Q = np.repeat(cur, nper, axis=0)
        d = rng.normal(size=(m * nper, 3)); d[:, :2] *= rad; d[:, 2] *= arad
        keep = rng.random(m * nper) < 0.3; d[keep, 2] = 0.0
        Q = Q + d
        spec = rng.random(m * nper) < 0.08
        Q[spec, 2] = rng.choice(specials, size=int(spec.sum()))
        Q = clip(Q); v = cov.mass(Q)
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


# ============================================================================== pool workers
def _worst30(args):
    cov, th, pitch, sh, box = args
    v, p = cov.scan(th, pitch, sh, box)
    if len(v) == 0: return np.zeros(0), np.zeros((0, 3))
    o = np.argsort(v)[:30]
    return v[o], p[o]


def _sep(args):
    cov, th, pitch, sh, perang, block, thr = args
    v, p = cov.scan(th, pitch, sh)
    if len(v) == 0: return np.zeros((0, 3)), np.zeros(0), 9e9
    bad = np.nonzero(v < thr)[0]
    if len(bad) == 0: return np.zeros((0, 3)), np.zeros(0), float(v.min())
    bad = bad[np.argsort(v[bad], kind='stable')]
    key = np.floor(p[bad, :2] / block).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
    _, first = np.unique(kk, return_index=True)
    sel = bad[np.sort(first)][:perang]
    return p[sel], v[sel], float(v.min())


def _pol(args):
    cov, seeds, seed = args
    return polish(cov, seeds, rounds=9, nper=64, rng=np.random.default_rng(seed))


def _mass(args):
    cov, Q = args
    return cov.mass(Q)


def tile_poses(s, pitch=0.0005, rad=0.015, degs=None):
    """near-tile family: every tile centre (i+1/2, j+1/2) +- rad on a pitch lattice, at small angles (degrees)."""
    if degs is None:
        degs = [0.0] + [math.degrees(10 ** e) for e in (-6, -5, -4.5, -4, -3.5, -3, -2.5)] + list(np.arange(0.01, 1.5 + 1e-9, 0.01))
    n = int(round(s)); g = np.arange(-rad, rad + 1e-12, pitch)
    out = []
    # D4 covers: the tiles with i <= j <= (n-1)/2 suffice (canonical octant), all others are images
    for i in range(n):
        for j in range(i, n):
            if j > (n - 1) / 2 + 1e-9: continue
            cx0, cy0 = i + 0.5, j + 0.5
            GX, GY = np.meshgrid(cx0 + g, cy0 + g, indexing='ij')
            for d in degs:
                th = math.radians(d)
                out.append(np.c_[GX.ravel(), GY.ravel(), np.full(GX.size, th)])
    Q = np.concatenate(out)
    return CS.clip_poses(Q, s)


def germ_poses(s, q, theta=1e-5, sub=1):
    """tile-germ blow-up family (LINE_COVER.md sec 2): c_x = i + 1/2 + theta * t on a half-integer line, c_y on the
    1/(sub q) grid, t on the 1/(sub q) grid in [-1/2, 1/2], at the small angle theta; and the transpose.  As theta -> 0
    the capture of a piecewise-constant line density is piecewise linear in (c_y, t) with vertices on the 1/q grid,
    so this finite family samples every vertex of the germ (up to O(theta) and the off-line points)."""
    n = int(round(s)); g = 1.0 / (q * sub)
    ts = np.arange(-0.5, 0.5 + 1e-12, g); w = (math.cos(theta) + math.sin(theta)) / 2
    free = np.arange(0, s + 1e-12, g); free = free[(free > w) & (free < s - w)]
    out = []
    for i in range(n):
        if i > (n - 1) / 2 + 1e-9: continue            # D4: half-integer lines up to the middle suffice
        T, V = np.meshgrid(ts, free, indexing='ij')
        cx = i + 0.5 + theta * T.ravel(); cy = V.ravel()
        out.append(np.c_[cx, cy, np.full(len(cx), theta)])
        out.append(np.c_[cy, cx, np.full(len(cx), theta)])
    # the 2-D blow-up at every tile centre (i + 1/2 + theta tx, j + 1/2 + theta ty): the vertical lines' chords depend
    # on tx, the horizontal lines' on ty (additively), so both have to sit at vertices at once
    T1, T2 = np.meshgrid(ts, ts, indexing='ij'); T1 = T1.ravel(); T2 = T2.ravel()
    for i in range(n):
        for j in range(i, n):
            if j > (n - 1) / 2 + 1e-9: continue
            out.append(np.c_[i + 0.5 + theta * T1, j + 0.5 + theta * T2, np.full(len(T1), theta)])
    return CS.clip_poses(np.concatenate(out), s)


def measure(cov, pool, pitch=0.002, nvisit=1000, nfam=300, family=True, tiles=True, tile_pitch=0.001, npar=28, germ_q=0):
    """strict protocol (hardscan at pitch + item-3 family) plus (tiles) the near-tile blow-up family."""
    t0 = time.time(); s = cov.s
    outs = pool.map(_worst30, [(cov, math.radians(d), pitch, (0.0, 0.0), None) for d in CS.hs_degs()])
    V = np.concatenate([o[0] for o in outs]); Q = np.concatenate([o[1] for o in outs])
    lat = float(V.min()); seeds = Q[np.argsort(V)[:1400]]
    res = pool.map(_pol, [(cov, ch, i) for i, ch in enumerate(np.array_split(seeds, 28))])
    fp = np.concatenate([r[0] for r in res]); fv = np.concatenate([r[1] for r in res])
    cur = np.concatenate([r[2] for r in res]); curv = np.concatenate([r[3] for r in res])
    j = int(curv.argmin()); vo = np.argsort(fv)[:nvisit]
    out = dict(lat=lat, pol=float(curv[j]), pose=cur[j].copy(), cur=cur, curv=curv, visit=fp[vo], visitv=fv[vo],
               t_hs=time.time() - t0)
    fmin, fpose, fbad, fbadv, fn = 9e9, np.zeros(3), np.zeros((0, 3)), np.zeros(0), 0
    if family:
        angles = [0.0] + [10 ** e for e in (-7, -6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5)]
        FQ = np.array(F.item3(s, 0.0005, angles))
        v = np.concatenate(pool.map(_mass, [(cov, ch) for ch in np.array_split(FQ, 2 * npar)]))
        k = int(v.argmin()); bad = np.nonzero(v < 1 - 1e-9)[0]; fn = len(bad); bad = bad[np.argsort(v[bad])[:nfam]]
        fmin, fpose, fbad, fbadv = float(v[k]), FQ[k], FQ[bad], v[bad]
    out.update(fam=fmin, fam_pose=fpose, fam_rows=fbad, fam_rowsv=fbadv, fam_n=fn)
    tmin, tpose, tbad, tbadv = 9e9, np.zeros(3), np.zeros((0, 3)), np.zeros(0)
    if tiles:
        TQ = tile_poses(s, tile_pitch)
        v = np.concatenate(pool.map(_mass, [(cov, ch) for ch in np.array_split(TQ, 2 * npar)]))
        k = int(v.argmin()); bad = np.nonzero(v < 1 - 1e-9)[0]; bad = bad[np.argsort(v[bad])[:nfam]]
        tmin, tpose, tbad, tbadv = float(v[k]), TQ[k], TQ[bad], v[bad]
    if germ_q:
        GQ = np.r_[germ_poses(s, germ_q, 1e-5, 2), germ_poses(s, germ_q, 1e-3, 1)]
        v = np.concatenate(pool.map(_mass, [(cov, ch) for ch in np.array_split(GQ, 2 * npar)]))
        k = int(v.argmin()); bad = np.nonzero(v < 1 - 1e-9)[0]; bad = bad[np.argsort(v[bad])[:nfam]]
        out.update(germ=float(v[k]), germ_pose=GQ[k])
        if float(v[k]) < tmin: tpose = GQ[k]
        tmin = min(tmin, float(v[k])); tbad = np.r_[tbad, GQ[bad]]; tbadv = np.r_[tbadv, v[bad]]
    out.update(tile=tmin, tile_pose=tpose, tile_rows=tbad, tile_rowsv=tbadv,
               hs=min(lat, out['pol']), strict=min(lat, out['pol'], fmin), all=min(lat, out['pol'], fmin, tmin),
               t_all=time.time() - t0)
    return out


# ============================================================================== LP model
def d4_pt(X, Y, K):
    return sorted(set([(X, Y), (K - X, Y), (X, K - Y), (K - X, K - Y), (Y, X), (K - Y, X), (Y, K - X), (K - Y, K - X)]))


def d4_seg(o, A, T0, T1, K):
    """images of the piece (orient o, line A, [T0, T1]) under D4 (integers)."""
    out = set()
    for (oo, aa, t0, t1) in [(o, A, T0, T1), (o, K - A, T0, T1), (o, A, K - T1, K - T0), (o, K - A, K - T1, K - T0)]:
        out.add((oo, aa, t0, t1)); out.add((1 - oo, aa, t0, t1))
    return sorted(out)


class MixedModel:
    """columns = D4 orbits of points or of segment pieces; rows = poses; A_rk = captured mass of column k
    (per unit column value = mass of each atom of the orbit) in the closed square of row r."""

    def __init__(self, s, D=1000):
        self.s = float(s); self.D = D; self.K = int(round(s * D)); self.sym = True
        self.sizes = []; self.kind = []; self.keys = {}
        self.P = np.zeros((0, 2)); self.own = np.zeros(0, dtype=np.int64)          # point atoms
        self.S = np.zeros((0, 4)); self.sown = np.zeros(0, dtype=np.int64)         # segment atoms [o, a, t0, t1] floats
        self.SI = []                                                               # segment atoms, integer tuples
        self.rows = []; self.rkey = set(); self.R = []; self.C = []; self.V = []
        self.hs = None
        self.ban = set()                     # canonical D4 keys of point orbits never used as columns (--ban)

    @property
    def orbits(self):                        # closed4.price uses len(m.orbits) only through m.rows / y
        return self.sizes

    # ---------------------------------------------------------------- columns
    def add_point(self, x, y, tag=''):
        X = min(max(int(round(x * self.D)), 0), self.K); Y = min(max(int(round(y * self.D)), 0), self.K)
        imgs = d4_pt(X, Y, self.K); key = ('p', imgs[0])
        if key in self.keys or imgs[0] in self.ban: return False
        k = len(self.sizes); self.keys[key] = k; self.sizes.append(float(len(imgs))); self.kind.append('p')
        of = np.array(imgs, dtype=float) / self.D
        self.P = np.concatenate([self.P, of]); self.own = np.concatenate([self.own, np.full(len(of), k)])
        if self.rows: self._col_block(np.arange(len(self.P) - len(of), len(self.P)), None)
        return True

    def add_seg(self, o, A, T0, T1):
        imgs = d4_seg(o, A, T0, T1, self.K); key = ('s', imgs[0])
        if key in self.keys: return False
        k = len(self.sizes); self.keys[key] = k; self.sizes.append(float(len(imgs))); self.kind.append('s')
        n0 = len(self.S)
        self.SI.extend(imgs)
        self.S = np.concatenate([self.S, np.array(imgs, dtype=float) * np.array([1.0, 1.0 / self.D, 1.0 / self.D, 1.0 / self.D])])
        self.sown = np.concatenate([self.sown, np.full(len(imgs), k)])
        if self.rows: self._col_block(None, np.arange(n0, len(self.S)))
        return True

    # ---------------------------------------------------------------- coefficients
    @staticmethod
    def _pt_coef(Q, P, own):
        """sparse (rows, atoms-> columns) via captured-style containment; returns COO triplets (r, col, count)."""
        R, Cc, V = [], [], []
        for i in range(0, len(Q), 256):
            B = Q[i:i + 256]; ct = np.cos(B[:, 2]); st = np.sin(B[:, 2])
            u0 = B[:, 0] * ct + B[:, 1] * st; u1 = -B[:, 0] * st + B[:, 1] * ct
            q0 = P[:, 0:1] * ct[None, :] + P[:, 1:2] * st[None, :]
            q1 = -P[:, 0:1] * st[None, :] + P[:, 1:2] * ct[None, :]
            ins = np.maximum(np.abs(q0 - u0[None, :]), np.abs(q1 - u1[None, :])) <= 0.5 + TOL
            a, r = np.nonzero(ins)
            if len(a):
                key = (r + i).astype(np.int64) * (1 << 31) + own[a]
                u, cnt = np.unique(key, return_counts=True)
                R.append(u // (1 << 31)); Cc.append(u % (1 << 31)); V.append(cnt.astype(float))
        if not R: return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
        return np.concatenate(R), np.concatenate(Cc), np.concatenate(V)

    @staticmethod
    def _seg_coef(Q, S, sown):
        """COO triplets (r, col, length fraction summed over the orbit's atoms)."""
        R, Cc, V = [], [], []
        if len(S) == 0 or len(Q) == 0: return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
        keys = sorted(set((int(o), round(float(a), 9)) for o, a in S[:, :2]))
        for (o, a) in keys:
            idx = np.nonzero((S[:, 0] == o) & (np.abs(S[:, 1] - a) < 1e-9))[0]
            t0 = S[idx, 2]; t1 = S[idx, 3]; L = t1 - t0
            for i in range(0, len(Q), 4096):
                B = Q[i:i + 4096]
                lo, hi = chord(o, a, B[:, 0], B[:, 1], B[:, 2])
                ok = np.nonzero(hi > lo)[0]
                if len(ok) == 0: continue
                ov = (np.minimum(hi[ok, None], t1[None, :]) - np.maximum(lo[ok, None], t0[None, :])) / L[None, :]
                rr, jj = np.nonzero(ov > 1e-12)
                if len(rr) == 0: continue
                key = (ok[rr] + i).astype(np.int64) * (1 << 31) + sown[idx[jj]]
                vals = np.minimum(ov[rr, jj], 1.0)
                u, inv = np.unique(key, return_inverse=True)
                R.append(u // (1 << 31)); Cc.append(u % (1 << 31)); V.append(np.bincount(inv, weights=vals))
        if not R: return np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0)
        return np.concatenate(R), np.concatenate(Cc), np.concatenate(V)

    def _col_block(self, pidx, sidx):
        Q = np.array(self.rows)
        if pidx is not None:
            r, c, v = self._pt_coef(Q, self.P[pidx], self.own[pidx])
        else:
            r, c, v = self._seg_coef(Q, self.S[sidx], self.sown[sidx])
        if len(r): self.R.append(r.astype(np.int32)); self.C.append(c.astype(np.int32)); self.V.append(v)
        self.hs = None

    def canon(self, cx, cy, th):
        s = self.s
        imgs = [(cx, cy), (s - cy, cx), (s - cx, s - cy), (cy, s - cx)]
        a, b = min((round(x, 7), round(y, 7)) for x, y in imgs)
        return (a, b, round(th, 9))

    def add_rows(self, poses):
        new = []
        for cx, cy, th in np.asarray(poses, dtype=float).reshape(-1, 3):
            key = self.canon(cx, cy, th)
            if key in self.rkey: continue
            self.rkey.add(key); new.append((float(cx), float(cy), float(th)))
        if not new: return 0
        r0 = len(self.rows); self.rows.extend(new); Q = np.array(new)
        for (r, c, v) in (self._pt_coef(Q, self.P, self.own) if len(self.P) else (np.zeros(0),) * 3,
                          self._seg_coef(Q, self.S, self.sown)):
            if len(r): self.R.append((r + r0).astype(np.int32)); self.C.append(c.astype(np.int32)); self.V.append(v)
        return len(new)

    def matrix(self):
        nr, nc = len(self.rows), len(self.sizes)
        if not self.R: return sp.csr_matrix((nr, nc))
        return sp.coo_matrix((np.concatenate(self.V), (np.concatenate(self.R), np.concatenate(self.C))), shape=(nr, nc)).tocsr()

    def prune(self, keep):
        idx = np.full(len(self.rows), -1, dtype=np.int64); idx[np.nonzero(keep)[0]] = np.arange(int(keep.sum()))
        R = np.concatenate(self.R); Cc = np.concatenate(self.C); V = np.concatenate(self.V)
        nr = idx[R]; sel = nr >= 0
        self.R = [nr[sel].astype(np.int32)]; self.C = [Cc[sel]]; self.V = [V[sel]]
        self.rows = [rw for rw, k in zip(self.rows, keep) if k]
        self.rkey = set(self.canon(*rw) for rw in self.rows); self.hs = None

    def solve(self, ipm=True):
        import highspy
        A = self.matrix().tocsr(); nr, nc = A.shape; inf = highspy.kHighsInf
        h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('random_seed', 0)
        if ipm:
            h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', 'off')
        cost = np.array(self.sizes, float)
        h.addVars(nc, np.zeros(nc), np.full(nc, inf)); h.changeColsCost(nc, np.arange(nc, dtype=np.int32), cost)
        A.sort_indices()
        h.addRows(nr, np.ones(nr), np.full(nr, inf), int(A.nnz), A.indptr[:-1].astype(np.int32), A.indices.astype(np.int32),
                  A.data.astype(float))
        h.run()
        if h.getModelStatus() != highspy.HighsModelStatus.kOptimal: return None
        sol = h.getSolution(); x = np.maximum(np.array(sol.col_value), 0.0); y = np.maximum(np.array(sol.row_dual), 0.0)
        if ipm: x[x < 1e-10] = 0.0; y[y < 1e-10] = 0.0
        return float(cost @ x), x, y, A @ x

    # ---------------------------------------------------------------- cover / export
    def cover(self, x):
        return Cover(self.s, self.P, x[self.own] if len(self.own) else np.zeros(0), self.S,
                     x[self.sown] if len(self.sown) else np.zeros(0))

    def export(self, x, path, W=10 ** 8, comment=None):
        pts, segs = [], []
        for (X, Y), k in zip(np.round(self.P * self.D).astype(np.int64), self.own):
            if x[k] > 1e-12:
                wq = int(math.ceil(x[k] * W - 1e-9))
                if wq > 0: pts.append((int(X), int(Y), wq))
        for (o, A, T0, T1), k in zip(self.SI, self.sown):
            if x[k] > 1e-12:
                wq = int(math.ceil(x[k] * W - 1e-9))
                if wq > 0: segs.append((A, T0, A, T1, wq) if o == 0 else (T0, A, T1, A, wq))
        cv = dict(s_num=self.K, s_den=self.D, D=self.D, W=W, points=pts, segments=segs, polygons=[])
        g = math.gcd(self.K, self.D); cv['s_num'] //= g; cv['s_den'] //= g
        MC.validate(cv); MC.write(path, cv, comment)
        return float(MC.total(cv)), len(pts), len(segs)


def scale_file(path, out, factor):
    cv = MC.load(path); f = Fraction(factor).limit_denominator(10 ** 9)
    up = lambda w: -((-w * f.numerator) // f.denominator)
    cv['points'] = [(x, y, up(w)) for x, y, w in cv['points']]
    cv['segments'] = [(a, b, c, d, up(w)) for a, b, c, d, w in cv['segments']]
    cv['polygons'] = [(up(w), vs) for w, vs in cv['polygons']]
    MC.validate(cv); MC.write(out, cv, f"scaled by {f} from {os.path.basename(path)}")
    return float(MC.total(cv))


# ============================================================================== the loop
def fmt(p):
    return f"({p[0]:.5f}, {p[1]:.5f}, {math.degrees(p[2]):.4f} deg)"


def build_model(a, log):
    s = float(a.s); m = MixedModel(s, D=1000); K = m.K; D = m.D
    for b in (a.ban or '').split(';'):       # e.g. --ban 1,1 (S20_LB.md: zmx2 cannot count a point at the corner germ's
        if b.strip():                        # second-order tangency (1,1); banning its orbit costs the LP almost nothing)
            bx, by = (float(v) for v in b.split(','))
            m.ban.add(d4_pt(int(round(bx * D)), int(round(by * D)), K)[0])
    if m.ban: log(f"banned point orbits: {sorted(m.ban)}")
    lines = [float(v) for v in a.lines.split(',')] if a.lines else list(range(1, int(round(s))))

    def online(X, Y):
        return any((X == int(round(l * D)) and 0 < Y < K) or (Y == int(round(l * D)) and 0 < X < K) for l in lines)
    npc = 0
    if a.col_pitch > 0:
        for x, y in C.column_lattice(s, a.col_pitch):
            if min(x, y) < 1e-9 or max(x, y) > s - 1e-9: continue
            X, Y = int(round(x * D)), int(round(y * D))
            if a.no_line_points and online(X, Y): continue
            npc += m.add_point(x, y, 'lattice')
    for f in (a.cols_from or '').split(','):
        if not f.strip(): continue
        if open(f.strip()).read(5) == 'mixed':
            cv = MC.load(f.strip()); Pp = [(x / cv['D'], y / cv['D']) for x, y, _ in cv['points']]
        else:
            Pp = C.load_points(f.strip())[1]
        for x, y in Pp:
            X, Y = int(round(x * D)), int(round(y * D))
            if a.no_line_points and online(X, Y): continue
            npc += m.add_point(x, y, 'cols')
    nsc = 0
    if a.q > 0:
        step = D // a.q; assert step * a.q == D, "q must divide D = 1000"
        for l in lines:
            A = int(round(l * D))
            for T in range(0, K, step):
                nsc += m.add_seg(0, A, T, T + step)
    log(f"columns: {npc} point orbits ({len(m.P)} atoms), {nsc} segment orbits ({len(m.S)} pieces, q = {a.q}, lines {lines}), "
        f"no_line_points={a.no_line_points}")
    return m


def loop(a):
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/lc_{a.tag}.log", 'a')

    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"line_cover.py loop {' '.join(sys.argv[2:])}   pid {os.getpid()}  {time.strftime('%F %T')}")
    rng = np.random.default_rng(a.seed); s = float(a.s)
    m = build_model(a, log)
    pool = mp.get_context('fork').Pool(a.nproc)
    prot = []

    def add(poses, protect):
        poses = np.asarray(poses, dtype=float).reshape(-1, 3)
        if protect and len(prot):
            idx = {m.canon(*rw): i for i, rw in enumerate(m.rows)}
            for q in poses:
                i = idx.get(m.canon(*q))
                if i is not None: prot[i] = True
        n = m.add_rows(poses); prot.extend([protect] * n); return n
    thetas = [math.radians(d) for d in C.angle_list(0.5)]
    # ---- initial rows
    if a.warm:
        cov0 = cover_from_file(a.warm)[0] if open(a.warm).read(5) == 'mixed' else Cover(s, *C.load_points(a.warm)[1:])
        log(f"warm start from {a.warm}: total {cov0.total():.6f}")
    else:
        cov0 = Cover(s, np.zeros((0, 2)), np.zeros(0))
    thr0 = 1.0 + a.init_thr
    outs = pool.map(_sep, [(cov0, th, 0.01 if a.warm else 0.05, (0.0, 0.0), 800 if a.warm else 150,
                            0.03 if a.warm else 0.05, thr0) for th in thetas])
    n0 = add(np.concatenate([o[0] for o in outs]), False)
    n0 += add(tile_poses(s, 0.004, 0.012, [0.0, 1e-4, 0.01, 0.03, 0.1, 0.3, 1.0]), False)
    if a.germ and a.q > 0:
        G = germ_poses(s, a.q, 1e-5, 1); ng = add(G, True); n0 += ng
        log(f"germ family (q = {a.q}, theta = 1e-5): {len(G)} poses, {ng} new permanent rows")
    for f in (a.rows_from or '').split(','):
        if f.strip():
            R = C.read_poses(f.strip(), s); n0 += add(R, True)
    log(f"initial rows: {n0} ({len(m.rows)})")
    t0 = time.time(); hist = []; r = 0; phase = 'A' if a.lp_rounds > 0 else 'B'; y = Ax = None; nA = 0
    lastLP = None
    while True:
        tl = time.time(); out = m.solve(ipm=not a.simplex); tlp = time.time() - tl
        if out is None: log("LP failed"); break
        val, x, y, Ax = out
        cov = m.cover(x); path = f"runs/lc_{a.tag}_r{r}.txt"
        tot, npt, nsg = m.export(x, path, comment=f"line_cover.py {a.tag} round {r} LP {val:.6f}")
        segw = float(cov.sw.sum())
        rec = dict(round=r, phase=phase, t=round(time.time() - t0), LP=val, total=tot, points=npt, segments=nsg, segw=segw,
                   rows=len(m.rows), perm=int(sum(prot)), cols=len(m.sizes), t_lp=round(tlp))
        if phase == 'A':
            # ---- cutting planes: 0.01 lattice, polish of the worst and the tight rows, point pricing
            outs = pool.map(_sep, [(cov, th, 0.01, tuple(rng.random(2)), 150, 0.05, 1 - 1e-9) for th in thetas])
            Lp = np.concatenate([o[0] for o in outs]); Lv = np.concatenate([o[1] for o in outs]); gmin = min(o[2] for o in outs)
            seeds = np.r_[Lp[np.argsort(Lv)[:400]], np.array(m.rows)[np.argsort(Ax)[:200]]]
            res = pool.map(_pol, [(cov, ch, 1000 * r + i) for i, ch in enumerate(np.array_split(seeds, 2 * a.nproc))])
            fp = np.concatenate([q[0] for q in res]); fv = np.concatenate([q[1] for q in res])
            pmin = min(float(q[3].min()) for q in res if len(q[3]))
            fp = fp[np.argsort(fv)[:3000]]
            tq = tile_poses(s, 0.002, 0.01, [0.0, 1e-5, 1e-3, 0.01, 0.03, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5])
            tv = cov.mass(tq); tb = tq[tv < 1 - 1e-9]; tb = tb[np.argsort(tv[tv < 1 - 1e-9])[:1000]]
            FQ = np.array(F.item3(s, 0.002, [0.0, 1e-6, 1e-4, 1e-3, 1e-2, 10 ** -1.5]))
            fqv = np.concatenate(pool.map(_mass, [(cov, ch) for ch in np.array_split(FQ, 2 * a.nproc)]))
            fb = FQ[fqv < 1 - 1e-9]; fb = fb[np.argsort(fqv[fqv < 1 - 1e-9])[:1000]]
            ncols = 0; dcov = None
            if not a.no_colgen:
                cand = C.price(m, y, pitch=0.01, want=a.cg_want)
                dcov = cand[0][0] if cand else 1.0
                for cv_, X, Y in cand:
                    if a.no_line_points and ((X % 1000 == 0 and 0 < Y < m.K) or (Y % 1000 == 0 and 0 < X < m.K)): continue
                    ncols += m.add_point(X / m.D, Y / m.D, 'priced')
            n1 = add(Lp, False) + add(fp, False) + add(tb, False) + add(fb, False)
            rec.update(latmin=gmin, polmin=pmin, tilemin=float(tv.min()), fammin=float(fqv.min()), newrows=n1, newcols=ncols, dualcov=dcov)
            log(f"[{a.tag}] A{r} LP={val:.6f} (segments {segw:.4f}) exported {tot:.6f} ({npt} pts, {nsg} segs) rows={len(m.rows)} "
                f"cols={len(m.sizes)} | lat {gmin:.5f} pol {pmin:.5f} tile {tv.min():.5f} fam {fqv.min():.5f} | +{n1} rows +{ncols} cols"
                + (f" dualcov {dcov:.4f}" if dcov else "") + f" | lp {tlp:.0f}s t={time.time() - t0:.0f}s")
            hist.append(rec); nA += 1
            if nA >= a.lp_rounds or (n1 == 0 and ncols == 0) or (lastLP is not None and abs(val - lastLP) < a.lp_tol * val and nA >= 4):
                phase = 'B'; log(f"[{a.tag}] -> phase B")
            lastLP = val
        else:
            ms = measure(cov, pool, pitch=a.hs_pitch, nvisit=a.nvisit, nfam=a.nfam, tiles=True, tile_pitch=a.tile_pitch,
                         germ_q=(a.q if a.germ else 0))
            smin = ms['strict']; worst = np.argsort(ms['curv'])
            rec.update(hs_lat=ms['lat'], hs_pol=ms['pol'], hs_pose=[float(ms['pose'][0]), float(ms['pose'][1]), math.degrees(ms['pose'][2])],
                       fam=ms['fam'], fam_n=ms['fam_n'], tile=ms['tile'],
                       tile_pose=[float(ms['tile_pose'][0]), float(ms['tile_pose'][1]), math.degrees(ms['tile_pose'][2])],
                       strict=smin, cost=tot / smin, cost_all=tot / ms['all'], t_meas=round(ms['t_all']),
                       worst=[(float(ms['curv'][k]), float(ms['cur'][k, 0]), float(ms['cur'][k, 1]), math.degrees(ms['cur'][k, 2]))
                              for k in worst[:8]])
            hist.append(rec)
            log(f"[{a.tag}] B{r} LP={val:.6f} (segments {segw:.4f}) total={tot:.6f} ({npt} pts, {nsg} segs) rows={len(m.rows)} "
                f"({rec['perm']} perm) | hs p{a.hs_pitch} lat {ms['lat']:.7f} pol {ms['pol']:.7f} at {fmt(ms['pose'])} | fam {ms['fam']:.7f} "
                f"| tiles {ms['tile']:.7f} at {fmt(ms['tile_pose'])} | strict {smin:.7f} cost {tot / smin:.4f} (with tiles {tot / ms['all']:.4f}) "
                f"| LP {tlp:.0f}s meas {ms['t_all']:.0f}s t={time.time() - t0:.0f}s")
            seen = [ms['pose']]
            for k in worst[1:]:
                q = ms['cur'][k]
                if all(np.abs(q - p_).max() > 0.01 for p_ in seen):
                    seen.append(q); log(f"      {ms['curv'][k]:.7f} at {fmt(q)}")
                if len(seen) >= 6: break
            hb = [h_ for h_ in hist if h_['phase'] == 'B']
            if len(hb) > a.stable:
                last = [h_['strict'] for h_ in hb[-(a.stable + 1):]]
                if all(abs(u - v) / v < a.stable_tol for u, v in zip(last[1:], last[:-1])):
                    log(f"[{a.tag}] STABLE: strict min over the last {a.stable + 1} rounds {', '.join(f'{v:.6f}' for v in last)}")
                    json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/lc_{a.tag}.json", 'w'), indent=1); break
            if smin >= 1 - 1e-7 and ms['all'] >= 1 - 1e-7: log(f"[{a.tag}] strict min >= 1"); break
            # ---- rows (close_shifted recipe)
            po = worst[:a.npol]; po = po[ms['curv'][po] < 1 - 1e-9]
            Dp = [ms['cur'][po], ms['visit'], ms['fam_rows'], ms['tile_rows']]
            Dv = [ms['curv'][po], ms['visitv'], ms['fam_rowsv'], ms['tile_rowsv']]
            sm = po[ms['curv'][po] < a.protect]
            st = CS.clip_poses(CS.stencil(ms['cur'][sm], ms['curv'][sm], a.nstencil, a.hs_pitch / 2, math.radians(0.01)), s)
            sh = rng.random(2)
            jit = [0.0] + list(NEAR + 0.02 * rng.random(len(NEAR))) + list(TILT[:-1] + 0.25 * rng.random(len(TILT) - 1)) + [45.0]
            outs = pool.map(_sep, [(cov, math.radians(d), a.hs_pitch, tuple(sh), a.perang, a.block, 1 - 1e-9) for d in jit])
            LP_ = np.concatenate([o[0] for o in outs]); LV = np.concatenate([o[1] for o in outs])
            lowm = np.concatenate([(np.arange(len(o[1])) < 4) & (o[1] < a.protect) for o in outs])
            Dp.append(LP_[lowm]); Dv.append(LV[lowm])
            Dp = np.concatenate(Dp).reshape(-1, 3); Dv = np.concatenate(Dv); pm = Dv < a.protect
            CS.write_poses(f"runs/lc_{a.tag}_dips.txt", Dp[pm], Dv[pm], note=f"round {r} strict {smin:.7f}")
            n1 = add(Dp[pm], True) + add(st, True)
            inj = f"runs/lc_{a.tag}_inject.txt"           # extra permanent rows dropped in by another process (S20_LB.md)
            if os.path.exists(inj):
                R_ = C.read_poses(inj, s); os.rename(inj, inj + f".used_r{r}")
                ni = add(np.asarray(R_).reshape(-1, 3), True); n1 += ni; log(f"   injected {ni} permanent rows from {inj}")
            n2 = add(Dp[~pm], False) + add(LP_[~lowm], False)
            outs = pool.map(_sep, [(cov, th, 0.01, (0.0, 0.0), 100, 0.08, 1 - 1e-7) for th in thetas])
            n4 = add(np.concatenate([o[0] for o in outs]), False)
            ncols = 0
            if a.colgen_b and not a.no_colgen:
                rows_all = m.rows; m.rows = rows_all[:len(y)]      # price against the rows y belongs to (rows were added above)
                try: cand = C.price(m, y, pitch=0.01, want=a.cg_want // 2)
                finally: m.rows = rows_all
                for cv_, X, Y in cand:
                    if a.no_line_points and ((X % 1000 == 0 and 0 < Y < m.K) or (Y % 1000 == 0 and 0 < X < m.K)): continue
                    ncols += m.add_point(X / m.D, Y / m.D, 'priced')
            log(f"   +rows {n1} permanent + {n2} other + {n4} 0.01-lattice; +{ncols} cols | rows {len(m.rows)} ({sum(prot)} perm)")
            if time.time() - t0 > a.time: log(f"[{a.tag}] time limit"); json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/lc_{a.tag}.json", 'w'), indent=1); break
            if len(hb) >= a.rounds: log(f"[{a.tag}] round limit"); json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/lc_{a.tag}.json", 'w'), indent=1); break
        json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/lc_{a.tag}.json", 'w'), indent=1)
        if a.prune_at and len(m.rows) > a.prune_at:
            keep = np.r_[(y > 1e-8) | (Ax <= 1 + a.prune_slack), np.zeros(len(m.rows) - len(y), dtype=bool)]
            keep[-a.prune_keep:] = True; keep |= np.array(prot, dtype=bool)
            m.prune(keep); prot = [p_ for p_, k in zip(prot, keep) if k]
            log(f"   pruned to {len(m.rows)} rows ({sum(prot)} permanent)")
        r += 1
    pool.close(); pool.join()
    hb = [h_ for h_ in hist if h_['phase'] == 'B']
    if hb:
        b = min(hb, key=lambda h_: h_['cost'])
        log(f"[{a.tag}] DONE; best strict cost {b['cost']:.4f} at round {b['round']} (strict {b['strict']:.7f}, total {b['total']:.6f}, "
            f"LP {b['LP']:.6f}); {time.strftime('%F %T')}")


# ============================================================================== eval / confirm
def eval_file(a):
    cov = cover_from_file(a.file)[0] if open(a.file).read(5) == 'mixed' else Cover(*C.load_points(a.file))
    pool = mp.get_context('fork').Pool(a.nproc)
    ms = measure(cov, pool, pitch=a.pitch, tiles=a.tiles, tile_pitch=a.tile_pitch, germ_q=a.germ_q)
    pool.close(); pool.join(); tot = cov.total()
    print(f"{a.file}: total {tot:.6f} ({len(cov.w)} points, {len(cov.sw)} segments, segment mass {cov.sw.sum():.4f}); "
          f"hardscan p{a.pitch}: lattice {ms['lat']:.7f}, polished {ms['pol']:.7f} at {fmt(ms['pose'])}; family {ms['fam']:.7f} "
          f"at {fmt(ms['fam_pose'])}; strict {ms['strict']:.7f} cost {tot / ms['strict']:.6f}"
          + (f"; tiles {ms['tile']:.7f} at {fmt(ms['tile_pose'])}, with tiles cost {tot / ms['all']:.6f}" if a.tiles else '')
          + f"  [{ms['t_all']:.0f}s]", flush=True)
    for k in np.argsort(ms['curv'])[:6]: print(f"   {ms['curv'][k]:.7f} at {fmt(ms['cur'][k])}")


def confirm(a):
    cov = cover_from_file(a.file)[0] if open(a.file).read(5) == 'mixed' else Cover(*C.load_points(a.file))
    s = cov.s; tot = cov.total(); t0 = time.time()
    pool = mp.get_context('fork').Pool(a.nproc)
    print(f"confirm {a.file}: total {tot:.6f}  {time.strftime('%F %T')}", flush=True)
    ms = measure(cov, pool, pitch=0.002, tiles=True, tile_pitch=0.0005, germ_q=a.germ_q)
    print(f"   hardscan p0.002 lat {ms['lat']:.7f} pol {ms['pol']:.7f} at {fmt(ms['pose'])}; family {ms['fam']:.7f}; "
          f"tiles p0.0005 {ms['tile']:.7f} at {fmt(ms['tile_pose'])} ({time.time() - t0:.0f}s)", flush=True)
    V = [np.array([ms['tile']])]; Q = [ms['tile_pose'][None]]
    Dd = [np.c_[ms['cur'], ms['curv']]]
    if a.dips and os.path.exists(a.dips): Dd.append(np.loadtxt(a.dips, comments='#').reshape(-1, 4))
    Dd = np.concatenate(Dd); Dd = Dd[np.argsort(Dd[:, 3])]; keep = []
    for d in Dd:
        if all(abs(d[0] - k[0]) > a.box or abs(d[1] - k[1]) > a.box or abs(math.degrees(d[2] - k[2])) > 0.3 for k in keep): keep.append(d)
        if len(keep) >= a.ndips: break
    jobs = []
    for d in keep:
        deg = math.degrees(d[2]); box = (d[0] - a.box, d[0] + a.box, d[1] - a.box, d[1] + a.box)
        angs = np.r_[0.0, np.arange(max(0.0, deg - 0.06), deg + 0.06 + 1e-12, 0.002)] if deg < 3.0 else \
            np.arange(deg - 0.3, min(45.0, deg + 0.3) + 1e-12, 0.01)
        jobs += [(cov, math.radians(t), a.pitch, (0.0, 0.0), box) for t in angs]
    outs = pool.map(_worst30, jobs, chunksize=4)
    Vb = np.concatenate([o[0] for o in outs]); Qb = np.concatenate([o[1] for o in outs])
    print(f"   {len(keep)} dip boxes (+-{a.box}), {len(jobs)} scans at pitch {a.pitch}: min {Vb.min():.7f} at {fmt(Qb[Vb.argmin()])} "
          f"({time.time() - t0:.0f}s)", flush=True)
    V.append(Vb); Q.append(Qb)
    for name, degs in (('hardscan', CS.hs_degs()), ('interleaved', list(np.arange(0.01, 3.0, 0.02)) + list(np.arange(3.125, 45.0, 0.25)))):
        outs = pool.map(_worst30, [(cov, math.radians(d), a.pitch, (0.0, 0.0), None) for d in degs])
        Vf = np.concatenate([o[0] for o in outs]); Qf = np.concatenate([o[1] for o in outs])
        print(f"   full container, {len(degs)} {name} angles, pitch {a.pitch}: min {Vf.min():.7f} at {fmt(Qf[Vf.argmin()])} "
              f"({time.time() - t0:.0f}s)", flush=True)
        V.append(Vf); Q.append(Qf)
    V = np.concatenate(V); Q = np.concatenate(Q)
    seeds = Q[np.argsort(V)[:1400]]
    res = pool.map(_pol, [(cov, ch, i) for i, ch in enumerate(np.array_split(seeds, 28))])
    pool.close(); pool.join()
    cur = np.concatenate([q[2] for q in res]); curv = np.concatenate([q[3] for q in res])
    j = int(curv.argmin()); mn = min(float(V.min()), float(curv[j]), ms['all'])
    print(f"{a.file}: confirm pitch {a.pitch}: scans min {V.min():.7f}; polished {curv[j]:.7f} at {fmt(cur[j])}; overall min "
          f"(incl. strict p0.002 + tiles) {mn:.7f}; cost total/min = {tot / mn:.6f}; {time.time() - t0:.0f}s", flush=True)
    for k in np.argsort(curv)[:8]: print(f"   {curv[k]:.7f} at {fmt(cur[k])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='mode', required=True)
    L = sub.add_parser('loop'); L.add_argument('tag'); L.add_argument('--s', type=float, required=True)
    L.add_argument('--q', type=int, default=50); L.add_argument('--lines', default=None)
    L.add_argument('--col-pitch', type=float, default=0.05); L.add_argument('--cols-from', default=None)
    L.add_argument('--no-line-points', action='store_true'); L.add_argument('--no-colgen', action='store_true')
    L.add_argument('--colgen-b', action='store_true'); L.add_argument('--cg-want', type=int, default=200)
    L.add_argument('--warm', default=None); L.add_argument('--rows-from', default=None); L.add_argument('--init-thr', type=float, default=0.005)
    L.add_argument('--lp-rounds', type=int, default=12); L.add_argument('--lp-tol', type=float, default=2e-4)
    L.add_argument('--hs-pitch', type=float, default=0.002); L.add_argument('--tile-pitch', type=float, default=0.001)
    L.add_argument('--perang', type=int, default=40); L.add_argument('--block', type=float, default=0.02)
    L.add_argument('--protect', type=float, default=0.995); L.add_argument('--nvisit', type=int, default=1000)
    L.add_argument('--npol', type=int, default=600); L.add_argument('--nfam', type=int, default=300); L.add_argument('--nstencil', type=int, default=40)
    L.add_argument('--prune-at', type=int, default=0); L.add_argument('--prune-keep', type=int, default=25000)
    L.add_argument('--prune-slack', type=float, default=0.003)
    L.add_argument('--rounds', type=int, default=40); L.add_argument('--time', type=float, default=36000)
    L.add_argument('--stable', type=int, default=3); L.add_argument('--stable-tol', type=float, default=0.001)
    L.add_argument('--germ', action='store_true', help='tile-germ blow-up family: permanent rows + measured each round')
    L.add_argument('--ban', default=None, help='point orbits never used as columns, "x,y;x,y"')
    L.add_argument('--simplex', action='store_true'); L.add_argument('--seed', type=int, default=1); L.add_argument('--nproc', type=int, default=20)
    E = sub.add_parser('eval'); E.add_argument('file'); E.add_argument('--pitch', type=float, default=0.004)
    E.add_argument('--germ-q', type=int, default=0)
    E.add_argument('--tiles', action='store_true'); E.add_argument('--tile-pitch', type=float, default=0.001); E.add_argument('--nproc', type=int, default=20)
    K = sub.add_parser('confirm'); K.add_argument('file'); K.add_argument('--dips', default=None)
    K.add_argument('--pitch', type=float, default=0.001); K.add_argument('--box', type=float, default=0.08)
    K.add_argument('--germ-q', type=int, default=50)
    K.add_argument('--ndips', type=int, default=80); K.add_argument('--nproc', type=int, default=20)
    S = sub.add_parser('scale'); S.add_argument('file'); S.add_argument('out'); S.add_argument('--factor', type=float, required=True)
    a = ap.parse_args()
    if a.mode == 'loop': loop(a)
    elif a.mode == 'eval': eval_file(a)
    elif a.mode == 'confirm': confirm(a)
    else: print(f"{a.out}: total {scale_file(a.file, a.out, a.factor):.9f}")


if __name__ == '__main__':
    main()
