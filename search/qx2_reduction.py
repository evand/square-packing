#!/usr/bin/env python3
"""qx2_reduction.py -- independent numeric double-check of the quadrant -> box reduction (HEURISTIC, float).

Given a quadrant solution (quadrant_lp `lp` mode, or qx2_lp `lp` mode), build the D4 box measure mu_k on [0,k]^2
ELEMENT BY ELEMENT with this file's own code (no quadrant_tools / BoxModel):
  * corner module nu = the model's corner elements (KA == 0; already the diagonal orbit) [+ Lebesgue [a,R]^2 for qx2];
  * profile pi = the model's base-cell atoms m.base, repeated at x = j + phase, clipped to [R, k-R] (closed), on the
    bottom wall [+ the Lebesgue layer y in [a,w] for qx2];
  * the box = the four rotations r^i (r(x,y) = (k-y, x)) of (nu + bottom band), plus Lebesgue on [w,k-w]^2 minus the
    four corner squares.
Then: total mass vs k^2 - 4D (D computed here from the solution), and a float oracle over the whole box with this
file's own evaluator (closed containment, tolerance TOL), cross-checked against quadrant_lp.Model.contrib.

usage: python3 qx2_reduction.py SOL.npz [--old] [--ks 6,7,...] [--nrand N] [--nproc P] [--xcheck N]
  --old: the npz has no 'args' (quadrant_lp lp mode); args --R --w --hc --hca --hp --hpa are used.
"""
import sys, os, math, json, argparse, time
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'): os.environ.setdefault(_v, '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import multiprocessing as mp

TOL = 1e-9       # closed containment tolerance (slight over-count, same convention as quadrant_lp)
WALL = 1e-9


# ============================================================================ load solution -> (nu, base profile)
def load(path, a):
    z = np.load(path)
    x = np.asarray(z['x'], float)
    if 'args' in z.files and not a.old:
        import qx2_lp
        args = argparse.Namespace(**json.loads(str(z['args'])))
        if args.mode != 'lp': raise SystemExit(f"{path}: mode {args.mode} is not a quadrant solution")
        m = qx2_lp.build(args); m.use_margin = False
        kind = 'qx2'; aa = m.a
    else:
        import quadrant_lp
        args = argparse.Namespace(mode='lp', R=a.R, w=a.w, hc=a.hc, hca=a.hca, hp=a.hp, hpa=a.hpa)
        m = quadrant_lp.build(args); kind = 'old'; aa = m.w
    assert len(x) == m.nvar, (len(x), m.nvar)
    R, w = m.R, m.w

    def mass(var, fm):
        var = int(var); return x[var] if var >= 0 else fm

    nu = dict(pts=[], hs=[], vs=[], rects=[])        # rects: (x0,x1,y0,y1,density)
    for arr_k, dst in (('pts', 'pts'), ('hs', 'hs'), ('vs', 'vs'), ('cells', 'rects')):
        A = m.A[arr_k]; K = m.KA[arr_k]
        for row, kk in zip(A, K):
            if kk != 0: continue
            ms = mass(row[-2], row[-1])
            if ms == 0: continue
            if dst == 'pts': nu['pts'].append((row[0], row[1], ms))
            elif dst == 'hs': nu['hs'].append((row[0], row[1], row[2], ms))       # y, x0, x1
            elif dst == 'vs': nu['vs'].append((row[0], row[1], row[2], ms))       # x, y0, y1
            else: nu['rects'].append((row[0], row[1], row[2], row[3], ms / ((row[1] - row[0]) * (row[3] - row[2]))))
    if kind == 'qx2' and aa < R: nu['rects'].append((aa, float(R), aa, float(R), 1.0))
    # profile base cell (phase in [0,1)), full mirror orbit
    base = []                                  # (typ, p0, p1, y0, y1, mass)
    for (typ, p0, p1, y0, y1, var, fm) in m.base:
        ms = mass(var, fm)
        if ms != 0: base.append((typ, p0, p1, y0, y1, ms))
    if kind == 'qx2' and aa < w: base.append(('c', 0.0, 1.0, aa, w, (w - aa)))
    MC = sum(t[2] for t in nu['pts']) + sum(t[3] for t in nu['hs']) + sum(t[3] for t in nu['vs']) + \
        sum(t[4] * (t[1] - t[0]) * (t[3] - t[2]) for t in nu['rects'])
    mv = sum(t[5] for t in base if t[0] in ('p', 'v') and abs(t[1]) < 1e-12)
    per = sum(t[5] for t in base)
    D = R * R - MC - mv
    info = dict(kind=kind, R=R, w=w, a=aa, MC=MC, mv=mv, sigma=w - per, D=D,
                D_npz=float(z['D']) if 'D' in z.files else float('nan'), D_model=float(m.D(x)))
    # symmetry checks of the input
    def key(t): return tuple(round(float(v), 9) for v in t)
    S = set(key(t) for t in nu['pts']); info['nu_pts_diag'] = all(key((b, a_, ms)) in S for (a_, b, ms) in nu['pts'])
    H = set(key(t) for t in nu['hs']); V = set(key(t) for t in nu['vs'])
    info['nu_seg_diag'] = H == V
    Rc = set(key(t) for t in nu['rects']); info['nu_rect_diag'] = all(key((c, d, a_, b, dd)) in Rc for (a_, b, c, d, dd) in nu['rects'])
    Bs = set((t[0],) + key(t[1:]) for t in base)
    def mir(t):
        typ, p0, p1, y0, y1, ms = t
        if typ in ('p', 'v'): q = (-p0) % 1.0; return (typ, q, q, y0, y1, ms)
        return (typ, 1.0 - p1, 1.0 - p0, y0, y1, ms)
    info['pi_mirror'] = all((mir(t)[0],) + key(mir(t)[1:]) in Bs for t in base)
    return m, x, nu, base, info


# ============================================================================ box construction
def build_box(nu, base, R, w, k):
    pts, hs, vs, rects = [], [], [], []
    # bottom band: x = j + phase, clipped to the closed interval [R, k-R]
    lo, hi = float(R), float(k - R); e = 1e-12
    band = dict(pts=[], hs=[], vs=[], rects=[])
    for (typ, p0, p1, y0, y1, ms) in base:
        for j in range(R - 1, k - R + 1):
            if typ == 'p':
                X = j + p0
                if lo - e <= X <= hi + e: band['pts'].append((X, y0, ms))
            elif typ == 'v':
                X = j + p0
                if lo - e <= X <= hi + e: band['vs'].append((X, y0, y1, ms))
            else:
                X0, X1 = max(j + p0, lo), min(j + p1, hi)
                if X1 - X0 <= e: continue
                fr = (X1 - X0) / (p1 - p0)
                if typ == 'h': band['hs'].append((y0, X0, X1, ms * fr))
                else: band['rects'].append((X0, X1, y0, y1, ms / ((p1 - p0) * (y1 - y0))))
    unit = dict(pts=nu['pts'] + band['pts'], hs=nu['hs'] + band['hs'], vs=nu['vs'] + band['vs'],
                rects=nu['rects'] + band['rects'])
    P = np.array(unit['pts'], float).reshape(-1, 3); Hs = np.array(unit['hs'], float).reshape(-1, 4)
    Vs = np.array(unit['vs'], float).reshape(-1, 4); Rc = np.array(unit['rects'], float).reshape(-1, 5)
    out = dict(pts=[], hs=[], vs=[], rects=[])
    for i in range(4):
        out['pts'].append(P.copy()); out['hs'].append(Hs.copy()); out['vs'].append(Vs.copy()); out['rects'].append(Rc.copy())
        # rotate by r(x,y) = (k-y, x)
        P = np.c_[k - P[:, 1], P[:, 0], P[:, 2]]
        nH = np.c_[Vs[:, 0], k - Vs[:, 2], k - Vs[:, 1], Vs[:, 3]]       # vseg x=X,[y0,y1] -> hseg y=X, x in [k-y1,k-y0]
        nV = np.c_[k - Hs[:, 0], Hs[:, 1], Hs[:, 2], Hs[:, 3]]           # hseg y=Y,[x0,x1] -> vseg x=k-Y, y in [x0,x1]
        Hs, Vs = nH, nV
        Rc = np.c_[k - Rc[:, 3], k - Rc[:, 2], Rc[:, 0], Rc[:, 1], Rc[:, 4]]
    box = {kk: np.concatenate(v) for kk, v in out.items()}
    # Lebesgue: [w,k-w]^2 minus the corner squares
    if w >= R - 1e-12: leb = [(w, k - w, w, k - w, 1.0)]
    else: leb = [(R, k - R, w, k - w, 1.0), (w, R, R, k - R, 1.0), (k - R, k - w, R, k - R, 1.0)]
    box['rects'] = np.r_[box['rects'], np.array(leb, float)]
    return box


def total(box):
    return box['pts'][:, 2].sum() + box['hs'][:, 3].sum() + box['vs'][:, 3].sum() + \
        (box['rects'][:, 4] * (box['rects'][:, 1] - box['rects'][:, 0]) * (box['rects'][:, 3] - box['rects'][:, 2])).sum()


# ============================================================================ evaluator (own code)
def _lin_interval(alpha, beta, h):
    """{t : |alpha t + beta| <= h} as (lo, hi); alpha, beta arrays."""
    big = 1e300
    nz = np.abs(alpha) > 1e-15
    al = np.where(nz, alpha, 1.0)
    t1 = (-h - beta) / al; t2 = (h - beta) / al
    lo = np.where(nz, np.minimum(t1, t2), np.where(np.abs(beta) <= h, -big, big))
    hi = np.where(nz, np.maximum(t1, t2), np.where(np.abs(beta) <= h, big, -big))
    return lo, hi


def seg_frac_h(Y, X0, X1, cx, cy, c, s, tol=TOL):
    """fraction of the segment y=Y, x in [X0,X1] inside the closed square (cx,cy,angle with cos c, sin s)."""
    h = 0.5 + tol; dy = Y - cy
    # u = (x-cx) c + dy s ; v = -(x-cx) s + dy c ; parameter t = x - cx
    l1, h1 = _lin_interval(c, dy * s, h); l2, h2 = _lin_interval(-s, dy * c, h)
    lo = cx + np.maximum(l1, l2); hi = cx + np.minimum(h1, h2)
    L = np.minimum(hi, X1) - np.maximum(lo, X0)
    return np.clip(L / (X1 - X0), 0.0, 1.0)


def seg_frac_v(X, Y0, Y1, cx, cy, c, s, tol=TOL):
    h = 0.5 + tol; dx = X - cx
    # u = dx c + (y-cy) s ; v = -dx s + (y-cy) c ; t = y - cy
    l1, h1 = _lin_interval(s, dx * c, h); l2, h2 = _lin_interval(c, -dx * s, h)
    lo = cy + np.maximum(l1, l2); hi = cy + np.minimum(h1, h2)
    L = np.minimum(hi, Y1) - np.maximum(lo, Y0)
    return np.clip(L / (Y1 - Y0), 0.0, 1.0)


def _clip_half(V, cnt, nx, ny, cc):
    """clip polygons V (n,8,2) with cnt vertices by nx*x+ny*y <= cc (scalars or (n,))."""
    n = len(V); M = V.shape[1]
    out = np.zeros_like(V); oc = np.zeros(n, int)
    f = V[:, :, 0] * nx + V[:, :, 1] * ny - np.asarray(cc).reshape(-1, 1)
    for i in range(M):
        act = i < cnt
        j = np.where(i + 1 < cnt, i + 1, 0)
        Pi = V[:, i]; Pj = V[np.arange(n), j]
        fi = f[:, i]; fj = f[np.arange(n), j]
        ini = act & (fi <= 0)
        idx = np.nonzero(ini)[0]
        out[idx, oc[idx]] = Pi[idx]; oc[idx] += 1
        crs = act & ((fi <= 0) != (fj <= 0))
        idx = np.nonzero(crs)[0]
        t = fi[idx] / (fi[idx] - fj[idx])
        out[idx, oc[idx]] = Pi[idx] + t[:, None] * (Pj[idx] - Pi[idx]); oc[idx] += 1
    return out, oc


def rect_area(cx, cy, c, s, x0, x1, y0, y1):
    n = len(cx)
    V = np.zeros((n, 8, 2))
    for i, (a, b) in enumerate(((-1, -1), (1, -1), (1, 1), (-1, 1))):
        V[:, i, 0] = cx + 0.5 * (a * c - b * s); V[:, i, 1] = cy + 0.5 * (a * s + b * c)
    cnt = np.full(n, 4)
    V, cnt = _clip_half(V, cnt, 1.0, 0.0, x1)
    V, cnt = _clip_half(V, cnt, -1.0, 0.0, -x0)
    V, cnt = _clip_half(V, cnt, 0.0, 1.0, y1)
    V, cnt = _clip_half(V, cnt, 0.0, -1.0, -y0)
    ar = np.zeros(n)
    for i in range(8):
        j = np.where(i + 1 < cnt, i + 1, 0)
        act = i < cnt
        a = V[:, i]; b = V[np.arange(n), j]
        ar += np.where(act, a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1], 0.0)
    return 0.5 * np.abs(ar)


class Evaluator:
    BIN = 0.5
    REACH = 0.7072

    def __init__(self, box):
        self.box = box
        b = box
        self.bb = dict(pts=np.c_[b['pts'][:, 0], b['pts'][:, 0], b['pts'][:, 1], b['pts'][:, 1]],
                       hs=np.c_[b['hs'][:, 1], b['hs'][:, 2], b['hs'][:, 0], b['hs'][:, 0]],
                       vs=np.c_[b['vs'][:, 0], b['vs'][:, 0], b['vs'][:, 1], b['vs'][:, 2]],
                       rects=b['rects'][:, :4])
        self.cache = {}

    def _cand(self, key):
        if key in self.cache: return self.cache[key]
        i, j = key; bx0 = i * self.BIN - self.REACH; bx1 = (i + 1) * self.BIN + self.REACH
        by0 = j * self.BIN - self.REACH; by1 = (j + 1) * self.BIN + self.REACH
        res = {}
        for kk, bb in self.bb.items():
            sel = (bb[:, 0] <= bx1) & (bb[:, 1] >= bx0) & (bb[:, 2] <= by1) & (bb[:, 3] >= by0)
            res[kk] = self.box[kk][sel]
        self.cache[key] = res
        return res

    def __call__(self, P):
        P = np.asarray(P, float); n = len(P); out = np.zeros(n)
        keys = np.floor(P[:, :2] / self.BIN).astype(int)
        kk = keys[:, 0] * 100003 + keys[:, 1]
        order = np.argsort(kk, kind='stable'); kks = kk[order]
        starts = np.r_[0, np.nonzero(np.diff(kks))[0] + 1, n]
        for a0, a1 in zip(starts[:-1], starts[1:]):
            idxall = order[a0:a1]
            key = tuple(keys[idxall[0]])
            E = self._cand(key)
            for c0 in range(0, len(idxall), 256):
                idx = idxall[c0:c0 + 256]
                out[idx] = self._eval(P[idx], E)
        return out

    @staticmethod
    def _eval(B, E):
        cx = B[:, 0:1]; cy = B[:, 1:2]; c = np.cos(B[:, 2:3]); s = np.sin(B[:, 2:3])
        tot = np.zeros(len(B))
        p = E['pts']
        if len(p):
            dx = p[None, :, 0] - cx; dy = p[None, :, 1] - cy
            u = dx * c + dy * s; v = -dx * s + dy * c
            ins = (np.abs(u) <= 0.5 + TOL) & (np.abs(v) <= 0.5 + TOL)
            tot += (ins * p[None, :, 2]).sum(1)
        hs = E['hs']
        if len(hs):
            fr = seg_frac_h(hs[None, :, 0], hs[None, :, 1], hs[None, :, 2], cx, cy, c, s)
            tot += (fr * hs[None, :, 3]).sum(1)
        vs = E['vs']
        if len(vs):
            fr = seg_frac_v(vs[None, :, 0], vs[None, :, 1], vs[None, :, 2], cx, cy, c, s)
            tot += (fr * vs[None, :, 3]).sum(1)
        rc = E['rects']
        if len(rc):
            nb, ne = len(B), len(rc)
            I = np.repeat(np.arange(nb), ne); J = np.tile(np.arange(ne), nb)
            # cheap bbox reject
            hw = 0.7072
            ok = (rc[J, 0] <= B[I, 0] + hw) & (rc[J, 1] >= B[I, 0] - hw) & (rc[J, 2] <= B[I, 1] + hw) & (rc[J, 3] >= B[I, 1] - hw)
            I = I[ok]; J = J[ok]
            if len(I):
                ar = rect_area(B[I, 0], B[I, 1], np.cos(B[I, 2]), np.sin(B[I, 2]), rc[J, 0], rc[J, 1], rc[J, 2], rc[J, 3])
                np.add.at(tot, I, ar * rc[J, 4])
        return tot


# ============================================================================ poses in the box
def hwid(th): return (np.abs(np.cos(th)) + np.abs(np.sin(th))) / 2


def clipbox(P, k):
    hw = hwid(P[:, 2]) + WALL
    P[:, 0] = np.clip(P[:, 0], hw, k - hw); P[:, 1] = np.clip(P[:, 1], hw, k - hw)
    P[:, 2] = np.mod(P[:, 2], math.pi / 2)
    return P


def random_poses(k, n, rng):
    P = np.c_[rng.uniform(0, k, n), rng.uniform(0, k, n), rng.uniform(0, math.pi / 2, n)]
    r = rng.random(n); sp = r < 0.25
    P[sp, 2] = rng.choice([0.0, 1e-6, 1e-4, 1e-2, math.pi / 4, math.pi / 2 - 1e-6], size=int(sp.sum()))
    P = clipbox(P, k); hw = hwid(P[:, 2]) + WALL
    r = rng.random(n); wall = rng.integers(0, 4, n)
    for wi in range(4):
        sel = (r < 0.3) & (wall == wi)
        if wi == 0: P[sel, 1] = hw[sel]
        elif wi == 1: P[sel, 1] = k - hw[sel]
        elif wi == 2: P[sel, 0] = hw[sel]
        else: P[sel, 0] = k - hw[sel]
    cor = rng.random(n) < 0.1; cx = rng.integers(0, 2, n); cy = rng.integers(0, 2, n)
    P[cor, 0] = np.where(cx[cor] == 0, hw[cor], k - hw[cor]); P[cor, 1] = np.where(cy[cor] == 0, hw[cor], k - hw[cor])
    return P


def near_lattice_poses(k, pitches=(0.1, 0.25), d=1e-7, ths=(0.0, 1e-6, -1e-6)):
    g = set()
    for h in pitches:
        for i in range(int(round((k - 1) / h)) + 1): g.add(round(i * h, 9))
    g = np.array(sorted(g))
    out = []
    for sx in (-1, 1):
        for sy in (-1, 1):
            X, Y = np.meshgrid(g + 0.5 + sx * d, g + 0.5 + sy * d, indexing='ij')
            for th in ths: out.append(np.c_[X.ravel(), Y.ravel(), np.full(X.size, th)])
    P = np.concatenate(out)
    hw = hwid(P[:, 2])
    keep = (P[:, 0] >= hw - 1e-12) & (P[:, 0] <= k - hw + 1e-12) & (P[:, 1] >= hw - 1e-12) & (P[:, 1] <= k - hw + 1e-12)
    P = P[keep]; P[:, 2] = np.mod(P[:, 2], math.pi / 2)
    return clipbox(P, k)


_W = {}


def _winit(box): _W['ev'] = Evaluator(box)


def _weval(P): return _W['ev'](P)


def _wpolish(args):
    P0, seed, rounds, k = args
    ev = _W['ev']; rng = np.random.default_rng(seed)
    cur = P0.copy(); curv = ev(cur); nper = 24
    for r in range(rounds):
        rad = 0.05 * 0.4 ** r; arad = math.radians(2.0) * 0.4 ** r
        Q = np.repeat(cur, nper, 0); d = rng.normal(size=Q.shape); d[:, :2] *= rad; d[:, 2] *= arad
        d[rng.random(len(Q)) < 0.3, 2] = 0
        Q = clipbox(Q + d, k)
        hw = hwid(Q[:, 2]) + WALL
        sel = rng.random(len(Q)) < 0.1; Q[sel, 1] = np.where(Q[sel, 1] < k / 2, hw[sel], k - hw[sel])
        sel = rng.random(len(Q)) < 0.1; Q[sel, 0] = np.where(Q[sel, 0] < k / 2, hw[sel], k - hw[sel])
        v = ev(Q).reshape(len(cur), nper); j = v.argmin(1)
        bv = v[np.arange(len(cur)), j]; better = bv < curv
        cur[better] = Q.reshape(len(cur), nper, 3)[np.arange(len(cur)), j][better]; curv[better] = bv[better]
    return cur, curv


def oracle(pool, k, rng, nrand, npol, nproc):
    P = np.r_[random_poses(k, nrand, rng), near_lattice_poses(k)]
    ch = np.array_split(P, nproc * 8)
    v = np.concatenate(pool.map(_weval, ch))
    o = np.argsort(v)
    key = np.floor(P[o, :2] / 0.05).astype(np.int64); kk = key[:, 0] * 100000 + key[:, 1]
    _, first = np.unique(kk, return_index=True)
    seeds = P[o[np.sort(first)][:npol]]
    res = pool.map(_wpolish, [(s, int(rng.integers(1 << 30)), 8, k) for s in np.array_split(seeds, nproc * 2)])
    cur = np.concatenate([r_[0] for r_ in res]); curv = np.concatenate([r_[1] for r_ in res])
    allp = np.r_[P, cur]; allv = np.r_[v, curv]; j = int(allv.argmin())
    nl = len(near_lattice_poses(k)); vl = v[nrand:]
    return dict(n=len(P), min=float(allv[j]), worst=allp[j].tolist(), min_random=float(v[:nrand].min()),
                min_lattice=float(vl.min()) if len(vl) else float('nan'), min_polish=float(curv.min()), nlat=nl)


# ============================================================================ cross-check vs quadrant_lp.Model.contrib
def xcheck(m, x, box, R, n, rng):
    ext = R + 3.0
    P = np.c_[rng.uniform(0, ext, n), rng.uniform(0, ext, n), rng.uniform(0, math.pi / 2, n)]
    sp = rng.random(n) < 0.3; P[sp, 2] = rng.choice([0.0, 1e-6, math.pi / 4], size=int(sp.sum()))
    hw = hwid(P[:, 2]) + 1e-7
    P[:, 0] = np.maximum(P[:, 0], hw); P[:, 1] = np.maximum(P[:, 1], hw)
    r = rng.random(n); P[r < 0.2, 1] = hw[r < 0.2]; P[r > 0.9, 0] = hw[r > 0.9]
    # near-lattice axis poses too
    g = np.round(rng.integers(0, int(ext / 0.05), n // 4) * 0.05, 9); g2 = np.round(rng.integers(0, int(ext / 0.05), n // 4) * 0.05, 9)
    PL = np.c_[g + 0.5 + rng.choice([-1e-7, 1e-7], n // 4), g2 + 0.5 + rng.choice([-1e-7, 1e-7], n // 4), np.zeros(n // 4)]
    PL[:, :2] = np.maximum(PL[:, :2], 0.5 + 1e-7)
    P = np.r_[P, PL]
    mine = Evaluator(box)(P)
    theirs = m.contrib(P, xval=x)
    d = np.abs(mine - theirs)
    j = int(d.argmax())
    return dict(n=len(P), maxdiff=float(d.max()), n_gt_1e6=int((d > 1e-6).sum()), worst=P[j].tolist(),
                mine=float(mine[j]), theirs=float(theirs[j]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sol'); ap.add_argument('--old', action='store_true')
    ap.add_argument('--R', type=int, default=2); ap.add_argument('--w', type=float, default=2.0)
    ap.add_argument('--hc', type=float, default=0.1); ap.add_argument('--hca', type=float, default=0.25)
    ap.add_argument('--hp', type=float, default=0.1); ap.add_argument('--hpa', type=float, default=0.25)
    ap.add_argument('--ks', default='6,7,8,9,10,11,12')
    ap.add_argument('--nrand', type=int, default=400000); ap.add_argument('--npol', type=int, default=3000)
    ap.add_argument('--nproc', type=int, default=4); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--xcheck', type=int, default=4000)
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    m, x, nu, base, info = load(a.sol, a)
    R, w = info['R'], info['w']
    print(f"== {a.sol}\n   model={info['kind']} R={R} w={w} a={info['a']}  nu(C)={info['MC']:.9f}  m_v={info['mv']:.9f}  "
          f"sigma={info['sigma']:.3e}  D={info['D']:.9f} (npz {info['D_npz']:.9f}, model {info['D_model']:.9f})  4D={4*info['D']:.9f}")
    print(f"   symmetry: nu points diag {info['nu_pts_diag']}, nu segs diag {info['nu_seg_diag']}, nu rects diag "
          f"{info['nu_rect_diag']}, pi mirror {info['pi_mirror']};  #nu atoms {sum(len(v) for v in nu.values())}, #base atoms {len(base)}", flush=True)
    res = dict(sol=a.sol, info=info, ks={})
    if a.xcheck:
        box = build_box(nu, base, R, w, 12)
        xc = xcheck(m, x, box, R, a.xcheck, rng); res['xcheck'] = xc
        print(f"   xcheck vs quadrant_lp.Model.contrib (box k=12, poses in [0,{R+3}]^2): n={xc['n']} max|diff|={xc['maxdiff']:.3e} "
              f"#>1e-6={xc['n_gt_1e6']} worst {xc['worst']} mine {xc['mine']:.9f} theirs {xc['theirs']:.9f}", flush=True)
    print(f"   {'k':>3} {'#elts':>7} {'total':>14} {'k^2-total':>12} {'4D':>12} {'diff':>10} {'#poses':>9} {'oracle min':>12}  worst pose (cx,cy,deg)")
    for k in [int(t) for t in a.ks.split(',')]:
        t0 = time.time()
        box = build_box(nu, base, R, w, k)
        tot = total(box); sav = k * k - tot; ne = sum(len(v) for v in box.values())
        pool = mp.Pool(a.nproc, initializer=_winit, initargs=(box,))
        o = oracle(pool, k, rng, a.nrand, a.npol, a.nproc); pool.terminate()
        wp = o['worst']
        print(f"   {k:>3} {ne:>7} {tot:>14.9f} {sav:>12.9f} {4*info['D']:>12.9f} {sav-4*info['D']:>10.2e} {o['n']:>9} {o['min']:>12.9f}  "
              f"({wp[0]:.6f},{wp[1]:.6f},{math.degrees(wp[2]):.5f})  [rand {o['min_random']:.7f} lat {o['min_lattice']:.7f} "
              f"pol {o['min_polish']:.7f}; {time.time()-t0:.0f}s]", flush=True)
        res['ks'][k] = dict(total=tot, saving=sav, fourD=4 * info['D'], **o)
    if a.out: json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
