"""Fill: add squares to a partial packing (10-20; one module for every constructor and operator).

A fill = placement (where the new squares go; overlap allowed) + optional relax (a named fq schedule,
pk/schedules.py).  Evan (10-20): exact fills undershoot in every useful case, so the defaults are inexact and random,
and a fill may end at a count other than its target (n +- 1 or 2): every result is a start for the n it reached.

Placers (registered with explicit parameter spaces, as pk/moves.py: unspecified parameters are drawn, every value
used comes back in info['params'] so outcomes can be regressed on them):

  bestofK   noisy greedy best fit (moved from shove/fill_greedy.py, CS1 agent A, 10-16): per square, raster the covered
            region, FFT-score overlap - LAM * contact per palette angle, draw K random poses, slide each to its local
            score minimum, keep the best, SAT-push it out of overlaps.  K = 1 random, large K near greedy (CS1: K = 4
            worst, K >= 16 similar; fill 1-3 s, the quench dominates, so large K is free).
  leastbad  least total penetration of K random poses (moved from shove/fill_relax.py, CS1 agent B): positions from
            uncovered raster points (biased to `region`), angles 0 / a neighbour's tilt / uniform.
  exact     fillfast exact (largest set of non-overlapping upright squares among fillfast's candidates; uprights only).
  greedy    fillfast greedy with overlap tolerance tau (overcounts by design; the quench repairs it).
  farthest  axis squares at the sampled points farthest from every centre (bulbgen's first top-up; a cheap baseline).
  uniform   centres U(0.5, s - 0.5)^2, angles U(0, 90) (pk.moves insert_k 'uniform').

Stopping (every placer): m = a count (place exactly m, or fewer if the placer runs out: exact / greedy), or
m = None with stop='full' (place until the placer's own "no room" test fires or m_max is reached).

  from pk.fill import fill
  X, info = fill(S, sq, random.Random(0), placer='bestofK', m=12, relax='round.15c', K=64)
  info: params, added, count, pen (total SAT penetration before relax), s (after relax; None = infeasible / no relax)

rng is a random.Random everywhere (as pk/moves.py).  The engines are the CS1 code verbatim (same draws, same output
for the same seed: checked against the 10-16 files on 10-20); shove/fill_greedy.py and fill_relax.py are adapters.
"""
from __future__ import annotations
import math, random, tempfile
import numpy as np


# ================================================================== geometry (moved from shove/harness.py)

def half_extent(a):
    r = math.radians(a % 90.0)
    return 0.5 * (abs(math.cos(r)) + abs(math.sin(r)))


def _inside(px, py, x, y, a):
    r = math.radians(a); c, s_ = math.cos(r), math.sin(r)
    u = (px - x) * c + (py - y) * s_; v = -(px - x) * s_ + (py - y) * c
    return (np.abs(u) <= 0.5) & (np.abs(v) <= 0.5)


def free_mask(S, sq, h=0.05):
    """(xs, ys, free) : free[j, i] = point (xs[i], ys[j]) lies in no square (raster of the uncovered region)"""
    xs = np.arange(h / 2, S, h); ys = np.arange(h / 2, S, h)
    X, Y = np.meshgrid(xs, ys)
    free = np.ones(X.shape, bool)
    for x, y, a in np.asarray(sq, float):
        i0, i1 = max(0, int((x - 0.75) / h)), min(len(xs), int((x + 0.75) / h) + 2)
        j0, j1 = max(0, int((y - 0.75) / h)), min(len(ys), int((y + 0.75) / h) + 2)
        free[j0:j1, i0:i1] &= ~_inside(X[j0:j1, i0:i1], Y[j0:j1, i0:i1], x, y, a)
    return xs, ys, free


def pen_vec(S, X, pose):
    """(depths vs each row of X, wall excess) for unit square `pose`; depth 0 = separated"""
    x, y, a = pose
    if len(X):
        dx, dy = X[:, 0] - x, X[:, 1] - y
        near = dx * dx + dy * dy < 2.0
        dep = np.zeros(len(X))
        if near.any():
            dxn, dyn, an = dx[near], dy[near], np.radians(X[near, 2])
            a0 = math.radians(a)
            best = np.full(len(dxn), np.inf)
            for phi in (np.full(len(dxn), a0), np.full(len(dxn), a0 + math.pi / 2), an, an + math.pi / 2):
                ux, uy = np.cos(phi), np.sin(phi)
                h1 = 0.5 * (np.abs(np.cos(a0 - phi)) + np.abs(np.sin(a0 - phi)))
                h2 = 0.5 * (np.abs(np.cos(an - phi)) + np.abs(np.sin(an - phi)))
                best = np.minimum(best, h1 + h2 - np.abs(dxn * ux + dyn * uy))
            dep[near] = np.maximum(best, 0.0)
    else:
        dep = np.zeros(0)
    h = half_extent(a)
    wall = max(0.0, h - x) + max(0.0, h - y) + max(0.0, x + h - S) + max(0.0, y + h - S)
    return dep, wall


def total_pen(S, X, pose):
    dep, wall = pen_vec(S, X, pose)
    return float(dep.sum()) + wall


def packing_pen(S, X):
    """total SAT penetration of a whole packing (pairs + walls): the 'how infeasible' number every placer reports"""
    X = np.asarray(X, float)
    return float(sum(total_pen(S, X[i + 1:], X[i]) for i in range(len(X))))     # each pair once, each wall once


# ================================================================== bestofK (moved from shove/fill_greedy.py)

RING = 0.08
PAD = 1.0
BESTOFK_DEFAULTS = dict(K=64, LAM=0.05, SNAP=1, SNAP_R=0.3, h=0.04, TILT=1, PUSH=1, RANDANG=0)


def _kernel(a, h, side):
    """raster footprint of a square of side `side` at angle a (deg), centred, step h"""
    R = side * 0.75
    n = int(math.ceil(R / h))
    t = (np.arange(-n, n + 1)) * h
    X, Y = np.meshgrid(t, t)
    r = math.radians(a); c, s = math.cos(r), math.sin(r)
    u = X * c + Y * s; v = -X * s + Y * c
    return ((np.abs(u) <= side / 2) & (np.abs(v) <= side / 2)).astype(float)


class Grid:
    def __init__(self, S, sq, h):
        self.S, self.h = S, h
        self.np_ = int(round(PAD / h))
        N = int(math.ceil(S / h)) + 2 * self.np_
        self.xs = (np.arange(N) - self.np_ + 0.5) * h          # pixel centres (x == y lattice)
        X, Y = np.meshgrid(self.xs, self.xs)
        self.X, self.Y = X, Y
        self.cov = ~((X > 0) & (X < S) & (Y > 0) & (Y < S))     # outside the box counts as covered
        for q in np.asarray(sq, float):
            self.paint(q)

    def paint(self, q):
        x, y, a = q
        h = self.h
        i0 = max(0, int((x - 0.75) / h) + self.np_ - 1); i1 = min(len(self.xs), int((x + 0.75) / h) + self.np_ + 2)
        j0 = max(0, int((y - 0.75) / h) + self.np_ - 1); j1 = min(len(self.xs), int((y + 0.75) / h) + self.np_ + 2)
        Xs, Ys = self.X[j0:j1, i0:i1], self.Y[j0:j1, i0:i1]
        r = math.radians(a); c, s = math.cos(r), math.sin(r)
        u = (Xs - x) * c + (Ys - y) * s; v = -(Xs - x) * s + (Ys - y) * c
        self.cov[j0:j1, i0:i1] |= (np.abs(u) <= 0.5) & (np.abs(v) <= 0.5)

    def score(self, a, lam):
        h = self.h
        K1 = _kernel(a, h, 1.0 - 1.5 * h)                     # shrunk: touching squares give no overlap pixels
        K2 = _kernel(a, h, 1.0 + 2 * RING)
        ring = np.clip(K2 - np.pad(K1, ((K2.shape[0] - K1.shape[0]) // 2,) * 2), 0, 1)
        C = self.cov.astype(float)
        ov = fftconvolve(C, K1, mode='same')
        ct = fftconvolve(C, ring, mode='same')
        ov = np.maximum(ov, 0); ov[ov < 0.5] = 0
        return ov * h * h, ov * h * h - lam * ct * h * h      # (overlap area, score)


def fftconvolve(*a, **k):
    from scipy.signal import fftconvolve as f
    return f(*a, **k)


def palette(S, sq, tilt):
    """angles to try: 0 plus up to 5 clustered signed tilts of the squares within 1.6 of the uncovered region"""
    if not tilt:
        return [0.0]
    sq = np.asarray(sq, float)
    xs, ys, free = free_mask(S, sq, h=0.1)
    J, I = np.nonzero(free)
    if not len(J):
        return [0.0]
    P = np.stack([xs[I], ys[J]], 1)
    near = []
    for x, y, a in sq:
        if np.min(np.hypot(P[:, 0] - x, P[:, 1] - y)) < 1.6:
            r = a % 90.0
            r = min(r, 90 - r) * (1 if (a % 90.0) <= 45 else -1)   # signed tilt in (-45, 45]
            if abs(r) > 1.0:
                near.append(r)
    angs = [0.0]
    for r in sorted(near, key=lambda r: -sum(abs(r - t) < 1.5 for t in near)):
        if all(abs(r - t) > 3.0 for t in angs) and len(angs) < 6:
            angs.append(r)
    return angs


def _corners(x, y, a):
    r = math.radians(a); c, s = math.cos(r), math.sin(r)
    off = np.array([[.5, .5], [-.5, .5], [-.5, -.5], [.5, -.5]])
    return np.array([x, y]) + off @ np.array([[c, s], [-s, c]]), (np.array([c, s]), np.array([-s, c]))


def _sep(P, Q, axP, axQ, cP, cQ):
    """SAT: (depth, unit vector pushing P away from Q); depth <= 0 means disjoint"""
    best = (math.inf, None)
    for ax in axP + axQ:
        a, b = P @ ax, Q @ ax
        ov = min(a.max() - b.min(), b.max() - a.min())
        if ov <= 0:
            return 0.0, None
        if ov < best[0]:
            sgn = 1.0 if (cP - cQ) @ ax >= 0 else -1.0
            best = (ov, sgn * ax)
    return best


def _push(S, sq, pose, iters=30):
    """continuous clean-up of a pixel-snapped pose: move it along the SAT separation vectors of its overlaps
    (neighbours and walls) until it fits or stops improving; returns the pose with the least max penetration"""
    x, y, a = pose
    he = half_extent(a)
    near = sq[np.hypot(sq[:, 0] - x, sq[:, 1] - y) < 2.2]
    NB = [_corners(*q) for q in near]
    best = None
    for _ in range(iters):
        P, axP = _corners(x, y, a)
        mv = np.zeros(2); worst = 0.0
        for (Q, axQ), q in zip(NB, near):
            d, v = _sep(P, Q, axP, axQ, np.array([x, y]), q[:2])
            if d > 0:
                mv += d * v; worst = max(worst, d)
        for dw, v in ((he - x, (1, 0)), (he - y, (0, 1)), (x + he - S, (-1, 0)), (y + he - S, (0, -1))):
            if dw > 0:
                mv += dw * np.array(v); worst = max(worst, dw)
        if best is None or worst < best[0]:
            best = (worst, (x, y, a))
        if worst < 1e-12:
            break
        x, y = x + 0.6 * mv[0], y + 0.6 * mv[1]
    return best[1]


def _descend(sc, j, i, w, iters=12):
    """slide (j, i) to a local minimum of sc by repeated argmin over a (2w+1)^2 window"""
    for _ in range(iters):
        j0, j1 = max(0, j - w), min(sc.shape[0], j + w + 1)
        i0, i1 = max(0, i - w), min(sc.shape[1], i + w + 1)
        win = sc[j0:j1, i0:i1]
        dj, di = np.unravel_index(np.argmin(win), win.shape)
        nj, ni = j0 + dj, i0 + di
        if (nj, ni) == (j, i) or win[dj, di] >= sc[j, i] - 1e-12:
            break
        j, i = nj, ni
    return j, i


def _region_mask(G, region, grow=0.6):
    """pixels of Grid G inside the union of rectangles `region` grown by `grow`"""
    m = np.zeros(G.cov.shape, bool)
    for x0, y0, x1, y1 in region:
        m |= (G.X > x0 - grow) & (G.X < x1 + grow) & (G.Y > y0 - grow) & (G.Y < y1 + grow)
    return m


def bestofk_candidate(S, sq, m, rng, kn, angs, region=None, p_region=0.75, ov_full=None):
    """the CS1 noisy greedy best fit: add m squares (m None: up to kn['m_max'], stopping when the best of K poses
    overlaps more than ov_full in area).  region: rectangles to draw centres from with probability p_region (None:
    every uncovered pixel, the CS1 behaviour).  -> (array of sq + added, number added)"""
    G = Grid(S, sq, kn['h'])
    add = []
    snapw = int(round(kn['SNAP_R'] / kn['h']))
    rmask = _region_mask(G, region) if region else None
    total = m if m is not None else kn.get('m_max', 10 ** 6)
    for _ in range(total):
        maps = {a: G.score(a, kn['LAM']) for a in angs}
        free = np.argwhere(~G.cov)
        if not len(free):
            if m is None:
                break
            free = np.argwhere(np.ones_like(G.cov))
        freeR = np.argwhere(~G.cov & rmask) if rmask is not None else None
        best = None
        for _k in range(max(1, kn['K'])):
            if freeR is not None and len(freeR) and rng.random() < p_region:
                j, i = freeR[rng.randrange(len(freeR))]
            else:
                j, i = free[rng.randrange(len(free))]
            a = angs[rng.randrange(len(angs))]
            sc = maps[a][1]
            if kn['SNAP']:
                j, i = _descend(sc, j, i, snapw)
            if best is None or sc[j, i] < best[0]:
                best = (sc[j, i], j, i, a)
        _, j, i, a = best
        if m is None and ov_full is not None and maps[a][0][j, i] > ov_full:
            break                                          # 'full': even the best of K overlaps too much
        pose = (float(G.xs[i]), float(G.xs[j]), float(a % 90.0))
        if kn['PUSH']:
            cur = np.vstack([np.asarray(sq, float)] + ([np.array(add)] if add else []))
            pose = tuple(map(float, _push(S, cur, pose)))
        add.append(pose)
        G.paint(pose)
    return np.vstack([np.asarray(sq, float)] + ([np.array(add)] if add else [])), len(add)


# ================================================================== leastbad (moved from shove/fill_relax.py)

def _sample_pts(S, X, strip, h=0.1):
    xs, ys, free = free_mask(S, X, h=h)
    J, I = np.nonzero(free)
    px, py = xs[I], ys[J]
    ins = np.zeros(len(px), bool)
    for x0, y0, x1, y1 in strip:
        ins |= (px > x0 - 0.6) & (px < x1 + 0.6) & (py > y0 - 0.6) & (py < y1 + 0.6)
    return np.c_[px, py], np.c_[px[ins], py[ins]]


def _angle(X, x, y, rng):
    u = rng.random()
    if u < 0.5:
        return 0.0
    if u < 0.85 and len(X):
        d = np.hypot(X[:, 0] - x, X[:, 1] - y)
        nb = np.nonzero((d < 2.5))[0]
        tilts = [X[j, 2] % 90 for j in nb if 1.0 < X[j, 2] % 90 < 89.0]
        if tilts:
            return (rng.choice(tilts) + rng.gauss(0, 2.0)) % 90
    return rng.uniform(0, 90)


def leastbad_place(S, X, strip, m, rng, K, resample_every=4, p_region=0.75, pen_full=None, m_max=10 ** 6):
    """add m squares to X one at a time, each the least-total-penetration of K random poses -> (X', total pen, added).
    strip: rectangles the centres are drawn from with probability p_region (CS1: the opened strip; [] = anywhere
    uncovered).  m None: until the best pose's penetration exceeds pen_full (or m_max squares)."""
    X = np.asarray(X, float)
    tot = 0.0
    allp, stp = _sample_pts(S, X, strip)
    added = 0
    for i in range(m if m is not None else m_max):
        if i and i % resample_every == 0:
            allp, stp = _sample_pts(S, X, strip)
        best = None
        for _ in range(K):
            if len(stp) and (rng.random() < p_region or not len(allp)):
                p = stp[rng.randrange(len(stp))]
            elif len(allp):
                p = allp[rng.randrange(len(allp))]
            else:                                          # no uncovered raster point at all: anywhere in the strip
                r = strip[rng.randrange(len(strip))] if strip else (0, 0, S, S)
                p = (rng.uniform(r[0], r[2]), rng.uniform(r[1], r[3]))
            x, y = float(p[0]) + rng.uniform(-0.05, 0.05), float(p[1]) + rng.uniform(-0.05, 0.05)
            a = _angle(X, x, y, rng)
            hh = half_extent(a)
            x, y = min(max(x, hh), S - hh), min(max(y, hh), S - hh)
            pv = total_pen(S, X, (x, y, a))
            if best is None or pv < best[0]:
                best = (pv, (x, y, a))
        if m is None and pen_full is not None and best[0] > pen_full:
            break
        tot += best[0]
        X = np.vstack([X, np.array(best[1])])
        added += 1
    return X, tot, added


# ================================================================== simple placers

def farthest_points(S, sq, m, rng, samples=4000):
    """m axis squares at the sampled points farthest from every centre (cheap proxy for the largest hole; bulbgen's
    10-19 top-up, here with a random.Random)"""
    X = [np.array(q[:2], float) for q in sq]
    out = []
    for _ in range(m):
        Pts = np.array([[rng.uniform(0.5, S - 0.5), rng.uniform(0.5, S - 0.5)] for _ in range(samples)])
        C = np.array(X) if X else np.full((1, 2), -1e9)
        d = np.min(np.hypot(Pts[:, None, 0] - C[None, :, 0], Pts[:, None, 1] - C[None, :, 1]), axis=1)
        p = Pts[int(np.argmax(d))]
        out.append((float(p[0]), float(p[1]), 0.0))
        X.append(p)
    return out


# ================================================================== registry

PLACERS = {}


def placer(name, **space):
    """register f(S, sq, m, rng, region, **params) -> (X, extra info) with its parameter space (pk.moves.draw specs)"""
    def deco(f):
        f.space, f.name = space, name
        PLACERS[name] = f
        return f
    return deco


def _draw(spec, rng):
    from .moves import draw
    return draw(spec, rng)


@placer('bestofK', K=('choice', [1, 4, 16, 64, 256]), LAM=0.05, SNAP=1, SNAP_R=0.3, h=0.04, TILT=1, PUSH=1,
        RANDANG=0, p_region=0.75, ov_full=0.02, m_max=10 ** 6)
def _p_bestofk(S, sq, m, rng, region, **P):
    """noisy greedy best fit: K random raster poses (palette angles), slide each to its local score minimum, keep the best, SAT-push"""
    kn = {k: P[k] for k in BESTOFK_DEFAULTS}
    kn['m_max'] = P['m_max']
    angs = palette(S, sq, kn['TILT']) + [rng.uniform(-45, 45) for _ in range(kn['RANDANG'])]
    X, added = bestofk_candidate(S, sq, m, rng, kn, angs, region=region, p_region=P['p_region'], ov_full=P['ov_full'])
    return X, dict(angles=[round(a, 2) for a in angs])


@placer('leastbad', K=('choice', [1, 4, 16, 64, 256]), p_region=0.75, pen_full=0.05, m_max=10 ** 6)
def _p_leastbad(S, sq, m, rng, region, **P):
    """least total penetration of K random poses (uncovered raster points, angles 0 / a neighbour tilt / uniform)"""
    X, pen, added = leastbad_place(S, sq, region or [], m, rng, P['K'], p_region=P['p_region'], pen_full=P['pen_full'],
                                   m_max=P['m_max'])
    return X, dict(pen_place=round(pen, 6))


def _fillfast(S, sq, m, rng, mode, tau, drop):
    import fillfast
    T = np.asarray(sq, float).reshape(-1, 3)
    X = fillfast.fill(S, T, mode, tau=tau)
    U = X[len(T):]
    if m is not None and len(U) > m:                      # more room than asked: keep m of them
        idx = list(range(len(U)))
        keep = sorted(rng.sample(idx, m)) if drop == 'random' else idx[:m]
        U = U[keep]
    return np.r_[T, U], dict(room=int(len(X) - len(T)))


@placer('exact', drop=('choice', ['random', 'first']))
def _p_exact(S, sq, m, rng, region, **P):
    """fillfast exact: largest set of non-overlapping upright squares among its candidates"""
    return _fillfast(S, sq, m, rng, 'exact', 0.0, P['drop'])


@placer('greedy', tau=('choice', [0.0, 0.01, 0.03, 0.06]), drop=('choice', ['random', 'first']))
def _p_greedy(S, sq, m, rng, region, **P):
    """fillfast greedy with overlap tolerance tau (overcounts by design; the relax repairs it)"""
    return _fillfast(S, sq, m, rng, 'greedy', P['tau'], P['drop'])


@placer('farthest', samples=4000)
def _p_farthest(S, sq, m, rng, region, **P):
    """upright squares at the sampled points farthest from every centre"""
    if m is None:
        raise ValueError("farthest needs a count (it never runs out of room)")
    add = farthest_points(S, sq, m, rng, P['samples'])
    return np.vstack([np.asarray(sq, float)] + ([np.array(add)] if add else [])), {}


@placer('uniform')
def _p_uniform(S, sq, m, rng, region, **P):
    """centres U(0.5, s - 0.5)^2, angles U(0, 90); ignores the region"""
    if m is None:
        raise ValueError("uniform needs a count")
    add = [(rng.uniform(0.5, S - 0.5), rng.uniform(0.5, S - 0.5), rng.uniform(0, 90)) for _ in range(m)]
    return np.vstack([np.asarray(sq, float)] + ([np.array(add)] if add else [])), {}


def blf_uprights(S, sq, rng, corner=0, order='yx', h=0.02, m=None):
    """Bottom-left-fill of upright unit squares around fixed squares (10-21, bulbgen multi-band): raster h (1/h an
    integer, so uprights abut exactly; against tilted squares the gap is < h), feasible lower-left cells = a free
    (1/h)^2 block (summed-area table), place the first in `order` ('yx': lowest then leftmost; 'xy': leftmost then
    lowest), seen from `corner` (0 BL, 1 BR, 2 TR, 3 TL: coordinates flipped so that corner is the origin).  Stops when
    no block is free or m are placed.  -> list of (x, y, 0)"""
    w = int(round(1 / h)); assert abs(w * h - 1) < 1e-12
    fx, fy = corner in (1, 2), corner in (2, 3)
    T = [((S - x) if fx else x, (S - y) if fy else y, (-a if fx != fy else a)) for x, y, a in np.asarray(sq, float)]
    N = int(math.floor(S / h + 1e-9))
    xs = (np.arange(N) + 0.5) * h
    X, Y = np.meshgrid(xs, xs)
    occ = np.zeros((N, N), bool)
    for x, y, a in T:
        i0, i1 = max(0, int((x - 0.75) / h)), min(N, int((x + 0.75) / h) + 2)
        j0, j1 = max(0, int((y - 0.75) / h)), min(N, int((y + 0.75) / h) + 2)
        occ[j0:j1, i0:i1] |= _inside(X[j0:j1, i0:i1], Y[j0:j1, i0:i1], x, y, a)
    out = []
    while m is None or len(out) < m:
        A = np.zeros((N + 1, N + 1), np.int32)
        A[1:, 1:] = occ.cumsum(0).cumsum(1)
        B = A[w:, w:] - A[:-w, w:] - A[w:, :-w] + A[:-w, :-w]          # occupied cells in the block at (j, i)
        F = B == 0
        if not F.any():
            break
        j, i = np.unravel_index(np.argmax(F if order == 'yx' else F.T), F.shape if order == 'yx' else F.T.shape)
        if order != 'yx':
            i, j = j, i
        occ[j:j + w, i:i + w] = True
        x, y = (i + w / 2) * h, (j + w / 2) * h
        out.append(((S - x) if fx else x, (S - y) if fy else y, 0.0))
    return out


@placer('blf', corner=('choice', [0, 1, 2, 3]), order=('choice', ['yx', 'xy']), h=0.02)
def _p_blf(S, sq, m, rng, region, **P):
    """bottom-left-fill uprights from one corner"""
    add = blf_uprights(S, sq, rng, P['corner'], P['order'], P['h'], m)
    return np.vstack([np.asarray(sq, float)] + ([np.array(add)] if add else [])), dict(room=len(add))


# ================================================================== relax + the one entry point

def relax(S, X, schedule, tmp=None, seed=None):
    """fq quench of (S, X) with a named schedule (pk/schedules.py) -> (s, array) or None if infeasible / failed"""
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import hop
    from .schedules import get
    opt = get(schedule) if seed is None else get(schedule).but(seed=seed)
    own = tmp is None
    td = tempfile.TemporaryDirectory(prefix='pkfill_') if own else None
    try:
        r = hop.quench(S, [tuple(map(float, q)) for q in np.asarray(X, float)], td.name if own else tmp, opt=opt)
    except Exception:                                      # fq divergence (inf in its JSON) etc.: a failed relax
        r = None
    finally:
        if own:
            td.cleanup()
    return None if r is None else (float(r[0]), np.asarray(r[1], float))


def fill(S, sq, rng=None, placer='bestofK', m=None, region=None, relax_with=None, tmp=None, **fixed):
    """Add squares to (S, sq).  m: count to add, or None = until full (placers that support it).  region: optional
    rectangles (x0, y0, x1, y1) to bias placement toward.  relax_with: a schedule name (pk/schedules.py) or None.
    Unspecified placer parameters are drawn from its space.  -> (X, info); X is the relaxed packing if relaxed and
    feasible, else the raw placement (info['s'] None tells which)."""
    rng = rng or random.Random()
    f = PLACERS[placer]
    P = {k: (fixed[k] if k in fixed else _draw(v, rng)) for k, v in f.space.items()}
    unknown = set(fixed) - set(f.space)
    if unknown:
        raise TypeError(f'{placer}: unknown parameters {sorted(unknown)}')
    n0 = len(sq)
    X, extra = f(S, np.asarray(sq, float).reshape(-1, 3), m, rng, region, **P)
    info = dict(placer=placer, params={k: v for k, v in P.items() if isinstance(v, (int, float, str))},
                m=m, added=int(len(X) - n0), count=int(len(X)), pen=round(packing_pen(S, X), 6), relax=relax_with, s=None,
                **extra)
    if relax_with:
        r = relax(S, X, relax_with, tmp=tmp)
        if r is not None:
            info['s'] = r[0]
            X = r[1]
    return X, info
