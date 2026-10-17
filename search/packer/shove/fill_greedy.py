"""Mechanism A: noisy greedy best fit (Evan's sketch, 10-16).

Each placement: rasterise the covered region (squares + everything outside the box) at step h; for every angle in a
small palette (0 plus the dominant tilts of the squares near the uncovered region) compute by FFT convolution
    overlap(x, y)  = covered pixels inside the unit square at (x, y, angle)
    contact(x, y)  = covered pixels in a thin ring (width RING) around it
and the cheap score  ov - LAM * contact  (low = tight fit: no overlap, hugging walls / squares).
Then draw K random poses (uniform over uncovered pixels as centres, uniform angle from the palette) and keep the best
by the score ("best of K random options": K = 1 is pure random, large K is near-greedy).  With SNAP on, the chosen
pose slides to the local score minimum within SNAP_R at its angle (pushes it against its neighbours).  Paint it,
repeat until m squares.  Small overlaps are allowed (the harness's fq quench relaxes them); when nothing fits the
least-bad pose is taken.  Returns C candidates (different random draws) per call.

Knobs (env vars, read at call time): SHOVE_K (comma list, cycled over the candidates; default '64,16,16'), SHOVE_C (3), SHOVE_LAM (0.05), SHOVE_SNAP (1),
SHOVE_SNAP_R (0.3), SHOVE_H (0.04), SHOVE_TILT (1: use neighbourhood tilts; 0: uprights only),
SHOVE_PUSH (1: continuous SAT push-out of each placed pose), SHOVE_RANDANG (0: extra random angles per candidate).
"""
from __future__ import annotations
import math, os
import numpy as np
from scipy.signal import fftconvolve
import harness as H

RING = 0.08
PAD = 1.0


def _knobs():
    e = os.environ.get
    return dict(K=[int(v) for v in e('SHOVE_K', '64,16,16').split(',')], C=int(e('SHOVE_C', 3)), LAM=float(e('SHOVE_LAM', 0.05)),
                SNAP=int(e('SHOVE_SNAP', 1)), SNAP_R=float(e('SHOVE_SNAP_R', 0.3)), h=float(e('SHOVE_H', 0.04)),
                TILT=int(e('SHOVE_TILT', 1)), PUSH=int(e('SHOVE_PUSH', 1)),
                RANDANG=int(e('SHOVE_RANDANG', 0)))


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


def palette(S, sq, meta, tilt):
    if not tilt:
        return [0.0]
    sq = np.asarray(sq, float)
    xs, ys, free = H.free_mask(S, sq, h=0.1)
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
    he = H.half_extent(a)
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


def one_candidate(S, sq, meta, m, rng, kn, angs):
    G = Grid(S, sq, kn['h'])
    add = []
    snapw = int(round(kn['SNAP_R'] / kn['h']))
    for _ in range(m):
        maps = {a: G.score(a, kn['LAM']) for a in angs}
        free = np.argwhere(~G.cov)
        if not len(free):
            free = np.argwhere(np.ones_like(G.cov))
        best = None
        for _k in range(max(1, kn['K'])):
            j, i = free[rng.randrange(len(free))]
            a = angs[rng.randrange(len(angs))]
            sc = maps[a][1]
            if kn['SNAP']:
                j, i = _descend(sc, j, i, snapw)
            if best is None or sc[j, i] < best[0]:
                best = (sc[j, i], j, i, a)
        _, j, i, a = best
        pose = (float(G.xs[i]), float(G.xs[j]), float(a % 90.0))
        if kn['PUSH']:
            cur = np.vstack([np.asarray(sq, float)] + ([np.array(add)] if add else []))
            pose = tuple(map(float, _push(S, cur, pose)))
        add.append(pose)
        G.paint(pose)
    return np.vstack([np.asarray(sq, float), np.array(add)])


def fill(S, sq, meta, m, rng):
    kn = _knobs()
    angs = palette(S, sq, meta, kn['TILT'])
    out = []
    for c in range(kn['C']):
        kc = dict(kn, K=kn['K'][c % len(kn['K'])])
        A = angs + [rng.uniform(-45, 45) for _ in range(kn['RANDANG'])]
        out.append(one_candidate(S, sq, meta, m, rng, kc, A))
    return out
