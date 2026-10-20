"""Corner-shove harness (10-16; Evan: "an L routed through the interior").

Operator: take a packing (n squares, side s), pick one corner's grid-aligned block (cornersig's corner groups), grow the
box by d and move the block with its corner while everything else stays with the opposite corner.  A strip of width d
opens between the block and the jumble (plus along the walls beyond the block); a FILL adds m squares (any poses,
small overlaps allowed: fq relaxes them), and EVALUATE quenches the result with the explorer's fq stages.

Fill interface (one module per mechanism, shove/fill_<name>.py):

    def fill(S, sq, meta, m, rng) -> list of candidate arrays
        S     new box side (s + d)
        sq    (n, 3) float array: x, y, angle in degrees, after the shove (box [0, S]^2)
        meta  dict: corner ('LL'|'LR'|'UL'|'UR'), d, block (indices into sq), s0 (old side), n0, k0 = floor(s0),
              strip (list of axis rectangles (x0, y0, x1, y1) covering the opened strip: a guide, not a constraint),
              cells (clean upright poses on the strip's unit lattices: usually well short of 2k + 1)
        m     squares to add
        rng   random.Random
      Each candidate is an (n + m', 3) array (m' may differ from m; existing squares may move).  The harness evaluates
      each candidate at its own count.

Helpers for fills: free_mask (raster of uncovered points), pen_one (max penetration of one pose against the others and
the walls), half_extent (half the bounding-box width of a tilted unit square).
"""
from __future__ import annotations
import math, os, sys, tempfile
from collections import deque
import numpy as np

PK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PK); sys.path.insert(0, os.path.join(PK, 'analysis'))

CORNERS = ('LL', 'LR', 'UL', 'UR')        # (y side, x side): 'LR' = lower right


from pk.fill import half_extent  # noqa: E402  (10-20: moved to pk/fill.py)

# ------------------------------------------------------------------ corner blocks

def corner_blocks(s, sq, tol=0.02):
    """{corner: indices of that corner's grid-aligned block} (cornersig v2's main component per corner, computed on
    rattler-snapped coordinates; indices refer to sq).  No L peeling: an L-addition counts as part of the blocks."""
    import cornersig as CS
    sq = np.asarray(sq, float)
    q = CS.snap_rattlers(float(s), sq, tol)
    ax = np.nonzero(CS._axis(q))[0]
    A = q[ax]
    xl, xr = CS._on(A[:, 0] - 0.5, tol), CS._on(s - A[:, 0] - 0.5, tol)
    yb, yt = CS._on(A[:, 1] - 0.5, tol), CS._on(s - A[:, 1] - 0.5, tol)
    out = {}
    for c in CORNERS:
        xe = A[:, 0] - 0.5 if c[1] == 'L' else s - A[:, 0] - 0.5
        ye = A[:, 1] - 0.5 if c[0] == 'L' else s - A[:, 1] - 0.5
        ok = (xl if c[1] == 'L' else xr) & (yb if c[0] == 'L' else yt)
        cells = {(int(round(x)), int(round(y))): i for i, (x, y) in enumerate(zip(xe, ye)) if ok[i]}
        comps = CS._components(cells)
        cand = [cp for cp in comps if any(i == 0 or j == 0 for i, j in cp)]
        main = next((cp for cp in comps if (0, 0) in cp), None) or (max(cand, key=len) if cand else set())
        out[c] = sorted(int(ax[cells[cell]]) for cell in main)
    return out


# ------------------------------------------------------------------ the shove

def pair_pen(sq):
    """(i, j, depth) for every pair with centre distance < sqrt 2 (SAT over the 4 axes; depth 0 = separated)"""
    sq = np.asarray(sq, float)
    x, y, a = sq[:, 0], sq[:, 1], np.radians(sq[:, 2])
    dx = x[None, :] - x[:, None]; dy = y[None, :] - y[:, None]
    I, J = np.nonzero(np.triu(dx * dx + dy * dy < 2.0, 1))
    if len(I) == 0:
        return I, J, np.zeros(0)
    ddx, ddy, ai, aj = dx[I, J], dy[I, J], a[I], a[J]
    best = np.full(len(I), np.inf)
    for phi in (ai, ai + math.pi / 2, aj, aj + math.pi / 2):
        ux, uy = np.cos(phi), np.sin(phi)
        h1 = 0.5 * (np.abs(np.cos(ai - phi)) + np.abs(np.sin(ai - phi)))
        h2 = 0.5 * (np.abs(np.cos(aj - phi)) + np.abs(np.sin(aj - phi)))
        best = np.minimum(best, h1 + h2 - np.abs(ddx * ux + ddy * uy))
    return I, J, np.maximum(best, 0.0)


def fix_levels(sq, t, b, o, groups=(), tol=1e-7, max_iter=10000):
    """10-18 (collision fix): per-square level t (0 = moves with the corner block by b, 1 = with the rest by o;
    offset = b + t (o - b)).  While some pair's penetration exceeds its source penetration by > tol and the two levels
    differ, the higher square drops to the lower level (a group, e.g. an island, drops as a whole).  t only decreases
    and all-equal levels move nothing, so this terminates with no new overlap.  -> (t, number of squares lowered)"""
    sq = np.asarray(sq, float); t = np.asarray(t, float).copy(); t0 = t.copy()
    b = np.asarray(b, float); o = np.asarray(o, float)
    I, J, D0 = pair_pen(sq)
    gid = -np.ones(len(sq), int)
    for g, G in enumerate(groups):
        gid[np.asarray(G, int)] = g
    for _ in range(max_iter):
        off = b[None, :] + t[:, None] * (o - b)[None, :]
        P = sq.copy(); P[:, :2] += off
        I1, J1, D1 = pair_pen(P)
        d0 = dict(zip(zip(I.tolist(), J.tolist()), D0.tolist()))
        bad = [(i, j) for i, j, dd in zip(I1.tolist(), J1.tolist(), D1.tolist())
               if dd > d0.get((i, j), 0.0) + tol and abs(t[i] - t[j]) > 1e-12]
        if not bad:
            break
        for i, j in bad:
            hi, lo = (i, j) if t[i] > t[j] else (j, i)
            if t[hi] <= t[lo]:
                continue
            if gid[hi] >= 0:
                m = gid == gid[hi]; t[m] = np.minimum(t[m], t[lo])
            else:
                t[hi] = t[lo]
    return t, int((t < t0 - 1e-12).sum())


def shove_fix_on():
    return os.environ.get('SHOVE_FIX', '0') not in ('0', '')


def shove(s, sq, corner, d=1.0, block=None, fix=None):
    """grow the box to s + d; the block moves with `corner`, everything else with the opposite corner.
    -> (S, sq2, meta)"""
    sq = np.asarray(sq, float).copy()
    if block is None:
        block = corner_blocks(s, sq)[corner]
    blk = np.zeros(len(sq), bool); blk[list(block)] = True
    bx = d if corner[1] == 'R' else 0.0          # block's offset (moves with its corner)
    by = d if corner[0] == 'U' else 0.0
    ox, oy = d - bx, d - by                      # everything else: with the opposite corner
    fix = shove_fix_on() if fix is None else fix
    nfix = 0
    if fix:                                      # 10-18: squares the shove would drive into the block move with it
        tl, nfix = fix_levels(sq, (~blk).astype(float), (bx, by), (ox, oy))
        blk = tl < 0.5
        block = np.nonzero(blk)[0]
    sq[blk, 0] += bx; sq[blk, 1] += by
    sq[~blk, 0] += ox; sq[~blk, 1] += oy
    S = s + d
    meta = dict(corner=corner, d=d, block=sorted(int(i) for i in block), s0=float(s), n0=len(sq),
                k0=int(math.floor(s + 1e-9)), strip=strip_rects(S, sq, blk, corner, d), fixed=nfix)
    meta['cells'] = clean_cells(S, sq, meta)
    return S, sq, meta


def strip_rects(S, sq, blk, corner, d):
    """axis rectangles covering the opened strip: the block's bounding box grown by d toward the jumble, minus the
    block's own bounding box, plus the two wall strips beyond it.  A guide for fills (cells of the block's lattice)."""
    if not blk.any():
        return []
    B = sq[blk]
    x0, x1 = B[:, 0].min() - 0.5, B[:, 0].max() + 0.5
    y0, y1 = B[:, 1].min() - 0.5, B[:, 1].max() + 0.5
    R = []
    if corner[1] == 'R':                                  # vertical strip on the block's inner (left) face, + wall strip
        R.append((x0 - d, 0.0 if corner[0] == 'L' else y0 - d, x0, y1 if corner[0] == 'L' else S))
        R.append((S - d, 0.0 if corner[0] == 'U' else y1, S, y0 - d if corner[0] == 'U' else S))
    else:
        R.append((x1, 0.0 if corner[0] == 'L' else y0 - d, x1 + d, y1 if corner[0] == 'L' else S))
        R.append((0.0, 0.0 if corner[0] == 'U' else y1, d, y0 - d if corner[0] == 'U' else S))
    if corner[0] == 'U':                                  # horizontal strip on the block's inner (lower) face, + wall
        R.append((0.0 if corner[1] == 'L' else x0 - d, y0 - d, x1 if corner[1] == 'L' else S, y0))
        R.append((0.0 if corner[1] == 'R' else x1, S - d, x0 - d if corner[1] == 'R' else S, S))
    else:
        R.append((0.0 if corner[1] == 'L' else x0 - d, y1, x1 + d if corner[1] == 'L' else S, y1 + d))
        R.append((0.0 if corner[1] == 'R' else x1, 0.0, x0 - d if corner[1] == 'R' else S, d))
    return [r for r in R if r[2] - r[0] > 1e-9 and r[3] - r[1] > 1e-9]


# ------------------------------------------------------------------ helpers for fills

from pk.fill import _inside, free_mask  # noqa: E402  (10-20: moved to pk/fill.py)


def pen_one(S, sq, pose, skip=None):
    """max SAT penetration of unit square `pose` (x, y, deg) against the squares in sq (except index skip) and the
    walls of [0, S]^2; 0 = fits"""
    from pk.packing import Packing
    x, y, a = pose
    h = half_extent(a)
    wall = max(0.0, h - x, h - y, x + h - S, y + h - S)
    sq = np.asarray(sq, float)
    dd = np.hypot(sq[:, 0] - x, sq[:, 1] - y)
    near = np.nonzero(dd < math.sqrt(2) + 1e-9)[0]
    worst = wall
    for j in near:
        if j == skip:
            continue
        worst = max(worst, Packing(100.0, np.array([pose, sq[j]]) + [40, 40, 0]).max_pen())
    return worst


# ------------------------------------------------------------------ evaluation (the explorer's fq stages)

def evaluate(S, sq, loosens=(1.0, 1.02, 1.05), polish=True):
    """fq SCREEN (ALM squeeze from the overlapping state) at each loosen, then POLISH the best feasible; grid check.
    -> dict(n, s, grid, sq) with s None if nothing feasible"""
    import hop
    from pk.engine import SCREEN, POLISH
    from layout import full_lines
    import mcmin
    sq = [tuple(map(float, q)) for q in np.asarray(sq, float)]
    n = len(sq); k = math.ceil(math.sqrt(n) - 1e-12)
    best = None
    with tempfile.TemporaryDirectory(prefix='shove_') as tmp:
        for L in loosens:
            try:                                          # 10-17 (agent A): fq can diverge (min_gap inf -> bad JSON)
                r = hop.quench(S, sq, tmp, opt=SCREEN.but(loosen=L))
            except Exception:
                r = None
            if r and (best is None or r[0] < best[0]):
                best = (r[0], r[1])
        if best is None:
            return dict(n=n, s=None, grid=None)
        if polish:
            try:
                r = hop.quench(best[0], best[1], tmp, opt=POLISH)
            except Exception:
                r = None
            if r and r[0] < best[0]:
                best = (r[0], r[1])
        p = os.path.join(tmp, 'e.txt'); mcmin.write_deg(p, best[0], best[1])
        g = full_lines(p, k)
    return dict(n=n, s=best[0], grid=bool(g) or best[0] >= k - 1e-9, sq=np.asarray(best[1], float))


# ------------------------------------------------------------------ baselines

def lstep(s, sq):
    """Friedman's L-step: n + 2k + 1 squares in side s + 1 (k = floor s): the packing in the lower left, a column of
    k + 1 upright squares on the right wall and a row of k on the top wall."""
    k = int(math.floor(s + 1e-9)); S = s + 1
    add = [(S - 0.5, 0.5 + j * (S - 1) / k, 0.0) for j in range(k + 1)] + \
          [(0.5 + i * (S - 2) / max(1, k - 1), S - 0.5, 0.0) for i in range(k)]
    return S, np.vstack([np.asarray(sq, float), np.array(add)])


def clean_cells(S, sq, meta):
    """upright unit-square poses that fit cleanly (penetration < 1e-9) on the unit lattices of the strip guide,
    anchored at either end of each rectangle; greedy, so the set is non-overlapping.  The opened region follows the
    block's staircase, so this is usually well short of the L-step's 2k + 1 (10-16 smoke test: 16 / 13 / 8 of
    21 / 23 / 25 at n = 101 / 126 / 151): fills must make up the rest (tilted squares, relaxation)."""
    sq = np.asarray(sq, float)
    add = []
    for x0, y0, x1, y1 in meta['strip']:
        nx, ny = max(1, int(math.floor(x1 - x0 + 1e-9))), max(1, int(math.floor(y1 - y0 + 1e-9)))
        for end in (0, 1):
            for i in range(nx):
                for j in range(ny):
                    pose = (x0 + 0.5 + i, y0 + 0.5 + j, 0.0) if end == 0 else (x1 - 0.5 - i, y1 - 0.5 - j, 0.0)
                    cur = np.vstack([sq] + ([np.array(add)] if add else []))
                    if H_pen(S, cur, pose) < 1e-9:
                        add.append(pose)
    return add


def H_pen(S, sq, pose):
    return pen_one(S, sq, pose)
