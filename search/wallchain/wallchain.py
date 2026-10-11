#!/usr/bin/env python3
"""Wall-chain elimination for n = 17 closed-cell states (float demonstrator).

The squares whose centres lie in the corner and side cells of jlevy/squares' exp-247 cover (24 cells, cap
U = 1169/250) form a ring along the container walls.  Along one wall the joint problem is a path: walk the wall,
keep for each square the function F(q, t) = least position p along the wall that square can reach with depth q
(distance of its centre from the wall) and angle t (mod 90 deg), given that every earlier square in the chain is in
its cell, inside the container and disjoint from its chain neighbours.

  F_B(qB, tB) = max(cell_lo_B, min_{qA, tA} F_A(qA, tA) + b(qB - qA, tA, tB))       (infeasible if > cell_hi_B)

where b is the least centre gap along the wall for which unit squares at the given relative depth and angles are
disjoint (separating-axis theorem, four candidate normals).  Valid because |qB - qA| < 1 for wall cells: the
Minkowski difference of two unit squares contains the disc of radius 1, so the forbidden gaps form an interval
containing 0 and the allowed gaps to the right are an up-set [b, oo).  Dropping non-neighbour pairs is a
relaxation, so a chain shown infeasible stays infeasible with them.  At a corner the square changes frame
(proper rotation, angles unchanged mod 90): its depth on the next wall is U - p, its position there is q.

FLOAT DEMONSTRATOR: F is sampled on a grid, so 'infeasible' here is evidence, not proof (the grid minimum can
exceed the true minimum).  A sound version replaces grid points by boxes and b by a lower bound over each box.
"""
import json, sys, os
import numpy as np
from fractions import Fraction as Fr

U = 1169 / 250
SIGMA = 1.0     # square side used by gap / wall_floor (set < 1 to give float grids room; see ringcond --shrink)
COVER = os.path.expanduser('~/math/_untrusted-third-party/jlevy-squares/packing/campaign/series/'
                           'series-000-smoke-and-calibration/results/exp-247-n17-unique-state-cover/run-001/receipt.json')

RING = ['corner-SW', 'side-S0', 'side-S1', 'side-S2', 'corner-SE', 'side-E0', 'side-E1', 'side-E2', 'corner-NE',
        'side-N2', 'side-N1', 'side-N0', 'corner-NW', 'side-W2', 'side-W1', 'side-W0']
WALL_OF = {'S': ['side-S0', 'side-S1', 'side-S2'], 'E': ['side-E0', 'side-E1', 'side-E2'],
           'N': ['side-N2', 'side-N1', 'side-N0'], 'W': ['side-W2', 'side-W1', 'side-W0']}
# frame maps (x, y) -> (p, q): proper rotations taking each wall to the bottom, traversed counter-clockwise
FRAME = {'S': lambda x, y: (x, y), 'E': lambda x, y: (y, U - x), 'N': lambda x, y: (U - x, U - y),
         'W': lambda x, y: (U - y, x)}
# corner -> (incoming wall, outgoing wall) in counter-clockwise order
CORNER = {'corner-SE': ('S', 'E'), 'corner-NE': ('E', 'N'), 'corner-NW': ('N', 'W'), 'corner-SW': ('W', 'S')}


def load_cells():
    r = json.load(open(COVER))
    out, names = {}, []
    for c in r['cells']:
        xs = [float(Fr(v[0])) for v in c['vertices']]
        ys = [float(Fr(v[1])) for v in c['vertices']]
        out[c['name']] = (min(xs), max(xs), min(ys), max(ys), c['kind'])
        names.append(c['name'])
    return out, names


def rect_in_frame(cell, wall):
    x0, x1, y0, y1, _ = cell
    pts = [FRAME[wall](x, y) for x in (x0, x1) for y in (y0, y1)]
    ps, qs = [a for a, _ in pts], [b for _, b in pts]
    return min(ps), max(ps), min(qs), max(qs)


def half_width(phi, t):
    d = phi - t
    return 0.5 * SIGMA * (np.abs(np.cos(d)) + np.abs(np.sin(d)))


def gap(dq, tA, tB):
    """least dp >= 0 with A (centre 0, angle tA) and B (centre (dp, dq), angle tB) disjoint; broadcasting."""
    best = np.full(np.broadcast(dq, tA, tB).shape, np.inf)
    for phi in (tA, tA + np.pi / 2, tB, tB + np.pi / 2):
        nx, ny = np.cos(phi), np.sin(phi)
        H = half_width(phi, tA) + half_width(phi, tB)
        s = np.sign(nx)
        with np.errstate(divide='ignore', invalid='ignore'):
            cand = np.where(np.abs(nx) > 1e-12, (H - s * ny * dq) / np.abs(nx), np.inf)
        best = np.minimum(best, cand)
    return best


class Grid:
    def __init__(self, nq_per_unit=200, nt=90):
        self.nqu = nq_per_unit
        self.t = np.arange(nt) * (np.pi / 2) / nt          # [0, 90deg)

    def qs(self, q0, q1):
        n = max(2, int(round((q1 - q0) * self.nqu)) + 1)
        return np.linspace(q0, q1, n)


def wall_floor(t):
    """centre-to-wall distance needed by a unit square at angle t (mod 90)"""
    return 0.5 * SIGMA * (np.abs(np.cos(t)) + np.abs(np.sin(t)))


def step(FA, qA, qB, t, cellB_p, chunk=8):
    """FB(qB, t) from FA(qA, t) along one wall.  FA shape (len(qA), nt)."""
    nt = len(t)
    dq = qB[:, None] - qA[None, :]                                  # (nB, nA)
    FB = np.full((len(qB), nt), np.inf)
    finite = np.isfinite(FA)
    for jA in range(nt):
        if not finite[:, jA].any():
            continue
        fa = FA[:, jA]
        # g[b, a, jB] = gap(qB-qA, tA, tB)
        g = gap(dq[:, :, None], t[jA], t[None, None, :])          # (nB, nA, nt)
        cand = (fa[None, :, None] + g).min(axis=1)                  # (nB, nt)
        FB = np.minimum(FB, cand)
    p0, p1 = cellB_p
    FB = np.maximum(FB, p0)
    FB = np.maximum(FB, wall_floor(t)[None, :])                     # side wall at p = 0 (matters near corners)
    FB[FB > p1] = np.inf
    FB[qB[:, None] < wall_floor(t)[None, :]] = np.inf               # own wall at q = 0
    return FB


def corner_turn(F1, q1, cell_p1_hi, q2grid, t):
    """corner square: F1(q1, t) = least p1 on the incoming wall.  Outgoing wall: q2 = U - p1, p2 = q1.
    G(q2, t) = least q1 with F1(q1, t) <= U - q2."""
    G = np.full((len(q2grid), len(t)), np.inf)
    for j in range(len(t)):
        f = F1[:, j]
        for i, q2 in enumerate(q2grid):
            ok = f <= U - q2 + 1e-12
            if ok.any():
                G[i, j] = q1[ok].min()
    G[q2grid[:, None] < wall_floor(t)[None, :]] = np.inf
    return G


def run_chain(order, cells, grid, start_F=None, verbose=True):
    """order: list of occupied ring cells in counter-clockwise order (a contiguous path).  Returns per-square
    feasibility summary; the chain is infeasible as soon as some F is identically inf."""
    t = grid.t
    log = []
    wall = None
    F = q = None
    for k, name in enumerate(order):
        if name.startswith('corner'):
            win, wout = CORNER[name]
        else:
            win = wout = name.split('-')[1][0]
        if k == 0:
            wall = wout
            p0, p1, q0, q1 = rect_in_frame(cells[name], wall)
            q = grid.qs(q0, q1)
            F = np.maximum(np.full((len(q), len(t)), p0), wall_floor(t)[None, :]) if start_F is None else start_F
            F[q[:, None] < wall_floor(t)[None, :]] = np.inf
            F[F > p1] = np.inf
        else:
            assert win == wall, (name, win, wall)
            p0, p1, q0, q1 = rect_in_frame(cells[name], wall)
            qB = grid.qs(q0, q1)
            F = step(F, q, qB, t, (p0, p1))
            q = qB
            if name.startswith('corner'):
                # turn: incoming frame wall, outgoing wout
                P0, P1, Q0, Q1 = rect_in_frame(cells[name], wout)
                q2 = grid.qs(Q0, Q1)
                F = corner_turn(F, q, p1, q2, t)
                F[F > P1] = np.inf
                F = np.maximum(F, P0)
                q = q2
                wall = wout
        feas = np.isfinite(F)
        slack = None
        if feas.any():
            p_hi = rect_in_frame(cells[name], wall)[1]
            slack = float(p_hi - F[feas].min())
        log.append(dict(square=name, wall=wall, feasible_states=int(feas.sum()), states=int(feas.size),
                        slack_to_cell_end=slack))
        if verbose:
            print(f'{k:2d} {name:10s} wall {wall}  feasible {feas.sum():6d}/{feas.size}  '
                  f'slack {slack if slack is None else round(slack, 4)}', flush=True)
        if not feas.any():
            break
    return log, F, q


def ring_path(mask, names, start=None):
    """occupied ring cells, counter-clockwise, as one path starting after a vacant ring cell (or at start)."""
    occ = {names[i] for i in range(24) if mask >> i & 1}
    ring = [c for c in RING if c in occ]
    vac = [i for i, c in enumerate(RING) if c not in occ]
    if start is None and vac:
        i0 = (vac[-1] + 1) % len(RING)
    else:
        i0 = RING.index(start or 'corner-SW')
    path = []
    for k in range(len(RING)):
        c = RING[(i0 + k) % len(RING)]
        if c in occ:
            path.append(c)
        elif path:
            path.append(None)      # gap marker
    return path


if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument('mask', type=int)
    a.add_argument('--nq', type=int, default=200, help='depth grid points per unit')
    a.add_argument('--nt', type=int, default=90, help='angle grid points on [0, 90deg)')
    o = a.parse_args()
    cells, names = load_cells()
    path = ring_path(o.mask, names)
    print('ring path:', path)
    # split at gaps into maximal open chains
    chains, cur = [], []
    for c in path:
        if c is None:
            if cur:
                chains.append(cur)
            cur = []
        else:
            cur.append(c)
    if cur:
        chains.append(cur)
    g = Grid(o.nq, o.nt)
    for ch in chains:
        print(f'--- chain of {len(ch)}: {ch[0]} .. {ch[-1]}')
        run_chain(ch, cells, g)
