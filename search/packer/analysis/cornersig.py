#!/usr/bin/env python3
"""Corner signature of a packing (Evan 10-15): how the axis-aligned squares organise around the four corners.

Each axis-aligned square (tilt < 1 deg) lies on the unit lattice anchored at the left or right wall in x (left edge at
an integer distance from that wall, to --tol) and at the bottom or top wall in y.  A corner's group = the squares on
that corner's lattice that are face-connected (lattice 4-neighbours, i.e. touching) to the square in the corner itself.
Squares in no corner group are "floating" (shifted against the tilted structure).  Before counting, L-additions are
peeled (a full column on one wall + a full row on an adjacent wall, s -> s - 1, as in Friedman's L-step); the number
peeled is the L depth.

Signature = (n, U = axis-aligned count, corner groups sorted descending, floating, L depth); the unsorted groups are kept
as (LL, LR, UL, UR).  Expected shapes (Evan): strips: two large + two small; bulbs / bent strips: big, medium, small or
zero; diamonds: four similar.

  python3 analysis/cornersig.py frontier [--nmin 100 --nmax 324]      # signatures of the store-best packings + family codes
  python3 analysis/cornersig.py file PACKING.txt
"""
import argparse, math, os, sys
from collections import deque
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK)


def _axis(sq, ang=1.0):
    a = sq[:, 2] % 90.0
    return np.minimum(a, 90 - a) < ang


def _on(v, tol):
    return np.abs(v - np.round(v)) < tol


def peel_L(s, sq, tol=0.02):
    """remove L-additions: a full column (>= floor(s) axis squares, one per unit, on the wall's lattice) on the left or
    right wall plus a full row on the bottom or top wall; -> (s', sq', depth)"""
    depth = 0
    while s > 2:
        k = int(math.floor(s + 1e-9))
        ax = _axis(sq)
        done = False
        for cx in ('L', 'R'):
            for cy in ('B', 'T'):
                xe = sq[:, 0] - 0.5 if cx == 'L' else s - sq[:, 0] - 0.5      # edge distance from the x wall
                ye = sq[:, 1] - 0.5 if cy == 'B' else s - sq[:, 1] - 0.5
                col = ax & (np.abs(xe) < tol)
                row = ax & (np.abs(ye) < tol)
                if col.sum() >= k and row.sum() >= k:
                    keep = ~(col | row)
                    q = sq[keep].copy()
                    if cx == 'L': q[:, 0] -= 1
                    if cy == 'B': q[:, 1] -= 1
                    r = np.radians(q[:, 2] % 90.0); h = 0.5 * (np.abs(np.cos(r)) + np.abs(np.sin(r)))
                    if ((q[:, 0] - h < -1e-6) | (q[:, 1] - h < -1e-6) | (q[:, 0] + h > s - 1 + 1e-6) |
                            (q[:, 1] + h > s - 1 + 1e-6)).any():
                        continue                         # the rest does not fit in s - 1: not an L-addition
                    sq, s, depth, done = q, s - 1, depth + 1, True
                    break
            if done:
                break
        if not done:
            break
    return s, sq, depth


def _components(cells):
    """4-connected components of a set of lattice cells -> list of sets"""
    left, comps = set(cells), []
    while left:
        c0 = left.pop(); comp, dq = {c0}, deque([c0])
        while dq:
            i, j = dq.popleft()
            for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                nb = (i + d[0], j + d[1])
                if nb in left:
                    left.discard(nb); comp.add(nb); dq.append(nb)
        comps.append(comp)
    return comps


def _pen_one(sq, i, cand, s):
    """max SAT penetration of square i placed at cand against the others and the walls (local check)"""
    from pk.packing import Packing
    x, y, a = cand
    r = math.radians(a % 90.0); h = 0.5 * (abs(math.cos(r)) + abs(math.sin(r)))
    if x - h < -1e-9 or y - h < -1e-9 or x + h > s + 1e-9 or y + h > s + 1e-9:
        return 1.0
    d = np.hypot(sq[:, 0] - x, sq[:, 1] - y)
    near = np.nonzero((d < math.sqrt(2) + 1e-9) & (np.arange(len(sq)) != i))[0]
    worst = 0.0
    for j in near:
        worst = max(worst, Packing(100.0, np.array([cand, sq[j]]) + [40, 40, 0]).max_pen())
        if worst > 1e-7:
            break
    return worst


def snap_rattlers(s, sq, tol=0.02, rounds=3):
    """rattler normalisation (10-15: copies of the same 110 record gave different signatures): move every axis-aligned
    square that is off its walls' lattices onto the nearest lattice position (x from the left or right wall, y from the
    bottom or top), one coordinate at a time, if the move (< 0.5) creates no overlap; repeat a few rounds."""
    sq = np.asarray(sq, float).copy()
    ax = _axis(sq)
    for _ in range(rounds):
        moved = False
        for i in np.nonzero(ax)[0]:
            for dim in (0, 1):
                v = sq[i, dim]
                opts = [math.floor(v - 0.5 + 0.5) + 0.5, s - (math.floor(s - v - 0.5 + 0.5) + 0.5)]
                if any(abs(v - o) < tol for o in opts):
                    continue
                for o in sorted(opts, key=lambda o: abs(o - v)):
                    if abs(o - v) >= 0.5:
                        continue
                    cand = sq[i].copy(); cand[dim] = o
                    if _pen_one(sq, i, cand, s) <= 1e-7:
                        sq[i] = cand; moved = True; break
        if not moved:
            break
    return sq


def signature(s, sq, tol=0.02, snap=True):
    """v2 (10-15, Evan's one-dimension question): per corner, the largest 4-connected component of axis squares on that
    corner's lattice in BOTH axes (and whether it reaches the corner cell: a pocket in the corner must not hide the
    block, as at 126); then the axis squares in no corner group, split by anchoring: slid rows (on a bottom / top
    lattice in y only: they slide along x), slid columns (on a left / right lattice in x only), cut-off lattice squares
    (both axes on some lattice, not in the corner's main block) and free (neither)."""
    sq = np.asarray(sq, float)
    n = len(sq)
    if snap:
        sq = snap_rattlers(float(s), sq, tol)
    s2, q, depth = peel_L(float(s), sq, tol)
    ax = _axis(q)
    A = q[ax]
    xl, xr = _on(A[:, 0] - 0.5, tol), _on(s2 - A[:, 0] - 0.5, tol)
    yb, yt = _on(A[:, 1] - 0.5, tol), _on(s2 - A[:, 1] - 0.5, tol)
    groups, reach, member = {}, {}, np.zeros(len(A), bool)
    for c in ('LL', 'LR', 'UL', 'UR'):
        xe = A[:, 0] - 0.5 if c[1] == 'L' else s2 - A[:, 0] - 0.5
        ye = A[:, 1] - 0.5 if c[0] == 'L' else s2 - A[:, 1] - 0.5
        ok = (xl if c[1] == 'L' else xr) & (yb if c[0] == 'L' else yt)
        cells = {(int(round(x)), int(round(y))): i for i, (x, y) in enumerate(zip(xe, ye)) if ok[i]}
        comps = _components(cells)
        # the corner's block: the component containing the corner cell, else the largest one adjacent to a wall of
        # this corner (a pocket may occupy the corner itself)
        cand = [cp for cp in comps if any(i == 0 or j == 0 for i, j in cp)]
        main = next((cp for cp in comps if (0, 0) in cp), None) or (max(cand, key=len) if cand else set())
        groups[c] = len(main); reach[c] = (0, 0) in main
        for cell in main:
            member[cells[cell]] = True
    rest = ~member
    xa, ya = xl | xr, yb | yt
    return dict(n=n, U=int(ax.sum()), T=int(len(q) - ax.sum()), groups=groups, reach=reach,
                sorted=sorted(groups.values(), reverse=True),
                slid_rows=int((rest & ya & ~xa).sum()), slid_cols=int((rest & xa & ~ya).sum()),
                cutoff=int((rest & xa & ya).sum()), free=int((rest & ~xa & ~ya).sum()),
                floating=int(rest.sum()), L=depth)


def cmd_frontier(a):
    import sqlite3, json
    from pk.store import Store
    from pk import eventbench as EB
    import familymap as FM
    db = sqlite3.connect(f"file:{os.path.join(PK, 'runs', 'store.sqlite')}?mode=ro", uri=True)
    st = Store()
    out = []
    for n in range(a.nmin, a.nmax + 1):
        T = EB.frontier_at(n, None)
        if T >= math.ceil(math.sqrt(n) - 1e-12) - 1e-9:
            continue                                   # grid
        r = db.execute("SELECT id FROM packing WHERE n=? AND abs(s-?)<2e-9 AND status!='infeasible' ORDER BY added LIMIT 1",
                       (n, T)).fetchone()
        if not r:
            continue
        p = st.get(r[0]); sq = np.asarray(p.sq, float)
        sig = signature(p.s, sq)
        code = FM.classify(n, p.s, sq)['code']
        k = math.ceil(math.sqrt(n)); c = k * k - n
        g = sig['groups']
        print(f"{n:3d} k{k} c{c:2d} {code}  U {sig['U']:3d} T {sig['T']:3d} L{sig['L']}  corners {sig['sorted']}  "
              f"slid rows {sig['slid_rows']:2d} cols {sig['slid_cols']:2d} cutoff {sig['cutoff']:2d} free {sig['free']}  "
              f"(LL {g['LL']} LR {g['LR']} UL {g['UL']} UR {g['UR']}{' corner-pocket' if not all(sig['reach'][c] or g[c] == 0 for c in g) else ''})")
        out.append(dict(code=code, **sig))
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=0)


def cmd_file(a):
    with open(a.path) as f:
        n, s = f.readline().split()[:2]
        sq = np.array([list(map(float, l.split()[:3])) for l in f if l.strip()])
    from fractions import Fraction
    print(signature(float(Fraction(s)), sq))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    f = sp.add_parser('frontier'); f.add_argument('--nmin', type=int, default=100); f.add_argument('--nmax', type=int, default=324)
    f.add_argument('--json')
    g = sp.add_parser('file'); g.add_argument('path')
    a = ap.parse_args()
    {'frontier': cmd_frontier, 'file': cmd_file}[a.cmd](a)
