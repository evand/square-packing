#!/usr/bin/env python3
"""Seed generator for packer: tilted squares + corner-anchored axis fill.

A seed = a set of tilted squares (given as blocks, or transplanted from a record packing) plus axis-aligned
squares placed greedily from four lattices, each flush with one corner of the container (records show exactly
this: axis regions anchored to different walls, meeting along seams).  Candidates are accepted in order of
distance from their anchor corner if they overlap nothing already placed by more than `tol` (the optimizer
repairs small overlaps).  The seed is written in packer's text format with exactly n squares when the fill
reaches n; extra axis squares are dropped from the far end of the fill order.

Usage as a library: seed(n, s, tilted, tol) -> list of (x, y, deg) or None.
"""
import math, sys, json

H = 0.5


def corners(x, y, deg):
    t = math.radians(deg)
    c, s = math.cos(t), math.sin(t)
    return [(x + c * u - s * v, y + s * u + c * v) for u, v in ((-H, -H), (H, -H), (H, H), (-H, H))]


def pen(a, b):
    """Separating-axis penetration depth of squares a, b = (x, y, deg); <= 0 when disjoint."""
    dx, dy = b[0] - a[0], b[1] - a[1]
    if dx * dx + dy * dy >= 2.0:
        return -1.0
    best = 1e9
    for own, oth in ((a, b), (b, a)):
        for k in (0, 90):
            phi = math.radians(own[2] + k)
            cp, sp = math.cos(phi), math.sin(phi)
            al = math.radians(own[2] + k - oth[2])
            o = H + H * (abs(math.cos(al)) + abs(math.sin(al))) - abs(dx * cp + dy * sp)
            if o <= 0:
                return o
            best = min(best, o)
    return best


def wall_viol(sq, s):
    t = math.radians(sq[2])
    w = H * (abs(math.cos(t)) + abs(math.sin(t)))
    return max(w - sq[0], sq[0] + w - s, w - sq[1], sq[1] + w - s)


def block(origin, deg, a, b):
    """a x b block of squares at angle deg; origin = centre of the (0,0) square."""
    t = math.radians(deg)
    u, v = (math.cos(t), math.sin(t)), (-math.sin(t), math.cos(t))
    return [(origin[0] + i * u[0] + j * v[0], origin[1] + i * u[1] + j * v[1], deg) for i in range(a) for j in range(b)]


def axis_fill(s, placed, tol=0.0, anchors=(0, 1, 2, 3)):
    """Greedy axis fill.  anchor 0 = lower-left, 1 = lower-right, 2 = upper-right, 3 = upper-left."""
    m = int(math.floor(s))
    cand = []
    for an in anchors:
        for i in range(m):
            for j in range(m):
                x = 0.5 + i if an in (0, 3) else s - 0.5 - i
                y = 0.5 + j if an in (0, 1) else s - 0.5 - j
                cand.append((max(i, j) + 0.01 * (i + j), an, (x, y, 0.0)))
    cand.sort(key=lambda c: c[0])
    out = []
    allsq = list(placed)
    for _, _, q in cand:
        if all(pen(q, p) <= tol for p in allsq):
            out.append(q)
            allsq.append(q)
    return out


def seed(n, s, tilted, tol=0.0, anchors=(0, 1, 2, 3)):
    tilted = [q for q in tilted if wall_viol(q, s) <= max(tol, 1e-9)]
    ax = axis_fill(s, tilted, tol, anchors)
    tot = len(tilted) + len(ax)
    if tot < n:
        return None, tot
    return (tilted + ax)[:n], tot


def write(path, s, sq):
    with open(path, 'w') as f:
        f.write(f"{len(sq)} {s!r}\n")
        for x, y, a in sq:
            f.write(f"{x!r} {y!r} {a % 90!r}\n")


def load_site(n):
    d = json.load(open(f'/home/evand/math/square-packing/public/site/www/data/p/square-{n}.json'))
    return float(d['s']), [tuple(map(float, q)) for q in d['squares']]


def is_axis(q, eps=1e-6):
    a = q[2] % 90
    return min(a, 90 - a) < eps


if __name__ == '__main__':
    print(__doc__)


def _line_fill(s, allsq, out, horiz, fromlow, pos, toward_low, tol, step):
    """Fill one row (horiz) or column at fixed coordinate pos, packing squares toward one end."""
    along = 0.5 if toward_low else s - 0.5
    d = 1.0 if toward_low else -1.0
    while 0.5 - 1e-12 <= along <= s - 0.5 + 1e-12:
        q = (along, pos, 0.0) if horiz else (pos, along, 0.0)
        worst = max((pen(q, p) for p in allsq), default=-1.0)
        if worst <= tol:
            out.append(q); allsq.append(q)
            along += d
        else:
            along += d * max(step, min(worst, 1.0) * 0.5)
    return out


def sweep_fill(s, placed, tol=0.0, order='rows', step=0.004):
    """Axis fill by lines.  order: 'rows' (rows from bottom and top alternately), 'cols', or 'mixed'.
    Lines alternate between the two anchor walls; within each line squares pack toward both ends
    (first toward the low end, then a second pass toward the high end fills the remaining gaps)."""
    m = int(math.floor(s))
    allsq = list(placed)
    out = []
    lines = []
    for j in range(m):
        for fromlow in (True, False):
            pos = 0.5 + j if fromlow else s - 0.5 - j
            if order in ('rows', 'mixed'):
                lines.append((j, True, fromlow, pos))
            if order in ('cols', 'mixed'):
                lines.append((j, False, fromlow, pos))
    for j, horiz, fromlow, pos in lines:
        _line_fill(s, allsq, out, horiz, fromlow, pos, True, tol, step)
        _line_fill(s, allsq, out, horiz, fromlow, pos, False, tol, step)
    return out


def best_seed(n, s, tilted, tol=0.0):
    """Try the fill variants; return (squares or None, best total, variant)."""
    tilted = [q for q in tilted if wall_viol(q, s) <= max(tol, 1e-9)]
    best = (None, -1, None)
    for v in ('corner', 'rows', 'cols', 'mixed'):
        ax = axis_fill(s, tilted, tol) if v == 'corner' else sweep_fill(s, tilted, tol, v)
        tot = len(tilted) + len(ax)
        if tot > best[1]:
            best = ((tilted + ax)[:n] if tot >= n else None, tot, v)
    return best


def candidates(s, placed, tol, extra_offsets=()):
    """Axis candidates: the four corner lattices (+ optional extra x/y offsets), minus those hitting placed squares."""
    m = int(math.floor(s))
    xs = sorted(set([0.5 + i for i in range(m)] + [s - 0.5 - i for i in range(m)] +
                    [o + i for o in extra_offsets for i in range(m + 1) if 0.5 <= o + i <= s - 0.5]))
    out = []
    for x in xs:
        for y in xs:
            q = (x, y, 0.0)
            if all(pen(q, p) <= tol for p in placed):
                out.append(q)
    return out


def mis_fill(s, placed, tol=1e-6, extra_offsets=(), time_limit=20):
    """Maximum set of mutually non-overlapping axis candidates (exact MILP, HiGHS)."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import lil_matrix
    C = candidates(s, placed, tol, extra_offsets)
    N = len(C)
    # axis-aligned unit squares overlap iff |dx| < 1 and |dy| < 1 (beyond tol)
    edges = [(a, b) for a in range(N) for b in range(a + 1, N)
             if abs(C[a][0] - C[b][0]) < 1 - tol and abs(C[a][1] - C[b][1]) < 1 - tol]
    if not edges:
        return C
    A = lil_matrix((len(edges), N))
    for r, (a, b) in enumerate(edges):
        A[r, a] = 1; A[r, b] = 1
    res = milp(c=-np.ones(N), constraints=LinearConstraint(A.tocsr(), -np.inf, 1),
               integrality=np.ones(N), bounds=Bounds(0, 1), options={'time_limit': time_limit})
    x = res.x if res.x is not None else np.zeros(N)
    return [C[i] for i in range(N) if x[i] > 0.5]


def _clip_extent(poly, lo, hi, axis):
    """Min/max of the other coordinate over poly ∩ {lo <= p[axis] <= hi} (None if empty)."""
    def clip(P, keep):
        out = []
        for k in range(len(P)):
            a, b = P[k], P[(k + 1) % len(P)]
            ia, ib = keep(a), keep(b)
            if ia:
                out.append(a)
            if ia != ib:
                fa = a[axis] - (lo if (a[axis] < lo) != (b[axis] < lo) else hi)
                fb = b[axis] - (lo if (a[axis] < lo) != (b[axis] < lo) else hi)
                t = fa / (fa - fb) if fa != fb else 0.0
                out.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
        return out
    P = clip(poly, lambda p: p[axis] >= lo)
    if P:
        P = clip(P, lambda p: p[axis] <= hi)
    if not P:
        return None
    o = 1 - axis
    return min(p[o] for p in P), max(p[o] for p in P)


def candidates2(s, placed, tol=1e-6, reach=4):
    """Corner-lattice candidates plus staircase candidates: on each row/column line anchored to a wall, axis squares
    touching a placed (tilted) square from either side, continued up to `reach` steps away from it."""
    m = int(math.floor(s))
    lat = [0.5 + i for i in range(m)] + [s - 0.5 - i for i in range(m)]
    C = set()
    for x in lat:
        for y in lat:
            C.add((round(x, 12), round(y, 12)))
    polys = [corners(*q) for q in placed]
    for axis in (1, 0):            # axis=1: rows (fixed y), axis=0: columns (fixed x)
        for c in lat:
            for P in polys:
                ext = _clip_extent(P, c - 0.5, c + 0.5, axis)
                if ext is None:
                    continue
                for start, d in ((ext[0] - 0.5, -1.0), (ext[1] + 0.5, 1.0)):
                    for i in range(reach):
                        u = start + d * i
                        if 0.5 - 1e-9 <= u <= s - 0.5 + 1e-9:
                            u = min(max(u, 0.5), s - 0.5)
                            C.add((round(u, 12), round(c, 12)) if axis == 1 else (round(c, 12), round(u, 12)))
    near = {}
    for p in placed:
        near.setdefault((int(math.floor(p[0])), int(math.floor(p[1]))), []).append(p)
    out = []
    for x, y in C:
        q = (x, y, 0.0)
        bx, by = int(math.floor(x)), int(math.floor(y))
        if all(pen(q, p) <= tol for ux in (bx - 2, bx - 1, bx, bx + 1, bx + 2) for uy in (by - 2, by - 1, by, by + 1, by + 2)
               for p in near.get((ux, uy), ())):
            out.append(q)
    return out


def mis_fill2(s, placed, tol=1e-6, reach=4, time_limit=30, axtol=1e-6):
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import coo_matrix
    C = candidates2(s, placed, tol, reach)
    N = len(C)
    buckets = {}
    for k, (x, y, _) in enumerate(C):
        buckets.setdefault((int(x), int(y)), []).append(k)
    rows, cols, r = [], [], 0
    for k, (x, y, _) in enumerate(C):
        bx, by = int(x), int(y)
        for ux in (bx - 1, bx, bx + 1):
            for uy in (by - 1, by, by + 1):
                for l in buckets.get((ux, uy), ()):
                    if l > k and abs(C[l][0] - x) < 1 - axtol and abs(C[l][1] - y) < 1 - axtol:
                        rows += [r, r]; cols += [k, l]; r += 1
    if r == 0:
        return C
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(r, N)).tocsr()
    res = milp(c=-np.ones(N), constraints=LinearConstraint(A, -np.inf, 1), integrality=np.ones(N),
               bounds=Bounds(0, 1), options={'time_limit': time_limit})
    x = res.x if res.x is not None else np.zeros(N)
    return [C[i] for i in range(N) if x[i] > 0.5]
