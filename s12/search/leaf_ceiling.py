#!/usr/bin/env python3
"""leaf_ceiling.py -- EXACT certification of a finite fractional packing measure in a level-2 leaf.

The no-go instrument for the t = 4 endgame.  Given a finite measure `mu` on poses of CLOSED unit
squares inside the CLOSED container `[0,t]^2` (t rational), it certifies, with every load-bearing
comparison in Python integers / `Fraction`s:

  (a) COVERAGE      cov(x) = mu({S : x in S}) <= 1 at EVERY point x of the container.
                    `cov` is piecewise constant on the arrangement of the closed squares and its
                    maximum is attained at an arrangement vertex (square corner or crossing of two
                    closed edges) -- the argument of `search/DUAL_EXACT.md`, reproduced in
                    `search/LEAF_CEILING.md` 1.  All vertices are enumerated in exact integers.

  (b) ANCHOR CLIQUES   mu(K(p, A)) <= 1 for every clique of the certifiable family
                    K(p, A) = {S : p in S and S meets A} u {S : A subset S}
                    (`notes/clique-family.md` Lemma 0), p any point, A any point or closed segment.
                    Three instruments, all exact (see `search/LEAF_CEILING.md` 2):
                      `exhaust`  the COMPLETE finite candidate set of the theorem in LEAF_CEILING.md
                                 2.3 (p an arrangement vertex; the anchor's supporting line through
                                 two arrangement vertices, both inside a common support square; its
                                 endpoints at the interval breakpoints on that line), with two
                                 sound prunes -- the exact max-weight clique through the primary
                                 container, and the clique ceiling of `clique`.  Decisive when it
                                 finishes inside the budget; on the t = 4 leaf measures it does
                                 not, and it reports `complete: false` (LEAF_CEILING.md 2.5).
                      `clique`   the exact maximum-weight clique of the closed-intersection graph
                                 (integer weights, branch and bound).  Every K(p, A) is such a
                                 clique, so this is a SOUND OVER-APPROXIMATION: `<= 1` proves (b).
                      `scan`     floats choose candidates, exact rationals evaluate them, and an
                                 exact local ascent over the candidate set of 2.3 polishes them.
                                 Gives a certified LOWER bound (a witness), never a proof of (b).

  (c) REGION MASSES   the mass centred in each closed region of a leaf (four corner boxes, eight
                    wall slots, interior; `search/level2_regions.py`) is exactly the leaf's count,
                    with a pose whose CENTRE lies on a region boundary free to be assigned to
                    either adjacent region -- the packing's choice, which is what makes the leaf a
                    partition.  Optionally the chord inequality mu(wall strip) <= 3.

A measure certified with (a), (b), (c) and mass >= 12 at side t is a THEOREM-GRADE NO-GO: by weak
LP duality no certificate built from point cliques, anchor cliques of this family and per-region
multipliers can have weight < 12 at any container side >= t, so that leaf does not close.

Measure file format (all exact):

    # t = 4/1                  (or a `t 4/1` line; `--t` overrides)
    # sym 8                    (8 = D4-symmetrised: mass/8 on each of the 8 dihedral images;
    #                           1 or absent = flat, one square per line with its own mass)
    pose p q cx cy mass        theta = 2 arctan(p/q), cos = (q^2-p^2)/(q^2+p^2), sin = 2pq/(q^2+p^2)

`search/dual_exact.py`'s and `search/clique_ceiling.py`'s support files are read directly (they
declare "D4-symmetrised" in the header).

Usage
    python3 search/leaf_ceiling.py check FILE [--t T] [--corners 1111] [--patterns 01010101]
            [--r 1] [--chord] [--anchor scan|exhaust|clique|all|none] [--procs 4] [--budget 3600]
    python3 search/leaf_ceiling.py snap FLOATFILE --t T --out FILE [--Q 100000] [--Dc 1000000]
            [--DM 1000000000] [--polish] [--corners 1111] [--patterns ........]
    python3 search/leaf_ceiling.py selftest
"""
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr
from math import gcd, lcm

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')

LOG = None


def log(msg=''):
    print(msg, flush=True)
    if LOG:
        LOG.write(str(msg) + '\n')
        LOG.flush()


# ============================================================================ exact squares
# A square record is the tuple
#     (CX, CY, Dk, a, b, r, corners, Dq, mass_int)
# with centre (CX/Dk, CY/Dk), cos = a/r, sin = b/r, and the four corners (X/Dq, Y/Dq) in cyclic
# order.  `mass_int` is the integer numerator of the mass over the global denominator DM.

def cos_sin(p, q):
    return q * q - p * p, 2 * p * q, q * q + p * p


def half_width(a, b, r):
    return Fr(abs(a) + abs(b), 2 * r)


def make_square(cx, cy, p, q, t, mass_int=0, tag=None):
    """exact square record; asserts the closed square lies in the closed container [0,t]^2"""
    a, b, r = cos_sin(p, q)
    Dk = lcm(cx.denominator, cy.denominator)
    CX = cx.numerator * (Dk // cx.denominator)
    CY = cy.numerator * (Dk // cy.denominator)
    Dq = 2 * r * Dk
    corners = []
    for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        X = 2 * r * CX + Dk * (a * sx - b * sy)
        Y = 2 * r * CY + Dk * (b * sx + a * sy)
        corners.append((X, Y))
    g = Dq
    for X, Y in corners:
        g = gcd(g, gcd(abs(X), abs(Y)))
    corners = [(X // g, Y // g) for X, Y in corners]
    Dq //= g
    tn, td = t.numerator, t.denominator
    for X, Y in corners:
        assert 0 <= X and td * X <= tn * Dq and 0 <= Y and td * Y <= tn * Dq, \
            ('square not admissible', tag, cx, cy, p, q)
    return (CX, CY, Dk, a, b, r, tuple(corners), Dq, mass_int)


def sq_bbox(s):
    """(Xlo, Xhi, Ylo, Yhi, Dq) exact bounding box of the square"""
    xs = [X for X, Y in s[6]]
    ys = [Y for X, Y in s[6]]
    return min(xs), max(xs), min(ys), max(ys), s[7]


def sq_contains(s, X, Y, D):
    """exact: is the point (X/D, Y/D) in the CLOSED square s?"""
    CX, CY, Dk, a, b, r = s[:6]
    dx = X * Dk - CX * D
    dy = Y * Dk - CY * D
    lim = r * D * Dk
    u = a * dx + b * dy
    if u < 0:
        u = -u
    if 2 * u > lim:
        return False
    v = a * dy - b * dx
    if v < 0:
        v = -v
    return 2 * v <= lim


def sq_meets_segment(s, P0, P1, Ds):
    """exact: does the CLOSED square s meet the closed segment (P0/Ds)-(P1/Ds)?

    Separating-axis test for two convex polygons: it suffices to test the edge normals of both.
    The square gives (a, b) and (-b, a); the segment, as a degenerate polygon with the two edges
    P0->P1 and P1->P0, gives its own normal n = (-(dy), dx) -- the axis whose omission is the near
    miss of `notes/clique-family.md` 7.  Three axes, all in integers."""
    corners, Dq = s[6], s[7]
    a, b = s[3], s[4]
    x0, y0 = P0
    x1, y1 = P1
    axes = ((a, b), (-b, a), (-(y1 - y0), x1 - x0))
    for ux, uy in axes:
        if ux == 0 and uy == 0:          # degenerate segment: the third axis vanishes
            continue
        d = [ux * X + uy * Y for X, Y in corners]
        lo1, hi1 = min(d), max(d)
        e0 = ux * x0 + uy * y0
        e1 = ux * x1 + uy * y1
        lo2, hi2 = (e0, e1) if e0 <= e1 else (e1, e0)
        # square projection is [lo1/Dq, hi1/Dq]; segment's is [lo2/Ds, hi2/Ds]
        if hi1 * Ds < lo2 * Dq or hi2 * Dq < lo1 * Ds:
            return False
    return True


def sq_meets_sq(si, sj):
    """exact: do two closed squares meet?  (4 edge normals, integers only)"""
    for s in (si, sj):
        a, b = s[3], s[4]
        for ux, uy in ((a, b), (-b, a)):
            di = [ux * X + uy * Y for X, Y in si[6]]
            dj = [ux * X + uy * Y for X, Y in sj[6]]
            if max(di) * sj[7] < min(dj) * si[7] or max(dj) * si[7] < min(di) * sj[7]:
                return False
    return True


# ============================================================================ measure files
def _parse_t(text):
    for tok in text.replace(',', ' ').replace(';', ' ').split():
        try:
            return Fr(tok)                      # accepts "399/100", "4", "3.99"
        except (ValueError, ZeroDivisionError):
            pass
    return None


def read_measure(path, t_override=None):
    """-> (t, sym, [(p, q, cx, cy, mass)])  with exact Fraction cx, cy, mass"""
    t = None
    sym = None
    poses = []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith('#'):
            low = s.lower()
            if 'd4-symmetrised' in low or 'd4 symmetrised' in low:
                sym = 8
            if t is None and 't =' in s:
                t = _parse_t(s.split('t =')[1].split()[0])
            if t is None and s.lstrip('# ').startswith('t '):
                t = _parse_t(s.lstrip('# ')[2:])
            if 'sym' in low:
                for i, tok in enumerate(low.split()):
                    if tok == 'sym' and i + 1 < len(low.split()):
                        try:
                            sym = int(low.split()[i + 1])
                        except ValueError:
                            pass
            continue
        q = s.split()
        if q[0] == 't':
            t = Fr(q[1]) if '/' in q[1] else Fr(int(q[1]))
        elif q[0] == 'sym':
            sym = int(q[1])
        elif q[0] == 'pose' and len(q) == 6:
            poses.append((int(q[1]), int(q[2]), Fr(q[3]), Fr(q[4]), Fr(q[5])))
        elif q[0] == 'pose':
            raise SystemExit(f'{path}: `pose` line with {len(q)-1} fields; this reader needs the '
                             f'exact form `pose p q cx cy mass` (use `leaf_ceiling.py snap` first)')
    if t_override is not None:
        t = t_override
    if t is None:
        raise SystemExit(f'{path}: no container side; pass --t')
    if sym is None:
        sym = 1
    return t, sym, poses


def d4_images(cx, cy, p, q, t):
    """the 8 dihedral images, `packing_dual.py` / `dual_exact.py` convention (-theta <-> -p)"""
    return [(cx, cy, p, q), (t - cx, cy, -p, q), (cx, t - cy, -p, q), (t - cx, t - cy, p, q),
            (cy, cx, -p, q), (t - cy, cx, p, q), (cy, t - cx, p, q), (t - cy, t - cx, -p, q)]


def build_squares(t, sym, poses, DM=None):
    """-> (squares, DM, pose_of).  Masses become integers over the common denominator DM;
    with sym = 8 each pose contributes 8 images of mass mass/8 (so DM absorbs the 8)."""
    if DM is None:
        DM = 1
        for (_, _, _, _, m) in poses:
            DM = lcm(DM, m.denominator)
        if sym == 8:
            DM *= 8
    squares = []
    pose_of = []
    for k, (p, q, cx, cy, m) in enumerate(poses):
        if m == 0:
            continue
        if sym == 8:
            share = m / 8
            ims = d4_images(cx, cy, p, q, t)
        else:
            share = m
            ims = [(cx, cy, p, q)]
        num = share * DM
        assert num.denominator == 1, f'mass {m} is not a multiple of 1/{DM}'
        num = int(num)
        for (x, y, pp, qq) in ims:
            squares.append(make_square(x, y, pp, qq, t, num, tag=k))
            pose_of.append(k)
    return squares, DM, pose_of


# ============================================================================ arrangement vertices
def reduce3(X, Y, D):
    g = gcd(gcd(abs(X), abs(Y)), D)
    return (X // g, Y // g, D // g)


def seg_intersect(e1, e2):
    """intersection point of two closed segments with rational endpoints, or None"""
    x1, y1, x2, y2, D1 = e1
    x3, y3, x4, y4, D2 = e2
    ax = x1 * D2
    ay = y1 * D2
    rx = x2 * D2 - ax
    ry = y2 * D2 - ay
    cx = x3 * D1
    cy = y3 * D1
    sx = x4 * D1 - cx
    sy = y4 * D1 - cy
    det = rx * sy - ry * sx
    if det == 0:
        return None
    qx = cx - ax
    qy = cy - ay
    tn = qx * sy - qy * sx
    un = qx * ry - qy * rx
    if det < 0:
        det = -det
        tn = -tn
        un = -un
    if tn < 0 or tn > det or un < 0 or un > det:
        return None
    return reduce3(ax * det + tn * rx, ay * det + tn * ry, D1 * D2 * det)


G_EDGE = 5          # grid resolution of the exact edge-pair prefilter
G_SQ = 10           # grid resolution of the exact vertex -> square prefilter


# globals shared with forked workers (the `dual_exact.py` pattern)
_EDGES = []
_EMIN = []
_EBUCKET = {}
_SQ = []
_SBUCKET = {}
_VERTS = []


def _pair_worker(cells):
    out = set()
    for (cxi, cyi) in cells:
        L = _EBUCKET.get((cxi, cyi))
        if not L:
            continue
        n = len(L)
        for i in range(n):
            ei = L[i]
            mxi, myi = _EMIN[ei]
            si = ei >> 2
            e1 = _EDGES[ei]
            for j in range(i + 1, n):
                ej = L[j]
                if (ej >> 2) == si:
                    continue
                mxj, myj = _EMIN[ej]
                if max(mxi, mxj) != cxi or max(myi, myj) != cyi:
                    continue                       # each pair is handled exactly once
                v = seg_intersect(e1, _EDGES[ej])
                if v is not None:
                    out.add(v)
    return out


def _cov_worker(rng):
    lo, hi = rng
    pairs = []
    for vi in range(lo, hi):
        X, Y, D = _VERTS[vi]
        for si in _SBUCKET.get(((X * G_SQ) // D, (Y * G_SQ) // D), ()):
            if sq_contains(_SQ[si], X, Y, D):
                pairs.append((si, vi))
    return pairs


def _chunks(lst, n):
    k = max(1, (len(lst) + n - 1) // n)
    return [lst[i:i + k] for i in range(0, len(lst), k)]


def enumerate_vertices(squares, t, verbose=True, procs=1):
    """every arrangement vertex of the closed squares, as reduced integer triples (X, Y, D).

    Complete by the argument of LEAF_CEILING.md 1: the maximum of any subcollection's coverage is
    attained at a corner of one of its squares or at a crossing of two closed edges."""
    global _EDGES, _EMIN, _EBUCKET, _VERTS
    t0 = time.time()
    _EDGES = []
    _EMIN = []
    _EBUCKET = {}
    for si, s in enumerate(squares):
        corners, Dq = s[6], s[7]
        for k in range(4):
            (x1, y1), (x2, y2) = corners[k], corners[(k + 1) % 4]
            ei = len(_EDGES)
            _EDGES.append((x1, y1, x2, y2, Dq))
            cx0 = (min(x1, x2) * G_EDGE) // Dq
            cx1 = (max(x1, x2) * G_EDGE) // Dq
            cy0 = (min(y1, y2) * G_EDGE) // Dq
            cy1 = (max(y1, y2) * G_EDGE) // Dq
            _EMIN.append((cx0, cy0))
            for cx in range(cx0, cx1 + 1):
                for cy in range(cy0, cy1 + 1):
                    _EBUCKET.setdefault((cx, cy), []).append(ei)
    cells = sorted(_EBUCKET.keys())
    if procs > 1:
        from multiprocessing import Pool
        with Pool(procs) as pool:
            parts = pool.map(_pair_worker, _chunks(cells, 4 * procs))
        verts = set()
        for p in parts:
            verts |= p
    else:
        verts = _pair_worker(cells)
    n_int = len(verts)
    for s in squares:
        for X, Y in s[6]:
            verts.add(reduce3(X, Y, s[7]))
    V = sorted(verts)
    _VERTS = V
    if verbose:
        log(f"  vertices: {n_int} distinct edge crossings + corners -> {len(V)} distinct "
            f"({len(_EDGES)} edges, {len(cells)} cells, {time.time() - t0:.1f} s)")
    return V


class Incidence:
    """for each square k, the (int32) indices of the arrangement vertices it CONTAINS -- the exact
    incidence relation, stored by column so that the coverage of an arbitrary SUBCOLLECTION costs
    only that subcollection's incidences.

    `cov(subset)` returns an exact int64 array over the vertices; the masses are integers over the
    global denominator DM, so every entry is exact."""

    def __init__(self, nverts, cols):
        self.nv = nverts
        self.cols = cols                            # list of np.int32 arrays
        self.nnz = int(sum(len(c) for c in cols))

    def cov(self, w, subset=None, out=None):
        acc = np.zeros(self.nv, dtype=np.int64) if out is None else out
        if out is not None:
            acc.fill(0)
        it = range(len(self.cols)) if subset is None else subset
        for k in it:
            if w[k]:
                acc[self.cols[k]] += w[k]
        return acc

    def to_csr(self, rows=None):
        """the same relation as a scipy CSR matrix of floats, restricted to the vertex indices
        `rows` (used only to CHOOSE masses; every number is re-checked exactly afterwards)"""
        import scipy.sparse as sp
        if rows is None:
            sel = None
            nr = self.nv
        else:
            sel = np.full(self.nv, -1, dtype=np.int64)
            sel[rows] = np.arange(len(rows))
            nr = len(rows)
        R = []
        C = []
        for k, c in enumerate(self.cols):
            if not len(c):
                continue
            if sel is None:
                r = c
            else:
                r = sel[c]
                r = r[r >= 0]
                if not len(r):
                    continue
            R.append(r)
            C.append(np.full(len(r), k, dtype=np.int32))
        if not R:
            return sp.csr_matrix((nr, len(self.cols)))
        R = np.concatenate(R)
        C = np.concatenate(C)
        return sp.csr_matrix((np.ones(len(R)), (R, C)), shape=(nr, len(self.cols)))


def incidences(squares, V, verbose=True, procs=1):
    """the exact incidence relation as an `Incidence` (vertex lists per square)"""
    global _SQ, _SBUCKET, _VERTS
    t0 = time.time()
    _SQ = squares
    _VERTS = V
    _SBUCKET = {}
    for si, s in enumerate(squares):
        Xlo, Xhi, Ylo, Yhi, Dq = sq_bbox(s)
        for cx in range((Xlo * G_SQ) // Dq, (Xhi * G_SQ) // Dq + 1):
            for cy in range((Ylo * G_SQ) // Dq, (Yhi * G_SQ) // Dq + 1):
                _SBUCKET.setdefault((cx, cy), []).append(si)
    n = len(V)
    rngs = [(i, min(i + 4000, n)) for i in range(0, n, 4000)]
    cols = [[] for _ in squares]
    if procs > 1:
        from multiprocessing import Pool
        with Pool(procs) as pool:
            for pairs in pool.imap_unordered(_cov_worker, rngs):
                for (si, vi) in pairs:
                    cols[si].append(vi)
    else:
        for rng in rngs:
            for (si, vi) in _cov_worker(rng):
                cols[si].append(vi)
    cols = [np.array(sorted(c), dtype=np.int32) for c in cols]
    Inc = Incidence(n, cols)
    if verbose:
        log(f"  incidences: {n} vertices, {Inc.nnz} (vertex, square) pairs, "
            f"{time.time() - t0:.1f} s")
    return Inc


def max_coverage(Inc, squares, DM):
    """exact maximum of cov over the container, as a Fraction, plus the argmax vertex index"""
    w = np.array([s[8] for s in squares], dtype=np.int64)
    assert w.sum() * 1 < 2 ** 62, 'integer masses too large for int64 accumulation'
    cov = Inc.cov(w)
    i = int(np.argmax(cov))
    return Fr(int(cov[i]), DM), i, cov


# ============================================================================ (c) region masses
def region_boxes(t, r):
    """the 13 closed regions of a level-2 leaf as (label, xlo, xhi, ylo, yhi), exact Fractions.
    Numbering follows `search/level2_regions.classify`."""
    h = t / 2
    lo, hi = Fr(r), t - Fr(r)
    B = []
    B.append(('C0', Fr(0), lo, Fr(0), lo))
    B.append(('C1', hi, t, Fr(0), lo))
    B.append(('C2', Fr(0), lo, hi, t))
    B.append(('C3', hi, t, hi, t))
    B.append(('W0', lo, h, Fr(0), lo))
    B.append(('W1', h, hi, Fr(0), lo))
    B.append(('W2', hi, t, lo, h))
    B.append(('W3', hi, t, h, hi))
    B.append(('W4', h, hi, hi, t))
    B.append(('W5', lo, h, hi, t))
    B.append(('W6', Fr(0), lo, h, hi))
    B.append(('W7', Fr(0), lo, lo, h))
    B.append(('I', lo, hi, lo, hi))
    return B


def region_targets(t, args):
    """{label: Fraction} from --corners / --patterns ('.' = the region is free)"""
    cs = getattr(args, 'corners', None)
    ps = getattr(args, 'patterns', None)
    if not cs and not ps:
        return None
    out = {}
    for i, ch in enumerate(cs or '....'):
        if ch.isdigit():
            out[f'C{i}'] = Fr(int(ch))
    for j, ch in enumerate(ps or '........'):
        if ch.isdigit():
            out[f'W{j}'] = Fr(int(ch))
    return out or None


def regions_of(cx, cy, boxes):
    """every closed region whose box contains the centre (exact).  More than one <=> the centre
    lies on a region boundary and the packing may declare it in either."""
    return [lab for (lab, x0, x1, y0, y1) in boxes if x0 <= cx <= x1 and y0 <= cy <= y1]


def region_report(squares_meta, boxes, target, DM):
    """squares_meta: [(cx, cy, mass_int)].  target: {label: Fraction} for the labels that are
    pinned.  Returns (ok, assigned_masses, ambiguous, note)."""
    forced = {lab: 0 for (lab, *_ ) in boxes}
    amb = []
    for (cx, cy, m) in squares_meta:
        labs = regions_of(cx, cy, boxes)
        assert labs, f'centre ({cx}, {cy}) is in no region'
        if len(labs) == 1:
            forced[labs[0]] += m
        else:
            amb.append((m, labs, cx, cy))
    if target is None:
        return True, forced, amb, 'no region equalities requested'
    tgt = {lab: int(v * DM) for lab, v in target.items()}
    for lab, v in target.items():
        assert Fr(tgt[lab], DM) == v, 'region target is not a multiple of 1/DM'
    # exact search over the assignment of the ambiguous poses (few, and each has <= 4 choices)
    labs_pinned = sorted(tgt)
    order = sorted(range(len(amb)), key=lambda i: -amb[i][0])
    best = [None]
    rest = [0] * (len(amb) + 1)
    for i in range(len(amb) - 1, -1, -1):
        rest[i] = rest[i + 1] + amb[order[i]][0]

    def rec(i, cur):
        if best[0] is not None:
            return
        if i == len(amb):
            if all(cur.get(lab, 0) == tgt[lab] for lab in labs_pinned):
                best[0] = dict(cur)
            return
        deficit = sum(max(0, tgt[lab] - cur.get(lab, 0)) for lab in labs_pinned)
        if deficit > rest[i]:
            return
        for lab in labs_pinned + [None]:
            m, labs, _, _ = amb[order[i]]
            if lab is None:
                for l2 in labs:
                    if l2 not in tgt:
                        cur[l2] = cur.get(l2, 0) + m
                        rec(i + 1, cur)
                        cur[l2] -= m
                        if best[0] is not None:
                            return
                        break
                continue
            if lab not in labs:
                continue
            if cur.get(lab, 0) + m > tgt[lab]:
                continue
            cur[lab] = cur.get(lab, 0) + m
            rec(i + 1, cur)
            cur[lab] -= m
            if best[0] is not None:
                return

    rec(0, dict(forced))
    if best[0] is None:
        return False, forced, amb, 'no assignment of the boundary poses meets the region targets'
    return True, best[0], amb, 'assignment found'


# ============================================================================ (b) anchor cliques
def anchor_members(squares, P0, P1, Ds, near=None):
    """exact: (contains_mask, meets_mask) for the closed segment (P0/Ds)-(P1/Ds).

    `contains` = the segment lies in the square (both endpoints, by convexity);
    `meets`    = the segment intersects the square.  Note contains => meets."""
    n = len(squares)
    cont = bytearray(n)
    meet = bytearray(n)
    sx0 = min(P0[0], P1[0])
    sx1 = max(P0[0], P1[0])
    sy0 = min(P0[1], P1[1])
    sy1 = max(P0[1], P1[1])
    idx = range(n) if near is None else near
    for i in idx:
        s = squares[i]
        Xlo, Xhi, Ylo, Yhi, Dq = sq_bbox(s)
        if Xhi * Ds < sx0 * Dq or sx1 * Dq < Xlo * Ds:
            continue
        if Yhi * Ds < sy0 * Dq or sy1 * Dq < Ylo * Ds:
            continue
        if not sq_meets_segment(s, P0, P1, Ds):
            continue
        meet[i] = 1
        if sq_contains(s, P0[0], P0[1], Ds) and sq_contains(s, P1[0], P1[1], Ds):
            cont[i] = 1
    return cont, meet


def anchor_value(squares, Inc, DM, P0, P1, Ds, w=None):
    """exact max over p of mu(K(p, A)) for the fixed anchor A = [P0/Ds, P1/Ds].

    mu(K(p,A)) = mu(Cont) + mu({i in Meet \\ Cont : p in S_i}); the second term is the coverage of
    the subcollection D = Meet \\ Cont, whose maximum over the container is attained at an
    arrangement vertex of D -- a subset of the vertices of the whole arrangement (LEAF_CEILING 2.2).
    Returns (value, best vertex index, |Cont|, |Meet|)."""
    if w is None:
        w = np.array([s[8] for s in squares], dtype=np.int64)
    cont, meet = anchor_members(squares, P0, P1, Ds)
    c = np.frombuffer(bytes(cont), dtype=np.uint8).astype(bool)
    m = np.frombuffer(bytes(meet), dtype=np.uint8).astype(bool)
    d = m & ~c
    base = int(w[c].sum())
    if not d.any():
        return Fr(base, DM), -1, int(c.sum()), int(m.sum())
    cov = Inc.cov(w, np.nonzero(d)[0])
    i = int(np.argmax(cov))
    return Fr(base + int(cov[i]), DM), i, int(c.sum()), int(m.sum())


def anchor_value_at(squares, DM, P0, P1, Ds, pX, pY, pD):
    """exact mu(K(p, A)) at a GIVEN point p -- the definition, for cross-checking"""
    cont, meet = anchor_members(squares, P0, P1, Ds)
    tot = 0
    for i, s in enumerate(squares):
        if cont[i] or (meet[i] and sq_contains(s, pX, pY, pD)):
            tot += s[8]
    return Fr(tot, DM)


# ---------------------------------------------------------------- the exact 1-D problem on a line
# the exact maximum coverage of the measure under test, as an integer over DM.  `best_on_line`
# uses it as the sound upper bound  max_p cov_D(p) <= M  in its pruning; it must never be smaller
# than the true maximum, so it is set from the certified value in `certify`.
MCAP = [None]

_SQF_CACHE = {}


def square_fractions(squares):
    """per-square rational constants (centre and the two unit axes), computed once: the inner loop
    of `line_intervals` is otherwise dominated by rebuilding these `Fraction`s"""
    key = id(squares)
    got = _SQF_CACHE.get(key)
    if got is not None and len(got) == len(squares):
        return got
    out = []
    for s in squares:
        CX, CY, Dk, a, b, r = s[:6]
        out.append((Fr(CX, Dk), Fr(CY, Dk), Fr(a, r), Fr(b, r), Fr(-b, r), Fr(a, r)))
    _SQF_CACHE.clear()
    _SQF_CACHE[key] = out
    return out


def line_intervals(squares, A, B, C, idx):
    """the closed interval [lo_k, hi_k] (as Fractions of the parameter u along the line
    x = A + u*B, in exact rationals) cut out of the line by each square in `idx`.

    The line is given by a rational point A = (ax, ay, ad) and a rational direction B = (bx, by).
    Returns [(k, lo, hi)] for the squares the line actually meets, and C is unused (kept for
    signature symmetry).  Exact: each square is the intersection of four half planes
    |a dx + b dy| <= r/2 * ..., clipped by Liang-Barsky in Fractions."""
    ax, ay, ad = A
    bx, by = B
    nx0, ny0 = -by, bx                      # a normal of the line (exact, not normalised)
    axF = Fr(ax, ad)
    ayF = Fr(ay, ad)
    SQF = square_fractions(squares)
    out = []
    for k in idx:
        cxF, cyF, e1x, e1y, e2x, e2y = SQF[k]
        # point on the line: (ax/ad + u bx, ay/ad + u by); frame coordinates relative to the centre
        ex = axF - cxF
        ey = ayF - cyF
        # EXACT one-line prefilter: the line meets the square iff the square's projection on the
        # line's normal straddles the line, i.e. |n.(c - A)| <= (|n.e1| + |n.e2|)/2 with e1, e2 the
        # square's unit axes.  A handful of rational operations instead of the full clip.
        d0 = nx0 * ex + ny0 * ey
        u1 = nx0 * e1x + ny0 * e1y
        u2 = nx0 * e2x + ny0 * e2y
        hw = (abs(u1) + abs(u2)) / 2
        if (d0 if d0 >= 0 else -d0) > hw:
            continue
        lo, hi = None, None
        ok = True
        for (nx, ny) in ((e1x, e1y), (e2x, e2y)):
            c0 = nx * ex + ny * ey
            c1 = nx * bx + ny * by
            for sgn in (1, -1):
                # sgn*(c0 + u c1) <= 1/2
                A_ = sgn * c1
                B_ = Fr(1, 2) - sgn * c0
                if A_ == 0:
                    if B_ < 0:
                        ok = False
                        break
                elif A_ > 0:
                    u = B_ / A_
                    hi = u if hi is None else min(hi, u)
                else:
                    u = B_ / A_
                    lo = u if lo is None else max(lo, u)
            if not ok:
                break
        if not ok or lo is None or hi is None or lo > hi:
            continue
        out.append((k, lo, hi))
    return out


def best_on_line(squares, Inc, DM, w, A, B, idx, incumbent=None, Mcap=None, maxeval=None):
    """exact maximum of mu(K(p, A')) over all segments A' on the line and all points p.

    The 1-D problem (LEAF_CEILING.md 2.2).  With I_k = [lo_k, hi_k] the interval the line cuts out
    of S_k, a segment [s, u] has
        Cont = {k : lo_k <= s, u <= hi_k},   Meet = {k : lo_k <= u, s <= hi_k},
    and the theorem says it is enough to take s and u among the finitely many lo_k, hi_k.  Ranking
    those breakpoints turns both conditions into comparisons of RANKS, so the mass of Cont and of
    Meet for ALL pairs (s, u) at once is a two-dimensional prefix sum in exact integers; the
    expensive step -- the exact maximum over p, a coverage evaluation of the subcollection
    Meet \\ Cont at every arrangement vertex -- then runs only on the pairs whose upper bound
    mu(Cont) + min(M, mu(Meet) - mu(Cont)) beats the incumbent.

    Returns (value, (s, u), cont_mask, meet_mask) or None."""
    iv = line_intervals(squares, A, B, None, idx)
    if not iv:
        return None
    pts = sorted(set([x for (_, lo, hi) in iv for x in (lo, hi)]))
    rank = {x: i for i, x in enumerate(pts)}
    P = len(pts)
    ks = np.array([k for (k, _, _) in iv], dtype=np.int64)
    lo_r = np.array([rank[lo] for (_, lo, _) in iv], dtype=np.int64)
    hi_r = np.array([rank[hi] for (_, _, hi) in iv], dtype=np.int64)
    mk = np.array([squares[k][8] for k in ks], dtype=np.int64)
    Amat = np.zeros((P, P), dtype=np.int64)          # mass by (rank of lo, rank of hi)
    np.add.at(Amat, (lo_r, hi_r), mk)
    cont = np.cumsum(Amat, axis=0)                   # a <= i
    cont = np.cumsum(cont[:, ::-1], axis=1)[:, ::-1]  # b >= j;  cont[i, j] = mu(Cont(pts[i], pts[j]))
    meet = cont.T                                    # meet[i, j] = cont[j, i] = mu(Meet)
    # cap = an upper bound on max_p cov_D(p); without a certified M, fall back to no cap
    cap = (int(MCAP[0]) if MCAP[0] is not None else (1 << 62)) if Mcap is None else int(Mcap * DM)
    iu = np.triu_indices(P)                          # i <= j
    c = cont[iu]
    m = meet[iu]
    ub = c + np.minimum(cap, m - c)
    inc = -1 if incumbent is None else int(incumbent * DM)
    sel = np.nonzero(ub > inc)[0]
    if not len(sel):
        return None
    sel = sel[np.argsort(-ub[sel])]
    best = None
    seen = set()
    for si in sel:
        if best is not None and ub[si] <= int(best[0] * DM):
            break
        if maxeval is not None and len(seen) >= maxeval:
            break                      # heuristic cap: used only by the ascent (a lower bound)
        i = int(iu[0][si])
        j = int(iu[1][si])
        cm = np.zeros(len(squares), dtype=bool)
        mm = np.zeros(len(squares), dtype=bool)
        cm[ks[(lo_r <= i) & (j <= hi_r)]] = True
        mm[ks[(lo_r <= j) & (i <= hi_r)]] = True
        key = (cm.tobytes(), mm.tobytes())
        if key in seen:
            continue
        seen.add(key)
        d = mm & ~cm
        val = int(c[si])
        if d.any():
            val += int(Inc.cov(w, np.nonzero(d)[0]).max())
        v = Fr(val, DM)
        if best is None or v > best[0]:
            best = (v, (pts[i], pts[j]), cm, mm)
    return best


# ---------------------------------------------------------------- exact maximum-weight clique
def max_clique_int(nbr, w, time_limit, lb=0):
    """maximum-weight clique with INTEGER weights, exact bound arithmetic (branch and bound).
    Returns (weight, members, complete).  Same algorithm as `search/clique_ceiling.py`.

    `lb` seeds the incumbent: the search then only looks for cliques of weight > lb, and a return
    value of exactly `lb` means "no clique heavier than lb exists" (a proof, when complete)."""
    n = len(w)
    best = [lb, []]
    t0 = time.time()
    timed_out = [False]

    def colour_bound(cands):
        classes = []
        cw = 0
        for v in sorted(cands, key=lambda v: -w[v]):
            for cl in classes:
                if not (nbr[v] & cl):
                    cl.add(v)
                    break
            else:
                classes.append({v})
                cw += w[v]
        return cw

    def expand(cands, cur, curw):
        if time.time() - t0 > time_limit:
            timed_out[0] = True
            return
        if not cands:
            if curw > best[0]:
                best[0] = curw
                best[1] = list(cur)
            return
        if curw + colour_bound(cands) <= best[0]:
            return
        cl = sorted(cands, key=lambda v: -w[v])
        while cl:
            v = cl.pop(0)
            expand([u for u in cl if u in nbr[v]], cur + [v], curw + w[v])
            if timed_out[0]:
                return
            if curw + colour_bound(cl) <= best[0]:
                return

    expand([v for v in range(n) if w[v] > 0], [], 0)
    return best[0], best[1], not timed_out[0]


def closed_graph(squares, verbose=True):
    """exact closed-intersection graph of the squares (integers only)"""
    t0 = time.time()
    n = len(squares)
    C = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2]))] for s in squares])
    nbr = [set() for _ in range(n)]
    npair = 0
    for i in range(n):
        d2 = (C[:, 0] - C[i, 0]) ** 2 + (C[:, 1] - C[i, 1]) ** 2
        for j in np.nonzero((d2 <= 2.0 + 1e-6) & (np.arange(n) > i))[0]:
            npair += 1
            if sq_meets_sq(squares[i], squares[int(j)]):
                nbr[i].add(int(j))
                nbr[int(j)].add(i)
    if verbose:
        log(f"  closed-intersection graph: {n} squares, {sum(map(len, nbr)) // 2} edges "
            f"({npair} candidate pairs, {time.time() - t0:.1f} s)")
    return nbr


# ---------------------------------------------------------------- float scan (candidate choice)
def float_squares(squares, DM):
    return np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2])),
                      math.atan2(s[4], s[3]), s[8] / DM] for s in squares])


def _contains_f(F, px, py):
    ct = np.cos(F[:, 2]); st = np.sin(F[:, 2])
    dx = px - F[:, 0]; dy = py - F[:, 1]
    u = np.abs(dx * ct + dy * st)
    v = np.abs(-dx * st + dy * ct)
    return (u <= 0.5) & (v <= 0.5)


def _meets_seg_f(F, P, Q):
    """float three-axis SAT, unit squares vs the segment PQ"""
    ct = np.cos(F[:, 2]); st = np.sin(F[:, 2])
    ok = np.ones(len(F), bool)
    for (ux, uy) in ((ct, st), (-st, ct)):
        c = F[:, 0] * ux + F[:, 1] * uy
        e = 0.5 * (np.abs(ux * ct + uy * st) + np.abs(-ux * st + uy * ct))
        p0 = P[0] * ux + P[1] * uy
        p1 = Q[0] * ux + Q[1] * uy
        lo = np.minimum(p0, p1); hi = np.maximum(p0, p1)
        ok &= (c + e >= lo) & (hi >= c - e)
    nx, ny = -(Q[1] - P[1]), Q[0] - P[0]
    nn = math.hypot(nx, ny)
    if nn > 0:
        nx /= nn; ny /= nn
        c = F[:, 0] * nx + F[:, 1] * ny
        e = 0.5 * (np.abs(nx * ct + ny * st) + np.abs(-nx * st + ny * ct))
        d = P[0] * nx + P[1] * ny
        ok &= (c + e >= d) & (d >= c - e)
    return ok


def scan_candidates(squares, F, t, Vf, cov, npts, ndir, neps, nrho, epsmax, covfrac=0.5,
                    verbose=True):
    """floats CHOOSE candidate (p, A): p over the heaviest arrangement vertices, A the segment
    perpendicular to a direction d at offset eps and half-length rho (`search/clique_scan.py`'s
    parameterisation, offset cap lifted).  Every hit is re-evaluated exactly by the caller."""
    # anchor points: arrangement vertices of high coverage, deduplicated on a 0.01 grid and then
    # spread evenly over the sorted list (`clique_scan.py`'s `tight_points` rule)
    keep = np.nonzero(cov >= covfrac * cov.max())[0]
    keep = keep[np.argsort(-cov[keep])]
    seen = set()
    ded = []
    for i in keep:
        k = (int(Vf[i, 0] * 100), int(Vf[i, 1] * 100))
        if k in seen:
            continue
        seen.add(k)
        ded.append(int(i))
    step = max(1, len(ded) // max(npts, 1))
    order = ded[::step][:npts]
    if verbose:
        log(f"    scan: {len(keep)} vertices with cov >= {covfrac} max, {len(ded)} distinct on a 0.01 "
            f"grid, {len(order)} anchor points")
    w = F[:, 3]
    dirs = [(math.cos(a), math.sin(a)) for a in np.linspace(0, 2 * math.pi, ndir, endpoint=False)]
    out = []
    t0 = time.time()
    for ii, vi in enumerate(order):
        px, py = Vf[vi]
        thru = _contains_f(F, px, py)
        mthru = w[thru].sum()
        for (dx, dy) in dirs:
            ex, ey = -dy, dx
            for ie in range(1, neps + 1):
                eps = epsmax * ie / neps
                mx, my = px + eps * dx, py + eps * dy
                for ir in range(1, nrho + 1):
                    rho = 0.49 * ir / nrho
                    P = (mx - rho * ex, my - rho * ey)
                    Q = (mx + rho * ex, my + rho * ey)
                    if min(P[0], Q[0]) < -1e-9 or max(P[0], Q[0]) > float(t) + 1e-9:
                        continue
                    if min(P[1], Q[1]) < -1e-9 or max(P[1], Q[1]) > float(t) + 1e-9:
                        continue
                    meet = _meets_seg_f(F, P, Q)
                    contA = _contains_f(F, P[0], P[1]) & _contains_f(F, Q[0], Q[1])
                    val = w[contA].sum() + w[thru & meet & ~contA].sum()
                    if val > 1.0 + 1e-9 or val > mthru + 1e-9:
                        out.append((val, (px, py), P, Q))
        if verbose and ii % 50 == 0:
            log(f"    scan {ii}/{len(order)} best {max((z[0] for z in out), default=0):.6f} "
                f"({time.time() - t0:.0f} s)")
    out.sort(key=lambda z: -z[0])
    return out


def rationalise_segment(P, Q, D):
    return ((int(round(P[0] * D)), int(round(P[1] * D))),
            (int(round(Q[0] * D)), int(round(Q[1] * D))), D)


# ---------------------------------------------------------------- exact local ascent
def local_ascent(squares, Inc, DM, w, V, seed, budget, nvert=24, maxeval=80, verbose=True):
    """exact hill-climb over the candidate set of LEAF_CEILING 2.3.

    From a seed anchor: (i) re-optimise the segment on its own supporting line by the exact 1-D
    problem; (ii) try every line through the current anchor's endpoint region -- lines through two
    arrangement vertices lying in a square that contains the current anchor -- and keep the best.
    Every value is exact.  This gives a certified lower bound (a witness), never a proof."""
    t0 = time.time()
    (P0, P1, Ds) = seed
    near = list(range(len(squares)))
    Vfa = np.array([[X / D, Y / D] for (X, Y, D) in V])
    best = anchor_value(squares, Inc, DM, P0, P1, Ds, w)
    bestseg = (P0, P1, Ds)
    improved = True
    rounds = 0
    while improved and time.time() - t0 < budget:
        improved = False
        rounds += 1
        (P0, P1, Ds) = bestseg
        # (i) the exact 1-D problem on the current line
        if P0 != P1:
            A = (P0[0], P0[1], Ds)
            B = (Fr(P1[0] - P0[0], Ds), Fr(P1[1] - P0[1], Ds))
            r = best_on_line(squares, Inc, DM, w, A, B, near, incumbent=best[0], maxeval=maxeval)
            if r is not None and r[0] > best[0]:
                (s0, u0) = r[1]
                nD = lcm(lcm(s0.denominator, u0.denominator), Ds)
                def pt(uu):
                    x = Fr(P0[0], Ds) + uu * B[0]
                    y = Fr(P0[1], Ds) + uu * B[1]
                    dd = lcm(x.denominator, y.denominator)
                    return (x, y, dd)
                x0, y0, d0 = pt(s0)
                x1, y1, d1 = pt(u0)
                dd = lcm(d0, d1)
                Q0 = (int(x0 * dd), int(y0 * dd))
                Q1 = (int(x1 * dd), int(y1 * dd))
                v = anchor_value(squares, Inc, DM, Q0, Q1, dd, w)
                if v[0] > best[0]:
                    best = v
                    bestseg = (Q0, Q1, dd)
                    improved = True
                    continue
        # (ii) replace one endpoint by a nearby arrangement vertex (a candidate line through it)
        px0, py0 = Fr(bestseg[0][0], bestseg[2]), Fr(bestseg[0][1], bestseg[2])
        px1, py1 = Fr(bestseg[1][0], bestseg[2]), Fr(bestseg[1][1], bestseg[2])
        cand = []
        for (ex, ey, which) in ((px0, py0, 0), (px1, py1, 1)):
            d2 = (Vfa[:, 0] - float(ex)) ** 2 + (Vfa[:, 1] - float(ey)) ** 2
            for vi in np.argsort(d2)[:nvert]:
                X, Y, D = V[int(vi)]
                vx, vy = Fr(X, D), Fr(Y, D)
                cand.append((vx, vy, px1, py1) if which == 0 else (px0, py0, vx, vy))
        for (xa, ya, xb, yb) in cand:
            dd = lcm(lcm(xa.denominator, ya.denominator), lcm(xb.denominator, yb.denominator))
            Q0 = (int(xa * dd), int(ya * dd))
            Q1 = (int(xb * dd), int(yb * dd))
            if Q0 == Q1:
                continue
            v = anchor_value(squares, Inc, DM, Q0, Q1, dd, w)
            if v[0] > best[0]:
                best = v
                bestseg = (Q0, Q1, dd)
                improved = True
        if time.time() - t0 > budget:
            break
    return best, bestseg, rounds


# ---------------------------------------------------------------- complete enumeration
def max_clique_through(nbr, w, j, time_limit, lb=0):
    """exact max-weight clique CONTAINING j, over the closed-intersection graph.

    Sound prune for the exhaustive anchor search: for an anchor A with j in Cont(A), the whole
    clique K(p, A) = Cont(A) u (Meet(A) n P(p)) lies in {j} u N[j] and contains j, so
    mu(K(p, A)) <= this number.  Returns (weight, complete)."""
    nb = sorted(nbr[j])
    pos = {v: i for i, v in enumerate(nb)}
    sub = [set() for _ in nb]
    for a, v in enumerate(nb):
        for u in nbr[v]:
            if u in pos:
                sub[a].add(pos[u])
    ww = [w[v] for v in nb]
    kw, _, comp = max_clique_int(sub, ww, time_limit, lb=max(0, lb - w[j]))
    return w[j] + kw, comp


def exhaust(squares, Inc, DM, w, V, t, budget, incumbent=None, verbose=True, nbr=None,
            clique_time=30.0, ceiling=None):
    """the COMPLETE candidate enumeration of LEAF_CEILING 2.3.

    For every square S_j (the anchor's primary container: the optimum has some j with A subset
    S_j, or else mu(K) = cov(p) <= M), for every pair of distinct arrangement vertices inside
    S_j (the supporting line of the optimal anchor passes through two such), solve the exact 1-D
    problem on that line.  Returns (value, witness, complete)."""
    t0 = time.time()
    n = len(squares)
    Vf = np.array([[X / D, Y / D] for (X, Y, D) in V])
    Cf = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2]))] for s in squares])
    best = (incumbent if incumbent is not None else Fr(0), None)
    complete = True
    npairs = 0
    nskip = 0
    wi = [int(x[8]) for x in squares]
    if nbr is not None:
        # visit the squares in decreasing neighbourhood mass: the incumbent then rises fast and
        # the clique bound prunes the rest outright
        hot = np.array([wi[k] + sum(wi[u] for u in nbr[k]) for k in range(n)])
    else:
        hot = np.array([s[8] for s in squares])
    order = np.argsort(-hot)
    for jj in order:
        j = int(jj)
        if ceiling is not None and best[0] >= ceiling:
            # the incumbent already attains the exact maximum-weight clique of the whole
            # closed-intersection graph, which upper-bounds every K(p, A): nothing is left to find
            if verbose:
                log(f"    exhaust: incumbent {float(best[0]):.9f} attains the exact clique "
                    f"ceiling {float(ceiling):.9f}; the maximum is proved")
            break
        if time.time() - t0 > budget:
            complete = False
            break
        s = squares[j]
        if nbr is not None:
            # SOUND PRUNE: Cont(A) u (Meet(A) n P(p)) is a clique containing j whenever j in
            # Cont(A), so no anchor with A subset S_j can beat the max-weight clique through j.
            lb = int(best[0] * DM)
            if sum(wi[k] for k in nbr[j]) + wi[j] <= lb:      # cheap: the whole neighbourhood
                nskip += 1
                continue
            kb, kcomp = max_clique_through(nbr, wi, j, clique_time, lb=lb)
            if kcomp and kb <= lb:
                nskip += 1
                continue
        # the arrangement vertices S_j contains -- already computed exactly by `incidences`
        inside = [int(i) for i in Inc.cols[j]]
        near = [k for k in range(n)
                if (Cf[k, 0] - Cf[j, 0]) ** 2 + (Cf[k, 1] - Cf[j, 1]) ** 2 <= 2.0 + 1e-6]
        npairs += len(inside) * (len(inside) - 1) // 2
        for ai in range(len(inside)):
            if time.time() - t0 > budget:
                complete = False
                break
            X0, Y0, D0 = V[inside[ai]]
            for bi in range(ai + 1, len(inside)):
                X1, Y1, D1 = V[inside[bi]]
                bx = Fr(X1, D1) - Fr(X0, D0)
                by = Fr(Y1, D1) - Fr(Y0, D0)
                if bx == 0 and by == 0:
                    continue
                r = best_on_line(squares, Inc, DM, w, (X0, Y0, D0), (bx, by), near,
                                 incumbent=best[0])
                if r is not None and r[0] > best[0]:
                    best = (r[0], (j, V[inside[ai]], V[inside[bi]], r[1]))
        if verbose:
            log(f"    exhaust j={j} ({len(inside)} vertices inside, {npairs} candidate lines so "
                f"far, {nskip} squares pruned by the clique bound) best {float(best[0]):.9f} "
                f"({time.time() - t0:.0f} s)")
    return best[0], best[1], complete


# ============================================================================ the certifier
def certify(args):
    global LOG
    t_all = time.time()
    t, sym, poses = read_measure(args.FILE, Fr(args.t) if args.t else None)
    squares, DM, pose_of = build_squares(t, sym, poses)
    mass = Fr(sum(s[8] for s in squares), DM)
    log(f"leaf_ceiling: {args.FILE}")
    log(f"  t = {t} = {float(t):.6f}; sym = {sym}; {len(poses)} poses -> {len(squares)} squares "
        f"with positive mass; mass denominator {DM}")
    log(f"  MASS = {mass} = {float(mass):.12f}")
    res = dict(file=args.FILE, t=str(t), sym=sym, poses=len(poses), squares=len(squares),
               DM=DM, mass=str(mass), mass_float=float(mass))

    # ---------------------------------------------------------------- (a) coverage
    V = enumerate_vertices(squares, t, procs=args.procs)
    Inc = incidences(squares, V, procs=args.procs)
    w = np.array([s[8] for s in squares], dtype=np.int64)
    M, mi, cov = max_coverage(Inc, squares, DM)
    MCAP[0] = int(cov.max())
    X, Y, D = V[mi]
    log(f"  (a) COVERAGE  M = max cov = {M} = {float(M):.15f}  "
        f"{'<= 1  OK' if M <= 1 else '> 1  FAIL'}   (attained at ({X}/{D}, {Y}/{D}) "
        f"~ ({X / D:.6f}, {Y / D:.6f}), {len(V)} vertices)")
    res.update(vertices=len(V), M=str(M), M_float=float(M), coverage_ok=bool(M <= 1))

    # ---------------------------------------------------------------- (c) regions
    meta = [(Fr(s[0], s[2]), Fr(s[1], s[2]), s[8]) for s in squares]
    boxes = region_boxes(t, Fr(args.r))
    target = region_targets(t, args)
    ok, assigned, amb, note = region_report(meta, boxes, target, DM)
    labs = [b[0] for b in boxes]
    log(f"  (c) REGIONS   " + "  ".join(
        f"{lab}={float(Fr(assigned.get(lab, 0), DM)):.6f}" for lab in labs))
    log(f"                {len(amb)} poses with the centre on a region boundary "
        f"(mass {float(Fr(sum(a[0] for a in amb), DM)):.6f}); {note}")
    if target:
        log(f"      targets {', '.join(f'{k}={v}' for k, v in sorted(target.items()))}: "
            f"{'OK' if ok else 'FAIL'}")
    res.update(regions={lab: str(Fr(assigned.get(lab, 0), DM)) for lab in labs},
               regions_ok=bool(ok), ambiguous=len(amb),
               region_target={k: str(v) for k, v in (target or {}).items()})
    chord_ok = None
    if args.chord:
        # mu(wall strip) for the four strips [0,t]x[0,r] and images
        strips = [('S0', Fr(0), t, Fr(0), Fr(args.r)), ('S1', Fr(0), t, t - Fr(args.r), t),
                  ('S2', Fr(0), Fr(args.r), Fr(0), t), ('S3', t - Fr(args.r), t, Fr(0), t)]
        chord_ok = True
        vals = []
        for (lab, x0, x1, y0, y1) in strips:
            v = Fr(sum(m for (cx, cy, m) in meta if x0 <= cx <= x1 and y0 <= cy <= y1), DM)
            vals.append((lab, v))
            if v > 3:
                chord_ok = False
        log("      chord: " + "  ".join(f"{lab}={float(v):.6f}" for lab, v in vals)
            + f"   {'<= 3 OK' if chord_ok else '> 3 FAIL'}")
        res.update(chord={lab: str(v) for lab, v in vals}, chord_ok=chord_ok)

    # ---------------------------------------------------------------- (b) anchor cliques
    anchor = dict()
    nbr = None
    mode = args.anchor
    if mode in ('clique', 'all'):
        nbr = closed_graph(squares)
        Kw, Kc, comp = max_clique_int(nbr, [int(s[8]) for s in squares], args.clique_time)
        K = Fr(Kw, DM)
        log(f"  (b) CLIQUE UPPER BOUND  max-weight closed-intersection clique = {K} = "
            f"{float(K):.12f} (size {len(Kc)}, B&B complete: {comp})"
            + ("  => every K(p,A) is <= 1: (b) PROVED" if comp and K <= 1 else ""))
        anchor.update(clique_bound=str(K), clique_bound_float=float(K),
                      clique_complete=bool(comp), clique_size=len(Kc))
    if mode in ('scan', 'all'):
        F = float_squares(squares, DM)
        Vf = np.array([[X / D, Y / D] for (X, Y, D) in V])
        cand = scan_candidates(squares, F, t, Vf, cov.astype(float) / DM, args.scan_pts,
                               args.scan_dir, args.scan_eps, args.scan_rho, args.epsmax,
                               covfrac=args.scan_covfrac)
        log(f"      scan: {len(cand)} float candidates over 1; best float {cand[0][0]:.6f}"
            if cand else "      scan: no float candidate over 1")
        best = (Fr(0), None, None)
        seen = 0
        for (val, p, P, Q) in cand[:args.scan_verify]:
            P0, P1, Ds = rationalise_segment(P, Q, args.D)
            if P0 == P1:
                continue
            v = anchor_value(squares, Inc, DM, P0, P1, Ds, w)
            seen += 1
            if v[0] > best[0]:
                best = (v[0], (P0, P1, Ds), v)
        if best[1] is not None:
            log(f"      exact on the {seen} best float candidates: {float(best[0]):.9f}")
            bv, bseg, rounds = local_ascent(squares, Inc, DM, w, V, best[1], args.ascent_time)
            log(f"      exact local ascent ({rounds} rounds): "
                f"{bv[0]} = {float(bv[0]):.9f}  (|Cont| = {bv[2]}, |Meet| = {bv[3]})")
            P0, P1, Ds = bseg
            vX, vY, vD = V[bv[1]] if bv[1] >= 0 else (0, 0, 1)
            log(f"      witness: p = ({vX}/{vD}, {vY}/{vD}) ~ ({vX / vD:.6f}, {vY / vD:.6f}); "
                f"A = ({P0[0]}/{Ds}, {P0[1]}/{Ds})-({P1[0]}/{Ds}, {P1[1]}/{Ds}) "
                f"~ ({P0[0] / Ds:.5f}, {P0[1] / Ds:.5f})-({P1[0] / Ds:.5f}, {P1[1] / Ds:.5f})")
            # INDEPENDENT re-check of the witness: rebuild K(p, A) from the definition and test
            # every unordered pair of members with the exact four-axis square-square SAT.  Lemma 0
            # says the answer must be 0 non-intersecting pairs; this is the check that catches a
            # membership bug (the missing-axis near miss of `notes/clique-family.md` 7).
            cont, meet = anchor_members(squares, P0, P1, Ds)
            mem = [i for i in range(len(squares))
                   if cont[i] or (meet[i] and sq_contains(squares[i], vX, vY, vD))]
            bad = 0
            for ii in range(len(mem)):
                for jj in range(ii + 1, len(mem)):
                    if not sq_meets_sq(squares[mem[ii]], squares[mem[jj]]):
                        bad += 1
            recomputed = Fr(sum(squares[i][8] for i in mem), DM)
            log(f"      witness re-check: {len(mem)} members, mass {recomputed} = "
                f"{float(recomputed):.9f}, {bad} non-intersecting pairs "
                f"({'OK' if bad == 0 and recomputed == bv[0] else 'MISMATCH'})")
            anchor.update(scan_max=str(bv[0]), scan_max_float=float(bv[0]),
                          witness_p=[vX, vY, vD], witness_A=[list(P0), list(P1), Ds],
                          n_contains=bv[2], n_meets=bv[3], witness_members=len(mem),
                          witness_badpairs=bad, witness_recomputed=str(recomputed))
        else:
            log("      scan found nothing over 1 (a lower bound only)")
            anchor.update(scan_max='0', scan_max_float=0.0)
    if mode in ('exhaust', 'all'):
        inc0 = None
        if 'scan_max' in anchor:
            inc0 = Fr(anchor['scan_max'])
        if nbr is None:
            nbr = closed_graph(squares)
        ceil0 = None
        if anchor.get('clique_complete'):
            ceil0 = Fr(anchor['clique_bound'])
        val, wit, comp = exhaust(squares, Inc, DM, w, V, t, args.budget, incumbent=inc0,
                                 nbr=nbr, clique_time=args.clique_time / 10.0, ceiling=ceil0)
        log(f"  (b) EXHAUST  max over the complete candidate set = {val} = {float(val):.12f} "
            f"(complete: {comp})")
        anchor.update(exhaust_max=str(val), exhaust_max_float=float(val), exhaust_complete=comp)
    res['anchor'] = anchor

    # ---------------------------------------------------------------- verdict
    proved_b = (anchor.get('clique_complete') and Fr(anchor.get('clique_bound', '99')) <= 1) or \
               (anchor.get('exhaust_complete') and Fr(anchor.get('exhaust_max', '99')) <= 1)
    lower_b = max([Fr(anchor.get(k, '0')) for k in ('scan_max', 'exhaust_max')] + [Fr(0)])
    log("")
    log(f"  VERDICT  mass = {float(mass):.9f}   coverage {'OK' if M <= 1 else 'FAIL'}   "
        f"regions {'OK' if ok else ('FAIL' if target else 'n/a')}   "
        f"anchor cliques {'PROVED <= 1' if proved_b else f'NOT PROVED (best found {float(lower_b):.9f})'}")
    if M <= 1 and ok and proved_b and mass >= 12:
        log(f"  ==> CERTIFIED NO-GO: a measure of mass {float(mass):.9f} >= 12 feasible for points, "
            f"anchor cliques and the region equalities at t = {t}.")
    res['proved_b'] = bool(proved_b)
    res['lower_b'] = str(lower_b)
    res['seconds'] = time.time() - t_all
    log(f"  ({time.time() - t_all:.1f} s)")
    if args.json:
        json.dump(res, open(args.json, 'w'), indent=1)
    return res


# ============================================================================ snap (floats -> exact)
def rational_angle(th_rad, Q):
    p = round(math.tan(th_rad / 2) * Q)
    q = Q
    g = gcd(abs(p), q)
    return p // g, q // g


def snap_pose(cx, cy, th_rad, t, Q, Dc):
    p, q = rational_angle(th_rad, Q)
    a, b, r = cos_sin(p, q)
    w2 = half_width(a, b, r)
    lo, hi = w2, t - w2
    assert lo <= hi
    x = min(max(Fr(round(cx * Dc), Dc), lo), hi)
    y = min(max(Fr(round(cy * Dc), Dc), lo), hi)
    return p, q, x, y


def read_float_poses(path):
    """`pose cx cy theta_deg mu` (t4screen / packing_dual support files)"""
    out = []
    t = None
    for line in open(path):
        s = line.strip()
        if s.startswith('#'):
            for tok in s.replace(',', ' ').split():
                if tok.startswith('t='):
                    try:
                        t = float(tok[2:])
                    except ValueError:
                        pass
            continue
        q = s.split()
        if q and q[0] == 'pose' and len(q) == 5:
            out.append((float(q[1]), float(q[2]), math.radians(float(q[3])), float(q[4])))
    return t, out


def cmd_snap(args):
    t = Fr(args.t)
    tf, src = read_float_poses(args.FILE)
    src = [z for z in src if z[3] > 1e-12]
    log(f"snap: {len(src)} float poses with positive mass from {args.FILE}; t = {t}; "
        f"Q = {args.Q}, Dc = {args.Dc}, DM = {args.DM}")
    poses = []
    dmax = 0.0
    for cx, cy, th, mu in src:
        p, q, x, y = snap_pose(cx, cy, th, t, args.Q, args.Dc)
        a, b, r = cos_sin(p, q)
        dmax = max(dmax, abs(math.atan2(b, a) - th), abs(float(x) - cx), abs(float(y) - cy))
        poses.append((p, q, x, y, mu))
    log(f"  max snap displacement (angle rad / centre): {dmax:.2e}")
    squares, _, _ = build_squares(t, 1, [(p, q, x, y, Fr(1)) for (p, q, x, y, _) in poses], DM=1)
    log(f"  {len(squares)} squares, all exactly admissible")
    V = enumerate_vertices(squares, t, procs=args.procs)
    Inc = incidences(squares, V, procs=args.procs)
    DM = args.DM
    MU = [int(math.floor(mu * DM)) for (_, _, _, _, mu) in poses]
    w = np.array(MU, dtype=np.int64)
    cov = Inc.cov(w)
    M = Fr(int(cov.max()), DM)
    log(f"  rounded down to /{DM}: mass {sum(MU) / DM:.9f}, exact M = {float(M):.12f}")
    if args.polish:
        # re-optimise the masses on the EXACT rows of this pose set.  The rows are the whole
        # arrangement of the snapped support, so they are complete for every sub-measure of it
        # (Lemma 1): the LP value is a value the poses really attain, not an over-estimate on a
        # sampled row set.  Floats only CHOOSE the masses; everything is re-checked exactly below.
        from scipy.optimize import linprog
        Aeq = beq = None
        tgt0 = region_targets(t, args)
        if tgt0:
            boxes0 = region_boxes(t, Fr(args.r))
            erows = []
            vals = []
            for lab, kv in sorted(tgt0.items()):
                erows.append(np.array([1.0 if lab in regions_of(x, y, boxes0) else 0.0
                                       for (p, q, x, y, _) in poses]))
                vals.append(float(kv))
            Aeq = np.array(erows)
            beq = np.array(vals)
            log(f"  polish: {len(vals)} region equalities "
                + ", ".join(f"{k}={v}" for k, v in sorted(tgt0.items())))
        # ROW GENERATION over the EXACT vertex set.  The vertex set is complete for every
        # sub-measure of this pose set (Lemma 1), so when no vertex is violated the LP value is a
        # value the poses really attain -- and the whole vertex set (10^5-10^6 rows, 10^8
        # incidences) is far too big to hand to the solver at once.
        rows = [np.argsort(-Inc.cov(np.array(MU, dtype=np.int64)))[:args.rows0]]
        for c in Inc.cols:                 # every column needs a row, or the LP is unbounded
            if len(c):
                rows.append(c[np.linspace(0, len(c) - 1, min(8, len(c))).astype(int)])
        rows = np.unique(np.concatenate(rows))
        for it in range(args.polish_rounds):
            A = Inc.to_csr(rows)
            res = linprog(-np.ones(len(poses)), A_ub=A,
                          b_ub=np.full(A.shape[0], 1.0 - args.margin),
                          A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
            if res.status != 0:
                log(f"  polish {it}: LP failed ({res.message}); keeping the previous masses")
                break
            MUf = np.maximum(res.x, 0.0)
            MU = [int(math.floor(x * DM)) for x in MUf]
            cov = Inc.cov(np.array(MU, dtype=np.int64))
            bad = np.nonzero(cov > DM)[0]
            log(f"  polish {it}: {A.shape[0]} rows, LP mass {-res.fun:.9f}, exact mass after "
                f"rounding {sum(MU) / DM:.9f}, exact M = {float(Fr(int(cov.max()), DM)):.12f}, "
                f"{len(bad)} violated vertices")
            if not len(bad):
                break
            add = bad[np.argsort(-cov[bad])[:args.row_cap]]
            rows = np.unique(np.concatenate([rows, add]))
        w = np.array(MU, dtype=np.int64)
        cov = Inc.cov(w)
        M = Fr(int(cov.max()), DM)
        log(f"  polished: mass {sum(MU) / DM:.9f}, exact M = {float(M):.12f}")
    if M > 1:
        MU = [(m * M.denominator) // M.numerator for m in MU]
        w = np.array(MU, dtype=np.int64)
        cov = Inc.cov(w)
        M = Fr(int(cov.max()), DM)
        log(f"  scaled by 1/M and rounded down: mass {sum(MU) / DM:.9f}, exact M = {float(M):.12f}")
    assert M <= 1

    # ---- region top-up: rounding DOWN leaves each pinned region a hair short of its leaf count.
    # Adding mass raises coverage, so the deficit is spent greedily on the poses with EXACT
    # coverage slack: for pose i, slack_i = min over the arrangement vertices of S_i of
    # (1 - cov(v)), and adding d <= slack_i to pose i keeps coverage <= 1 exactly.
    target = region_targets(t, args)
    if target:
        boxes = region_boxes(t, Fr(args.r))
        cap = int(DM)                       # 1 unit of mass, in integer units
        for lab, kv in sorted(target.items()):
            want = int(kv * DM)
            assert Fr(want, DM) == kv
            idx = [i for i, (p, q, x, y, _) in enumerate(poses)
                   if lab in regions_of(x, y, boxes)]
            have = sum(MU[i] for i in idx)
            d = want - have
            if d == 0:
                log(f"  region {lab}: exactly {kv} already")
                continue
            if d < 0:
                # too much: take it off the heaviest pose that is only in this region
                only = [i for i in idx if regions_of(poses[i][2], poses[i][3], boxes) == [lab]]
                for i in sorted(only, key=lambda i: -MU[i]):
                    take = min(MU[i], -d)
                    MU[i] -= take
                    d += take
                    if d == 0:
                        break
                log(f"  region {lab}: removed {(want - have - d) / DM:+.3e} to reach {kv}"
                    if d == 0 else f"  region {lab}: STILL {d / DM:+.3e} off {kv}")
                w = np.array(MU, dtype=np.int64)
                cov = Inc.cov(w)
                continue
            for _ in range(200):
                if d == 0:
                    break
                w = np.array(MU, dtype=np.int64)
                cov = Inc.cov(w)
                slack = {i: int(cap - cov[Inc.cols[i]].max()) if len(Inc.cols[i]) else cap
                         for i in idx}
                i = max(idx, key=lambda i: slack[i])
                if slack[i] <= 0:
                    break
                add = min(d, slack[i])
                MU[i] += add
                d -= add
            w = np.array(MU, dtype=np.int64)
            cov = Inc.cov(w)
            M = Fr(int(cov.max()), DM)
            log(f"  region {lab}: topped up to {Fr(sum(MU[i] for i in idx), DM)} "
                f"(target {kv}, residual {d / DM:+.3e}); M now {float(M):.12f}")
        w = np.array(MU, dtype=np.int64)
        cov = Inc.cov(w)
        M = Fr(int(cov.max()), DM)
        log(f"  after region top-up: mass {sum(MU) / DM:.9f}, exact M = {float(M):.12f}")
        assert M <= 1

    with open(args.out, 'w') as f:
        f.write(f"# leaf_ceiling measure; t = {t}; sym 1; snapped from "
                f"{os.path.basename(args.FILE)} (Q={args.Q}, Dc={args.Dc}, DM={DM})\n")
        f.write(f"# mass = {Fr(sum(MU), DM)} = {sum(MU) / DM:.12f}; M = {M}\n")
        f.write("# pose p q cx cy mass : theta = 2 arctan(p/q)\n")
        for (p, q, x, y, _), m in zip(poses, MU):
            if m > 0:
                f.write(f"pose {p} {q} {x} {y} {Fr(m, DM)}\n")
    log(f"  wrote {args.out}: {sum(1 for m in MU if m > 0)} poses, mass {sum(MU) / DM:.9f}")


# ============================================================================ selftest
def _sq(cx, cy, deg, mass, t):
    """a test pose: rational rotation nearest to `deg` with denominator 10^5"""
    p, q = rational_angle(math.radians(deg), 10 ** 5)
    return (p, q, Fr(cx), Fr(cy), Fr(mass))


def selftest():
    ok = True

    def check(name, cond, extra=''):
        nonlocal ok
        log(f"  {'PASS' if cond else 'FAIL'}  {name}  {extra}")
        ok = ok and cond

    t = Fr(4)
    # --- two disjoint axis-aligned squares, mass 1 each: coverage 1, no violated clique
    poses = [(0, 1, Fr(1, 2), Fr(1, 2), Fr(1)), (0, 1, Fr(3), Fr(3), Fr(1))]
    sq, DM, _ = build_squares(t, 1, poses)
    V = enumerate_vertices(sq, t, verbose=False)
    Inc = incidences(sq, V, verbose=False)
    M, _, _ = max_coverage(Inc, sq, DM)
    check('two disjoint squares: mass 2', Fr(sum(s[8] for s in sq), DM) == 2)
    check('two disjoint squares: M = 1', M == 1, f'M = {M}')
    nbr = closed_graph(sq, verbose=False)
    Kw, _, comp = max_clique_int(nbr, [int(s[8]) for s in sq], 10)
    check('two disjoint squares: max clique = 1', Fr(Kw, DM) == 1 and comp, f'K = {Fr(Kw, DM)}')

    # --- two touching squares (share the edge x = 1): still a packing, but a clique of mass 2
    poses = [(0, 1, Fr(1, 2), Fr(1, 2), Fr(1)), (0, 1, Fr(3, 2), Fr(1, 2), Fr(1))]
    sq, DM, _ = build_squares(t, 1, poses)
    V = enumerate_vertices(sq, t, verbose=False)
    Inc = incidences(sq, V, verbose=False)
    M, _, _ = max_coverage(Inc, sq, DM)
    check('two touching squares: M = 2 (the shared edge is covered twice)', M == 2, f'M = {M}')

    # --- two overlapping squares of mass 1: coverage 2 > 1, must NOT certify
    poses = [(0, 1, Fr(1, 2), Fr(1, 2), Fr(1)), (0, 1, Fr(1), Fr(1), Fr(1))]
    sq, DM, _ = build_squares(t, 1, poses)
    V = enumerate_vertices(sq, t, verbose=False)
    Inc = incidences(sq, V, verbose=False)
    M, _, _ = max_coverage(Inc, sq, DM)
    check('two overlapping squares: M = 2 > 1 (does not certify)', M == 2, f'M = {M}')

    # --- four disjoint unit squares of mass 1: everything certifies, mass 4
    poses = [(0, 1, Fr(1, 2), Fr(1, 2), Fr(1)), (0, 1, Fr(7, 2), Fr(1, 2), Fr(1)),
             (0, 1, Fr(1, 2), Fr(7, 2), Fr(1)), (0, 1, Fr(7, 2), Fr(7, 2), Fr(1))]
    sq, DM, _ = build_squares(t, 1, poses)
    V = enumerate_vertices(sq, t, verbose=False)
    Inc = incidences(sq, V, verbose=False)
    M, _, _ = max_coverage(Inc, sq, DM)
    w = np.array([s[8] for s in sq], dtype=np.int64)
    val, _, comp2 = exhaust(sq, Inc, DM, w, V, t, 120, verbose=False)
    check('four disjoint corner squares: M = 1, exhaust max = 1, complete',
          M == 1 and val == 1 and comp2, f'M = {M}, anchor max = {val}, complete = {comp2}')

    # --- the PINWHEEL: three unit squares at 0, 30, 60 degrees, centres at distance 3/5 from
    #     (2,2) in the directions 90, 210, 330 degrees.  They pairwise closed-intersect and have
    #     EMPTY triple intersection (3-fold symmetric and convex, so a common point would force
    #     the centroid in), so with mass 1/2 each the coverage is exactly 1 -- and yet the anchor
    #     clique K(p, A) with p in S1 n S2 and A a segment of S3 meeting both has mass 3/2.
    #     This is the discriminating case: coverage feasible, anchor-clique INfeasible.
    R = Fr(3, 5)
    cs = []
    for k in range(3):
        phi = math.radians(90 + 120 * k)
        cs.append((Fr(2) + Fr(round(float(R) * math.cos(phi) * 10 ** 6), 10 ** 6),
                   Fr(2) + Fr(round(float(R) * math.sin(phi) * 10 ** 6), 10 ** 6),
                   math.radians(30 * k)))
    poses = []
    for (cx, cy, th) in cs:
        p, q = rational_angle(th, 10 ** 5)
        poses.append((p, q, cx, cy, Fr(1, 2)))
    sq, DM, _ = build_squares(t, 1, poses)
    V = enumerate_vertices(sq, t, verbose=False)
    Inc = incidences(sq, V, verbose=False)
    M, _, _ = max_coverage(Inc, sq, DM)
    w = np.array([s[8] for s in sq], dtype=np.int64)
    check('pinwheel: M = 1 (coverage feasible)', M == 1, f'M = {M}')
    nbr = closed_graph(sq, verbose=False)
    Kw, _, comp = max_clique_int(nbr, [int(s[8]) for s in sq], 10)
    check('pinwheel: the three squares pairwise meet, max clique = 3/2',
          Fr(Kw, DM) == Fr(3, 2) and comp, f'K = {Fr(Kw, DM)}')
    val, wit, comp2 = exhaust(sq, Inc, DM, w, V, t, 300, verbose=False)
    check('pinwheel: the complete anchor enumeration finds 3/2 and finishes',
          val == Fr(3, 2) and comp2, f'anchor max = {val}, complete = {comp2}')

    # --- three squares in a row: overlapping, coverage 3/2 at the double crossing
    poses = [(0, 1, Fr(1, 2), Fr(1, 2), Fr(1, 2)), (0, 1, Fr(1), Fr(1, 2), Fr(1, 2)),
             (0, 1, Fr(3, 2), Fr(1, 2), Fr(1, 2))]
    sq, DM, _ = build_squares(t, 1, poses)
    V = enumerate_vertices(sq, t, verbose=False)
    Inc = incidences(sq, V, verbose=False)
    M, _, _ = max_coverage(Inc, sq, DM)
    w = np.array([s[8] for s in sq], dtype=np.int64)
    check('three-in-a-row: M = 3/2 > 1 (does not certify)', M == Fr(3, 2), f'M = {M}')
    P0 = (3 * 10 ** 6 // 4, 10 ** 6 // 2)
    P1 = (5 * 10 ** 6 // 4, 10 ** 6 // 2)
    v = anchor_value(sq, Inc, DM, P0, P1, 10 ** 6, w)
    check('three-in-a-row: the anchor [0.75,0.5]-[1.25,0.5] gives mass 3/2',
          v[0] == Fr(3, 2), f'{v[0]}')
    v2 = anchor_value_at(sq, DM, P0, P1, 10 ** 6, 10 ** 6, 10 ** 6 // 2, 10 ** 6)
    check('three-in-a-row: the direct definition agrees at p = (1, 1/2)',
          v2 == Fr(3, 2), f'{v2}')

    # --- exact/float cross-check of the segment SAT on random instances
    rng = np.random.default_rng(7)
    bad = 0
    for _ in range(400):
        cx, cy = rng.uniform(0.8, 3.2, 2)
        deg = rng.uniform(0, 90)
        p, q = rational_angle(math.radians(deg), 10 ** 5)
        s = make_square(Fr(round(cx * 10 ** 6), 10 ** 6), Fr(round(cy * 10 ** 6), 10 ** 6), p, q,
                        t, 1)
        P = rng.uniform(0.5, 3.5, 2)
        Q = P + rng.uniform(-0.9, 0.9, 2)
        D = 10 ** 6
        P0 = (int(round(P[0] * D)), int(round(P[1] * D)))
        P1 = (int(round(Q[0] * D)), int(round(Q[1] * D)))
        e = sq_meets_segment(s, P0, P1, D)
        F = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2])), math.atan2(s[4], s[3]), 1.0]])
        f = bool(_meets_seg_f(F, (P0[0] / D, P0[1] / D), (P1[0] / D, P1[1] / D))[0])
        if e != f:
            bad += 1
    check('segment SAT: exact vs float agree on 400 random instances', bad == 0, f'{bad} mismatches')

    # --- a mid-wall pose is counted in either slot
    boxes = region_boxes(Fr(4), Fr(1))
    check('mid-wall pose (2, 1/2) is in both slots W0 and W1',
          set(regions_of(Fr(2), Fr(1, 2), boxes)) == {'W0', 'W1'},
          str(regions_of(Fr(2), Fr(1, 2), boxes)))
    check('grid pose (3/2, 1/2) is in slot W0 only',
          regions_of(Fr(3, 2), Fr(1, 2), boxes) == ['W0'])
    check('corner pose (1/2, 1/2) is in corner box C0 only',
          regions_of(Fr(1, 2), Fr(1, 2), boxes) == ['C0'])
    check('pose (1, 1) is in C0, W0, W7 and I',
          set(regions_of(Fr(1), Fr(1), boxes)) == {'C0', 'W0', 'W7', 'I'},
          str(regions_of(Fr(1), Fr(1), boxes)))
    # region feasibility with a mid-wall pose: mass 1 at (2, 1/2) can satisfy W0 = 1, W1 = 0
    meta = [(Fr(2), Fr(1, 2), 1)]
    okA, asg, amb, _ = region_report(meta, boxes, {'W0': Fr(1), 'W1': Fr(0)}, 1)
    check('mid-wall mass 1 can be declared in W0', okA and asg['W0'] == 1)
    okB, asg, amb, _ = region_report(meta, boxes, {'W0': Fr(0), 'W1': Fr(1)}, 1)
    check('mid-wall mass 1 can be declared in W1', okB and asg['W1'] == 1)
    okC, *_ = region_report(meta, boxes, {'W0': Fr(1), 'W1': Fr(1)}, 1)
    check('mid-wall mass 1 cannot satisfy W0 = W1 = 1', not okC)

    log(f"  selftest: {'ALL PASS' if ok else 'FAILURES'}")
    return 0 if ok else 1


# ============================================================================ main
def main():
    global LOG
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    c = sub.add_parser('check')
    c.add_argument('FILE')
    c.add_argument('--t', default=None)
    c.add_argument('--r', default='1')
    c.add_argument('--corners', default=None, help="e.g. 1111 ('.' = free)")
    c.add_argument('--patterns', default=None, help="e.g. 01010101 ('.' = free)")
    c.add_argument('--chord', action='store_true', help='also check mu(wall strip) <= 3')
    c.add_argument('--anchor', default='all', choices=['none', 'scan', 'exhaust', 'clique', 'all'])
    c.add_argument('--budget', type=float, default=1800.0, help='seconds for the exhaustive search')
    c.add_argument('--clique-time', type=float, default=600.0)
    c.add_argument('--ascent-time', type=float, default=600.0)
    c.add_argument('--scan-pts', type=int, default=200)
    c.add_argument('--scan-dir', type=int, default=24)
    c.add_argument('--scan-eps', type=int, default=12)
    c.add_argument('--scan-rho', type=int, default=10)
    c.add_argument('--scan-verify', type=int, default=60)
    c.add_argument('--epsmax', type=float, default=0.9)
    c.add_argument('--scan-covfrac', type=float, default=0.5,
                   help='scan anchor points: arrangement vertices with cov >= this fraction of M')
    c.add_argument('--D', type=int, default=10 ** 6, help='denominator of the scanned anchors')
    c.add_argument('--procs', type=int, default=4)
    c.add_argument('--json', default=None)
    c.add_argument('--log', default=None)

    s = sub.add_parser('snap')
    s.add_argument('FILE')
    s.add_argument('--t', required=True)
    s.add_argument('--out', required=True)
    s.add_argument('--Q', type=int, default=10 ** 5)
    s.add_argument('--Dc', type=int, default=10 ** 6)
    s.add_argument('--DM', type=int, default=10 ** 9)
    s.add_argument('--polish', action='store_true')
    s.add_argument('--margin', type=float, default=0.0)
    s.add_argument('--procs', type=int, default=4)
    s.add_argument('--r', default='1')
    s.add_argument('--corners', default=None, help="region top-up target, e.g. 1111")
    s.add_argument('--rows0', type=int, default=20000, help='initial exact rows for --polish')
    s.add_argument('--row-cap', type=int, default=20000, help='rows added per polish round')
    s.add_argument('--polish-rounds', type=int, default=30)
    s.add_argument('--patterns', default=None, help="region top-up target, e.g. 01010101")
    s.add_argument('--log', default=None)

    sub.add_parser('selftest')

    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    if getattr(a, 'log', None):
        LOG = open(a.log, 'w')
    if a.cmd == 'check':
        certify(a)
    elif a.cmd == 'snap':
        cmd_snap(a)
    else:
        sys.exit(selftest())


if __name__ == '__main__':
    main()
