#!/usr/bin/env python3
"""EXACT (rational, integer-only) certification of the packing lower bound  L(t) <= nu_f(t).

The container side is  t = --t  (an exact fraction; default 399/100, the certificate of
search/DUAL_EXACT.md -- `--t 4` is the certificate of search/COVER4.md).

Everything load-bearing is Python integers / Fractions:

  * every pose is a closed unit square with a RATIONAL rotation  theta = 2*arctan(p/q),
        cos = (q^2-p^2)/(q^2+p^2),  sin = 2pq/(q^2+p^2),
    and a rational centre; the D4 symmetrisation (8 dihedral images in [0,t]^2, t = 399/100)
    is exact, and every image is checked to lie in the closed container (4 rational corners);
  * the arrangement vertices (square corners, and intersection points of every pair of
    non-parallel closed edges) are enumerated exactly as integer triples (X, Y, D) = (X/D, Y/D);
  * the coverage  cov(v) = sum_k mu_k/8 * #{images of pose k containing v}  is evaluated at every
    vertex with the exact containment test  2|a dx + b dy| <= r D Dk,  2|a dy - b dx| <= r D Dk;
  * masses are rationals (denominator DM);  M = max cov over all vertices is an exact Fraction
    and  L = sum(mu) / M  is an exact Fraction (mu/M is a feasible packing measure).

Why the vertex set suffices (closed squares): let p maximise cov and let S be the set of support
squares containing p.  P = intersection of S is a nonempty closed convex polygon (all of whose
points have coverage >= cov(p), since they lie in every square of S).  Every vertex of P is either
a corner of a square of S or the intersection point of two non-parallel edges (closed segments) of
two squares of S -- all such points are enumerated.  So max cov = max over the enumerated set.
Container corners are not even needed (P lies inside the squares).  Since the support is exactly
D4-invariant, cov(g p) = cov(p) for g in D4 and the arrangement is D4-invariant, so it suffices to
evaluate the vertices in the fundamental domain F = {0 <= x <= y <= t/2}  (--full checks all).

Floating point is used ONLY to choose the masses (a HiGHS LP over the exact incidence rows); they
are then rounded DOWN to rationals and everything is re-checked exactly.

Usage
    python3 search/dual_exact.py build [--t 399/100] [--tag TAG] [--src FILE]
                                       [--Q 100000] [--Dc 1000000] [--procs 4]
        reads a float support (default runs/dual_PA2_support.txt), snaps, enumerates, polishes,
        certifies, writes runs/dual_exact_<TAG>_support.txt (+ .json, .log); TAG defaults to t
    python3 search/dual_exact.py check search/dual_exact_3.99_support.txt [--t T] [--full] [--procs 4]
                                       [--n N] [--stream] [--start-method fork|spawn|forkserver]
        independent exact re-certification from the exact support file (no LP); t is read from the
        file's "# t = p/q" header unless --t is given (a header that disagrees with --t is an error).
        With --n N the exit status ASSERTS the strict inequality  L > N  (exact Fractions): exit 0 iff
        it holds, exit 1 if it fails.  A malformed support (a non-comment line that is not
        `pose p q cx cy mass` with integers p, q > 0 and fractions cx, cy, mass >= 0) exits 2.
        Without --n the log line still compares L with 12 for display, and the exit status only
        reports that the check ran.

Multiprocessing: workers receive the arrangement state through an explicit Pool initializer, so the
checker is correct under every start method (fork, spawn, forkserver; --start-method, default the
platform's); --procs 1 runs everything in-process with no pool.

The angle snap Q is the binding approximation: the LP optimum is degenerate on tight vertices, so
Q = 10^5 can cost 0.3 of mass on a support the float search polished.  Q = 10^7 reproduces the
float value to 8 decimals at t = 4 (COVER4.md) at ~3x the integer cost.
"""
import sys, os, re, math, time, json, argparse
from fractions import Fraction as Fr
from math import gcd, lcm
from collections import Counter
import multiprocessing as mp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
TN, TD = 399, 100                     # container side t = TN/TD (default 3.99; --t sets it)
T = Fr(TN, TD)
TAG = '3.99'                          # file tag for runs/dual_exact_<TAG>_*

def set_t(frac, tag=None):
    """set the container side t = TN/TD exactly (a Fraction) and the output file tag"""
    global TN, TD, T, TAG
    T = Fr(frac); TN, TD = T.numerator, T.denominator
    TAG = tag if tag is not None else (str(TN) if TD == 1 else f"{float(T):g}")
G_EDGE = 5                            # grid resolution for the edge-pair prefilter (exact floors)
G_SQ = 10                             # grid resolution for the vertex -> square prefilter

LOG = None
def log(msg):
    print(msg, flush=True)
    if LOG: LOG.write(msg + '\n'); LOG.flush()


# ============================================================================== poses (exact)
def rational_angle(th_rad, Q):
    """p/q ~ tan(theta/2) with denominator Q -> exact (a, b, r): cos = a/r, sin = b/r"""
    p = round(math.tan(th_rad / 2) * Q); q = Q
    g = gcd(abs(p), q); p //= g; q //= g
    return p, q

def cos_sin(p, q):
    return q * q - p * p, 2 * p * q, q * q + p * p

def half_width(a, b, r):
    return Fr(abs(a) + abs(b), 2 * r)

def snap_pose(cx, cy, th_rad, Q, Dc):
    """float pose -> rational pose (p, q, cx, cy) with the centre clamped EXACTLY into the container"""
    p, q = rational_angle(th_rad, Q)
    a, b, r = cos_sin(p, q)
    w2 = half_width(a, b, r)
    lo, hi = w2, T - w2
    assert lo <= hi
    cx = min(max(Fr(round(cx * Dc), Dc), lo), hi)
    cy = min(max(Fr(round(cy * Dc), Dc), lo), hi)
    return p, q, cx, cy

def images(p, q, cx, cy):
    """the 8 dihedral images (cx, cy, p, q) -- angle -theta <-> -p (same convention as packing_dual)"""
    return [(cx, cy, p, q), (T - cx, cy, -p, q), (cx, T - cy, -p, q), (T - cx, T - cy, p, q),
            (cy, cx, -p, q), (T - cy, cx, p, q), (cy, T - cx, p, q), (T - cy, T - cx, -p, q)]

def make_square(cx, cy, p, q, pose):
    """exact square record: (CX, CY, Dk, a, b, r, corners[(X,Y)x4], Dq, pose)"""
    a, b, r = cos_sin(p, q)
    Dk = lcm(cx.denominator, cy.denominator)
    CX = cx.numerator * (Dk // cx.denominator); CY = cy.numerator * (Dk // cy.denominator)
    Dq = 2 * r * Dk
    corners = []
    for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):        # cyclic order
        X = 2 * r * CX + Dk * (a * sx - b * sy)
        Y = 2 * r * CY + Dk * (b * sx + a * sy)
        corners.append((X, Y))
    g = Dq
    for X, Y in corners: g = gcd(g, gcd(abs(X), abs(Y)))
    corners = [(X // g, Y // g) for X, Y in corners]; Dq //= g
    for X, Y in corners:                                       # admissibility: closed square in closed container
        if not (0 <= X and TD * X <= TN * Dq and 0 <= Y and TD * Y <= TN * Dq):     # explicit: survives python -O
            raise SupportError(f'not admissible: pose {pose}, centre ({cx}, {cy}), p/q = {p}/{q}')
    return (CX, CY, Dk, a, b, r, corners, Dq, pose)


# ============================================================================== arrangement vertices
def reduce3(X, Y, D):
    g = gcd(gcd(abs(X), abs(Y)), D)
    return (X // g, Y // g, D // g)

def intersect(e1, e2):
    """intersection point of two closed segments with rational endpoints, or None (parallel / disjoint)"""
    x1, y1, x2, y2, D1 = e1; x3, y3, x4, y4, D2 = e2
    ax = x1 * D2; ay = y1 * D2; rx = x2 * D2 - ax; ry = y2 * D2 - ay
    cx = x3 * D1; cy = y3 * D1; sx = x4 * D1 - cx; sy = y4 * D1 - cy
    det = rx * sy - ry * sx
    if det == 0: return None
    qx = cx - ax; qy = cy - ay
    tn = qx * sy - qy * sx
    un = qx * ry - qy * rx
    if det < 0: det = -det; tn = -tn; un = -un
    if tn < 0 or tn > det or un < 0 or un > det: return None
    return reduce3(ax * det + tn * rx, ay * det + tn * ry, D1 * D2 * det)

def in_F(v):
    X, Y, D = v
    return X >= 0 and X <= Y and 2 * TD * Y <= TN * D

# module globals read by the workers; a worker process gets them from _init_worker (pool initializer),
# never by inheritance, so every start method (fork / spawn / forkserver) sees the same state
EDGES = []; EMIN = []; EBUCKET = {}; SQ = []; SBUCKET = {}; VERTS = []; REGION_F = True

def pair_worker(cells):
    out = []
    for (cxi, cyi) in cells:
        L = EBUCKET.get((cxi, cyi))
        if not L: continue
        n = len(L)
        for i in range(n):
            ei = L[i]; mxi, myi = EMIN[ei]; si = ei >> 2; e1 = EDGES[ei]
            for j in range(i + 1, n):
                ej = L[j]
                if (ej >> 2) == si: continue
                mxj, myj = EMIN[ej]
                # each pair is handled once: in the cell (max of the min-cells)
                if (mxi if mxi > mxj else mxj) != cxi or (myi if myi > myj else myj) != cyi: continue
                v = intersect(e1, EDGES[ej])
                if v is not None and (not REGION_F or in_F(v)): out.append(v)
    return out

def cov_worker(rng):
    lo, hi = rng; res = []
    for vi in range(lo, hi):
        X, Y, D = VERTS[vi]
        inc = []
        for si in SBUCKET.get(((X * G_SQ) // D, (Y * G_SQ) // D), ()):
            CX, CY, Dk, a, b, r = SQ[si][:6]
            dx = X * Dk - CX * D; dy = Y * Dk - CY * D
            lim = r * D * Dk
            u = a * dx + b * dy
            if u < 0: u = -u
            if 2 * u > lim: continue
            v = a * dy - b * dx
            if v < 0: v = -v
            if 2 * v > lim: continue
            inc.append(si)
        res.append(inc)
    return lo, res

SQW = []          # per-square integer weight (check mode, streaming)
START_METHOD = None   # multiprocessing start method (None: platform default)

def _init_worker(state):
    """pool initializer: install the arrangement state (and the container side) in this worker"""
    global T
    globals().update(state)
    T = Fr(state['TN'], state['TD'])

def _pool_map(procs, func, items, state_keys, unordered=False):
    """map func over items: in-process when procs <= 1, else in a pool whose workers are initialised
    with the named module globals (plus TN, TD) -- independent of the start method"""
    if procs <= 1:
        for it in items: yield func(it)
        return
    g = globals()
    state = {k: g[k] for k in ('TN', 'TD') + tuple(state_keys)}
    ctx = mp.get_context(START_METHOD)
    with ctx.Pool(procs, initializer=_init_worker, initargs=(state,)) as pool:
        yield from (pool.imap_unordered(func, items) if unordered else pool.imap(func, items))

def covmax_worker(rng):
    """streaming variant of cov_worker for `check`: the max weighted coverage over the range, without
    keeping the incidence lists (same exact containment test, via cov_worker)"""
    lo, res = cov_worker(rng)
    best = -1; bi = -1; npairs = 0
    for off, inc in enumerate(res):
        npairs += len(inc)
        num = sum(SQW[si] for si in inc)
        if num > best: best = num; bi = lo + off
    return best, bi, npairs

def max_coverage(procs, weights):
    """(best numerator, vertex index, number of pairs) over all VERTS, weights[si] per square image"""
    global SQW
    SQW = weights
    t0 = time.time()
    n = len(VERTS); rngs = [(i, min(i + 2000, n)) for i in range(0, n, 2000)]
    best = -1; bi = -1; npairs = 0
    for b, i, np_ in _pool_map(procs, covmax_worker, rngs, ('VERTS', 'SQ', 'SBUCKET', 'SQW'), unordered=True):
        npairs += np_
        if b > best or (b == best and i < bi): best = b; bi = i
    log(f"  incidences (streamed): {n} vertices, {npairs} (vertex, square) pairs, {time.time() - t0:.1f} s")
    return best, bi

def chunks(lst, n):
    k = max(1, (len(lst) + n - 1) // n)
    return [lst[i:i + k] for i in range(0, len(lst), k)]

def enumerate_vertices(squares, full, procs):
    """all arrangement vertices in F (or the whole container if full) as reduced (X, Y, D)"""
    global EDGES, EMIN, EBUCKET, SQ, SBUCKET, VERTS, REGION_F
    REGION_F = not full
    SQ = squares
    K = (TN * G_EDGE) // TD if full else (TN * G_EDGE) // (2 * TD)      # largest cell index of the region
    EDGES = []; EMIN = []; EBUCKET = {}
    for si, s in enumerate(squares):
        corners, Dq = s[6], s[7]
        for k in range(4):
            (x1, y1), (x2, y2) = corners[k], corners[(k + 1) % 4]
            ei = len(EDGES); EDGES.append((x1, y1, x2, y2, Dq))
            cx0 = (min(x1, x2) * G_EDGE) // Dq; cx1 = (max(x1, x2) * G_EDGE) // Dq
            cy0 = (min(y1, y2) * G_EDGE) // Dq; cy1 = (max(y1, y2) * G_EDGE) // Dq
            EMIN.append((cx0, cy0))
            for cx in range(cx0, min(cx1, K) + 1):
                for cy in range(cy0, min(cy1, K) + 1):
                    EBUCKET.setdefault((cx, cy), []).append(ei)
    cells = sorted(EBUCKET.keys())
    t0 = time.time()
    parts = list(_pool_map(procs, pair_worker, chunks(cells, 8 * procs), ('EDGES', 'EMIN', 'EBUCKET', 'REGION_F')))
    verts = set()
    for part in parts: verts.update(part)
    n_int = len(verts)
    n_corner = 0
    for s in squares:
        Dq = s[7]
        for X, Y in s[6]:
            v = reduce3(X, Y, Dq)
            if not REGION_F or in_F(v): verts.add(v); n_corner += 1
    log(f"  vertices: {n_int} distinct edge intersections + {n_corner} corners -> {len(verts)} distinct "
        f"({len(EDGES)} edges, {len(cells)} cells, {time.time() - t0:.1f} s)")
    VERTS = sorted(verts)
    # square buckets for the coverage evaluation
    SBUCKET = {}
    for si, s in enumerate(squares):
        corners, Dq = s[6], s[7]
        xs = [X for X, Y in corners]; ys = [Y for X, Y in corners]
        for cx in range((min(xs) * G_SQ) // Dq, (max(xs) * G_SQ) // Dq + 1):
            for cy in range((min(ys) * G_SQ) // Dq, (max(ys) * G_SQ) // Dq + 1):
                SBUCKET.setdefault((cx, cy), []).append(si)
    return VERTS

def incidences(procs):
    """for every vertex, the list of square (image) indices containing it -- exact"""
    t0 = time.time()
    n = len(VERTS); rngs = [(i, min(i + 2000, n)) for i in range(0, n, 2000)]
    out = [None] * n
    for lo, res in _pool_map(procs, cov_worker, rngs, ('VERTS', 'SQ', 'SBUCKET'), unordered=True):
        out[lo:lo + len(res)] = res
    log(f"  incidences: {n} vertices, {sum(map(len, out))} (vertex, square) pairs, {time.time() - t0:.1f} s")
    return out

def patterns(inc, pose_of):
    """distinct incidence patterns: tuple(sorted((pose, count))) -> list of vertex indices"""
    pat = {}
    for vi, lst in enumerate(inc):
        key = tuple(sorted(Counter(pose_of[s] for s in lst).items()))
        pat.setdefault(key, []).append(vi)
    return pat

def exact_max(pat, MU, DM):
    """M = max_v cov(v) as a Fraction, with masses MU[k]/DM; returns (M, key attaining it)"""
    best = None; bkey = None
    for key in pat:
        num = sum(c * MU[k] for k, c in key)
        if best is None or num > best: best = num; bkey = key
    return Fr(best, 8 * DM), bkey


# ============================================================================== support files
def read_float_support(path):
    poses = []
    for line in open(path):
        if line.startswith('#') or not line.strip(): continue
        q = line.split()
        if q[0] == 'pose': poses.append((float(q[1]), float(q[2]), math.radians(float(q[3])), float(q[4])))
    return poses

def support_t(path):
    """container side recorded in the header of an exact support file, or None"""
    for line in open(path):
        m = re.match(r'#\s*t\s*=\s*(-?\d+(?:/\d+)?)', line)
        if m: return Fr(m.group(1))
    return None

class SupportError(ValueError):
    pass

_INT = re.compile(r'[+-]?\d+$')
_FRAC = re.compile(r'[+-]?\d+(/\d+)?$')

def read_exact_support(path):
    """strict parser: every non-blank, non-comment line must be `pose p q cx cy mass` with integers
    p, q (q > 0) and exact fractions cx, cy, mass (mass >= 0); anything else raises SupportError"""
    poses = []
    for ln, line in enumerate(open(path), 1):
        if line.startswith('#') or not line.strip(): continue
        q = line.split()
        where = f"{path}:{ln}"
        if q[0] != 'pose' or len(q) != 6:
            raise SupportError(f"{where}: expected 'pose p q cx cy mass', got {line.strip()!r}")
        if not (_INT.match(q[1]) and _INT.match(q[2])) or not all(_FRAC.match(x) for x in q[3:6]):
            raise SupportError(f"{where}: p, q must be integers and cx, cy, mass fractions: {line.strip()!r}")
        try:
            cx, cy, m = Fr(q[3]), Fr(q[4]), Fr(q[5])
        except ZeroDivisionError:
            raise SupportError(f"{where}: zero denominator: {line.strip()!r}")
        p_, q_ = int(q[1]), int(q[2])
        if q_ <= 0: raise SupportError(f"{where}: q must be positive: {line.strip()!r}")
        if m < 0: raise SupportError(f"{where}: negative mass: {line.strip()!r}")
        poses.append((p_, q_, cx, cy, m))
    if not poses: raise SupportError(f"{path}: no pose lines")
    return poses

def write_exact_support(path, poses, MU, DM, M, L, note):
    with open(path, 'w') as f:
        f.write(f"# t = {TN}/{TD} exact; D4-symmetrised measure of closed unit squares; {note}\n")
        f.write(f"# mass = {sum(MU)}/{DM} = {sum(MU) / DM:.12f};  M = max coverage = {M} = {float(M):.15f};  "
                f"L = mass/M = {L} = {float(L):.12f}\n")
        f.write("# pose p q cx cy mass : rotation theta = 2*arctan(p/q) (cos = (q^2-p^2)/(q^2+p^2), sin = 2pq/(q^2+p^2)),\n"
                "#   rational centre (cx, cy); mass/8 on each of the 8 dihedral images "
                "(cx,cy,+), (t-cx,cy,-), (cx,t-cy,-), (t-cx,t-cy,+), (cy,cx,-), (t-cy,cx,+), (cy,t-cx,+), (t-cy,t-cx,-)\n")
        for (p, q, cx, cy), m in zip(poses, MU):
            if m > 0: f.write(f"pose {p} {q} {cx} {cy} {Fr(m, DM)}\n")

def build_squares(poses):
    squares = []; pose_of = []
    for k, (p, q, cx, cy) in enumerate(poses):
        for (x, y, pp, qq) in images(p, q, cx, cy):
            squares.append(make_square(x, y, pp, qq, k)); pose_of.append(k)
    return squares, pose_of


# ============================================================================== modes
def polish_lp(pat, poses, DM):
    """LP (floats, HiGHS) on the exact rows -- only CHOOSES the masses; returns (MU, exact M, worst key)"""
    import numpy as np, scipy.sparse as sp
    from scipy.optimize import linprog
    keys = list(pat.keys()); nP = len(poses)
    rows = []; cols = []; vals = []
    for i, key in enumerate(keys):
        for k, c in key: rows.append(i); cols.append(k); vals.append(c / 8.0)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(len(keys), nP))
    t0 = time.time()
    res = linprog(-np.ones(nP), A_ub=A, b_ub=np.ones(len(keys)), bounds=(0, None), method='highs',
                  options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
    assert res.status == 0, res.message
    log(f"  LP: {len(keys)} rows x {nP} cols, {A.nnz} nz, mass {-res.fun:.9f}, {time.time() - t0:.1f} s")
    mu = np.maximum(res.x, 0.0)
    MU = [int(math.floor(float(m) * DM)) for m in mu]
    M, key = exact_max(pat, MU, DM)
    log(f"  rounded down to /{DM}: mass {sum(MU)/DM:.9f}, exact M = {M} = {float(M):.12f}")
    return MU, M, key


def cmd_build(args):
    global LOG
    LOG = open(os.path.join(RUNS, f'dual_exact_{TAG}.log'), 'w')
    t_all = time.time()
    src = read_float_support(args.src)
    log(f"build: {len(src)} float poses from {args.src}; snap Q={args.Q} (angle), Dc={args.Dc} (centre); t={TN}/{TD}")
    poses = []; dmax = 0.0
    for cx, cy, th, mu in src:
        p, q, x, y = snap_pose(cx, cy, th, args.Q, args.Dc)
        a, b, r = cos_sin(p, q)
        dth = abs(math.atan2(b, a) - th); dmax = max(dmax, dth, abs(float(x) - cx), abs(float(y) - cy))
        poses.append((p, q, x, y))
    log(f"  max snap displacement (angle rad / centre): {dmax:.2e}")
    squares, pose_of = build_squares(poses)
    log(f"  {len(squares)} images, all admissible (exact)")
    verts = enumerate_vertices(squares, full=False, procs=args.procs)
    inc = incidences(args.procs)
    pat = patterns(inc, pose_of)
    log(f"  {len(pat)} distinct incidence patterns (LP rows)")

    # ---- masses of the float support, rounded down: exact M before polishing
    DM = args.DM
    MU0 = [int(math.floor(mu * DM)) for (_, _, _, mu) in src]
    M0, _ = exact_max(pat, MU0, DM)
    log(f"  before polish: mass {sum(MU0)/DM:.9f}, exact M = {float(M0):.9f}, L = {sum(MU0)/DM/float(M0):.9f}")

    if args.no_lp:
        # keep the float masses (rounded down); skips the LP, which dominates memory on large supports
        MU = MU0; M, key = exact_max(pat, MU, DM)
        log(f"  --no-lp: float masses rounded down to /{DM}: mass {sum(MU)/DM:.9f}, exact M = {M} = {float(M):.12f}")
    else:
        MU, M, key = polish_lp(pat, poses, DM)
    if M > 1:
        MU = [(m * M.denominator) // M.numerator for m in MU]      # floor(mu / M): coverage <= 1 exactly
        M, key = exact_max(pat, MU, DM)
        log(f"  scaled by 1/M and rounded down: mass {sum(MU)/DM:.9f}, exact M = {M} = {float(M):.12f}")
    assert M <= 1
    mass = Fr(sum(MU), DM); L = mass / M
    log(f"  EXACT: mass = {mass} = {float(mass):.12f};  M = {M};  L = mass/M = {L} = {float(L):.12f}  "
        f"({'>=' if L >= args.n else '<'} {args.n}; L - {args.n} = {float(L - args.n):+.3e})")
    log(f"  support: {sum(1 for m in MU if m > 0)} poses with positive mass; worst vertex pattern {key}")
    out = os.path.join(RUNS, f'dual_exact_{TAG}_support.txt')
    write_exact_support(out, poses, MU, DM, M, L, f"snapped from {os.path.basename(args.src)} with Q={args.Q}, Dc={args.Dc}")
    js = dict(t=f"{TN}/{TD}", Q=args.Q, Dc=args.Dc, DM=DM, poses=len(poses), poses_positive=sum(1 for m in MU if m > 0),
              images=len(squares), vertices_F=len(verts), incidence_pairs=sum(map(len, inc)), lp_rows=len(pat),
              mass=str(mass), mass_float=float(mass), M=str(M), M_float=float(M), L=str(L), L_float=float(L),
              L_minus_12=str(L - 12), M_before_polish=str(M0), seconds=time.time() - t_all)
    json.dump(js, open(os.path.join(RUNS, f'dual_exact_{TAG}.json'), 'w'), indent=1)
    log(f"  wrote {out}; total {time.time() - t_all:.1f} s")
    return L

def cmd_check(args, sup):
    global LOG
    tag = 'full' if args.full else 'F'
    os.makedirs(RUNS, exist_ok=True)
    LOG = open(os.path.join(RUNS, f'dual_exact_{TAG}_check_{tag}.log'), 'w')
    t_all = time.time()
    poses = [(p, q, cx, cy) for (p, q, cx, cy, m) in sup]
    masses = [m for (_, _, _, _, m) in sup]
    DM = 1
    for m in masses: DM = lcm(DM, m.denominator)
    MU = [m.numerator * (DM // m.denominator) for m in masses]
    log(f"check ({'whole container' if args.full else 'fundamental domain'}): {len(poses)} exact poses from {args.support}, "
        f"t = {TN}/{TD}, mass denominator {DM}")
    squares, pose_of = build_squares(poses)
    log(f"  {len(squares)} images, all admissible (exact: 4 rational corners in the closed container)")
    verts = enumerate_vertices(squares, full=args.full, procs=args.procs)
    if args.stream:
        best, bi = max_coverage(args.procs, [MU[k] for k in pose_of])
    else:
        inc = incidences(args.procs)
        # exact coverage at every vertex
        best = -1; bi = -1
        for vi, lst in enumerate(inc):
            num = sum(MU[pose_of[s]] for s in lst)
            if num > best: best = num; bi = vi
    mass = Fr(sum(MU), DM)
    if best <= 0 or mass <= 0:
        raise SupportError(f"degenerate support: total mass {mass}, max coverage numerator {best}")
    M = Fr(best, 8 * DM); L = mass / M
    X, Y, D = verts[bi]
    n = args.n if args.n is not None else 12            # display target when nothing is asserted
    log(f"  EXACT: {len(verts)} vertices checked; mass = {mass} = {float(mass):.12f};  M = {M} = {float(M):.15f} "
        f"(attained at ({X}/{D}, {Y}/{D}) ~ ({X/D:.6f}, {Y/D:.6f}));  L = mass/M = {L} = {float(L):.12f}  "
        f"({'>=' if L >= n else '<'} {n});  {time.time() - t_all:.1f} s")
    if args.n is None:
        log("  no --n given: nothing asserted (exit status 0 only means the check ran)")
        return 0
    ok = L > args.n
    log(f"  ASSERT L > {args.n}: {'PASS' if ok else 'FAIL'} (L - {args.n} = {L - args.n} = "
        f"{float(L - args.n):+.6e}, exact)")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    b = sub.add_parser('build')
    b.add_argument('--t', default='399/100', help='container side (exact fraction, e.g. 4 or 399/100)')
    b.add_argument('--tag', default=None, help='output file tag (default: from t)')
    b.add_argument('--src', default=os.path.join(RUNS, 'dual_PA2_support.txt'))
    b.add_argument('--Q', type=int, default=100000, help='denominator of p/q ~ tan(theta/2)')
    b.add_argument('--Dc', type=int, default=1000000, help='denominator of the centre coordinates')
    b.add_argument('--DM', type=int, default=10 ** 9, help='denominator of the masses')
    b.add_argument('--procs', type=int, default=4)
    b.add_argument('--no-lp', action='store_true', help='keep the float masses (rounded down, scaled by 1/M if M > 1); no LP polish')
    b.add_argument('--n', type=int, default=12, help='target the log line compares L with (display only)')
    c = sub.add_parser('check')
    c.add_argument('support')
    c.add_argument('--t', default=None, help='container side (default: read from the support file header)')
    c.add_argument('--tag', default=None)
    c.add_argument('--full', action='store_true', help='whole container instead of the fundamental domain')
    c.add_argument('--stream', action='store_true', help='do not keep the incidence lists (low memory; same exact test)')
    c.add_argument('--procs', type=int, default=4)
    c.add_argument('--n', type=Fr, default=None,
                   help='ASSERT the strict exact inequality L > N: exit 0 iff it holds, 1 otherwise '
                        '(without --n nothing is asserted; the log line compares with 12 for display)')
    for q in (b, c):
        q.add_argument('--start-method', choices=('fork', 'spawn', 'forkserver'), default=None,
                       help='multiprocessing start method (default: platform default; all are supported)')
    args = ap.parse_args()
    global START_METHOD
    START_METHOD = args.start_method
    sys.set_int_max_str_digits(0)
    if args.cmd == 'build':
        set_t(Fr(args.t), args.tag)
        cmd_build(args)
    else:
        try:
            th = support_t(args.support)
            t = Fr(args.t) if args.t else th
            if t is None: raise SupportError('no "# t = p/q" header in the support file; pass --t')
            if th is not None and th != t:
                raise SupportError(f'--t {t} disagrees with the support header t = {th}')
            if t <= 0: raise SupportError(f't = {t} must be positive')
            set_t(t, args.tag)
            sup = read_exact_support(args.support)
            rc = cmd_check(args, sup)
        except (SupportError, OSError, ValueError) as e:
            print(f"error: {e}", file=sys.stderr); sys.exit(2)
        sys.exit(rc)

if __name__ == '__main__':
    main()
