#!/usr/bin/env python3
"""unavoid13: hitting-set cutting-plane loop for pure unavoidable sets of closed unit squares
in [0,m]^2 (tasks/unavoid13).  Library; drivers are unavoid13_loop.py / unavoid13_recheck.py.

Semantics (search/ZEROMARGIN.md s1): pose (cx, cy, u = tan(theta/2)), closed unit square
Q = c + R_theta [-1/2,1/2]^2, admissible iff cx, cy in [w/2, m - w/2], w = |cos| + |sin|; a point
on the boundary of Q counts as inside.

Lower bound: for a finite family F of admissible squares, h(F) = min |P| with P hitting every
Q in F is a lower bound on the pure unavoidable number.  Candidate points = the vertices of the
arrangement of F (square corners + intersections of non-parallel edges of two squares): for any
point p, the polygon  P(p) = intersection of the squares of F containing p  is compact convex and
contains p; each of its vertices is such an arrangement vertex and lies in every square containing
p, so every point is dominated by an arrangement vertex and nothing is lost.

Float incidence is LENIENT (tol): the computed incidence is a superset of the true one, so the
float IP value is <= the true h(F) -- the safe direction for a lower bound -- provided the vertex
list is complete.  The theorem-grade path (unavoid13_recheck.py) recomputes vertices and incidence
exactly in Fractions from the dumped family.
"""
import math, os, sys, time, json
from fractions import Fraction as F
import numpy as np

HALF = F(1, 2)
SQ2 = math.sqrt(2.0)

# ----------------------------------------------------------------------------------------------
# exact poses
# ----------------------------------------------------------------------------------------------
def trig(u):
    """cos, sin of theta = 2 atan u, exact (u a Fraction)."""
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d


def norm_u(u):
    """Normalise u = tan(theta/2) into [0, 1) (theta in [0, 90 deg)); the square is 90-deg symmetric."""
    u = F(u)
    while u < 0:
        u = (1 + u) / (1 - u)
    while u >= 1:
        u = (u - 1) / (u + 1)
    return u


def pose_key(u, cx, cy):
    return (norm_u(u), F(cx), F(cy))


def admissible_exact(m, u, cx, cy):
    c, s = trig(u)
    w = abs(c) + abs(s)
    return (w / 2 <= cx <= m - w / 2) and (w / 2 <= cy <= m - w / 2)


def d4_images(m, u, cx, cy):
    """The 8 dihedral images of the pose as normalised keys (deduped, admissible ones are all
    admissible since the container is D4-invariant)."""
    m = F(m); u = norm_u(u); cx = F(cx); cy = F(cy)
    out = set()
    # reflections x -> m - x send theta -> -theta
    for (x, y, uu) in ((cx, cy, u), (m - cx, cy, norm_u(-u))):
        # rotations by 90 deg about the container centre: (x, y) -> (m - y, x)
        for _ in range(4):
            out.add(pose_key(uu, x, y))
            x, y = m - y, x
    return sorted(out)


def d4_closure(m, poses):
    out = set()
    for (u, cx, cy) in poses:
        out.update(d4_images(m, u, cx, cy))
    return sorted(out)


def square_corners_exact(u, cx, cy):
    c, s = trig(u)
    e1 = (c, s); e2 = (-s, c)
    out = []
    for a, b in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)):
        out.append((cx + a * e1[0] + b * e2[0], cy + a * e1[1] + b * e2[1]))
    return out


def inside_exact(u, cx, cy, px, py):
    c, s = trig(u)
    dx, dy = px - cx, py - cy
    x = dx * c + dy * s; y = -dx * s + dy * c
    return -HALF <= x <= HALF and -HALF <= y <= HALF


def depth_exact(u, cx, cy, px, py):
    """1/2 - max(|x'|, |y'|): >= 0 iff inside (closed)."""
    c, s = trig(u)
    dx, dy = px - cx, py - cy
    x = dx * c + dy * s; y = -dx * s + dy * c
    return HALF - max(abs(x), abs(y))


# ----------------------------------------------------------------------------------------------
# family I/O
# ----------------------------------------------------------------------------------------------
def write_family(path, m, poses, header=''):
    with open(path, 'w') as f:
        f.write(f"# m = {m}; pose u cx cy (Fractions; u = tan(theta/2), theta in [0,90deg)); "
                f"closed unit square, closed containment\n")
        if header:
            for line in header.split('\n'):
                f.write('# ' + line + '\n')
        for (u, cx, cy) in poses:
            f.write(f"pose {u} {cx} {cy}\n")


def read_family(path):
    poses = []; m = None
    for line in open(path):
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('# m =') and m is None:
                m = F(line.split('=')[1].split(';')[0].strip())
            continue
        if not line: continue
        tok = line.split()
        if tok[0] == 'pose':
            poses.append(pose_key(F(tok[1]), F(tok[2]), F(tok[3])))
    return m, poses


def read_support(path, m):
    """search/cover4_exact_support.txt style: 'pose p q cx cy mass' (theta = 2 atan(p/q))."""
    out = []
    for line in open(path):
        tok = line.split()
        if not tok or tok[0] != 'pose': continue
        p, q = int(tok[1]), int(tok[2])
        u = F(p, q); cx = F(tok[3]); cy = F(tok[4]); mass = F(tok[5])
        out.append((pose_key(u, cx, cy), mass))
    return out


# ----------------------------------------------------------------------------------------------
# float squares, arrangement vertices, incidence
# ----------------------------------------------------------------------------------------------
class Squares:
    def __init__(self, poses):
        self.poses = list(poses)
        n = len(self.poses)
        self.cx = np.zeros(n); self.cy = np.zeros(n); self.c = np.zeros(n); self.s = np.zeros(n)
        for i, (u, cx, cy) in enumerate(self.poses):
            c, s = trig(u)
            self.cx[i] = float(cx); self.cy[i] = float(cy); self.c[i] = float(c); self.s[i] = float(s)
        # corners (n, 4, 2), counter-clockwise
        e1 = np.stack([self.c, self.s], 1); e2 = np.stack([-self.s, self.c], 1)
        ctr = np.stack([self.cx, self.cy], 1)
        self.corners = np.stack([ctr - .5 * e1 - .5 * e2, ctr + .5 * e1 - .5 * e2,
                                 ctr + .5 * e1 + .5 * e2, ctr - .5 * e1 + .5 * e2], 1)

    def __len__(self):
        return len(self.poses)

    def depth(self, X, chunk=None):
        """X (M,2) points -> (M, n) depths (1/2 - max(|x'|,|y'|)), float."""
        dx = X[:, 0, None] - self.cx[None, :]
        dy = X[:, 1, None] - self.cy[None, :]
        xp = dx * self.c[None, :] + dy * self.s[None, :]
        yp = -dx * self.s[None, :] + dy * self.c[None, :]
        return 0.5 - np.maximum(np.abs(xp), np.abs(yp))


def arrangement_vertices(sq, tol=1e-9, par_tol=1e-11):
    """Square corners + intersections of every pair of non-parallel edges of two different squares
    (float).  Returns (V, 2).  par_tol: |cross(r, q)| below this counts as parallel; the exact
    recheck does not use this cut-off."""
    n = len(sq)
    C = sq.corners
    ctr = np.stack([sq.cx, sq.cy], 1)
    # candidate pairs: centres within sqrt2 (two unit squares meet only if centres <= sqrt2 apart)
    I, J = np.triu_indices(n, 1)
    d2 = ((ctr[I] - ctr[J]) ** 2).sum(1)
    keep = d2 <= 2.0 + 1e-9
    I, J = I[keep], J[keep]
    pts = [C.reshape(-1, 2)]
    for a in range(4):
        P = C[I, a]; R = C[I, (a + 1) % 4] - P
        for b in range(4):
            Q = C[J, b]; S = C[J, (b + 1) % 4] - Q
            den = R[:, 0] * S[:, 1] - R[:, 1] * S[:, 0]
            ok = np.abs(den) > par_tol
            if not ok.any(): continue
            PQ = Q - P
            t = (PQ[:, 0] * S[:, 1] - PQ[:, 1] * S[:, 0])
            s = (PQ[:, 0] * R[:, 1] - PQ[:, 1] * R[:, 0])
            with np.errstate(divide='ignore', invalid='ignore'):
                t = np.where(ok, t / den, -1.0); s = np.where(ok, s / den, -1.0)
            ok &= (t >= -tol) & (t <= 1 + tol) & (s >= -tol) & (s <= 1 + tol)
            if ok.any():
                pts.append(P[ok] + t[ok, None] * R[ok])
    return np.concatenate(pts, 0)


def incidence_packed(sq, V, tol=1e-7, chunk=4000):
    """Lenient float incidence of vertices V in squares: packed bit rows (M, ceil(n/8))."""
    n = len(sq); M = len(V)
    out = np.zeros((M, (n + 7) // 8), dtype=np.uint8)
    for i0 in range(0, M, chunk):
        D = sq.depth(V[i0:i0 + chunk])
        out[i0:i0 + chunk] = np.packbits(D >= -tol, axis=1)
    return out


def dedupe_rows(packed):
    """Unique packed rows; returns (unique_rows, first_index)."""
    uniq, idx = np.unique(packed, axis=0, return_index=True)
    return uniq, idx


def unpack_rows(packed, n):
    return np.unpackbits(packed, axis=1)[:, :n].astype(bool)


def drop_dominated(B):
    """B (K, n) bool incidence.  Remove rows that are subsets of another row (dominated
    candidates).  O(K^2) bit ops in numpy blocks; use only for K up to a few 10^4."""
    K, n = B.shape
    if K == 0: return np.zeros(0, dtype=np.int64)
    P = np.packbits(B, axis=1)
    cnt = B.sum(1)
    order = np.argsort(-cnt, kind='stable')
    P = P[order]; cnt = cnt[order]
    keep = np.ones(K, dtype=bool)
    # a row i is dominated if some row j with cnt[j] > cnt[i] has P[j] & P[i] == P[i]
    blk = 2000
    for i0 in range(0, K, blk):
        Pi = P[i0:i0 + blk]                         # (b, bytes)
        ci = cnt[i0:i0 + blk]
        # only rows with strictly more bits can dominate: those are earlier in the order
        jmax = np.searchsorted(-cnt, -ci.min(), side='left')   # rows with cnt > min(ci)
        if jmax == 0: continue
        Pj = P[:jmax]; cj_all = cnt[:jmax]
        dom = np.zeros(len(Pi), dtype=bool)
        for j0 in range(0, jmax, 4000):
            Bj = Pj[j0:j0 + 4000]; cj = cj_all[j0:j0 + 4000]
            # Pi is a subset of Bj  iff  (Pi & ~Bj) == 0 on every byte
            sub = ((Pi[:, None, :] & ~Bj[None, :, :]) == 0).all(2)      # (b, bj)
            # exclude equal-count rows (equal count + superset = equal row; rows are unique)
            sub &= (cj[None, :] > ci[:, None])
            dom |= sub.any(1)
        keep[i0:i0 + blk] &= ~dom
    return np.sort(order[keep])


# ----------------------------------------------------------------------------------------------
# the hitting-set IP
# ----------------------------------------------------------------------------------------------
def build_candidates(m, poses, tol=1e-7, dominance_max=2000000, verbose=True):
    """poses -> (sq, V, B, rep): candidate incidence matrix B (K, n) bool with representative
    vertex coordinates rep (K, 2)."""
    t0 = time.time()
    sq = Squares(poses)
    V = arrangement_vertices(sq)
    # keep vertices inside the container (all are, up to float noise) -- clamp
    V = np.clip(V, 0.0, float(m))
    packed = incidence_packed(sq, V, tol=tol)
    uniq, idx = dedupe_rows(packed)
    B = unpack_rows(uniq, len(sq))
    nz = B.any(1)
    B = B[nz]; idx = idx[nz]
    ndom = 0
    if len(B) <= dominance_max:
        keep = drop_dominated_fast(B)
        ndom = len(B) - len(keep)
        B = B[keep]; idx = idx[keep]
    rep = V[idx]
    if verbose:
        print(f"  candidates: {len(sq)} squares, {len(V)} vertices, {len(uniq)} distinct incidence "
              f"sets, {ndom} dominated dropped, {len(B)} candidates, {time.time()-t0:.1f}s", flush=True)
    return sq, V, B, rep


def solve_hitting_set(B, threads=4, time_limit=None, verbose=False, ub=None, warm=None,
                      extra_rows=None, symmetry_off=False):
    """min sum x  s.t.  B^T x >= 1 (every square hit), x binary.  B (K, n) bool: K candidates.
    Returns dict(status, obj, x, lp) with lp = LP relaxation value (root, from a separate LP solve)."""
    import highspy
    K, n = B.shape
    # column-wise sparse: column k has rows where B[k] true
    starts = np.zeros(K + 1, dtype=np.int64)
    cnt = B.sum(1)
    starts[1:] = np.cumsum(cnt)
    index = np.concatenate([np.nonzero(B[k])[0] for k in range(K)]).astype(np.int32) if K else np.zeros(0, np.int32)
    value = np.ones(len(index))

    def make(integer):
        h = highspy.Highs()
        h.setOptionValue('output_flag', verbose)
        h.setOptionValue('threads', int(threads))
        h.setOptionValue('mip_rel_gap', 0.0)
        h.setOptionValue('mip_abs_gap', 0.0)
        h.setOptionValue('mip_feasibility_tolerance', 1e-9)
        h.setOptionValue('primal_feasibility_tolerance', 1e-9)
        if symmetry_off:
            h.setOptionValue('mip_detect_symmetry', False)
        if time_limit:
            h.setOptionValue('time_limit', float(time_limit))
        lp = highspy.HighsLp()
        lp.num_col_ = K; lp.num_row_ = n
        lp.col_cost_ = np.ones(K); lp.col_lower_ = np.zeros(K); lp.col_upper_ = np.ones(K)
        lp.row_lower_ = np.ones(n); lp.row_upper_ = np.full(n, highspy.kHighsInf)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = value
        if integer:
            lp.integrality_ = [highspy.HighsVarType.kInteger] * K
        h.passModel(lp)
        return h

    # LP relaxation
    hl = make(False); hl.run()
    lpval = hl.getInfo().objective_function_value
    # IP
    h = make(True)
    if ub is not None:
        # objective cutoff: any solution with objective >= ub is not wanted (proves >= ub if infeasible)
        h.setOptionValue('objective_bound', float(ub) - 0.5)
    if warm is not None:
        sol = highspy.HighsSolution()
        sol.col_value = list(map(float, warm)); sol.value_valid = True
        h.setSolution(sol)
    h.run()
    st = h.getModelStatus()
    info = h.getInfo()
    res = dict(status=h.modelStatusToString(st), lp=lpval, obj=None, x=None,
               bound=info.mip_dual_bound, gap=info.mip_gap)
    if st == highspy.HighsModelStatus.kOptimal:
        x = np.array(h.getSolution().col_value)
        res['obj'] = int(round(info.objective_function_value)); res['x'] = x
    return res


def solve_hitting_set_scipy(B, time_limit=None):
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import csc_matrix
    K, n = B.shape
    A = csc_matrix(B.T.astype(float))
    c = np.ones(K)
    res = milp(c, constraints=LinearConstraint(A, lb=np.ones(n), ub=np.full(n, np.inf)),
               integrality=np.ones(K), bounds=Bounds(0, 1),
               options=dict(mip_rel_gap=0.0, time_limit=time_limit, disp=False))
    return res


# ----------------------------------------------------------------------------------------------
# violation search (float) + exact confirmation
# ----------------------------------------------------------------------------------------------
def w_of_theta(th):
    return abs(math.cos(th)) + abs(math.sin(th))


def maxdepth_grid(P, m, dth=math.radians(0.5), pitch=0.01, th_max=math.pi / 2):
    """Grid scan of f(pose) = max_p depth(p, pose) over admissible poses.  Returns arrays
    (cx, cy, th, f) of all grid poses (can be large)."""
    P = np.asarray(P, float)
    out = []
    ths = np.arange(0.0, th_max, dth)
    for th in ths:
        w = w_of_theta(th)
        lo, hi = w / 2, m - w / 2
        if hi < lo: continue
        g = np.arange(lo, hi + 1e-12, pitch)
        if len(g) == 0 or g[-1] < hi - 1e-12: g = np.append(g, hi)
        CX, CY = np.meshgrid(g, g, indexing='ij')
        cx = CX.ravel(); cy = CY.ravel()
        c, s = math.cos(th), math.sin(th)
        dx = P[None, :, 0] - cx[:, None]; dy = P[None, :, 1] - cy[:, None]
        xp = dx * c + dy * s; yp = -dx * s + dy * c
        f = (0.5 - np.maximum(np.abs(xp), np.abs(yp))).max(1)
        out.append(np.stack([cx, cy, np.full_like(cx, th), f], 1))
    return np.concatenate(out, 0)


def f_pose(P, m, cx, cy, th, clamp=True):
    w = w_of_theta(th)
    lo, hi = w / 2, m - w / 2
    if clamp:
        cx = min(max(cx, lo), hi); cy = min(max(cy, lo), hi)
    c, s = math.cos(th), math.sin(th)
    dx = P[:, 0] - cx; dy = P[:, 1] - cy
    xp = dx * c + dy * s; yp = -dx * s + dy * c
    return (0.5 - np.maximum(np.abs(xp), np.abs(yp))).max(), cx, cy


def refine_pose(P, m, cx, cy, th, iters=3):
    """Local minimisation of f over (cx, cy, theta) with clamping to the admissible range."""
    from scipy.optimize import minimize
    P = np.asarray(P, float)
    x0 = np.array([cx, cy, th])
    best = (f_pose(P, m, cx, cy, th)[0], cx, cy, th)
    for _ in range(iters):
        r = minimize(lambda z: f_pose(P, m, z[0], z[1], z[2])[0], x0, method='Nelder-Mead',
                     options=dict(xatol=1e-12, fatol=1e-13, maxiter=4000, initial_simplex=None))
        f, cxx, cyy = f_pose(P, m, r.x[0], r.x[1], r.x[2])
        if f < best[0] - 1e-15:
            best = (f, cxx, cyy, r.x[2]); x0 = np.array([cxx, cyy, r.x[2]])
        else:
            break
    return best


def find_violations(P, m, dth_deg=0.5, pitch=0.01, top=400, min_sep=0.05, seeds=None,
                    cutoff=-1e-9, want=12, th_max=math.pi / 2, verbose=True):
    """Float search for admissible poses with f < cutoff (no point of P inside, with margin).
    Returns list of (f, cx, cy, th) local minima, sorted by f, separated by min_sep in pose."""
    P = np.asarray(P, float)
    G = maxdepth_grid(P, m, math.radians(dth_deg), pitch, th_max=th_max)
    order = np.argsort(G[:, 3])
    cand = []
    for k in order[:top * 5]:
        cx, cy, th, f = G[k]
        if all(abs(cx - a) + abs(cy - b) + abs(th - t) > min_sep for (_, a, b, t) in cand):
            cand.append((f, cx, cy, th))
        if len(cand) >= top: break
    if seeds:
        for (cx, cy, th) in seeds:
            cand.append((f_pose(P, m, cx, cy, th)[0], cx, cy, th))
    res = []
    for (f0, cx, cy, th) in cand:
        f, cxx, cyy, thh = refine_pose(P, m, cx, cy, th)
        if f < cutoff:
            if all(abs(cxx - a) + abs(cyy - b) + abs(thh - t) > min_sep for (_, a, b, t) in res):
                res.append((f, cxx, cyy, thh))
    res.sort()
    if verbose:
        print(f"  violation search: grid min f = {G[:,3].min():+.3e}; {len(res)} distinct local "
              f"minima below {cutoff:g}" + (f"; worst {res[0][0]:+.3e}" if res else ''), flush=True)
    return res[:want], float(G[:, 3].min())


def rationalise_pose(m, cx, cy, th, den=10 ** 6):
    """Snap a float pose to a rational one: u = tan(theta/2) -> Fraction with denominator <= den,
    centre likewise, then clamp EXACTLY into the admissible range."""
    m = F(m)
    u = norm_u(F(math.tan(th / 2)).limit_denominator(den))
    c, s = trig(u); w = abs(c) + abs(s)
    lo, hi = w / 2, m - w / 2
    cxr = F(cx).limit_denominator(den); cyr = F(cy).limit_denominator(den)
    cxr = min(max(cxr, lo), hi); cyr = min(max(cyr, lo), hi)
    return pose_key(u, cxr, cyr)


def exact_violation(m, Pex, u, cx, cy):
    """Exact test: is the (rational) pose admissible and does it avoid every point of Pex?
    Returns (admissible, hit_count, min_over_points_of_max(|x'|,|y'|) - 1/2)."""
    m = F(m)
    if not admissible_exact(m, u, cx, cy):
        return False, None, None
    c, s = trig(u)
    hits = 0; worst = None
    for (px, py) in Pex:
        dx, dy = px - cx, py - cy
        x = dx * c + dy * s; y = -dx * s + dy * c
        d = HALF - max(abs(x), abs(y))
        if d >= 0: hits += 1
        worst = d if worst is None else max(worst, d)
    return True, hits, worst


def snap_points(P, den):
    return [(F(round(x * den), den), F(round(y * den), den)) for x, y in np.asarray(P, float)]


def write_cert(path, m, Pex, header=None):
    """certificates/FORMAT.md, unit weights (W = 1).  Pex: exact rational points."""
    D = 1
    for x, y in Pex:
        D = D * x.denominator // math.gcd(D, x.denominator)
        D = D * y.denominator // math.gcd(D, y.denominator)
    m = F(m)
    with open(path, 'w') as f:
        f.write(f"{m.numerator} {m.denominator}\n{D}\n1\n{len(Pex)}\n")
        for x, y in Pex:
            f.write(f"{int(x * D)} {int(y * D)} 1\n")
    return D


# ----------------------------------------------------------------------------------------------
# positioning polish: maximise the minimum margin for a fixed assignment (an LP), alternating
# with re-assignment
# ----------------------------------------------------------------------------------------------
def polish_positions(P, m, sq_list, rounds=6, delta=0.03, trust=0.1, move_pen=0.02, verbose=True):
    """P (n,2) float points; sq_list: Squares objects (families the points must hit).  Alternate:
    assign each square to its deepest point; LP: max t s.t. depth(p_a(Q), Q) >= t (4 linear
    constraints per square), 0 <= p <= m, |p - p0| <= trust; then a second LP with t fixed at t*
    maximising sum_Q min(slack_Q, delta) - move_pen * sum |p - p0|_1 to pull every non-forced
    constraint away from zero without wandering.  Returns (P_new, t*)."""
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix, csr_matrix
    P = np.asarray(P, float).copy(); n = len(P)
    tstar = None
    for it in range(rounds):
        P0 = P.copy()
        cons = []          # (point index, cx, cy, c, s)
        for sq in sq_list:
            D = sq.depth(P)                    # (n, N)
            a = D.argmax(0)
            for j in range(len(sq)):
                cons.append((int(a[j]), sq.cx[j], sq.cy[j], sq.c[j], sq.s[j]))
        N = len(cons)
        # variables: P (2n), t (1), slack u_Q (N), movement d (2n): |p - p0| <= d
        nv = 2 * n + 1 + N + 2 * n
        A = lil_matrix((4 * N + 4 * n, nv)); b = np.zeros(4 * N + 4 * n)
        for k, (i, cx, cy, c, s) in enumerate(cons):
            # x' = (px-cx) c + (py-cy) s ;  y' = -(px-cx) s + (py-cy) c ;  |x'|, |y'| <= 1/2 - t - u
            off_x = cx * c + cy * s; off_y = -cx * s + cy * c
            for r, (ax, ay, off) in enumerate(((c, s, off_x), (-s, c, off_y))):
                row = 4 * k + 2 * r
                A[row, 2 * i] = ax; A[row, 2 * i + 1] = ay; A[row, 2 * n] = 1; A[row, 2 * n + 1 + k] = 1
                b[row] = 0.5 + off
                A[row + 1, 2 * i] = -ax; A[row + 1, 2 * i + 1] = -ay; A[row + 1, 2 * n] = 1; A[row + 1, 2 * n + 1 + k] = 1
                b[row + 1] = 0.5 - off
        base = 4 * N; dof = 2 * n + 1 + N
        for v in range(2 * n):
            A[base + 2 * v, v] = 1; A[base + 2 * v, dof + v] = -1; b[base + 2 * v] = P0.ravel()[v]
            A[base + 2 * v + 1, v] = -1; A[base + 2 * v + 1, dof + v] = -1; b[base + 2 * v + 1] = -P0.ravel()[v]
        A = csr_matrix(A)
        lo = np.concatenate([np.maximum(P.ravel() - trust, 0.0), [-1.0], np.zeros(N), np.zeros(2 * n)])
        hi = np.concatenate([np.minimum(P.ravel() + trust, float(m)), [1.0], np.zeros(N), np.full(2 * n, trust)])
        c1 = np.zeros(nv); c1[2 * n] = -1.0
        r1 = linprog(c1, A_ub=A, b_ub=b, bounds=list(zip(lo, hi)), method='highs')
        if not r1.success:
            if verbose: print("  polish: stage-1 LP failed", r1.message); break
        t = -r1.fun
        lo2 = lo.copy(); hi2 = hi.copy()
        lo2[2 * n] = t - 1e-12; hi2[2 * n] = t - 1e-12
        hi2[2 * n + 1:dof] = delta
        c2 = np.zeros(nv); c2[2 * n + 1:dof] = -1.0; c2[dof:] = move_pen
        r2 = linprog(c2, A_ub=A, b_ub=b, bounds=list(zip(lo2, hi2)), method='highs')
        Pn = (r2.x if r2.success else r1.x)[:2 * n].reshape(n, 2)
        moved = np.abs(Pn - P).max()
        P = Pn
        if verbose:
            print(f"  polish {it}: t* = {t:+.3e}, moved {moved:.2e}, stage-2 obj = {-r2.fun if r2.success else float('nan'):.4f}", flush=True)
        if tstar is not None and abs(t - tstar) < 1e-12 and moved < 1e-9:
            tstar = t; break
        tstar = t
    return P, tstar


# ----------------------------------------------------------------------------------------------
# dumps for the independent re-solve
# ----------------------------------------------------------------------------------------------
def dump_instance(path_prefix, m, poses, B, rep, extra=None):
    write_family(path_prefix + '_family.txt', m, poses)
    np.save(path_prefix + '_B.npy', np.packbits(B, axis=1))
    np.save(path_prefix + '_rep.npy', rep)
    meta = dict(m=str(m), n_squares=len(poses), n_cand=int(B.shape[0]))
    if extra: meta.update(extra)
    with open(path_prefix + '_meta.json', 'w') as f:
        json.dump(meta, f, indent=1)


def load_instance(path_prefix):
    m, poses = read_family(path_prefix + '_family.txt')
    packed = np.load(path_prefix + '_B.npy')
    B = np.unpackbits(packed, axis=1)[:, :len(poses)].astype(bool)
    rep = np.load(path_prefix + '_rep.npy')
    return m, poses, B, rep


def drop_dominated_fast(B, verbose=False):
    """Dominance filter by posting lists: candidate i (row of the bool matrix B, K x n) is dropped
    iff some other row is a strict superset.  For each row, only rows sharing its rarest square are
    tested (a superset must contain that square), so the cost is sum_i |posting(rarest(i))|.
    Rows must be unique.  Returns the indices to keep (sorted)."""
    K, n = B.shape
    if K == 0: return np.zeros(0, dtype=np.int64)
    P = np.packbits(B, axis=1)
    cnt = B.sum(1).astype(np.int64)
    colcnt = B.sum(0)
    post = [np.nonzero(B[:, j])[0] for j in range(n)]
    keep = np.ones(K, dtype=bool)
    # rarest square of each row
    colrank = colcnt.astype(float) + 1e-9 * np.arange(n)
    for i in range(K):
        row = np.nonzero(B[i])[0]
        j = row[np.argmin(colrank[row])]
        cand = post[j]
        cand = cand[cnt[cand] > cnt[i]]
        if len(cand) == 0: continue
        sub = ((P[i][None, :] & ~P[cand]) == 0).all(1)
        if sub.any(): keep[i] = False
    if verbose:
        print(f"  dominance (fast): {K} -> {int(keep.sum())}", flush=True)
    return np.nonzero(keep)[0]


def solve_feasible(B, k, threads=4, time_limit=None, verbose=False, warm=None):
    """Feasibility IP: is there x in {0,1}^K with B^T x >= 1 and sum x <= k?  Returns
    dict(status, feasible (True/False/None), x, lp).  Infeasible <=> h(F) >= k + 1.  This is what
    the loop needs ("is there a 13-point hitting set?"), and it avoids proving that no 12-point
    set exists, which is where an optimisation IP spends most of its time."""
    import highspy
    K, n = B.shape
    starts = np.zeros(K + 1, dtype=np.int64); starts[1:] = np.cumsum(B.sum(1))
    index = np.concatenate([np.nonzero(B[kk])[0] for kk in range(K)]).astype(np.int32) if K else np.zeros(0, np.int32)
    value = np.ones(len(index))
    # add the cardinality row as row n
    starts2 = starts + np.arange(K + 1)
    index2 = np.zeros(len(index) + K, dtype=np.int32); value2 = np.ones(len(index) + K)
    for kk in range(K):
        index2[starts2[kk]:starts2[kk + 1] - 1] = index[starts[kk]:starts[kk + 1]]
        index2[starts2[kk + 1] - 1] = n
    h = highspy.Highs()
    h.setOptionValue('output_flag', verbose)
    h.setOptionValue('threads', int(threads))
    h.setOptionValue('mip_feasibility_tolerance', 1e-9)
    h.setOptionValue('primal_feasibility_tolerance', 1e-9)
    if time_limit: h.setOptionValue('time_limit', float(time_limit))
    lp = highspy.HighsLp()
    lp.num_col_ = K; lp.num_row_ = n + 1
    lp.col_cost_ = np.zeros(K); lp.col_lower_ = np.zeros(K); lp.col_upper_ = np.ones(K)
    lp.row_lower_ = np.concatenate([np.ones(n), [-highspy.kHighsInf]])
    lp.row_upper_ = np.concatenate([np.full(n, highspy.kHighsInf), [float(k)]])
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts2; lp.a_matrix_.index_ = index2; lp.a_matrix_.value_ = value2
    lp.integrality_ = [highspy.HighsVarType.kInteger] * K
    h.passModel(lp)
    if warm is not None:
        sol = highspy.HighsSolution(); sol.col_value = list(map(float, warm)); sol.value_valid = True
        h.setSolution(sol)
    h.run()
    st = h.getModelStatus()
    res = dict(status=h.modelStatusToString(st), feasible=None, x=None)
    if st == highspy.HighsModelStatus.kOptimal:
        res['feasible'] = True; res['x'] = np.array(h.getSolution().col_value)
    elif st == highspy.HighsModelStatus.kInfeasible:
        res['feasible'] = False
    return res


def lp_relaxation(B, threads=4):
    import highspy
    K, n = B.shape
    starts = np.zeros(K + 1, dtype=np.int64); starts[1:] = np.cumsum(B.sum(1))
    index = np.concatenate([np.nonzero(B[kk])[0] for kk in range(K)]).astype(np.int32)
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', int(threads))
    lp = highspy.HighsLp(); lp.num_col_ = K; lp.num_row_ = n
    lp.col_cost_ = np.ones(K); lp.col_lower_ = np.zeros(K); lp.col_upper_ = np.ones(K)
    lp.row_lower_ = np.ones(n); lp.row_upper_ = np.full(n, highspy.kHighsInf)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = np.ones(len(index))
    h.passModel(lp); h.run()
    return h.getInfo().objective_function_value
