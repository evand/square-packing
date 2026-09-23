#!/usr/bin/env python3
"""unavoid13-no: library for the lower-bound push "no 13-point pure unavoidable set for [0,4]^2"
(tasks/unavoid13-no).  Imports search/unavoid13_lib.py (never modified).

Columns ("candidate regions") instead of arrangement vertices.  A column is a closed region R of
the container, either a BOX [x0,x1] x [y0,y1] (a cell of a D4-symmetric quadtree over [0,m]^2) or a
POINT (an exact local cell, x0 = x1, y0 = y1).  Its incidence with a closed square Q is

    inc(R, Q) = 1   iff   R meets Q   (closed sets; for a point: the point lies in Q),

computed LENIENTLY in floats (superset of the truth) and exactly in Fractions for the certificate.

Soundness of the relaxation [proved, two lines].  The regions cover the container.  If P is any
point set hitting every square of F, map each p in P to a region R(p) containing it; then
inc(R(p), Q) >= [p in Q] for every Q, so the column set {R(p)} hits every row of F with <= |P|
columns.  Hence: if no k columns hit F, no k points hit F, i.e. h(F) >= k + 1.  A box is a
relaxation of every point inside it; refining a box (4 children, or the exact cells inside it)
keeps the cover-the-container property and can only make the IP harder to satisfy.

Exact cells of a box b [proved, unavoid13_lib docstring adapted].  For p in b let P(p) be the
intersection of the squares of F containing p; P(p) ∩ b is a compact convex polygon containing p
whose vertices are arrangement vertices of F inside b, intersections of square edges with the
boundary of b, or corners of b, and each such vertex lies in every square containing p.  So the
finite set of those points, with EXACT incidence, dominates every point of b.

D4 [proved].  F is D4-closed (the container is D4-invariant), the quadtree is D4-invariant by
construction (every split is applied to a whole D4-orbit of boxes) and exact cells of an orbit are
the images of one representative's cells; so g maps columns to columns and inc(gR, gQ) = inc(R, Q)
(exactly; in floats up to rounding far below the lenient tolerance).  Hence if a column set X hits
F, so does gX for every g.  ROOT CUT: let Q0 be a D4-invariant square of F (the tile-centred one at
(m/2, m/2), axis-parallel).  Every hitting set X contains a column meeting Q0; those columns form a
D4-invariant set C0; let R0 be a set of orbit representatives of C0.  Some g maps the column of X
meeting Q0 into R0, so gX satisfies  sum_{R in R0} x_R >= 1.  Adding that row preserves
feasibility, and infeasibility WITH the row implies infeasibility without it.
"""
import math, os, sys, time, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L

TOL = 1e-9          # lenient incidence tolerance (float)


# ----------------------------------------------------------------------------------------------
# columns
# ----------------------------------------------------------------------------------------------
class Columns:
    """Regions as float arrays x0, x1, y0, y1 (K,).  kind: 0 = box, 1 = point.
    Boxes carry an integer key (lvl, i, j) in a D4-invariant quadtree of [0,m]^2 with N0 boxes per
    side at level 0 (pitch m / (N0 2^lvl)); points carry key ('pt', tag) with a Fraction coordinate
    pair for exactness."""
    def __init__(self, m, N0):
        self.m = float(m); self.mF = F(m); self.N0 = N0
        self.x0 = np.zeros(0); self.x1 = np.zeros(0); self.y0 = np.zeros(0); self.y1 = np.zeros(0)
        self.kind = np.zeros(0, dtype=np.int8)
        self.keys = []          # list of keys
        self.exact = []         # Fraction coordinates (x0, x1, y0, y1)
        self.index = {}         # key -> column index

    def __len__(self): return len(self.keys)

    def _add(self, key, ex):
        x0, x1, y0, y1 = ex
        self.keys.append(key); self.exact.append(ex); self.index[key] = len(self.keys) - 1
        self.x0 = np.append(self.x0, float(x0)); self.x1 = np.append(self.x1, float(x1))
        self.y0 = np.append(self.y0, float(y0)); self.y1 = np.append(self.y1, float(y1))
        self.kind = np.append(self.kind, 0 if key[0] != 'pt' else 1)

    def add_many(self, keys, exs):
        for k, e in zip(keys, exs):
            if k not in self.index: self._add(k, e)

    def box_exact(self, lvl, i, j):
        n = self.N0 << lvl; d = self.mF / n
        return (i * d, (i + 1) * d, j * d, (j + 1) * d)

    def init_grid(self):
        for i in range(self.N0):
            for j in range(self.N0):
                self._add((0, i, j), self.box_exact(0, i, j))

    def rebuild(self, keep):
        """Keep the columns with indices in keep (bool mask or index array)."""
        keep = np.asarray(keep)
        if keep.dtype == bool: keep = np.nonzero(keep)[0]
        keys = [self.keys[k] for k in keep]; exs = [self.exact[k] for k in keep]
        self.x0 = np.zeros(0); self.x1 = np.zeros(0); self.y0 = np.zeros(0); self.y1 = np.zeros(0)
        self.kind = np.zeros(0, dtype=np.int8); self.keys = []; self.exact = []; self.index = {}
        # vectorised re-add
        self.keys = keys; self.exact = exs; self.index = {k: i for i, k in enumerate(keys)}
        self.x0 = np.array([float(e[0]) for e in exs]); self.x1 = np.array([float(e[1]) for e in exs])
        self.y0 = np.array([float(e[2]) for e in exs]); self.y1 = np.array([float(e[3]) for e in exs])
        self.kind = np.array([0 if k[0] != 'pt' else 1 for k in keys], dtype=np.int8)

    # ---- D4 on keys ----
    def d4_key(self, key, g):
        """g in 0..7: g = r + 4 f, r rotations by 90deg about the centre after an optional
        reflection x -> m - x."""
        if key[0] == 'pt':
            return ('pt', d4_point_exact(self.mF, key[1], g))
        lvl, i, j = key; n = self.N0 << lvl
        if g >= 4: i = n - 1 - i
        for _ in range(g % 4):
            i, j = n - 1 - j, i
        return (lvl, i, j)

    def d4_perm(self, g):
        """Column permutation for g (asserts closure)."""
        out = np.empty(len(self.keys), dtype=np.int64)
        for k, key in enumerate(self.keys):
            out[k] = self.index[self.d4_key(key, g)]
        return out

    def orbit_of(self, key):
        return sorted({self.d4_key(key, g) for g in range(8)})

    def reps(self):
        """Canonical representative (min key in orbit) for each column and a bool mask 'is rep'."""
        mask = np.zeros(len(self.keys), dtype=bool)
        for k, key in enumerate(self.keys):
            if key == min(self.orbit_of(key)): mask[k] = True
        return mask

    def split(self, keys):
        """Replace each box in keys (whole D4 orbits) by its 4 children.  Returns the new keys."""
        todo = set()
        for key in keys:
            if key[0] == 'pt': continue
            todo.update(self.orbit_of(key))
        drop = np.ones(len(self.keys), dtype=bool)
        newk = []; newe = []
        for key in todo:
            if key not in self.index: continue
            drop[self.index[key]] = False
            lvl, i, j = key
            for di in (0, 1):
                for dj in (0, 1):
                    kk = (lvl + 1, 2 * i + di, 2 * j + dj)
                    newk.append(kk); newe.append(self.box_exact(*kk))
        self.rebuild(drop)
        self.add_many(newk, newe)
        return newk

    def exactify(self, keys, sqF, tol=1e-7):
        """Replace each box in keys (whole D4 orbits) by its exact local cells (float coordinates,
        computed for one representative per orbit and mapped by D4 exactly as Fractions of the float
        values).  Returns the number of point columns added."""
        todo = {}
        for key in keys:
            if key[0] == 'pt': continue
            orb = self.orbit_of(key)
            todo[orb[0]] = orb
        drop = np.ones(len(self.keys), dtype=bool)
        newk = []; newe = []
        for rep, orb in todo.items():
            if rep not in self.index: continue
            for key in orb: drop[self.index[key]] = False
            pts = local_cells(self.mF, self.box_exact(*rep), sqF, tol=tol)
            # map to the orbit exactly
            for g in range(8):
                if self.d4_key(rep, g) not in orb: continue
                for (x, y) in pts:
                    ex = d4_point_exact(self.mF, (F(x), F(y)), g)
                    kk = ('pt', ex)
                    newk.append(kk); newe.append((ex[0], ex[0], ex[1], ex[1]))
        self.rebuild(drop)
        n0 = len(self.keys)
        self.add_many(newk, newe)
        return len(self.keys) - n0

    def rep_points(self, idx):
        """Representative point of each column (box centre / the point)."""
        return np.stack([(self.x0[idx] + self.x1[idx]) / 2, (self.y0[idx] + self.y1[idx]) / 2], 1)

    def save(self, path):
        with open(path, 'w') as f:
            f.write(f"# m = {self.mF}; N0 = {self.N0}; columns: 'box lvl i j' or 'pt x y' (Fractions)\n")
            for key in self.keys:
                if key[0] == 'pt': f.write(f"pt {key[1][0]} {key[1][1]}\n")
                else: f.write(f"box {key[0]} {key[1]} {key[2]}\n")

    @classmethod
    def load(cls, path):
        m = None; N0 = None; keys = []; exs = []
        for line in open(path):
            if line.startswith('#'):
                if 'm =' in line:
                    m = F(line.split('m =')[1].split(';')[0].strip()); N0 = int(line.split('N0 =')[1].split(';')[0].strip())
                continue
            t = line.split()
            if not t: continue
            if t[0] == 'box':
                keys.append((int(t[1]), int(t[2]), int(t[3])))
            else:
                keys.append(('pt', (F(t[1]), F(t[2]))))
        C = cls(m, N0)
        for key in keys:
            if key[0] == 'pt':
                exs.append((key[1][0], key[1][0], key[1][1], key[1][1]))
            else:
                exs.append(C.box_exact(*key))
        C.add_many(keys, exs)
        return C


def d4_point_exact(m, pt, g):
    x, y = pt
    if g >= 4: x = m - x
    for _ in range(g % 4):
        x, y = m - y, x
    return (x, y)


def d4_pose_perm(m, poses):
    """(8, n) permutations of the family under D4 (asserts closure)."""
    idx = {p: i for i, p in enumerate(poses)}
    out = np.empty((8, len(poses)), dtype=np.int64)
    for g in range(8):
        for i, (u, cx, cy) in enumerate(poses):
            uu = u; x, y = cx, cy
            if g >= 4: x = m - x; uu = L.norm_u(-uu)
            for _ in range(g % 4):
                x, y = m - y, x
            out[g, i] = idx[L.pose_key(uu, x, y)]
    return out


# ----------------------------------------------------------------------------------------------
# incidence
# ----------------------------------------------------------------------------------------------
def incidence(cols, sq, idx=None, tol=TOL):
    """Lenient float incidence (K, n) bool of the columns (boxes or points) with the squares sq
    (L.Squares).  SAT test on 4 axes for boxes; closed containment with tolerance."""
    if idx is None: idx = np.arange(len(cols))
    x0 = cols.x0[idx]; x1 = cols.x1[idx]; y0 = cols.y0[idx]; y1 = cols.y1[idx]
    mx = (x0 + x1) / 2; my = (y0 + y1) / 2; hx = (x1 - x0) / 2; hy = (y1 - y0) / 2
    c = sq.c; s = sq.s; cx = sq.cx; cy = sq.cy
    ac = np.abs(c); asn = np.abs(s); w2 = (ac + asn) / 2
    out = np.ones((len(idx), len(sq)), dtype=bool)
    chunk = max(1, 20000000 // max(1, len(sq)))
    for i0 in range(0, len(idx), chunk):
        sl = slice(i0, i0 + chunk)
        MX = mx[sl, None]; MY = my[sl, None]; HX = hx[sl, None]; HY = hy[sl, None]
        ok = np.abs(MX - cx[None, :]) <= w2[None, :] + HX + tol
        ok &= np.abs(MY - cy[None, :]) <= w2[None, :] + HY + tol
        p1 = MX * c[None, :] + MY * s[None, :] - (cx * c + cy * s)[None, :]
        ok &= np.abs(p1) <= 0.5 + HX * ac[None, :] + HY * asn[None, :] + tol
        p2 = -MX * s[None, :] + MY * c[None, :] - (-cx * s + cy * c)[None, :]
        ok &= np.abs(p2) <= 0.5 + HX * asn[None, :] + HY * ac[None, :] + tol
        out[sl] = ok
    return out


def incidence_exact_one(ex, pose):
    """Exact (Fraction) test: does the closed box/point ex = (x0, x1, y0, y1) meet the closed unit
    square pose = (u, cx, cy)?"""
    x0, x1, y0, y1 = ex
    u, cx, cy = pose
    c, s = L.trig(u)
    hx = (x1 - x0) / 2; hy = (y1 - y0) / 2; mx = (x0 + x1) / 2; my = (y0 + y1) / 2
    w2 = (abs(c) + abs(s)) / 2
    if abs(mx - cx) > w2 + hx: return False
    if abs(my - cy) > w2 + hy: return False
    if abs(mx * c + my * s - (cx * c + cy * s)) > F(1, 2) + hx * abs(c) + hy * abs(s): return False
    if abs(-mx * s + my * c - (-cx * s + cy * c)) > F(1, 2) + hx * abs(s) + hy * abs(c): return False
    return True


def incidence_exact(cols, poses, sq=None, verbose=True):
    """Exact incidence (K, n): float decides pairs far (> 1e-6 in every SAT slack) from the
    boundary, Fractions decide the rest."""
    t0 = time.time()
    if sq is None: sq = L.Squares(poses)
    K = len(cols); n = len(poses)
    Bl = incidence(cols, sq, tol=1e-6)         # lenient: superset
    Bs = incidence(cols, sq, tol=-1e-6)        # strict: subset
    B = Bs.copy()
    amb = np.nonzero(Bl & ~Bs)
    for k, j in zip(*amb):
        B[k, j] = incidence_exact_one(cols.exact[k], poses[j])
    if verbose:
        print(f"  exact incidence: {K} x {n}, {len(amb[0])} borderline pairs decided in Fractions, "
              f"{int(B.sum())} incidences, {time.time()-t0:.1f}s", flush=True)
    return B


# ----------------------------------------------------------------------------------------------
# exact local cells of a box (float coordinates, lenient incidence for the loop)
# ----------------------------------------------------------------------------------------------
def local_cells(m, box, sq, tol=1e-7):
    """Points dominating every point of the closed box: arrangement vertices of the squares meeting
    the box that lie in it, intersections of square edges with the box boundary, the 4 corners.
    Returns a list of (x, y) floats after dedupe of lenient incidence sets and dominance."""
    x0, x1, y0, y1 = map(float, box)
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2; hx, hy = (x1 - x0) / 2, (y1 - y0) / 2
    # squares meeting the box (lenient)
    cols1 = _OneBox(x0, x1, y0, y1)
    inc = incidence(cols1, sq, tol=1e-7)[0]
    js = np.nonzero(inc)[0]
    pts = [np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])]
    if len(js):
        sub = L.Squares([sq.poses[j] for j in js])
        V = L.arrangement_vertices(sub)
        inb = (V[:, 0] >= x0 - 1e-9) & (V[:, 0] <= x1 + 1e-9) & (V[:, 1] >= y0 - 1e-9) & (V[:, 1] <= y1 + 1e-9)
        pts.append(V[inb])
        # edge x box-boundary intersections
        C = sub.corners
        for a in range(4):
            P = C[:, a]; R = C[:, (a + 1) % 4] - P
            for (X, axis) in ((x0, 0), (x1, 0), (y0, 1), (y1, 1)):
                den = R[:, axis]
                ok = np.abs(den) > 1e-12
                t = np.where(ok, (X - P[:, axis]) / np.where(ok, den, 1.0), -1.0)
                ok &= (t >= -1e-9) & (t <= 1 + 1e-9)
                if ok.any():
                    Q = P[ok] + t[ok, None] * R[ok]
                    o = 1 - axis
                    lo, hi = (y0, y1) if axis == 0 else (x0, x1)
                    ok2 = (Q[:, o] >= lo - 1e-9) & (Q[:, o] <= hi + 1e-9)
                    pts.append(Q[ok2])
    V = np.concatenate(pts, 0)
    V[:, 0] = np.clip(V[:, 0], x0, x1); V[:, 1] = np.clip(V[:, 1], y0, y1)
    # lenient incidence with ALL squares, dedupe, dominance
    packed = L.incidence_packed(sq, V, tol=tol)
    uniq, idx = L.dedupe_rows(packed)
    B = L.unpack_rows(uniq, len(sq))
    keep = L.drop_dominated_fast(B) if len(B) > 1 else np.arange(len(B))
    return [tuple(V[idx[k]]) for k in keep]


class _OneBox:
    def __init__(self, x0, x1, y0, y1):
        self.x0 = np.array([x0]); self.x1 = np.array([x1]); self.y0 = np.array([y0]); self.y1 = np.array([y1])
    def __len__(self): return 1


# ----------------------------------------------------------------------------------------------
# the IP
# ----------------------------------------------------------------------------------------------
def _colwise(B, extra_rows=()):
    """B (K, n) bool -> (starts, index, value) column-wise with extra dense-ish rows appended:
    extra_rows: list of (row_coeffs (K,) float) appended after the n cover rows."""
    K, n = B.shape
    Bi = B.astype(np.int8)
    cnt = B.sum(1).astype(np.int64)
    ex = [np.asarray(r, float) for r in extra_rows]
    excnt = np.zeros(K, dtype=np.int64)
    for r in ex: excnt += (r != 0)
    starts = np.zeros(K + 1, dtype=np.int64); starts[1:] = np.cumsum(cnt + excnt)
    index = np.zeros(starts[-1], dtype=np.int32); value = np.zeros(starts[-1])
    # cover part: nonzero positions row by row
    rows_k, cols_j = np.nonzero(B)
    # positions: for candidate k, its j's are contiguous in np.nonzero order (row-major)
    pos = starts[rows_k] + (np.arange(len(rows_k)) - np.searchsorted(rows_k, rows_k, side='left'))
    index[pos] = cols_j; value[pos] = 1.0
    for ri, r in enumerate(ex):
        ks = np.nonzero(r)[0]
        # place after the cover entries and the previous extra rows
        off = starts[ks] + cnt[ks]
        for prev in ex[:ri]:
            off = off + (prev[ks] != 0)
        index[off] = n + ri; value[off] = r[ks]
    return starts, index, value


def solve_ip(B, k, threads=4, time_limit=None, verbose=False, warm=None, extra_rows=(),
             feasibility=True, seed=None, heuristic_effort=None, symmetry=None, presolve=None):
    """Feasibility IP: x in {0,1}^K, B^T x >= 1, sum x <= k, extra rows (coeffs, lo, hi).
    feasibility=True: zero objective; else minimise sum x with objective cutoff k + 0.5.
    Returns dict(status, feasible (True/False/None), x, nodes, time)."""
    import highspy
    K, n = B.shape
    rows = [r[0] for r in extra_rows]
    starts, index, value = _colwise(B, rows)
    # cardinality row
    starts2 = starts + np.arange(K + 1)
    index2 = np.zeros(len(index) + K, dtype=np.int32); value2 = np.ones(len(index) + K)
    nrow0 = n + len(rows)
    for kk in range(K):
        index2[starts2[kk]:starts2[kk + 1] - 1] = index[starts[kk]:starts[kk + 1]]
        value2[starts2[kk]:starts2[kk + 1] - 1] = value[starts[kk]:starts[kk + 1]]
        index2[starts2[kk + 1] - 1] = nrow0
    h = highspy.Highs()
    h.setOptionValue('output_flag', verbose)
    h.setOptionValue('threads', int(threads))
    h.setOptionValue('mip_feasibility_tolerance', 1e-9)
    h.setOptionValue('primal_feasibility_tolerance', 1e-9)
    if seed is not None: h.setOptionValue('random_seed', int(seed))
    if heuristic_effort is not None: h.setOptionValue('mip_heuristic_effort', float(heuristic_effort))
    if symmetry is not None: h.setOptionValue('mip_detect_symmetry', bool(symmetry))
    if presolve is not None: h.setOptionValue('presolve', presolve)
    if time_limit: h.setOptionValue('time_limit', float(time_limit))
    lp = highspy.HighsLp()
    lp.num_col_ = K; lp.num_row_ = nrow0 + 1
    lp.col_cost_ = np.zeros(K) if feasibility else np.ones(K)
    lp.col_lower_ = np.zeros(K); lp.col_upper_ = np.ones(K)
    rl = [np.ones(n)] + [[r[1]] for r in extra_rows] + [[-highspy.kHighsInf]]
    ru = [np.full(n, highspy.kHighsInf)] + [[r[2]] for r in extra_rows] + [[float(k)]]
    lp.row_lower_ = np.concatenate([np.asarray(v, float) for v in rl])
    lp.row_upper_ = np.concatenate([np.asarray(v, float) for v in ru])
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts2; lp.a_matrix_.index_ = index2; lp.a_matrix_.value_ = value2
    lp.integrality_ = [highspy.HighsVarType.kInteger] * K
    h.passModel(lp)
    if not feasibility:
        h.setOptionValue('objective_bound', float(k) + 0.5)
    if warm is not None:
        sol = highspy.HighsSolution(); sol.col_value = list(map(float, warm)); sol.value_valid = True
        h.setSolution(sol)
    t0 = time.time()
    h.run()
    st = h.getModelStatus()
    info = h.getInfo()
    res = dict(status=h.modelStatusToString(st), feasible=None, x=None, time=time.time() - t0,
               nodes=int(info.mip_node_count), bound=float(info.mip_dual_bound))
    if st == highspy.HighsModelStatus.kOptimal:
        res['feasible'] = True; res['x'] = np.array(h.getSolution().col_value)
    elif st in (highspy.HighsModelStatus.kInfeasible,):
        res['feasible'] = False
    elif not feasibility and st == highspy.HighsModelStatus.kObjectiveBound:
        res['feasible'] = False
    return res


def lp_value(B, threads=4, extra_rows=()):
    import highspy
    K, n = B.shape
    rows = [r[0] for r in extra_rows]
    starts, index, value = _colwise(B, rows)
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('threads', int(threads))
    lp = highspy.HighsLp(); lp.num_col_ = K; lp.num_row_ = n + len(rows)
    lp.col_cost_ = np.ones(K); lp.col_lower_ = np.zeros(K); lp.col_upper_ = np.ones(K)
    lp.row_lower_ = np.concatenate([np.ones(n), [r[1] for r in extra_rows]])
    lp.row_upper_ = np.concatenate([np.full(n, highspy.kHighsInf), [r[2] for r in extra_rows]])
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts; lp.a_matrix_.index_ = index; lp.a_matrix_.value_ = value
    h.passModel(lp); h.run()
    return h.getInfo().objective_function_value, np.array(h.getSolution().row_dual)


def root_cut(cols, B, j0):
    """Symmetry row: sum over D4-orbit representatives of the columns meeting square j0 >= 1.
    Returns (coeffs (K,), 1.0, inf).  Requires B[:, j0] to be a D4-invariant column set (checked)."""
    K = B.shape[0]
    meet = np.nonzero(B[:, j0])[0]
    meetset = set(int(k) for k in meet)
    coef = np.zeros(K)
    for k in meet:
        key = cols.keys[k]
        orb = cols.orbit_of(key)
        for o in orb:
            assert cols.index[o] in meetset, "columns meeting the invariant square are not D4-closed"
        if key == orb[0]: coef[k] = 1.0
    return (coef, 1.0, float('inf'))


def reduce_columns(B):
    """Dedupe identical incidence rows and drop dominated ones (rows that are strict subsets).
    Returns kept indices (sorted).  Sound for hitting sets: a dropped column's incidence is
    contained in a kept column's."""
    packed = np.packbits(B, axis=1)
    _, idx = np.unique(packed, axis=0, return_index=True)
    idx = np.sort(idx)
    uniq = B[idx]
    nz = uniq.any(1); uniq = uniq[nz]; idx = idx[nz]
    keep = L.drop_dominated_fast(uniq)
    return idx[keep]


# ----------------------------------------------------------------------------------------------
# families
# ----------------------------------------------------------------------------------------------
def support_poses(m, path='search/cover4_exact_support.txt'):
    sup = L.read_support(path, m)
    return L.d4_closure(m, [p for p, _ in sup])


def family_union(m, lists):
    seen = set(); out = []
    for lst in lists:
        for p in lst:
            if p not in seen:
                seen.add(p); out.append(p)
    return out


def hits_all(P, sq, tol=1e-9):
    D = sq.depth(np.asarray(P, float)).max(0)
    return D, np.nonzero(D < -tol)[0]


# ----------------------------------------------------------------------------------------------
# exact arrangement cells via the C helper (float, lenient, maximal cells only)
# ----------------------------------------------------------------------------------------------
_CELLS = None


def _cells_lib():
    global _CELLS
    if _CELLS is None:
        import ctypes, subprocess
        here = os.path.dirname(os.path.abspath(__file__))
        so = os.path.join(here, 'unavoid13no_cells.so'); src = os.path.join(here, 'unavoid13no_cells.c')
        if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
            subprocess.check_call(['gcc', '-O2', '-fopenmp', '-shared', '-fPIC', '-o', so, src, '-lm'])
        lib = ctypes.CDLL(so)
        lib.maximal_cells.restype = ctypes.c_long
        lib.maximal_cells.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                      ctypes.c_double, ctypes.c_double, ctypes.c_int,
                                      ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),
                                      ctypes.POINTER(ctypes.POINTER(ctypes.c_uint64)), ctypes.POINTER(ctypes.c_long)]
        lib.free_buf.argtypes = [ctypes.c_void_p]
        _CELLS = (lib, ctypes)
    return _CELLS


def maximal_cells(m, poses, tol=1e-11, threads=4, dominance=True, verbose=True):
    """Maximal cells of the arrangement of the family (float, lenient incidence with tol).
    Returns (sq, V (K,2) representative vertices, B (K,n) bool incidence), deduped, and (if
    dominance) passed through the exact posting-list dominance filter as well."""
    lib, ctypes = _cells_lib()
    t0 = time.time()
    sq = L.Squares(poses)
    n = len(sq)
    cx = np.ascontiguousarray(sq.cx); cy = np.ascontiguousarray(sq.cy)
    cc = np.ascontiguousarray(sq.c); ss = np.ascontiguousarray(sq.s)
    pc = ctypes.POINTER(ctypes.c_double); pb = ctypes.POINTER(ctypes.c_uint64)
    coords = pc(); bits = pb(); nvert = ctypes.c_long(0)
    total = lib.maximal_cells(n, cx.ctypes.data_as(pc), cy.ctypes.data_as(pc), cc.ctypes.data_as(pc),
                              ss.ctypes.data_as(pc), float(m), float(tol), int(threads),
                              ctypes.byref(coords), ctypes.byref(bits), ctypes.byref(nvert))
    NW = (n + 63) // 64
    V = np.ctypeslib.as_array(coords, shape=(max(total, 1), 2)).copy()[:total]
    W = np.ctypeslib.as_array(bits, shape=(max(total, 1), NW)).copy()[:total]
    lib.free_buf(ctypes.cast(coords, ctypes.c_void_p)); lib.free_buf(ctypes.cast(bits, ctypes.c_void_p))
    packed = np.ascontiguousarray(W).view(np.uint8).reshape(total, NW * 8)
    _, idx = np.unique(packed, axis=0, return_index=True)
    idx = np.sort(idx)
    B = np.unpackbits(packed[idx], axis=1, bitorder='little')[:, :n].astype(bool)
    V = V[idx]
    nz = B.any(1); B = B[nz]; V = V[nz]
    k1 = len(B)
    if dominance:
        keep = L.drop_dominated_fast(B); B = B[keep]; V = V[keep]
    if verbose:
        print(f"  cells: {n} squares, {nvert.value} vertices streamed, {total} survivors, {k1} distinct, "
              f"{len(B)} after dominance, {time.time()-t0:.1f}s", flush=True)
    return sq, V, B


def cell_d4_perms(m, poses, B):
    """(8, K) permutations of the cells under D4, matched through incidence sets: the image of a
    cell with incidence S is the cell with incidence g.S.  Asserts closure (float asymmetry would
    show up as a missing row)."""
    perm = d4_pose_perm(m, poses)
    K, n = B.shape
    packed = np.packbits(B, axis=1)
    index = {packed[k].tobytes(): k for k in range(K)}
    out = np.empty((8, K), dtype=np.int64)
    for g in range(8):
        Bg = np.zeros_like(B)
        Bg[:, perm[g]] = B           # column i of B goes to column perm[g][i]
        pg = np.packbits(Bg, axis=1)
        for k in range(K):
            key = pg[k].tobytes()
            assert key in index, f"cell {k} has no D4 image under g={g}"
            out[g, k] = index[key]
    return out


def root_cut_cells(m, V, B, j0, poses=None, perms=None):
    """Symmetry row over D4-orbit representatives of the cells meeting square j0."""
    if perms is None:
        perms = cell_d4_perms(m, poses, B)
    K = B.shape[0]
    meet = np.nonzero(B[:, j0])[0]
    coef = np.zeros(K)
    for k in meet:
        orb = perms[:, k]
        if k == orb.min(): coef[k] = 1.0
    return (coef, 1.0, float('inf'))


def d4_close_cells(m, poses, V, B, verbose=True):
    """Make the cell set exactly D4-closed by ADDING every missing image (sound: more candidates).
    Returns (V, B, perms (8, K))."""
    perm = d4_pose_perm(m, poses)
    mF = F(m)
    K, n = B.shape
    packed = np.packbits(B, axis=1)
    index = {packed[k].tobytes(): k for k in range(K)}
    addV = []; addB = []
    rows = list(B); vs = list(V)
    for g in range(1, 8):
        for k in range(len(rows)):
            row = rows[k]
            img = np.zeros(n, dtype=bool); img[perm[g]] = row
            key = np.packbits(img).tobytes()
            if key not in index:
                x, y = d4_point_exact(mF, (F(float(vs[k][0])), F(float(vs[k][1]))), g)
                index[key] = len(rows); rows.append(img); vs.append((float(x), float(y)))
    if len(rows) > K:
        B = np.array(rows); V = np.array(vs)
    if verbose:
        print(f"  D4-closing cells: {K} -> {len(rows)} (+{len(rows)-K} images added)", flush=True)
    perms = cell_d4_perms(m, poses, B)
    return V, B, perms


# ----------------------------------------------------------------------------------------------
# local search for k-cell hitting sets (heuristic; a found set is a true hitting set of the rows)
# ----------------------------------------------------------------------------------------------
class Cover:
    def __init__(self, B):
        self.B = B; self.K, self.n = B.shape
        self.rows_of = [np.nonzero(B[c])[0] for c in range(self.K)]
        self.cells_of = [np.nonzero(B[:, r])[0] for r in range(self.n)]


def sa_search(cov, k, start=None, iters=200000, seed=0, T0=1.5, T1=0.05, sample=40, tabu_len=15,
              time_limit=None, verbose=False):
    """Simulated annealing on the number of unhit rows.  Returns (best_cells, best_cost)."""
    rng = np.random.default_rng(seed)
    K, n = cov.K, cov.n
    if start is None or len(start) != k:
        chosen = list(rng.choice(K, k, replace=False))
    else:
        chosen = [int(c) for c in start]
    cnt = np.zeros(n, dtype=np.int32)
    for c in chosen: cnt[cov.rows_of[c]] += 1
    cost = int((cnt == 0).sum())
    best = (list(chosen), cost)
    tabu = {}
    t0 = time.time()
    for it in range(iters):
        if cost == 0: break
        if time_limit and (it & 1023) == 0 and time.time() - t0 > time_limit: break
        T = T0 * (T1 / T0) ** (it / max(1, iters))
        unhit = np.nonzero(cnt == 0)[0]
        r = int(unhit[rng.integers(len(unhit))])
        cand = cov.cells_of[r]
        if len(cand) > sample: cand = cand[rng.choice(len(cand), sample, replace=False)]
        # gain of each candidate: rows with cnt == 0 among its rows
        gains = np.array([int((cnt[cov.rows_of[c]] == 0).sum()) for c in cand])
        # loss of each chosen: rows with cnt == 1 among its rows
        losses = np.array([int((cnt[cov.rows_of[c]] == 1).sum()) for c in chosen])
        # best pair (approximate: ignore overlap between removed and added)
        ci = int(gains.argmax()); c_add = int(cand[ci])
        if c_add in tabu and tabu[c_add] > it and rng.random() < 0.9:
            # pick another
            order = np.argsort(-gains)
            for ci in order:
                if int(cand[ci]) not in tabu or tabu[int(cand[ci])] <= it: c_add = int(cand[ci]); break
        li = int(losses.argmin())
        # randomise the removal a little
        if rng.random() < 0.3: li = int(rng.integers(k))
        c_rem = chosen[li]
        if c_rem == c_add: continue
        # exact delta
        cnt[cov.rows_of[c_rem]] -= 1
        cnt[cov.rows_of[c_add]] += 1
        new_cost = int((cnt == 0).sum())
        delta = new_cost - cost
        if delta <= 0 or rng.random() < math.exp(-delta / T):
            chosen[li] = c_add; cost = new_cost
            tabu[c_rem] = it + tabu_len
            if cost < best[1]: best = (list(chosen), cost)
        else:
            cnt[cov.rows_of[c_add]] -= 1
            cnt[cov.rows_of[c_rem]] += 1
        if verbose and it % 20000 == 0:
            print(f"    sa it {it}: cost {cost}, best {best[1]}, T {T:.3f}", flush=True)
    return best


def multi_sa(B, k, starts=(), restarts=4, iters=200000, seed=0, time_limit=None, verbose=False):
    """Try SA from each start and from random starts; return the first exact hitting set found
    (index array) or None."""
    cov = Cover(B)
    t0 = time.time()
    trials = list(starts) + [None] * restarts
    for i, st in enumerate(trials):
        if time_limit and time.time() - t0 > time_limit: break
        cells, cost = sa_search(cov, k, start=st, iters=iters, seed=seed + i, verbose=verbose,
                                time_limit=None if not time_limit else max(1.0, time_limit - (time.time() - t0)))
        if verbose: print(f"    sa trial {i}: cost {cost}, {time.time()-t0:.0f}s", flush=True)
        if cost == 0:
            return np.array(cells)
    return None


# ----------------------------------------------------------------------------------------------
# portfolio finder: several HiGHS configurations in parallel processes, first verdict wins
# ----------------------------------------------------------------------------------------------
def _race_worker(packed, n, k, extra_rows, cfg, time_limit, q):
    try:
        B = np.unpackbits(packed, axis=1)[:, :n].astype(bool)
        extra_rows = [(np.asarray(c), lo, hi) for (c, lo, hi) in extra_rows]
        r = solve_ip(B, k, threads=1, time_limit=time_limit, extra_rows=extra_rows, **cfg)
        q.put((cfg, r['status'], r['feasible'], None if r['x'] is None else np.nonzero(r['x'] > 0.5)[0].tolist(), r['time'], r['nodes']))
    except Exception as e:
        q.put((cfg, f'worker error: {e!r}', None, None, 0.0, 0))


def solve_ip_race(B, k, extra_rows=(), nproc=4, time_limit=3600, seeds=None):
    """Run nproc HiGHS solves (alternating feasibility / optimisation form, different seeds) as
    separate processes; return the first verdict (feasible with a solution, or infeasible) and kill
    the rest.  A verdict from any single run is a verdict for the instance (same model)."""
    import multiprocessing as mp
    # 'spawn', not 'fork': a process forked from a parent that already ran HiGHS (its thread pool)
    # dies silently in HiGHS
    ctx = mp.get_context('spawn')
    q = ctx.Queue()
    seeds = seeds or list(range(1, nproc + 1))
    cfgs = [dict(feasibility=(i % 2 == 0), seed=int(seeds[i])) for i in range(nproc)]
    packed = np.packbits(B, axis=1); n = B.shape[1]
    ex = [(np.asarray(c, float), float(lo), float(hi)) for (c, lo, hi) in extra_rows]
    procs = [ctx.Process(target=_race_worker, args=(packed, n, k, ex, cfg, time_limit, q)) for cfg in cfgs]
    for p in procs: p.start()
    t0 = time.time(); result = None
    while time.time() - t0 < time_limit + 60:
        try:
            cfg, status, feas, sol, dt, nodes = q.get(timeout=5)
        except Exception:
            if not any(p.is_alive() for p in procs): break
            continue
        if feas is True or feas is False:
            result = dict(cfg=cfg, status=status, feasible=feas, sol=sol, time=dt, nodes=nodes); break
    for p in procs:
        if p.is_alive(): p.terminate()
    for p in procs: p.join(timeout=5)
    return result
