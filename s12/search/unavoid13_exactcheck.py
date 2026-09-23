#!/usr/bin/env python3
"""Independent exact checker for the unavoid13 hitting-set infeasibility certificates.

    python3 search/unavoid13_exactcheck.py FAMILY.txt CERT.txt [--threads 4]

Prints `VERIFIED: ...` or `NOT VERIFIED: ...` (exit 0 either way); malformed input -> `ERROR:` and
a non-zero exit.  Shares no code with the other unavoid13* modules: everything below is written
from the definitions.

CLAIM `plain`, bound k.  No set of at most k points meets every closed unit square of FAMILY.
CLAIM `C2`,    bound k.  No set of at most k points that is invariant under the half-turn
                         s(p) = (m, m) - p meets every closed unit square of FAMILY.

Objects.  A pose `u cx cy` (rationals) is the closed unit square with centre (cx, cy) and edge
directions e1 = (c, s), e2 = (-s, c), c = (1 - u^2)/(1 + u^2), s = 2u/(1 + u^2).  p is in it iff
|<p - ctr, e1>| <= 1/2 and |<p - ctr, e2>| <= 1/2.  Admissible in [0, m]^2 iff both centre
coordinates lie in [w/2, m - w/2], w = |c| + |s|  (checked exactly for every square).

What is checked, all in exact rational / integer arithmetic:

 (i)   Vertices.  E = the family (plain) or the family closed under the half-turn (C2; the image of
       pose (u, cx, cy) is (u, m - cx, m - cy)).  V = all corners of squares of E, plus every
       intersection point of two non-parallel closed edges of two different squares of E (pairs
       whose centres are more than sqrt 2 apart cannot meet and are skipped; the test is the exact
       |ctr_i - ctr_j|^2 > 2), plus (C2) the centre (m/2, m/2).
 (ii)  Incidence of every vertex and every certificate candidate in every square of E, exact.
       By default every (point, square) pair is decided in Fractions (after the exact test
       |p - ctr| > 3/4 in a coordinate => outside).  `--float-filter` (faster,
       same verdict) lets a float pre-filter decide a pair when both frame coordinates are more than
       1e-9 away from +-1/2: every input is a correctly rounded double (relative error 2^-53) of
       magnitude <= m + 1 <= 5, and <p - ctr, e> is 3 roundings + 2 products + 1 sum, so its float
       error is < 64 * 5 * 2^-53 < 4e-14; the remaining pairs are decided in Fractions.
 (iii) Rows = the squares of FAMILY in file order (0-based).  Column of a point p: plain,
       {i : p in Q_i}; C2, {i : p in Q_i or s(p) in Q_i}, with cost 1 if p is the centre, else 2
       (the size of its orbit).  Candidates = the points listed in the certificate.
       Domination: for every vertex v, some candidate's column contains v's column (C2: with cost
       <= cost(v)).  Soundness (notes/unavoid13.md §1, §4): for any point p the polygon
       P(p) = intersection of the squares of E containing p has a vertex in V lying in all of
       them; so any hitting set (C2: any symmetric one, orbit by orbit, the centre kept) maps to
       candidate columns of no larger total cost that still cover every row.
 (iv)  The tree.  Preorder; `B j` branches on candidate j (which must be free at that node):
       first subtree x_j = 1, then subtree x_j = 0.  The two children partition the node's 0/1
       box, so by induction the leaves partition {0,1}^N.  A leaf is closed by
         `U i`   row i has no candidate that is not fixed to 0 (no x in the leaf covers it), or
         `D den i:num ...`  y_i = num/den >= 0 (rows not listed: y_i = 0).  For every x in the
                 leaf's box with A x >= 1:  c.x >= L := sum_i y_i + sum_{j fixed 1} r_j +
                 sum_{j free} min(0, r_j),  r = c - A^T y  (weak LP duality, valid for any y >= 0).
                 Every integer x in the leaf has c.x in {C1 + g t : t >= 0 integer}, C1 = cost of
                 the fixed-1 candidates, g = gcd of the free candidates' costs (no free candidate:
                 c.x = C1).  The leaf is closed iff the least such value >= L exceeds k.
                 (plain: g = 1, i.e. the rounding is `ceil(L) > k`, i.e. L > k; C2 with the centre
                 fixed: g = 2, so L > k - 1 suffices when C1 and k have different parity.)
       Every leaf closed and every subtree present  =>  no 0/1 x with A x >= 1 has c.x <= k.
"""
import sys, os, time, math, argparse, hashlib
from fractions import Fraction as Fr
from multiprocessing import Pool
import numpy as np

HALF = Fr(1, 2)
EPS = 1e-9
Q34 = Fr(3, 4)


class Bad(Exception):
    pass


# --------------------------------------------------------------------------------------------
# geometry (exact)
# --------------------------------------------------------------------------------------------
def cos_sin(u):
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d


def frame(pose):
    u, cx, cy = pose
    c, s = cos_sin(u)
    return cx, cy, c, s


def corners_of(pose):
    cx, cy, c, s = frame(pose)
    out = []
    for a, b in ((-HALF, -HALF), (HALF, -HALF), (HALF, HALF), (-HALF, HALF)):
        out.append((cx + a * c - b * s, cy + a * s + b * c))
    return out


def inside_exact(fr, p):
    cx, cy, c, s = fr
    dx = p[0] - cx; dy = p[1] - cy
    x = dx * c + dy * s; y = -dx * s + dy * c
    return -HALF <= x <= HALF and -HALF <= y <= HALF


def admissible(pose, m):
    _, cx, cy, c, s = (None,) + frame(pose)
    h = (abs(c) + abs(s)) / 2
    return h <= cx <= m - h and h <= cy <= m - h


def edge_intersections(args):
    """All intersection points of non-parallel closed edges of square pairs (i, j) in the chunk."""
    pairs, corners = args
    out = set()
    for i, j in pairs:
        Ci = corners[i]; Cj = corners[j]
        for a in range(4):
            P = Ci[a]; P2 = Ci[(a + 1) % 4]
            R0 = P2[0] - P[0]; R1 = P2[1] - P[1]
            for b in range(4):
                Q = Cj[b]; Q2 = Cj[(b + 1) % 4]
                S0 = Q2[0] - Q[0]; S1 = Q2[1] - Q[1]
                den = R0 * S1 - R1 * S0
                if den == 0:
                    continue
                W0 = Q[0] - P[0]; W1 = Q[1] - P[1]
                tn = W0 * S1 - W1 * S0
                sn = W0 * R1 - W1 * R0
                if den > 0:
                    if tn < 0 or tn > den or sn < 0 or sn > den: continue
                else:
                    if tn > 0 or tn < den or sn > 0 or sn < den: continue
                t = tn / den
                out.add((P[0] + t * R0, P[1] + t * R1))
    return out


# incidence worker state (set before the pool forks)
_G = {}


def incidence_chunk(pts):
    """pts: list of exact points -> list of Python-int bitmasks over the squares of E."""
    fr = _G['frames']; ffl = _G['ffl']; pure = _G['pure']
    cx, cy, c, s = ffl
    out = []
    nexact = 0
    for p in pts:
        mask = 0
        if pure:
            # exact box pre-test: every point of a unit square lies within sqrt(2)/2 < 3/4 of its
            # centre in each coordinate, so |p - ctr| > 3/4 in a coordinate means "outside"
            px, py = p
            for j in range(len(fr)):
                f = fr[j]
                if abs(px - f[0]) > Q34 or abs(py - f[1]) > Q34: continue
                if inside_exact(f, p): mask |= 1 << j
            out.append(mask); continue
        px = float(p[0]); py = float(p[1])
        dx = px - cx; dy = py - cy
        x = np.abs(dx * c + dy * s); y = np.abs(-dx * s + dy * c)
        sure_in = (x < 0.5 - EPS) & (y < 0.5 - EPS)
        sure_out = (x > 0.5 + EPS) | (y > 0.5 + EPS)
        for j in np.nonzero(sure_in)[0]:
            mask |= 1 << int(j)
        for j in np.nonzero(~sure_in & ~sure_out)[0]:
            nexact += 1
            if inside_exact(fr[int(j)], p): mask |= 1 << int(j)
        out.append(mask)
    return out, nexact


def incidence(points, threads, chunk=500):
    jobs = [points[i:i + chunk] for i in range(0, len(points), chunk)]
    if threads > 1 and len(jobs) > 1:
        with Pool(threads) as pool:
            res = pool.map(incidence_chunk, jobs)
    else:
        res = [incidence_chunk(j) for j in jobs]
    masks = []; nexact = 0
    for r, ne in res:
        masks.extend(r); nexact += ne
    return masks, nexact


# --------------------------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------------------------
def parse_family(path):
    m = None; poses = []
    for ln, line in enumerate(open(path), 1):
        t = line.split()
        if not t: continue
        if t[0] == '#':
            if len(t) >= 4 and t[1] == 'm' and t[2] == '=' and m is None:
                m = Fr(t[3].rstrip(';'))
            continue
        if t[0].startswith('#'): continue
        if t[0] != 'pose' or len(t) != 4:
            raise Bad(f"family line {ln}: expected 'pose u cx cy'")
        try:
            poses.append(tuple(Fr(x) for x in t[1:]))
        except (ValueError, ZeroDivisionError):
            raise Bad(f"family line {ln}: not rationals")
    if m is None: raise Bad("family: no '# m = ...' header")
    if not poses: raise Bad("family: no squares")
    return m, poses


def parse_cert(path):
    hdr = {}; cands = []; tree = []
    section = 'hdr'
    for ln, line in enumerate(open(path), 1):
        t = line.split()
        if not t or t[0].startswith('#'): continue
        try:
            if section == 'hdr':
                if t[0] == 'cand':
                    cands.append((Fr(t[1]), Fr(t[2])))
                elif t[0] == 'tree':
                    section = 'tree'
                else:
                    hdr[t[0]] = t[1:]
            elif section == 'tree':
                if t[0] == 'end':
                    section = 'end'
                elif t[0] == 'B' and len(t) == 2:
                    tree.append(('B', int(t[1])))
                elif t[0] == 'U' and len(t) == 2:
                    tree.append(('U', int(t[1])))
                elif t[0] == 'D' and len(t) >= 2:
                    den = int(t[1])
                    ys = []
                    for tok in t[2:]:
                        i, num = tok.split(':')
                        ys.append((int(i), int(num)))
                    tree.append(('D', den, ys))
                else:
                    raise Bad(f"cert line {ln}: bad tree record")
            else:
                raise Bad(f"cert line {ln}: content after 'end'")
        except (ValueError, IndexError, ZeroDivisionError):
            raise Bad(f"cert line {ln}: malformed")
    if section != 'end': raise Bad("cert: missing 'tree' section or 'end' record (truncated?)")
    for key in ('claim', 'm', 'k'):
        if key not in hdr: raise Bad(f"cert: missing '{key}'")
    claim = hdr['claim'][0]
    if claim not in ('plain', 'C2'): raise Bad(f"cert: unknown claim {claim}")
    return claim, Fr(hdr['m'][0]), int(hdr['k'][0]), hdr, cands, tree


# --------------------------------------------------------------------------------------------
# the tree
# --------------------------------------------------------------------------------------------
def check_tree(tree, A, cost, k, log):
    """A: (N, R) int64 0/1 matrix (candidates x rows); cost: list of ints.  Iterative preorder."""
    N, R = A.shape
    colsets = [np.nonzero(A[:, i])[0] for i in range(R)]     # candidates covering row i
    state = np.zeros(N, dtype=np.int8)                       # 0 free, 1 fixed to 1, -1 fixed to 0
    stats = dict(branch=0, U=0, D=0, maxdepth=0, minslack=None)
    pos = 0
    # stack of pending actions: ('node',) visit next record; ('set0', j) switch j to 0 then visit;
    # ('undo', j)
    stack = [('node',)]
    depth = 0
    costs = np.array(cost, dtype=np.int64)
    while stack:
        act = stack.pop()
        if act[0] == 'undo':
            state[act[1]] = 0; depth -= 1; continue
        if act[0] == 'set0':
            state[act[1]] = -1; stack.append(('node',)); continue
        if pos >= len(tree): raise Bad("tree ends early: a subtree is missing")
        rec = tree[pos]; pos += 1
        if rec[0] == 'B':
            j = rec[1]
            if not (0 <= j < N): raise Bad(f"tree record {pos}: candidate {j} out of range")
            if state[j] != 0: raise Bad(f"tree record {pos}: candidate {j} already fixed")
            stats['branch'] += 1
            depth += 1; stats['maxdepth'] = max(stats['maxdepth'], depth)
            state[j] = 1
            # after the x_j = 1 subtree: set x_j = 0, do that subtree, then undo
            stack.append(('undo', j)); stack.append(('set0', j)); stack.append(('node',))
        elif rec[0] == 'U':
            i = rec[1]
            if not (0 <= i < R): raise Bad(f"tree record {pos}: row {i} out of range")
            if (state[colsets[i]] != -1).any():
                return False, f"leaf {pos} (U {i}): row {i} is still coverable", stats
            stats['U'] += 1
        else:
            _, den, ys = rec
            if den <= 0: return False, f"leaf {pos}: denominator {den} <= 0", stats
            y = np.zeros(R, dtype=np.int64)
            tot = 0
            for i, num in ys:
                if not (0 <= i < R): raise Bad(f"tree record {pos}: row {i} out of range")
                if num < 0: return False, f"leaf {pos}: negative dual y_{i}", stats
                if y[i] != 0: raise Bad(f"tree record {pos}: row {i} listed twice")
                if num >= 2 ** 40: return False, f"leaf {pos}: numerator too large", stats
                y[i] = num; tot += num
            if den >= 2 ** 40 or tot >= 2 ** 60:
                return False, f"leaf {pos}: numbers too large for the int64 path", stats
            # rD = den * r = den * c - A y, exact in int64: every partial sum of A y is <= tot < 2^60
            # and |den * c_j| <= 2 * 2^40 * max cost; the sums over candidates are Python ints
            rD = costs * den - A @ y
            fixed1 = state == 1; free = state == 0
            LD = (tot + sum(rD[fixed1].tolist())
                  + sum(min(v, 0) for v in rD[free].tolist()))                     # den * L
            L = Fr(LD, den)
            C1 = int(costs[fixed1].sum())
            if free.any():
                g = 0
                for cj in set(costs[free].tolist()): g = math.gcd(g, int(cj))
                t = max(0, math.ceil((L - C1) / g))
                least = C1 + g * t
            else:
                least = C1
                if L > C1:        # no free candidate and the LP bound exceeds the only value:
                    least = math.inf  # the box is a single point, infeasible for A x >= 1
            if not least > k:
                return False, (f"leaf {pos}: bound L = {float(L):.6f} (C1 = {C1}) does not "
                               f"exclude cost <= {k}"), stats
            sl = float(L) - k
            stats['minslack'] = sl if stats['minslack'] is None else min(stats['minslack'], sl)
            stats['D'] += 1
    if pos != len(tree): raise Bad(f"tree has {len(tree) - pos} records after the root is closed")
    return True, "", stats


# --------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family'); ap.add_argument('cert')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--float-filter', action='store_true',
                    help='let a float pre-filter with a proven 1e-9 margin decide clear incidences')
    a = ap.parse_args()
    t0 = time.time()
    try:
        m, poses = parse_family(a.family)
        claim, mc, k, hdr, cands, tree = parse_cert(a.cert)
    except (Bad, OSError) as e:
        print(f"ERROR: {e}"); sys.exit(2)
    sha = hashlib.sha256(open(a.family, 'rb').read()).hexdigest()
    print(f"family {a.family}: m = {m}, {len(poses)} squares, sha256 {sha[:16]}...")
    print(f"certificate {a.cert}: claim {claim}, k = {k}, {len(cands)} candidates, {len(tree)} tree records")
    if 'family_sha256' in hdr and hdr['family_sha256'][0] != sha:
        print("  note: the family's sha256 differs from the one recorded in the certificate")
    if mc != m:
        print(f"NOT VERIFIED: certificate m = {mc} but family m = {m}"); return
    R = len(poses)
    if len(set(poses)) != R:
        print("  note: the family has repeated poses (harmless: repeated rows)")
    for i, p in enumerate(poses):
        if not admissible(p, m):
            print(f"NOT VERIFIED: square {i} (pose {' '.join(map(str, p))}) is not admissible in [0,{m}]^2"); return
    if not (0 < k): print("NOT VERIFIED: k must be positive"); return
    # E: the family (+ half-turn images for C2); rows are E[:R]
    E = list(poses)
    if claim == 'C2':
        have = set(E)
        for (u, cx, cy) in poses:
            q = (u, m - cx, m - cy)
            if q not in have: have.add(q); E.append(q)
        sig = {p: i for i, p in enumerate(E)}
        partner = [sig[(u, m - cx, m - cy)] for (u, cx, cy) in E]
    print(f"  E: {len(E)} squares ({len(E) - R} half-turn images added)")
    # (i) vertices
    corners = [corners_of(p) for p in E]
    V = set()
    for C in corners: V.update(C)
    pairs = []
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            dx = E[i][1] - E[j][1]; dy = E[i][2] - E[j][2]
            if dx * dx + dy * dy <= 2: pairs.append((i, j))
    nch = max(1, a.threads * 8)
    jobs = [(pairs[i::nch], corners) for i in range(nch)]
    if a.threads > 1:
        with Pool(a.threads) as pool:
            for s in pool.map(edge_intersections, jobs): V |= s
    else:
        for jb in jobs: V |= edge_intersections(jb)
    centre = (m / 2, m / 2)
    if claim == 'C2': V.add(centre)
    V = sorted(V)
    print(f"  (i) {len(V)} exact arrangement vertices ({len(pairs)} close pairs), {time.time()-t0:.1f}s")
    # (ii) incidence
    _G['frames'] = [frame(p) for p in E]
    _G['ffl'] = tuple(np.array([float(fr[q]) for fr in _G['frames']]) for q in range(4))
    _G['pure'] = not a.float_filter
    vm, ne1 = incidence(V, a.threads)
    cm, ne2 = incidence(cands, a.threads)
    print(f"  (ii) incidence: {len(V) + len(cands)} points x {len(E)} squares, " +
          (f"{ne1 + ne2} near-boundary pairs decided in Fractions" if a.float_filter else "all pairs decided in Fractions")
          + f", {time.time()-t0:.1f}s")
    rowmask = (1 << R) - 1

    def column(mask):
        if claim == 'plain':
            return mask & rowmask
        img = 0
        mm = mask
        while mm:
            low = mm & -mm; j = low.bit_length() - 1; mm ^= low
            img |= 1 << partner[j]
        return (mask | img) & rowmask     # p in Q_i or p in s(Q_i)  <=>  orbit meets Q_i

    def pcost(p):
        return 1 if (claim == 'C2' and p == centre) else (2 if claim == 'C2' else 1)

    ccol = [column(x) for x in cm]
    ccost = [pcost(p) for p in cands]
    # (iii) domination
    vcols = {}
    for p, x in zip(V, vm):
        col = column(x)
        if col == 0: continue
        c = pcost(p)
        if vcols.get(col, 99) > c: vcols[col] = c
    post = [[] for _ in range(R)]
    for idx, col in enumerate(ccol):
        mm = col
        while mm:
            low = mm & -mm; post[low.bit_length() - 1].append(idx); mm ^= low
    postlen = [len(x) for x in post]
    bad = 0; example = None
    for col, c in vcols.items():
        # a dominating candidate must cover col's rarest row
        mm = col; best = None
        while mm:
            low = mm & -mm; i = low.bit_length() - 1; mm ^= low
            if best is None or postlen[i] < postlen[best]: best = i
        if not any((ccol[j] & col) == col and ccost[j] <= c for j in post[best]):
            bad += 1
            if example is None: example = (col, c)
    print(f"  (iii) {len(vcols)} distinct vertex columns; each dominated by a candidate: "
          f"{'yes' if bad == 0 else 'NO (' + str(bad) + ' not dominated)'}, {time.time()-t0:.1f}s")
    if bad:
        rows = [i for i in range(R) if example[0] >> i & 1]
        print(f"NOT VERIFIED: {bad} vertex columns are not dominated by any candidate "
              f"(e.g. rows {rows[:12]}{'...' if len(rows) > 12 else ''}, cost {example[1]})")
        return
    # (iv) tree
    A = np.zeros((len(cands), R), dtype=np.int64)
    for j, col in enumerate(ccol):
        mm = col
        while mm:
            low = mm & -mm; A[j, low.bit_length() - 1] = 1; mm ^= low
    try:
        ok, why, st = check_tree(tree, A, ccost, k, None)
    except Bad as e:
        print(f"NOT VERIFIED: malformed tree: {e}"); return
    if not ok:
        print(f"NOT VERIFIED: {why}"); return
    nleaf = st['U'] + st['D']
    print(f"  (iv) tree: {st['branch']} branchings, {nleaf} leaves (U {st['U']}, D {st['D']}), "
          f"depth {st['maxdepth']}, min leaf margin L - k = {st['minslack']:.3g}, {time.time()-t0:.1f}s")
    what = (f"no set of at most {k} points meets every one of the {R} closed unit squares"
            if claim == 'plain' else
            f"no half-turn-symmetric set of at most {k} points meets every one of the {R} closed unit squares")
    print(f"VERIFIED: {what} in {a.family} (container [0,{m}]^2), {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
