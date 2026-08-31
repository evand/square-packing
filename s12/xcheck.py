#!/usr/bin/env python3
"""INDEPENDENT exact re-check of a certificate, written from scratch in Python.

Deliberately a separate implementation from verify/ (different language, different code path,
re-derived from the definition in certificates/FORMAT.md rather than from the Rust), so that
agreement between the two is meaningful evidence.  Everything is exact: `fractions.Fraction`
for the derivation of each angle bin, Python big integers (all quantities of a bin scaled by one
common denominator) for the sweep itself.  No floating point is used for any decision.

WHAT IS CHECKED
    For every CLOSED unit square Q contained in the container C = [0,s]^2, at every centre and
    every angle, the total weight of the certificate points lying in Q is >= 1.  With `--n n`
    the total weight is also required to be < n, which is what makes the file a proof of
    s(n) >= s (see FORMAT.md).

HOW  (all steps exact)
  1. Angles.  theta_k = 2*arctan(k/N), k = 0, 1, ..., so cos theta_k = (N^2-k^2)/(N^2+k^2) and
     sin theta_k = 2kN/(N^2+k^2) are rational.  Bin k is the closed angle interval
     [theta_k, theta_{k+1}].  Bins are taken to cover [0, 90 deg) (k = 0..N-1), or only
     [0, 45 deg] when the point multiset is invariant under the reflection (x,y) -> (y,x)
     (checked exactly, together with the two axis reflections, i.e. full D4 invariance);
     that reflection preserves the container and the points and maps a square at angle theta
     to one at angle 90 deg - theta, so [45, 90) follows from [0, 45].
  2. Shrink.  Let delta = theta_{k+1} - theta_k (its cos and sin are rational by the addition
     formulas; the code asserts 0 <= delta <= 45 deg).  A square of side sigma at angle
     theta_k, concentric with a unit square at angle theta, has bounding-box half-width
     (sigma/2)(cos d + sin d) in the unit square's frame, where d = theta - theta_k.  Since
     cos d + sin d is increasing on [0, 45 deg], sigma_k = 1/(cos delta + sin delta) makes the
     sigma_k-square lie inside the unit square for EVERY theta in the bin.  Hence it suffices
     that every admissible placement of the closed sigma_k-square at angle theta_k covers >= 1.
     sigma_k is used exactly (not rounded).
  3. Centres.  A unit square at angle theta lies in C iff its centre lies in [w/2, s-w/2]^2 with
     w = cos theta + sin theta.  Over a bin, w is minimised at an endpoint (it is sqrt2 *
     sin(theta+45deg), unimodal on [0, 90deg]), so the centre box B_k = [w_min/2, s-w_min/2]^2
     with w_min = min(w(theta_k), w(theta_{k+1})) contains every admissible centre of the bin.
  4. Sweep.  Rotate everything by -theta_k, so the sigma_k-square is axis-aligned: point p is
     covered by the square centred at u iff |q(p) - u|_inf <= h, h = sigma_k/2, where q is the
     rotated position.  The covered weight f(u) is constant on each open cell of the arrangement
     of the lines u0 = q0 +- h, u1 = q1 +- h, and f is upper semicontinuous (closed squares), so
         min over the rotated box P = q(B_k) of f  =  min over open cells meeting P of f.
     (Every point of P is a limit of interior points of P, which lie in open cells meeting P.)
     The code enumerates x-strips between consecutive breakpoints, computes the exact y-range of
     the convex polygon P over each strip, and within it the y-cells of the atoms active in x;
     a cell reaching beyond all active atoms means some admissible centre covers weight 0.
     Cell values come from prefix sums over the atoms sorted by rotated y.

MODES
    python3 xcheck.py CERT N                 exhaustive: every bin k
    python3 xcheck.py CERT N --all           same
    python3 xcheck.py CERT N STRIDE          fast/sampled: only bins k = 0, STRIDE, 2*STRIDE, ...
                                             (STRIDE = 1 is exhaustive)
    options: --n n   also require total weight < n
             -v      print the exact minimum of every bin checked
             -j J    worker processes (default: all CPUs)
The program always prints the overall minimum covered weight as an exact fraction and a final
line VERIFIED or NOT VERIFIED.

BRANCH CERTIFICATES (FORMAT.md, "Branch certificates").  A trailer `region corner r_num r_den /
lambda L / k K` after the points changes the claim to: squares whose centre lies in one of the
four closed corner boxes [0,r]^2, [s-r,s]x[0,r], ... capture >= 1 + L/W, all others >= 1, and
total - (L/W)*K < n.  Here every cell of the sweep is classified EXACTLY (rational arithmetic,
no padding, unlike the verifier's float-with-padding test): the cell is a convex quadrilateral in
container coordinates, so it lies inside a box iff its bounding box does and can meet a box only
if its bounding box does; a cell that may meet a box must reach 1 + lambda, a cell that may leave
the boxes must reach 1, and the reported minimum is min over cells of (covered - required extra).
The condition (2r-1)^2 < 2 (one square per box) is checked exactly.  A sampled run can only say NOT VERIFIED with certainty; its
VERIFIED is qualified as "(sampled bins only)".
"""
from fractions import Fraction as F
from bisect import bisect_left, bisect_right
from math import lcm, degrees, atan
import sys, os, time, argparse, multiprocessing

# ----------------------------------------------------------------------------------------------
def load(path):
    t = open(path).read().split()
    sn, sd, D, WD, m = (int(v) for v in t[:5])
    assert sn > 0 and sd > 0 and D > 0 and WD > 0 and m >= 0
    A = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    assert all(w >= 0 for _, _, w in A), "weights must be non-negative"
    rest = t[5 + 3*m:]
    region = None; cliques = None
    if rest and rest[0] == "cliques":
        # `cliques N Q c`, then per clique `w b` and b box lines `k U0LO U0HI U1LO U1HI`
        N_cl, Q, c = int(rest[1]), int(rest[2]), int(rest[3])
        assert N_cl >= 3 and Q > 0 and c >= 0, "clique block: N >= 3, Q > 0, c >= 0 required"
        pos = 4; ws = []; boxes = []
        for ci in range(c):
            w, nb = int(rest[pos]), int(rest[pos + 1]); pos += 2
            assert w >= 0, "clique weights must be non-negative"
            assert nb >= 1, "a clique needs at least one box"
            ws.append(w)
            for _ in range(nb):
                k, lo0, hi0, lo1, hi1 = (int(v) for v in rest[pos:pos + 5]); pos += 5
                assert 0 <= k < N_cl, "box bin k must satisfy 0 <= k < N"
                assert lo0 <= hi0 and lo1 <= hi1, "box rectangle has LO > HI"
                boxes.append((ci, k, lo0, hi0, lo1, hi1))
        cliques = dict(N=N_cl, Q=Q, w=ws, boxes=boxes)
        rest = rest[pos:]
    if rest:
        # `region corner r_num r_den / lambda L [L2 L3 L4] / k K [K2 K3 K4]`
        assert len(rest) >= 8 and rest[0] == "region" and rest[1] == "corner" and rest[4] == "lambda" and "k" in rest[5:], \
            "trailing data is not a `region corner r_num r_den / lambda ... / k ...` trailer"
        ik = rest.index("k", 5)
        lams = [int(v) for v in rest[5:ik]]; ks = [int(v) for v in rest[ik+1:]]
        r = F(int(rest[2]), int(rest[3]))
        assert r > 0 and (2*r - 1)**2 < 2, "region trailer: r > 0 and (2r-1)^2 < 2 required (one square per corner box)"
        if len(lams) == 1 and len(ks) == 1:
            assert 0 <= ks[0] <= 4, "region trailer: 0 <= k <= 4 required"
            region = dict(r=r, lam=[lams[0]]*4, k=[ks[0]]*4, kdot=lams[0]*ks[0], single=True, ktot=ks[0])
        elif len(lams) == 4 and len(ks) == 4:
            assert all(v in (0, 1) for v in ks), "region trailer: per-box k must be 0 or 1"
            region = dict(r=r, lam=lams, k=ks, kdot=sum(l*k for l, k in zip(lams, ks)), single=False, ktot=sum(ks))
        else:
            raise AssertionError("region trailer: 1 or 4 lambda values and as many k values expected")
    return F(sn, sd), D, WD, A, region, cliques

def box_core(box, Q, N):
    """EXACT core of a box (FORMAT.md, 'Clique certificates'): in the frame of its bin k, every
    unit square with a pose in the box contains the axis-parallel rectangle
    [hi0 - h, lo0 + h] x [hi1 - h, lo1 + h], h = sigma_k / 2 with the exact (unrounded) sigma_k.
    Returns ((lo0, hi0, lo1, hi1) as Fractions, (cos, sin) of theta_k)."""
    _, k, lo0, hi0, lo1, hi1 = box
    c, sn, sigma, _ = bin_params(k, N)
    h = sigma / 2
    core = (F(hi0, Q) - h, F(lo0, Q) + h, F(hi1, Q) - h, F(lo1, Q) + h)
    return core, (c, sn)

def cores_meet(ca, ra, cb, rb):
    """do two cores (closed rectangles, each axis-parallel in its own bin's frame) intersect?
    Exact separating-axis test: two convex polygons are disjoint iff an edge normal separates
    them; here the four edge normals are the axes of the two frames."""
    def project_onto(core, rot, rot_axes):
        # corners of `core` (frame `rot`) expressed in the frame `rot_axes`: u' = R(theta - theta') u
        (c, s), (c2, s2) = rot, rot_axes
        cc = c * c2 + s * s2; ss = s * c2 - c * s2          # cos, sin of (theta - theta')
        p0 = []; p1 = []
        for u0 in (core[0], core[1]):
            for u1 in (core[2], core[3]):
                p0.append(cc * u0 - ss * u1); p1.append(ss * u0 + cc * u1)
        return (min(p0), max(p0)), (min(p1), max(p1))
    for (x, rx), (y, ry) in (((ca, ra), (cb, rb)), ((cb, rb), (ca, ra))):
        (q0lo, q0hi), (q1lo, q1hi) = project_onto(x, rx, ry)
        if q0hi < y[0] or y[1] < q0lo or q1hi < y[2] or y[3] < q1lo:
            return False
    return True

def check_cliques(cl, N):
    """refuse (AssertionError) any clique not certified by pairwise core intersection"""
    assert cl['N'] == N, f"clique block is defined on the net N={cl['N']} but the check runs with N={N}"
    cores = [box_core(b, cl['Q'], N) for b in cl['boxes']]
    for i, (core, _) in enumerate(cores):
        assert core[0] <= core[1] and core[2] <= core[3], f"box {i+1} of clique {cl['boxes'][i][0]+1} has an empty core: not a clique"
    by_cl = {}
    for i, b in enumerate(cl['boxes']): by_cl.setdefault(b[0], []).append(i)
    for cid, idx in by_cl.items():
        for a in range(len(idx)):
            for b in range(a, len(idx)):
                i, j = idx[a], idx[b]
                assert cores_meet(cores[i][0], cores[i][1], cores[j][0], cores[j][1]), \
                    f"clique {cid+1}: cores of boxes {a+1} and {b+1} do not meet: not certified to be a clique"

def d4_symmetric(A, sD):
    """exact: is the weighted point multiset invariant under the dihedral group of [0,s]^2 ?"""
    S = sorted(A)
    gens = [lambda x, y: (sD - x, y), lambda x, y: (x, sD - y), lambda x, y: (y, x)]
    return all(sorted((*g(x, y), w) for x, y, w in A) == S for g in gens)

def rot(k, N):
    """(cos, sin) of theta_k = 2*arctan(k/N), exactly"""
    g = N*N + k*k
    return F(N*N - k*k, g), F(2*k*N, g)

def bin_params(k, N):
    """(cos theta_k, sin theta_k, sigma_k, w_min) for bin k = [theta_k, theta_{k+1}]"""
    c0, s0 = rot(k, N); c1, s1 = rot(k+1, N)
    cd = c0*c1 + s0*s1          # cos(theta_{k+1} - theta_k)
    sd = c0*s1 - s0*c1          # sin(theta_{k+1} - theta_k)
    assert sd >= 0 and cd >= sd, "bin wider than 45 deg: N too small"
    sigma = 1 / (cd + sd)
    wmin = min(c0 + s0, c1 + s1)
    return c0, s0, sigma, wmin

def yrange_over_strip(V, a, b):
    """exact [ymin, ymax] of the convex polygon V (integer vertices) over the strip a<=x<=b"""
    ys = []
    n = len(V)
    for i in range(n):
        (px, py), (qx, qy) = V[i], V[(i+1) % n]
        if a <= px <= b:
            ys.append(F(py))
        if px != qx:
            for xc in (a, b):
                if min(px, qx) <= xc <= max(px, qx):
                    ys.append(py + F((qy - py) * (xc - px), qx - px))
    assert ys
    return min(ys), max(ys)

# state shared with worker processes (set once in main, inherited by fork)
G = {}

def min_cover_bin(k):
    """exact minimum, over every admissible centre, of the weight covered by the closed
    sigma_k-square at angle theta_k, minus the region's required extra (0 without a trailer).
    Returns (k, min_numerator or None, witness)."""
    N, S, D, A, region = G['N'], G['s'], G['D'], G['A'], G['region']
    WD = G['WD']; cliques = G['cliques']
    c, sn, sigma, wmin = bin_params(k, N)
    # clique boxes of this bin: (clique id, lo0, hi0, lo1, hi1) over Q, in the frame of bin k
    cboxes = [(b[0], b[2], b[3], b[4], b[5]) for b in cliques['boxes'] if b[1] == k] if cliques else []
    Qc = cliques['Q'] if cliques else 1; cw = cliques['w'] if cliques else []
    def clique_credit(u0a, u0b, u1a, u1b):
        """weight of the cliques one of whose boxes contains the closed u-space cell (Fractions), each once"""
        if not cboxes: return 0
        credit = 0; seen = set()
        for cid, lo0, hi0, lo1, hi1 in cboxes:
            if cid in seen: continue
            if F(lo0, Qc) <= u0a and u0b <= F(hi0, Qc) and F(lo1, Qc) <= u1a and u1b <= F(hi1, Qc):
                credit += cw[cid]; seen.add(cid)
        return credit
    h = sigma / 2
    L = wmin / 2; U = S - wmin / 2
    if U < L:
        return k, None, None        # no unit square at these angles fits in the container
    # rotate by -theta_k: atoms and the centre box
    Q = [(c*F(x, D) + sn*F(y, D), -sn*F(x, D) + c*F(y, D), w) for x, y, w in A]
    poly = [(c*a + sn*b, -sn*a + c*b) for a, b in ((L, L), (U, L), (U, U), (L, U))]
    # corner boxes (container coordinates, closed) and their u-space images
    boxes = []; bpolys = []; lams = [0, 0, 0, 0]
    if region is not None:
        r = region['r']; lams = region['lam']
        for x0, y0 in ((0, 0), (S - r, 0), (0, S - r), (S - r, S - r)):
            boxes.append((x0, x0 + r, y0, y0 + r))
            bpolys.append([(c*a + sn*b, -sn*a + c*b) for a, b in ((x0, y0), (x0 + r, y0), (x0 + r, y0 + r), (x0, y0 + r))])
    def cell_flags(u0a, u0b, u1a, u1b):
        """(list of boxes the cell may meet, index of a box it certainly lies inside or None) for
        the u-space cell [u0a,u0b]x[u1a,u1b] (Fractions, container units), EXACTLY: the cell's image
        is the convex hull of its four rotated corners, and a convex set lies inside / meets an
        axis-parallel box iff its bounding box does / does."""
        xs = []; ys = []
        for u0, u1 in ((u0a, u1a), (u0b, u1a), (u0a, u1b), (u0b, u1b)):
            xs.append(c*u0 - sn*u1); ys.append(sn*u0 + c*u1)
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        may = [j for j, (bx0, bx1, by0, by1) in enumerate(boxes) if x1 >= bx0 and x0 <= bx1 and y1 >= by0 and y0 <= by1]
        inside = next((j for j, (bx0, bx1, by0, by1) in enumerate(boxes) if x0 >= bx0 and x1 <= bx1 and y0 >= by0 and y1 <= by1), None)
        return may, inside
    def required(may, inside):
        # no contact -> 0; entirely inside box j -> lam_j; straddling -> max over the boxes met of max(lam_j, 0)
        if not may: return 0
        if inside is not None: return lams[inside]
        return max(max(lams[j], 0) for j in may)
    if L == U:                      # a single admissible centre: evaluate f there directly
        u0, u1 = poly[0]
        tot = sum(w for q0, q1, w in Q if abs(q0 - u0) <= h and abs(q1 - u1) <= h) + clique_credit(u0, u0, u1, u1)
        may, inside = cell_flags(u0, u0, u1, u1)
        return k, tot - required(may, inside), (u0, u1)
    # one common denominator for the whole bin -> integer sweep
    den = 1
    for v in [h] + [q[0] for q in Q] + [q[1] for q in Q] + [p[0] for p in poly] + [p[0] for bp in bpolys for p in bp] + [p[1] for bp in bpolys for p in bp]:
        den = lcm(den, v.denominator)
    for v in [p[1] for p in poly]:
        den = lcm(den, v.denominator)
    def I(v):
        v = v * den; assert v.denominator == 1; return v.numerator
    H = I(h)
    atoms = sorted((I(q0), I(q1), w) for q0, q1, w in Q)
    qx = [t[0] for t in atoms]
    V = [(I(px), I(py)) for px, py in poly]
    BV = [[(I(px), I(py)) for px, py in bp] for bp in bpolys]
    px0 = min(p[0] for p in V); px1 = max(p[0] for p in V)
    bx = sorted(set([x - H for x in qx] + [x + H for x in qx] + [px0, px1]))
    best = None; wit = None
    for a, b in zip(bx, bx[1:]):
        if b <= px0 or a >= px1:
            continue                                  # strip does not meet P
        # atoms whose x-window [q0-h, q0+h] contains the whole strip (CLOSED squares)
        i0 = bisect_left(qx, b - H); i1 = bisect_right(qx, a + H)
        ylo, yhi = yrange_over_strip(V, a, b)
        assert ylo < yhi
        # u1-bands where the strip meets a corner box (integer units); cells outside all bands
        # cannot meet the region
        bands = []
        for bv in BV:
            bpx0 = min(p[0] for p in bv); bpx1 = max(p[0] for p in bv)
            if b <= bpx0 or a >= bpx1: continue
            lo, hi = yrange_over_strip(bv, max(a, bpx0), min(b, bpx1))
            bands.append((lo, hi))
        def flags(cl, dl):
            """flags of the cell [a,b] x [max(cl,ylo), min(dl,yhi)] (a superset of its admissible part)"""
            cl = max(cl, ylo); dl = min(dl, yhi)
            if not any(dl >= lo and cl <= hi for lo, hi in bands): return [], None
            return cell_flags(F(a, den), F(b, den), F(cl, den), F(dl, den))
        if i1 <= i0:
            v = clique_credit(F(a, den), F(b, den), F(ylo, den), F(yhi, den)) - required(*flags(ylo, yhi))
            if v < WD:
                return k, v, (F(a + b, 2*den), (ylo + yhi) / (2*den))   # nothing covered in this strip
            if best is None or v < best: best = v; wit = (F(a + b, 2*den), (ylo + yhi) / (2*den))
            continue
        ys = sorted((atoms[i][1], atoms[i][2]) for i in range(i0, i1))
        yv = [t[0] for t in ys]
        pre = [0]
        for _, w in ys: pre.append(pre[-1] + w)
        by = sorted(set([y - H for y in yv] + [y + H for y in yv]))
        # y-cells (by[j], by[j+1]), j = -1..len(by)-1 with by[-1] = -inf, by[len] = +inf;
        # the cell meets P iff by[j+1] > ylo and by[j] < yhi
        j_start = bisect_right(by, ylo) - 1
        j_end = bisect_left(by, yhi) - 1
        for j in range(j_start, j_end + 1):
            if j < 0 or j >= len(by) - 1:
                cl = ylo if j < 0 else by[j]; dl = by[0] if j < 0 else yhi
                v = clique_credit(F(a, den), F(b, den), F(cl, den), F(dl, den)) - required(*flags(cl, dl))
                if v < WD:
                    return k, v, (F(a + b, 2*den), (ylo + yhi) / (2*den))   # beyond every atom's window
                if best is None or v < best: best = v; wit = (F(a + b, 2*den), (ylo + yhi) / (2*den))
                continue
            cl, dl = by[j], by[j+1]
            k0 = bisect_left(yv, dl - H); k1 = bisect_right(yv, cl + H)
            tot = (pre[k1] - pre[k0] if k1 > k0 else 0) + clique_credit(F(a, den), F(b, den), F(cl, den), F(dl, den))
            v = tot - required(*flags(cl, dl))
            if best is None or v < best:
                best = v
                wit = (F(a + b, 2*den), F(max(cl, ylo) + min(dl, yhi), 2*den))
    assert best is not None
    return k, best, wit

def witness_xy(k, N, wit):
    """rotate a witness centre back to container coordinates (for display only)"""
    c, sn = rot(k, N)
    u0, u1 = wit
    return float(c*u0 - sn*u1), float(sn*u0 + c*u1)

# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="independent exact re-check of a certificate")
    ap.add_argument("cert"); ap.add_argument("N", type=int)
    ap.add_argument("stride", type=int, nargs="?", default=1, help="check bins k = 0, stride, 2*stride, ... (1 = all)")
    ap.add_argument("--all", action="store_true", help="exhaustive (same as stride 1)")
    ap.add_argument("--n", type=int, default=None, help="also require total weight < n")
    ap.add_argument("-v", "--verbose", action="store_true", help="print every bin's exact minimum")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count() or 1)
    args = ap.parse_args()
    stride = 1 if args.all else args.stride
    assert stride >= 1 and args.N >= 3

    s, D, WD, A, region, cliques = load(args.cert)
    tot = sum(a[2] for a in A)
    if cliques:
        check_cliques(cliques, args.N)
        tot_cl = sum(cliques['w'])
        print(f"CLIQUE certificate: {len(cliques['w'])} cliques, {len(cliques['boxes'])} boxes, net N={cliques['N']}, Q={cliques['Q']}; "
              f"clique weight {F(tot_cl, WD)}={tot_cl/WD:.7f}; pairwise core intersection checked exactly (exact sigma_k)")
        tot += tot_cl
    print(f"certificate {args.cert}:  s={s}={float(s):.6f}  atoms={len(A)}  total weight{' (points + cliques)' if cliques else ''}={F(tot, WD)}={tot/WD:.7f}")
    kdot = region['kdot'] if region else 0
    kk = region['ktot'] if region else 0
    if region:
        if region['single']:
            print(f"BRANCH certificate: corner boxes [0,{region['r']}]^2 and images, lambda={F(region['lam'][0], WD)}={region['lam'][0]/WD:+.7f}, k={kk}; "
                  f"total - lambda*k = {F(tot - kdot, WD)} = {(tot - kdot)/WD:.7f}")
        else:
            print(f"BRANCH certificate: corner boxes [0,{region['r']}]^2 and images, PER BOX lambda={[l/WD for l in region['lam']]}, k={region['k']}; "
                  f"total - sum lambda_j k_j = {F(tot - kdot, WD)} = {(tot - kdot)/WD:.7f}")
    weight_ok = True
    if args.n is not None:
        weight_ok = tot - kdot < args.n * WD
        print(f"total{' - lambda.k' if region else ''} < n={args.n}: {'yes' if weight_ok else 'NO'}")
    sD = s * D; assert sD.denominator == 1, "s*D must be an integer"
    sym = d4_symmetric(A, int(sD)) and (region is None or len(set(region['lam'])) == 1)   # unequal per-box lambdas break the symmetry of the claim
    sym = sym and cliques is None            # clique boxes are tied to the net: no representable images under the symmetries
    N = args.N
    if sym:
        K = 0
        while (K + N)**2 < 2*N*N: K += 1     # theta_K >= 45 deg
        print(f"D4 symmetry of the atom set: OK  -> angles reduce to [0,45deg]: bins k=0..{K-1} (N={N})")
    else:
        K = N                                # theta_N = 90 deg
        print(f"{'clique block present' if cliques else 'atom set is NOT D4-symmetric'}  -> full range [0,90deg): bins k=0..{K-1} (N={N})")
    ks = list(range(0, K, stride))
    exhaustive = (stride == 1)
    print(f"{'EXHAUSTIVE' if exhaustive else 'SAMPLED'}: checking {len(ks)} of {K} bins with {args.jobs} worker(s)")
    G.update(N=N, s=s, D=D, A=A, WD=WD, region=region, cliques=cliques)

    t0 = time.time()
    res = []
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(args.jobs) as pool:
        done = 0; nxt = 0.1
        for r in pool.imap_unordered(min_cover_bin, ks, chunksize=max(1, len(ks)//(8*args.jobs))):
            res.append(r); done += 1
            if done >= nxt * len(ks) and not args.verbose and sys.stderr.isatty():
                print(f"  ... {done}/{len(ks)} bins, {time.time()-t0:.0f}s", file=sys.stderr); nxt += 0.1
    res.sort()
    worst = None; worst_k = None; nfail = 0
    for k, v, wit in res:
        if v is None:
            if args.verbose: print(f"  bin k={k:5d}  theta={degrees(2*atan(k/N)):7.3f}deg   no unit square fits (vacuous)")
            continue
        if args.verbose:
            print(f"  bin k={k:5d}  theta={degrees(2*atan(k/N)):7.3f}deg   exact min covered weight = {F(v, WD)} = {v/WD:.7f}  {'OK' if v >= WD else '*** FAIL ***'}")
        if v < WD:
            nfail += 1
            if not args.verbose and nfail <= 8:
                cx, cy = witness_xy(k, N, wit)
                print(f"  FAIL at bin k={k} (theta={degrees(2*atan(k/N)):.3f}deg): covered {F(v, WD)} = {v/WD:.7f}  near centre ({cx:.4f},{cy:.4f})")
        if worst is None or v < worst: worst, worst_k = v, k
    print(f"wall time {time.time()-t0:.1f}s")
    if worst is None:
        print("no admissible placements at all (container smaller than a unit square)")
        worst = 0
    print(f"minimum covered weight{' minus region threshold' if region else ''} over {'ALL' if exhaustive else 'SAMPLED'} bins = {F(worst, WD)} = {worst/WD:.7f}   (at bin k={worst_k})")
    ok = worst >= WD and weight_ok
    if ok:
        qual = "" if exhaustive else "  (sampled bins only -- not a proof)"
        nn = f", and total weight < {args.n}" if args.n is not None else " (total weight vs n not checked: pass --n)"
        if region:
            ktxt = str(kk) if region['single'] else "".join(str(v) for v in region['k'])
            print(f"VERIFIED: (branch k={ktxt}) every closed unit square inside [0,{s}]^2 centred in corner box j covers weight >= 1 + lambda_j, every other one >= 1{nn.replace('total weight', 'total - lambda.k')}{qual}")
        else:
            print(f"VERIFIED: every closed unit square inside [0,{s}]^2 covers weight >= 1{' (points plus cliques)' if cliques else ''}{nn}{qual}")
    else:
        print("NOT VERIFIED")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
