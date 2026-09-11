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
line VERIFIED or NOT VERIFIED.  A sampled run can only say NOT VERIFIED with certainty; its
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
    assert len(t) == 5 + 3*m, "wrong number of integers in certificate"
    A = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    assert all(w >= 0 for _, _, w in A), "weights must be non-negative"
    return F(sn, sd), D, WD, A

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
    sigma_k-square at angle theta_k.  Returns (k, min_weight_numerator or None, witness)."""
    N, S, D, A = G['N'], G['s'], G['D'], G['A']
    c, sn, sigma, wmin = bin_params(k, N)
    h = sigma / 2
    L = wmin / 2; U = S - wmin / 2
    if U < L:
        return k, None, None        # no unit square at these angles fits in the container
    # rotate by -theta_k: atoms and the centre box
    Q = [(c*F(x, D) + sn*F(y, D), -sn*F(x, D) + c*F(y, D), w) for x, y, w in A]
    poly = [(c*a + sn*b, -sn*a + c*b) for a, b in ((L, L), (U, L), (U, U), (L, U))]
    if L == U:                      # a single admissible centre: evaluate f there directly
        u0, u1 = poly[0]
        tot = sum(w for q0, q1, w in Q if abs(q0 - u0) <= h and abs(q1 - u1) <= h)
        return k, tot, (u0, u1)
    # one common denominator for the whole bin -> integer sweep
    den = 1
    for v in [h] + [q[0] for q in Q] + [q[1] for q in Q] + [p[0] for p in poly] + [p[1] for p in poly]:
        den = lcm(den, v.denominator)
    def I(v):
        v = v * den; assert v.denominator == 1; return v.numerator
    H = I(h)
    atoms = sorted((I(q0), I(q1), w) for q0, q1, w in Q)
    qx = [t[0] for t in atoms]
    V = [(I(px), I(py)) for px, py in poly]
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
        if i1 <= i0:
            return k, 0, (F(a + b, 2*den), (ylo + yhi) / (2*den))   # nothing covered in this strip
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
                return k, 0, (F(a + b, 2*den), (ylo + yhi) / (2*den))   # beyond every atom's window
            cl, dl = by[j], by[j+1]
            k0 = bisect_left(yv, dl - H); k1 = bisect_right(yv, cl + H)
            tot = pre[k1] - pre[k0] if k1 > k0 else 0
            if best is None or tot < best:
                best = tot
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

    s, D, WD, A = load(args.cert)
    tot = sum(a[2] for a in A)
    print(f"certificate {args.cert}:  s={s}={float(s):.6f}  atoms={len(A)}  total weight={F(tot, WD)}={tot/WD:.7f}")
    weight_ok = True
    if args.n is not None:
        weight_ok = tot < args.n * WD
        print(f"total weight < n={args.n}: {'yes' if weight_ok else 'NO'}")
    sD = s * D; assert sD.denominator == 1, "s*D must be an integer"
    sym = d4_symmetric(A, int(sD))
    N = args.N
    if sym:
        K = 0
        while (K + N)**2 < 2*N*N: K += 1     # theta_K >= 45 deg
        print(f"D4 symmetry of the atom set: OK  -> angles reduce to [0,45deg]: bins k=0..{K-1} (N={N})")
    else:
        K = N                                # theta_N = 90 deg
        print(f"atom set is NOT D4-symmetric  -> full range [0,90deg): bins k=0..{K-1} (N={N})")
    ks = list(range(0, K, stride))
    exhaustive = (stride == 1)
    print(f"{'EXHAUSTIVE' if exhaustive else 'SAMPLED'}: checking {len(ks)} of {K} bins with {args.jobs} worker(s)")
    G.update(N=N, s=s, D=D, A=A)

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
    print(f"minimum covered weight over {'ALL' if exhaustive else 'SAMPLED'} bins = {F(worst, WD)} = {worst/WD:.7f}   (at bin k={worst_k})")
    ok = worst >= WD and weight_ok
    if ok:
        qual = "" if exhaustive else "  (sampled bins only -- not a proof)"
        nn = f", and total weight < {args.n}" if args.n is not None else " (total weight vs n not checked: pass --n)"
        print(f"VERIFIED: every closed unit square inside [0,{s}]^2 covers weight >= 1{nn}{qual}")
    else:
        print("NOT VERIFIED")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
