#!/usr/bin/env python3
"""Independent stress test of the two anchor-clique predicates (task I).

The verifier credits a cell of the sweep a clique's weight only if one of its pieces holds EVERY
pose of the cell, which it decides with two exact predicates (`certificates/FORMAT.md`, "Anchor
cliques"): `contains` by the exact bin core and `meets` by the hexagon `A + [-h,h]^2`.  Both are
*sufficient* conditions, so the only thing that can go wrong is that one of them says yes when some
pose of the cell does not in fact contain / meet the anchor.  This script hunts for exactly that:

    python3 search/anchorstress.py [--n 20000] [--N 2000] [--samples 40] [--seed 0]

It draws a random bin, a random cell of centres and a random anchor, asks `xcheck.py`'s predicates
(the exact rational ones, which are the *inflated* version of the Rust's -- anything the Rust
credits, they credit), and where they say yes it samples poses of the cell and tests the geometry
directly, from the definition, with code that shares nothing with either checker: an explicit
rotation for `contains` and Liang-Barsky clipping of the segment against the square for `meets`.
Any failure is a soundness bug.  It also reports how often the predicates are conservative (all
sampled poses fine but the predicate said no), which is the price of the sufficient condition.
"""
import sys, os, math, random, argparse
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import xcheck as X


def theta(k, N):
    return 2.0 * math.atan(k / N)


def contains_pose(z, cx, cy, th):
    """does the closed UNIT square at (cx, cy, th) contain the point z?  (direct, floats)"""
    dx, dy = z[0] - cx, z[1] - cy
    c, s = math.cos(th), math.sin(th)
    return abs(dx * c + dy * s) <= 0.5 and abs(-dx * s + dy * c) <= 0.5


def meets_pose(p, q, cx, cy, th):
    """does the closed unit square at (cx, cy, th) meet the closed segment pq?  Liang-Barsky
    clipping of the segment against the axis-parallel square in the square's own frame."""
    c, s = math.cos(th), math.sin(th)
    def loc(z):
        dx, dy = z[0] - cx, z[1] - cy
        return (dx * c + dy * s, -dx * s + dy * c)
    (x0, y0), (x1, y1) = loc(p), loc(q)
    dx, dy = x1 - x0, y1 - y0
    t0, t1 = 0.0, 1.0
    for (pp, qq) in ((-dx, x0 + 0.5), (dx, 0.5 - x0), (-dy, y0 + 0.5), (0 + dy, 0.5 - y0)):
        if pp == 0.0:
            if qq < 0: return False
        else:
            r = qq / pp
            if pp < 0:
                if r > t1: return False
                if r > t0: t0 = r
            else:
                if r < t0: return False
                if r < t1: t1 = r
    return t0 <= t1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=20000); ap.add_argument('--N', type=int, default=2000)
    ap.add_argument('--samples', type=int, default=40); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--s', default='199/50'); ap.add_argument('--D', type=int, default=784000)
    a = ap.parse_args()
    random.seed(a.seed)
    s = F(a.s); D = a.D; N = a.N
    ncon = nmeet = 0; failc = failm = 0; consc = consm = 0
    for it in range(a.n):
        k = random.randrange(0, N)
        c, sn, sigma, wmin, cd, sd = X.bin_params(k, N, full=True)
        h = sigma / 2
        th0, th1 = theta(k, N), theta(k + 1, N)
        # a random anchor: a point, or a segment of random direction and length
        px = F(random.randrange(0, int(s * D) + 1), D); py = F(random.randrange(0, int(s * D) + 1), D)
        if random.random() < 0.5:
            anc = ((px, py), (px, py))
        else:
            L = random.choice([1, 3, 10, 40, 200, 1000]) * random.randrange(1, 200)
            dxn = random.randrange(-L, L + 1); dyn = random.choice([-1, 1]) * int(math.isqrt(max(L * L - dxn * dxn, 0)))
            qx = min(max(px + F(dxn, D), F(0)), s); qy = min(max(py + F(dyn, D), F(0)), s)
            anc = ((px, py), (qx, qy))
        # a random cell of centres in the frame of the bin, near the anchor
        au = [(c * z[0] + sn * z[1], -sn * z[0] + c * z[1]) for z in anc]
        w = F(random.choice([1, 10, 100, 1000, 20000]), D)
        u0a = au[0][0] + F(random.randrange(-700000, 700001), 10 ** 6)
        u1a = au[0][1] + F(random.randrange(-700000, 700001), 10 ** 6)
        u0b, u1b = u0a + w, u1a + w
        # the two predicates, exactly as xcheck.py evaluates them for a cell
        okc = all(X.in_core(e[0] - u0, e[1] - u1, cd, sd)
                  for e in au for u0 in (u0a, u0b) for u1 in (u1a, u1b))
        nx, ny = anc[0][1] - anc[1][1], anc[1][0] - anc[0][0]
        n0, n1 = c * nx + sn * ny, -sn * nx + c * ny
        okm = (u0b <= max(au[0][0], au[1][0]) + h and u0a >= min(au[0][0], au[1][0]) - h
               and u1b <= max(au[0][1], au[1][1]) + h and u1a >= min(au[0][1], au[1][1]) - h)
        if okm and (n0 or n1):
            P = au[0]
            t0 = (n0 * (u0a - P[0]), n0 * (u0b - P[0])); t1 = (n1 * (u1a - P[1]), n1 * (u1b - P[1]))
            bnd = h * (abs(n0) + abs(n1))
            okm = not (min(t0) + min(t1) < -bnd or max(t0) + max(t1) > bnd)
        # sample poses of the cell and test the geometry directly
        cf = float(c); snf = float(sn)
        badc = badm = False
        for _ in range(a.samples):
            t = random.random(); u0 = float(u0a) + t * float(w)
            t = random.random(); u1 = float(u1a) + t * float(w)
            for u0v, u1v in [(u0, u1), (float(u0a), float(u1a)), (float(u0b), float(u1b)),
                             (float(u0a), float(u1b)), (float(u0b), float(u1a))]:
                cx = cf * u0v - snf * u1v; cy = snf * u0v + cf * u1v
                for th in (th0, th1, th0 + random.random() * (th1 - th0)):
                    if not all(contains_pose(z, cx, cy, th) for z in [(float(z[0]), float(z[1])) for z in anc]):
                        badc = True
                    if not meets_pose((float(anc[0][0]), float(anc[0][1])), (float(anc[1][0]), float(anc[1][1])), cx, cy, th):
                        badm = True
        if okc:
            ncon += 1
            if badc: failc += 1; print(f"  CONTAINS FAILURE at it={it} k={k} anchor={anc} cell u0=[{float(u0a)},{float(u0b)}] u1=[{float(u1a)},{float(u1b)}]")
        elif not badc: consc += 1
        if okm:
            nmeet += 1
            if badm: failm += 1; print(f"  MEETS FAILURE at it={it} k={k} anchor={anc} cell u0=[{float(u0a)},{float(u0b)}] u1=[{float(u1a)},{float(u1b)}]")
        elif not badm: consm += 1
    print(f"{a.n} random (bin, cell, anchor) cases at N={N}, {a.samples} pose samples each "
          f"(x5 cell points x3 angles):")
    print(f"  contains: {ncon} credited, {failc} FAILURES; {consc} conservative (no pose violated but the predicate said no)")
    print(f"  meets:    {nmeet} credited, {failm} FAILURES; {consm} conservative")
    print("VERDICT:", "OK" if failc == 0 and failm == 0 else "*** SOUNDNESS FAILURE ***")
    sys.exit(0 if failc == 0 and failm == 0 else 1)


if __name__ == '__main__':
    main()
