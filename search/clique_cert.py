#!/usr/bin/env python3
"""EXACT certification of the anchor-clique property  (task G, part 1).

The object.  For a point p and a compact convex anchor A (here a closed segment), the family

    K(p, A) = { S admissible : p in S }  u  { S admissible : A subset S }

is a clique -- every two of its members closed-intersect -- as soon as

    (T)   every admissible closed unit square S with p in S meets A.

(Two members of the first part share p; two of the second share A; a mixed pair S, S' has
S n S' contained... contains S n A, nonempty by (T).)  So certifying the clique reduces to the
three-dimensional statement (T), which this script proves by exact rational subdivision of pose
space -- the same machinery as search/zeromargin.py, with the same conventions:

    pose (cx, cy, u), u = tan(theta/2), cos = (1-u^2)/(1+u^2), sin = 2u/(1+u^2) (exact
    rationals at the box corners), admissible iff cx, cy in [w/2, t - w/2], w = cos + sin
    (theta in [0, 90 deg], where cos, sin >= 0).

A box is a leaf when one of
    EMPTY  no admissible pose in the box has p in S;
    MEET   every pose in the box (after clipping to cx, cy in [w_lo/2, t - w_lo/2], which keeps
           every admissible pose) has S meeting A.
MEET is decided by the separating-axis theorem: a convex polygon and a segment are disjoint iff
one of three axes separates them -- the segment's normal d, and the square's two edge normals.
The test refutes all six (axis, side) separations at once over the box with interval arithmetic
in exact rationals.

Termination.  If rho > rho*(p, eps) strictly (Lemma 2 of notes/clique-family.md) the statement
holds with a uniform positive margin and the subdivision terminates.  At rho = rho*(p, eps)
exactly -- the maximal family -- it does NOT: the pose  theta = arccos(p_x),  cx = w(theta)/2
touches A in a single point and every box around it contains inadmissible poses that separate.
That configuration is what Lemma 2 handles in closed form; there is no way to reach it by
subdivision alone.  `--rho-star` reports the exact rho*(p, eps) as a rational bound.

Usage
    python3 search/clique_cert.py cert --px 99/100 --py 2 --eps 1/250 --rho 3/50   [--depth 26]
    python3 search/clique_cert.py stress  [--n 200000]
"""
import sys, os, math, time, json, argparse
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')

T_DEFAULT = Fr(399, 100)


# ------------------------------------------------------------------ rational interval arithmetic
class Iv:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi

    def __add__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __mul__(self, o):
        o = o if isinstance(o, Iv) else Iv(o)
        v = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(v), max(v))

    __rmul__ = __mul__
    __radd__ = __add__

    def __repr__(self):
        return f"[{float(self.lo):.6g},{float(self.hi):.6g}]"


def cs_interval(u0, u1):
    """exact cos, sin intervals for theta with tan(theta/2) in [u0, u1] subset [0, 1]"""
    c = lambda u: Fr(1 - u * u, 1) / (1 + u * u)
    s = lambda u: Fr(2 * u, 1) / (1 + u * u)
    return Iv(c(u1), c(u0)), Iv(s(u0), s(u1))          # cos decreasing, sin increasing on [0,1]


def w_lo_exact(u0, u1):
    """min of w = cos + sin over the bin (w is unimodal on [0, 90 deg], max at u = sqrt2 - 1)"""
    wf = lambda u: (1 - u * u + 2 * u) / (1 + u * u)
    return min(wf(u0), wf(u1))


# ------------------------------------------------------------------ the box test
class Problem:
    def __init__(self, t, p, a0, a1):
        self.t = t
        self.p = p
        self.a0 = a0
        self.a1 = a1
        dx = a1[0] - a0[0]
        dy = a1[1] - a0[1]
        self.d = (-dy, dx)                       # a normal of the segment (unnormalised)
        self.md0 = a0[0] * self.d[0] + a0[1] * self.d[1]
        assert self.md0 == a1[0] * self.d[0] + a1[1] * self.d[1]
        self.dabs = (abs(self.d[0]), abs(self.d[1]))

    def test(self, X0, X1, Y0, Y1, u0, u1):
        """returns 'EMPTY', 'MEET', or None (subdivide)"""
        t = self.t
        C, S = cs_interval(u0, u1)
        wlo = w_lo_exact(u0, u1)
        # clip to the admissible slab (keeps every admissible pose of the box)
        lo = Fr(wlo, 2)
        hi = t - lo
        X0 = max(X0, lo)
        X1 = min(X1, hi)
        Y0 = max(Y0, lo)
        Y1 = min(Y1, hi)
        if X0 > X1 or Y0 > Y1:
            return 'EMPTY'
        CX = Iv(X0, X1)
        CY = Iv(Y0, Y1)
        px, py = self.p
        # p in S ?
        e1p = (Iv(px) - CX) * C + (Iv(py) - CY) * S
        e2p = (Iv(px) - CX) * (-S) + (Iv(py) - CY) * C
        if e1p.lo > Fr(1, 2) or e1p.hi < Fr(-1, 2) or e2p.lo > Fr(1, 2) or e2p.hi < Fr(-1, 2):
            return 'EMPTY'
        # ---- SAT axis 0: the segment's normal d.  Separated iff
        #      max_S <x,d> < <A,d>  or  min_S <x,d> > <A,d>
        wd = Iv(self.dabs[0]) * C + Iv(self.dabs[1]) * S          # |<e1,d>|, both cos,sin >= 0
        wd = wd + (Iv(self.dabs[0]) * S + Iv(self.dabs[1]) * C)   # + |<e2,d>|
        cd = CX * self.d[0] + CY * self.d[1]
        hiS = cd + Fr(1, 2) * wd
        loS = cd - Fr(1, 2) * wd
        if not (hiS.lo >= self.md0 and loS.hi <= self.md0):
            return None
        # ---- SAT axes 1, 2: the square's edge normals
        for which in (1, 2):
            okhi = False
            oklo = False
            for a in (self.a0, self.a1):
                if which == 1:
                    v = (Iv(a[0]) - CX) * C + (Iv(a[1]) - CY) * S
                else:
                    v = (Iv(a[0]) - CX) * (-S) + (Iv(a[1]) - CY) * C
                if v.hi <= Fr(1, 2):
                    okhi = True
                if v.lo >= Fr(-1, 2):
                    oklo = True
            if not (okhi and oklo):
                return None
        return 'MEET'


def certify(prob, root_pitch=Fr(1, 10), ubins=8, depth=26, verbose=True, cap=4_000_000):
    """adaptive subdivision; returns a dict of counts and the uncertified boxes"""
    t = prob.t
    stack = []
    nx = int(math.ceil(float(t) / float(root_pitch)))
    for i in range(nx):
        for j in range(nx):
            for k in range(ubins):
                stack.append((min(i * root_pitch, t), min((i + 1) * root_pitch, t),
                              min(j * root_pitch, t), min((j + 1) * root_pitch, t),
                              Fr(k, ubins), Fr(k + 1, ubins), 0))
    cnt = dict(EMPTY=0, MEET=0, boxes=0)
    bad = []
    maxdepth = 0
    while stack:
        X0, X1, Y0, Y1, u0, u1, dep = stack.pop()
        cnt['boxes'] += 1
        maxdepth = max(maxdepth, dep)
        if cnt['boxes'] > cap:
            bad.append(('CAP', X0, X1, Y0, Y1, u0, u1, dep))
            break
        r = prob.test(X0, X1, Y0, Y1, u0, u1)
        if r is not None:
            cnt[r] += 1
            continue
        if dep >= depth:
            bad.append(('DEPTH', X0, X1, Y0, Y1, u0, u1, dep))
            continue
        # split the longest scaled side (u scaled by 2 -- one unit of u is ~2 radians)
        lx = X1 - X0
        ly = Y1 - Y0
        lu = 2 * (u1 - u0)
        if lx >= ly and lx >= lu:
            m = (X0 + X1) / 2
            stack.append((X0, m, Y0, Y1, u0, u1, dep + 1))
            stack.append((m, X1, Y0, Y1, u0, u1, dep + 1))
        elif ly >= lu:
            m = (Y0 + Y1) / 2
            stack.append((X0, X1, Y0, m, u0, u1, dep + 1))
            stack.append((X0, X1, m, Y1, u0, u1, dep + 1))
        else:
            m = (u0 + u1) / 2
            stack.append((X0, X1, Y0, Y1, u0, m, dep + 1))
            stack.append((X0, X1, Y0, Y1, m, u1, dep + 1))
    cnt['depth'] = maxdepth
    cnt['uncertified'] = len(bad)
    return cnt, bad


# ------------------------------------------------------------------ rho* (exact, Lemma 2)
def rho_star_exact(px):
    """rho*/eps = px / sqrt(1 - px^2) -- returned as an exact rational UPPER bound"""
    from fractions import Fraction
    v = float(px) / math.sqrt(1 - float(px) ** 2)
    return Fraction(v).limit_denominator(10 ** 9) + Fraction(1, 10 ** 9)


# ------------------------------------------------------------------ independent float stress
def stress(n=200000, seed=7, t=3.99):
    """random (p, eps, rho, pose): check Lemma 2 and the transversal claim directly"""
    import random
    rng = random.Random(seed)
    bad_lemma = 0
    bad_meet = 0
    tested = 0
    for _ in range(n):
        px = rng.uniform(0.30, 0.999)
        py = rng.uniform(1.5, 2.5)
        kappa = px / math.sqrt(1 - px * px)
        eps = rng.uniform(1e-4, 0.999) * (1 - px)
        slack = rng.choice([1.001, 1.01, 1.1, 1.5])
        rho = kappa * eps * slack
        if rho >= 0.5:
            continue
        x0 = px + eps
        a0 = (x0, py - rho)
        a1 = (x0, py + rho)
        # a random admissible pose through p
        th = rng.uniform(0.0, math.pi / 2)
        w2 = (abs(math.cos(th)) + abs(math.sin(th))) / 2
        for _try in range(60):
            cx = rng.uniform(max(w2, px - 0.71), min(t - w2, px + 0.71))
            cy = rng.uniform(max(w2, py - 0.71), min(t - w2, py + 0.71))
            ct, st = math.cos(th), math.sin(th)
            if abs((px - cx) * ct + (py - cy) * st) <= 0.5 and abs(-(px - cx) * st + (py - cy) * ct) <= 0.5:
                break
        else:
            continue
        tested += 1
        # independent test: does the segment meet the square?
        if not seg_meets(a0, a1, (cx, cy), th):
            bad_meet += 1
        # sharpness: with rho' slightly BELOW kappa*eps the witness pose of the proof
        # (theta0 = arccos px, p at the corner with s1 = s2 = 1/2) must be admissible,
        # contain p, and MISS the segment -- provided py >= w(theta0)
        th0 = math.acos(px)
        c0, s0 = px, math.sqrt(1 - px * px)
        w0 = c0 + s0
        if py >= w0 and py <= t - w0:
            cx0 = px - 0.5 * c0 + 0.5 * s0
            cy0 = py - 0.5 * w0
            rho2 = kappa * eps * 0.999
            b0 = (x0, py - rho2)
            b1 = (x0, py + rho2)
            inside = (abs((px - cx0) * c0 + (py - cy0) * s0) <= 0.5 + 1e-12 and
                      abs(-(px - cx0) * s0 + (py - cy0) * c0) <= 0.5 + 1e-12)
            adm = (cx0 >= w0 / 2 - 1e-12 and cy0 >= w0 / 2 - 1e-12
                   and cx0 <= t - w0 / 2 + 1e-12 and cy0 <= t - w0 / 2 + 1e-12)
            if inside and adm and seg_meets(b0, b1, (cx0, cy0), th0):
                bad_lemma += 1
    return tested, bad_meet, bad_lemma


def seg_meets(a0, a1, c, th):
    ct, st = math.cos(th), math.sin(th)
    u0 = (a0[0] - c[0]) * ct + (a0[1] - c[1]) * st
    v0 = -(a0[0] - c[0]) * st + (a0[1] - c[1]) * ct
    u1 = (a1[0] - c[0]) * ct + (a1[1] - c[1]) * st
    v1 = -(a1[0] - c[0]) * st + (a1[1] - c[1]) * ct
    lo, hi = 0.0, 1.0
    for (q0, q1) in ((u0, u1), (v0, v1)):
        dq = q1 - q0
        for (sgn, bnd) in ((1.0, 0.5), (-1.0, 0.5)):
            num = bnd - sgn * q0
            den = sgn * dq
            if abs(den) < 1e-300:
                if num < 0:
                    return False
                continue
            r = num / den
            if den > 0:
                hi = min(hi, r)
            else:
                lo = max(lo, r)
            if lo > hi:
                return False
    return True


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('cert')
    c.add_argument('--t', default='399/100')
    c.add_argument('--px', default='99/100')
    c.add_argument('--py', default='2')
    c.add_argument('--eps', default='1/250')
    c.add_argument('--rho', default=None, help='default: 1.05 * rho*(px) * eps')
    c.add_argument('--depth', type=int, default=26)
    c.add_argument('--pitch', default='1/10')
    c.add_argument('--ubins', type=int, default=8)
    c.add_argument('--cap', type=int, default=4000000)
    s = sub.add_parser('stress')
    s.add_argument('--n', type=int, default=200000)
    args = ap.parse_args()
    if args.cmd == 'stress':
        t0 = time.time()
        n, bad, badl = stress(args.n)
        print(f"stress: {n} random (p, eps, rho > rho*, admissible pose through p): "
              f"{bad} failures of 'S meets A'; {badl} failures of the sharpness witness at "
              f"rho = 0.999 rho*  ({time.time()-t0:.1f} s)")
        return
    t = Fr(args.t)
    px = Fr(args.px)
    py = Fr(args.py)
    eps = Fr(args.eps)
    if args.rho:
        rho = Fr(args.rho)
    else:
        rho = Fr(21, 20) * rho_star_exact(px) * eps
        rho = Fr(rho).limit_denominator(10 ** 7)
    p = (px, py)
    a0 = (px + eps, py - rho)
    a1 = (px + eps, py + rho)
    print(f"# t = {t}, p = ({px}, {py}), A = vertical segment x = {px + eps}, "
          f"y in [{py - rho}, {py + rho}]")
    print(f"#   eps = {eps} = {float(eps):.6f} (< 1 - px = {float(1-px):.6f}), "
          f"rho = {rho} = {float(rho):.6f}, rho* = {float(rho_star_exact(px)*eps):.6f}, "
          f"slack = {float(rho/(rho_star_exact(px)*eps)):.4f}")
    prob = Problem(t, p, a0, a1)
    t0 = time.time()
    cnt, bad = certify(prob, root_pitch=Fr(args.pitch), ubins=args.ubins, depth=args.depth,
                       cap=args.cap)
    print(f"# boxes {cnt['boxes']}, depth {cnt['depth']}, EMPTY {cnt['EMPTY']}, "
          f"MEET {cnt['MEET']}, uncertified {cnt['uncertified']}  ({time.time()-t0:.1f} s)")
    if bad:
        print("# uncertified boxes (first 10):")
        for b in bad[:10]:
            print("   ", b[0], [float(z) for z in b[1:7]], "depth", b[7])
        print("VERDICT: NOT CERTIFIED")
    else:
        print("VERDICT: CERTIFIED -- K(p, A) is a clique (every admissible S through p meets A)")


if __name__ == '__main__':
    main()
