#!/usr/bin/env python3
"""unavoid13: zeromargin.py's exact checker plus one more primitive, SEG (unit-segment lemma),
without modifying search/zeromargin.py (imported, subclassed).

    python3 search/unavoid13_check.py cert FILE [--tri] [--seg] [--full] [--depth D] [--nproc N]
                                       [--dump OUT] [--oracle OUT] [--pitch 1/10] [--ubins 8]

SEG lemma [proved].  Let p, q be points with |q - p| <= 1 and let X be a set of poses (a box) such
that for every pose in X, in the square's own frame (x' along e1 = (cos t, sin t), y' along
e2 = (-sin t, cos t), origin at the centre):  |x'_p| <= 1/2, |x'_q| <= 1/2, y'_p <= 1/2 and
y'_q >= -1/2.  Then every pose in X contains p or q.
Proof: if p is not in Q then, since |x'_p| <= 1/2 and y'_p <= 1/2, y'_p < -1/2; then
y'_q = y'_p + <q - p, e2> < -1/2 + |q - p| <= 1/2, and y'_q >= -1/2 with |x'_q| <= 1/2 puts q in Q.
The four variants (p <-> q, x' <-> y') are all tried.  This is the two-point sibling of TRI and is
what certifies a wall pose whose witness switches between two points on the same vertical/horizontal
unit segment (e.g. Kearney-Shiu's (a,1), (a,2) at the wall square (1/2, 3/2, theta -> 0)), where every
fixed-witness primitive fails at every depth (RUNG2.md's one-sided poses).

Exactness of the box test: for v = p - c, c in the centre rectangle, theta in [t0, t1]:
max <v, e(theta)> over the box is attained at a corner of the v-rectangle (linear in v) and, for
fixed v, at a bin endpoint unless dir(v) lies in the bin's cone, where it equals |v|; so
"max <= 1/2" is: endpoint dots <= 1/2 at all four corners and |v|^2 <= 1/4 whenever the cone test
(two exact cross products with the rational endpoint directions) passes.  "min >= -1/2" is the same
test for -v.
"""
import sys, os, math, time, argparse, itertools
from fractions import Fraction as F
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zeromargin as Z

HALF = F(1, 2)


def _le_half(v, E0, E1):
    """max over theta in bin of <v, e(theta)> <= 1/2, exact; E0 = e(t0), E1 = e(t1) as rational vectors."""
    vx, vy = v
    if vx * E0[0] + vy * E0[1] > HALF: return False
    if vx * E1[0] + vy * E1[1] > HALF: return False
    # dir(v) in the closed cone from E0 to E1 (cone angle < 180 deg)?
    if (E0[0] * vy - E0[1] * vx) >= 0 and (vx * E1[1] - vy * E1[0]) >= 0:
        return vx * vx + vy * vy <= F(1, 4)
    return True


def _le_half_f(v, E0, E1, tol=1e-9):
    vx, vy = v
    if vx * E0[0] + vy * E0[1] > 0.5 + tol: return False
    if vx * E1[0] + vy * E1[1] > 0.5 + tol: return False
    if (E0[0] * vy - E0[1] * vx) >= -tol and (vx * E1[1] - vy * E1[0]) >= -tol:
        return vx * vx + vy * vy <= 0.25 + tol
    return True


class SegChecker(Z.Checker):
    def __init__(self, *a, use_seg=True, **k):
        super().__init__(*a, **k)
        self.use_seg = use_seg
        self.pairs = []
        if use_seg:
            n = len(self.P)
            for i, j in itertools.combinations(range(n), 2):
                dx = self.P[i][0] - self.P[j][0]; dy = self.P[i][1] - self.P[j][1]
                if dx * dx + dy * dy <= 1 and self.W[i] >= 1 and self.W[j] >= 1:
                    self.pairs.append((i, j))

    def _seg_ok(self, box, p, q, axis, le, exact):
        """axis 0: slab along e1 (x'), free coordinate y' along e2; axis 1: the reverse.
        Checks (b) both in the slab, (c) 'free(p) <= 1/2' and 'free(q) >= -1/2'."""
        cx0, cx1, cy0, cy1, u0, u1 = box
        if exact:
            c0, s0 = Z.trig(u0); c1, s1 = Z.trig(u1)
            E1a, E1b = (c0, s0), (c1, s1)          # e1 at the bin ends
            E2a, E2b = (-s0, c0), (-s1, c1)        # e2 at the bin ends
            neg = lambda v: (-v[0], -v[1])
        else:
            c0, s0 = Z.trig(u0); c1, s1 = Z.trig(u1)
            c0, s0, c1, s1 = float(c0), float(s0), float(c1), float(s1)
            E1a, E1b = (c0, s0), (c1, s1); E2a, E2b = (-s0, c0), (-s1, c1)
            neg = lambda v: (-v[0], -v[1])
            cx0, cx1, cy0, cy1 = float(cx0), float(cx1), float(cy0), float(cy1)
            p = (float(p[0]), float(p[1])); q = (float(q[0]), float(q[1]))
        slabA, slabB = (E1a, E1b) if axis == 0 else (E2a, E2b)     # the coordinate that must be in [-1/2, 1/2]
        freeA, freeB = (E2a, E2b) if axis == 0 else (E1a, E1b)     # the coordinate the segment spans
        for pt, side in ((p, 'le'), (q, 'ge')):
            for ax in (cx0, cx1):
                for ay in (cy0, cy1):
                    v = (pt[0] - ax, pt[1] - ay)
                    if not le(v, slabA, slabB): return False
                    if not le(neg(v), slabA, slabB): return False
                    if side == 'le':
                        if not le(v, freeA, freeB): return False          # free(p) <= 1/2
                    else:
                        if not le(neg(v), freeA, freeB): return False     # free(q) >= -1/2
        return True

    def cert_seg(self, box):
        for (i, j) in self.pairs:
            for (a, b) in ((i, j), (j, i)):
                for axis in (0, 1):
                    if self._seg_ok(box, self.P[a], self.P[b], axis, _le_half_f, exact=False):
                        if self._seg_ok(box, self.P[a], self.P[b], axis, _le_half, exact=True):
                            return ('TRI', ('SEG', a, b, axis))
        return (None, None)

    def cert_tri(self, box):
        kind, wit = super().cert_tri(box) if self.use_tri else (None, None)
        if kind is None and self.use_seg:
            kind, wit = self.cert_seg(box)
        return kind, wit


def _worker(args):
    chk, root = args
    return chk.run_box(root)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what', help='cert')
    ap.add_argument('path')
    ap.add_argument('--tri', action='store_true')
    ap.add_argument('--seg', action='store_true')
    ap.add_argument('--disj', action='store_true')
    ap.add_argument('--depth', type=int, default=16)
    ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--pitch', type=str, default='1/10')
    ap.add_argument('--ubins', type=int, default=8)
    ap.add_argument('--full', action='store_true')
    ap.add_argument('--dump', type=str, default=None)
    ap.add_argument('--oracle', type=str, default=None)
    a = ap.parse_args()
    m, pts, ws = Z.read_cert(a.path)
    # cert_tri is only called when use_tri is set in run_box: force it on and let the subclass
    # decide whether the TRI lemma itself is used.
    chk = SegChecker(m, pts, ws, use_tri=True, max_depth=a.depth, dump=a.dump, use_adm=True,
                     theta_bias=4, use_chain=a.disj, chain_from=0, clip=True, use_seg=a.seg)
    if not a.tri:
        chk.tris = []
    sym = chk.symmetric()
    print(f"container [0,{m}]^2, {len(pts)} points, total weight {float(sum(chk.W)):.6f}, "
          f"symmetric under x->m-x and y->m-y: {sym}, triangles: {len(chk.tris)}, unit pairs: {len(chk.pairs)}")
    if not sym and not a.full:
        print("point set not symmetric: use --full"); sys.exit(2)
    pitch = F(a.pitch)
    if a.full:
        R = Z.roots(m, pitch, a.ubins * 2, cy_max=m)
        R = [(x0, x1, y0, y1, u0 * 2, u1 * 2) for (x0, x1, y0, y1, u0, u1) in R]
    else:
        R = Z.roots(m, pitch, a.ubins)
    print(f"{len(R)} root boxes, depth limit {a.depth}, pitch {pitch}, u-bins {a.ubins}")
    t0 = time.time()
    tot = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'EMPTY': 0, 'UNCERT': 0,
           'boxes': 0, 'maxdepth': 0}
    unc_all = []; leaves_all = []; nseg = 0
    with Pool(a.nproc) as pool:
        for i, (st, unc, leaves) in enumerate(pool.imap_unordered(_worker, [(chk, r) for r in R], chunksize=4)):
            for k in tot:
                tot[k] = max(tot[k], st[k]) if k == 'maxdepth' else tot[k] + st[k]
            unc_all += unc; leaves_all += leaves
            if (i + 1) % 1000 == 0:
                print(f"  {i+1}/{len(R)} roots, {tot['boxes']} boxes, uncert {tot['UNCERT']}, {time.time()-t0:.0f}s", flush=True)
    nseg = sum(1 for (_, kind, wit) in leaves_all if kind == 'TRI' and isinstance(wit, tuple) and wit and wit[0] == 'SEG')
    print(f"done in {time.time()-t0:.0f}s: boxes {tot['boxes']}, max depth {tot['maxdepth']}")
    print(f"  leaves: ADM {tot['ADM']}  CORE {tot['CORE']}  P1 {tot['P1']}  MIX {tot['MIX']}  "
          f"CHAIN {tot['CHAIN']}  TRI+SEG {tot['TRI']} (of which SEG {nseg if a.dump else 'n/a without --dump'})  "
          f"EMPTY {tot['EMPTY']}  UNCERTIFIED {tot['UNCERT']}")
    if unc_all:
        print("uncertified boxes (cx0 cx1 cy0 cy1 -> theta0 theta1 deg):")
        for b in sorted(unc_all)[:30]:
            print("  ", *[f"{float(v):.6f}" for v in b[:4]],
                  f"{math.degrees(2*math.atan(float(b[4]))):.3f} {math.degrees(2*math.atan(float(b[5]))):.3f}")
        if len(unc_all) > 30: print(f"  ... {len(unc_all)} in total")
    if a.dump:
        with open(a.dump, 'w') as f:
            f.write(f"# container {m}; box: cx0 cx1 cy0 cy1 u0 u1 ; kind ; witness (TRI leaves with witness ('SEG', i, j, axis) are SEG leaves)\n")
            for box, kind, wit in leaves_all:
                f.write(" ".join(str(v) for v in box) + f" ; {kind} ; {wit}\n")
        print(f"leaves written to {a.dump}")
    if a.oracle:
        with open(a.oracle, 'w') as f:
            f.write(f"# container {m}; uncertified box worst-case poses: cx cy theta_rad\n")
            for box in unc_all:
                bx0, bx1, by0, by1, bu0, bu1 = box
                th0, th1 = 2 * math.atan(float(bu0)), 2 * math.atan(float(bu1))
                for cxv in sorted({float(bx0), float(bx1), float((bx0 + bx1) / 2)}):
                    for cyv in sorted({float(by0), float(by1), float((by0 + by1) / 2)}):
                        for thv in (th0, th1, (th0 + th1) / 2):
                            f.write(f"{cxv!r} {cyv!r} {thv!r}\n")
        print(f"oracle rows written to {a.oracle}: {len(unc_all)} boxes")
    print("VERIFIED" if tot['UNCERT'] == 0 else "NOT VERIFIED")


if __name__ == '__main__':
    main()
