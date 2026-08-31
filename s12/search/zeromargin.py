#!/usr/bin/env python3
"""Exact zero-margin checker for closed covers at the container [0,m]^2 (task E, rung 1).

Statement checked (unit weights):  every closed unit square Q, at every angle, contained in the
closed square [0,m]^2, contains at least one point of the set P.

Method: adaptive subdivision of pose space (cx, cy, u = tan(theta/2)) into boxes with rational
endpoints; a box is certified by one of

  CORE  the point p lies in the exact core of the box: for every centre c of the box's centre
        rectangle and every angle in the box's angle bin, p is in the closed unit square (c, theta).
        The core over an angle bin [t0, t1] at a fixed centre is the intersection of the rotated
        squares R_t Q, which equals  R_t0 Q  ∩  R_t1 Q  ∩  {x : dir(x) mod 90deg in [t0,t1] => |x| <= 1/2}
        (support-function argument; see search/ZEROMARGIN.md).  It is convex, so testing the four
        corners of the rectangle p - rect suffices.  All arithmetic is in Fractions.
  P1    the "2x2 box lemma": if |c - p|_inf <= 1 - w(theta)/2 with w = |cos| + |sin| (i.e. the
        square lies in p + [-1,1]^2), then p in S.  Proof: in the rotated frame the coordinate of
        p - c is at most (1 - w/2) w = w - w^2/2 <= 1/2 since (w-1)^2 >= 0.  Applied to the
        admissible poses of the box only (a wall makes the constraint on that side automatic).
        Exact at the wall/corner tight poses, where CORE never terminates (quadratic margin).
  TRI   (optional) the triangle lemma: a triangle with vertices in P and sides <= 1 containing the
        box's centre rectangle certifies it for every angle (Friedman/Stromquist).
  EMPTY the box contains no admissible pose.

Only admissible poses matter (cx, cy in [w/2, m - w/2]); the centre rectangle is clipped to the
admissible range at the bin's smallest w before CORE, which over-tests inadmissible poses in the
sliver (sound), and P1 handles the slivers at the tight wall poses.

Symmetry: if P is invariant under x -> m-x and y -> m-y (checked exactly), the fundamental domain
is theta in [0, 45deg] (x-reflection sends theta to -theta) and cy <= m/2 (the 180deg rotation).
The root domain uses u in [0, 1/2] (theta up to 53deg) to keep endpoints rational.

Usage:  python3 search/zeromargin.py friedman14 [--tri] [--depth 16] [--nproc 4] [--dump FILE]
        python3 search/zeromargin.py cert certificates/xxx.txt ...   (weighted: sum of certified weights >= 1)
"""
import sys, os, time, argparse, itertools
from fractions import Fraction as F
from multiprocessing import Pool

HALF = F(1, 2)

def trig(u):
    """cos, sin of theta = 2 atan u, exact."""
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d

def bin_data(u0, u1):
    c0, s0 = trig(u0); c1, s1 = trig(u1)
    cD = c1 * c0 + s1 * s0          # cos(t1 - t0)
    sD = s1 * c0 - c1 * s0          # sin(t1 - t0)  (>= 0)
    w0, w1 = c0 + s0, c1 + s1
    # max of w = cos+sin over the bin: at 45deg if the bin contains it (u = sqrt2-1: u^2+2u-1 = 0)
    if u0 * u0 + 2 * u0 - 1 < 0 < u1 * u1 + 2 * u1 - 1:
        whi = F(14143, 10000)       # > sqrt 2
    else:
        whi = max(w0, w1)
    wlo = min(w0, w1)
    return dict(c0=c0, s0=s0, c1=c1, s1=s1, cD=cD, sD=sD, whi=whi, wlo=wlo)

def in_rot_square(a, b, c, s):
    """(a,b) in R_theta Q  (closed unit square centred at 0, angle theta with cos c, sin s)."""
    x = a * c + b * s
    y = -a * s + b * c
    return -HALF <= x <= HALF and -HALF <= y <= HALF

def in_core(a, b, B):
    """(a,b) in the intersection of R_t Q over t in [t0, t1]."""
    if not in_rot_square(a, b, B['c0'], B['s0']): return False
    if not in_rot_square(a, b, B['c1'], B['s1']): return False
    # direction of (a,b) relative to t0, folded mod 90deg into the first quadrant
    x = a * B['c0'] + b * B['s0']
    y = -a * B['s0'] + b * B['c0']
    if x >= 0 and y >= 0: p, q = x, y
    elif x < 0 and y >= 0: p, q = y, -x
    elif x < 0 and y < 0: p, q = -x, -y
    else: p, q = -y, x
    in_sector = (p == 0) or (q * B['cD'] <= p * B['sD'])
    if in_sector:
        return a * a + b * b <= F(1, 4)
    return True

def in_core_f(a, b, Bf, tol=1e-9):
    """float version, lenient by tol (a superset of the exact test)."""
    c0, s0, c1, s1, cD, sD = Bf
    h = 0.5 + tol
    x = a * c0 + b * s0; y = -a * s0 + b * c0
    if not (-h <= x <= h and -h <= y <= h): return False
    x1 = a * c1 + b * s1; y1 = -a * s1 + b * c1
    if not (-h <= x1 <= h and -h <= y1 <= h): return False
    if x >= 0 and y >= 0: p, q = x, y
    elif x < 0 and y >= 0: p, q = y, -x
    elif x < 0 and y < 0: p, q = -x, -y
    else: p, q = -y, x
    if p <= tol or q * cD <= p * sD + tol:
        return a * a + b * b <= 0.25 + tol
    return True

class Checker:
    def __init__(self, m, points, weights=None, use_tri=False, max_depth=16, dump=None):
        self.m = F(m)
        self.P = [(F(x), F(y)) for x, y in points]
        self.W = [F(w) for w in weights] if weights else [F(1)] * len(self.P)
        self.Pf = [(float(x), float(y)) for x, y in self.P]
        self.use_tri = use_tri
        self.max_depth = max_depth
        self.tris = self._triangles() if use_tri else []
        self.dump = dump

    def _triangles(self):
        n = len(self.P); T = []
        for i, j, k in itertools.combinations(range(n), 3):
            ok = True
            for a, b in ((i, j), (j, k), (i, k)):
                dx = self.P[a][0] - self.P[b][0]; dy = self.P[a][1] - self.P[b][1]
                if dx * dx + dy * dy > 1: ok = False; break
            if ok: T.append((self.P[i], self.P[j], self.P[k]))
        return T

    # ---- symmetry -------------------------------------------------------------------------
    def symmetric(self):
        S = set(self.P)
        m = self.m
        return all((m - x, y) in S for x, y in S) and all((x, m - y) in S for x, y in S)

    # ---- certification tests ----------------------------------------------------------------
    def cert_core(self, box, B, Bf):
        cx0, cx1, cy0, cy1 = box[:4]
        # clip to the widest admissible range over the bin (w minimal)
        lo = B['wlo'] / 2; hi = self.m - lo
        cx0 = max(cx0, lo); cx1 = min(cx1, hi); cy0 = max(cy0, lo); cy1 = min(cy1, hi)
        if cx0 > cx1 or cy0 > cy1: return ('EMPTY', None)
        cxf = (float(cx0), float(cx1)); cyf = (float(cy0), float(cy1))
        total = F(0); used = []
        for k, ((px, py), (pxf, pyf)) in enumerate(zip(self.P, self.Pf)):
            okf = True
            for a in (pxf - cxf[0], pxf - cxf[1]):
                for b in (pyf - cyf[0], pyf - cyf[1]):
                    if not in_core_f(a, b, Bf): okf = False; break
                if not okf: break
            if not okf: continue
            ok = True
            for a in (px - cx0, px - cx1):
                for b in (py - cy0, py - cy1):
                    if not in_core(a, b, B): ok = False; break
                if not ok: break
            if ok:
                total += self.W[k]; used.append(k)
                if total >= 1: return ('CORE', used)
        return (None, None)

    def cert_p1(self, box, B):
        cx0, cx1, cy0, cy1 = box[:4]
        t = 1 - B['whi'] / 2
        m = self.m
        total = F(0); used = []
        for k, (px, py) in enumerate(self.P):
            if not (px - 1 <= 0 or cx0 >= px - t): continue
            if not (px + 1 >= m or cx1 <= px + t): continue
            if not (py - 1 <= 0 or cy0 >= py - t): continue
            if not (py + 1 >= m or cy1 <= py + t): continue
            total += self.W[k]; used.append(k)
            if total >= 1: return ('P1', used)
        return (None, None)

    def cert_tri(self, box):
        cx0, cx1, cy0, cy1 = box[:4]
        corners = [(cx0, cy0), (cx1, cy0), (cx0, cy1), (cx1, cy1)]
        for (A, Bp, C) in self.tris:
            if all(self._in_tri(c, A, Bp, C) for c in corners):
                wmin = min(self.W[self.P.index(v)] for v in (A, Bp, C))
                if wmin >= 1: return ('TRI', (A, Bp, C))
        return (None, None)

    @staticmethod
    def _in_tri(c, A, B, C):
        def cross(o, p, q): return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])
        d1 = cross(A, B, c); d2 = cross(B, C, c); d3 = cross(C, A, c)
        return (d1 >= 0 and d2 >= 0 and d3 >= 0) or (d1 <= 0 and d2 <= 0 and d3 <= 0)

    # ---- recursion --------------------------------------------------------------------------
    def run_box(self, root):
        """Certify one root box; returns stats dict and the list of uncertified boxes."""
        stats = {'CORE': 0, 'P1': 0, 'TRI': 0, 'EMPTY': 0, 'UNCERT': 0, 'boxes': 0, 'maxdepth': 0}
        unc = []; leaves = []
        stack = [(root, 0)]
        bincache = {}
        while stack:
            box, depth = stack.pop()
            stats['boxes'] += 1
            stats['maxdepth'] = max(stats['maxdepth'], depth)
            cx0, cx1, cy0, cy1, u0, u1 = box
            key = (u0, u1)
            if key not in bincache:
                B = bin_data(u0, u1)
                Bf = tuple(float(B[k]) for k in ('c0', 's0', 'c1', 's1', 'cD', 'sD'))
                bincache[key] = (B, Bf)
            B, Bf = bincache[key]
            # no admissible pose at all?
            lo = B['wlo'] / 2
            if cx1 < lo or cx0 > self.m - lo or cy1 < lo or cy0 > self.m - lo:
                stats['EMPTY'] += 1; leaves.append((box, 'EMPTY', None)); continue
            kind, wit = self.cert_core(box, B, Bf)
            if kind is None:
                kind, wit = self.cert_p1(box, B)
            if kind is None and self.use_tri:
                kind, wit = self.cert_tri(box)
            if kind is not None:
                stats[kind] += 1
                if self.dump: leaves.append((box, kind, wit))
                continue
            if depth >= self.max_depth:
                stats['UNCERT'] += 1; unc.append(box); continue
            # split the longest dimension (theta measured in radians ~ 2 du)
            dx, dy, du = cx1 - cx0, cy1 - cy0, 2 * (u1 - u0)
            if dx >= dy and dx >= du:
                mid = (cx0 + cx1) / 2
                stack.append(((cx0, mid, cy0, cy1, u0, u1), depth + 1))
                stack.append(((mid, cx1, cy0, cy1, u0, u1), depth + 1))
            elif dy >= du:
                mid = (cy0 + cy1) / 2
                stack.append(((cx0, cx1, cy0, mid, u0, u1), depth + 1))
                stack.append(((cx0, cx1, mid, cy1, u0, u1), depth + 1))
            else:
                mid = (u0 + u1) / 2
                stack.append(((cx0, cx1, cy0, cy1, u0, mid), depth + 1))
                stack.append(((cx0, cx1, cy0, cy1, mid, u1), depth + 1))
        return stats, unc, leaves

def _worker(args):
    chk, root = args
    return chk.run_box(root)

def roots(m, pitch=F(1, 10), ubins=8, cy_max=None):
    """Root boxes aligned to the pitch grid (so tight poses sit on box boundaries)."""
    m = F(m); cy_max = m / 2 if cy_max is None else F(cy_max)
    nx = int(m / pitch); ny = int(cy_max / pitch)
    R = []
    for i in range(nx):
        for j in range(ny):
            for k in range(ubins):
                R.append((i * pitch, (i + 1) * pitch, j * pitch, (j + 1) * pitch,
                          F(k, 2 * ubins), F(k + 1, 2 * ubins)))
    return R

FRIEDMAN14 = [(1, 1), (F(8, 5), 1), (F(12, 5), 1), (3, 1),
              (1, F(9, 5)), (2, F(9, 5)), (3, F(9, 5)),
              (1, F(11, 5)), (2, F(11, 5)), (3, F(11, 5)),
              (1, 3), (F(8, 5), 3), (F(12, 5), 3), (3, 3)]

def read_cert(path):
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    pts, ws = [], []
    for i in range(n):
        X, Y, w = map(int, tok[5 + 3 * i: 8 + 3 * i])
        pts.append((F(X, D), F(Y, D))); ws.append(F(w, W))
    return F(sn, sd), pts, ws

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what', help='friedman14 | cert')
    ap.add_argument('path', nargs='?')
    ap.add_argument('--tri', action='store_true')
    ap.add_argument('--depth', type=int, default=16)
    ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--pitch', type=str, default='1/10')
    ap.add_argument('--ubins', type=int, default=8)
    ap.add_argument('--dump', type=str, default=None)
    ap.add_argument('--full', action='store_true', help='no symmetry reduction (cy up to m, u up to 1)')
    a = ap.parse_args()
    if a.what == 'friedman14':
        m, pts, ws = 4, FRIEDMAN14, None
    else:
        m, pts, ws = read_cert(a.path)
    chk = Checker(m, pts, ws, use_tri=a.tri, max_depth=a.depth, dump=a.dump)
    sym = chk.symmetric()
    print(f"container [0,{m}]^2, {len(pts)} points, total weight {float(sum(chk.W)):.6f}, "
          f"symmetric under x->m-x and y->m-y: {sym}, triangles: {len(chk.tris) if a.tri else 'off'}")
    if not sym and not a.full:
        print("point set not symmetric: use --full"); sys.exit(2)
    pitch = F(a.pitch)
    if a.full:
        R = roots(m, pitch, a.ubins * 2, cy_max=m)
        R = [(x0, x1, y0, y1, u0 * 2, u1 * 2) for (x0, x1, y0, y1, u0, u1) in R]  # u in [0,1]: theta to 90deg
    else:
        R = roots(m, pitch, a.ubins)
    print(f"{len(R)} root boxes, depth limit {a.depth}, pitch {pitch}, u-bins {a.ubins}")
    t0 = time.time()
    tot = {'CORE': 0, 'P1': 0, 'TRI': 0, 'EMPTY': 0, 'UNCERT': 0, 'boxes': 0, 'maxdepth': 0}
    unc_all = []; leaves_all = []
    with Pool(a.nproc) as pool:
        for i, (st, unc, leaves) in enumerate(pool.imap_unordered(_worker, [(chk, r) for r in R], chunksize=4)):
            for k in tot:
                tot[k] = max(tot[k], st[k]) if k == 'maxdepth' else tot[k] + st[k]
            unc_all += unc; leaves_all += leaves
            if (i + 1) % 500 == 0:
                print(f"  {i+1}/{len(R)} roots, {tot['boxes']} boxes, uncert {tot['UNCERT']}, {time.time()-t0:.0f}s", flush=True)
    print(f"done in {time.time()-t0:.0f}s: boxes {tot['boxes']}, max depth {tot['maxdepth']}")
    print(f"  leaves: CORE {tot['CORE']}  P1 {tot['P1']}  TRI {tot['TRI']}  EMPTY {tot['EMPTY']}  UNCERTIFIED {tot['UNCERT']}")
    if unc_all:
        print("uncertified boxes (cx0 cx1 cy0 cy1 u0 u1 -> theta0 theta1 deg):")
        import math
        for b in sorted(unc_all)[:40]:
            print("  ", *[f"{float(v):.6f}" for v in b[:4]],
                  f"{math.degrees(2*math.atan(float(b[4]))):.3f} {math.degrees(2*math.atan(float(b[5]))):.3f}")
        if len(unc_all) > 40: print(f"  ... {len(unc_all)} in total")
    if a.dump:
        with open(a.dump, 'w') as f:
            f.write(f"# container {m}; box: cx0 cx1 cy0 cy1 u0 u1 ; kind ; witness\n")
            for box, kind, wit in leaves_all:
                f.write(" ".join(str(v) for v in box) + f" ; {kind} ; {wit}\n")
        print(f"leaves written to {a.dump}")
    print("VERIFIED" if tot['UNCERT'] == 0 else "NOT VERIFIED")

if __name__ == '__main__':
    main()
