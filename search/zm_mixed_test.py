#!/usr/bin/env python3
"""Tests for search/zm_mixed.py (search/ZM_MIXED.md sec 4).

  python3 search/zm_mixed_test.py selftest [--n 300]     randomised tests of Lemma S (lines, polygons), Lemma T,
                                                          the whole PIECE bound, and exact vs float mass
  python3 search/zm_mixed_test.py toys OUTDIR            write the m = 4 toy covers (valid one + rejection variants)
  python3 search/zm_mixed_test.py holes FILE [...]       float search for the worst poses of a cover (where a
                                                          rejection test must fail); prints exact mass there
  python3 search/zm_mixed.py stress FILE LEAFDUMP        re-check the certified leaves of a dump by sampling

Every check of the form "bound <= true mass" is made with EXACT (Fraction) masses at rational poses; the float
evaluator of mixed_cover.py is only compared against, never trusted.
"""
import sys, os, math, random, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np
import zeromargin as zm
import mixed_cover as MC
import zm_mixed as Z


def rnd_frac(lo, hi, den=10 ** 6):
    return F(random.randint(int(lo * den), int(hi * den)), den)


def rand_pose_in(box, m, tries=50, special=True):
    """a random ADMISSIBLE rational pose in the box (u rational), biased to the box's faces/corners and to the
    pivot poses of Lemma T; None if none found."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    for _ in range(tries):
        r = random.random()
        def pick(a, b):
            q = random.random()
            if special and q < 0.2: return a
            if special and q < 0.4: return b
            return a + (b - a) * F(random.randint(0, 10 ** 6), 10 ** 6)
        u = pick(u0, u1)
        if special and random.random() < 0.15 and u0 == 0:
            u = F(1, random.choice([10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]))
            if u > u1: u = u0
        cx = pick(cx0, cx1); cy = pick(cy0, cy1)
        if Z.admissible(m, cx, cy, u): return cx, cy, u
        # push into admissibility: clamp to [w/2, m - w/2]
        c, s = zm.trig(u); w = c + s
        cx = min(max(cx, w / 2), m - w / 2); cy = min(max(cy, w / 2), m - w / 2)
        if cx0 <= cx <= cx1 and cy0 <= cy <= cy1: return cx, cy, u
    return None


def pivot_poses(box, m, xi=None, eta=None, n=6):
    """poses of the box near the Lemma T pivots: theta tiny, centre shifted by ~theta^2/4 off the tile line."""
    out = []
    cx0, cx1, cy0, cy1, u0, u1 = box
    for _ in range(n):
        u = F(1, random.choice([50, 200, 1000, 10 ** 4]))
        if not (u0 <= u <= u1): u = u0 + (u1 - u0) * F(random.randint(0, 1000), 1000)
        c, s = zm.trig(u)
        cx = (xi + 1 / (2 * c)) if xi is not None else rnd_frac(float(cx0), float(cx1))
        cy = (eta + 1 / (2 * c)) if eta is not None else rnd_frac(float(cy0), float(cy1))
        cx += F(random.randint(-3, 3), 10 ** 6); cy += F(random.randint(-3, 3), 10 ** 6)
        if cx0 <= cx <= cx1 and cy0 <= cy <= cy1 and Z.admissible(m, cx, cy, u): out.append((cx, cy, u))
    return out


def seg_mass_exact(cov, cx, cy, u):
    """exact mass of the pieces only (segments + polygons)."""
    save = cov.points
    cov.points = []
    try:
        return Z.exact_mass(cov, cx, cy, u)
    finally:
        cov.points = save


# ---------------------------------------------------------------------------------------------- random covers
def random_mixed(m=4, D=1000, W=10 ** 6, nseg=40, npoly=4, npts=10, axis_frac=0.7, grid=True):
    S = m * D
    pts, segs, polys = [], [], []
    for _ in range(npts):
        pts.append((random.randint(0, S), random.randint(0, S), random.randint(1, W // 10)))
    for _ in range(nseg):
        r = random.random()
        if r < axis_frac:
            k = random.randint(0, m) * D if grid else random.randint(0, S)
            a = random.randint(0, S - 1); b = min(S, a + random.randint(1, D))
            if random.random() < 0.5: segs.append((k, a, k, b, random.randint(1, W // 10)))
            else: segs.append((a, k, b, k, random.randint(1, W // 10)))
        else:
            x0, y0 = random.randint(0, S), random.randint(0, S)
            x1 = min(S, max(0, x0 + random.randint(-D, D))); y1 = min(S, max(0, y0 + random.randint(-D, D)))
            if (x0, y0) != (x1, y1): segs.append((x0, y0, x1, y1, random.randint(1, W // 10)))
    for _ in range(npoly):
        cx, cy = random.randint(D // 2, S - D // 2), random.randint(D // 2, S - D // 2)
        k = random.randint(3, 6)
        ang = sorted(random.random() * 2 * math.pi for _ in range(k))
        rr = random.randint(D // 20, D)
        vs = [(min(S, max(0, int(cx + rr * math.cos(t)))), min(S, max(0, int(cy + rr * math.sin(t))))) for t in ang]
        H = Z.convex_hull(vs)
        if len(H) >= 3: polys.append((random.randint(1, W // 10), [tuple(map(int, v)) for v in H]))
    return dict(s_num=m, s_den=1, s=F(m), D=D, W=W, points=pts, segments=segs, polygons=polys)


def random_box(m, germ_frac=0.5, size=None):
    """a random rational pose box; with probability germ_frac a box touching a tile germ (theta bin at 0)."""
    size = size or random.choice([F(1, 10), F(1, 20), F(1, 80), F(1, 640)])
    if random.random() < germ_frac:
        xi = random.randint(0, int(m) - 1); eta = random.randint(0, int(m) - 1)
        cx0 = xi + HALF - size * random.randint(0, 1); cy0 = eta + HALF - size * random.randint(0, 1)
        u0 = F(0); u1 = random.choice([F(1, 16), F(1, 128), F(1, 1024), F(1, 8192)])
    else:
        cx0 = F(random.randint(0, int(m / size) - 1)) * size
        cy0 = F(random.randint(0, int(m / size) - 1)) * size
        k = random.randint(0, 7); u0 = F(k, 16); u1 = F(k + 1, 16)
        if random.random() < 0.5:
            u1 = u0 + (u1 - u0) / random.choice([1, 4, 32])
    return (cx0, cx0 + size, cy0, cy0 + size, u0, u1)


HALF = F(1, 2)


def selftest(n=300, seed=7):
    random.seed(seed)
    fails = 0
    t0 = time.time()
    # ---- (1) exact vs float mass (also a check of mixed_cover.mass_in_square_float)
    worst = 0.0
    for it in range(40):
        cv = random_mixed()
        cov = Z.Cover(cv)
        for _ in range(10):
            box = random_box(cov.m, germ_frac=0.2)
            p = rand_pose_in(box, cov.m)
            if p is None: continue
            cx, cy, u = p
            ex = float(Z.exact_mass(cov, cx, cy, u))
            fl = MC.mass_in_square_float(cv, float(cx), float(cy), 2 * math.atan(float(u)))
            worst = max(worst, abs(ex - fl))
    ok = worst < 1e-6
    print(f"[1] exact vs mixed_cover float mass, 400 poses: max |diff| = {worst:.3e}  {'ok' if ok else 'FAIL'}")
    fails += (not ok)
    # (the float evaluator has a tolerance 1e-12 and a single piece on an edge can flip; the random poses are
    #  generic, so the two agree to rounding)

    # ---- (2) Lemma S (lines): every sampled point of J4 lies in Q at sampled admissible poses; every cond too
    nchk = nbad = 0
    chk_dummy = None
    for it in range(n):
        cv = random_mixed(nseg=6, npoly=0, npts=0)
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        box = random_box(cov.m)
        cu1 = zm.clip_bin(box, cov.m)
        box = box[:5] + (cu1,)
        B = zm.bin_data(box[4], box[5])
        specs = mc.zc._adm_specs(box, B)
        for L in cov.lines.values():
            I = [Z.line_cond_iv(L['P0'], L['d'], specs, c, box[4], box[5] - box[4], cov.m) for c in range(4)]
            for c in range(5):
                J = I[c] if c < 4 else Z.iv_and(Z.iv_and(I[0], I[1]), Z.iv_and(I[2], I[3]))
                if J is None: continue
                lo = J[0] if J[0] is not None else F(-3); hi = J[1] if J[1] is not None else F(8)
                if lo > hi: continue
                for _ in range(4):
                    t = lo + (hi - lo) * F(random.choice([0, 1000, random.randint(0, 1000)]), 1000)
                    px = L['P0'][0] + t * L['d'][0]; py = L['P0'][1] + t * L['d'][1]
                    for _ in range(3):
                        p = rand_pose_in(box, cov.m)
                        if p is None: continue
                        cx, cy, u = p
                        cs, sn = zm.trig(u)
                        X = (px - cx) * cs + (py - cy) * sn; Y = -(px - cx) * sn + (py - cy) * cs
                        okc = [X <= HALF, X >= -HALF, Y <= HALF, Y >= -HALF]
                        good = okc[c] if c < 4 else all(okc)
                        nchk += 1
                        if not good:
                            nbad += 1
                            if nbad <= 3: print("   Lemma S violation:", box, L['key'], c, t, p)
    print(f"[2] Lemma S (lines): {nchk} (box, line, t, pose) checks, {nbad} violations  {'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (3) Lemma S (polygons): K's vertices in Q; polygon bound <= exact polygon mass
    nchk = nbad = 0
    for it in range(max(n // 3, 20)):
        cv = random_mixed(nseg=0, npoly=6, npts=0)
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        box = random_box(cov.m, germ_frac=0.3)
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        L, _ = mc.piece_bound(box, B)
        for _ in range(8):
            p = rand_pose_in(box, cov.m)
            if p is None: continue
            v = seg_mass_exact(cov, *p)
            nchk += 1
            if L > v:
                nbad += 1
                if nbad <= 3: print("   polygon bound violation:", box, p, float(L), float(v))
    print(f"[3] Lemma S (polygons): {nchk} (box, pose) checks of bound <= mass, {nbad} violations  "
          f"{'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (4) Lemma T and the whole PIECE bound at germ boxes: bound <= exact mass at poses incl. pivots
    nchk = nbad = 0; ngain = 0
    for it in range(n):
        m = 4
        cv = random_mixed(nseg=0, npoly=0, npts=0)
        # dense grid-line segments near one tile, random densities, random gaps
        D, W = cv['D'], cv['W']
        xi, eta = random.randint(0, 3), random.randint(0, 3)
        segs = []
        q = random.choice([4, 10, 20])
        for line in (xi, xi + 1):
            for k in range(-q, 2 * q):
                if random.random() < 0.15: continue
                a = eta * D + k * D // q; b = a + D // q
                if a < 0 or b > m * D: continue
                segs.append((line * D, a, line * D, b, random.randint(1, W // 20)))
        for line in (eta, eta + 1):
            for k in range(-q, 2 * q):
                if random.random() < 0.15: continue
                a = xi * D + k * D // q; b = a + D // q
                if a < 0 or b > m * D: continue
                segs.append((a, line * D, b, line * D, random.randint(1, W // 20)))
        cv['segments'] = segs
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        size = random.choice([F(1, 10), F(1, 40), F(1, 320)])
        sx, sy = random.randint(-1, 0), random.randint(-1, 0)
        box = (xi + HALF + sx * size, xi + HALF + (sx + 1) * size, eta + HALF + sy * size,
               eta + HALF + (sy + 1) * size, F(0), random.choice([F(1, 16), F(1, 256), F(1, 4096)]))
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        Lb, thr = mc.piece_bound(box, B)
        ngain += (thr == 'T')
        poses = [p for p in (rand_pose_in(box, cov.m) for _ in range(12)) if p]
        poses += pivot_poses(box, cov.m, xi=xi, n=4) + pivot_poses(box, cov.m, eta=eta, n=4)
        poses += pivot_poses(box, cov.m, xi=xi, eta=eta, n=4)
        for p in poses:
            v = seg_mass_exact(cov, *p)
            nchk += 1
            if Lb > v:
                nbad += 1
                if nbad <= 3: print("   Lemma T / PIECE violation:", box, p, float(Lb), float(v))
    print(f"[4] Lemma T at tile germs: {nchk} (box, pose) checks of bound <= exact mass, {nbad} violations, "
          f"Lemma T raised the bound in {ngain}/{n} boxes  {'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (4c) Lemma T with line points: bound(with_pts) <= exact mass of pieces + the moved points
    nchk = nbad = 0; nmv = 0
    for it in range(n):
        m = 4; D, W = 1000, 10 ** 6
        xi, eta = random.randint(0, 3), random.randint(0, 3)
        q = random.choice([4, 10])
        segs = []; pts = []
        for line in (xi, xi + 1):
            for k in range(-q, 2 * q):
                a = eta * D + k * D // q; b = a + D // q
                if a < 0 or b > m * D: continue
                if random.random() < 0.7: segs.append((line * D, a, line * D, b, random.randint(1, W // 40)))
                for _ in range(random.randint(0, 3)):
                    pts.append((line * D, random.randint(a, b), random.randint(1, W // 60)))
        for line in (eta, eta + 1):
            for k in range(-q, 2 * q):
                a = xi * D + k * D // q; b = a + D // q
                if a < 0 or b > m * D: continue
                if random.random() < 0.7: segs.append((a, line * D, b, line * D, random.randint(1, W // 40)))
                for _ in range(random.randint(0, 3)):
                    pts.append((random.randint(a, b), line * D, random.randint(1, W // 60)))
        pts += [(xi * D, eta * D, 5000), ((xi + 1) * D, (eta + 1) * D, 5000)]      # grid crossings
        cv = dict(s_num=m, s_den=1, s=F(m), D=D, W=W, points=pts, segments=segs, polygons=[])
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        size = random.choice([F(1, 10), F(1, 40), F(1, 320)])
        sx, sy = random.randint(-1, 0), random.randint(-1, 0)
        box = (xi + HALF + sx * size, xi + HALF + (sx + 1) * size, eta + HALF + sy * size,
               eta + HALF + (sy + 1) * size, F(0), random.choice([F(1, 16), F(1, 256), F(1, 4096)]))
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        moved = []
        Lb, _ = mc.piece_bound(box, B, with_pts=True, moved=moved)
        assert len(set(moved)) == len(moved)
        nmv += len(moved)
        poses = [p for p in (rand_pose_in(box, cov.m) for _ in range(12)) if p]
        poses += pivot_poses(box, cov.m, xi=xi, n=4) + pivot_poses(box, cov.m, eta=eta, n=4)
        poses += pivot_poses(box, cov.m, xi=xi, eta=eta, n=4)
        for p in poses:
            cxp, cyp, up = p
            cs, sn = zm.trig(up)
            v = seg_mass_exact(cov, *p) + sum((cov.points[k][2] for k in moved
                                               if zm.in_rot_square(cov.points[k][0] - cxp, cov.points[k][1] - cyp, cs, sn)), F(0))
            nchk += 1
            if Lb > v:
                nbad += 1
                if nbad <= 3: print("   Lemma T+points violation:", box, p, float(Lb), float(v))
    print(f"[4c] Lemma T with line points: {nchk} checks, {nbad} violations ({nmv} points moved)  "
          f"{'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (4b) Lemma L: dense axis-line covers (random densities), random small boxes anywhere, bound <= mass
    nchk = nbad = 0; nl = 0
    for it in range(n):
        m = 4
        D, W = 1000, 10 ** 6
        q = random.choice([10, 20, 50])
        segs = []
        for line in range(0, m + 1):
            for k in range(m * q):
                if random.random() < 0.05: continue
                a = k * D // q; b = a + D // q
                segs.append((line * D, a, line * D, b, random.randint(1, W // (2 * q))))
                segs.append((a, line * D, b, line * D, random.randint(1, W // (2 * q))))
        cv = dict(s_num=m, s_den=1, s=F(m), D=D, W=W, points=[], segments=segs, polygons=[])
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        box = random_box(cov.m, germ_frac=0.2, size=random.choice([F(1, 40), F(1, 160), F(1, 640)]))
        if it % 2:        # a wall / corner box: centre within the admissible band of a wall
            sz = box[1] - box[0]
            wx = random.choice([F(1, 2), F(7072, 10000) - sz, F(3, 5)]); wy = random.choice([wx, F(3, 2), F(2)])
            if random.random() < .5: wx = cov.m - wx - sz
            box = (wx, wx + sz, wy, wy + sz, box[4], box[5])
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        Lb, thr = mc.piece_bound(box, B)
        nl += (thr == 'L')
        for _ in range(10):
            p = rand_pose_in(box, cov.m)
            if p is None: continue
            v = seg_mass_exact(cov, *p)
            nchk += 1
            if Lb > v:
                nbad += 1
                if nbad <= 3: print("   Lemma L / PIECE violation:", box, p, float(Lb), float(v))
    print(f"[4b] Lemma L, dense axis covers: {nchk} (box, pose) checks, {nbad} violations, Lemma L used in "
          f"{nl}/{n} boxes  {'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (4d) Lemma L / L' on near-tight grid covers (T1-like, densities jittered per segment), random small boxes
    nchk = nbad = 0; nl = 0
    base_cv = grid_cover(rho=F(60507, 100000))
    for it in range(n):
        cv = dict(base_cv)
        cv['segments'] = [sg[:4] + (sg[4] + random.randint(-sg[4] // 20, sg[4] // 20),) for sg in base_cv['segments']]
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        sz = random.choice([F(1, 80), F(1, 320), F(1, 1280)])
        cx0 = F(random.randint(0, int(2 / sz) - 1)) * sz; cy0 = F(random.randint(0, int(2 / sz) - 1)) * sz
        k = random.randint(0, 63); u0 = F(k, 128); u1 = F(k + 1, 128)
        box = (cx0, cx0 + sz, cy0, cy0 + sz, u0, u1)
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        lo = B['wlo'] / 2
        if box[1] < lo or box[3] < lo: continue
        Lb, thr = mc.piece_bound(box, B)
        nl += (thr == 'L')
        for _ in range(8):
            p = rand_pose_in(box, cov.m)
            if p is None: continue
            v = seg_mass_exact(cov, *p)
            nchk += 1
            if Lb > v:
                nbad += 1
                if nbad <= 3: print("   Lemma L violation (near-tight):", box, p, float(Lb), float(v))
    print(f"[4d] Lemma L on near-tight grid covers: {nchk} checks, {nbad} violations, Lemma L used in {nl}/{n} boxes  "
          f"{'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (4e) PIECE bound on a real LP cover (agent A's r7, if present): random small boxes, bound <= exact mass
    for rp in [os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', f)
               for f in ('lc_m4mixA_r7.txt', 'line-cover_m5_candidate.txt')]:
      if os.path.exists(rp):
          cvr = MC.load(rp); cvr = dict(cvr, points=[])
          cov = Z.Cover(cvr); mc = Z.MixedChecker(cov)
          nchk = nbad = 0; nl = 0
          vcount = [0]; _vs = Z.vertex_split
          def _vcount(*a, **k):
              r = _vs(*a, **k)
              if r is not None: vcount[0] += 1
              return r
          Z.vertex_split = _vcount
          for it in range(n):
              sz = random.choice([F(1, 80), F(1, 320), F(1, 1280)])
              cx0 = F(random.randint(0, int(2 / sz) - 1)) * sz; cy0 = F(random.randint(0, int(2 / sz) - 1)) * sz
              k = random.randint(0, 63); u0 = F(k, 128); u1 = u0 + F(1, 128) / random.choice([1, 4, 16])
              if it % 3 == 0:        # the m5 wall / near-vertex area of ZM_MIXED.md sec 7.4
                  cx0 = F(random.randint(3300, 3580), 6400); cy0 = F(random.randint(9200, 9470), 6400); sz = F(1, 640)
                  u0 = F(random.randint(200, 640), 10240); u1 = u0 + F(1, 10240 * random.choice([1, 4, 16]))
              box = (cx0, cx0 + sz, cy0, cy0 + sz, u0, u1)
              box = box[:5] + (zm.clip_bin(box, cov.m),)
              B = zm.bin_data(box[4], box[5])
              if box[1] < B['wlo'] / 2 or box[3] < B['wlo'] / 2: continue
              Lb, thr = mc.piece_bound(box, B)
              nl += (thr == 'L')
              for _ in range(8):
                  p = rand_pose_in(box, cov.m)
                  if p is None: continue
                  v = seg_mass_exact(cov, *p)
                  nchk += 1
                  if Lb > v:
                      nbad += 1
                      if nbad <= 3: print("   PIECE violation (r7):", box, p, float(Lb), float(v))
          Z.vertex_split = _vs
          print(f"[4e] PIECE bound on {os.path.basename(rp)} segments: {nchk} checks, {nbad} violations, Lemma L used in {nl}/{n} boxes, Lemma V evaluated {vcount[0]} times  "
                f"{'ok' if nbad == 0 else 'FAIL'}")
          fails += (nbad > 0)

    # ---- (6) Lemma R / SPLIT: per region, every claimed point is in Q and the region piece bound <= piece mass
    nchk = nbad = 0; nreg = 0; nok = 0
    base_cv = grid_cover(rho=F(45, 100))
    for it in range(n):
        cv = dict(base_cv)
        cv['segments'] = [sg[:4] + (sg[4] + random.randint(-sg[4] // 10, sg[4] // 10),) for sg in base_cv['segments']]
        cx0 = F(random.randint(10, 150), 100); cy0 = F(random.randint(10, 150), 100)
        sz = random.choice([F(1, 40), F(1, 160), F(1, 640)])
        pts = []
        for _ in range(random.randint(20, 80)):      # points around the box, many near the reach of its squares
            ang = random.random() * 2 * math.pi; rr = random.uniform(0.3, 0.75)
            X = int(1000 * (float(cx0) + rr * math.cos(ang))); Y = int(1000 * (float(cy0) + rr * math.sin(ang)))
            if 0 <= X <= 4000 and 0 <= Y <= 4000: pts.append((X, Y, random.randint(1000, 60000)))
        cv['points'] = pts
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov, use_chain=True)
        k = random.randint(1, 60); u0 = F(k, 128); u1 = u0 + F(1, 128) / random.choice([1, 8, 64])
        box = (cx0, cx0 + sz, cy0, cy0 + sz, u0, u1)
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        if box[1] < B['wlo'] / 2 or box[3] < B['wlo'] / 2: continue
        Lb, _ = mc.piece_bound(box, B)
        diag = []
        mc.cert_split(box, B, None, Lb, mc._lparts, diag=diag)
        for dg in diag:
            nok += dg['ok']
            chain = dg['chain']; kd = dg['kind']
            poses = [rand_pose_in(box, cov.m) for _ in range(6)]
            for _ in range(6):             # poses ON (and just off) a pivot surface G_q = 0: solve for cx
                if not chain: break
                q = random.choice(chain); p0 = rand_pose_in(box, cov.m)
                if p0 is None: continue
                _, cyp, up = p0
                # G = g0 + g1 u + g2 u^2 with a = px - cx affine: G(cx) = G(0) + cx * dG/dcx
                def Gx(cx):
                    g = zm.Checker._gcoef(kd, cov.points[q][0] - cx, cov.points[q][1] - cyp)
                    return g[0] + g[1] * up + g[2] * up * up
                s0 = Gx(F(0)); s1 = Gx(F(1)) - s0
                if s1 == 0: continue
                cxs = -s0 / s1 + random.choice([0, 0, F(1, 10 ** 9), -F(1, 10 ** 9)])
                if box[0] <= cxs <= box[1] and Z.admissible(cov.m, cxs, cyp, up): poses.append((cxs, cyp, up))
            for p in poses:
                if p is None: continue
                cxp, cyp, up = p
                cs, sn = zm.trig(up)
                def G(k):
                    g = zm.Checker._gcoef(kd, cov.points[k][0] - cxp, cov.points[k][1] - cyp)
                    return g[0] + g[1] * up + g[2] * up * up
                r = 0
                for j, q in enumerate(chain):
                    if G(q) <= 0: r = j + 1
                reg = dg['regions'][r]
                nreg += 1
                inq = all(zm.in_rot_square(cov.points[k][0] - cxp, cov.points[k][1] - cyp, cs, sn) for k in reg['pts'])
                pm = seg_mass_exact(cov, *p)
                nchk += 1
                if not inq or reg['pieces'] > pm or reg['empty']:
                    nbad += 1
                    if nbad <= 3: print("   Lemma R violation:", box, p, r, inq, float(reg['pieces']), float(pm), reg['empty'])
    print(f"[6] Lemma R (SPLIT regions): {nchk} (box, pose) checks, {nbad} violations; {nok} chains fully certified  "
          f"{'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)

    # ---- (5) the whole PIECE bound on random mixed covers and random boxes
    nchk = nbad = 0
    for it in range(n):
        cv = random_mixed(nseg=30, npoly=3, npts=0)
        cov = Z.Cover(cv)
        mc = Z.MixedChecker(cov)
        box = random_box(cov.m, germ_frac=0.4)
        box = box[:5] + (zm.clip_bin(box, cov.m),)
        B = zm.bin_data(box[4], box[5])
        Lb, _ = mc.piece_bound(box, B)
        for _ in range(6):
            p = rand_pose_in(box, cov.m)
            if p is None: continue
            v = seg_mass_exact(cov, *p)
            nchk += 1
            if Lb > v:
                nbad += 1
                if nbad <= 3: print("   PIECE violation:", box, p, float(Lb), float(v))
    print(f"[5] PIECE bound, random covers/boxes: {nchk} checks, {nbad} violations  {'ok' if nbad == 0 else 'FAIL'}")
    fails += (nbad > 0)
    print(f"selftest: {'PASS' if fails == 0 else 'FAIL (%d)' % fails} in {time.time()-t0:.0f}s")
    return 1 if fails else 0


# ---------------------------------------------------------------------------------------------- leaf stress
def stress(cov, dump, per=20, seed=1):
    """For every certified leaf of a zm_mixed dump: `per` random admissible rational poses (faces, corners and
    pivot poses favoured); float mass by mixed_cover.MixedEval, exact mass (Fractions) for any float value
    < 1 + 1e-6.  A certified leaf with an exact mass < 1 would be a soundness failure."""
    random.seed(seed)
    ev = MC.MixedEval.from_cover(cov.raw)
    covs = {'-': cov, 'swap:': None}
    nleaf = npose = 0; worst = (9.0, None); bad = 0
    for line in open(dump):
        if line.startswith('#'): continue
        parts = [p.strip() for p in line.split(';')]
        lab, boxs, kind = parts[0], parts[1], parts[2]
        if kind in ('UNCERT', 'EMPTY'): continue
        if lab != '-':
            if covs['swap:'] is None: covs['swap:'] = Z.Cover(Z.swapped(cov.raw))
        cc = covs[lab]
        e = ev if lab == '-' else MC.MixedEval.from_cover(cc.raw)
        box = tuple(F(t) for t in boxs.split())
        nleaf += 1
        poses = [p for p in (rand_pose_in(box, cc.m) for _ in range(per)) if p]
        if box[4] == 0:
            for xi in range(int(cc.m)):
                if box[0] <= xi + HALF + F(1, 100) and box[1] >= xi + HALF - F(1, 100):
                    poses += pivot_poses(box, cc.m, xi=xi, n=3)
            for eta in range(int(cc.m)):
                if box[2] <= eta + HALF + F(1, 100) and box[3] >= eta + HALF - F(1, 100):
                    poses += pivot_poses(box, cc.m, eta=eta, n=3)
        if not poses: continue
        cx = np.array([float(p[0]) for p in poses]); cy = np.array([float(p[1]) for p in poses])
        th = np.array([2 * math.atan(float(p[2])) for p in poses])
        fm = e.mass(cx, cy, th)
        for p, v in zip(poses, fm):
            npose += 1
            if v < 1 + 1e-6:
                ex = Z.exact_mass(cc, *p)
                if ex < worst[0]: worst = (ex, (lab, p, box, kind))
                if ex < 1:
                    bad += 1
                    print("*** SOUNDNESS FAILURE: certified leaf", lab, box, kind, "pose", p, "exact", ex)
            elif v < worst[0]:
                worst = (v, (lab, p, box, kind))
    print(f"stress: {nleaf} certified leaves, {npose} poses, min mass {float(worst[0]):.9f} at {worst[1]}; "
          f"{bad} exact violations  {'ok' if bad == 0 else 'FAIL'}")
    return 1 if bad else 0


# ---------------------------------------------------------------------------------------------- toy covers
def grid_cover(m=4, rho=F(5, 8), q=10, D=1000, W=10 ** 6, pts=()):
    """the interior grid lines x, y = 1..m-1 at uniform density rho per unit length, cut into segments of
    length 1/q (all of equal mass rho/q)."""
    segs = []
    wseg = rho / q * W
    assert wseg.denominator == 1
    wseg = int(wseg)
    for k in range(1, m):
        for j in range(m * q):
            a, b = j * D // q, (j + 1) * D // q
            segs.append((k * D, a, k * D, b, wseg))
            segs.append((a, k * D, b, k * D, wseg))
    return dict(s_num=m, s_den=1, s=F(m), D=D, W=W, points=list(pts), segments=segs, polygons=[])


def toys(outdir):
    os.makedirs(outdir, exist_ok=True)
    out = {}
    D, W, q = 1000, 10 ** 6, 10
    # T1: uniform grid lines at rho = 5/8 (total 15): the minimum over poses is at theta near 45 deg (vertical +
    # horizontal chords 2 * 0.8284 * rho = 1.0355), every tile germ has 2 rho = 1.25 one-sided.
    T1 = grid_cover(rho=F(5, 8)); out['T1_grid'] = T1
    # T1 with a centre point on every tile centre (a mixed point + line cover; D4-symmetric)
    pts = [(D // 2 + i * D, D // 2 + j * D, 20000) for i in range(4) for j in range(4)]
    out['T1p_grid_pts'] = grid_cover(rho=F(5, 8), pts=pts)
    # R_lighten: all weights x 0.95 (min 0.984 < 1)
    Tl = grid_cover(rho=F(5, 8) * F(95, 100)); out['R_lighten'] = Tl
    # R_remove: T1 without the four symmetric images of one 1/10 segment on x = 2 (y in [1.2, 1.3])
    def drop(cv, segset):
        c = dict(cv); c['segments'] = [s for s in cv['segments'] if tuple(s[:4]) not in segset]; return c
    S = 4 * D
    def d4set(X0, Y0, X1, Y1):
        imgs = set()
        for g in (lambda x, y: (x, y), lambda x, y: (S - x, y), lambda x, y: (x, S - y), lambda x, y: (S - x, S - y),
                  lambda x, y: (y, x), lambda x, y: (S - y, x), lambda x, y: (y, S - x), lambda x, y: (S - y, S - x)):
            a, b = g(X0, Y0); c, d = g(X1, Y1)
            imgs.add((a, b, c, d)); imgs.add((c, d, a, b))
        return imgs
    out['R_remove'] = drop(T1, d4set(2000, 1200, 2000, 1300))
    # R_pivot (germ-specific): rho = 0.81 and at the tile [1,2]x[1,2] the pieces x = 1, y in [1.5, 1.9] and
    # x = 2, y in [1.1, 1.5] removed (D4 images too).  At theta = 0 every square keeps >= 1.296; at theta = 0+
    # with the pivot T = 1.5 (centre (1 + 1/(2 cos), 1.5)) the pair x = 1 / x = 2 loses 0.8 of its unit.
    Tp = grid_cover(rho=F(81, 100))
    rem = set()
    for (y0, y1, x) in ((1500, 1900, 1000), (1100, 1500, 2000)):
        for j in range((y1 - y0) * q // D):
            rem |= d4set(x, y0 + j * D // q, x, y0 + (j + 1) * D // q)
    out['R_pivot'] = drop(Tp, rem)
    # R_shift: the same pieces SHIFTED off the grid lines by 1/100 instead of removed: x = 1 -> 0.99, and
    # x = 2 -> 2.01 and 1.99 at half weight each (so that the cover stays D4-symmetric); mass unchanged
    Ts = drop(Tp, rem)
    wseg = int(F(81, 100) / q * W)
    canon = {}
    for (y0, y1, x, shifts) in ((1500, 1900, 1000, ((990, 1),)), (1100, 1500, 2000, ((2010, 2), (1990, 2)))):
        for j in range((y1 - y0) * q // D):
            for x2, div in shifts:
                for a, b, c, d in d4set(x2, y0 + j * D // q, x2, y0 + (j + 1) * D // q):
                    canon[min((a, b, c, d), (c, d, a, b))] = wseg // div
    Ts = dict(Ts); Ts['segments'] = Ts['segments'] + [k + (w,) for k, w in sorted(canon.items())]
    out['R_shift'] = Ts
    # P_T1: the POINT analogue of T1 (the same lines as points at spacing 1/100, weight rho/100 each), for the
    # performance comparison with zeromargin.py (a plain-format file: ns = npg = 0)
    pp = []
    for k in range(1, 4):
        for j in range(0, 401):
            pp.append((k * D, j * 10, 6250)); pp.append((j * 10, k * D, 6250))
    out['P_T1pts'] = dict(s_num=4, s_den=1, s=F(4), D=D, W=W, points=pp, segments=[], polygons=[])
    # the valid version of R_pivot's density: rho = 0.81 uniform (sanity: must certify)
    out['T2_grid81'] = Tp
    for name, cv in out.items():
        MC.validate(cv)
        path = os.path.join(outdir, f'zmm_{name}.txt')
        MC.write(path, cv, comment=f'zm_mixed toy cover {name} (search/zm_mixed_test.py toys)')
        cov = Z.Cover(MC.load(path))
        print(f"{path}: {len(cv['points'])} points, {len(cv['segments'])} segments, total {float(cov.total):.6f}, "
              f"D4 {cov.symmetric_d4()}")


def holes(path, pitch=0.01, thetas=None):
    """float scan for the minimum of mu over poses (theta grid incl. tiny angles), then exact mass at the
    argmin snapped to rationals.  Heuristic locator for where a rejection must happen."""
    cv = MC.load(path); cov = Z.Cover(cv)
    ev = MC.MixedEval.from_cover(cv)
    m = float(cov.m)
    if thetas is None:
        thetas = [0.0, 1e-5, 1e-4, 1e-3, 0.003, 0.01, 0.03] + list(np.radians(np.arange(1, 46, 1.0)))
    best = []
    for th in thetas:
        w = math.cos(th) + math.sin(th)
        xs = np.arange(w / 2, m - w / 2 + 1e-12, pitch)
        X, Y = np.meshgrid(xs, xs); X = X.ravel(); Y = Y.ravel()
        v = ev.mass(X, Y, np.full(len(X), th))
        i = int(np.argmin(v))
        best.append((float(v[i]), float(X[i]), float(Y[i]), th))
    # pivot family around every tile: centre (xi + 1/(2cos), eta + t), theta tiny
    for th in (1e-4, 1e-3, 1e-2):
        c = math.cos(th)
        for xi in range(int(m)):
            for eta in range(int(m)):
                ts = np.linspace(0.3, 0.7, 81)
                X = np.full(len(ts), xi + 1 / (2 * c)); Y = eta + ts
                w = c + math.sin(th)
                ok = (X >= w / 2) & (X <= m - w / 2) & (Y >= w / 2) & (Y <= m - w / 2)
                if not ok.any(): continue
                X, Y = X[ok], Y[ok]
                v = ev.mass(X, Y, np.full(len(X), th)); i = int(np.argmin(v))
                best.append((float(v[i]), float(X[i]), float(Y[i]), th))
    best.sort()
    print(f"{path}: float minimum over the scan {best[0][0]:.6f} at cx={best[0][1]:.6f} cy={best[0][2]:.6f} "
          f"theta={math.degrees(best[0][3]):.5f} deg")
    for v, x, y, th in best[:5]:
        u = F(math.tan(th / 2)).limit_denominator(10 ** 7); cx = F(x).limit_denominator(10 ** 6)
        cy = F(y).limit_denominator(10 ** 6)
        if not Z.admissible(cov.m, cx, cy, u): continue
        ex = Z.exact_mass(cov, cx, cy, u)
        print(f"   float {v:.6f}  exact {float(ex):.9f}  at cx={cx} cy={cy} u={u} ({math.degrees(th):.5f} deg)")
    return best


def locate(dump, cx, cy, u):
    """which leaves of a dump contain the pose (cx, cy, u)?  A rejection test passes when every box containing
    a pose of exact mass < 1 is UNCERT (a certified one would be a soundness failure)."""
    hits = []
    for line in open(dump):
        if line.startswith('#'): continue
        parts = [p.strip() for p in line.split(';')]
        b = tuple(F(t) for t in parts[1].split())
        if b[0] <= cx <= b[1] and b[2] <= cy <= b[3] and b[4] <= u <= b[5]: hits.append((parts[2], b))
    return hits


def unc_diag(path, dump, per=3000, seed=3):
    """float sampling of every UNCERT box of a dump: is the COVER short there (min < 1: a genuine hole, the
    rejection is correct) or is the checker short (min >= 1)?  Heuristic."""
    rng = np.random.default_rng(seed)
    cv = MC.load(path); cov = Z.Cover(cv); ev = MC.MixedEval.from_cover(cv)
    m = float(cov.m); res = []
    for line in open(dump):
        parts = [p.strip() for p in line.split(';')]
        if line.startswith('#') or parts[2] != 'UNCERT': continue
        b = [float(F(t)) for t in parts[1].split()]
        u = rng.uniform(b[4], b[5], per); u[:per // 10] = b[4]; u[per // 10: per // 5] = b[5]
        th = 2 * np.arctan(u); w = np.cos(th) + np.sin(th)
        cx = rng.uniform(b[0], b[1], per); cy = rng.uniform(b[2], b[3], per)
        cx = np.clip(cx, np.maximum(b[0], w / 2), np.minimum(b[1], m - w / 2))
        cy = np.clip(cy, np.maximum(b[2], w / 2), np.minimum(b[3], m - w / 2))
        ok = (cx >= w / 2 - 1e-15) & (cx <= m - w / 2 + 1e-15) & (cy >= w / 2 - 1e-15) & (cy <= m - w / 2 + 1e-15)
        if not ok.any(): continue
        v = ev.mass(cx[ok], cy[ok], th[ok]); i = int(np.argmin(v))
        res.append((float(v[i]), b, (cx[ok][i], cy[ok][i], math.degrees(th[ok][i]))))
    res.sort(key=lambda r: r[0])
    nh = sum(1 for r in res if r[0] < 1)
    print(f"{len(res)} uncertified boxes: {nh} with a sampled pose of mass < 1 (cover short), "
          f"{len(res) - nh} with sampled min >= 1 (checker short); overall min {res[0][0] if res else None}")
    for r in res[:5]: print("  ", f"{r[0]:.6f}", [round(x, 6) for x in r[1][:4]], r[2])
    if res and res[-1][0] >= 1: print("   largest sampled min among them:", f"{res[-1][0]:.6f}", [round(x, 6) for x in res[-1][1][:4]])
    return res


def scale_cover(path, out, f):
    """exact rational scaling of every mass of a mixed file by f (a Fraction): W -> W * den, w -> w * num."""
    f = F(f); cv = MC.load(path); a, b = f.numerator, f.denominator
    cv2 = dict(cv, W=cv['W'] * b, points=[(x, y, w * a) for x, y, w in cv['points']],
               segments=[sg[:4] + (sg[4] * a,) for sg in cv['segments']], polygons=[(w * a, v) for w, v in cv['polygons']])
    MC.write(out, cv2, comment=f'{path} x {f} (checker test)')
    return MC.total(cv2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what')
    ap.add_argument('args', nargs='*')
    ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()
    if a.what == 'selftest': return selftest(a.n, a.seed)
    if a.what == 'toys': return toys(a.args[0])
    if a.what == 'locate':
        dump, cx, cy, u = a.args[0], F(a.args[1]), F(a.args[2]), F(a.args[3])
        hits = locate(dump, cx, cy, u)
        for k, b in hits: print(k, [str(v) for v in b])
        cert = [k for k, _ in hits if k not in ('UNCERT', 'EMPTY')]
        print(f"pose in {len(hits)} leaves; certified leaves containing it: {len(cert)} "
              f"({'REJECTED here: ok' if hits and not cert else 'NOT rejected here'})")
        return 0 if hits and not cert else 1
    if a.what == 'scale':
        print(float(scale_cover(a.args[0], a.args[1], F(a.args[2])))); return 0
    if a.what == 'unc':
        unc_diag(a.args[0], a.args[1]); return 0
    if a.what == 'holes':
        for p in a.args: holes(p)
        return 0
    print(__doc__); return 2


if __name__ == '__main__':
    sys.exit(main() or 0)
