#!/usr/bin/env python3
"""qx2_zm_test.py -- tests of the new exact primitives in qx2_zm.py (task quadrant-exact-w2).

  opts            the chord-end option formulas against a direct exact chord computation
  poly            Bernstein positivity against dense exact evaluation
  exact COVER     Lemma E on random boxes: (i) soundness at tau = sampled minimum + delta (must never certify
                  when some sampled pose is below tau), (ii) completeness at tau = sampled min - delta.
  capleb COVER    Lemmas U and K against exact masses at random admissible poses of certified boxes
"""
import sys, os, math, random, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import qx2_zm as X
import zm_mixed as ZM
import zeromargin as zm
import mixed_cover as MC
import qx2_exact as E


def t_opts(n=3000, seed=1):
    rng = random.Random(seed); bad = 0
    for _ in range(n):
        cx = F(rng.randint(0, 4000), 1000); cy = F(rng.randint(0, 4000), 1000)
        u = F(rng.randint(1, 499), 1000)
        C, S = E.trig(u)
        for o in ('H', 'V'):
            p = F(rng.randint(0, 4000), 1000)
            if o == 'H': P0, d = (F(0), p), (F(1), F(0))
            else: P0, d = (p, F(0)), (F(0), F(1))
            iv = E.chord_iv(P0, d, cx, cy, C, S)
            opts = X.options(o, p)
            vals = []
            for (typ, Fk, a0, ax, ay) in opts:
                Fv = X.peval(X.PC if Fk == 'C' else X.PS, u)
                vals.append((typ, (X.peval(a0, u) + X.peval(ax, u) * cx + X.peval(ay, u) * cy) / Fv))
            lo = max(v for t, v in vals if t == 'lo'); hi = min(v for t, v in vals if t == 'up')
            if iv is None:
                if lo <= hi: bad += 1
            elif (lo, hi) != iv: bad += 1
    print(f"[opts] {n} poses x 2 lines: {bad} mismatches  {'ok' if bad == 0 else 'FAIL'}")
    return bad == 0


def t_poly(n=400, seed=2):
    rng = random.Random(seed); bad = 0; cert = 0
    for _ in range(n):
        deg = rng.randint(1, 7)
        p = [F(rng.randint(-20, 20), rng.randint(1, 5)) for _ in range(deg + 1)]
        u0 = F(rng.randint(0, 50), 100); u1 = u0 + F(rng.randint(1, 50), 100)
        if X.nonneg(p, u0, u1):
            cert += 1
            for k in range(201):
                if X.peval(p, u0 + (u1 - u0) * k / 200) < 0: bad += 1; break
    print(f"[poly] {n} random polys: {cert} certified >= 0, {bad} counterexamples  {'ok' if bad == 0 else 'FAIL'}")
    return bad == 0


def sample_poses(box, rng, k=60):
    cx0, cx1, cy0, cy1, u0, u1 = box
    out = []
    for i in range(k):
        # favour corners, faces and tiny u
        def pick(a, b):
            r = rng.random()
            if r < 0.25: return a
            if r < 0.5: return b
            return a + (b - a) * F(rng.randint(0, 10 ** 6), 10 ** 6)
        u = pick(u0, u1)
        if u == 0: u = u1 * F(1, rng.choice([10, 1000, 10 ** 6]))
        c, s = zm.trig(u); hw = (c + s) / 2
        cx = max(pick(cx0, cx1), hw); cy = max(pick(cy0, cy1), hw)
        if cx > cx1 or cy > cy1: continue
        out.append((cx, cy, u))
    return out


def t_exact(path, n=40, seed=3, region=None, sizes=(F(1, 40), F(1, 160)), umax=F(1, 50)):
    cov = ZM.Cover(MC.load(path))
    a, b = X.u_square(cov)
    ex = X.Exact(cov, a, b)
    rng = random.Random(seed)
    unsound = 0; certs = 0; comp_ok = 0; comp_tries = 0
    t0 = time.time()
    for it in range(n):
        h = rng.choice(sizes)
        if region:
            x0 = F(rng.randint(int(region[0] * 1000), int(region[1] * 1000)), 1000)
            y0 = F(rng.randint(int(region[2] * 1000), int(region[3] * 1000)), 1000)
        else:
            x0 = F(rng.randint(500, 3400), 1000); y0 = F(rng.randint(500, 3400), 1000)
        du = umax * F(rng.choice([1, 2, 4, 8, 32]), 32)
        u0 = rng.choice([F(0), F(0), du * rng.randint(1, 4)])
        box = (x0, x0 + h, y0, y0 + h, u0, u0 + du)
        P = sample_poses(box, rng)
        if not P: continue
        vals = [ZM.exact_mass(cov, *p) for p in P]
        mn = min(vals)
        # soundness: tau just above the sampled minimum must not certify
        tau = mn + F(1, 10 ** 9)
        ok, why = ex.certify(box, tau)
        if ok: unsound += 1; print("  UNSOUND at", [float(v) for v in box], float(mn), why)
        # completeness probe: tau a bit below
        tau2 = mn - F(1, 10 ** 4)
        ok2, why2 = ex.certify(box, tau2)
        comp_tries += 1; comp_ok += ok2
        if ok2: certs += 1
    print(f"[exact] {n} boxes: unsound certifications {unsound}; certified at tau = min - 1e-4: {comp_ok}/{comp_tries} "
          f"[{time.time()-t0:.0f}s] {'ok' if unsound == 0 else 'FAIL'}")
    return unsound == 0


def t_capleb(path, n=3000, seed=4):
    cov = ZM.Cover(MC.load(path)); a, b = X.u_square(cov)
    rng = random.Random(seed); bad = 0; nc = {'LEB': 0, 'CAP': 0}; checks = 0
    for it in range(n):
        h = F(1, rng.choice([10, 40, 160, 640]))
        x0 = F(rng.randint(1500, 3400), 1000); y0 = F(rng.randint(1500, 3400), 1000)
        du = F(1, rng.choice([2, 8, 32, 256])) / 2
        u0 = F(rng.randint(0, 3), 8) if rng.random() < 0.5 else F(0)
        box = (x0, x0 + h, y0, y0 + h, u0, min(u0 + du, F(1, 2)))
        B = zm.bin_data(box[4], box[5])
        kind = 'LEB' if X.cert_leb(box, B, a, b) else ('CAP' if X.cert_cap(box, B, cov, a, b) else None)
        if kind is None: continue
        nc[kind] += 1
        for p in sample_poses(box, rng, 12):
            checks += 1
            v = ZM.exact_mass(cov, *p)
            if v < 1: bad += 1; print("  VIOLATION", kind, [float(t) for t in box], [float(t) for t in p], float(v))
    print(f"[capleb] certified boxes {nc}, {checks} exact pose checks, {bad} below 1  {'ok' if bad == 0 else 'FAIL'}")
    return bad == 0


def float_mass(ex, cx, cy, u):
    """float mass of the lines + the Lebesgue square (independent of the exact path except the option formulas,
    which t_opts checks against a direct chord computation)."""
    tot = 0.0
    for L in ex.L:
        vals = [X.opt_float(o, cx, cy, u) for o in L['opts']]
        up = min(v for o, v in zip(L['opts'], vals) if o[0] == 'up')
        lo = max(v for o, v in zip(L['opts'], vals) if o[0] == 'lo')
        if up > lo: tot += L['prof'].Gf(up) - L['prof'].Gf(lo)
    # area of Q ∩ [a,b]^2 by clipping (float)
    th = 2 * math.atan(u); c, s = math.cos(th), math.sin(th)
    V = [(cx + (sx * c - sy * s) / 2, cy + (sx * s + sy * c) / 2) for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    a, b = float(ex.Ua), float(ex.Ub)
    for (nx, ny, cc) in ((-1, 0, -a), (0, -1, -a), (1, 0, b), (0, 1, b)):
        out = []
        for i in range(len(V)):
            P, Q = V[i], V[(i + 1) % len(V)]
            sp = nx * P[0] + ny * P[1] - cc; sq = nx * Q[0] + ny * Q[1] - cc
            if sp <= 0: out.append(P)
            if (sp < 0 < sq) or (sq < 0 < sp):
                t = sp / (sp - sq); out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
        V = out
        if not V: break
    ar = 0.0
    for i in range(len(V)):
        ar += V[i][0] * V[(i + 1) % len(V)][1] - V[(i + 1) % len(V)][0] * V[i][1]
    return tot + abs(ar) / 2


def random_cover(rng, m=4, nlines=10, a=None):
    """random axis-parallel piecewise-uniform lines + Lebesgue square [a, m-a]^2, D4-symmetrised, as a FORMAT dict"""
    D = 20; W = 1000
    if a is None: a = rng.choice([F(7, 5), F(6, 5), F(3, 2)])
    segs = []
    for _ in range(nlines):
        o = rng.choice('HV'); p = F(rng.randint(1, D * m - 1), D)
        x0 = F(rng.randint(0, D * m - 4), D); x1 = x0 + F(rng.randint(1, 20), D)
        x1 = min(x1, F(m))
        w = F(rng.randint(1, 400), W)
        for (X0, Y0, X1, Y1) in ((x0, p, x1, p), (p, x0, p, x1)) if o == 'H' else ((p, x0, p, x1),):
            for fx in (lambda x, y: (x, y), lambda x, y: (m - x, y), lambda x, y: (x, m - y), lambda x, y: (m - x, m - y),
                       lambda x, y: (y, x), lambda x, y: (m - y, x), lambda x, y: (y, m - x), lambda x, y: (m - y, m - x)):
                A = fx(X0, Y0); B = fx(X1, Y1)
                segs.append((int(A[0] * D), int(A[1] * D), int(B[0] * D), int(B[1] * D), int(w * W)))
    L = (m - 2 * a) ** 2
    cv = dict(s_num=m, s_den=1, D=D, W=W * 1, points=[], segments=segs,
              polygons=[(int(L * W), [(int(a * D), int(a * D)), (int((m - a) * D), int(a * D)),
                                       (int((m - a) * D), int((m - a) * D)), (int(a * D), int((m - a) * D))])])
    return cv


def t_adversarial(n=200, seed=11):
    """random covers, random boxes; tau = (dense float minimum over the box) + 1e-7: certification = soundness bug."""
    rng = random.Random(seed); unsound = 0; tried = 0; comp = 0; reasons = {}
    t0 = time.time()
    for it in range(n):
        mode = it % 3
        if mode == 0:
            cv = random_cover(rng, m=5, a=F(11, 5))          # U far from the boxes: regime 'none'
        else:
            cv = random_cover(rng, m=6, a=F(3, 2), nlines=14)   # boxes straddling y = a (1) or near the corner (2)
        cov = ZM.Cover(cv)
        a, b = X.u_square(cov)
        ex = X.Exact(cov, a, b)
        h = F(1, rng.choice([8, 20, 40]))
        if mode == 0: x0 = F(rng.randint(40, 130), 100); y0 = F(rng.randint(40, 130), 100)
        elif mode == 1: x0 = F(rng.randint(230, 330), 100); y0 = F(rng.randint(80, 190), 100)
        else: x0 = F(rng.randint(95, 160), 100); y0 = F(rng.randint(95, 160), 100)
        du = F(1, rng.choice([4, 10, 40, 200]))
        u0 = rng.choice([F(0), F(1, 20), F(1, 10)])
        box = (x0, x0 + h, y0, y0 + h, u0, min(u0 + du, F(39, 100)))
        # dense float min over admissible poses
        best = None
        for _ in range(3000):
            u = float(box[4]) + (float(box[5]) - float(box[4])) * rng.random() ** 2
            if u == 0: continue
            th = 2 * math.atan(u); hw = (math.cos(th) + math.sin(th)) / 2
            cx = float(box[0]) + float(h) * rng.random(); cy = float(box[2]) + float(h) * rng.random()
            if rng.random() < 0.3: cx = max(float(box[0]), hw)
            if rng.random() < 0.3: cy = max(float(box[2]), hw)
            if cx < hw or cy < hw or cx > float(box[1]) or cy > float(box[3]): continue
            if cx + hw > float(cov.m) or cy + hw > float(cov.m): continue
            v = float_mass(ex, cx, cy, u)
            if best is None or v < best[0]: best = (v, cx, cy, u)
        if best is None: continue
        # local polish
        v, cx, cy, u = best
        for r in range(400):
            sc = 0.3 ** (r // 100)
            c2 = cx + (rng.random() - .5) * float(h) * sc; d2 = cy + (rng.random() - .5) * float(h) * sc
            u2 = u + (rng.random() - .5) * float(box[5] - box[4]) * sc
            if not (float(box[4]) < u2 <= float(box[5])): continue
            th = 2 * math.atan(u2); hw = (math.cos(th) + math.sin(th)) / 2
            c2 = max(c2, hw); d2 = max(d2, hw)
            if not (float(box[0]) <= c2 <= float(box[1]) and float(box[2]) <= d2 <= float(box[3])): continue
            w2 = float_mass(ex, c2, d2, u2)
            if w2 < v: v, cx, cy, u = w2, c2, d2, u2
        tried += 1
        tau = F(v + 1e-7).limit_denominator(10 ** 12)
        ok, why = ex.certify(box, tau)
        if ok:
            unsound += 1
            print("  UNSOUND", [float(t) for t in box], 'float min', v, 'at', (cx, cy, u), why)
        ok2, why2 = ex.certify(box, F(v - 1e-5).limit_denominator(10 ** 12))
        comp += ok2
        reasons[why2 if not ok2 else 'ok'] = reasons.get(why2 if not ok2 else 'ok', 0) + 1
    print(f"[adversarial] {tried} random covers/boxes: unsound {unsound}; certified at min - 1e-5: {comp}/{tried} "
          f"[{time.time()-t0:.0f}s]  {'ok' if unsound == 0 else 'FAIL'}  reasons {reasons}")
    return unsound == 0


def t_gap(n=30, seed=21):
    """completeness: the largest tau EXACT certifies (bisection) against the float minimum over the box."""
    rng = random.Random(seed); gaps = []; unsound = 0
    for it in range(n):
        mode = it % 2
        cv = random_cover(rng, m=5, a=F(11, 5)) if mode == 0 else random_cover(rng, m=6, a=F(3, 2), nlines=14)
        cov = ZM.Cover(cv); a, b = X.u_square(cov); ex = X.Exact(cov, a, b)
        h = F(1, rng.choice([20, 40, 80]))
        if mode == 0: x0 = F(rng.randint(60, 130), 100); y0 = F(rng.randint(60, 130), 100)
        else: x0 = F(rng.randint(230, 330), 100); y0 = F(rng.randint(80, 190), 100)
        du = F(1, rng.choice([40, 200, 1000])); u0 = rng.choice([F(0), F(1, 20)])
        box = (x0, x0 + h, y0, y0 + h, u0, u0 + du)
        if ex.certify(box, F(0))[1] in ('rect', 'U-regime'): continue
        best = None
        for _ in range(4000):
            u = float(box[4]) + float(du) * rng.random()
            if u == 0: continue
            th = 2 * math.atan(u); hw = (math.cos(th) + math.sin(th)) / 2
            cx = float(box[0]) + float(h) * rng.choice([0, 1, rng.random()])
            cy = float(box[2]) + float(h) * rng.choice([0, 1, rng.random()])
            cx = max(cx, hw); cy = max(cy, hw)
            if cx > float(box[1]) or cy > float(box[3]): continue
            v = float_mass(ex, cx, cy, u)
            if best is None or v < best: best = v
        lo, hi = F(0), F(best).limit_denominator(10 ** 9) + F(1, 10 ** 7)
        if ex.certify(box, hi)[0]: unsound += 1; print("  UNSOUND", [float(t) for t in box]); continue
        if not ex.certify(box, lo)[0]: gaps.append(best); continue
        for _ in range(14):
            mid = (lo + hi) / 2
            if ex.certify(box, mid)[0]: lo = mid
            else: hi = mid
        gaps.append(best - float(lo))
    gaps.sort()
    print(f"[gap] {len(gaps)} boxes: float min - certified tau: median {gaps[len(gaps)//2]:.2e}, "
          f"80% {gaps[int(len(gaps)*0.8)]:.2e}, max {gaps[-1]:.2e}; unsound {unsound}")


def make_holes(path, outdir):
    """rejection-test covers built from the real one:
       R_scale : every segment mass x (1 - 1e-6)            (tight families fall to 1 - O(1e-6): theta = 0 and u -> 0)
       R_piece : one piece of the line y = 1 in the bottom band (x in [3, 3.2]) lighter by 1e-4 (and its D4 images)
       R_seam  : the seam piece x = 3, y in [0.8, 1.0] lighter by 1e-4 (and its D4 images)"""
    cv = MC.load(path); os.makedirs(outdir, exist_ok=True)
    D, W, m = cv['D'], cv['W'], cv['s_num'] // cv['s_den']
    K = 10 ** 6
    def rescale(cv2, f):
        c = dict(cv2); c['W'] = cv2['W'] * K
        c['segments'] = [(a, b, cc, d, f(s) ) for s in cv2['segments'] for (a, b, cc, d) in [s[:4]]]
        c['polygons'] = [(w * K, vs) for w, vs in cv2['polygons']]
        return c
    out = {}
    c = rescale(cv, lambda s: s[4] * (K - 1)); MC.write(os.path.join(outdir, 'R_scale.txt'), c, 'x (1 - 1e-6)'); out['R_scale'] = c
    def d4imgs(X0, Y0, X1, Y1):
        S = m * D
        fs = [lambda x, y: (x, y), lambda x, y: (S - x, y), lambda x, y: (x, S - y), lambda x, y: (S - x, S - y),
              lambda x, y: (y, x), lambda x, y: (S - y, x), lambda x, y: (y, S - x), lambda x, y: (S - y, S - x)]
        res = set()
        for f in fs:
            A = f(X0, Y0); B = f(X1, Y1); res.add(tuple(sorted([A, B])))
        return res
    for name, (X0, Y0, X1, Y1) in (('R_piece', (3 * D, 1 * D, 16 * D // 5, 1 * D)), ('R_seam', (3 * D, 4 * D // 5, 3 * D, 1 * D))):
        tgt = d4imgs(X0, Y0, X1, Y1)
        def f(s, tgt=tgt):
            key = tuple(sorted([(s[0], s[1]), (s[2], s[3])]))
            return s[4] * K - (s[4] * K // 10 ** 4 if key in tgt else 0)
        c = rescale(cv, f); hit = sum(1 for s in cv['segments'] if tuple(sorted([(s[0], s[1]), (s[2], s[3])])) in tgt)
        MC.write(os.path.join(outdir, name + '.txt'), c, name); out[name] = c
        print(name, 'pieces lightened:', hit)
    return out


def t_slope(path, lam_hi=F(1, 2), lam_lo=F(1, 20), n=40, seed=31, region=(2.5, 3.5, 0.5, 0.52), umax=F(1, 200)):
    """zero-margin test of Lemma E at u -> 0: on boxes touching u = 0 along a tight family, the claim
    mass >= 1 + lam u must be REFUSED when lam exceeds the float growth rate somewhere in the box (checked by float
    search) and should be certified for small lam."""
    cov = ZM.Cover(MC.load(path)); a, b = X.u_square(cov); ex = X.Exact(cov, a, b)
    rng = random.Random(seed); bad = 0; acc_lo = 0; tot = 0
    for it in range(n):
        h = F(1, rng.choice([80, 160, 320]))
        x0 = F(rng.randint(int(region[0] * 1000), int(region[1] * 1000)), 1000)
        y0 = F(region[2]).limit_denominator(1000) if rng.random() < 0.7 else F(rng.randint(int(region[2] * 1000), int(region[3] * 1000)), 1000)
        du = umax * F(1, rng.choice([1, 4, 16]))
        box = (x0, x0 + h, y0, y0 + h, F(0), du)
        # float minimum of (f - 1)/u over the box
        g = None
        for _ in range(3000):
            u = float(du) * 10 ** rng.uniform(-4, 0)
            th = 2 * math.atan(u); hw = (math.cos(th) + math.sin(th)) / 2
            cx = max(float(x0) + float(h) * rng.random(), hw)
            cy = max(float(y0) + float(h) * rng.choice([0, 0, rng.random()]), hw)
            if cx > float(box[1]) or cy > float(box[3]): continue
            v = (float_mass(ex, cx, cy, u) - 1) / u
            if g is None or v < g: g = v
        if g is None: continue
        tot += 1
        lam_hi = F(g + 0.02).limit_denominator(10 ** 6); lam_lo = F(max(0.0, g - 0.05)).limit_denominator(10 ** 6)
        okhi = ex.certify(box, F(1), lam_hi)[0]
        if okhi: bad += 1; print("  UNSOUND slope", [float(t) for t in box], 'float growth', g)
        oklo = ex.certify(box, F(1), lam_lo)[0]
        acc_lo += oklo
        if it < 8: print(f"   box {[round(float(t),5) for t in box]} float min growth {g:.4f}: lam=g+0.02 refused={not okhi}, lam=g-0.05 certified={oklo}")
    print(f"[slope] {tot} boxes at u -> 0: unsound {bad}; certified at lam = {float(lam_lo)}: {acc_lo}/{tot}  {'ok' if bad == 0 else 'FAIL'}")


def t_corner(n=3000, seed=81):
    """Lemma E'': at random poses of random boxes satisfying corner_ok, (i) the quadrilateral formula equals the exact
    area of Q ∩ {x < a, y < a} (polygon clipping), (ii) the McCormick/tangent minorant used by corner_term is <= it."""
    import qx2_exact as E2
    rng = random.Random(seed)
    cv = random_cover(rng, m=6, a=F(3, 2), nlines=4)
    cov = ZM.Cover(cv); a, b = X.u_square(cov); ex = X.Exact(cov, a, b)
    ok_boxes = 0; checks = 0; bad_formula = 0; bad_minor = 0; worst_loss = 0
    for it in range(n):
        h = F(1, rng.choice([10, 40, 160])); du = F(1, rng.choice([20, 100, 1000]))
        x0 = F(rng.randint(90, 170), 100); y0 = F(rng.randint(90, 170), 100)
        u0 = rng.choice([F(0), F(1, 50), F(1, 10)])
        box = (x0, x0 + h, y0, y0 + h, u0, min(u0 + du, F(39, 100)))
        if not ex.corner_ok(box): continue
        ok_boxes += 1
        for _ in range(4):
            cx = x0 + h * F(rng.randint(0, 1000), 1000); cy = y0 + h * F(rng.randint(0, 1000), 1000)
            u = box[4] + (box[5] - box[4]) * F(rng.randint(1, 1000), 1000)
            C, S = E2.trig(u); N = 1  # trig returns cos, sin
            Qv = E2.square_vertices(cx, cy, C, S)
            LL = E2.poly_area(E2.clip_hp(E2.clip_hp(Qv, 1, 0, a), 0, 1, a))      # x <= a, y <= a
            uu = F(u); Cp, Sp, Np = 1 - uu * uu, 2 * uu, 1 + uu * uu
            xBL = cx + (-C + S) / 2; yBL = cy + (-S - C) / 2
            al = (a - xBL) * Np; be = (a - yBL) * Np
            form = (2 * Cp * al * be + Sp * (be * be - al * al)) / (2 * Cp * Np * Np)
            # the minorant, with the box's alpha_lo, beta_lo, beta*
            anc = rng.choice(['lo', 'hi'])
            if anc == 'lo': al_lo = (a - box[1]) * Np + (Cp - Sp) / 2; be_lo = (a - box[3]) * Np + (Cp + Sp) / 2
            else: al_lo = (a - box[0]) * Np + (Cp - Sp) / 2; be_lo = (a - box[2]) * Np + (Cp + Sp) / 2
            bst = (a - (box[2] + box[3]) / 2) * Np + (Cp + Sp) / 2
            minor = (2 * Cp * (al_lo * be + be_lo * al - al_lo * be_lo) + Sp * (2 * bst * be - bst * bst) - Sp * al * al) / (2 * Cp * Np * Np)
            # and the RF path of corner_term at a degenerate 'candidate' (D = 1, X = cx, Y = cy)
            rf = ex.corner_term(box, [F(1)], [F(cx)], [F(cy)], anc)
            val = X.peval(rf.num, uu) / (Cp ** rf.e[0] * Sp ** rf.e[1] * Np ** rf.e[2])
            checks += 1
            if form != LL: bad_formula += 1
            if minor > LL or val != minor: bad_minor += 1
            worst_loss = max(worst_loss, float(LL - minor))
    print(f"[corner] {ok_boxes} boxes in the corner regime, {checks} exact poses: formula mismatches {bad_formula}, "
          f"minorant violations / RF mismatches {bad_minor}; max loss {worst_loss:.2e}  "
          f"{'ok' if bad_formula == 0 and bad_minor == 0 else 'FAIL'}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('what'); ap.add_argument('path', nargs='?')
    ap.add_argument('--n', type=int, default=40); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--region', type=str, default='')
    a = ap.parse_args()
    reg = [float(t) for t in a.region.split(',')] if (a.region and a.what != 'holes') else None
    if a.what == 'opts': t_opts()
    elif a.what == 'poly': t_poly()
    elif a.what == 'exact': t_exact(a.path, a.n, a.seed, reg)
    elif a.what == 'capleb': t_capleb(a.path, a.n, a.seed)
    elif a.what == 'adv': t_adversarial(a.n, a.seed)
    elif a.what == 'gap': t_gap(a.n, a.seed)
    elif a.what == 'holes': make_holes(a.path, a.region or 'holes')
    elif a.what == 'slope': t_slope(a.path, n=a.n, seed=a.seed, region=reg or (2.5, 3.5, 0.5, 0.52))
    elif a.what == 'corner': t_corner(a.n, a.seed)
