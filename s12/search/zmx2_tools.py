#!/usr/bin/env python3
"""zmx2_tools.py -- test harness for the independent mixed-cover checker `zmx2` (search/ZMX2.md).

Written independently of search/zm_mixed.py (never read).  Own parser, own exact evaluator.

  load / write                 mixed v1 files (tasks/line-cover/FORMAT.md), exact integers
  mu_exact(cv, x, y, u)        exact mu(Q) at a rational pose (Fractions; any segment direction)
  admissible(cv, x, y, u)      exact admissibility of the pose (Q inside [0,s]^2)

Subcommands:
  sound FILE [--n N] [--poses P] [--seed S] [--refl] [--near X,Y,U] [--zflags f1,f2] [--area] [--small-u]
        soundness harness: random small pose boxes (biased to germs, walls, grid lines), the
        bound printed by `zmx2 boxes`, and exact mu at P random admissible rational poses of
        each box (plus corners); FAIL if any mu < bound.
  sound0 FILE [--n N] [--poses P] [--seed S]   theta = 0 harness for `zmx2 boxes0` (ZMX2_AREA.md)
  mu FILE x y u                exact mu at a rational pose (x, y, u as p/q)
  toy OUT --s S --kind K [...] toy covers for the certify / reject tests (see ZMX2.md sec 8)
  fmin FILE [--pitch P] [--ubins N]   float minimum over a grid of admissible poses (all angles)
  perturb IN OUT --op ...      rejection-test covers (drop mass in a region; break symmetry)
  uncert LOG FILE [--max N]    exact mu at the corners of the UNCERT boxes of a cert log (min, pose)
"""
import math
import os
import random
import subprocess
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ZMX2 = os.environ.get('ZMX2') or os.path.join(HERE, '..', 'verify2', 'target', 'release', 'zmx2')


def load(path):
    toks = []
    with open(path) as f:
        for line in f:
            line = line.split('#', 1)[0]
            toks.extend(line.split())
    i = 0
    mixed = toks[0] == 'mixed'
    if mixed:
        assert toks[1] == '1'
        i = 2

    def nx():
        nonlocal i
        v = int(toks[i])
        i += 1
        return v
    s_num, s_den, D, W = nx(), nx(), nx(), nx()
    pts = [(nx(), nx(), nx()) for _ in range(nx())]
    segs = []
    polys = []
    if mixed:
        segs = [(nx(), nx(), nx(), nx(), nx()) for _ in range(nx())]
        for _ in range(nx()):   # polygons: k w X1 Y1 ... Xk Yk (convex, ccw; mass w spread by area)
            k = nx()
            w = nx()
            polys.append((w, [(nx(), nx()) for _ in range(k)]))
    assert i == len(toks)
    return dict(s=Fr(s_num, s_den), s_num=s_num, s_den=s_den, D=D, W=W, points=pts, segments=segs, polygons=polys)


def write(path, cv, comment=None):
    with open(path, 'w') as f:
        f.write('mixed 1\n')
        if comment:
            for l in comment.splitlines():
                f.write('# %s\n' % l)
        f.write('%d %d\n%d\n%d\n' % (cv['s_num'], cv['s_den'], cv['D'], cv['W']))
        f.write('%d\n' % len(cv['points']))
        for p in cv['points']:
            f.write('%d %d %d\n' % p)
        f.write('%d\n' % len(cv['segments']))
        for sg in cv['segments']:
            f.write('%d %d %d %d %d\n' % sg)
        pg = cv.get('polygons', [])
        f.write('%d\n' % len(pg))
        for (w, vs) in pg:
            f.write('%d %d %s\n' % (len(vs), w, ' '.join('%d %d' % v for v in vs)))


def total(cv):
    return Fr(sum(p[2] for p in cv['points']) + sum(s[4] for s in cv['segments'])
              + sum(p[0] for p in cv.get('polygons', [])), cv['W'])


def poly_area(vs):
    """exact signed area (shoelace) of a polygon given as a list of Fraction/int pairs"""
    n = len(vs)
    return sum(vs[k][0] * vs[(k + 1) % n][1] - vs[(k + 1) % n][0] * vs[k][1] for k in range(n)) / Fr(2)


def clip_halfplane(vs, a, b, c):
    """Sutherland-Hodgman: the part of the convex polygon vs with a x + b y <= c (exact)"""
    out = []
    n = len(vs)
    for k in range(n):
        p, q = vs[k], vs[(k + 1) % n]
        fp = a * p[0] + b * p[1] - c
        fq = a * q[0] + b * q[1] - c
        if fp <= 0:
            out.append(p)
        if (fp < 0 < fq) or (fq < 0 < fp):
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def square_corners(x, y, u):
    """the closed unit square Q(x, y, u) as a ccw list of exact corners"""
    c, s = trig(u)
    h = Fr(1, 2)
    return [(x + c * a - s * b, y + s * a + c * b) for (a, b) in ((-h, -h), (h, -h), (h, h), (-h, h))]


def area_in_poly(x, y, u, P):
    """exact area of Q(x, y, u) cap P, P a convex ccw polygon of Fractions"""
    vs = square_corners(x, y, u)
    n = len(P)
    for k in range(n):
        (px, py), (qx, qy) = P[k], P[(k + 1) % n]
        # inside of a ccw edge p->q: cross(q - p, z - p) >= 0  <=>  a zx + b zy <= c
        a, b = (qy - py), -(qx - px)
        c = a * px + b * py
        vs = clip_halfplane(vs, a, b, c)
        if len(vs) < 3:
            return Fr(0)
    return poly_area(vs)


def trig(u):
    n = 1 + u * u
    return (1 - u * u) / n, 2 * u / n


def admissible(cv, x, y, u):
    c, s = trig(u)
    w = c + s  # u in [0,1]: theta in [0,90]
    S = cv['s']
    return w / 2 <= x <= S - w / 2 and w / 2 <= y <= S - w / 2


def _tint(a, b, lo, hi, cur):
    """restrict the t-interval cur=(t0,t1) to {t : lo <= a + b t <= hi}"""
    t0, t1 = cur
    if t0 > t1:
        return cur
    if b == 0:
        return cur if lo <= a <= hi else (Fr(1), Fr(0))
    ta, tb = (lo - a) / b, (hi - a) / b
    if b < 0:
        ta, tb = tb, ta
    return (max(t0, ta), min(t1, tb))


def mu_exact(cv, x, y, u, prefilter=True):
    """exact mu(Q(x,y,u)); x, y, u Fractions, u in [0,1)."""
    c, s = trig(u)
    D, W = cv['D'], cv['W']
    h = Fr(1, 2)
    xf, yf = float(x), float(y)
    tot = Fr(0)
    for (X, Y, w) in cv['points']:
        if w == 0:
            continue
        if prefilter and (abs(X / D - xf) > 0.72 or abs(Y / D - yf) > 0.72):
            continue
        a, b = Fr(X, D) - x, Fr(Y, D) - y
        if abs(a * c + b * s) <= h and abs(-a * s + b * c) <= h:
            tot += w
    for (X0, Y0, X1, Y1, w) in cv['segments']:
        if w == 0:
            continue
        if prefilter:
            if min(X0, X1) / D > xf + 0.72 or max(X0, X1) / D < xf - 0.72:
                continue
            if min(Y0, Y1) / D > yf + 0.72 or max(Y0, Y1) / D < yf - 0.72:
                continue
        a0, b0 = Fr(X0, D) - x, Fr(Y0, D) - y
        da, db = Fr(X1 - X0, D), Fr(Y1 - Y0, D)
        # X(t) = (a0 + t da) c + (b0 + t db) s ; Y(t) = -(a0 + t da) s + (b0 + t db) c
        cur = (Fr(0), Fr(1))
        cur = _tint(a0 * c + b0 * s, da * c + db * s, -h, h, cur)
        cur = _tint(-a0 * s + b0 * c, -da * s + db * c, -h, h, cur)
        if cur[1] > cur[0]:
            tot += w * (cur[1] - cur[0])
    for (w, vs) in cv.get('polygons', []):
        if w == 0:
            continue
        P = [(Fr(X, D), Fr(Y, D)) for (X, Y) in vs]
        if prefilter and (min(p[0] for p in P) > xf + 0.72 or max(p[0] for p in P) < xf - 0.72
                          or min(p[1] for p in P) > yf + 0.72 or max(p[1] for p in P) < yf - 0.72):
            continue
        tot += w * area_in_poly(x, y, u, P) / poly_area(P)
    return tot / W


def reflect_y(cv):
    S = cv['s'] * cv['D']
    assert S.denominator == 1
    S = S.numerator
    out = dict(cv)
    out['points'] = [(X, S - Y, w) for (X, Y, w) in cv['points']]
    out['segments'] = [(a, S - b, c, S - d, w) for (a, b, c, d, w) in cv['segments']]
    # reflection reverses orientation: reverse the vertex order to stay ccw
    out['polygons'] = [(w, [(X, S - Y) for (X, Y) in vs][::-1]) for (w, vs) in cv.get('polygons', [])]
    return out


def parse_fr(t):
    return Fr(t)


# ----------------------------------------------------------------------------- soundness harness

def rand_box(rng, S):
    """random pose box with centre denominators 10*2^a and u denominators 8*2^b."""
    a = rng.randint(2, 11)
    b = rng.randint(1, 11)
    cd, ud = 10 << a, 8 << b
    kind = rng.random()
    Sf = float(S)
    if kind < 0.3:      # near a tile centre (germ), small angle
        cx = rng.randint(0, int(Sf) - 1) + 0.5 + rng.uniform(-0.02, 0.02)
        cy = rng.randint(0, int(Sf) - 1) + 0.5 + rng.uniform(-0.02, 0.02)
        uc = rng.choice([0.0, 0.0, rng.uniform(0, 0.02)])
    elif kind < 0.5:    # near a wall
        cx = rng.choice([0.5, Sf - 0.5]) + rng.uniform(-0.02, 0.05) * (1 if rng.random() < .5 else -1)
        cy = rng.uniform(0.4, Sf - 0.4)
        uc = rng.choice([0.0, rng.uniform(0, 0.5)])
        if rng.random() < 0.5:
            cx, cy = cy, cx
    else:               # anywhere
        cx, cy = rng.uniform(0.4, Sf - 0.4), rng.uniform(0.4, Sf - 0.4)
        uc = rng.uniform(0, 0.5)
    i = int(cx * cd)
    j = int(cy * cd)
    k = min(int(uc * ud), ud // 2 - 1)
    wx = rng.choice([1, 1, 2, 4])
    wy = rng.choice([1, 1, 2, 4])
    wu = rng.choice([1, 1, 2])
    return (Fr(i, cd), Fr(i + wx, cd), Fr(j, cd), Fr(j + wy, cd), Fr(k, ud), Fr(min(k + wu, ud // 2), ud))


def near_box(rng, x, y, u, r=0.002, ru=0.001):
    """random small pose box near the pose (x, y, u): centre within r, u within ru (u >= 0), sides
    1/(10*2^a) and 1/(8*2^b) for a in 7..16, b in 6..16 (the scale of a germ's uncertified leaves)."""
    a = rng.randint(7, 16)
    b = rng.randint(6, 16)
    cd, ud = 10 << a, 8 << b
    cx = x + rng.uniform(-r, r)
    cy = y + rng.uniform(-r, r)
    uc = max(0.0, u + rng.uniform(-ru, ru)) if rng.random() < 0.85 else 0.0
    i, j, k = int(cx * cd), int(cy * cd), int(uc * ud)
    wx, wy, wu = rng.choice([1, 1, 2]), rng.choice([1, 1, 2]), rng.choice([1, 1, 2])
    return (Fr(i, cd), Fr(i + wx, cd), Fr(j, cd), Fr(j + wy, cd), Fr(k, ud), Fr(k + wu, ud))


def cmd_sound(argv):
    path = argv[0]
    n = int(opt(argv, '--n', 200))
    npose = int(opt(argv, '--poses', 25))
    seed = int(opt(argv, '--seed', 1))
    refl = '--refl' in argv
    rng = random.Random(seed)
    small_u = '--small-u' in argv
    cv = load(path)
    cve = reflect_y(cv) if refl else cv
    near = opt(argv, '--near', None)
    if near is None and '--area' in argv:
        boxes = [rand_box_area(rng, cve) for _ in range(n)]
    elif near is None:
        boxes = [rand_box(rng, cv['s']) for _ in range(n)]
    else:
        boxes = [near_box(rng, *(float(t) for t in near.split(','))) for _ in range(n)]
    inp = '\n'.join(','.join('%d/%d' % (f.numerator, f.denominator) for f in bx) for bx in boxes) + '\n'
    zflags = ['--' + t for t in opt(argv, '--zflags', '').split(',') if t]  # e.g. pair-points,sym-atoms
    cmd = [ZMX2, 'boxes', path] + (['--refl'] if refl else []) + zflags
    out = subprocess.run(cmd, input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
    worst = None
    nfail = 0
    nposes = 0
    tight = []
    for bx, line in zip(boxes, out):
        bound, unit, empty, _ = (int(t) for t in line.split())
        B = Fr(bound, unit)
        x0, x1, y0, y1, u0, u1 = bx
        poses = []
        for cx in (x0, x1):
            for cy in (y0, y1):
                for cu in (u0, u1):
                    poses.append((cx, cy, cu))
        for _ in range(npose):
            m = 1009
            fx, fy, fu = (Fr(rng.randint(0, m), m) for _ in range(3))
            if rng.random() < 0.25:
                fu = Fr(0)
            elif small_u and rng.random() < 0.6:
                # tiny angles relative to the box: explores the germ variables (offset / theta) of u0 = 0 boxes
                fu = Fr(1, 2 ** rng.randint(1, 40))
            poses.append((x0 + (x1 - x0) * fx, y0 + (y1 - y0) * fy, u0 + (u1 - u0) * fu))
        mmin = None
        for (px, py, pu) in poses:
            if not admissible(cve, px, py, pu):
                continue
            nposes += 1
            if empty:
                print('FAIL: box reported EMPTY has admissible pose', bx, (px, py, pu))
                nfail += 1
                continue
            mu = mu_exact(cve, px, py, pu)
            if mmin is None or mu < mmin:
                mmin = mu
            if mu < B:
                nfail += 1
                print('FAIL: mu %s < bound %s at pose %s in box %s' % (float(mu), float(B), (px, py, pu), bx))
        if mmin is not None and not empty:
            gap = float(mmin - B)
            tight.append(gap)
            if worst is None or gap < worst[0]:
                worst = (gap, bx)
    tight.sort()
    print('sound: %d boxes, %d admissible poses checked, %d FAIL; min(mu)-bound: min %.3g, median %.3g%s' % (
        n, nposes, nfail, tight[0] if tight else float('nan'), tight[len(tight) // 2] if tight else float('nan'),
        '' if not refl else ' (reflected cover)'))
    return 1 if nfail else 0


def rand_box_area(rng, cv):
    """random small pose box near the boundary of a rectangle of the cover (area densities), or inside it:
    centre within 0.75 of a side (both sides of it), any angle, small angles favoured."""
    a = rng.randint(3, 12)
    b = rng.randint(2, 12)
    cd, ud = 10 << a, 8 << b
    D = cv['D']
    w, vs = rng.choice(cv['polygons'])
    xs = [X / D for (X, Y) in vs]
    ys = [Y / D for (X, Y) in vs]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    kind = rng.random()
    if kind < 0.6:      # near a side
        side = rng.randint(0, 3)
        t = rng.uniform(-0.75, 0.75)
        if side < 2:
            cx = (x0 if side == 0 else x1) + t
            cy = rng.uniform(y0 - 0.5, y1 + 0.5)
        else:
            cy = (y0 if side == 2 else y1) + t
            cx = rng.uniform(x0 - 0.5, x1 + 0.5)
    elif kind < 0.8:    # near a corner
        cx = rng.choice([x0, x1]) + rng.uniform(-0.75, 0.75)
        cy = rng.choice([y0, y1]) + rng.uniform(-0.75, 0.75)
    else:               # inside, touching from inside: centre at w/2 from a side
        uc = rng.uniform(0, 0.5)
        c, s = trig(Fr(uc).limit_denominator(1000))
        hw = float(c + s) / 2
        cx = x0 + hw + rng.uniform(-0.01, 0.01)
        cy = rng.uniform(y0 + 0.8, y1 - 0.8)
        if rng.random() < 0.5:
            cx, cy = cy, cx
    uc = rng.choice([0.0, rng.uniform(0, 0.03), rng.uniform(0, 0.5)])
    S = float(cv['s'])
    cx = min(max(cx, 0.45), S - 0.45)
    cy = min(max(cy, 0.45), S - 0.45)
    i, j = int(cx * cd), int(cy * cd)
    k = min(int(uc * ud), ud // 2 - 1)
    wx, wy, wu = rng.choice([1, 1, 2, 4]), rng.choice([1, 1, 2, 4]), rng.choice([1, 1, 2])
    return (Fr(i, cd), Fr(i + wx, cd), Fr(j, cd), Fr(j + wy, cd), Fr(k, ud), Fr(min(k + wu, ud // 2), ud))


def cmd_sound0(argv):
    """sound0 FILE [--n N] [--poses P] [--seed S]: theta = 0 harness for `zmx2 boxes0` (Lemma Z0): random
    centre boxes (biased to walls, rectangle sides, and the half-integer offsets of lines), exact mu of the
    axis-parallel square at random rational centres of each (and its corners); FAIL if any mu < bound."""
    path = argv[0]
    n = int(opt(argv, '--n', 300))
    npose = int(opt(argv, '--poses', 20))
    rng = random.Random(int(opt(argv, '--seed', 1)))
    cv = load(path)
    S = cv['s']
    Sf = float(S)
    D = cv['D']
    special = [Fr(1, 2), S - Fr(1, 2)]
    for (X0, Y0, X1, Y1, w) in cv['segments']:
        for v in (X0, Y0):
            special += [Fr(v, D) - Fr(1, 2), Fr(v, D) + Fr(1, 2)]
    for (w, vs) in cv.get('polygons', []):
        for (X, Y) in vs:
            for v in (X, Y):
                special += [Fr(v, D) - Fr(1, 2), Fr(v, D) + Fr(1, 2), Fr(v, D)]
    special = sorted(set(v for v in special if Fr(1, 2) <= v <= S - Fr(1, 2)))
    boxes = []
    for _ in range(n):
        a = rng.randint(2, 14)
        cd = 10 << a
        cs = []
        for _ in range(2):
            if rng.random() < 0.7:
                v = float(rng.choice(special)) + rng.choice([0, 0, -1, 1]) * rng.random() * 2.0 / cd
            else:
                v = rng.uniform(0.5, Sf - 0.5)
            i = int(v * cd)
            wdt = rng.choice([1, 1, 2, 3])
            lo = max(Fr(i, cd), Fr(1, 2))
            hi = min(Fr(i + wdt, cd), S - Fr(1, 2))
            if lo >= hi:
                lo, hi = Fr(cd // 2, cd), Fr(cd // 2 + 1, cd)
            cs.append((lo, hi))
        boxes.append((cs[0][0], cs[0][1], cs[1][0], cs[1][1]))
    inp = '\n'.join(','.join('%d/%d' % (f.numerator, f.denominator) for f in bx) for bx in boxes) + '\n'
    out = subprocess.run([ZMX2, 'boxes0', path], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
    nfail = nposes = 0
    gaps = []
    for bx, line in zip(boxes, out):
        bound, unit = (int(t) for t in line.split())
        B = Fr(bound, unit)
        x0, x1, y0, y1 = bx
        poses = [(x, y) for x in (x0, x1) for y in (y0, y1)]
        for _ in range(npose):
            m = 1009
            poses.append((x0 + (x1 - x0) * Fr(rng.randint(0, m), m), y0 + (y1 - y0) * Fr(rng.randint(0, m), m)))
        mn = None
        for (px, py) in poses:
            if not admissible(cv, px, py, Fr(0)):
                continue
            nposes += 1
            mu = mu_exact(cv, px, py, Fr(0))
            mn = mu if mn is None or mu < mn else mn
            if mu < B:
                nfail += 1
                print('FAIL: mu %s < bound %s at (%s, %s) in box %s' % (float(mu), float(B), px, py, bx))
        if mn is not None:
            gaps.append(float(mn - B))
    gaps.sort()
    print('sound0: %d boxes, %d admissible poses checked, %d FAIL; min(mu)-bound: min %.3g, median %.3g' % (
        n, nposes, nfail, gaps[0] if gaps else float('nan'), gaps[len(gaps) // 2] if gaps else float('nan')))
    return 1 if nfail else 0


def cmd_tight(argv):
    """tight DUMP FILE [--poses P] [--max N] [--seed S]: the certified leaves with the smallest
    margin (zmx2 cert --tight T --dump-tight DUMP); exact mu at random admissible rational poses of
    each must be >= its bound (and hence >= 1)."""
    dump, path = argv[0], argv[1]
    npose = int(opt(argv, '--poses', 20))
    mx = int(opt(argv, '--max', 400))
    rng = random.Random(int(opt(argv, '--seed', 7)))
    cv = load(path)
    cvs = [cv, reflect_y(cv)]
    rows = []
    with open(dump) as f:
        for line in f:
            pas, bx, bd = line.split()
            rows.append((int(pas), [Fr(t) for t in bx.split(',')], int(bd)))
    rng.shuffle(rows)
    rows = rows[:mx]
    # unit = W * Lc * 2^30 ; recompute Lc
    Lc = 1
    for (a, b, c, d, w) in cv['segments']:
        if w:
            Lc = Lc * (abs(c - a) + abs(d - b)) // math.gcd(Lc, abs(c - a) + abs(d - b))
    unit = cv['W'] * Lc * 2 ** 30
    nfail = nch = 0
    gaps = []
    for pas, bx, bd in rows:
        B = Fr(bd, unit)
        x0, x1, y0, y1, u0, u1 = bx
        cve = cvs[pas]
        mn = None
        for k in range(npose + 8):
            if k < 8:
                fx, fy, fu = (k >> 2) & 1, (k >> 1) & 1, k & 1
                fx, fy, fu = Fr(fx), Fr(fy), Fr(fu)
            else:
                m = 1013
                fx, fy, fu = (Fr(rng.randint(0, m), m) for _ in range(3))
            px, py, pu = x0 + (x1 - x0) * fx, y0 + (y1 - y0) * fy, u0 + (u1 - u0) * fu
            if not admissible(cve, px, py, pu):
                continue
            nch += 1
            mu = mu_exact(cve, px, py, pu)
            mn = mu if mn is None or mu < mn else mn
            if mu < B:
                nfail += 1
                print('FAIL: mu %.12f < bound %.12f at %s in %s' % (float(mu), float(B), (px, py, pu), bx))
        if mn is not None:
            gaps.append(float(mn - B))
    gaps.sort()
    print('tight: %d leaves, %d admissible poses, %d FAIL; min(mu)-bound: min %.3g median %.3g; bounds in [%.6f, %.6f]' % (
        len(rows), nch, nfail, gaps[0] if gaps else float('nan'), gaps[len(gaps) // 2] if gaps else float('nan'),
        min(float(Fr(r[2], unit)) for r in rows), max(float(Fr(r[2], unit)) for r in rows)))
    return 1 if nfail else 0


def cmd_uncert(argv):
    """uncert LOG FILE [--max N]: exact mu at the 8 corners and the centre of the first N (default 20) UNCERT
    boxes of a `zmx2 cert` log/output (pass 1 boxes on the reflected cover); prints the smallest admissible value
    and its pose.  For rejection tests: a refusal is a true one if some uncertified box contains a pose with mu < 1."""
    log, path = argv[0], argv[1]
    mx = int(opt(argv, '--max', 20))
    cv = load(path)
    cvs = [cv, reflect_y(cv)]
    best = None
    n = 0
    for line in open(log):
        if not line.startswith('UNCERT ') or n >= mx:
            continue
        f = line.split()
        pas = int(f[3])
        x0, x1, y0, y1, u0, u1 = (Fr(t) for t in f[5].split(','))
        n += 1
        poses = [(x, y, u) for x in (x0, x1) for y in (y0, y1) for u in (u0, u1)]
        poses.append(((x0 + x1) / 2, (y0 + y1) / 2, (u0 + u1) / 2))
        for (x, y, u) in poses:
            if admissible(cvs[pas], x, y, u):
                m = mu_exact(cvs[pas], x, y, u)
                if best is None or m < best[0]:
                    best = (m, pas, x, y, u)
    if best is None:
        print('uncert: %d boxes, no admissible corner' % n)
        return 1
    m, pas, x, y, u = best
    print('uncert: %d boxes; min exact mu %.9f at pass %d (%s, %s, u %s) = (%.6f, %.6f, %.4f deg)' % (
        n, float(m), pas, x, y, u, float(x), float(y), 2 * math.degrees(math.atan(float(u)))))
    return 0


def opt(argv, key, dflt):
    if key in argv:
        return argv[argv.index(key) + 1]
    return dflt


# ----------------------------------------------------------------------------- float scan

def fmass(cv, x, y, th):
    c, s = math.cos(th), math.sin(th)
    D, W = cv['D'], cv['W']
    tot = 0.0
    h = 0.5 + 1e-12
    for (X, Y, w) in cv['points']:
        a, b = X / D - x, Y / D - y
        if abs(a * c + b * s) <= h and abs(-a * s + b * c) <= h:
            tot += w
    for (X0, Y0, X1, Y1, w) in cv['segments']:
        a0, b0 = X0 / D - x, Y0 / D - y
        da, db = (X1 - X0) / D, (Y1 - Y0) / D
        t0, t1 = 0.0, 1.0
        for (a, b) in ((a0 * c + b0 * s, da * c + db * s), (-a0 * s + b0 * c, -da * s + db * c)):
            if abs(b) < 1e-15:
                if abs(a) > h:
                    t0, t1 = 1, 0
                continue
            ta, tb = (-h - a) / b, (h - a) / b
            if b < 0:
                ta, tb = tb, ta
            t0, t1 = max(t0, ta), min(t1, tb)
        if t1 > t0:
            tot += w * (t1 - t0)
    if cv.get('polygons'):
        h = 0.5
        corners = [(x + c * a - s * b, y + s * a + c * b) for (a, b) in ((-h, -h), (h, -h), (h, h), (-h, h))]
        for (w, vs) in cv['polygons']:
            P = [(X / D, Y / D) for (X, Y) in vs]
            q = corners
            for k in range(len(P)):
                (px, py), (qx, qy) = P[k], P[(k + 1) % len(P)]
                a, b = (qy - py), -(qx - px)
                q = clip_halfplane(q, a, b, a * px + b * py)
                if len(q) < 3:
                    break
            if len(q) >= 3:
                tot += w * float(poly_area(q)) / float(poly_area(P))
    return tot / W


def fmin(cv, pitch=0.05, nth=46, thmax=90.0):
    S = float(cv['s'])
    best = (9e9, None)
    n = int(round(S / pitch))
    for k in range(nth + 1):
        th = math.radians(thmax * k / nth)
        w = math.cos(th) + math.sin(th)
        for i in range(n + 1):
            x = i * pitch
            if x < w / 2 - 1e-12 or x > S - w / 2 + 1e-12:
                continue
            for j in range(n + 1):
                y = j * pitch
                if y < w / 2 - 1e-12 or y > S - w / 2 + 1e-12:
                    continue
                m = fmass(cv, x, y, th)
                if m < best[0]:
                    best = (m, (x, y, math.degrees(th)))
    return best


# ----------------------------------------------------------------------------- toy covers

def toy(s, kind, D=1000, W=10 ** 9, rho=None, q=10, extra_pts=()):
    """Toy mixed covers of [0,s]^2.

    kind 'lines':  uniform density rho (mass per unit length, rational) on every interior grid
                   line x = k, y = k (k = 1..s-1), as pieces of length 1/q.
    kind 'wave':   same lines, piece i carries rho * (1 + (i mod 3)/4) (non-uniform, D4-symmetric
                   pieces), exercising the pair lemma with unequal densities.
    extra_pts:     (x, y, w) floats -> points with weight w (converted to integers)
    """
    segs = []
    L = D // q
    for k in range(1, s):
        for i in range(s * q):
            if kind == 'lines':
                f = Fr(1)
            else:
                ii = min(i, s * q - 1 - i)  # symmetric about the line's midpoint
                f = 1 + Fr(ii % 3, 4)
            w = rho * f * Fr(1, q) * W
            wn = -(-w.numerator // w.denominator)  # round up
            segs.append((k * D, i * L, k * D, (i + 1) * L, wn))
            segs.append((i * L, k * D, (i + 1) * L, k * D, wn))
    pts = [(int(round(x * D)), int(round(y * D)), int(math.ceil(w * W))) for (x, y, w) in extra_pts]
    return dict(s=Fr(s), s_num=s, s_den=1, D=D, W=W, points=pts, segments=segs)


def cmd_toy(argv):
    out = argv[0]
    s = int(opt(argv, '--s', 3))
    kind = opt(argv, '--kind', 'lines')
    if kind == 'lebesgue':
        # area density rho on the whole container [0,s]^2 (ZMX2_AREA.md sec 10): mu(Q) = rho for every
        # admissible square; rho = 1 is valid with zero margin everywhere
        rho = Fr(opt(argv, '--rho', '1'))
        D, W = 10, 10 ** 9
        w = rho * s * s * W
        assert w.denominator == 1
        cv = dict(s=Fr(s), s_num=s, s_den=1, D=D, W=W, points=[], segments=[],
                  polygons=[(int(w), [(0, 0), (s * D, 0), (s * D, s * D), (0, s * D)])])
        write(out, cv, comment='toy lebesgue s=%d rho=%s' % (s, rho))
        print('wrote %s: total %.9f' % (out, float(total(cv))))
        return 0
    rho = Fr(opt(argv, '--rho', '1'))
    q = int(opt(argv, '--q', 10))
    pts = []
    if '--centre-pts' in argv:  # a point at every tile centre
        wc = float(opt(argv, '--centre-w', 0.1))
        pts = [(i + 0.5, j + 0.5, wc) for i in range(s) for j in range(s)]
    cv = toy(s, kind, rho=rho, q=q, extra_pts=pts)
    write(out, cv, comment='toy %s s=%d rho=%s q=%d centre-pts=%s' % (kind, s, rho, q, bool(pts)))
    print('wrote %s: total %.6f' % (out, float(total(cv))))


def cmd_perturb(argv):
    """perturb IN OUT --op zero-seg --line x|y --at K --lo A --hi B   (zero the pieces of line
    x=K (or y=K) inside [A,B]);  --op scale --f F (all masses x F, rounded down);
    --op pt-weight --at X,Y --w W (set the weight of the point at X,Y, file units);
    --op bump (add 1 to the first point's weight: breaks D4);
    --op rect-scale --f F;  --op rect-delta --dw DW;  --op seg-delta --at X0,Y0,X1,Y1 --dw DW."""
    cv = load(argv[0])
    op = opt(argv, '--op', '')
    D = cv['D']
    if op == 'zero-seg':
        ori = opt(argv, '--line', 'x')
        K = Fr(opt(argv, '--at', '1')) * D
        lo, hi = Fr(opt(argv, '--lo', '0')) * D, Fr(opt(argv, '--hi', '0')) * D
        segs = []
        for (a, b, c, d, w) in cv['segments']:
            if ori == 'x' and a == c == K and lo <= min(b, d) and max(b, d) <= hi:
                w = 0
            if ori == 'y' and b == d == K and lo <= min(a, c) and max(a, c) <= hi:
                w = 0
            segs.append((a, b, c, d, w))
        cv['segments'] = segs
    elif op == 'scale':
        f = Fr(opt(argv, '--f', '1'))
        cv['points'] = [(x, y, int(w * f)) for (x, y, w) in cv['points']]
        cv['segments'] = [(a, b, c, d, int(w * f)) for (a, b, c, d, w) in cv['segments']]
    elif op == 'pt-weight':
        # set the weight of the point(s) at X,Y (file units) to W (file units)
        X, Y = (int(t) for t in opt(argv, '--at', '').split(','))
        W = int(opt(argv, '--w', '0'))
        hit = [k for k, p in enumerate(cv['points']) if (p[0], p[1]) == (X, Y)]
        if not hit:
            raise SystemExit('no point at %d,%d' % (X, Y))
        for n, k in enumerate(hit):
            cv['points'][k] = (X, Y, W if n == 0 else 0)
    elif op == 'rect-scale':
        # multiply every polygon weight by F (rounded down)
        f = Fr(opt(argv, '--f', '1'))
        cv['polygons'] = [(int(w * f), vs) for (w, vs) in cv['polygons']]
    elif op == 'rect-delta':
        # add DW (file units, may be negative) to the weight of the first polygon
        dw = int(opt(argv, '--dw', '0'))
        w, vs = cv['polygons'][0]
        cv['polygons'][0] = (w + dw, vs)
    elif op == 'seg-delta':
        # add DW (file units) to the weight of the segment X0,Y0,X1,Y1 (either orientation)
        a = tuple(int(t) for t in opt(argv, '--at', '').split(','))
        dw = int(opt(argv, '--dw', '0'))
        hit = 0
        segs = []
        for sg in cv['segments']:
            if sg[:4] == a or (sg[2], sg[3], sg[0], sg[1]) == a:
                sg = sg[:4] + (sg[4] + dw,)
                hit += 1
            segs.append(sg)
        if hit != 1:
            raise SystemExit('segment %s found %d times' % (a, hit))
        cv['segments'] = segs
    elif op == 'bump':
        p = cv['points'][0]
        cv['points'][0] = (p[0], p[1], p[2] + 1)
    else:
        raise SystemExit('unknown op')
    write(argv[1], cv, comment='perturbed from %s: %s' % (argv[0], ' '.join(argv[2:])))
    print('wrote %s: total %.9f' % (argv[1], float(total(cv))))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    cmd, argv = sys.argv[1], sys.argv[2:]
    if cmd == 'sound':
        return cmd_sound(argv)
    if cmd == 'tight':
        return cmd_tight(argv)
    if cmd == 'sound0':
        return cmd_sound0(argv)
    if cmd == 'mu':
        cv = load(argv[0])
        x, y, u = (Fr(t) for t in argv[1:4])
        m = mu_exact(cv, x, y, u)
        print('mu = %s = %.12f  admissible %s' % (m, float(m), admissible(cv, x, y, u)))
        return 0
    if cmd == 'toy':
        return cmd_toy(argv)
    if cmd == 'perturb':
        return cmd_perturb(argv)
    if cmd == 'uncert':
        return cmd_uncert(argv)
    if cmd == 'fmin':
        cv = load(argv[0])
        print(fmin(cv, float(opt(argv, '--pitch', 0.05)), int(opt(argv, '--nth', 46))))
        return 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main())
