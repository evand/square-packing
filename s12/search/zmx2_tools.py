#!/usr/bin/env python3
"""zmx2_tools.py -- test harness for the independent mixed-cover checker `zmx2` (search/ZMX2.md).

Written independently of search/zm_mixed.py (never read).  Own parser, own exact evaluator.

  load / write                 mixed v1 files (tasks/line-cover/FORMAT.md), exact integers
  mu_exact(cv, x, y, u)        exact mu(Q) at a rational pose (Fractions; any segment direction)
  admissible(cv, x, y, u)      exact admissibility of the pose (Q inside [0,s]^2)

Subcommands:
  sound FILE [--n N] [--poses P] [--seed S] [--refl]
        soundness harness: random small pose boxes (biased to germs, walls, grid lines), the
        bound printed by `zmx2 boxes`, and exact mu at P random admissible rational poses of
        each box (plus corners); FAIL if any mu < bound.
  mu FILE x y u                exact mu at a rational pose (x, y, u as p/q)
  toy OUT --s S --kind K [...] toy covers for the certify / reject tests (see ZMX2.md sec 8)
  fmin FILE [--pitch P] [--ubins N]   float minimum over a grid of admissible poses (all angles)
  perturb IN OUT --op ...      rejection-test covers (drop mass in a region; break symmetry)
"""
import math
import os
import random
import subprocess
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ZMX2 = os.path.join(HERE, '..', 'verify2', 'target', 'release', 'zmx2')


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
    if mixed:
        segs = [(nx(), nx(), nx(), nx(), nx()) for _ in range(nx())]
        assert nx() == 0, 'polygons unsupported'
    assert i == len(toks)
    return dict(s=Fr(s_num, s_den), s_num=s_num, s_den=s_den, D=D, W=W, points=pts, segments=segs)


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
        f.write('0\n')


def total(cv):
    return Fr(sum(p[2] for p in cv['points']) + sum(s[4] for s in cv['segments']), cv['W'])


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
    return tot / W


def reflect_y(cv):
    S = cv['s'] * cv['D']
    assert S.denominator == 1
    S = S.numerator
    out = dict(cv)
    out['points'] = [(X, S - Y, w) for (X, Y, w) in cv['points']]
    out['segments'] = [(a, S - b, c, S - d, w) for (a, b, c, d, w) in cv['segments']]
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


def cmd_sound(argv):
    path = argv[0]
    n = int(opt(argv, '--n', 200))
    npose = int(opt(argv, '--poses', 25))
    seed = int(opt(argv, '--seed', 1))
    refl = '--refl' in argv
    rng = random.Random(seed)
    cv = load(path)
    cve = reflect_y(cv) if refl else cv
    boxes = [rand_box(rng, cv['s']) for _ in range(n)]
    inp = '\n'.join(','.join('%d/%d' % (f.numerator, f.denominator) for f in bx) for bx in boxes) + '\n'
    cmd = [ZMX2, 'boxes', path] + (['--refl'] if refl else [])
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
    --op bump (add 1 to the first point's weight: breaks D4)."""
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
    if cmd == 'fmin':
        cv = load(argv[0])
        print(fmin(cv, float(opt(argv, '--pitch', 0.05)), int(opt(argv, '--nth', 46))))
        return 0
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main())
