"""End-to-end test of the MATHEMATICS of Lemma E at a fixed rational u (independent implementation).

For a random box of centres, a random set of axis lines with random piecewise-uniform profiles, and the Lebesgue
square U = [a, b]^2, build the arrangement L(u) exactly as the lemma statement prescribes, the cell function
F_K = sum_l mu_l + Phi_U^K (E' std/tan, E'' corner3/xcut with McCormick, both anchors), and check exactly:
  T1  mu(Q(c)) >= F_K(c)                                (validity of the minorant)
  T2  F_K(c) >= min over the vertices of the closed cell K of F_K(v)   (the lemma's conclusion, cell-wise)
  T3  F_K(midpoint) >= average  for pairs in the same open cell       (concavity)
"""
import random, sys, itertools, os
MUT = os.environ.get('MUT', '')
from fractions import Fraction as F
from geom import *

a, b = F(9, 5), F(26, 5)
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
ntrials = int(sys.argv[2]) if len(sys.argv) > 2 else 50
rng = random.Random(seed)


def opts(orient, p, u):
    """(type, F, a0, ax, ay): F t = a0 + ax cx + ay cy  (QUADRANT_EXACT 4.3, re-derived by hand)"""
    C, S, N = 1 - u * u, 2 * u, 1 + u * u
    if orient == 'H':
        eta = p
        return [('up', C, N / 2 - eta * S, C, S), ('lo', C, -N / 2 - eta * S, C, S),
                ('lo', S, eta * C - N / 2, S, -C), ('up', S, eta * C + N / 2, S, -C)]
    xi = p
    return [('up', S, N / 2 - xi * C, C, S), ('lo', S, -N / 2 - xi * C, C, S),
            ('up', C, N / 2 + xi * S, -S, C), ('lo', C, -N / 2 + xi * S, -S, C)]


def ghat(d, s, c):
    if d <= 0: return F(0)
    if d <= s: return d * d / (2 * s * c)
    return (d - s / 2) / c


def run_trial(u, box, lines, anchor):
    cx0, cx1, cy0, cy1 = box
    C, S, N = 1 - u * u, 2 * u, 1 + u * u
    c, s = trig(u); w = c + s
    corners = [(x, y) for x in (cx0, cx1) for y in (cy0, cy1)]
    Ls = []            # (A, B, K): A cx + B cy + K = 0
    Ls += [(F(1), F(0), -cx0), (F(1), F(0), -cx1), (F(0), F(1), -cy0), (F(0), F(1), -cy1)]

    zholder = []

    def add(l, isz=False):
        vals = [l[0] * x + l[1] * y + l[2] for x, y in corners]
        if all(v > 0 for v in vals) or all(v < 0 for v in vals): return
        if l[0] == 0 and l[1] == 0: return
        if isz: zholder.append(len(Ls))
        Ls.append(l)
    for (orient, pos, prof) in lines:
        op = opts(orient, pos, u)
        for (typ, Fk, a0, ax, ay) in op:
            for i, bp in enumerate(prof.b):
                left, right = prof.rho[i], prof.rho[i + 1]
                isbad = (typ == 'up' and right > left) or (typ == 'lo' and right < left)
                if MUT == 'nobad': isbad = False
                if MUT == 'flip': isbad = not isbad
                if isbad:
                    add((ax, ay, a0 - Fk * bp))
        for o1 in op:
            for o2 in op:
                if o1[0] == 'up' and o2[0] == 'lo' and MUT != 'nopair':
                    # t1 = t2:  (a0_1 + ...)/F1 = (a0_2 + ...)/F2
                    add((o1[3] * o2[1] - o2[3] * o1[1], o1[4] * o2[1] - o2[4] * o1[1], o1[2] * o2[1] - o2[2] * o1[1]))
    # ---- U regime (at this u; the box is the whole R(u))
    if cx1 + w / 2 <= a or cy1 + w / 2 <= a: reg = 'none'
    else:
        assert cx1 + w / 2 <= b and cy1 + w / 2 <= b
        caps = []
        for axis, lo, hi in ((0, cx0, cx1), (1, cy0, cy1)):
            dmax = a - lo + w / 2
            if dmax <= 0: continue
            dmin = a - hi + w / 2
            if dmin >= c: caps.append((axis, 'tan', min(F(1), (dmin + dmax) / 2)))
            else: caps.append((axis, 'std', None))
        reg = caps
        if len(caps) == 2 and all(m == 'std' for _, m, _ in caps):
            xBL_max = cx1 - (c - s) / 2; xBR_min = cx0 + w / 2; yTL_min = cy0 + (c - s) / 2
            Xaa = ((a - cx0) * C + (a - cy0) * S) / N
            if xBL_max <= a <= xBR_min and yTL_min >= a and Xaa <= HALF: reg = caps + [(2, 'corner3', None)]
            elif cx0 - (c - s) / 2 >= a and yTL_min >= a: reg = [caps[1], (0, 'xcut', None)]
    if reg != 'none':
        for axis, mode, _ in reg:
            if mode in ('std', 'xcut') and axis in (0, 1):
                for dd in (F(0), s):          # d = a - coord + w/2 = dd
                    add((F(1), F(0), -(a + w / 2 - dd)) if axis == 0 else (F(0), F(1), -(a + w / 2 - dd)))
            if mode == 'corner3':
                # z = beta C - alpha S, alpha = (a - x_BL) N, x_BL = cx - (c - s)/2 ; beta = (a - y_BL) N, y_BL = cy - w/2
                # z = N [ C (a - cy + w/2) - S (a - cx + (c-s)/2) ]
                (MUT != 'noz') and add((S * N, -C * N, N * (C * (a + w / 2) - S * (a + (c - s) / 2))), True)
            if mode == 'xcut':
                # z' = alpha' C + beta'' S, alpha' = (a - x_TL) N, x_TL = cx - w/2 ; beta'' = (a - y_TL) N, y_TL = cy + (c-s)/2
                (MUT != 'noz') and add((-C * N, -S * N, N * (C * (a + w / 2) + S * (a - (c - s) / 2))), True)
    modes = [m for _, m, _ in reg] if reg != 'none' else []

    def zval(cx_, cy_):
        if 'corner3' in modes:
            al = (a - cx_ + (c - s) / 2) * N; be = (a - cy_ + w / 2) * N
            return be * C - al * S
        if 'xcut' in modes:
            al = (a - cx_ + w / 2) * N; be = (a - cy_ - (c - s) / 2) * N
            return al * C + be * S
        return None

    def mcq(cx_, cy_, kind):
        if kind == 'corner3':
            fa = lambda x: (a - x + (c - s) / 2) * N; fb = lambda y: (a - y + w / 2) * N
        else:
            fa = lambda x: (a - x + w / 2) * N; fb = lambda y: (a - y - (c - s) / 2) * N
        al, be = fa(cx_), fb(cy_)
        A_, B_ = (fa(cx0), fb(cy0)) if anchor == 'hi' else (fa(cx1), fb(cy1))
        bst = fb((cy0 + cy1) / 2)
        ab = A_ * be + B_ * al - A_ * B_
        return (2 * C * ab + S * (2 * bst * be - bst * bst) - S * al * al) / (2 * C * N * N)

    def phiU(cx_, cy_, side):
        """cell function of the Lebesgue part; side = sign of z on the cell (+1/-1), 0 = min of both"""
        if reg == 'none': return F(0)
        if side == 0 and ('corner3' in modes or 'xcut' in modes):
            return min(phiU(cx_, cy_, 1), phiU(cx_, cy_, -1))
        dx = a - cx_ + w / 2; dy = a - cy_ + w / 2
        tot = F(0); ncap = 0
        for axis, mode, ds in reg:
            d = dx if axis == 0 else dy
            if mode == 'std': tot += 1 - ghat(d, s, c); ncap += 1
            elif mode == 'tan':
                k = s + c - ds; tot += k * k / (2 * s * c) - k / (s * c) * (d - ds); ncap += 1
        tot -= (ncap - 1)
        if 'corner3' in modes:
            tot += ghat(dy, s, c) if side < 0 else mcq(cx_, cy_, 'corner3')
        if 'xcut' in modes:
            tot += -ghat(dx, s, c) if side < 0 else mcq(cx_, cy_, 'xcut')
        return tot

    def lines_mass(cx_, cy_):
        P = square(cx_, cy_, u)
        return sum((mass_line(P, o, p, pr) for (o, p, pr) in lines), F(0))

    def true_mass(cx_, cy_):
        P = square(cx_, cy_, u)
        lam = area(clip(clip(clip(clip(P, F(-1), F(0), a), F(0), F(-1), a), F(1), F(0), -b), F(0), F(1), -b))
        return lines_mass(cx_, cy_) + lam

    def sgnvec(p):
        return tuple((1 if v > 0 else (-1 if v < 0 else 0)) for v in (l[0] * p[0] + l[1] * p[1] + l[2] for l in Ls))

    # vertices
    V = []
    for l1, l2 in itertools.combinations(Ls, 2):
        D = l1[0] * l2[1] - l2[0] * l1[1]
        if D == 0: continue
        x = (l1[1] * l2[2] - l2[1] * l1[2]) / D; y = (l2[0] * l1[2] - l1[0] * l2[2]) / D
        if cx0 <= x <= cx1 and cy0 <= y <= cy1: V.append((x, y))
    V = list(set(V))
    Vs = [(v, sgnvec(v), lines_mass(*v)) for v in V]
    zi = zholder[0] if zholder else None
    zconst = None
    if zi is None and ('corner3' in modes or 'xcut' in modes):
        zc_ = zval(cx0, cy0); zconst = 1 if zc_ > 0 else -1   # z has constant strict sign on the box (line not crossing)
    fails = []
    npts = 0
    for it in range(12):
        p = (cx0 + (cx1 - cx0) * F(rng.randint(1, 999), 1000), cy0 + (cy1 - cy0) * F(rng.randint(1, 999), 1000))
        sv = sgnvec(p)
        if 0 in sv: continue
        npts += 1
        side = sv[zi] if zi is not None else (zconst or 1)
        if zi is not None and zval(*p) is not None:
            zs = zval(*p); assert (zs > 0) == (side > 0) or zs == 0
        Fp = lines_mass(*p) + phiU(*p, side)
        tm = true_mass(*p)
        if tm < Fp: fails.append(('T1', p, float(tm - Fp)))
        # T1 on the z line's other side too is not needed; T2:
        vK = [(v, lm) for (v, svv, lm) in Vs if all(x == 0 or x == y for x, y in zip(svv, sv))]
        if not vK: fails.append(('T2-novertex', p)); continue
        mv = min(lm + phiU(v[0], v[1], side) for v, lm in vK)
        if Fp < mv: fails.append(('T2', p, float(Fp - mv)))
        # T3 midpoint concavity
        for jt in range(3):
            q = (p[0] + (cx1 - cx0) * F(rng.randint(-50, 50), 1000), p[1] + (cy1 - cy0) * F(rng.randint(-50, 50), 1000))
            if not (cx0 < q[0] < cx1 and cy0 < q[1] < cy1) or sgnvec(q) != sv: continue
            m = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
            Fq = lines_mass(*q) + phiU(*q, side); Fm = lines_mass(*m) + phiU(*m, side)
            if Fm < (Fp + Fq) / 2: fails.append(('T3', p, q, float(Fm - (Fp + Fq) / 2)))
    # also T1 at the vertices with the 'min' rule (side 0 at z = 0 vertices)
    for v, svv, lm in Vs:
        side = 0 if (zi is not None and svv[zi] == 0) else (svv[zi] if zi is not None else (zconst or 1))
        if true_mass(*v) < lm + phiU(v[0], v[1], side): fails.append(('T1v', v))
    return fails, npts, len(Ls), len(V), modes


def rand_setup():
    k = rng.random()
    if k < 0.2: u = F(rng.randint(1, 2000), 10 ** 6)
    elif k < 0.35: u = F(41421356, 10 ** 8) - F(rng.randint(0, 10 ** 5), 10 ** 8)
    else: u = F(rng.randint(1, 4142), 10 ** 4)
    h = F(rng.choice([1, 2, 5, 10, 20, 40]), 100)
    where = rng.random()
    if where < 0.5:     # near the corner of U
        cx0 = a + F(rng.randint(-800, 600), 1000); cy0 = a + F(rng.randint(-800, 600), 1000)
    elif where < 0.8:   # near the wall x = a
        cx0 = a + F(rng.randint(-800, 600), 1000); cy0 = F(rng.randint(2600, 4000), 1000)
    else:
        cx0 = F(rng.randint(500, 4000), 1000); cy0 = F(rng.randint(500, 1500), 1000)
    box = (cx0, cx0 + h * F(rng.randint(5, 10), 10), cy0, cy0 + h * F(rng.randint(5, 10), 10))
    lines = []
    for _ in range(rng.randint(2, 5)):
        o = rng.choice('HV')
        if o == 'H': pos = F(rng.randint(int((box[2] - 1) * 20), int((box[3] + 1) * 20)), 20); lo, hi = box[0] - 1, box[1] + 1
        else: pos = F(rng.randint(int((box[0] - 1) * 20), int((box[1] + 1) * 20)), 20); lo, hi = box[2] - 1, box[3] + 1
        lines.append((o, pos, rand_prof(rng, lo, hi, rng.randint(2, 6))))
    return u, box, lines


tot = dict(trials=0, pts=0, fails=0)
modecount = {}
for t in range(ntrials):
    u, box, lines = rand_setup()
    for anchor in ('hi', 'lo'):
        fails, npts, nL, nV, modes = run_trial(u, box, lines, anchor)
        tot['trials'] += 1; tot['pts'] += npts; tot['fails'] += len(fails)
        modecount[tuple(modes)] = modecount.get(tuple(modes), 0) + 1
        for f_ in fails[:3]: print("FAIL", t, anchor, float(u), [float(x) for x in box], f_)
print(tot, modecount)
