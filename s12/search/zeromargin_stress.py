#!/usr/bin/env python3
"""Float stress test of a zeromargin.py leaf dump: sample poses inside every certified leaf and
check its claimed witness directly (ADM/CORE/P1/MIX: some listed point is in the square; TRI: some
vertex is).  Also random tests of the four primitives themselves.  Independent of the exact code
paths (plain floats, no Fractions, no shared helper).

usage: python3 search/zeromargin_stress.py LEAVES.txt [samples_per_leaf] [--cert CERT.txt]

With no --cert the point set is Friedman's 14 (rung 1).  For a weighted cover pass the same
certificate file that produced the dump.
"""
import sys, math, random, ast
from fractions import Fraction as F
random.seed(1)

def inside(p, c, th, tol=1e-9):
    dx, dy = p[0] - c[0], p[1] - c[1]
    x = dx * math.cos(th) + dy * math.sin(th); y = -dx * math.sin(th) + dy * math.cos(th)
    return abs(x) <= 0.5 + tol and abs(y) <= 0.5 + tol

FRIED = [(1, 1), (1.6, 1), (2.4, 1), (3, 1), (1, 1.8), (2, 1.8), (3, 1.8), (1, 2.2), (2, 2.2),
         (3, 2.2), (1, 3), (1.6, 3), (2.4, 3), (3, 3)]

def read_cert(path):
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    P = []; WT = []
    for i in range(n):
        X, Y, w = map(int, tok[5 + 3 * i: 8 + 3 * i])
        P.append((X / D, Y / D)); WT.append(w / W)
    return sn / sd, P, WT

args = [a for a in sys.argv[1:] if not a.startswith('--')]
path = args[0]; nsamp = int(args[1]) if len(args) > 1 else 20
m = 4.0; P = FRIED; WT = [1.0] * len(FRIED)
if '--cert' in sys.argv:
    m, P, WT = read_cert(sys.argv[sys.argv.index('--cert') + 1])

bad = 0; n = {'ADM': 0, 'CORE': 0, 'P1': 0, 'MIX': 0, 'CHAIN': 0, 'TRI': 0, 'EMPTY': 0}
for line in open(path):
    if line.startswith('#'): continue
    boxs, kind, wit = [t.strip() for t in line.split(';')]
    box = [float(F(t)) for t in boxs.split()]
    n[kind] += 1
    cx0, cx1, cy0, cy1, u0, u1 = box
    th0, th1 = 2 * math.atan(u0), 2 * math.atan(u1)
    wit = eval(wit, {'Fraction': F, 'F': F}) if kind != 'EMPTY' else None
    for _ in range(nsamp):
        cx = random.choice([cx0, cx1, random.uniform(cx0, cx1)])
        cy = random.choice([cy0, cy1, random.uniform(cy0, cy1)])
        th = random.choice([th0, th1, random.uniform(th0, th1)])
        w = math.cos(th) + math.sin(th)
        adm = (w / 2 - 1e-12 <= cx <= m - w / 2 + 1e-12) and (w / 2 - 1e-12 <= cy <= m - w / 2 + 1e-12)
        if kind == 'EMPTY':
            if adm: bad += 1; print("EMPTY box has admissible pose", box, cx, cy, th)
            continue
        if kind in ('P1', 'ADM', 'MIX', 'CHAIN') and not adm: continue  # admissible poses only
        if kind == 'CORE':
            # CORE claims every pose of the box clipped to the widest admissible range of the bin
            lo = min(math.cos(th0) + math.sin(th0), math.cos(th1) + math.sin(th1)) / 2
            if not (lo - 1e-12 <= cx <= m - lo + 1e-12 and lo - 1e-12 <= cy <= m - lo + 1e-12): continue
        if kind == 'CHAIN':
            # the witness list of a CHAIN leaf is the UNION over its regions; the claim it makes
            # is that at every admissible pose of the box the captured weight drawn from that list
            # reaches 1 (a different subset does it in each region).  That is what is checked here,
            # with no knowledge of the regions themselves.
            ok = sum(WT[k] for k in wit if inside(P[k], (cx, cy), th)) >= 1.0 - 1e-9
        elif kind in ('CORE', 'P1', 'ADM', 'MIX'):
            ok = any(inside(P[k], (cx, cy), th) for k in wit)
        else:
            ok = any(inside((float(v[0]), float(v[1])), (cx, cy), th) for v in wit)
        if not ok:
            bad += 1
            if bad <= 10: print("FAIL", kind, box, (cx, cy, math.degrees(th)), wit)
print("leaves:", n, "samples per leaf:", nsamp, "failures:", bad)

# --- the primitives themselves on random inputs ----------------------------------------------
fails = 0; counts = {'P1': 0, 'core': 0, 'tri': 0, 'adm': 0}
for _ in range(200000):
    th = random.uniform(0, math.pi / 2); w = math.cos(th) + math.sin(th)
    # P1: |c-p|_inf <= 1 - w/2  =>  p in S
    p = (random.uniform(0, 4), random.uniform(0, 4)); t = 1 - w / 2
    c = (p[0] + random.uniform(-t, t), p[1] + random.uniform(-t, t))
    counts['P1'] += 1
    if not inside(p, c, th): fails += 1; print("P1 fail", p, c, th)
    # core lemma: x in R_t0 Q and R_t1 Q, and (dir(x) mod 90 not in [t0,t1] or |x| <= 1/2)  =>  x in R_t Q on the bin
    t0 = random.uniform(0, math.pi / 4); t1 = t0 + random.uniform(0, math.pi / 4)
    x = (random.uniform(-0.8, 0.8), random.uniform(-0.8, 0.8))
    if inside(x, (0, 0), t0) and inside(x, (0, 0), t1):
        phi = math.atan2(x[1], x[0]) % (math.pi / 2); rel = (phi - t0) % (math.pi / 2)
        sector = rel <= t1 - t0 + 1e-12
        if (not sector) or math.hypot(*x) <= 0.5:
            counts['core'] += 1
            for tt in [t0 + (t1 - t0) * k / 20 for k in range(21)]:
                if not inside(x, (0, 0), tt): fails += 1; print("core fail", x, t0, t1, tt); break
    # triangle lemma: sides <= 1, centre in T  =>  some vertex in S
    A = (random.uniform(0, 1), random.uniform(0, 1))
    B = (A[0] + random.uniform(-1, 1), A[1] + random.uniform(-1, 1))
    C = (A[0] + random.uniform(-1, 1), A[1] + random.uniform(-1, 1))
    if all(math.dist(u, v) <= 1 for u, v in ((A, B), (B, C), (A, C))):
        a, b = random.random(), random.random()
        if a + b > 1: a, b = 1 - a, 1 - b
        cc = (A[0] + a * (B[0] - A[0]) + b * (C[0] - A[0]), A[1] + a * (B[1] - A[1]) + b * (C[1] - A[1]))
        counts['tri'] += 1
        if not any(inside(v, cc, th) for v in (A, B, C)): fails += 1; print("TRI fail", A, B, C, cc, th)
print("primitive random tests:", counts, "failures", fails)

# --- ADM (search/RUNG2.md Lemma A+B) on random boxes, checked against brute-force sampling -----
# Lemma A/B claim: if the four corner inequalities hold for every theta of the bin -- with the
# centre bounds A_x = max(cx0, w/2), B_x = min(cx1, m - w/2) etc. -- then p is in Q(c, theta) at
# EVERY admissible pose of the box.  Here the four inequalities are evaluated directly in floats on
# a fine theta grid (no polynomials, no Bernstein) and the conclusion is tested by sampling poses.
mfl = 4.0; afail = 0; atests = 0; atrue = 0
rng = random.Random(7)
for _ in range(40000):
    cx0 = rng.uniform(0, mfl); cx1 = cx0 + rng.choice([0.0, 0.001, 0.01, 0.1]) * rng.random()
    cy0 = rng.uniform(0, mfl); cy1 = cy0 + rng.choice([0.0, 0.001, 0.01, 0.1]) * rng.random()
    t0 = rng.uniform(0, math.pi / 2); t1 = min(math.pi / 2, t0 + rng.choice([1e-4, 1e-2, 0.1, 0.4]) * rng.random())
    px = rng.uniform(0, mfl); py = rng.uniform(0, mfl)
    TH = [t0 + (t1 - t0) * k / 200 for k in range(201)]
    ok = True
    for th in TH:
        w = math.cos(th) + math.sin(th); c, s = math.cos(th), math.sin(th)
        Ax = max(cx0, w / 2); Bx = min(cx1, mfl - w / 2)
        Ay = max(cy0, w / 2); By = min(cy1, mfl - w / 2)
        if (px - Ax) * c + (py - Ay) * s > 0.5 + 1e-12: ok = False; break
        if (px - Bx) * c + (py - By) * s < -0.5 - 1e-12: ok = False; break
        if -(px - Bx) * s + (py - Ay) * c > 0.5 + 1e-12: ok = False; break
        if -(px - Ax) * s + (py - By) * c < -0.5 - 1e-12: ok = False; break
    if not ok: continue
    atrue += 1
    for _ in range(40):                      # the claim: p is captured at every admissible pose
        th = rng.choice([t0, t1, rng.uniform(t0, t1)])
        w = math.cos(th) + math.sin(th)
        lo, hi = w / 2, mfl - w / 2
        a, b = max(cx0, lo), min(cx1, hi)
        cc, dd = max(cy0, lo), min(cy1, hi)
        if a > b or cc > dd: continue
        cxs = rng.choice([a, b, rng.uniform(a, b)]); cys = rng.choice([cc, dd, rng.uniform(cc, dd)])
        atests += 1
        if not inside((px, py), (cxs, cys), th, tol=1e-9):
            afail += 1
            if afail <= 5: print("ADM fail", (cx0, cx1, cy0, cy1, t0, t1), (px, py), (cxs, cys, th))
print(f"ADM lemma random tests: {atrue} certifying (point, box) pairs, {atests} pose samples, failures {afail}")

# --- the CHAIN algebra (RUNG2.md Lemmas E-G): the violation polynomials and their exact maximum ---
# G_{p,k} > 0 must mean "p violates the k-th containment inequality", and _gmax must be an upper
# bound for any nonnegative combination of them over the whole pose box.  Both are re-derived here
# from the geometry, in floats, with no reference to zeromargin.py's Fraction code path.
import importlib.util as _ilu
_sp = _ilu.spec_from_file_location('zm', __file__.rsplit('/', 1)[0] + '/zeromargin.py')
_zm = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_zm)
from fractions import Fraction as _F

def gfloat(kind, px, py, cx, cy, u):
    a, b = px - cx, py - cy
    C, S, N = 1 - u * u, 2 * u, 1 + u * u
    if kind == 0: return 2 * a * C + 2 * b * S - N
    if kind == 1: return -2 * a * C - 2 * b * S - N
    if kind == 2: return -2 * a * S + 2 * b * C - N
    return 2 * a * S - 2 * b * C - N

rng2 = random.Random(11)
sgn_fail = 0; sgn_n = 0
for _ in range(200000):
    px, py = rng2.uniform(0, 4), rng2.uniform(0, 4)
    cx, cy = rng2.uniform(0, 4), rng2.uniform(0, 4)
    u = rng2.uniform(0, 1); th = 2 * math.atan(u)
    g = [gfloat(k, px, py, cx, cy, u) for k in range(4)]
    ins = inside((px, py), (cx, cy), th, tol=0.0)
    if max(g) > 1e-9 and ins: sgn_fail += 1
    if max(g) < -1e-9 and not ins: sgn_fail += 1
    sgn_n += 1
print(f"CHAIN violation polynomials: {sgn_n} random (point, pose) pairs, "
      f"sign disagreements with the geometry: {sgn_fail}")

chk = _zm.Checker(4, [(0, 0)], None)
gm_fail = 0; gm_n = 0
for _ in range(20000):
    pts = [(_F(rng2.randrange(0, 4001), 1000), _F(rng2.randrange(0, 4001), 1000)) for _ in range(3)]
    chk.P = pts
    cx0 = _F(rng2.randrange(0, 4001), 1000); cx1 = cx0 + _F(rng2.randrange(0, 300), 1000)
    cy0 = _F(rng2.randrange(0, 4001), 1000); cy1 = cy0 + _F(rng2.randrange(0, 300), 1000)
    u0 = _F(rng2.randrange(0, 1000), 1000); u1 = u0 + _F(rng2.randrange(0, 200), 1000)
    box = (cx0, cx1, cy0, cy1, u0, u1)
    terms = [(_F(rng2.choice([1, 1, 2, 1]), rng2.choice([1, 2, 1])), i, rng2.randrange(4))
             for i in range(rng2.randrange(1, 3))]
    bnd = float(chk._gmax(terms, box))
    gm_n += 1
    for _ in range(30):
        cx = rng2.uniform(float(cx0), float(cx1)); cy = rng2.uniform(float(cy0), float(cy1))
        u = rng2.uniform(float(u0), float(u1))
        v = sum(float(lam) * gfloat(kd, float(pts[i][0]), float(pts[i][1]), cx, cy, u)
                for lam, i, kd in terms)
        if v > bnd + 1e-9:
            gm_fail += 1
            if gm_fail <= 5: print("  _gmax fail", box, terms, v, bnd)
            break
print(f"CHAIN _gmax enclosure: {gm_n} random (box, combination) pairs x 30 interior poses, "
      f"violations of the exact bound: {gm_fail}")
