#!/usr/bin/env python3
"""Float stress test of a zeromargin.py leaf dump: sample poses inside every certified leaf and
check its claimed witness directly (CORE/P1: the point is in the square; TRI: some vertex is).
Also random tests of the three primitives themselves.  Independent of the exact code paths.

usage: python3 search/zeromargin_stress.py runs/zeromargin_friedman14_leaves.txt [samples_per_leaf]
"""
import sys, math, random, ast
from fractions import Fraction as F
random.seed(1)

def inside(p, c, th, tol=1e-9):
    dx, dy = p[0] - c[0], p[1] - c[1]
    x = dx * math.cos(th) + dy * math.sin(th); y = -dx * math.sin(th) + dy * math.cos(th)
    return abs(x) <= 0.5 + tol and abs(y) <= 0.5 + tol

P = [(1, 1), (1.6, 1), (2.4, 1), (3, 1), (1, 1.8), (2, 1.8), (3, 1.8), (1, 2.2), (2, 2.2), (3, 2.2),
     (1, 3), (1.6, 3), (2.4, 3), (3, 3)]
m = 4.0
path = sys.argv[1]; nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 20
bad = 0; n = {'CORE': 0, 'P1': 0, 'TRI': 0, 'EMPTY': 0}
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
        if kind == 'P1' and not adm: continue          # P1 only claims admissible poses
        if kind == 'CORE':
            # CORE claims every pose of the box clipped to the widest admissible range of the bin
            lo = min(math.cos(th0) + math.sin(th0), math.cos(th1) + math.sin(th1)) / 2
            if not (lo - 1e-12 <= cx <= m - lo + 1e-12 and lo - 1e-12 <= cy <= m - lo + 1e-12): continue
        if kind in ('CORE', 'P1'):
            ok = any(inside(P[k], (cx, cy), th) for k in wit)
        else:
            ok = any(inside((float(v[0]), float(v[1])), (cx, cy), th) for v in wit)
        if not ok:
            bad += 1
            if bad <= 10: print("FAIL", kind, box, (cx, cy, math.degrees(th)), wit)
print("leaves:", n, "samples per leaf:", nsamp, "failures:", bad)

# --- the primitives themselves on random inputs ----------------------------------------------
fails = 0; counts = {'P1': 0, 'core': 0, 'tri': 0}
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
