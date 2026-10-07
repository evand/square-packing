#!/usr/bin/env python3
"""Independent data-integrity check: family file <-> box file <-> exact.json.  Exact Fractions only."""
import json, hashlib, sys
from fractions import Fraction as F
from collections import Counter
B = '/home/evand/math/square-packing/public/s12/search/qx2_data/'
k = 7; R = 2; w = 2; a = F(9, 5)

prof = []; mod = []
for line in open(B + 'L4_k02_family.txt'):
    if line.startswith('#') or not line.strip(): continue
    t = line.split(); prof.append(t) if t[0] in 'hv' else mod.append(t)
prof = [(t[0],) + tuple(F(x) for x in t[1:]) for t in prof]
mod = [(t[0],) + tuple(F(x) for x in t[1:]) for t in mod]
print(len(prof), 'profile pieces;', len(mod), 'module pieces')
allm = [p[-1] for p in prof] + [p[-1] for p in mod]
assert all(m > 0 for m in allm), 'nonpositive mass'
# supports inside expected regions
for p in prof:
    if p[0] == 'h': _, y, p0, p1, m = p; assert 0 <= p0 < p1 <= 1 and 0 <= y <= w
    else: _, ph, y0, y1, m = p; assert 0 <= ph < 1 and 0 <= y0 < y1 <= w
for p in mod:
    _, c, t0, t1, m = p; assert 0 <= c <= R and 0 <= t0 < t1 <= R
# mirror symmetry of pi: x -> -x mod 1
def mir(p):
    if p[0] == 'h': return ('h', p[1], 1 - p[3], 1 - p[2], p[4])
    return ('v', (-p[1]) % 1, p[2], p[3], p[4])
print('profile mirror symmetric:', Counter(prof) == Counter(mir(p) for p in prof))
dg = lambda p: ('V' if p[0] == 'H' else 'H',) + p[1:]
print('module diagonal symmetric:', Counter(mod) == Counter(dg(p) for p in mod))
period = sum(p[-1] for p in prof) + (w - a) * 1
sigma = w - period
nu = sum(p[-1] for p in mod) + (R - a) ** 2
mv = sum(p[-1] for p in prof if p[0] == 'v' and p[1] == 0)
D = R * R - nu - mv
Dst = F(423621306389, 500000000000)
print('sigma =', sigma, '; D =', D, '; equals stated:', D == Dst, '; D > 3/4:', D > F(3, 4))
# exact.json
ej = json.load(open(B + 'L4_k02_exact.json'))
print('exact.json D matches:', F(ej['D']) == D, 'sigma', ej['sigma'])
xs = Counter(F(v) for v in ej['x'] if F(v) != 0)
fm = Counter(allm)
print('exact.json nonzero x values:', len(xs), 'distinct', sum(xs.values()), 'total; family masses distinct', len(fm),
      '; same value set:', set(xs) == set(fm), '; x all >= 0:', all(F(v) >= 0 for v in ej['x']))
# ---- build mu_7 myself
segs = []
def seg(x0, y0, x1, y1, m):
    P, Q = (x0, y0), (x1, y1)
    segs.append((min(P, Q), max(P, Q), m))
refl = [lambda x, y: (x, y), lambda x, y: (k - x, y), lambda x, y: (x, k - y), lambda x, y: (k - x, k - y)]
for (o, c, t0, t1, m) in mod:
    A, Bp = ((t0, c), (t1, c)) if o == 'H' else ((c, t0), (c, t1))
    for f in refl: seg(*f(*A), *f(*Bp), m)
rot = [lambda x, y: (x, y), lambda x, y: (k - y, x), lambda x, y: (k - x, k - y), lambda x, y: (y, k - x)]
for j in range(-1, k + 1):
    for p in prof:
        if p[0] == 'h': A, Bp = (j + p[2], p[1]), (j + p[3], p[1])
        else: A, Bp = (j + p[1], p[2]), (j + p[1], p[3])
        if min(A[0], Bp[0]) >= R and max(A[0], Bp[0]) <= k - R:
            for f in rot: seg(*f(*A), *f(*Bp), p[-1])
leb = (k - 2 * a) ** 2
tot = sum(s[2] for s in segs) + leb
print('my mu_7:', len(segs), 'segments; total', tot, '= 49 - 4D:', tot == 49 - 4 * D, '; < 46:', tot < 46)
# ---- box file
tok = []
for line in open(B + 'L4_k02_box7.txt'): tok += line.split('#', 1)[0].split()
it = iter(tok); assert next(it) == 'mixed' and next(it) == '1'
sn, sd, Dd, W = (int(next(it)) for _ in range(4)); assert F(sn, sd) == k
npt = int(next(it)); pts = [tuple(int(next(it)) for _ in range(3)) for _ in range(npt)]
ns = int(next(it)); bs = [tuple(int(next(it)) for _ in range(5)) for _ in range(ns)]
npg = int(next(it)); pgs = []
for _ in range(npg):
    kk = int(next(it)); wm = int(next(it)); V = [(F(int(next(it)), Dd), F(int(next(it)), Dd)) for _ in range(kk)]; pgs.append((wm, V))
assert list(it) == [], 'trailing tokens'
print('box: points', npt, 'segments', ns, 'polygons', npg, [ (F(wm, W), V) for wm, V in pgs])
bsegs = []
for (X0, Y0, X1, Y1, wm) in bs:
    assert wm >= 0
    P, Q = (F(X0, Dd), F(Y0, Dd)), (F(X1, Dd), F(Y1, Dd)); assert P != Q
    for c in P + Q: assert 0 <= c <= k
    bsegs.append((min(P, Q), max(P, Q), F(wm, W)))
print('box segments == my segments (multiset):', Counter(bsegs) == Counter(segs))
assert npg == 1 and pgs[0][1] == [(a, a), (k - a, a), (k - a, k - a), (a, k - a)] and F(pgs[0][0], W) == leb
btot = sum(s[2] for s in bsegs) + sum(F(wm, W) for wm, _ in pgs) + sum(F(p[2], W) for p in pts)
print('box total', btot, '== 49-4D:', btot == 49 - 4 * D)
# D4 symmetry of box cover
maps = [lambda x, y: (x, y), lambda x, y: (k - y, x), lambda x, y: (k - x, k - y), lambda x, y: (y, k - x),
        lambda x, y: (k - x, y), lambda x, y: (x, k - y), lambda x, y: (y, x), lambda x, y: (k - y, k - x)]
C0 = Counter(bsegs)
for i, f in enumerate(maps):
    Ci = Counter((min(f(*P), f(*Q)), max(f(*P), f(*Q)), m) for P, Q, m in bsegs)
    assert Ci == C0, i
print('box D4 invariant: yes (8 maps); polygon [a,k-a]^2 is D4 invariant')
for fn, pre in [('L4_k02_family.txt', '3f4771d3'), ('L4_k02_box7.txt', 'c0a67509'), ('L4_k02_exact.json', '8931e71d')]:
    h = hashlib.sha256(open(B + fn, 'rb').read()).hexdigest(); print(fn, h, h.startswith(pre))
