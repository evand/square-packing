# Independent Lemma Z enumeration (theta = 0) from the box file. Exact Fractions.
import sys
from fractions import Fraction as F
from collections import Counter
fn = sys.argv[1]
toks = []
for line in open(fn):
    line = line.split('#')[0]
    toks += line.split()
assert toks[0] == 'mixed' and toks[1] == '1'
it = iter(toks[2:])
nx = lambda: int(next(it))
s = F(nx(), nx()); D = nx(); W = nx()
npnt = nx(); assert npnt == 0
ns = nx()
segs = []
for _ in range(ns):
    X0, Y0, X1, Y1, w = [nx() for _ in range(5)]
    segs.append((F(X0, D), F(Y0, D), F(X1, D), F(Y1, D), F(w, W)))
npg = nx()
polys = []
for _ in range(npg):
    k = nx(); w = F(nx(), W)
    pts = [(F(nx(), D), F(nx(), D)) for _ in range(k)]
    polys.append((w, pts))
rest = list(it); assert rest == [], rest
print('s', s, 'segs', len(segs), 'polys', polys)
# all axis parallel?
H = [(min(x0, x1), max(x0, x1), y0, w) for x0, y0, x1, y1, w in segs if y0 == y1]
V = [(min(y0, y1), max(y0, y1), x0, w) for x0, y0, x1, y1, w in segs if x0 == x1]
assert len(H) + len(V) == len(segs)
assert all(p < q for p, q, _, _ in H + V)
(pw, ppts), = polys
xs = sorted(set(p[0] for p in ppts)); ys = sorted(set(p[1] for p in ppts))
assert len(xs) == 2 and len(ys) == 2 and len(ppts) == 4
a, b = xs; assert ys == xs
assert pw == (b - a) ** 2
print('U =', a, b, 'H', len(H), 'V', len(V))
tot = sum(w for *_, w in segs) + pw
print('total', tot, float(tot))
# closed axis-parallel square mass
def ov(p, q, x0, x1):
    return max(F(0), min(q, x1) - max(p, x0))
def mass(cx, cy):
    x0, x1, y0, y1 = cx - F(1, 2), cx + F(1, 2), cy - F(1, 2), cy + F(1, 2)
    m = ov(a, b, x0, x1) * ov(a, b, y0, y1)
    for p, q, e, w in H:
        if y0 <= e <= y1:
            m += w * ov(p, q, x0, x1) / (q - p)
    for p, q, e, w in V:
        if x0 <= e <= x1:
            m += w * ov(p, q, y0, y1) / (q - p)
    return m
vals = set([a, b])
for p, q, e, w in H + V:
    vals |= {p, q, e}
G = sorted(set(v + d for v in vals for d in (F(1, 2), F(-1, 2))))
lo, hi = F(1, 2), s - F(1, 2)
G = sorted(set([g for g in G if lo <= g <= hi] + [lo, hi]))
print('Gamma size', len(G), 'denominators', Counter(g.denominator for g in G))
gaps = min(G[i + 1] - G[i] for i in range(len(G) - 1))
dl = gaps / 10
best = None; cnt = 0; tight = 0
res = []
for gx in G:
    for gy in G:
        vc = mass(gx, gy)
        for sx in (-1, 1):
            for sy in (-1, 1):
                pts = []
                ok = True
                for k in (1, 2, 3):
                    cx, cy = gx + sx * k * dl, gy + sy * k * dl
                    if not (lo <= cx <= hi and lo <= cy <= hi):
                        ok = False; break
                    pts.append(mass(cx, cy))
                if not ok:
                    continue
                # quadratic in k through k=1,2,3 -> value at k=0
                f1, f2, f3 = pts
                lim = 3 * f1 - 3 * f2 + f3
                # sanity: multilinearity along diagonal => quadratic; check a 4th point
                f4 = mass(gx + sx * 4 * dl, gy + sy * 4 * dl) if (lo <= gx + sx*4*dl <= hi and lo <= gy + sy*4*dl <= hi) else None
                if f4 is not None:
                    assert f4 == f1 - 3 * f2 + 3 * f3 + 0 or f4 == 4*f3 - 6*f2 + 4*f1 - lim, (gx, gy)
                assert vc >= lim, ('usc fails', gx, gy, sx, sy, vc, lim)
                cnt += 1
                res.append((lim, gx, gy, sx, sy))
                if lim == 1: tight += 1
res.sort()
print('limits', cnt, 'min', res[0][0], 'tight', tight)
for r in res[:15]:
    print(float(r[0]), r[0], r[1:])
print('count < 1:', sum(1 for r in res if r[0] < 1))
