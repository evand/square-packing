#!/usr/bin/env python3
"""Independent theta = 0 check on the box file, full container (no D4), exact one-sided limits at all cell corners."""
from fractions import Fraction as F
B = '/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt'
tok = []
for line in open(B): tok += line.split('#', 1)[0].split()
it = iter(tok); next(it); next(it)
sn, sd, Dd, W = (int(next(it)) for _ in range(4)); k = F(sn, sd)
assert int(next(it)) == 0
segs = [tuple(int(next(it)) for _ in range(5)) for _ in range(int(next(it)))]
segs = [(F(a, Dd), F(b, Dd), F(c, Dd), F(d, Dd), F(w, W)) for a, b, c, d, w in segs]
a = b = None
npg = int(next(it)); kk = int(next(it)); wm = int(next(it)); V = [(F(int(next(it)), Dd), F(int(next(it)), Dd)) for _ in range(kk)]
a, b = min(v[0] for v in V), max(v[0] for v in V); assert F(wm, W) == (b - a) ** 2
H = [(y0, min(x0, x1), max(x0, x1), w) for x0, y0, x1, y1, w in segs if y0 == y1]
Vv = [(x0, min(y0, y1), max(y0, y1), w) for x0, y0, x1, y1, w in segs if x0 == x1]
assert len(H) + len(Vv) == len(segs)
h = F(1, 2)
pts = set([a, b])
for p, t0, t1, w in H + Vv: pts.update([p, t0, t1])
grid = sorted(set(v + s for v in pts for s in (-h, h) if h <= v + s <= k - h) | {h, k - h})
def inside(p, c, s):
    lo, hi = c - h, c + h
    return (lo < p <= hi) if s > 0 else (lo <= p < hi)
def ov(a0, a1, b0, b1): return max(F(0), min(a1, b1) - max(a0, b0))
best = None; n = 0; tight = 0
for cx in grid:
    for cy in grid:
        for sx in (-1, 1):
            if (cx == h and sx < 0) or (cx == k - h and sx > 0): continue
            for sy in (-1, 1):
                if (cy == h and sy < 0) or (cy == k - h and sy > 0): continue
                v = ov(a, b, cx - h, cx + h) * ov(a, b, cy - h, cy + h)
                for p, t0, t1, w in H:
                    if inside(p, cy, sy): v += w * ov(t0, t1, cx - h, cx + h) / (t1 - t0)
                for p, t0, t1, w in Vv:
                    if inside(p, cx, sx): v += w * ov(t0, t1, cy - h, cy + h) / (t1 - t0)
                n += 1; tight += v == 1
                if best is None or v < best[0]: best = (v, cx, cy, sx, sy)
print('grid', len(grid), 'corners', n, 'min', best[0], [str(t) for t in best[1:]], 'tight', tight)
