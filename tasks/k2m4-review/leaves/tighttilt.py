"""At every exactly-tight theta = 0 one-sided limit corner of the full face (own Lemma Z grid, no D4), probe the
germ: centre offsets (sx e, sy e') into the cell and theta of both signs, exact masses.  Also report the minimum
first-order growth (mass - 1)/|u| at the smallest tilt."""
from fractions import Fraction as F
from mass import mass, M

g = sorted(set([F(k, 10) for k in range(5, 86) if k % 2 == 1]))
L = {F(1, 4): lambda s: (F(3, 4) - s) * 2, F(3, 4): lambda s: (s - F(1, 4)) * 2}
tight = []
for i in range(len(g) - 1):
    for j in range(len(g) - 1):
        x0, x1, y0, y1 = g[i], g[i + 1], g[j], g[j + 1]
        hx, hy = x1 - x0, y1 - y0
        pts = {(a, b): mass(x0 + a * hx, y0 + b * hy, F(0)) for a in L for b in L}
        for s in (0, 1):
            for t in (0, 1):
                v = sum(pts[(a, b)] * L[a](F(s)) * L[b](F(t)) for a in L for b in L)
                if v == 1: tight.append(((x0 + s * hx, y0 + t * hy), (1 - 2 * s, 1 - 2 * t)))
print('tight one-sided corners (full face):', len(tight))
mn = None; below = []; slopes = []
for (cx, cy), (sx, sy) in tight:
    for e in (F(0), F(1, 10 ** 9), F(1, 10 ** 5), F(1, 10 ** 3)):
        for e2 in (F(0), F(1, 10 ** 9), F(1, 10 ** 5), F(1, 10 ** 3)):
            x, y = cx + sx * e, cy + sy * e2
            for u in (F(1, 10 ** 12), F(-1, 10 ** 12), F(1, 10 ** 7), F(-1, 10 ** 7), F(1, 10 ** 4), F(-1, 10 ** 4), F(1, 200), F(-1, 200)):
                c, s = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
                hh = (c + abs(s)) / 2
                if not (hh <= x <= M - hh and hh <= y <= M - hh): continue
                v = mass(x, y, u)
                if v < 1: below.append((x, y, u, v))
                if mn is None or v < mn[0]: mn = (v, x, y, u)
                if e == 0 and e2 == 0 and abs(u) == F(1, 10 ** 12): slopes.append(((v - 1) / abs(u), cx, cy, sx, sy, u > 0))
print('poses below 1:', len(below), below[:5])
print('min mass', float(mn[0]), [str(t) for t in mn[1:]])
slopes.sort(key=lambda r: r[0])
print('smallest (mu - 1)/|u| at |u| = 1e-12 (centre exactly at the corner):', [(float(r[0]), str(r[1]), str(r[2]), r[5]) for r in slopes[:8]])
