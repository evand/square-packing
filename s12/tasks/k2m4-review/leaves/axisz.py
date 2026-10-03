"""Own Lemma Z on the whole theta = 0 face (no D4): centres [1/2, 17/2]^2.  On each open cell of the odd-tenths grid the
mass is multilinear; evaluate at 4 interior points per cell, extrapolate to the 4 corners (one-sided limits).
Also evaluate the actual (closed) value at every grid point and check it is >= every adjacent limit (usc)."""
from fractions import Fraction as F
from mass import mass, M
import itertools

g = sorted(set([F(1, 2), F(17, 2)] + [F(k, 10) for k in range(5, 86) if k % 2 == 1]))
print('grid', len(g), g[0], g[-1])
best = None; nlim = 0; tight = 0; usc_bad = 0
vals = {}
for i in range(len(g) - 1):
    for j in range(len(g) - 1):
        x0, x1, y0, y1 = g[i], g[i + 1], g[j], g[j + 1]
        hx, hy = x1 - x0, y1 - y0
        pts = {}
        for a in (F(1, 4), F(3, 4)):
            for b in (F(1, 4), F(3, 4)):
                pts[(a, b)] = mass(x0 + a * hx, y0 + b * hy, F(0))
        # bilinear through (1/4,3/4)^2: f(s,t) = sum f_ab * L_a(s) L_b(t), L_{1/4}(s) = (3/4 - s)*2, L_{3/4}(s) = (s - 1/4)*2
        L = {F(1, 4): lambda s: (F(3, 4) - s) * 2, F(3, 4): lambda s: (s - F(1, 4)) * 2}
        # check multilinearity at the cell centre
        fc = sum(pts[(a, b)] * L[a](F(1, 2)) * L[b](F(1, 2)) for a in L for b in L)
        assert fc == mass(x0 + hx / 2, y0 + hy / 2, F(0)), ('not bilinear', x0, y0)
        for s in (0, 1):
            for t in (0, 1):
                v = sum(pts[(a, b)] * L[a](F(s)) * L[b](F(t)) for a in L for b in L)
                nlim += 1
                if v == 1: tight += 1
                if best is None or v < best[0]: best = (v, (x0 + s * hx, y0 + t * hy), (i, j))
                key = (x0 + s * hx, y0 + t * hy)
                vals.setdefault(key, []).append(v)
for key, lims in vals.items():
    v = mass(key[0], key[1], F(0))
    if v < max(lims): usc_bad += 1
print('limits', nlim, 'min', best[0], float(best[0]), 'at', [str(t) for t in best[1]], 'tight', tight,
      'below 1:', sum(1 for l in vals.values() for v in l if v < 1), 'usc violations', usc_bad)
