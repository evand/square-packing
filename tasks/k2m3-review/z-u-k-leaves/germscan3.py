# tiny-tilt exact scan restricted to squares partly overlapping U (0 < lam < 1), fine ratio grid
import sys
from ev import *
c = Cov('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
u = F(sys.argv[1]); th = 2 * u; K = int(sys.argv[2]); den = int(sys.argv[3])
lo, hi = F(1, 2), F(7, 2)
G = [F(2 * k + 5, 10) for k in range(-10, 40)]
G = [g for g in G if lo <= g <= hi and (F(11, 10) <= g <= F(27, 10))]
Gy = [F(2 * k + 5, 10) for k in range(-10, 40) if lo <= F(2 * k + 5, 10) <= hi]
res = []
for gx in G:
    for gy in Gy:
        for i in range(-K, K + 1):
            for j in range(-K, K + 1):
                cx, cy = gx + th * F(i, den), gy + th * F(j, den)
                if not c.admissible(cx, cy, u): continue
                l, s = c.parts(cx, cy, u)
                if l == 1 or l == 0: continue
                res.append((l + s - 1, gx, gy, i, j, l))
res.sort(key=lambda r: r[0])
print('u', u, 'n (0<lam<1, gx near a)', len(res))
for r in res[:20]:
    print('excess', float(r[0]), 'excess/theta', float(r[0] / th), 'gamma', r[1], r[2], 'ij', r[3], r[4], '/', den, 'lam', float(r[5]))
print('below 1:', sum(1 for r in res if r[0] < 0))
