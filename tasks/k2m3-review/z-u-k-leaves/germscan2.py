import sys, time
from ev import *
c = Cov('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
u = F(sys.argv[1]); th = 2 * u
lo, hi = F(1, 2), F(7, 2)
G = [F(2 * k + 5, 10) for k in range(-10, 40)]
G = [g for g in G if lo <= g <= hi]
xr = (F(sys.argv[3]), F(sys.argv[4])) if len(sys.argv) > 4 else (lo, hi)
K = int(sys.argv[2])
res = []
for gx in G:
    if not (xr[0] <= gx <= xr[1]): continue
    for gy in G:
        for i in range(-K, K + 1):
            for j in range(-K, K + 1):
                cx, cy = gx + th * F(i, 20), gy + th * F(j, 20)
                if not c.admissible(cx, cy, u): continue
                l, s = c.parts(cx, cy, u)
                if l == 1: continue
                res.append(((l + s - 1) / th, l + s - 1, gx, gy, i, j))
res.sort(key=lambda r: r[1])
print('u', u, 'n (not inside U)', len(res))
for r in res[:25]:
    print('excess', float(r[1]), 'excess/theta', float(r[0]), 'gamma', r[2], r[3], 'ij/20', r[4], r[5])
print('below 1:', sum(1 for r in res if r[1] < 0))
