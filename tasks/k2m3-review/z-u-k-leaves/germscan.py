import sys, time
from ev import *
c = Cov('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
u = F(sys.argv[1]); th = 2 * u
lo, hi = F(1, 2), F(7, 2)
G = [F(2 * k + 5, 10) for k in range(-10, 40)]
G = [g for g in G if lo <= g <= hi]
K = int(sys.argv[2]) if len(sys.argv) > 2 else 14
res = []; t0 = time.time()
for gx in G:
    for gy in G:
        for i in range(-K, K + 1):
            for j in range(-K, K + 1):
                cx, cy = gx + th * F(i, 20), gy + th * F(j, 20)
                if not c.admissible(cx, cy, u): continue
                m = c.mass(cx, cy, u)
                res.append((m, gx, gy, i, j))
res.sort()
print('u', u, 'n', len(res), 'time', time.time() - t0)
for r in res[:12]:
    print(float(r[0] - 1), r[1:], )
print('below 1:', sum(1 for r in res if r[0] < 1))
