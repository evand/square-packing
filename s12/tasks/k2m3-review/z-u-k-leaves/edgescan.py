import sys
from ev import *
c = Cov('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
a = c.a
us = [F(sys.argv[1]), -F(sys.argv[1])]
res = []
for u in us:
    th = 2 * abs(u)
    for i in range(-20, 21):
        for j in range(50, 351):
            cx, cy = a + HALF + th * F(i, 20), F(j, 100)
            if not c.admissible(cx, cy, u): continue
            l, s = c.parts(cx, cy, u)
            if l == 1: continue
            res.append((l + s - 1, float(u), i, j))
res.sort()
for r in res[:8]: print(float(r[0]), r[1:])
print('below', sum(1 for r in res if r[0] < 0), 'n', len(res))
