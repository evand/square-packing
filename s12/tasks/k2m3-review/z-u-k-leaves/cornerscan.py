import sys
from ev import *
c = Cov('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
us = [F(1, 10**6), F(1, 10**4), F(1, 1000), F(1, 200), F(1, 100), F(3, 100), F(1, 20), F(1, 10), F(1, 5), F(2, 5)]
k0, k1 = int(sys.argv[1]), int(sys.argv[2])
res = []
for u in us:
    for i in range(k0, k1):
        for j in range(-40, 41):
            cx, cy = F(23, 10) + F(i, 200), F(23, 10) + F(j, 200)
            l, s = c.parts(cx, cy, u)
            res.append((l + s - 1, float(u), float(cx), float(cy), float(l)))
res.sort()
for r in res[:10]: print(r[0].__float__(), r[1:])
print('below:', sum(1 for r in res if r[0] < 0), 'n', len(res))
