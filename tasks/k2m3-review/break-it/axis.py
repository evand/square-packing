import sys
from fractions import Fraction as F
import ev
COV=ev.Cover(sys.argv[1] if '/' in sys.argv[1] else '/home/evand/math/square-packing/public/s12/search/qx2_data/'+sys.argv[1])
eps=F(1,10**15); G=[F(2*i+5,10) for i in range(0,31)]
res=[]
for x in G:
  for y in G:
    for sx in (-1,1):
      for sy in (-1,1):
        cx,cy=x+sx*eps,y+sy*eps
        if not COV.admissible(cx,cy,0): continue
        res.append((COV.mass(cx,cy,0),str(x),str(y),sx,sy))
# also cell midpoints
for x in G[:-1]:
  for y in G[:-1]:
    res.append((COV.mass(x+F(1,10),y+F(1,10),0),'mid',str(x),str(y),0))
res.sort(key=lambda r:r[0])
print('n',len(res),'below1',sum(r[0]<1 for r in res),'eq1(limit≈1, <1+1e-12)',sum(r[0]<1+F(1,10**12) for r in res))
for r in res[:8]: print(float(r[0]-1), r[1:])
