import sys, multiprocessing as mp
from fractions import Fraction as F
import ev, random
COV=ev.Cover('/home/evand/math/square-packing/public/s12/search/qx2_data/'+sys.argv[1])
random.seed(0)
CX=[F(50+i,100) for i in range(0,301,1)]+[F(random.randint(500000,3500000),10**6) for _ in range(100)]
def work(cx):
    best=(9,None)
    for base in (F(1,2),F(3,2),F(1,1),F(2,1)):
        for k in (3,6,9,12):
            for sg in (1,-1):
                u=sg*F(1,10**k); c,s=ev.cs(u); e=(abs(c)+abs(s))/2
                for ty in [F(i,4) for i in range(-8,13)]:
                    cy=base+abs(u)*ty
                    if base==F(1,2): cy=max(cy,e)
                    if not COV.admissible(cx,cy,u): continue
                    m=COV.mass(cx,cy,u)
                    if m<best[0]: best=(m,(str(cx),str(cy),str(u)))
    return best
with mp.Pool(2) as p: res=p.map(work,CX)
res.sort(key=lambda r:r[0])
print('below1',sum(r[0]<1 for r in res),'n',len(res))
for r in res[:10]: print(float(r[0]-1),r[1])
