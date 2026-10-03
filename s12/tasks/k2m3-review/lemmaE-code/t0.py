import time, random
from fractions import Fraction as F
from mass import *
import qx2_zm as QZ
cov,a,b=load(); M=Mass(cov,a,b)
ex=QZ.Exact(cov,a,b)
# sanity: wall pose
for u in (F(1,10**6),F(1,1000),F(1,50)):
    print(u, float(M.exact(F(3),F(1,2)+F(1,10**6),u)), M.flt(3.0,0.5+1e-6,float(u)))
random.seed(1); mn=9
for i in range(3000):
    cx=F(random.randint(700,3500),1000); cy=F(random.randint(700,3500),1000); u=F(random.randint(1,400),1000)
    c,s=QZ.zm.trig(u); w=c+s
    if cx<w/2 or cy<w/2: continue
    m=M.exact(cx,cy,u); mf=M.flt(float(cx),float(cy),float(u))
    assert abs(float(m)-mf)<1e-9,(cx,cy,u,m,mf)
    mn=min(mn,m)
print('min',float(mn))
for box in [(F(3),F(3)+F(1,80),F(1,2),F(1,2)+F(1,80),F(0),F(1,500)),
            (F(23,10),F(23,10)+F(1,160),F(23,10),F(23,10)+F(1,160),F(0),F(1,1000)),
            (F(3,2),F(3,2)+F(1,80),F(3,2),F(3,2)+F(1,80),F(0),F(1,500))]:
    t=time.time(); r=ex.certify(box); print(r, time.time()-t, ex.nbrk)
    t=time.time(); r=ex.certify(box,tau=F(1)+F(1,10**9)); print('tau+',r, time.time()-t)
