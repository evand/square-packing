import random
from fractions import Fraction as F
import ev
COV=ev.Cover('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
random.seed(11)
def rnd(sc): return F(random.randint(-10**6,10**6),10**6)*sc
# base poses: tight germs (corner, t, sign) and the 45-degree diamond
bases=[]
for (cx,cy,sg,tx,ty) in [('3/2','5/2',-1,'1','3/5'),('3/2','29/10',1,'1','3/5'),('3/2','3/2',1,'1','1/5'),('3/2','23/10',1,'1','3/5'),
                         ('1/2','29/10',-1,'1','-13/20'),('1/2','5/2',1,'1','13/20'),('1/2','27/10',1,'1','-3/20'),('3/2','5/2',1,'1','19/20')]:
    bases.append(('germ',F(cx),F(cy),sg,F(tx),F(ty)))
low=(9,None); cnt=0; below=0
for b in bases:
    _,cx,cy,sg,tx,ty=b
    for k in range(4,15):
        u0=F(1,10**k)
        for _ in range(60):
            u=sg*u0*(1+rnd(F(1,2)))
            x=cx+abs(u)*tx+rnd(F(1,10**k)); y=cy+abs(u)*ty+rnd(F(1,10**k))
            if not COV.admissible(x,y,u): continue
            m=COV.mass(x,y,u); cnt+=1; below+=m<1
            r=(m-1)/abs(u)
            if r<low[0]: low=(r,(str(x),str(y),str(u)),float(m-1))
print('germ perturbations',cnt,'below1',below,'min (m-1)/|u|',float(low[0]),low[1:])
# 45 degree diamond near (2.9, 1.7+sqrt2/2)
import math
low=(9,None); cnt=0
for k in range(3,15):
    for _ in range(150):
        th=math.pi/4+float(rnd(F(1,10**k)))
        u=F(math.tan(th/2)).limit_denominator(10**15); c,s=ev.cs(u); e=(abs(c)+abs(s))/2
        x=F(29,10)+rnd(F(1,10**k)); y=F(17,10)+e+rnd(F(1,10**k))
        m=COV.mass(x,y,u); cnt+=1
        if m<low[0]: low=(m,(str(x),str(y),str(u)))
print('45deg perturbations',cnt,'min m-1',float(low[0]-1),low[1])
