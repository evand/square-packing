# Own exact germ scan: limits theta->0+/- of mass for squares converging to [x0,x0+1]x[y0,y0+1].
from cover import *
from collections import defaultdict
import itertools, json
G=F(1,5)
H=defaultdict(list); V=defaultdict(list)   # line -> list of (lo,hi,w)
for x0,y0,x1,y1,w in SEGS:
    if y0==y1: H[y0].append((min(x0,x1),max(x0,x1),w))
    else: V[x0].append((min(y0,y1),max(y0,y1),w))
def linemass(lst,lo,hi):
    m=F(0)
    for a,b,w in lst:
        l=max(a,lo);r=min(b,hi)
        if r>l: m+=w*(r-l)/(b-a)
    return m
def closed_mass(x0,y0):
    return mass_exact(x0+F(1,2),y0+F(1,2),1,0)
res=[]
grid=[G*i for i in range(0,41)]
gmin=None
for x0 in grid:
  for y0 in grid:
    M=closed_mass(x0,y0)
    xs=[x0+G*i for i in range(6)]; ys=[y0+G*i for i in range(6)]
    for sign in (1,-1):
      for X in xs:
        # sign +: bottom keeps x<=X, top keeps x>=X. drop the complement
        if sign==1: dropH=linemass(H[y0],X,x0+1)+linemass(H[y0+1],x0,X)
        else: dropH=linemass(H[y0],x0,X)+linemass(H[y0+1],X,x0+1)
        # wall feasibility
        if y0==0 and not (X==(x0 if sign==1 else x0+1)): pass  # still allowed? handled below
        for Y in ys:
          if sign==1: dropV=linemass(V[x0],y0,Y)+linemass(V[x0+1],Y,y0+1)
          else: dropV=linemass(V[x0],Y,y0+1)+linemass(V[x0+1],y0,Y)
          # wall constraints (sign + : left wall => Y=top, bottom wall => X=left, right wall => Y=bottom, top wall => X=right)
          ok=True
          if sign==1:
            if x0==0 and Y!=y0+1: ok=False
            if x0+1==S and Y!=y0: ok=False
            if y0==0 and X!=x0: ok=False
            if y0+1==S and X!=x0+1: ok=False
          else:
            if x0==0 and Y!=y0: ok=False
            if x0+1==S and Y!=y0+1: ok=False
            if y0==0 and X!=x0+1: ok=False
            if y0+1==S and X!=x0: ok=False
          if not ok: continue
          g=M-dropH-dropV
          res.append((g,x0,y0,sign,X,Y))
res.sort()
print('n',len(res),'min',float(res[0][0]),res[0][0])
tight=[r for r in res if r[0]<=1]
print('germs <=1:',len(tight),'  <1:',sum(1 for r in res if r[0]<1))
for r in res[:30]: print(float(r[0]),[str(v) for v in r[1:]])
json.dump([[str(v) for v in r] for r in res if r[0]<=F(1001,1000)],open('germs_near.json','w'))
# also theta=0 face min over grid
fm=min((closed_mass(x0,y0),x0,y0) for x0 in grid for y0 in grid)
print('theta0 grid min',float(fm[0]),fm)
