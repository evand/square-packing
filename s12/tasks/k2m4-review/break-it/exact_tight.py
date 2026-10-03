# exact evaluation at tiny rational tilts around tight germs on/near the doubled end lines x,y in {3,6} and corner/band junction
import json, sys, time
from fractions import Fraction as F
import cover
from cover import rot, inside_box, clip_area, corners, PX0,PX1,PY0,PY1,PDENS
ALL=cover.SEGS
def mass_local(cx,cy,c,s):
    h=F(1,2); m=F(0)
    for (x0,y0,x1,y1,w) in ALL:
        if abs(x0-cx)>1 or abs(y0-cy)>1: continue
        a0=(x0-cx)*c+(y0-cy)*s; b0=-(x0-cx)*s+(y0-cy)*c
        a1=(x1-cx)*c+(y1-cy)*s; b1=-(x1-cx)*s+(y1-cy)*c
        lo,hi=F(0),F(1)
        for (p0,p1) in ((a0,a1),(b0,b1)):
            d=p1-p0
            if d==0:
                if not (-h<=p0<=h): lo,hi=F(1),F(0)
            else:
                l1=(-h-p0)/d; l2=(h-p0)/d
                if l1>l2: l1,l2=l2,l1
                lo=max(lo,l1); hi=min(hi,l2)
        if hi>lo: m+=w*(hi-lo)
    return m+PDENS*clip_area(corners(cx,cy,c,s),PX0,PX1,PY0,PY1)
part,nparts=int(sys.argv[1]),int(sys.argv[2])
rows=[r for r in json.load(open('germs_near.json')) if r[0]=='1']
special={F(3),F(6),F(2),F(14,5),F(31,5),F(7),F(1),F(8)}
sel=[]
for r in rows:
    x0,y0,sg,X,Y=[F(v) for v in r[1:]]
    if {x0,x0+1,y0,y0+1}&{F(3),F(6)} or {x0,x0+1,y0,y0+1}&{F(14,5),F(31,5)}:
        sel.append((x0,y0,int(sg),X,Y))
sel=sel[part::nparts]
print('germs',len(sel)); t0=time.time(); minrel=None; nv=0; n=0
for (x0,y0,sg,X,Y) in sel:
  for tt in (F(1,10**6),F(1,10**4),F(1,300)):
    t=sg*tt; c,s=rot(t); th_approx=2*t   # theta ~ 2t
    for dX in (F(-1,10),F(-1,50),F(0),F(1,50),F(1,10)):
      for dY in (F(-1,10),F(-1,50),F(0),F(1,50),F(1,10)):
        cx0=x0+F(1,2); cy0=y0+F(1,2)
        dy=cx0-(X+dX); dx=(Y+dY)-cy0
        cx=cx0+th_approx*dx; cy=cy0+th_approx*dy
        if not inside_box(cx,cy,c,s): continue
        m=mass_local(cx,cy,c,s); n+=1
        rel=(m-1)/tt
        if minrel is None or rel<minrel[0]: minrel=(rel,m,(x0,y0,sg,X,Y),tt,dX,dY)
        if m<1: nv+=1; print('VIOLATION',m,cx,cy,t)
print('evaluated',n,'violations',nv,'time',time.time()-t0)
print('min (m-1)/t',float(minrel[0]),'m',float(minrel[1]),[str(v) for v in minrel[2]],minrel[3:])
