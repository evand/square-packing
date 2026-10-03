import numpy as np
from cover import SEGS, PX0,PX1,PY0,PY1, PDENS, S
SA=np.array([[float(v) for v in sg] for sg in SEGS])
X0,Y0,X1,Y1,WW=SA.T
bx0,bx1,by0,by1=map(float,(PX0,PX1,PY0,PY1)); dens=float(PDENS)
BE=np.array([[bx0,by0,bx1,by0],[bx1,by0,bx1,by1],[bx1,by1,bx0,by1],[bx0,by1,bx0,by0]])
def clip_ab(p0a,p0b,p1a,p1b):
    lo=np.zeros(np.broadcast(p0a,p1a).shape);hi=np.ones_like(lo)
    for p0,p1 in ((p0a,p1a),(p0b,p1b)):
        d=p1-p0
        par=np.abs(d)<1e-300
        dd=np.where(par,1.0,d)
        l1=(-0.5-p0)/dd; l2=(0.5-p0)/dd
        inn=np.abs(p0)<=0.5
        lmin=np.where(par,np.where(inn,-np.inf,np.inf),np.minimum(l1,l2))
        lmax=np.where(par,np.where(inn,np.inf,-np.inf),np.maximum(l1,l2))
        lo=np.maximum(lo,lmin);hi=np.minimum(hi,lmax)
    return lo,hi
def seg_mass(cx,cy,c,s):
    cx=cx[:,None];cy=cy[:,None];c=c[:,None];s=s[:,None]
    a0=(X0-cx)*c+(Y0-cy)*s; b0=-(X0-cx)*s+(Y0-cy)*c
    a1=(X1-cx)*c+(Y1-cy)*s; b1=-(X1-cx)*s+(Y1-cy)*c
    lo,hi=clip_ab(a0,b0,a1,b1)
    return (WW*np.clip(hi-lo,0,None)).sum(1)
def leb_area(cx,cy,c,s):
    A=np.zeros(len(cx))
    # edges of Q (ccw) clipped to box
    cor=[(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]
    P=[(cx+a*c-b*s,cy+a*s+b*c) for a,b in cor]
    for k in range(4):
        (px,py),(rx,ry)=P[k],P[(k+1)%4]
        lo=np.zeros(len(cx));hi=np.ones(len(cx))
        for p0,p1,L,U in ((px,rx,bx0,bx1),(py,ry,by0,by1)):
            d=p1-p0; par=np.abs(d)<1e-300; dd=np.where(par,1.0,d)
            l1=(L-p0)/dd;l2=(U-p0)/dd
            inn=(p0>=L)&(p0<=U)
            lo=np.maximum(lo,np.where(par,np.where(inn,-np.inf,np.inf),np.minimum(l1,l2)))
            hi=np.minimum(hi,np.where(par,np.where(inn,np.inf,-np.inf),np.maximum(l1,l2)))
        ok=hi>lo
        A+=np.where(ok,(ry-py)*(px*(hi-lo)+(rx-px)*(hi**2-lo**2)/2),0)
    for (ex0,ey0,ex1,ey1) in BE:
        a0=(ex0-cx)*c+(ey0-cy)*s; b0=-(ex0-cx)*s+(ey0-cy)*c
        a1=(ex1-cx)*c+(ey1-cy)*s; b1=-(ex1-cx)*s+(ey1-cy)*c
        lo,hi=clip_ab(a0,b0,a1,b1)
        ok=hi>lo
        A+=np.where(ok,(ey1-ey0)*(ex0*(hi-lo)+(ex1-ex0)*(hi**2-lo**2)/2),0)
    return A
def mass(cx,cy,th):
    cx=np.asarray(cx,float);cy=np.asarray(cy,float);th=np.asarray(th,float)
    c=np.cos(th);s=np.sin(th)
    return seg_mass(cx,cy,c,s)+dens*leb_area(cx,cy,c,s)
def feasible(cx,cy,th):
    c=np.abs(np.cos(th));s=np.abs(np.sin(th)); h=(c+s)/2
    return (cx-h>=0)&(cx+h<=float(S))&(cy-h>=0)&(cy+h<=float(S))
