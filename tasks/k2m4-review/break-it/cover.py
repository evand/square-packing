# Independent parser + exact/float evaluators for the mixed box file (break-it review).
from fractions import Fraction as F
import numpy as np, sys
BOX='/home/evand/math/square-packing/public/s12/search/qx2_data/K4_k008_box9.txt'
def load(path=BOX):
    toks=[]
    for line in open(path):
        line=line.split('#')[0]
        toks+=line.split()
    assert toks[0]=='mixed' and toks[1]=='1'
    it=iter(toks[2:]); nx=lambda:int(next(it))
    sn,sd=nx(),nx(); D=nx(); W=nx()
    npt=nx(); pts=[(F(nx(),D),F(nx(),D),F(nx(),W)) for _ in range(npt)]
    ns=nx(); segs=[(F(nx(),D),F(nx(),D),F(nx(),D),F(nx(),D),F(nx(),W)) for _ in range(ns)]
    npg=nx(); polys=[]
    for _ in range(npg):
        k=nx(); w=F(nx(),W); vs=[(F(nx(),D),F(nx(),D)) for _ in range(k)]
        polys.append((vs,w))
    rest=list(it); assert not rest, rest
    return F(sn,sd),pts,segs,polys
S,PTS,SEGS,POLYS=load()
assert not PTS
assert len(POLYS)==1
PV,PW=POLYS[0]
PX0=min(v[0] for v in PV);PX1=max(v[0] for v in PV);PY0=min(v[1] for v in PV);PY1=max(v[1] for v in PV)
assert sorted(PV)==sorted([(PX0,PY0),(PX1,PY0),(PX1,PY1),(PX0,PY1)])
PDENS=PW/((PX1-PX0)*(PY1-PY0))
TOTAL=sum(s[4] for s in SEGS)+PW

def rot(t):
    t=F(t); d=1+t*t
    return (1-t*t)/d, 2*t/d   # cos, sin of theta=2 atan t

def clip_area(poly, x0,x1,y0,y1):
    def clip(P, inside, inter):
        out=[]
        n=len(P)
        for i in range(n):
            A=P[i];B=P[(i+1)%n]
            ia,ib=inside(A),inside(B)
            if ia: out.append(A)
            if ia!=ib: out.append(inter(A,B))
        return out
    P=poly
    for (axis,val,sgn) in ((0,x0,1),(0,x1,-1),(1,y0,1),(1,y1,-1)):
        if not P: return F(0)
        def inside(p,axis=axis,val=val,sgn=sgn): return sgn*(p[axis]-val)>=0
        def inter(A,B,axis=axis,val=val):
            l=(val-A[axis])/(B[axis]-A[axis])
            return (A[0]+l*(B[0]-A[0]),A[1]+l*(B[1]-A[1]))
        P=clip(P,inside,inter)
    if len(P)<3: return F(0)
    a=F(0)
    for i in range(len(P)):
        a+=P[i][0]*P[(i+1)%len(P)][1]-P[(i+1)%len(P)][0]*P[i][1]
    return abs(a)/2

def corners(cx,cy,c,s):
    h=F(1,2)
    return [(cx+a*c-b*s, cy+a*s+b*c) for (a,b) in ((-h,-h),(h,-h),(h,h),(-h,h))]

def inside_box(cx,cy,c,s):
    return all(0<=x<=S and 0<=y<=S for x,y in corners(cx,cy,c,s))

def mass_exact(cx,cy,c,s,detail=False):
    cx=F(cx);cy=F(cy);c=F(c);s=F(s)
    assert c*c+s*s==1
    h=F(1,2); m=F(0); parts=[]
    for (x0,y0,x1,y1,w) in SEGS:
        # a = (p-c).u, b=(p-c).v, u=(c,s), v=(-s,c)
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
        if hi>lo:
            m+=w*(hi-lo); 
            if detail: parts.append(((x0,y0,x1,y1,w),hi-lo))
    A=clip_area(corners(cx,cy,c,s),PX0,PX1,PY0,PY1)
    m+=PDENS*A
    return (m,parts,A) if detail else m

# float vectorised over many poses
SA=np.array([[float(v) for v in sg] for sg in SEGS])
X0,Y0,X1,Y1,WW=SA.T
def mass_float(cx,cy,th):
    """cx,cy,th arrays (n,) -> masses (n,). Lebesgue via polygon clip numerically."""
    cx=np.atleast_1d(cx)[:,None];cy=np.atleast_1d(cy)[:,None];th=np.atleast_1d(th)[:,None]
    c=np.cos(th);s=np.sin(th)
    a0=(X0-cx)*c+(Y0-cy)*s; b0=-(X0-cx)*s+(Y0-cy)*c
    a1=(X1-cx)*c+(Y1-cy)*s; b1=-(X1-cx)*s+(Y1-cy)*c
    lo=np.zeros_like(a0);hi=np.ones_like(a0)
    for p0,p1 in ((a0,a1),(b0,b1)):
        d=p1-p0
        with np.errstate(divide='ignore',invalid='ignore'):
            l1=(-0.5-p0)/d; l2=(0.5-p0)/d
        par=np.abs(d)<1e-15
        lmin=np.where(par,np.where(np.abs(p0)<=0.5+1e-12,-np.inf,np.inf),np.minimum(l1,l2))
        lmax=np.where(par,np.where(np.abs(p0)<=0.5+1e-12,np.inf,-np.inf),np.maximum(l1,l2))
        lo=np.maximum(lo,lmin);hi=np.minimum(hi,lmax)
    m=(WW*np.clip(hi-lo,0,None)).sum(1)
    m+=float(PDENS)*leb_area_float(cx[:,0],cy[:,0],th[:,0])
    return m
def leb_area_float(cx,cy,th):
    out=np.empty(len(cx))
    x0,x1,y0,y1=map(float,(PX0,PX1,PY0,PY1))
    for i in range(len(cx)):
        c=np.cos(th[i]);s=np.sin(th[i])
        P=[(cx[i]+a*c-b*s,cy[i]+a*s+b*c) for a,b in ((-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5))]
        # quick reject / full-inside
        xs=[p[0] for p in P];ys=[p[1] for p in P]
        if max(xs)<=x0 or min(xs)>=x1 or max(ys)<=y0 or min(ys)>=y1: out[i]=0;continue
        if min(xs)>=x0 and max(xs)<=x1 and min(ys)>=y0 and max(ys)<=y1: out[i]=1;continue
        for ax,val,sg in ((0,x0,1),(0,x1,-1),(1,y0,1),(1,y1,-1)):
            Q=[];n=len(P)
            for k in range(n):
                A=P[k];B=P[(k+1)%n]
                ia=sg*(A[ax]-val)>=0;ib=sg*(B[ax]-val)>=0
                if ia:Q.append(A)
                if ia!=ib:
                    l=(val-A[ax])/(B[ax]-A[ax]);Q.append((A[0]+l*(B[0]-A[0]),A[1]+l*(B[1]-A[1])))
            P=Q
            if not P:break
        if len(P)<3: out[i]=0;continue
        out[i]=abs(sum(P[k][0]*P[(k+1)%len(P)][1]-P[(k+1)%len(P)][0]*P[k][1] for k in range(len(P))))/2
    return out
