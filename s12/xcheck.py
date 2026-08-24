"""INDEPENDENT exact re-check of a certificate, written from scratch in Python Fractions.

For a sample of angle bins it computes, in exact rational arithmetic, the minimum over ALL
admissible centres of the weight covered by the closed sigma_k-square at angle theta_k.
Deliberately a separate implementation from verify/ (different language, different code path)
so that agreement is meaningful evidence.
"""
from fractions import Fraction as F
import sys, math

def load(path):
    t=open(path).read().split()
    sn,sd,D,WD,m=int(t[0]),int(t[1]),int(t[2]),int(t[3]),int(t[4])
    A=[(int(t[5+3*i]),int(t[6+3*i]),int(t[7+3*i])) for i in range(m)]
    return F(sn,sd),D,WD,A

def yrange_over_strip(poly,a,b):
    """exact [ymin,ymax] of the convex polygon over the vertical strip [a,b] (empty -> None)"""
    ys=[]
    n=len(poly)
    for i in range(n):
        p=poly[i]; q=poly[(i+1)%n]
        if a<=p[0]<=b: ys.append(p[1])
        for xc in (a,b):
            if (p[0]-xc)*(q[0]-xc)<=0 and p[0]!=q[0]:
                ys.append(p[1]+(q[1]-p[1])*(xc-p[0])/(q[0]-p[0]))
    if not ys: return None
    return min(ys),max(ys)

def check_angle(s,D,WD,A,k,N):
    c0=F(N*N-k*k,N*N+k*k); s0=F(2*k*N,N*N+k*k)
    k2=k+1; c1=F(N*N-k2*k2,N*N+k2*k2); s1=F(2*k2*N,N*N+k2*k2)
    cd=c0*c1+s0*s1; sd_=c0*s1-s0*c1
    sig=F(1,1)/(cd+sd_)                       # sigma_k = 1/(cos d + sin d)
    sig=F(int(sig*10**6),10**6)               # round DOWN (sound: smaller square)
    h=sig/2
    wm=c0+s0; wm=F(int(wm*10**6),10**6)       # round DOWN (sound: larger centre box)
    lo=wm/2; hi=s-wm/2
    if hi<=lo: return None
    Q=[(c0*F(x,D)+s0*F(y,D), -s0*F(x,D)+c0*F(y,D), w) for x,y,w in A]
    poly=[(c0*a+s0*b, -s0*a+c0*b) for a,b in [(lo,lo),(hi,lo),(hi,hi),(lo,hi)]]
    px0=min(p[0] for p in poly); px1=max(p[0] for p in poly)
    bx=sorted(set([q[0]-h for q in Q]+[q[0]+h for q in Q]+[px0,px1]))
    best=None
    for i in range(len(bx)-1):
        a,b=bx[i],bx[i+1]
        if b<=px0 or a>=px1: continue
        yr=yrange_over_strip(poly,a,b)
        if yr is None: continue
        ylo,yhi=yr
        act=[q for q in Q if q[0]-h<=a and b<=q[0]+h]      # CLOSED convention
        if not act: return 0
        by=sorted(set([q[1]-h for q in act]+[q[1]+h for q in act]+[ylo,yhi]))
        for j in range(len(by)-1):
            c_,d_=by[j],by[j+1]
            if d_<=ylo or c_>=yhi: continue
            tot=sum(q[2] for q in act if q[1]-h<=c_ and d_<=q[1]+h)
            if best is None or tot<best: best=tot
    return best

if __name__=="__main__":
    path=sys.argv[1]; N=int(sys.argv[2]); nsample=int(sys.argv[3])
    s,D,WD,A=load(path)
    tot=sum(a[2] for a in A)
    print(f"certificate {path}:  s={s}={float(s):.6f}  atoms={len(A)}  total weight={F(tot,WD)}={tot/WD:.7f}")
    sD=s*D; assert sD.denominator==1; sD=int(sD)
    S=sorted(A)
    for g in [lambda x,y:(sD-x,y), lambda x,y:(x,sD-y), lambda x,y:(y,x)]:
        assert sorted((*g(x,y),w) for x,y,w in A)==S, "atom set is NOT D4-symmetric"
    print("D4 symmetry of the atom set: OK  (so angles in [45,90) reduce to [0,45])")
    K=0
    while (K+N)**2 < 2*N*N: K+=1
    ks=sorted(set(int(i*K/(nsample-1)) for i in range(nsample)))
    worst=None
    for k in ks:
        v=check_angle(s,D,WD,A,k,N)
        if v is None: continue
        print(f"  bin k={k:4d}  theta={math.degrees(2*math.atan(k/N)):7.3f}deg   exact min covered weight = {F(v,WD)} = {v/WD:.7f}  {'OK' if v>=WD else '*** FAIL ***'}")
        if worst is None or v<worst: worst=v
    print(f"minimum over sampled bins = {worst/WD:.7f}   ->  {'all >= 1  (consistent with verify/)' if worst>=WD else 'VIOLATION'}")
