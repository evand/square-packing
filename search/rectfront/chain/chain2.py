# Corrected chain: signed vertical offsets, joint middle height y2.
# Pair (left a, right b): need dx >= rw(M_ab, dy), dy = y_b - y_a, M = A0 + B0 (centrally symmetric, not mirror symmetric).
import numpy as np, math, sys
from scipy.spatial import ConvexHull
def sqv(t):
    c,s=math.cos(t),math.sin(t)
    return np.array([[c*x-s*y,s*x+c*y] for x,y in [(.5,.5),(-.5,.5),(-.5,-.5),(.5,-.5)]])
def rw(P,Y):
    best=0.0; n=len(P)
    for i in range(n):
        a,b=P[i],P[(i+1)%n]
        if (a[1]-Y)*(b[1]-Y)<=0 and a[1]!=b[1]:
            best=max(best,a[0]+(Y-a[1])*(b[0]-a[0])/(b[1]-a[1]))
    return best
step=float(sys.argv[1]) if len(sys.argv)>1 else 3.0
deg=np.arange(0,90,step); T=np.radians(deg); K=(np.cos(T)+np.sin(T))/2; n=len(T)
P={(i,j):(lambda p:p[ConvexHull(p).vertices])(np.array([a+b for a in sqv(T[i]) for b in sqv(T[j])])) for i in range(n) for j in range(n)}
r2=math.sqrt(2)
for w in [float(x) for x in sys.argv[2].split(",")]:
    best=(9,None)
    for b in range(n):
        if 2*K[b]>w: continue
        for y2 in np.linspace(K[b],w-K[b],25):
            # left square a at y1, dy=y2-y1 ; right square c at y3, dy=y3-y2. min over endpoints of y-range
            Lb=min(K[a]+min(rw(P[a,b],y2-y1) for y1 in (K[a],w-K[a])) for a in range(n) if 2*K[a]<=w)
            Rb=min(K[c]+min(rw(P[b,c],y3-y2) for y3 in (K[c],w-K[c])) for c in range(n) if 2*K[c]<=w)
            if Lb+Rb<best[0]: best=(Lb+Rb,(deg[b],round(y2,3)))
    t=min(3,4+2*r2-2*w)
    print(f"w={w} chain h>={best[0]:.4f} need {t:.4f} slack {best[0]-t:+.4f} at middle tilt,y2 {best[1]}",flush=True)
