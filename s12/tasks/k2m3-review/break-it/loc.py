import sys, numpy as np
from fractions import Fraction as F
from scipy.optimize import minimize
import ev
COV=ev.Cover('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
def f(x):
    cx,cy,th=x
    ee=(abs(np.cos(th))+abs(np.sin(th)))/2
    pen=max(0,ee-cx)+max(0,ee-cy)
    return COV.fmass([max(cx,ee)],[max(cy,ee)],[th])[0]+10*pen
def exact(x):
    cx,cy,th=x
    u=F(np.tan(th/2)).limit_denominator(10**12); c,s=ev.cs(u); E=(abs(c)+abs(s))/2
    X=max(F(cx).limit_denominator(10**13),E); Y=max(F(cy).limit_denominator(10**13),E)
    return COV.mass(X,Y,u),X,Y,u
starts=[eval(a) for a in sys.argv[1:]]
for st in starts:
    best=None
    for d in ([0,0,0],[0.01,0,0],[0,0.01,0],[0,0,-0.02],[0.02,-0.01,-0.05],[-0.02,0.01,0.03]):
        x0=np.array(st)+np.array(d); x0[2]=np.radians(st[2])+d[2]
        r=minimize(f,x0,method='Nelder-Mead',options=dict(xatol=1e-13,fatol=1e-16,maxiter=4000))
        m,X,Y,u=exact(r.x)
        if best is None or m<best[0]: best=(m,X,Y,u,r.x)
    m,X,Y,u,x=best
    print(st,'->',float(m-1), float(X),float(Y),np.degrees(x[2]), 'exact pose',X,Y,u)
