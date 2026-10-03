import sys, numpy as np, json
from fractions import Fraction as F
from scipy.optimize import minimize
import ev
P='/home/evand/math/square-packing/public/s12/search/qx2_data/'
fn=sys.argv[1]; ths=[float(x) for x in sys.argv[2].split(',')]; N=int(sys.argv[3]); seed=int(sys.argv[4])
DL=float(sys.argv[5]); COV=ev.Cover(P+fn); rng=np.random.default_rng(seed)
out=[]
for deg in ths:
    th=np.radians(deg); ee=(abs(np.cos(th))+abs(np.sin(th)))/2
    def f(x):
        cx,cy=x
        pen=max(0,ee-cx)+max(0,ee-cy)+max(0,cx-3.5)+max(0,cy-3.5)+max(0,min(cx,cy)-(1.8+ee)+DL)
        cxc=min(max(cx,ee),3.5); cyc=min(max(cy,ee),3.5)
        return COV.fmass([cxc],[cyc],[th])[0]+10*pen
    a=ee+rng.uniform(0,1,N)*(3.5-ee); b=ee+rng.uniform(0,1,N)*(1.8-DL)
    h=N//3; b[:h]=ee+rng.exponential(1e-2,h)
    sw=rng.uniform(size=N)<0.5; cx=np.where(sw,b,a); cy=np.where(sw,a,b)
    v=np.concatenate([COV.fmass(cx[i:i+2000],cy[i:i+2000],np.full(min(2000,N-i),th)) for i in range(0,N,2000)])
    best=[]
    for i in np.argsort(v)[:25]:
        r=minimize(f,[cx[i],cy[i]],method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-15,maxiter=2000))
        u=F(np.tan(th/2)).limit_denominator(10**9); c,s=ev.cs(u); E=(abs(c)+abs(s))/2
        X=min(max(F(r.x[0]).limit_denominator(10**12),E),7-E); Y=min(max(F(r.x[1]).limit_denominator(10**12),E),7-E)
        m=COV.mass(X,Y,u); best.append((float(m-1),float(X),float(Y),str(X),str(Y),str(u)))
    best.sort(); out.append((deg,best[:5]))
    print(deg, best[0][:3], best[1][:3], flush=True)
json.dump(out,open(f'fixth_{seed}_{DL}.json','w'))
