import sys, numpy as np, json
from fractions import Fraction as F
from scipy.optimize import minimize
import ev
P='/home/evand/math/square-packing/public/s12/search/qx2_data/'
fn=sys.argv[1]; N=int(sys.argv[2]); seed=int(sys.argv[3]); thmax=float(sys.argv[4]) if len(sys.argv)>4 else np.pi/4
COV=ev.Cover(P+fn)
rng=np.random.default_rng(seed)
def e(th): return (abs(np.cos(th))+abs(np.sin(th)))/2
def f1(x):
    cx,cy,th=x; ee=e(th)
    pen=max(0,ee-cx)+max(0,ee-cy)+max(0,cx-3.5)+max(0,cy-3.5)+max(0,min(cx,cy)-(1.8+ee)+1e-9)+max(0,abs(th)-np.pi/4)
    return COV.fmass([min(max(cx,ee),3.5)],[min(max(cy,ee),3.5)],[th])[0]+10*pen
th=rng.uniform(-thmax,thmax,N); ee=np.array([e(t) for t in th])
cx=ee+rng.uniform(0,1,N)*(3.5-ee); cy=ee+rng.uniform(0,1,N)*(1.8-ee)
sw=rng.uniform(size=N)<0.5; cx,cy=np.where(sw,cy,cx),np.where(sw,cx,cy)
# bias: half the samples hug a wall (cy = e + tiny)
h=N//2; cy[:h]=np.minimum(cy[:h],3.5); cy[:h]=ee[:h]+rng.exponential(1e-3,h)*np.where(rng.uniform(size=h)<0.5,1,100)
vals=np.concatenate([COV.fmass(cx[i:i+2000],cy[i:i+2000],th[i:i+2000]) for i in range(0,N,2000)])
idx=np.argsort(vals)[:60]
res=[]
for i in idx:
    r=minimize(f1,[cx[i],cy[i],th[i]],method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-14,maxiter=3000))
    x=r.x; t=x[2]; u=F(np.tan(t/2)).limit_denominator(10**12)
    c,s=ev.cs(u); E=(abs(c)+abs(s))/2
    X=max(F(x[0]).limit_denominator(10**12),E); Y=max(F(x[1]).limit_denominator(10**12),E)
    X=min(X,7-E); Y=min(Y,7-E)
    m=COV.mass(X,Y,u)
    inU=all(F(9,5)<=a<=F(26,5) and F(9,5)<=b<=F(26,5) for a,b in ev.square(X,Y,u))
    res.append((float(m-1),inU,float(X),float(Y),float(np.degrees(t)),vals[i],str(X),str(Y),str(u)))
res.sort()
for r in res[:15]: print(r[:6])
json.dump(res,open(f'fs_{fn[:6]}_{seed}.json','w'))
