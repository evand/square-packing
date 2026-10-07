import numpy as np
from fastf import mass, feasible
from scipy.optimize import minimize
from fractions import Fraction as Fr
from cover import mass_exact, rot, inside_box
rng=np.random.default_rng(31)
def f(x):
    if not feasible(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]: return 10
    return mass(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]
# dense sample: square protruding below y=2.8, cx in [3.2,4.5], angle 30..60
n=400000; th=rng.uniform(np.radians(25),np.radians(65),n); hh=(np.cos(th)+np.sin(th))/2
p=rng.uniform(0.0,0.35,n); cy=2.8+hh-p; cx=rng.uniform(2.5,4.5,n)
m=np.concatenate([mass(cx[i:i+20000],cy[i:i+20000],th[i:i+20000]) for i in range(0,n,20000)])
idx=np.argsort(m)
# skip flat ones (p tiny): choose lowest with p>0.01
cand=[i for i in idx if p[i]>0.01][:40]
out=[]
for i in cand:
    r=minimize(f,[cx[i],cy[i],th[i]],method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-15,maxiter=5000))
    hhh=(abs(np.cos(r.x[2]))+abs(np.sin(r.x[2])))/2
    if 2.8-(r.x[1]-hhh)>0.005: out.append((r.fun,*r.x))
out.sort()
for o in out[:6]:
    t=Fr(np.tan(o[3]/2)); c,s=rot(t); me=mass_exact(Fr(o[1]),Fr(o[2]),c,s)
    hhh=(abs(np.cos(o[3]))+abs(np.sin(o[3])))/2
    print('exact %.12f inside %s cx %.6f cy %.6f deg %.4f protr %.4f'%(float(me),inside_box(Fr(o[1]),Fr(o[2]),c,s),o[1],o[2],np.degrees(o[3]),2.8-(o[2]-hhh)))
