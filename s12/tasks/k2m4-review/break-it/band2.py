import numpy as np
from fastf import mass, feasible
from scipy.optimize import minimize
rng=np.random.default_rng(8); best=[]
for it in range(30):
    n=20000
    th=rng.uniform(0,np.pi/2,n); hh=(np.cos(th)+np.sin(th))/2
    p=rng.uniform(0,0.6,n); cx=2.8+hh-p; cy=rng.uniform(hh,4.5,n)
    m=mass(cx,cy,th); ok=feasible(cx,cy,th); m=np.where(ok,m,9)
    for i in np.argsort(m)[:5]: best.append((m[i],cx[i],cy[i],th[i],p[i]))
best.sort(); print(best[:5])
def f(x):
    if not feasible(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]: return 10
    hh=(abs(np.cos(x[2]))+abs(np.sin(x[2])))/2
    if x[0]-hh>=2.75 and x[1]-hh>=2.75: return 10
    return mass(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]
out=[]
for b in best[:40]:
    r=minimize(f,b[1:4],method='Nelder-Mead',options=dict(xatol=1e-11,fatol=1e-15,maxiter=4000)); out.append((r.fun,*r.x))
out.sort()
for o in out[:10]:
    hh=(abs(np.cos(o[3]))+abs(np.sin(o[3])))/2
    print('%.12f cx %.6f cy %.6f deg %.4f protr_x %.4f protr_y %.4f'%(o[0],o[1],o[2],np.degrees(o[3]),2.8-(o[1]-hh),2.8-(o[2]-hh)))
