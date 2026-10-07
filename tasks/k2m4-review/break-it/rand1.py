# global random sampling, D4-reduced (quadrant), all angles; report lowest masses
import numpy as np, sys, time
from fastf import mass, feasible
rng=np.random.default_rng(int(sys.argv[1]))
N=int(sys.argv[2]); best=[]
t0=time.time()
for it in range(N//20000):
    n=20000
    th=rng.uniform(0,np.pi/2,n)
    # bias half to tiny angles
    k=n//2; th[:k]=10**rng.uniform(-7,-1,k)*rng.choice([1,-1],k)%(np.pi/2)
    h=(np.abs(np.cos(th))+np.abs(np.sin(th)))/2
    cx=rng.uniform(h,4.5);cy=rng.uniform(h,4.5)
    m=mass(cx,cy,th)
    idx=np.argsort(m)[:5]
    best+= [(m[i],cx[i],cy[i],th[i]) for i in idx]
best.sort(); print('time',time.time()-t0)
for b in best[:15]: print('%.12f %.9f %.9f %.3e'%b)
