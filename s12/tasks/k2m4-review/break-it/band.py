# poses protruding slightly out of the Lebesgue box: check mass-1 vs protrusion
import numpy as np
from fastf import mass, feasible
rng=np.random.default_rng(7); worst=[];nv=0;N=0
for it in range(40):
    n=20000
    th=np.where(rng.random(n)<0.3,10**rng.uniform(-6,-1,n),rng.uniform(0,np.pi/2,n))
    hh=(np.cos(th)+np.sin(th))/2
    p=10**rng.uniform(-7,-1,n)  # protrusion on left side
    cx=2.8+hh-p
    cy=rng.uniform(hh,4.5,n)
    m=mass(cx,cy,th); N+=n
    nv+=int((m<1-1e-12).sum())
    r=(m-1)/p; i=np.argmin(r); worst.append((r[i],m[i],cx[i],cy[i],th[i],p[i]))
worst.sort(); print('N',N,'viol',nv)
for w in worst[:8]: print('(m-1)/p %.5f m %.15f cx %.9f cy %.9f th %.3e p %.2e'%w)
