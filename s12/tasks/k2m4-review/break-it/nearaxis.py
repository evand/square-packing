# near-axis multi-scale search: centre = grid-aligned square centre + independent multi-scale offsets, tilt multi-scale
import numpy as np, sys, time
from fastf import mass, feasible
from fractions import Fraction as Fr
from cover import mass_exact, rot, inside_box
seed=int(sys.argv[1]); N=int(sys.argv[2])
rng=np.random.default_rng(seed); nv=0; worst=[]; t0=time.time(); tot=0
for it in range(N//20000):
    n=20000
    i=rng.integers(0,23,n); j=rng.integers(0,23,n)   # x0=i/5 in [0,4.4]
    sgn=lambda: rng.choice([-1,1],n)
    def off():
        r=rng.random(n)
        return np.where(r<0.2,0,np.where(r<0.6,sgn()*10**rng.uniform(-8,-1,n),rng.uniform(-0.1,0.1,n)))
    th=sgn()*10**rng.uniform(-8,-1,n)
    # option: offsets proportional to theta
    prop=rng.random(n)<0.4
    ox=np.where(prop,th*rng.uniform(-0.7,0.7,n),off()); oy=np.where(prop,th*rng.uniform(-0.7,0.7,n),off())
    cx=i/5+0.5+ox; cy=j/5+0.5+oy
    ok=feasible(cx,cy,th)
    m=mass(cx[ok],cy[ok],th[ok]); tot+=ok.sum()
    v=m<1-1e-12; nv+=int(v.sum())
    for k in np.where(v)[0][:5]: worst.append((m[k],cx[ok][k],cy[ok][k],th[ok][k]))
    # also keep lowest nonflat for info
print('N',tot,'float-viol',nv,'time',time.time()-t0)
for w in sorted(worst)[:20]:
    m,cx,cy,th=w; t=Fr(np.tan(th/2)); c,s=rot(t)
    me=mass_exact(Fr(cx),Fr(cy),c,s); ins=inside_box(Fr(cx),Fr(cy),c,s)
    print('float %.15f exact %.15f inside %s pose %r %r %r'%(m,float(me),ins,cx,cy,th))
