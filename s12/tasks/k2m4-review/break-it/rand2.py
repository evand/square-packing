# global random sampling + Nelder-Mead refinement (quadrant-reduced), report violations and lowest non-flat poses
import numpy as np, sys, time
from fastf import mass, feasible
from scipy.optimize import minimize
rng=np.random.default_rng(int(sys.argv[1]))
N=int(sys.argv[2]); mode=sys.argv[3]
def inside_leb(cx,cy,th):
    h=(np.abs(np.cos(th))+np.abs(np.sin(th)))/2
    return (cx-h>=2.8)&(cx+h<=6.2)&(cy-h>=2.8)&(cy+h<=6.2)
cands=[];viol=0;t0=time.time()
for it in range(N//20000):
    n=20000
    if mode=='all':
        th=rng.uniform(0,np.pi/2,n)
    else:
        th=10**rng.uniform(-6,-0.5,n)*rng.choice([1,-1],n)%(np.pi/2)
    h=(np.abs(np.cos(th))+np.abs(np.sin(th)))/2
    cx=rng.uniform(h,4.5);cy=rng.uniform(h,4.5)
    m=mass(cx,cy,th)
    viol+=int((m<1-1e-12).sum())
    m=np.where(inside_leb(cx,cy,th),9,m)
    idx=np.argsort(m)[:3]
    cands+=[(m[i],cx[i],cy[i],th[i]) for i in idx]
cands.sort()
print('sampled',N,'violations',viol,'time',time.time()-t0)
def f(p):
    cx,cy,th=p
    if not feasible(np.array([cx]),np.array([cy]),np.array([th]))[0]: return 10
    return mass(np.array([cx]),np.array([cy]),np.array([th]))[0]
out=[]
for m,cx,cy,th in cands[:60]:
    r=minimize(f,[cx,cy,th],method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-15,maxiter=3000,initial_simplex=None))
    out.append((r.fun,*r.x,m))
out.sort()
for o in out[:12]: print('refined %.13f  at %.9f %.9f th=%.6e (deg %.6f) from %.6f'%(o[0],o[1],o[2],o[3],np.degrees(o[3]),o[4]))
from fractions import Fraction as Fr
from cover import mass_exact, rot, inside_box
print('exact checks of refined poses with float mass < 1+1e-9:')
for o in out:
    if o[0]<1+1e-9:
        cx,cy,th=o[1],o[2],o[3]
        t=Fr(np.tan(th/2)); c,s=rot(t)
        ins=inside_box(Fr(cx),Fr(cy),c,s)
        me=mass_exact(Fr(cx),Fr(cy),c,s)
        print('  float %.15f exact %.15f inside %s  pose %r %r t=%r'%(o[0],float(me),ins,cx,cy,float(t)), 'VIOLATION' if ins and me<1 else '')
