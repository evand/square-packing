# stratified grid over quadrant x angle, then NM refinement from lowest non-flat local candidates; exact confirm
import numpy as np, sys, time
from fastf import mass, feasible
from scipy.optimize import minimize
from fractions import Fraction as Fr
from cover import mass_exact, rot, inside_box
part,nparts=int(sys.argv[1]),int(sys.argv[2])
step=0.05
ths=np.radians(np.arange(0.25,90,0.5))[part::nparts]
res=[];t0=time.time();nv=0
for th0 in ths:
    h=(np.cos(th0)+np.sin(th0))/2
    xs=np.arange(h,4.5+1e-9,step)
    CX,CY=np.meshgrid(xs,xs);CX=CX.ravel();CY=CY.ravel();TH=np.full(CX.size,th0)
    m=mass(CX,CY,TH)
    nv+=int((m<1-1e-12).sum())
    hh=(np.cos(th0)+np.sin(th0))/2; nonflat=~((CX-hh>=2.7)&(CY-hh>=2.7)&(CX+hh<=6.3)&(CY+hh<=6.3))
    idx=np.where(nonflat)[0]
    j=idx[np.argsort(m[idx])[:4]]
    res+=[(m[i],CX[i],CY[i],TH[i]) for i in j]
print('grid done',time.time()-t0,'viol',nv)
res.sort()
def f(p):
    cx,cy,th=p
    if not feasible(np.array([cx]),np.array([cy]),np.array([th]))[0]: return 10
    hh=(abs(np.cos(th))+abs(np.sin(th)))/2
    if cx-hh>=2.7 and cy-hh>=2.7 and cx+hh<=6.3 and cy+hh<=6.3: return 10
    return mass(np.array([cx]),np.array([cy]),np.array([th]))[0]
out=[]
for m,cx,cy,th in res[:120]:
    r=minimize(f,[cx,cy,th],method='Nelder-Mead',options=dict(xatol=1e-11,fatol=1e-15,maxiter=4000))
    out.append((r.fun,*r.x,m))
out.sort()
for o in out[:10]:
    cx,cy,th=o[1],o[2],o[3]
    t=Fr(np.tan(th/2)); c,s=rot(t); me=mass_exact(Fr(cx),Fr(cy),c,s); ins=inside_box(Fr(cx),Fr(cy),c,s)
    print('refined float %.13f exact %.13f inside %s at %.9f %.9f deg=%.6f from %.6f'%(o[0],float(me),ins,cx,cy,np.degrees(th),o[4]))
print('time',time.time()-t0)
