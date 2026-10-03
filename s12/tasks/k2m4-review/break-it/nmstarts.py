# many NM local minimisations from random starts outside the Lebesgue core; report lowest local minima (excluding theta~0 germs)
import numpy as np, sys
from fastf import mass, feasible
from scipy.optimize import minimize
from fractions import Fraction as Fr
from cover import mass_exact, rot, inside_box
rng=np.random.default_rng(int(sys.argv[1])); K=int(sys.argv[2])
def f(x):
    if not feasible(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]: return 10
    hh=(abs(np.cos(x[2]))+abs(np.sin(x[2])))/2
    if x[0]-hh>=2.7 and x[1]-hh>=2.7: return 10
    return mass(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]
out=[]
for k in range(K):
    th=rng.uniform(0.002,np.pi/2-0.002); hh=(np.cos(th)+np.sin(th))/2
    x=[rng.uniform(hh,4.5),rng.uniform(hh,4.5),th]
    if f(x)>=10: continue
    r=minimize(f,x,method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-14,maxiter=3000))
    out.append((r.fun,*r.x))
out.sort(key=lambda o:o[0])
print('starts',len(out))
gen=[o for o in out if min(abs(o[3]%(np.pi/2)),np.pi/2-abs(o[3]%(np.pi/2)))>1e-3]
print('lowest generic-angle (|theta mod 90|>0.057deg) minima:')
for o in gen[:8]:
    t=Fr(np.tan(o[3]/2)); c,s=rot(t); me=mass_exact(Fr(o[1]),Fr(o[2]),c,s)
    print(' float %.10f exact %.10f inside %s cx %.6f cy %.6f deg %.4f'%(o[0],float(me),inside_box(Fr(o[1]),Fr(o[2]),c,s),o[1],o[2],np.degrees(o[3])))
print('overall lowest:',['%.10f'%o[0] for o in out[:5]])
