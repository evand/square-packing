import numpy as np
from fastf import mass, feasible
from scipy.optimize import minimize
from fractions import Fraction as Fr
from cover import mass_exact, rot, inside_box
def f(x):
    if not feasible(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]: return 10
    return mass(np.array([x[0]]),np.array([x[1]]),np.array([x[2]]))[0]
rng=np.random.default_rng(5); out=[]
for st in ([4.904328,3.407081,np.radians(44.5094)],[3.890290,3.406979,np.radians(46.0892)]):
  for k in range(30):
    x=np.array(st)+rng.normal(0,[0.02,0.01,0.02])
    r=minimize(f,x,method='Nelder-Mead',options=dict(xatol=1e-13,fatol=1e-16,maxiter=8000)); out.append((r.fun,*r.x))
out.sort()
for o in out[:4]:
    t=Fr(np.tan(o[3]/2)); c,s=rot(t); me=mass_exact(Fr(o[1]),Fr(o[2]),c,s)
    print('exact %.12f inside %s cx %.9f cy %.9f deg %.6f'%(float(me),inside_box(Fr(o[1]),Fr(o[2]),c,s),o[1],o[2],np.degrees(o[3])))
