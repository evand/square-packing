import numpy as np
from fractions import Fraction as F
from scipy.optimize import minimize
from fastf import mass, feasible
from cover import mass_exact, rot, inside_box
def f(p):
    cx,cy,th=p
    if not feasible(np.array([cx]),np.array([cy]),np.array([th]))[0]: return 10
    return mass(np.array([cx]),np.array([cy]),np.array([th]))[0]
for st in ([3.98,1.5,1e-5],[3.98,1.5001,1e-4],[3.9795,1.5,1e-6]):
    r=minimize(f,st,method='Nelder-Mead',options=dict(xatol=1e-14,fatol=1e-16,maxiter=5000))
    cx,cy,th=r.x
    print(r.fun, repr(cx),repr(cy),repr(th))
    t=F(np.tan(th/2)); c,s=rot(t)
    print('  exact', float(mass_exact(F(cx),F(cy),c,s)), 'inside',inside_box(F(cx),F(cy),c,s),'dy/th',(cy-1.5)/th)
