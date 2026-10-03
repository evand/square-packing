import numpy as np
from fractions import Fraction as F
from cover import mass_exact, rot
from fastf import mass
rng=np.random.default_rng(1)
n=300
cx=rng.uniform(.8,8.2,n);cy=rng.uniform(.8,8.2,n);t=rng.uniform(-.5,.5,n)
mx=0
for i in range(n):
    tt=F(t[i]).limit_denominator(10**6); c,s=rot(tt)
    ex=mass_exact(F(cx[i]),F(cy[i]),c,s)
    th=2*np.arctan(float(tt))
    fl=mass(np.array([cx[i]]),np.array([cy[i]]),np.array([th]))[0]
    mx=max(mx,abs(fl-float(ex)))
print('max err',mx)
