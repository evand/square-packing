import numpy as np
from fractions import Fraction as F
from scipy.optimize import minimize
import ev
COV=ev.Cover('/home/evand/math/square-packing/public/s12/search/qx2_data/L4_k02_box7.txt')
for cx0,cy0 in ((2.9,2.4071),(3.1,2.4071)):
  for deg in [36,38,40,41,42,43,44,44.5,44.9,44.99,45]:
    th=np.radians(deg)
    f=lambda x: COV.fmass([x[0]],[x[1]],[th])[0]+10*max(0,x[1]-2.45)
    r=minimize(f,[cx0,cy0],method='Nelder-Mead',options=dict(xatol=1e-13,fatol=1e-16))
    u=F(np.tan(th/2)).limit_denominator(10**9)
    m=COV.mass(F(r.x[0]).limit_denominator(10**12),F(r.x[1]).limit_denominator(10**12),u)
    print(cx0,deg, r.x, float(m-1))
