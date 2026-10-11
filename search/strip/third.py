import numpy as np
from scipy.optimize import minimize
from rectpack import penalty, gap, half
v = np.load('leaf_w1.99_s2.npy'); n = 12
x, y, t = v[:4], v[n:n+4], v[2*n:2*n+4]
print("strip at L=4:"); [print('  x=%.3f y=%.3f th=%.1f'%(x[k],y[k],np.degrees(t[k])%90)) for k in range(4)]
rng = np.random.default_rng(0)
for W in (1.99, 1.97, 1.95):
  for H in (3.99, 3.95):
    res = None
    for s in range(60):
        z = np.r_[x + rng.normal(0,0.03,4), y*H/4 + rng.normal(0,0.03,4), t + rng.normal(0,0.03,4)]
        z = minimize(lambda q: penalty(q, 4, W, H), z, method='L-BFGS-B', options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
        if gap(z, 4, W, H) > -1e-10: res = z; break
    if res is None: print(f"W={W} H={H}: not found"); continue
    a = np.degrees(res[8:]) % 90
    print(f"W={W} H={H}: fits;", ' '.join('(%.2f,%.2f,%.1f°)'%(res[k],res[4+k],min(a[k],90-a[k])) for k in np.argsort(res[4:8])))
