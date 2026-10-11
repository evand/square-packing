import numpy as np, sys
from scipy.optimize import minimize
from rectpack import penalty, gap
W, H = float(sys.argv[1]), float(sys.argv[2])
rng = np.random.default_rng(int(sys.argv[3])); found = []
for s in range(150):
    v = np.r_[np.array([0.5, W-0.5, 0.5, W-0.5]) + rng.normal(0, 0.05, 4), (np.arange(4)+0.5)*H/4 + rng.normal(0, 0.03, 4), rng.normal(0, np.radians(6), 4)]
    v = minimize(lambda vv: penalty(vv, 4, W, H), v, method='L-BFGS-B', options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
    g = gap(v, 4, W, H)
    if g > -1e-10:
        a = np.degrees(v[8:]) % 90; a = np.minimum(a, 90-a)
        found.append(a.max())
        if len(found) == 1: print("feasible:", np.round(v[:4],3), np.round(v[4:8],3), np.round(a,2))
print(f"W={W} H={H}: {len(found)} feasible of 150; max tilt among feasible: min {min(found) if found else None}")
