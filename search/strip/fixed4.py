"""Given 4 fixed strip squares in [0,H]^2, can 8 more unit squares be added (all disjoint, inside)?"""
import numpy as np, sys
from scipy.optimize import minimize
from rectpack import penalty, gap
H = float(sys.argv[1]); W = float(sys.argv[2])
# third-class strip, re-solved at (W,H)
v = np.load('leaf_w1.99_s2.npy'); n = 12
x, y, t = v[:4], v[n:n+4], v[2*n:2*n+4]
rng = np.random.default_rng(0)
for s in range(60):
    z = np.r_[x + rng.normal(0,0.03,4), y*H/4 + rng.normal(0,0.03,4), t + rng.normal(0,0.03,4)]
    z = minimize(lambda q: penalty(q, 4, W, H), z, method='L-BFGS-B', options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
    if gap(z, 4, W, H) > -1e-10: break
fx, fy, ft = z[:4], z[4:8], z[8:]
base = np.load('cfg_n8_H3.99.npy'); bx, by, bt = base[:8], base[8:16], base[16:]
def full(q): return np.r_[fx, q[:8], fy, q[8:16], ft, q[16:]]
best = -9
for s in range(int(sys.argv[3])):
    if rng.random() < 0.6:
        fl = rng.random() < 0.5
        qy = H - by if fl else by.copy(); qt = -bt if fl else bt.copy()
        q = np.r_[bx + (H - 2.8625) + rng.normal(0, 0.1, 8), qy + rng.normal(0, 0.1, 8), qt + rng.normal(0, 0.1, 8)]
    else:
        q = np.r_[rng.uniform(1.0, H-0.5, 8), rng.uniform(0.5, H-0.5, 8), rng.uniform(0, np.pi/2, 8)]
    for sc in (1.08, 1.03, 1.0):
        q = minimize(lambda qq: penalty(full(qq), 12, H*sc, H*sc), q, method='L-BFGS-B', options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
    g = gap(full(q), 12, H, H); best = max(best, g)
    if g > -1e-10:
        print(f"H={H} strip W={W}: 8 MORE FIT (start {s})"); np.save(f"fixed4_fit_H{H}_W{W}.npy", full(q)); sys.exit()
print(f"H={H} strip W={W}: no fit in {sys.argv[3]} starts, best gap {best:.3e}")
