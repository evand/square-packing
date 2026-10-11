import numpy as np, sys
from scipy.optimize import minimize
from notch import P, ok
n=8; H, c, W, starts = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
base = np.load('cfg_n8_H3.99.npy'); bx, by, bt = base[:8], base[8:16], base[16:]
bx = bx - (bx - 0.5*(np.abs(np.cos(bt))+np.abs(np.sin(bt)))).min()
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv)>5 else 0)
best = None
for s in range(starts):
    fx, fy = rng.random() < 0.5, rng.random() < 0.5
    x = (bx.max()+bx.min()) - bx if fx else bx.copy(); y = H - by*(H/3.99) if fy else by*(H/3.99)
    t = -bt if fx ^ fy else bt.copy()
    x = x * (W / 2.8625)  # stretch to the width
    x += rng.normal(0, 0.15, 8); y += rng.normal(0, 0.15, 8); t = t + rng.normal(0, 0.15, 8)
    v = np.r_[x, y, t]
    for s_ in (1.1, 1.04, 1.0):
        v = minimize(lambda vv: P(vv, n, W*s_, H*s_, c*s_) if s_ > 1 else P(vv, n, W, H, c), v, method='L-BFGS-B', options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
    if ok(v, n, W, H, c):
        print(f"H={H} c={c} W={W}: FITS (start {s})", flush=True); np.save(f"notchfit_H{H}_c{c}_W{W}.npy", v); sys.exit()
print(f"H={H} c={c} W={W}: no fit in {starts} seeded starts", flush=True)
