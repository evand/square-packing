"""Leaf slack: min side L of 12 unit squares in [0,L]^2 such that squares 0..3 lie wholly in [0,w] x [0,L].
Float penalty multistart; result is an upper bound on the leaf's true minimum."""
import numpy as np, sys
from scipy.optimize import minimize
from rectpack import penalty, gap, half
n = 12
def strip_pen(v, w):
    x, t = v[:4], v[2*n:2*n+4]
    return np.sum(np.maximum(0, x + half(t) - w)**2)
def run(w, rng):
    L = 4.6
    x = np.r_[rng.uniform(0.5, w-0.5, 4), rng.uniform(w, L-0.5, 8)]
    y = rng.uniform(0.5, L-0.5, n)
    t = rng.uniform(0, np.pi/2, n)
    if rng.random() < 0.5: t = np.where(rng.random(n) < 0.6, 0.0, t)
    z = np.r_[x, y, t, L]
    for mu in (1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7):
        f = lambda z: z[-1] + mu*(penalty(z[:-1], n, z[-1], z[-1]) + strip_pen(z[:-1], w))
        z = minimize(f, z, method='L-BFGS-B', options={'maxiter': 20000}).x
    v, L = z[:-1], z[-1]
    for k in range(400):
        r = minimize(lambda vv: penalty(vv, n, L, L) + strip_pen(vv, w), v, method='L-BFGS-B', options={'ftol':1e-22,'gtol':1e-15})
        v = r.x
        if gap(v, n, L, L) > -1e-9 and strip_pen(v, w) < 1e-18: return L, v
        L += 5e-4
    return np.inf, v
w, starts, seed = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed); best = (np.inf, None)
for s in range(starts):
    L, v = run(w, rng)
    if L < best[0]: best = (L, v); np.save(f"leaf_w{w}_s{seed}.npy", v)
print(f"w={w} seed={seed} starts={starts} best L={best[0]:.4f}", flush=True)
