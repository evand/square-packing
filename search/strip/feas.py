"""Feasibility search: can n unit squares fit strictly in W x H?  Multistart penalty minimisation."""
import numpy as np, sys, os
from scipy.optimize import minimize
from rectpack import penalty, gap, unpack

def attempt(n, W, H, rng):
    v = np.concatenate([rng.uniform(0.5, W-0.5, n), rng.uniform(0.5, H-0.5, n), rng.uniform(0, np.pi/2, n)])
    if rng.random() < 0.4:
        v[2*n:] = rng.choice([0.0, 0.0, np.pi/4, rng.uniform(0, np.pi/2)], n) + rng.normal(0, 0.03, n)
    # shrink-then-target: solve in an inflated box first, then the real one
    for s in (1.15, 1.07, 1.03, 1.0):
        v = minimize(lambda vv: penalty(vv, n, W*s, H*s) if s == 1.0 else penalty(vv*np.r_[np.ones(2*n)*1, np.ones(n)], n, W*s, H*s),
                     v, method='L-BFGS-B', options={'ftol':1e-22,'gtol':1e-15,'maxiter':5000}).x
        if s != 1.0:
            v[:2*n] /= 1.0  # keep
    return v

if __name__ == '__main__':
    n, W, H, starts, seed = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    rng = np.random.default_rng(seed)
    best = -np.inf; bv = None
    for s in range(starts):
        v = attempt(n, W, H, rng)
        g = gap(v, n, W, H)
        if g > best:
            best, bv = g, v
        if g > -1e-10:
            break
    print(f"n={n} W={W} H={H} bestgap={best:.3e} tries={s+1} {'FEASIBLE' if best > -1e-10 else 'not found'}")
    if best > -1e-10:
        x, y, t = unpack(bv, n)
        for k in range(n):
            print(f"  {x[k]:.5f} {y[k]:.5f} {np.degrees(t[k])%90:.3f}")
