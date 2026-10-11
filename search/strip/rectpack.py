"""Float search: min width W of a W x H rectangle holding n unit squares (any angle).
Penalty method; feasibility = penetration penalty driven to ~0.  Upper bounds only."""
import numpy as np, sys
from scipy.optimize import minimize

def unpack(v, n):
    return v[:n], v[n:2*n], v[2*n:3*n]

def half(th):
    return 0.5*(np.abs(np.cos(th))+np.abs(np.sin(th)))

def penalty(v, n, W, H, fixed=None):
    x, y, t = unpack(v, n)
    hx = half(t)  # half-extent along x and y are equal for squares: (|c|+|s|)/2
    P = np.sum(np.maximum(0, hx - x)**2) + np.sum(np.maximum(0, x + hx - W)**2)
    P += np.sum(np.maximum(0, hx - y)**2) + np.sum(np.maximum(0, y + hx - H)**2)
    i, j = np.triu_indices(n, 1)
    dx, dy = x[j]-x[i], y[j]-y[i]
    best = np.full(len(i), np.inf)
    for phis in (t[i], t[j]):
        for off in (0.0, np.pi/2):
            ph = phis + off
            u = np.abs(np.cos(ph)*dx + np.sin(ph)*dy)
            r = half(t[i]-ph) + half(t[j]-ph)
            best = np.minimum(best, r - u)
    P += np.sum(np.maximum(0, best)**2)
    return P

def gap(v, n, W, H):
    """min over pairs of separation (SAT: max over axes of u - r) and walls; >0 = strictly disjoint."""
    x, y, t = unpack(v, n)
    hx = half(t)
    g = min(np.min(x-hx), np.min(W-x-hx), np.min(y-hx), np.min(H-y-hx))
    i, j = np.triu_indices(n, 1)
    dx, dy = x[j]-x[i], y[j]-y[i]
    best = np.full(len(i), -np.inf)
    for phis in (t[i], t[j]):
        for off in (0.0, np.pi/2):
            ph = phis + off
            u = np.abs(np.cos(ph)*dx + np.sin(ph)*dy)
            r = half(t[i]-ph) + half(t[j]-ph)
            best = np.maximum(best, u - r)
    return min(g, np.min(best)) if len(i) else g

def minwidth(n, H, rng, W0):
    """penalty ramp on (config, W): minimise W + mu*P"""
    v = np.concatenate([rng.uniform(0.5, W0-0.5, n), rng.uniform(0.5, H-0.5, n), rng.uniform(0, np.pi/2, n)])
    # occasionally axis-parallel start
    if rng.random() < 0.3:
        v[2*n:] = rng.choice([0.0, 0.0, 0.0, np.pi/4], n) + rng.normal(0, 0.05, n)
    W = W0
    z = np.concatenate([v, [W]])
    for mu in (1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8):
        f = lambda z: z[-1] + mu*penalty(z[:-1], n, z[-1], H)
        z = minimize(f, z, method='L-BFGS-B').x
    v, W = z[:-1], z[-1]
    # repair: grow W until strictly feasible after a pure-penalty polish
    for k in range(60):
        r = minimize(lambda vv: penalty(vv, n, W, H), v, method='L-BFGS-B', options={'ftol':1e-20,'gtol':1e-14})
        v = r.x
        if gap(v, n, W, H) > -1e-9:
            return W, v
        W += 2e-4
    return np.inf, v

if __name__ == '__main__':
    n, H, starts, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    W0 = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for s in range(starts):
        W, v = minwidth(n, H, rng, W0)
        if W < best[0]:
            best = (W, v)
    W, v = best
    print(f"n={n} H={H} Wmin~{W:.6f} gap={gap(v,n,W,H):.2e}")
    x, y, t = unpack(v, n)
    for k in range(n):
        print(f"  {x[k]:.5f} {y[k]:.5f} {np.degrees(t[k])%90:.3f}")
