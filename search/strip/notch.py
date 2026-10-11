"""Min width W of region [0,W]x[0,H] minus notch [0,c]x[H-1,H] holding n unit squares (float, upper bound)."""
import numpy as np, sys
from scipy.optimize import minimize
from rectpack import penalty, gap, half

def obst_pen(v, n, c, H, sign=1):
    x, y, t = v[:n], v[n:2*n], v[2*n:]
    ox, oy, hx, hy = c/2, H-0.5, c/2, 0.5
    dx, dy = x-ox, y-oy
    best = np.full(n, np.inf)
    for ph in (np.zeros(n), np.full(n, np.pi/2), t, t+np.pi/2):
        u = np.abs(np.cos(ph)*dx + np.sin(ph)*dy)
        r = half(t-ph) + hx*np.abs(np.cos(ph)) + hy*np.abs(np.sin(ph))
        best = np.minimum(best, r-u)
    return best  # >0 = overlap depth

def P(v, n, W, H, c):
    return penalty(v, n, W, H) + np.sum(np.maximum(0, obst_pen(v, n, c, H))**2)

def ok(v, n, W, H, c):
    return gap(v, n, W, H) > -1e-10 and np.max(obst_pen(v, n, c, H)) < 1e-10

def feasible(n, W, H, c, starts, rng):
    for s in range(starts):
        v = np.r_[rng.uniform(0.5, W-0.5, n), rng.uniform(0.5, H-0.5, n), rng.uniform(0, np.pi/2, n)]
        if rng.random() < 0.4: v[2*n:] = np.where(rng.random(n) < 0.5, 0.0, rng.uniform(0, np.pi/2, n))
        for s_ in (1.12, 1.05, 1.0):
            v = minimize(lambda vv: P(vv, n, W*s_, H*s_, c*s_) if s_ > 1 else P(vv, n, W, H, c), v, method='L-BFGS-B',
                         options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
        if ok(v, n, W, H, c): return True, v
    return False, None

if __name__ == '__main__':
    n, H, c, lo, hi, starts = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
    rng = np.random.default_rng(7)
    f, v = feasible(n, hi, H, c, starts, rng)
    if not f: print(f"n={n} H={H} c={c}: none even at W={hi}", flush=True); sys.exit()
    while hi-lo > 4e-3:
        m = (lo+hi)/2; f, vv = feasible(n, m, H, c, starts, rng)
        if f: hi, v = m, vv
        else: lo = m
    np.save(f"notch_n{n}_H{H}_c{c}.npy", v)
    print(f"n={n} H={H} c={c}: W_min <= {hi:.4f} (none at {lo:.4f});  kill threshold H-1 = {H-1:.3f}  => {'NO KILL (fits)' if hi <= H-1 else 'kill if search is right'}", flush=True)
