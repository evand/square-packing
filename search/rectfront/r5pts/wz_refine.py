import numpy as np
from scipy.optimize import minimize
from wz_check import *
def f(p, lem):
    a, th, x, y = p
    a = np.clip(a, 1+R2, 4-R2); th = np.clip(th, 0, np.pi/2); b = 4+R2-a; k = (np.cos(th)+np.sin(th))/2
    lo, hi = (k, 1) if lem == 'W' else (1, a/2)
    x2 = np.clip(x, lo, hi); y2 = np.clip(y, max(1, k), b/2)
    pen = abs(x-x2) + abs(y-y2) + abs(p[0]-a) + abs(p[1]-th)
    return -(val(lem, a, th, x2, y2) - 10*pen)
rng = np.random.default_rng(1)
for lem in 'WZ':
    a, th, x, y = sample(lem, 2000000, rng); v = val(lem, a, th, x, y)
    idx = np.argsort(-v)[:200]
    # also spread: random 200
    idx = np.concatenate([idx, rng.integers(0, len(v), 200)])
    best = -9; bp = None
    for i in idx:
        r = minimize(f, [a[i], th[i], x[i], y[i]], args=(lem,), method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
        if -r.fun > best: best, bp = -r.fun, r.x
    print(lem, 'refined worst', best, 'a=%.5f th=%.4fdeg x=%.5f y=%.5f' % (bp[0], np.degrees(bp[1]), bp[2], bp[3]))
