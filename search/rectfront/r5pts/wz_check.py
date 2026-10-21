"""Numerical check of lemmas W and Z (cover C, c0 = 0.91): worst case of min over options of max gauge - 1/2."""
import numpy as np, sys
from scipy.optimize import minimize
R2 = np.sqrt(2); C0 = 0.91
def g(px, py, x, y, c, s):
    dx, dy = px - x, py - y
    return np.maximum(abs(dx*c + dy*s), abs(-dx*s + dy*c)) - 0.5
def val(lem, a, th, x, y):
    b = 4 + R2 - a; c, s = np.cos(th), np.sin(th)
    r = g(C0, 1, x, y, c, s); S = g(1, 1, x, y, c, s); m = g(1, b/2, x, y, c, s)
    rh = g(C0, b-1, x, y, c, s); sh = g(1, b-1, x, y, c, s); mp = g(a-1, b/2, x, y, c, s)
    if lem == 'W':
        opts = [np.maximum(r, S), np.maximum(r, m), np.maximum(rh, m)]
    else:
        opts = [np.maximum(r, S), np.maximum(np.maximum(S, m), mp), np.maximum(np.maximum(sh, m), mp)]
    return np.minimum.reduce(opts)
def sample(lem, n, rng):
    a = rng.uniform(1+R2, 4-R2, n); b = 4 + R2 - a
    th = rng.uniform(0, np.pi/2, n); k = (np.cos(th)+np.sin(th))/2
    if lem == 'W':
        x = rng.uniform(k, 1, n)
    else:
        x = rng.uniform(1, a/2, n)
    y = rng.uniform(np.maximum(1, k), b/2, n)
    return a, th, x, y
if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for lem in 'WZ':
        worst = -9; arg = None
        for rep in range(40):
            a, th, x, y = sample(lem, 500000, rng)
            v = val(lem, a, th, x, y)
            i = np.argmax(v)
            if v[i] > worst: worst, arg = v[i], (a[i], th[i], x[i], y[i])
        print(lem, 'worst sampled', worst, 'at a=%.4f th=%.2fdeg x=%.4f y=%.4f' % (arg[0], np.degrees(arg[1]), arg[2], arg[3]))
