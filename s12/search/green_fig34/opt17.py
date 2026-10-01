"""Float exploration: optimise the 4 shape parameters of the Figure-34 n=17 family to maximise the
critical side t* (smax is 1-homogeneous under scaling, so t* = T17/smax at side T17).  Slow (Nelder-Mead on a
nonsmooth objective, ~2 s per evaluation); one 400-step round from the pixel-read start reached
t* = 4.43935 at (a,b,s,v) = (0.98004, 0.91707, 0.99170, 0.86840), close to Hstar (GREEN_FIG34.md sec. 3)."""
import sys, os, math, time, functools
print = functools.partial(print, flush=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
from smax import smax
from family import green17, T17

def obj(x, nth=45):
    a, b, s, v = x
    P = green17(T17, a * T17, b * T17, s * T17, v * T17)
    val, th, _ = smax(P, T17, nth=nth, nloc=2)
    return val

if __name__ == '__main__':
    t = T17
    # start: a=(t-2.5)/2, s=1, b~0.915, v~0.866, normalised by t
    x0 = np.array([(t - 2.5) / 2, 0.915, 1.0, 0.866]) / t
    print('start t* =', T17 / obj(x0))
    best = None
    for it in range(3):
        r = minimize(obj, x0, method='Nelder-Mead',
                     options=dict(xatol=1e-7, fatol=1e-10, maxiter=400, initial_simplex=None))
        x0 = r.x
        print(it, 't* =', T17 / r.fun, 'params a,b,s,v (unit = 1) =', list(r.x * T17 / r.fun))
    print('dense check t* =', T17 / obj(x0, nth=360))
