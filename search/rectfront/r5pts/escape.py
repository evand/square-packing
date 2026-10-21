"""Escape search: max over unit squares (centre m, angle th) inside [0,a]x[0,b] of
   esc = min_points gauge(pt - m) - 1/2   (esc <= 0  <=>  every such square contains a point, closed)."""
import numpy as np
from scipy.optimize import minimize

def esc_grid(pts, a, b, nx=60, ny=60, nth=91):
    ths = np.linspace(0, np.pi/2, nth, endpoint=False)
    best = (-9, None)
    out = []
    for th in ths:
        c, s = np.cos(th), np.sin(th); k = (abs(c)+abs(s))/2
        if 2*k > min(a, b): continue
        xs = np.linspace(k, a-k, nx); ys = np.linspace(k, b-k, ny)
        Xg, Yg = np.meshgrid(xs, ys, indexing='ij')
        m = np.full(Xg.shape, np.inf)
        for px, py in pts:
            dx, dy = px-Xg, py-Yg
            m = np.minimum(m, np.maximum(abs(dx*c+dy*s), abs(-dx*s+dy*c)))
        i = np.unravel_index(np.argmax(m), m.shape)
        out.append((m[i]-0.5, Xg[i], Yg[i], th))
    out.sort(key=lambda t: -t[0])
    return out

def esc_fun(p, pts, a, b):
    th = p[2]; c, s = np.cos(th), np.sin(th); k = (abs(c)+abs(s))/2
    # map p0,p1 in R to [k,a-k] via clip; penalise outside
    x = np.clip(p[0], k, a-k); y = np.clip(p[1], k, b-k)
    pen = abs(p[0]-x)+abs(p[1]-y)
    g = min(max(abs((px-x)*c+(py-y)*s), abs(-(px-x)*s+(py-y)*c)) for px, py in pts)
    return g - 0.5 - pen

def escape(pts, a, b, nstart=12, **kw):
    cand = esc_grid(pts, a, b, **kw)[:nstart]
    q = np.pi/4; h = np.sqrt(2)/2
    cand += [(0, a/2, b/2, q), (0, h, b/2, q), (0, a-h, b/2, q), (0, a/2, h, q), (0, a/2, b-h, q)]
    best = (-9, None)
    for e0, x, y, th in cand:
        r = minimize(lambda p: -esc_fun(p, pts, a, b), [x, y, th], method='Nelder-Mead',
                     options=dict(xatol=1e-11, fatol=1e-13, maxiter=3000))
        if -r.fun > best[0]: best = (-r.fun, r.x)
    return best

def rect_pts(a, b, p=1.0, q=1.0):
    return [(p, q), (a-p, q), (p, b-q), (a-p, b-q)]
