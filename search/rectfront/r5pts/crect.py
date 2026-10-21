"""Centre-rectangle lemma, numerics.
G(X,Y) = max over centre m in [0,X]x[0,Y], theta of  min over the 4 corners of gauge_theta(corner - m),
gauge = max(|u|,|v|) in the square's frame.  Lemma <=> G <= 1/2 (closed square contains a corner)."""
import numpy as np, sys
from scipy.optimize import minimize

def gauge(d, th):
    c, s = np.cos(th), np.sin(th)
    u = d[..., 0]*c + d[..., 1]*s; v = -d[..., 0]*s + d[..., 1]*c
    return np.maximum(abs(u), abs(v))

def mingauge(X, Y, x, y, th):
    cs = np.array([[0, 0], [X, 0], [0, Y], [X, Y]], float)
    return min(gauge(c - np.array([x, y]), th) for c in cs)

def G(X, Y, ngrid=41, nth=91, rng=np.random.default_rng(0)):
    xs = np.linspace(0, X, ngrid); ys = np.linspace(0, Y, ngrid); ths = np.linspace(0, np.pi/2, nth)
    Xg, Yg, Tg = np.meshgrid(xs, ys, ths, indexing='ij')
    cs = np.array([[0, 0], [X, 0], [0, Y], [X, Y]], float)
    c, s = np.cos(Tg), np.sin(Tg)
    m = None
    for cx, cy in cs:
        dx, dy = cx - Xg, cy - Yg
        g = np.maximum(abs(dx*c + dy*s), abs(-dx*s + dy*c))
        m = g if m is None else np.minimum(m, g)
    idx = np.argsort(m.ravel())[-20:]
    best = -1; arg = None
    for i in idx:
        p0 = np.array([Xg.ravel()[i], Yg.ravel()[i], Tg.ravel()[i]])
        f = lambda p: -mingauge(X, Y, np.clip(p[0], 0, X), np.clip(p[1], 0, Y), p[2])
        r = minimize(f, p0, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
        if -r.fun > best: best = -r.fun; arg = (np.clip(r.x[0], 0, X), np.clip(r.x[1], 0, Y), np.degrees(r.x[2]) % 90)
    return best, arg

if __name__ == '__main__':
    r2 = np.sqrt(2)
    print("X  Y  G  arg(x,y,deg)   [lemma: G<=1/2 iff X+Y<=sqrt2, X,Y<=1]")
    for X in np.linspace(0.0, 1.0, 11):
        Ybd = min(1.0, r2 - X)
        for Y in [Ybd - 0.01, Ybd, Ybd + 0.01]:
            if Y < 0: continue
            g, a = G(X, Y)
            print(f"{X:.3f} {Y:.3f} {g:.10f} {'IN ' if (X+Y<=r2+1e-12 and max(X,Y)<=1+1e-12) else 'OUT'} {'ok' if (g<=0.5+1e-9)==(X+Y<=r2+1e-12 and max(X,Y)<=1+1e-12) else 'MISMATCH'}  {a[0]:.4f} {a[1]:.4f} {a[2]:.3f}")
