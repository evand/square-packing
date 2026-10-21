#!/usr/bin/env python3
"""R_5 strip, plan B: X(c) >= min( edge contacts (lines along c's edges), vertical line,
vertex contacts Lambda(y_V) - ox_V over c's three left-facing vertices ), Lambda(y) = min over lines of the
line's x at height y (both branches of b's bound, both signs).  Compare with the true X and the final slack."""
import math
import numpy as np
from scipy.optimize import minimize_scalar
from xsplit5 import E, X, target, osb

R2 = math.sqrt(2)

def line_x(t, eps, y, H):
    """x where the b|c line (unit direction t, sign eps) at its least offset meets height y."""
    C, S = math.cos(t), math.sin(t); v = C if eps else -C
    return (osb(C, S, H) + (C * H if eps else 0) - v * y) / S

def Lam(y, H):
    best = math.inf
    for eps in (False, True):
        ts = np.radians(np.linspace(0.02, 89.98, 4497))
        vals = np.array([line_x(t, eps, y, H) for t in ts])
        i = int(np.argmin(vals))
        r = minimize_scalar(lambda t: line_x(t, eps, y, H), bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                            method='bounded', options={'xatol': 1e-12})
        best = min(best, r.fun, vals.min())
    return best

def Xlb(th, y, H):
    P, Q = math.cos(th), math.sin(th)
    vals = []
    # edge contacts: lines along c's edges (directions th and th + 90), both signs; and the vertical line
    for d in (th, th + math.pi / 2):
        d = d % math.pi
        t, eps = (d, False) if d < math.pi / 2 else (math.pi - d, True)
        if 1e-9 < t < math.pi / 2 - 1e-12: vals.append(E(t, eps, th, y, H))
    for eps in (False, True): vals.append(E(math.pi / 2 - 1e-9, eps, th, y, H))
    # vertex contacts: vertices c + (+-P/2 -+ Q/2, ...) ; the left-facing ones: all but the rightmost
    verts = [(sx * P / 2 - sy * Q / 2, sx * Q / 2 + sy * P / 2) for sx in (-1, 1) for sy in (-1, 1)]
    verts.sort(key=lambda o: o[0])
    for ox, oy in verts[:3]:
        vals.append(Lam(y + oy, H) - ox)
    return min(vals)

if __name__ == '__main__':
    for H in (1.8, 1.8713, 1.9, 1.95, 1.99):
        worst = (math.inf, None); gap = 0
        for wall in ('floor', 'ceil'):
            for thd in np.linspace(0, 90, 46):
                th = math.radians(thd); k = (math.cos(th) + math.sin(th)) / 2
                if 2 * k > H: continue
                y = k if wall == 'floor' else H - k
                a, b = Xlb(th, y, H), Xlb(math.pi / 2 - th, y, H)
                xt = X(th, y, H)[0]
                gap = max(gap, xt - a)
                s = a + b - target(H)
                if s < worst[0]: worst = (s, thd, wall)
        print(f'H={H}: worst slack {worst[0]:+.2e} at theta={worst[1]} {worst[2]};  max (true X - bound) = {gap:.3f}', flush=True)
