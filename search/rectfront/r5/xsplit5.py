#!/usr/bin/env python3
"""R_5 strip, plan B (2026-10-21): split at the middle square c.
X(theta, y, H) = min over the b|c line (unit direction t, sign eps) of
    [OSB bound for b (+ shift) + h_c(n)/2 - v y] / S         (lower bound for c's centre x)
L >= X(c) + X(mirror c); the sum is concave in y, so c on the floor (y = k) or ceiling (y = H - k).
Report X on the floor / ceiling, the minimising line (aligned with c or not), and the sum vs the target."""
import math
import numpy as np
from scipy.optimize import minimize_scalar

R2 = math.sqrt(2)
def osb(C, S, H):
    up = 2 * S - C * (H - 1)
    psi = (C * (S - C) + 1 + C + S - H) / C if C > 1e-14 else math.inf
    return min(up, psi)

def E(t, eps, th, y, H):
    C, S = math.cos(t), math.sin(t)
    v = C if eps else -C
    n = (S, v)
    P, Q = math.cos(th), math.sin(th)
    h = abs(n[0] * P + n[1] * Q) + abs(-n[0] * Q + n[1] * P)
    return (osb(C, S, H) + (C * H if eps else 0) + h / 2 - v * y) / S

def X(th, y, H):
    best = (math.inf, None)
    for eps in (False, True):
        ts = np.radians(np.linspace(0.05, 89.95, 1799))
        vals = [E(t, eps, th, y, H) for t in ts]
        i = int(np.argmin(vals))
        r = minimize_scalar(lambda t: E(t, eps, th, y, H), bounds=(ts[max(i - 1, 0)], ts[min(i + 1, len(ts) - 1)]),
                            method='bounded', options={'xatol': 1e-12})
        if r.fun < best[0]: best = (r.fun, math.degrees(r.x), eps)
    return best

def target(H): return min(5.0, 4 + 6 * R2 - 4 * H)

if __name__ == '__main__':
    for H in (1.8, 1.8713, 1.9, 1.95, 1.99):
        worst = (math.inf, None)
        rows = []
        for thd in np.linspace(0, 90, 91):
            th = math.radians(thd); k = (math.cos(th) + math.sin(th)) / 2
            if 2 * k > H: continue
            for wall, y in (('floor', k), ('ceil', H - k)):
                xl = X(th, y, H)
                xr = X(math.pi / 2 - th, y, H)          # mirror: tilt 90 - theta, same height
                s = xl[0] + xr[0] - target(H)
                rows.append((thd, wall, xl, xr, s))
                if s < worst[0]: worst = (s, thd, wall)
        print(f'H={H}: worst slack {worst[0]:+.2e} at theta={worst[1]} {worst[2]}')
        for thd, wall, xl, xr, s in rows:
            if thd in (0, 10, 20, 30, 40, 45, 50) and wall == 'floor':
                print(f'   th={thd:4.0f} floor  X={xl[0]:.4f} (line {xl[1]:.2f} eps {int(xl[2])})  Xmir={xr[0]:.4f} (line {xr[1]:.2f})  slack {s:+.3f}')
