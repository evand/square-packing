#!/usr/bin/env python3
"""R_5 strip, plan B, explicit lower bound for X(c) (2026-10-21).
X(c) >= min( upright-branch candidates (lines along c's edges, vertical line),
             frontier-branch relaxed (h_c >= 1): chord in y of m(y) = min_t [Psi(t) + shift + 1/2 - v y]/S ).
Check the final one-variable inequality  Xlb(theta) + Xlb(mirror) >= target, c on the floor and the ceiling."""
import math
import numpy as np
from scipy.optimize import minimize_scalar
from xsplit5 import E, X, target

R2 = math.sqrt(2)

def EU(t, eps, th, y, H):
    C, S = math.cos(t), math.sin(t); v = C if eps else -C
    P, Q = math.cos(th), math.sin(th); n = (S, v)
    h = abs(n[0] * P + n[1] * Q) + abs(-n[0] * Q + n[1] * P)
    return (2 * S - C * (H - 1) + (C * H if eps else 0) + h / 2 - v * y) / S

def EPsi_relaxed(t, eps, y, H):
    C, S = math.cos(t), math.sin(t); v = C if eps else -C
    psi = (C * (S - C) + 1 + C + S - H) / C
    return (psi + (C * H if eps else 0) + 0.5 - v * y) / S

def up_cands(th, y, H):
    out = []
    for eps in (False, True):
        # lines along c's edges: direction angle th or th + 90 -> (t, eps) with m = -cot t (eps False) or +cot t
        for d in (th, th + math.pi / 2):
            d = d % math.pi
            t, e = (d, False) if d < math.pi / 2 else (math.pi - d, True)
            if 1e-9 < t < math.pi / 2 - 1e-12 and e == eps: out.append(EU(t, eps, th, y, H))
        out.append(EU(math.pi / 2 - 1e-9, eps, th, y, H))       # vertical line
    return min(out)

def m_psi(y, H):
    best = math.inf
    for eps in (False, True):
        r = minimize_scalar(lambda t: EPsi_relaxed(t, eps, y, H), bounds=(1e-4, math.pi / 2 - 1e-4), method='bounded',
                            options={'xatol': 1e-12})
        ts = np.radians(np.linspace(0.05, 89.95, 1799))
        best = min(best, r.fun, min(EPsi_relaxed(t, eps, y, H) for t in ts))
    return best

def Xlb(th, y, H, ylo, yhi, mlo, mhi):
    chord = mlo + (mhi - mlo) * (y - ylo) / (yhi - ylo)
    return min(up_cands(th, y, H), chord)

if __name__ == '__main__':
    for H in (1.6, 1.8, 1.8713, 1.9, 1.95, 1.99):
        worst = (math.inf, None)
        for wall in ('floor', 'ceil'):
            ylo, yhi = (0.5, R2 / 2) if wall == 'floor' else (H - R2 / 2, H - 0.5)
            mlo, mhi = m_psi(ylo, H), m_psi(yhi, H)
            for thd in np.linspace(0, 90, 181):
                th = math.radians(thd); k = (math.cos(th) + math.sin(th)) / 2
                if 2 * k > H: continue
                y = k if wall == 'floor' else H - k
                a = Xlb(th, y, H, ylo, yhi, mlo, mhi)
                b = Xlb(math.pi / 2 - th, y, H, ylo, yhi, mlo, mhi)
                s = a + b - target(H)
                if s < worst[0]: worst = (s, thd, wall, a, b, X(th, y, H)[0])
        print(f'H={H}: worst slack {worst[0]:+.2e} at theta={worst[1]} {worst[2]}  (Xlb {worst[3]:.4f}, mirror {worst[4]:.4f}, true X {worst[5]:.4f})', flush=True)
