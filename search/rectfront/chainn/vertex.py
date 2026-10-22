#!/usr/bin/env python3
"""Exact-in-slopes check of the per-square E5 inequalities: for fixed tilt t and height H everything is piecewise
linear in (m, m'), so the minimum over a polyhedral region is at a vertex of the arrangement (kinks m, m' in
{0, tan t, -cot t}, the diagonal m = m', and the region's boundary lines).  Grid over (t, H) only.
Reports the minimum and the binding vertex for each statement."""
import math, itertools, sys
import numpy as np

def kap(t, H, m, mp):
    C, S = math.cos(t), math.sin(t); k = (C + S) / 2
    g = lambda m: abs(C + m * S) + abs(m * C - S)
    return (g(m) + g(mp)) / 2 - abs(m - mp) * (H / 2 - k)

def lines_for(t, extra):
    """Lines a*m + b*m' = c of the arrangement."""
    ks = [0.0, math.tan(t)] + ([-1 / math.tan(t)] if t > 1e-9 else [])
    L = [(1, 0, x) for x in ks] + [(0, 1, x) for x in ks] + [(1, -1, 0.0)] + extra
    return L

def vertices(L):
    V = []
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(L, 2):
        d = a1 * b2 - a2 * b1
        if abs(d) < 1e-12: continue
        V.append(((c1 * b2 - c2 * b1) / d, (a1 * c2 - a2 * c1) / d))
    return V

def stmt(name, H, t):
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2
    if name == "A1":   # |m|+|m'| >= tau  =>  c >= b(2 - |m| - |m'|)
        tau = 2 - a / b
        extra = [(1, 1, tau), (1, 1, -tau), (1, -1, tau), (1, -1, -tau)]
        reg = lambda m, mp: abs(m) + abs(mp) >= tau - 1e-12
        f = lambda m, mp: kap(t, H, m, mp) - lam - b * (2 - abs(m) - abs(mp))
    elif name == "A2":  # |m|+|m'| <= tau  =>  kap >= 1
        tau = 2 - a / b
        extra = [(1, 1, tau), (1, 1, -tau), (1, -1, tau), (1, -1, -tau)]
        reg = lambda m, mp: abs(m) + abs(mp) <= tau + 1e-12
        f = lambda m, mp: kap(t, H, m, mp) - 1
    elif name == "S1a":  # |m'| >= tau'  =>  c >= b(2 - |m| - |m'|)
        tau = 1 - a / (2 * b)
        extra = [(0, 1, tau), (0, 1, -tau)]
        reg = lambda m, mp: abs(mp) >= tau - 1e-12
        f = lambda m, mp: kap(t, H, m, mp) - lam - b * (2 - abs(m) - abs(mp))
    elif name == "S1b":  # |m'| <= tau'  =>  c >= a + b(|m'| - |m|)
        tau = 1 - a / (2 * b)
        extra = [(0, 1, tau), (0, 1, -tau)]
        reg = lambda m, mp: abs(mp) <= tau + 1e-12
        f = lambda m, mp: kap(t, H, m, mp) - lam - a - b * (abs(mp) - abs(m))
    elif name == "Cp":   # c >= b(|m'| - |m|) - a
        extra = []
        reg = lambda m, mp: True
        f = lambda m, mp: kap(t, H, m, mp) - lam - (b * (abs(mp) - abs(m)) - a)
    V = [v for v in vertices(lines_for(t, extra)) if reg(*v)]
    if not V: return math.inf, None
    vals = [(f(*v), v) for v in V]
    return min(vals)

if __name__ == '__main__':
    Hlo = float(sys.argv[1]) if len(sys.argv) > 1 else 1.87
    for name in ("A1", "A2", "S1a", "S1b", "Cp"):
        best = (math.inf, None)
        for H in np.linspace(Hlo, 1.9999, 60):
            for t in np.radians(np.linspace(0, 89.99, 900)):
                if math.cos(t) + math.sin(t) > H: continue
                v, x = stmt(name, H, t)
                if v < best[0]:
                    best = (v, (round(math.degrees(t), 2), round(H, 4),
                                None if x is None else tuple(round(math.degrees(math.atan(y)), 2) for y in x)))
        print(f"{name}: min {best[0]:+.3e} at (t, H, (phi, phi')) = {best[1]}", flush=True)
