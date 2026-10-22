#!/usr/bin/env python3
"""Unified (S): kap >= min(1 + b||m| - |m'||, Z - b(|m| + |m'|)),  Z = 2 + 2sqrt2 - 2H, b = 1 - H/2.
Vertex check (exact in slopes), split by slope signs, with the binding vertex per (sign case, branch)."""
import math, itertools, sys
import numpy as np
from vertex import kap, vertices

def lines_S(t, H):
    b = 1 - H / 2; Z = 2 + 2 * math.sqrt(2) - 2 * H
    ks = [0.0, math.tan(t)] + ([-1 / math.tan(t)] if t > 1e-9 else [])
    L = [(1, 0, x) for x in ks] + [(0, 1, x) for x in ks] + [(1, -1, 0.0), (1, 1, 0.0)]
    # branch-equality lines 1 + b||m|-|m'|| = Z - b(|m|+|m'|): in each sign quadrant, |m|, |m'| linear
    for s1, s2 in itertools.product((1, -1), repeat=2):
        # |m| = s1 m, |m'| = s2 mp;  case |m| >= |m'|:  1 + b(s1 m - s2 mp) = Z - b(s1 m + s2 mp) -> 2 b s1 m = Z - 1
        if b > 1e-12:
            L.append((1, 0, s1 * (Z - 1) / (2 * b))); L.append((0, 1, s2 * (Z - 1) / (2 * b)))
    L += [(1, -1, 0.0), (1, 1, 0.0)]
    return L

def f(t, H, m, mp):
    b = 1 - H / 2; Z = 2 + 2 * math.sqrt(2) - 2 * H
    return kap(t, H, m, mp) - min(1 + b * abs(abs(m) - abs(mp)), Z - b * (abs(m) + abs(mp)))

Hlo = float(sys.argv[1]) if len(sys.argv) > 1 else 1.87
best = {}
for H in np.linspace(Hlo, 1.9999, 50):
    for t in np.radians(np.linspace(0, 89.99, 900)):
        if math.cos(t) + math.sin(t) > H: continue
        for v in vertices(lines_S(t, H)):
            m, mp = v
            if abs(m) > 50 or abs(mp) > 50: continue
            val = f(t, H, m, mp)
            case = ('opp' if m * mp < 0 else 'same')
            b = 1 - H / 2; Z = 2 + 2 * math.sqrt(2) - 2 * H
            br = 'grid' if 1 + b * abs(abs(m) - abs(mp)) <= Z - b * (abs(m) + abs(mp)) else 'zig'
            key = (case, br)
            if val < best.get(key, (math.inf,))[0]:
                best[key] = (val, round(math.degrees(t), 2), round(H, 4), round(math.degrees(math.atan(m)), 2),
                             round(math.degrees(math.atan(mp)), 2))
for k, v in sorted(best.items()): print(k, f"min {v[0]:+.3e}", "t, H, phi, phi' =", v[1:])
