"""(S), opposite slopes m >= 0 >= m' = -y: decoupled per-line checks.
 zig  (m >= mu*):  min_{m >= mu*} T'(m) + min_{y >= 0} T'(-y) >= Z
 grid (m <= mu*):  min_{0 <= m <= mu*} T(m)  + min_{y >= 0} T'(-y) >= 1
T(x) = g/2 - |x|(1-k),  T'(x) = g/2 - |x|(H-1-k),  mu* = (Z-1)/(2b)."""
import math, numpy as np
r2 = math.sqrt(2)
worst = {'zig': (9, None), 'grid': (9, None), 'zig_coupled': (9, None)}
for H in np.linspace(2 * r2 - 1, 1.9999, 80):
    b = 1 - H / 2; Z = 2 + 2 * r2 - 2 * H; mu = (Z - 1) / (2 * b)
    for t in np.linspace(1e-4, math.pi / 2 - 1e-4, 1500):
        C, S = math.cos(t), math.sin(t); k = (C + S) / 2
        if C + S > H: continue
        g = lambda x: abs(C + x * S) + abs(x * C - S)
        T = lambda x: g(x) / 2 - abs(x) * (1 - k)
        Tp = lambda x: g(x) / 2 - abs(x) * (H - 1 - k)
        ylow = min(Tp(0.0), Tp(-C / S))                    # convex PL in y >= 0: breakpoints 0, cot t
        # zig: m >= mu (mu may be < 0 -> m >= 0)
        lo = max(mu, 0.0)
        cands = [lo] + ([S / C] if S / C >= lo else [])
        mz = min(Tp(x) for x in cands)
        v = mz + ylow - Z
        if v < worst['zig'][0]: worst['zig'] = (v, (round(math.degrees(t), 2), round(H, 4)))
        if mu > 0:
            cands = [0.0, mu] + ([S / C] if S / C <= mu else [])
            mg = min(T(x) for x in cands)
            v = mg + ylow - 1
            if v < worst['grid'][0]: worst['grid'] = (v, (round(math.degrees(t), 2), round(H, 4)))
for k_, v in worst.items(): print(k_, f"{v[0]:+.4e}", v[1])
