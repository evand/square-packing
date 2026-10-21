#!/usr/bin/env python3
"""R_5 strip, plan B: frontier branch with csc t >= (1 + cot t)/sqrt2, sec t >= (1 + tan t)/sqrt2 (tangent at 45).
X_relax(c) = min over lines of min(upright branch, relaxed frontier branch) + c's exact support.
Check the final slack and where the minimising line is (edge contact / vertical / interior)."""
import math
import numpy as np
from xsplit5 import target

R2 = math.sqrt(2)
def E_relax(t, eps, th, y, H):
    C, S = math.cos(t), math.sin(t); v = C if eps else -C
    cot, tan = C / S, S / C
    up = 2 * S - C * (H - 1)
    # frontier branch, divided by S:  Psi/S = 1 - cot + (1 - H)(tan + cot) + csc + sec   (Psi = (C(S-C) + 1 + C + S - H)/C)
    psiS = 1 - cot + (1 - H) * (tan + cot) + (1 + cot) / R2 + (1 + tan) / R2
    P, Q = math.cos(th), math.sin(th); n = (S, v)
    h = abs(n[0] * P + n[1] * Q) + abs(-n[0] * Q + n[1] * P)
    rest = ((C * H if eps else 0) + h / 2 - v * y) / S
    return min(up / S, psiS) + rest

def X_relax(th, y, H):
    best = (math.inf, None)
    ts = np.radians(np.linspace(0.02, 89.98, 8997))
    for eps in (False, True):
        vals = np.array([E_relax(t, eps, th, y, H) for t in ts])
        i = int(np.argmin(vals))
        if vals[i] < best[0]: best = (vals[i], math.degrees(ts[i]), eps)
    return best

if __name__ == '__main__':
    for H in (1.0, 1.3, 1.6, 1.71, 1.8, 1.8713, 1.9, 1.95, 1.99):
        worst = (math.inf, None); interior = []
        for wall in ('floor', 'ceil'):
            for thd in np.linspace(0, 90, 91):
                th = math.radians(thd); k = (math.cos(th) + math.sin(th)) / 2
                if 2 * k > H: continue
                y = k if wall == 'floor' else H - k
                a, b = X_relax(th, y, H), X_relax(math.pi / 2 - th, y, H)
                for (val, tdeg, eps), tt in ((a, thd), (b, 90 - thd)):
                    kinks = {round(tt % 90, 2), round((90 - tt) % 90, 2), 89.98}
                    if not any(abs(tdeg - kk) < 0.03 for kk in kinks): interior.append((round(thd), wall, round(tdeg, 2)))
                s = a[0] + b[0] - target(H)
                if s < worst[0]: worst = (s, thd, wall)
        print(f'H={H}: worst slack {worst[0]:+.2e} at theta={worst[1]} {worst[2]};  interior minimisers: {len(interior)} {interior[:4]}', flush=True)
