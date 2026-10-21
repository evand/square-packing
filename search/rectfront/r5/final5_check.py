#!/usr/bin/env python3
"""Numerical check of the Lean statement Final5 (2026-10-21), with exactly its conventions:
alFr(eps, C, S) = (S, C) if eps else (C, S);  sgnC(eps, C) = C if eps else -C;
OSB: 2S - C(H-1) <= X  or  C(S-C) + 1 + C + S - H <= C X.
Least L allowed by the hypotheses, versus the conclusion 5 <= L or 4 + 6 sqrt2 <= L + 4H."""
import math, sys
import numpy as np

def osb_min(C, S, H):
    up = 2 * S - C * (H - 1)
    psi = (C * (S - C) + 1 + C + S - H) / C if C > 0 else math.inf
    return min(up, psi)

def cr(a, b): return a[0] * b[1] - a[1] * b[0]
def dt(a, b): return a[0] * b[0] + a[1] * b[1]

def midW(mu, c1, c2, a1, a2, x, y):
    return mu * y + c1 * (abs(cr(a1, x)) + dt(a1, x)) + c2 * (abs(cr(a2, x)) + dt(a2, x))

def least_L(e1, e2, t1, t2, H, xi):
    C1, S1, C2, S2 = math.cos(t1), math.sin(t1), math.cos(t2), math.sin(t2)
    a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
    v1 = C1 if e1 else -C1; v2 = C2 if e2 else -C2
    s1 = osb_min(C1, S1, H) + (C1 * H if e1 else 0)
    x = [(1.0, 0.0), a1, a2, (0.0, 1.0)][xi]
    mu = v2 / S2 - v1 / S1
    y = (x[0] + x[1]) / 2 if mu >= 0 else H - (x[0] + x[1]) / 2
    s2 = S2 * (s1 / S1 + midW(mu, 1 / (2 * S1), 1 / (2 * S2), a1, a2, x, y))
    return (s2 + osb_min(C2, S2, H) + (0 if e2 else C2 * H)) / S2

worst = (math.inf, None)
ts = np.radians(np.concatenate([np.linspace(0.01, 89.99, 400), [45.0, 90 - 1e-6]]))
for H in np.concatenate([np.linspace(1.0, 1.999, 60), [1.8713, 1.8714]]):
    for e1 in (False, True):
        for e2 in (False, True):
            for xi in range(4):
                for t1 in ts:
                    for t2 in ts[::4]:
                        L = least_L(e1, e2, t1, t2, H, xi)
                        slack = max(L - 5, L + 4 * H - 4 - 6 * math.sqrt(2))
                        if slack < worst[0]:
                            worst = (slack, (round(H, 4), e1, e2, xi, round(math.degrees(t1), 3), round(math.degrees(t2), 3)))
print('min slack of the conclusion over the grid:', worst)
