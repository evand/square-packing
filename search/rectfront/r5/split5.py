#!/usr/bin/env python3
"""R_5 strip: does  L >= B_b(line1) + M_c(line1, line2) + B_d(line2)  reach the target? (2026-10-21)
B from the Lean one-square lemma (min of the upright and Psi branches, all four orientations);
M = min over c (tilt, height) of max_c (x + m2 y) - min_c (x + m1 y)."""
import math, sys
import numpy as np
from scipy.optimize import minimize_scalar

def branch(C, S, H):
    """min(Up, Psi) of one_square for a unit direction (C, S); returns the bound on s (Psi = Psic / C)."""
    Up = 2 * S - C * (H - 1)
    Psi = (C * (S - C) + 1 + C + S - H) / C if C > 1e-12 else math.inf
    return min(Up, Psi)

def Bb(t, sgn, H):
    """lower bound on sigma1 = s1/S1, line1 slope m1 = sgn * C/S (sgn=-1: osl_set, +1: osl_set_up)."""
    C, S = math.cos(t), math.sin(t)
    b = branch(C, S, H)
    return b / S if sgn < 0 else (b + C * H) / S

def Bd(t, sgn, H):
    """lower bound on L - sigma2, line2 slope m2 = sgn * C/S (sgn=+1: osl_set_mir, -1: osl_set_rot)."""
    C, S = math.cos(t), math.sin(t)
    b = branch(C, S, H)
    return b / S if sgn > 0 else (b + C * H) / S

def width(th, m1, m2, H):
    P, Q = math.cos(th), math.sin(th); k = (P + Q) / 2
    if 2 * k > H: return math.inf
    h = lambda m: abs(P + m * Q) + abs(-Q + m * P)
    return min((m2 - m1) * y for y in (k, H - k)) + (h(m1) + h(m2)) / 2

def M(m1, m2, H):
    ths = np.linspace(0, math.pi / 2, 721)
    v = min(width(th, m1, m2, H) for th in ths)
    return v

def target(H): return min(5.0, 4 + 6 * math.sqrt(2) - 4 * H)

if __name__ == '__main__':
    worst = {}
    ts = np.radians(np.linspace(1, 89, 89))
    for H in [1.5, 1.8, 1.87, 1.9, 1.95, 1.99]:
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                best = (math.inf, None)
                for t1 in ts:
                    m1 = s1 * math.cos(t1) / math.sin(t1); b1 = Bb(t1, s1, H)
                    for t2 in ts:
                        m2 = s2 * math.cos(t2) / math.sin(t2)
                        g = b1 + M(m1, m2, H) + Bd(t2, s2, H)
                        if g < best[0]: best = (g, (round(math.degrees(t1), 1), round(math.degrees(t2), 1)))
                print(f'H={H} signs=({s1:+d},{s2:+d}) min G={best[0]:.6f} target={target(H):.6f} gap={best[0] - target(H):+.2e} at t1,t2={best[1]}', flush=True)
