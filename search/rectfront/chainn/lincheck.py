#!/usr/bin/env python3
"""Do piecewise-linear potentials (Lean style) suffice for odd n down to ~1.849?  (2026-10-24)"""
import sys, math, numpy as np
from kwalk import Kmat
phis = np.arange(-89, 89 + 1e-9, 0.25)
ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0; am = np.abs(ms)

def run(H, variants):
    K, _, _ = Kmat(H, ms, 3601)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2; g = a + b; c = K - lam; z = int(np.argmin(am))
    cheap = c < a - 1e-12
    for name, (u0, u1) in variants(a, b).items():
        r = {'W': (c[z] - u1).min(), 'S10': (u1[:, None] + c - u0[None, :]).min(), 'S01': (u0[:, None] + c - u1[None, :]).min(),
             'C0': np.where(cheap, u0[:, None] + c + u0[None, :], np.inf).min() - 2 * g,
             'C1': np.where(cheap, u1[:, None] + c + u1[None, :], np.inf).min() - 2 * g}
        print(f"H={H} {name:34s} " + " ".join(f"{k} {v:+.1e}" for k, v in r.items())
              + ("" if min(r.values()) > -1e-12 else "  FAIL"), flush=True)

def variants(a, b):
    uE = np.minimum(2 * a + b * am, a + b * (2 - am))
    return {'uEven; u1=min(uOdd, 2a+b(2-|m|))': (uE, np.minimum(a + b * am, 2 * a + b * (2 - am)))}

if __name__ == '__main__':
    for H in [float(x) for x in sys.argv[1:]] or [1.849, 1.86, 1.9]:
        run(H, variants)


def Efun(H, x):
    """E as a function of x = |m| > 0 (pcert.py, sympy form)."""
    return -H * x / 2 + 2 * H - H / x + np.sqrt(x * x + 1) - 4 * math.sqrt(2) + 1 + np.sqrt(x * x + 1) / x + 1 / x

def tangents(H, pts):
    """Max of tangent lines of E at |m| = pts (piecewise linear in |m|)."""
    h = 1e-6
    lines = [(Efun(H, p), (Efun(H, p + h) - Efun(H, p - h)) / (2 * h), p) for p in pts]
    return lambda x: np.max([v + s * (x - p) for v, s, p in lines], axis=0)

def run_tan(H, pts):
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2
    L = tangents(H, pts)(am)
    run(H, lambda a_, b_: {f'tangents {pts}': (np.minimum(2 * a + b * am, L), np.minimum(a + b * am, a + L))})


def exact_pots(H):
    """Lean candidate: L(x) = max(a + beta(2 - x), T(x)), T the tangent of E at x = 4/3 (exact, affine in H)."""
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; b = 1 - H / 2
    T = 7 * H / 12 + 14 / 3 - 4 * math.sqrt(2) + (H / 16 - 1 / 10) * (am - 4 / 3)
    L = np.maximum(a + b * (2 - am), T)
    return np.minimum(2 * a + b * am, L), np.minimum(a + b * am, a + L)
