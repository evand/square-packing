#!/usr/bin/env python3
"""band_dec with H-dependent branch splits y1 = 1 + a/(2b) (u1: a + by below, a + L above), y2 = 1 - a/(2b)
(u0: A0 below, L above): the points where the branches cross (2026-10-25)."""
import sys, math
import numpy as np
from band_dec import run, consts

if __name__ == '__main__':
    xs = np.concatenate([np.linspace(0, 3, 601), np.linspace(3.05, 60, 60)])
    for H in [float(h) for h in sys.argv[1:]] or [1.849, 1.86, 1.8699]:
        lam, a, b, g, T0, T1 = consts(H)
        y1 = 1 + a / (2 * b); y2 = 1 - a / (2 * b)
        w = run(H, y1, y2, nt=361, xs=np.unique(np.concatenate([xs, [y1, y2, 0.5, 1.0]])))
        bad = {k: v for k, v in w.items() if v[0] < -1e-12}
        print(f"H={H} y1={y1:.4f} y2={y2:.4f}: " + "; ".join(
            f"{k} {v[0]:+.1e}@t{v[1]:.1f},x{v[2]:.2f},y{v[3]:.2f}"
            for k, v in sorted(w.items(), key=lambda kv: kv[1][0])[:5]) + (f"  FAIL {list(bad)}" if bad else ""),
            flush=True)
