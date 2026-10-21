#!/usr/bin/env python3
"""Check the one-square lemma (2026-10-21): phi(m, H) >= min(2 + m H / 2, Phi(t_m, H)), m = -cot t_m <= 0."""
import math, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from phi import phi

def Phi(t, H):
    C, S = math.cos(t), math.sin(t); k = (C + S) / 2
    return 1 + (0.5 - (H - k - 1) * S) / C + (1 - 2 * (H - k) * C) / (2 * S)

def job(a):
    t, H = a
    m = -1 / math.tan(t)
    p = phi(m, H)
    lb = min(2 + m * H / 2, Phi(t, H))
    return math.degrees(t), H, p, 2 + m * H / 2, Phi(t, H), p - lb

if __name__ == '__main__':
    ts = np.radians([3, 10, 20, 30, 35, 40, 45, 50, 55, 60, 70, 80, 87])
    Hs = [1.0, 1.3, 1.6, 1.8, 1.88, 1.93, 1.97, 1.999]
    with ProcessPoolExecutor(8) as ex:
        worst = math.inf
        for td, H, p, a, b, gap in ex.map(job, [(t, H) for H in Hs for t in ts]):
            worst = min(worst, gap)
            tag = 'TIGHT' if abs(gap) < 1e-9 else ('FAIL' if gap < -1e-9 else '')
            print(f'H={H:.3f} t={td:5.1f} phi={p:.9f} 2+mH/2={a:.6f} Phi={b:.6f} slack={gap:+.2e} {tag}', flush=True)
        print('worst slack', worst)
