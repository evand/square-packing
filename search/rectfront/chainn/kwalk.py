#!/usr/bin/env python3
"""Slope-walk form of the chain relaxation (2026-10-22).
Separating lines x + m_i y = s_i (m_0 = m_n = 0 are the walls).  Square i between lines i-1 and i:
s_i - s_{i-1} >= (g_t(m_{i-1}) + g_t(m_i))/2 + (m_i - m_{i-1}) y,  y in [k_t, H - k_t]; the H/2 (m_i - m_{i-1}) parts
telescope away, leaving L >= sum_i K(m_{i-1}, m_i),
  K(m, m') = min_t [(g_t(m) + g_t(m'))/2 - |m - m'| (H/2 - k_t)],  g_t(m) = |cos t + m sin t| + |m cos t - sin t|."""
import math, sys
import numpy as np
from targets import target

def Kmat(H, ms, nt=3601):
    ts = np.linspace(0, math.pi / 2, nt)
    C, S = np.cos(ts), np.sin(ts); k = (C + S) / 2
    ok = 2 * k <= H + 1e-15
    C, S, k = C[ok], S[ok], k[ok]
    g = np.abs(C[None, :] + ms[:, None] * S[None, :]) + np.abs(ms[:, None] * C[None, :] - S[None, :])   # (M, T)
    K = np.full((len(ms), len(ms)), np.inf); arg = np.zeros((len(ms), len(ms)), int)
    for i in range(len(ms)):
        v = (g[i][None, :] + g) / 2 - np.abs(ms[i] - ms)[:, None] * (H / 2 - k)[None, :]
        K[i] = v.min(1); arg[i] = v.argmin(1)
    return K, ts[ok], arg

def walk(H, nmax, ms):
    K, ts, arg = Kmat(H, ms)
    z = int(np.argmin(np.abs(ms)))
    V = K[z].copy(); out = {1: K[z, z]}
    for n in range(2, nmax + 1):
        out[n] = float(V[z] if False else (V + K[:, z]).min())
        V = (V[:, None] + K).min(0)
    return out

if __name__ == '__main__':
    phis = np.radians(np.linspace(-85, 85, 681))      # 0.25 deg; includes 0, +-45
    ms = np.tan(phis); ms[np.abs(ms) < 1e-12] = 0
    for H in [float(a) for a in sys.argv[1:]] or [1.7, 1.83, 1.85, 1.87, 1.9, 1.93, 1.96, 1.99]:
        out = walk(H, 30, ms)
        g = {n: out[n] - target(n, H) for n in range(3, 31)}
        print(f"H={H:.3f} min {min(g.values()):+.2e} (n={min(g, key=g.get)}) max {max(g.values()):+.2e} | "
              + " ".join(f"{n}:{g[n]:+.4f}" for n in (3, 4, 5, 6, 7, 8, 9, 12, 20, 30)), flush=True)
