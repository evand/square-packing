#!/usr/bin/env python3
"""n = 7 strip, band H*(7) = (10√2 - 3)/6 <= H < 1.87 (2026-10-24).

Below 1.87 `KapC` fails, so `strip_odd_lower` stops.  For n = 7 only this band is missing (below H*(7) monotonicity
gives L >= 7).  Questions, on the slope grid of kwalk (0.25 deg, fine tilt grid):
 1. how far KapC fails (min over steps of kap - (2 lam - 1 + beta(|m'| - |m|)));
 2. does the 7-step walk still reach F_7 = 5 lam + 4 - H (relaxation exact)?
 3. exact prefix distances d_j(m) (j steps from the wall, excess c = K - lam) vs the Lean potentials
    uOdd = a + beta|m| (odd j), uEven = min(2a + beta|m|, a + beta(2 - |m|)) (even j): where is d_j < u_j?
 4. index-specific check for n = 7: for a cheap step between lines j, j+1, d_j(m) + c + d_{6-j}(m') >= 2g?
    (the whole walk argument with exact per-index potentials; if this holds, potentials indexed by j <= 3
    suffice, by the mirror)."""
import math, sys
import numpy as np
from kwalk import Kmat

def run(H, phis):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, ts, arg = Kmat(H, ms)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; beta = 1 - H / 2; g = a + beta
    c = K - lam; z = int(np.argmin(np.abs(ms))); am = np.abs(ms)
    # 1. KapC margin (K is min over tilts, so this is the worst tilt)
    kapc = (K - (2 * lam - 1 + beta * (am[None, :] - am[:, None]))).min()
    # 2. 7-step walk
    d = [None, c[z].copy()]
    for j in range(2, 7):
        d.append((d[-1][:, None] + c).min(0))
    walk7 = (d[6] + c[:, z]).min() + 7 * lam
    F7 = 5 * lam + 4 - H
    # 3. prefix distances vs Lean potentials
    uO = a + beta * am
    uE = np.minimum(2 * a + beta * am, a + beta * (2 - am))
    gaps = {j: (d[j] - (uO if j % 2 else uE)) for j in range(1, 7)}
    # 4. per-index cheap-step check (lines j, j+1 interior, j = 1..5; suffix length 6 - j)
    cheap = c < a - 1e-12
    worst4 = min(np.where(cheap, d[j][:, None] + c + d[6 - j][None, :], np.inf).min() for j in range(1, 6)) - 2 * g
    return ms, kapc, walk7 - F7, gaps, worst4

if __name__ == '__main__':
    phis = np.linspace(-85, 85, 681)
    H7 = (10 * math.sqrt(2) - 3) / 6
    for H in [float(x) for x in sys.argv[1:]] or [H7, 1.86, 1.865, 1.8685, 1.87]:
        ms, kapc, w7, gaps, w4 = run(H, phis)
        print(f"H={H:.5f}: KapC {kapc:+.2e}  walk7-F7 {w7:+.2e}  per-index cheap {w4:+.2e}", flush=True)
        for j, gp in gaps.items():
            i = int(np.argmin(gp))
            print(f"    d_{j} - u: min {gp[i]:+.3e} at phi {math.degrees(math.atan(ms[i])):+.2f}", flush=True)


# --- 2026-10-24 (later): closed-form line-3 potential and the n = 7 walk conditions ---------------------------
def kap(C, S, H, m, mp):
    """Lean `kap` (StripKernel): cost of a square with frame (C, S) between lines of slopes m, m'."""
    return (abs(C + m * S) + abs(m * C - S) + abs(C + mp * S) + abs(mp * C - S)) / 2 - abs(m - mp) * (H / 2 - (C + S) / 2)

def u3_aligned(H, ms):
    """Route upright, upright, aligned: line 2 at angle phi = psi +- 90 deg (perpendicular to line 3, slope tan psi),
    square 3 at tilt phi, edges along both lines.  Value uEven(tan phi) + kap - lam."""
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; beta = 1 - H / 2
    out = np.full(len(ms), np.inf)
    for i, m in enumerate(ms):
        psi = math.atan(m)
        if abs(m) < 1e-12: continue
        phi = psi - math.copysign(math.pi / 2, psi)          # perpendicular, in (-90, 90)
        m2 = math.tan(phi); t = phi % (math.pi / 2)
        C, S = math.cos(t), math.sin(t)
        if C + S > H: continue
        uE2 = min(2 * a + beta * abs(m2), a + beta * (2 - abs(m2)))
        out[i] = uE2 + kap(C, S, H, m2, m) - lam
    return out

def conds7(H, phis, u4=None):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, ts, arg = Kmat(H, ms)
    lam = 2 * math.sqrt(2) - H; a = 1 - lam; beta = 1 - H / 2; g = a + beta
    c = K - lam; z = int(np.argmin(np.abs(ms))); am = np.abs(ms)
    uO = a + beta * am; uE = np.minimum(2 * a + beta * am, a + beta * (2 - am))
    u3 = np.minimum(uO, u3_aligned(H, ms))
    d = [None, c[z].copy()]
    for j in range(2, 7): d.append((d[-1][:, None] + c).min(0))
    U = {1: uO, 2: uE, 3: u3, 4: uE if u4 is None else d[4], 5: uO, 6: uE if u4 is None else d[6]}
    res = {'wall': (c[z] - U[1]).min()}
    for j in range(1, 6):
        res[f'step{j}->{j+1}'] = (U[j][:, None] + c - U[j + 1][None, :]).min()
    cheap = c < a - 1e-12
    for j in range(1, 6):   # cheap step between lines j, j+1: prefix u_j, suffix u_{6-j} (mirror)
        res[f'cheap({j},{6-j})'] = np.where(cheap, U[j][:, None] + c + U[6 - j][None, :], np.inf).min() - 2 * g
    res['u3 vs d3 (<=0 ok)'] = (U[3] - d[3]).max()
    return res
