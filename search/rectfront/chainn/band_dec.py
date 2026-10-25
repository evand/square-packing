#!/usr/bin/env python3
"""Band theorem (odd n, 1.849 <= H < 1.87): the per-square conditions as separable convex sub-inequalities (2026-10-24).

Opposite slopes x = m >= 0, y = -m' >= 0, frame (C, S) >= 0, p = C + S <= H:
    kap(x, -y) = p + (S - H/2) x + (C - H/2) y + (Cx - S)^+ + (Sy - C)^+          (closed form; checked below)
Potentials (x = |m|): A0 = 2a + bx, Lz = a + b(2 - x), T = T0 + T1 x (tangent of E at 4/3), L = max(Lz, T),
    u0 = min(A0, L), u1 = min(a + bx, a + L).
A min on the + side splits into both branches (each must hold everywhere); on the - side a branch is chosen per region
(split points y1 for -u1, y2 for -u0; L = Lz below yLT, T above).  Sub-inequalities (c = kap - lam):
    Sa  c >= u1(y) - u0(x):  x-branch in {A0, L};  y: [0, y1] with a + by,  [y1, inf) with a + Lz and with a + T
    Sb  c >= u0(y) - u1(x):  x-branch in {a + bx, a + L};  y: [0, y2] with A0,  [y2, inf) with Lz and with T
    C0  x, y >= 1/2: c + u0(x) + u0(y) >= 2g:  branches {A0, L} x {A0, L}
    C1  x, y >= 1/2: c + u1(x) + u1(y) >= 2g:  branches {a + bx, a + L} x {a + bx, a + L}
    N   x <= 1/2: kap >= 1
(the x <-> y swapped versions follow by kap_mirror; same-sign slopes by kap_same.)"""
import math, sys
import numpy as np

R2 = math.sqrt(2)

def consts(H):
    lam = 2 * R2 - H; a = 1 - lam; b = 1 - H / 2; g = a + b
    T1 = H / 16 - 1 / 10; T0 = 7 * H / 12 + 14 / 3 - 4 * R2 - T1 * 4 / 3
    return lam, a, b, g, T0, T1

def kap_closed(C, S, H, x, y):
    return C + S + (S - H / 2) * x + (C - H / 2) * y + np.maximum(C * x - S, 0) + np.maximum(S * y - C, 0)

def kap_def(C, S, H, m, mp):
    return (abs(C + m * S) + abs(m * C - S) + abs(C + mp * S) + abs(mp * C - S)) / 2 - abs(m - mp) * (H / 2 - (C + S) / 2)

def run(H, y1, y2, nt=721, xs=None):
    lam, a, b, g, T0, T1 = consts(H)
    xs = np.concatenate([np.linspace(0, 3, 1201), np.linspace(3, 60, 200)]) if xs is None else xs
    X = xs[:, None]; Y = xs[None, :]
    A0 = lambda x: 2 * a + b * x
    Lz = lambda x: a + b * (2 - x)
    T = lambda x: T0 + T1 * x
    L = lambda x: np.maximum(Lz(x), T(x))
    O = lambda x: a + b * x
    yLT = (a + 2 * b - T0) / (b + T1)   # Lz = T: L = Lz below, T above
    worst = {}
    for t in np.linspace(0, math.pi / 2, nt):
        C, S = math.cos(t), math.sin(t)
        if C + S > H: continue
        c = kap_closed(C, S, H, X, Y) - lam
        subs = {}
        for xb, P in (('A0', A0), ('L', L)):
            subs[f'Sa {xb}|O[0,y1]'] = np.where(Y <= y1, c + P(X) - O(Y), np.inf)
            subs[f'Sa {xb}|aLz[y1,yLT]'] = np.where((Y >= y1) & (Y <= max(y1, yLT)), c + P(X) - a - Lz(Y), np.inf)
            subs[f'Sa {xb}|aT[yLT,)'] = np.where(Y >= max(y1, yLT), c + P(X) - a - T(Y), np.inf)
        for xb, P in (('O', O), ('aL', lambda x: a + L(x))):
            subs[f'Sb {xb}|A0[0,y2]'] = np.where(Y <= y2, c + P(X) - A0(Y), np.inf)
            subs[f'Sb {xb}|Lz[y2,yLT]'] = np.where((Y >= y2) & (Y <= max(y2, yLT)), c + P(X) - Lz(Y), np.inf)
            subs[f'Sb {xb}|T[yLT,)'] = np.where(Y >= max(y2, yLT), c + P(X) - T(Y), np.inf)
        big = (X >= 0.5) & (Y >= 0.5)
        for xb, P in (('A0', A0), ('L', L)):
            for yb, Q in (('A0', A0), ('L', L)):
                subs[f'C0 {xb}|{yb}'] = np.where(big, c + P(X) + Q(Y) - 2 * g, np.inf)
        for xb, P in (('O', O), ('aL', lambda x: a + L(x))):
            for yb, Q in (('O', O), ('aL', lambda x: a + L(x))):
                subs[f'C1 {xb}|{yb}'] = np.where(big, c + P(X) + Q(Y) - 2 * g, np.inf)
        subs['N'] = np.where(X <= 0.5, c + lam - 1, np.inf)
        for k, M in subs.items():
            i, j = np.unravel_index(np.argmin(M), M.shape)
            v = M[i, j]
            if k not in worst or v < worst[k][0]:
                worst[k] = (v, math.degrees(t), xs[i], xs[j])
    return worst

if __name__ == '__main__':
    # closed form vs definition
    rng = np.random.default_rng(1); err = 0
    for _ in range(20000):
        t = rng.uniform(0, math.pi / 2); C, S = math.cos(t), math.sin(t); H = rng.uniform(1.8, 2)
        x, y = rng.exponential(1.5), rng.exponential(1.5)
        err = max(err, abs(kap_closed(C, S, H, x, y) - kap_def(C, S, H, x, -y)))
    print(f"closed form vs kap: max err {err:.1e}")
    y1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.125
    y2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.875
    for H in (1.849, 1.86, 1.8699):
        w = run(H, y1, y2)
        bad = {k: v for k, v in w.items() if v[0] < -1e-12}
        print(f"H={H} y1={y1} y2={y2}: worst " + "; ".join(f"{k} {v[0]:+.1e}@t{v[1]:.1f},x{v[2]:.2f},y{v[3]:.2f}"
                                                            for k, v in sorted(w.items(), key=lambda kv: kv[1][0])[:6])
              + (f"  FAIL {list(bad)}" if bad else ""), flush=True)
