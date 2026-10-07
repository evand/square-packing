"""Condensed-seam profile (FRIEDMAN.md section 8): rho = Leb off the strips (-g/2, g/2) + Z, plus g * delta_0 (the seam line).
Float check of e_needed(g) = 1 - min_{th, cx} F(th, cx);  heuristic prediction e ~ g^1.5 / 4 (worst tilt ~ sqrt(g)).
F(th, cx) = 1 - int_{-g/2}^{g/2} G_th(x - cx) dx + g G_th(-cx),  G = periodised vertical chord V_th."""
import sys, numpy as np
def V(c, s, d):
    d = np.abs(d); a = (c + s) / 2
    if s < 1e-15: return (d <= 0.5).astype(float)
    return np.clip(np.minimum(1 / c, (a - d) / (s * c)), 0, None)
def Vint(c, s, lo, hi):          # exact integral of V over [lo, hi] (piecewise linear: integrate on breakpoints)
    a = (c + s) / 2; b = (c - s) / 2
    pts = np.unique(np.r_[lo, hi, [p for p in (-a, -b, b, a) if lo < p < hi]])
    v = V(c, s, pts); return float(np.sum((v[1:] + v[:-1]) / 2 * np.diff(pts)))
def F(th, cx, g):
    c, s = np.cos(th), np.sin(th); tot = 1.0
    for j in (-1, 0, 1):
        tot += -Vint(c, s, j - g / 2 - cx, j + g / 2 - cx) + g * float(V(c, s, np.array([j - cx]))[0])
    return tot
def worst(g, nth=400, ncx=400):
    best = (9, 0, 0)
    ths = np.unique(np.r_[np.geomspace(1e-5, np.pi / 4, nth), np.sqrt(g) * np.linspace(0.2, 3, 60)])
    for th in ths:
        c, s = np.cos(th), np.sin(th); a = (c + s) / 2; b = (c - s) / 2
        # critical cx: kinks of V aligned with strip ends / seam images, plus a grid
        crit = np.r_[np.linspace(0, 1, ncx, endpoint=False),
                     [(sgn * k + off) % 1 for k in (a, b) for sgn in (1, -1) for off in (0, g / 2, -g / 2, 1e-12, -1e-12)]]
        for cx in crit:
            f = F(th, cx, g)
            if f < best[0]: best = (f, th, cx)
    # local refine
    f0, th0, cx0 = best
    for _ in range(3):
        for th in np.clip(th0 * np.linspace(0.9, 1.1, 41), 1e-9, np.pi / 4):
            for cx in cx0 + np.linspace(-2e-3, 2e-3, 81):
                f = F(th, cx, g)
                if f < f0: f0, th0, cx0 = f, th, cx
    return 1 - f0, th0, cx0
for g in [float(t) for t in (sys.argv[1].split(',') if len(sys.argv) > 1 else ['0.3', '0.1', '0.03', '0.01', '0.003'])]:
    e, th, cx = worst(g)
    print(f"g={g:<7} e_needed={e:.3e}  e/g^1.5={e / g**1.5:.4f}  e/g={e / g:.4f}  worst th={np.degrees(th):.3f}deg (th/sqrt g={th / np.sqrt(g):.2f}) cx={cx % 1:.4f}", flush=True)
