#!/usr/bin/env python3
"""seam_w1_check.py -- FLOAT check of the phase-free sufficient inequality for the w = 1 seam profile (SEAM_W1.md sec 2).

Profile A (band [0,1], 1-periodic, Lebesgue on y > 1):
    seam segments {n} x [1/2, 1] at line density 3/2 (m_v = 3/4),  the line y = 1 at density 1/4.
For a closed unit square S (tilt th in (0, 45 deg], lowest vertex at height y0 = 1 - h, slice width W):
    mu(S) >= 1  <=  F := (1/4) W(at y=1) + (3/2) |{y in [1/2, 1] : W(y) >= 1}| - area(S cap {y <= 1})  >= 0,
using only "a closed slice of width >= 1 contains an integer" (so the horizontal phase drops out).
th = 0 is done by hand in SEAM_W1.md; th in [45, 90) is the mirror image.
Closed-form area G(h) and width; grid over (th, h).  Not a proof: the exact casework is in SEAM_W1.md sec 2.3 (TODO).
"""
import numpy as np

def F_grid(nth=4501, nh=10001):
    th = np.radians(np.linspace(1e-6, 45, nth))[:, None]
    s, c = np.sin(th), np.cos(th); H = s + c; k = 1 / (s * c)
    h = np.linspace(0, 1, nh)[None, :]                     # depth of the line y = 1 above the lowest vertex
    G = np.where(h <= s, k * h**2 / 2, np.where(h <= c, s / (2 * c) + (h - s) / c, 1 - k * np.clip(H - h, 0, None)**2 / 2))
    Wt = np.where(h <= s, k * h, np.where(h <= c, 1 / c, k * np.clip(H - h, 0, None)))
    y0 = 1 - h
    lo = np.maximum(0.5, y0 + s * c); hi = np.minimum(1, y0 + H - s * c)   # W >= 1 on depths [sc, H - sc]
    F = 0.25 * Wt + 1.5 * np.maximum(0, hi - lo) - G
    return np.degrees(th[:, 0]), y0[0], F

if __name__ == '__main__':
    ths, y0s, F = F_grid()
    i = np.unravel_index(F.argmin(), F.shape)
    print('min F = %.3e at theta = %.4f deg, y0 = %.4f' % (F[i], ths[i[0]], y0s[i[1]]))
    for d, y in ((45, 0.0), (45, 0.5)):
        a = np.abs(ths - d).argmin(); b = np.abs(y0s - y).argmin()
        print('  F(theta=%g, y0=%g) = %.2e   (tight family)' % (d, y, F[a, b]))
