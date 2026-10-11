"""Worst-case seam over strips with slack: strip of 4 in [0,w] x [0,H], with its top square pinned to the
top-right corner contact (touching right wall and top) and its lowest square touching the left wall and bottom
(the A/B/F pattern); a fifth square Q disjoint from all.  Minimise Q's leftmost x jointly over strip and Q."""
import numpy as np, sys
from scipy.optimize import minimize
from rectpack import penalty, gap, half
from seam_common import load_strip
import os
H = float(os.environ.get("HH", "3.99")); TH0 = np.radians(float(os.environ.get("TH0", "0"))); sg = sys.argv[1].replace("_", "|"); w = float(sys.argv[2]); N = int(sys.argv[3]); YLO, YHI = float(sys.argv[4]), float(sys.argv[5])
x0, y0, t0, W0 = load_strip(sg, 3.99); y0 = y0 * (H / 3.99)
x0 = x0 * (w / W0)   # stretch the representative to width w as a seed
rng = np.random.default_rng(1)
def split(z): return z[:5], z[5:10], z[10:15]
def contacts(z):
    x, y, t = split(z); hs = half(t[:4]); kt = np.argmax(y[:4]); kb = np.argmin(y[:4])
    om = np.argsort(y[:4])[1:3]; am = np.abs(((t[om] + np.pi/4) % (np.pi/2)) - np.pi/4)
    return (np.sum(np.maximum(0, TH0 - am)**2) + np.maximum(0, YLO - z[9])**2 + np.maximum(0, z[9] - YHI)**2 + (w - x[kt] - hs[kt])**2 + (H - y[kt] - hs[kt])**2 + (x[kb] - hs[kb])**2 + (y[kb] - hs[kb])**2)
def strip_box(z):
    x, y, t = split(z); hs = half(t[:4])
    return np.sum(np.maximum(0, x[:4] + hs - w)**2)
best = (9, None)
for s in range(N):
    q = np.array([rng.uniform(1.0, w + 0.8), rng.uniform(0.5, H - 0.5), rng.uniform(0, np.pi/2)])
    z = np.r_[np.r_[x0 + rng.normal(0, .02, 4), q[0]], np.r_[y0 + rng.normal(0, .02, 4), q[1]], np.r_[t0 + rng.normal(0, .02, 4), q[2]]]
    for mu in (1e2, 1e4, 1e6, 1e8):
        f = lambda z: (z[4] - half(z[14])) + mu*(penalty(z, 5, H, H) + strip_box(z) + contacts(z))
        z = minimize(f, z, method='L-BFGS-B', options={'maxiter': 20000}).x
    if gap(z, 5, H, H) > -1e-7 and strip_box(z) < 1e-12 and contacts(z) < 1e-10:
        m = z[4] - half(z[14])
        if m < best[0]: best = (m, z)
m, z = best
if z is None: print(f"{sg} H={H} tilt>={os.environ.get('TH0','0')} w={w}: no valid run"); sys.exit()
x, y, t = split(z); a = np.degrees(t) % 90; a = np.minimum(a, 90 - a)
print(f"{sg} H={H} tilt>={os.environ.get('TH0','0')} w={w}: worst leftmost reach {m:.4f} ({'SeamOK x0=1 holds, margin %.3f' % (m-1) if m >= 1 else 'x0=1 FAILS by %.3f' % (1-m)}); strip tilts {np.round(np.sort(a[:4]),1)}, Q tilt {a[4]:.1f} at y={y[4]:.2f}", flush=True)
np.save(f"seamband_H{H}_th{os.environ.get("TH0","0")}_{sys.argv[1]}_{w}_{YLO}.npy", z)
