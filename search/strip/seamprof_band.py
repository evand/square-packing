"""Margin of the top-square block: shift the strip's top square right by delta (giving back room).
8 squares must avoid column [0,1]x[0,H] and the shifted top square.  Finds whether 8 fit."""
import numpy as np, sys, json, glob
from scipy.optimize import minimize
from rectpack import penalty, gap, half
from seam_common import load_top, rect_pen
import os
H = float(os.environ.get('HH', '3.99')); sigwant, delta, starts, X0, M = sys.argv[1].replace('_', '|'), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
tx, ty, tt0, W = load_top(sigwant, 3.99); ty = ty + (H - 3.99); tx -= delta  # move the top square LEFT, giving room to the 8
OFF = 10.0
def parts(q):
    return q[:8], q[8:16], q[16:]
HB = 1.17
def band(fy):
    ys = np.sort(fy)
    return max(0, HB - ys[2])**2 + max(0, ys[-3] - (H - HB))**2
def P(q):
    fx, fy, ft = parts(q)
    p = band(fy) + penalty(q, 8, H, H) + np.sum(np.maximum(0, rect_pen(fx, fy, ft, X0/2, H/2, X0/2, H/2))**2) + np.sum(np.maximum(0, rect_pen(fx, fy, ft, (1+M)/2, H/2, (1+M)/2, H/2-1))**2)
    v = np.r_[np.r_[fx, tx]+OFF, np.r_[fy, ty]+OFF, np.r_[ft, tt0]]
    return p + penalty(v, 9, 100, 100)
def ok(q):
    fx, fy, ft = parts(q)
    v = np.r_[np.r_[fx, tx]+OFF, np.r_[fy, ty]+OFF, np.r_[ft, tt0]]
    return band(fy) < 1e-14 and gap(q, 8, H, H) > -1e-10 and gap(v, 9, 100, 100) > -1e-10 and rect_pen(fx, fy, ft, X0/2, H/2, X0/2, H/2).max() < 1e-10 and rect_pen(fx, fy, ft, (1+M)/2, H/2, (1+M)/2, H/2-1).max() < 1e-10
base = np.load('cfg_n8_H3.99.npy'); bx, by, bt = base[:8], base[8:16], base[16:]
rng = np.random.default_rng(0)
for s in range(starts):
    if rng.random() < 0.6:
        fl = rng.random() < 0.5
        q = np.r_[bx + (H - 2.8625) + rng.normal(0, .12, 8), (H - by if fl else by) + rng.normal(0, .12, 8), (-bt if fl else bt) + rng.normal(0, .12, 8)]
    else:
        q = np.r_[rng.uniform(1.5, H-.5, 8), rng.uniform(.5, H-.5, 8), rng.uniform(0, np.pi/2, 8)]
    q = minimize(P, q, method='L-BFGS-B', options={'ftol':1e-24,'gtol':1e-16,'maxiter':20000}).x
    if ok(q):
        np.save(f"spb_fit_{H}_{X0}.npy", q); print(f"[band] H={H} x0={X0}: 8 FIT (start {s})", flush=True); sys.exit()
print(f"[band<=2 of 8 within 1.17 of top/bottom] H={H} {sigwant} x0={X0} m={M} delta={delta}: no fit in {starts} starts", flush=True)
