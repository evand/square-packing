"""Optimistic seam with the real top strip square: 8 squares in [0,H]^2 avoiding the column [0,1]x[0,H]
and the strip's top square (oriented to touch the right wall and the top).  Float search."""
import numpy as np, sys, json, glob
from scipy.optimize import minimize
from rectpack import penalty, gap, half
H = 3.99; sigwant = sys.argv[1]; starts = int(sys.argv[2])
rows = [r for f in glob.glob('census4*_H3.99_s*.json') for r in json.load(open(f)) if r['sig'] == sigwant]
r = min(rows, key=lambda r: r['W']); z = np.array(r['z']); W = r['W']
x, y, t = z[:4], z[4:8], z[8:]
best = None
for fx in (0, 1):
    for fy in (0, 1):
        xx = W - x if fx else x; yy = H - y if fy else y; tt = -t if fx ^ fy else t
        k = np.argmax(yy)
        if W - xx[k] - half(tt[k]) < 2e-3 and H - yy[k] - half(tt[k]) < 2e-3: best = (xx[k], yy[k], tt[k]); break
    if best: break
assert best, "no orientation with top square touching right and top"
tx, ty, tt0 = best
print(f"{sigwant}: strip W={W:.4f}, top square centre ({tx:.3f},{ty:.3f}) tilt {np.degrees(tt0)%90:.1f}")
def rect_pen(x, y, t, cx, cy, hx, hy):
    dx, dy = x - cx, y - cy; best = np.full(len(x), np.inf)
    for ph in (np.zeros(len(x)), np.full(len(x), np.pi/2), t, t + np.pi/2):
        u = np.abs(np.cos(ph)*dx + np.sin(ph)*dy)
        r = half(t - ph) + hx*np.abs(np.cos(ph)) + hy*np.abs(np.sin(ph))
        best = np.minimum(best, r - u)
    return best
def P(q):
    fx, fy, ft = q[:8], q[8:16], q[16:]
    v = np.r_[np.r_[fx, tx], np.r_[fy, ty], np.r_[ft, tt0]]
    return penalty(v, 9, H, H) + np.sum(np.maximum(0, rect_pen(fx, fy, ft, 0.5, H/2, 0.5, H/2))**2)
def ok(q):
    fx, fy, ft = q[:8], q[8:16], q[16:]
    v = np.r_[np.r_[fx, tx], np.r_[fy, ty], np.r_[ft, tt0]]
    return gap(v, 9, H, H) > -1e-10 and rect_pen(fx, fy, ft, 0.5, H/2, 0.5, H/2).max() < 1e-10
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
        print(f"{sigwant}: 8 FIT (start {s})"); np.save(f'seam_fit_{sigwant.replace("|","_")}.npy', q); sys.exit()
print(f"{sigwant}: no fit in {starts} starts")
