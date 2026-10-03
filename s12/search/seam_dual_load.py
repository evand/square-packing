"""Load function of the seam dual: l(p) = sum_S y_S #{j : p + j e1 in S} on the band cell.
Weak duality (derived): for valid mu, m_v <= int_{[0,1)x[0,w]} (Lam - l) dLeb, provided l <= Lam everywhere
in the band and l <= Lam - 1 on the seam {x in Z}.  Lam = -lambda."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import quadrant_lp as Q
w = float(sys.argv[1]); z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', 'seam_dual', 'w%s.npz' % sys.argv[1]))
y = z['ydual']; keep = y > 1e-9; P = z['P'][keep]; y = y[keep]; Lam = -float(z['lam'])
print('w', w, 'Lam', Lam, 'nposes', len(y), 'obj', z['obj'])
def load(px, py, tol=1e-9):
    l = np.zeros(px.shape)
    for (cx, cy, th), wt in zip(P, y):
        c, s = np.cos(th), np.sin(th)
        for j in range(-2, 3):
            dx = px + j - cx; dy = py - cy
            u = dx * c + dy * s; v = -dx * s + dy * c
            l += wt * (np.maximum(abs(u), abs(v)) <= 0.5 + tol)
    return l
n = 1000
xs = (np.arange(n) + 0.5) / n; ys = (np.arange(int(n * w)) + 0.5) / n
X, Y = np.meshgrid(xs, ys); L = load(X, Y)
print('generic load: max %.4f  min %.4f' % (L.max(), L.min()))
print('bound int(Lam - l) = %.4f   (LP says %.4f)' % (((Lam - L).mean() * w), -z['obj']))
# seam line, closed (tol) and points just off it
yy = np.linspace(0, w, 20001)
for x0 in (0.0, 1e-6, -1e-6, 1e-3):
    ls = load(np.full_like(yy, x0), yy)
    print('x=%g: max %.4f (need <= %.4f there if seam)' % (x0, ls.max(), Lam - 1))
# also horizontal lines: closed-edge maxima off-seam
for y0 in np.linspace(0, w, 21):
    xx = np.linspace(0, 1, 20001)
    m = load(xx, np.full_like(xx, y0)).max()
    if m > Lam + 1e-6: print('  excess on y=%.2f: %.4f' % (y0, m))
pass
