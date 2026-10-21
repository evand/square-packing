"""Task 4a: squares avoiding the rectangular 4-point set at a = 1+sqrt2."""
import numpy as np
r2 = np.sqrt(2); a = 1 + r2
for b in [2.95, 2.99, 2.999]:
    pts = [(1, 1), (a-1, 1), (1, b-1), (a-1, b-1)]
    rows = []
    for th in np.linspace(0, np.pi/2, 1801, endpoint=False):
        c, s = np.cos(th), np.sin(th); k = (c+s)/2
        xs = np.linspace(k, a-k, 400); ys = np.linspace(k, b-k, 400)
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        m = np.full(X.shape, np.inf)
        for px, py in pts:
            dx, dy = px-X, py-Y
            m = np.minimum(m, np.maximum(abs(dx*c+dy*s), abs(-dx*s+dy*c)))
        e = m - .5
        if (e > 1e-9).any():
            sel = e > 1e-9
            rows.append((np.degrees(th), e.max(), X[sel].min(), X[sel].max(), Y[sel].min(), Y[sel].max()))
    rows = np.array(rows)
    print(f"b={b}: escaping tilts {rows[:,0].min():.2f}..{rows[:,0].max():.2f} deg, max esc {rows[:,1].max():.4e} at {rows[rows[:,1].argmax(),0]:.2f} deg")
    print(f"   centre x in [{rows[:,2].min():.4f},{rows[:,3].max():.4f}] (a={a:.4f}), y in [{rows[:,4].min():.4f},{rows[:,5].max():.4f}] (b/2={b/2})")
    left = rows  # report per-angle x-range for a few angles
    for t in [rows[0,0], 40, 45, 50, rows[-1,0]]:
        r = rows[np.argmin(abs(rows[:,0]-t))]
        print(f"   {r[0]:.2f} deg: x {r[2]:.4f}..{r[3]:.4f}, y {r[4]:.4f}..{r[5]:.4f}, esc {r[1]:.2e}")
