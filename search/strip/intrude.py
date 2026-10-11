"""Given a strip config (4 squares in [0,W]x[0,H]), how far left can one more unit square reach
(min x of its vertices) while disjoint from the 4 and inside [0,H]^2?  Reports depth c = W - minx."""
import numpy as np, sys
from scipy.optimize import minimize
from rectpack import penalty, gap, half
f, H = sys.argv[1], float(sys.argv[2])
v = np.load(f); n = len(v)//3
x, y, t = v[:n], v[n:2*n], v[2*n:]
W = max(x + half(t))
# normalise: reflect so the top flat square is flush right (user's orientation is irrelevant for depth: try both reflections)
rng = np.random.default_rng(0)
for refl in (False, True):
    xs = W - x if refl else x
    ts = -t if refl else t
    best = (-1, None)
    for s in range(400):
        q = np.array([rng.uniform(W-0.6, W+0.6), rng.uniform(0.5, H-0.5), rng.uniform(0, np.pi/2)])
        def obj(q, mu):
            vv = np.r_[np.r_[xs, q[0]], np.r_[y, q[1]], np.r_[ts, q[2]]]
            return (q[0] - half(q[2])) + mu*penalty(vv, n+1, H, H)
        for mu in (1e2, 1e4, 1e6, 1e8):
            q = minimize(obj, q, args=(mu,), method='L-BFGS-B').x
        vv = np.r_[np.r_[xs, q[0]], np.r_[y, q[1]], np.r_[ts, q[2]]]
        if gap(vv, n+1, H, H) > -1e-6:
            c = W - (q[0] - half(q[2]))
            if c > best[0]: best = (c, q)
    c, q = best
    print(f"{f} refl={refl} W={W:.4f} a={W-1:.4f} max intrusion c={c:.4f} at y={q[1]:.3f} theta={np.degrees(q[2])%90:.2f}  (a+b+... needed: c > {2*1.2071-3+W-1+(4-H):.3f} with h*, {2*1.19-3+W-1+(4-H):.3f} with h=1.19)")
