#!/usr/bin/env python3
"""plus_pack.py -- heuristic packing of m unit squares in the plus region P_t = [0,t]^2 minus its
four corner unit squares (interior-disjoint from [0,1]^2 and images).  Penalty minimisation with
random restarts; bisection on t for the smallest container in which m squares were found to fit.
The fractional packing number of P_t is 8 at t = 3.98 and 3.99 (corner leaf k = 4 of BRANCH.md);
this measures the integral side.

    python3 search/plus_pack.py --m 8 --t 3.98 --restarts 400      # best penetration at fixed t
    python3 search/plus_pack.py --m 7 --bisect 3.5 4.0 --restarts 150
"""
import sys, math, argparse, time
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 0)


def corners(cx, cy, th, h=0.5):
    c, s = math.cos(th), math.sin(th)
    return np.array([[cx + h * (c * sx - s * sy), cy + h * (s * sx + c * sy)] for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1))])


def penetration(A, B):
    """SAT penetration depth of two convex polygons (0 if interiors disjoint)"""
    best = float('inf')
    for poly in (A, B):
        for k in range(len(poly)):
            e = poly[(k + 1) % len(poly)] - poly[k]; n = np.array([-e[1], e[0]]); n /= np.linalg.norm(n)
            pa = A @ n; pb = B @ n
            ov = min(pa.max(), pb.max()) - max(pa.min(), pb.min())
            if ov <= 0: return 0.0
            best = min(best, ov)
    return best


def loss(v, m, t, fixed):
    v = v.reshape(m, 3); polys = [corners(*v[i]) for i in range(m)]; L = 0.0
    for i in range(m):
        p = polys[i]
        L += (np.maximum(0, -p[:, 0]) ** 2).sum() + (np.maximum(0, p[:, 0] - t) ** 2).sum()
        L += (np.maximum(0, -p[:, 1]) ** 2).sum() + (np.maximum(0, p[:, 1] - t) ** 2).sum()
        for F in fixed: L += penetration(p, F) ** 2
        for j in range(i + 1, m): L += penetration(p, polys[j]) ** 2
    return L


def solve(m, t, restarts, x0=None, iters=400):
    fixed = [corners(0.5, 0.5, 0), corners(t - 0.5, 0.5, 0), corners(0.5, t - 0.5, 0), corners(t - 0.5, t - 0.5, 0)]
    best = (float('inf'), None)
    for r in range(restarts):
        if x0 is not None and r == 0: v = x0.copy()
        else:
            v = np.c_[rng.uniform(0.5, t - 0.5, m), rng.uniform(0.5, t - 0.5, m), rng.uniform(0, math.pi / 2, m)].ravel()
            if x0 is not None and r < restarts // 3: v = x0 + rng.normal(0, 0.05, x0.shape)
        res = minimize(loss, v, args=(m, t, fixed), method='L-BFGS-B', options=dict(maxiter=iters))
        if res.fun < best[0]: best = (res.fun, res.x)
        if best[0] < 1e-12: break
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--m', type=int, default=8); ap.add_argument('--t', type=float, default=3.98)
    ap.add_argument('--restarts', type=int, default=200); ap.add_argument('--bisect', type=float, nargs=2, default=None)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--steps', type=int, default=10)
    a = ap.parse_args(); t0 = time.time()
    if a.bisect is None:
        L, x = solve(a.m, a.t, a.restarts)
        print(f"m={a.m} t={a.t}: best loss {L:.3e} ({'FITS' if L < 1e-12 else 'does not fit (heuristic)'}) in {time.time()-t0:.0f}s")
        for i in range(a.m): print(f"   ({x[3*i]:.4f}, {x[3*i+1]:.4f}) theta={math.degrees(x[3*i+2]) % 90:.2f}")
        return
    lo, hi = a.bisect; xhi = None
    L, xhi = solve(a.m, hi, a.restarts)
    print(f"m={a.m} t={hi}: loss {L:.3e} {'FITS' if L < 1e-12 else 'NO'} ({time.time()-t0:.0f}s)", flush=True)
    if L >= 1e-12: print("does not fit at the upper end; stop"); return
    for _ in range(a.steps):
        mid = (lo + hi) / 2
        # warm start: scale the last feasible configuration into the smaller container
        x0 = xhi.copy().reshape(-1, 3); x0[:, :2] *= mid / hi; x0 = x0.ravel()
        L, x = solve(a.m, mid, a.restarts, x0=x0)
        print(f"m={a.m} t={mid:.5f}: loss {L:.3e} {'FITS' if L < 1e-12 else 'NO'} ({time.time()-t0:.0f}s)", flush=True)
        if L < 1e-12: hi, xhi = mid, x
        else: lo = mid
    print(f"m={a.m}: fits at t={hi:.5f}, not found at t={lo:.5f}")
    for i in range(a.m): print(f"   ({xhi[3*i]:.4f}, {xhi[3*i+1]:.4f}) theta={math.degrees(xhi[3*i+2]) % 90:.2f}")


if __name__ == '__main__':
    main()
