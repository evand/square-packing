"""Are R5's strip-boundary packings rigid?  Five unit squares in [0, L] x [0, H] on the boundary curve, maximise
the common slack t (SLSQP, many random starts); report the distinct solutions with t >= -1e-7.  A continuum of
packings at the boundary means averaging them gives a fractional packing of value 5 and depth ~1, so no point
cover of total < 5 exists there.  Usage: OMP_NUM_THREADS=1 python3 strip_flex.py H [nstarts]"""
import sys, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from stage2 import gap
from strip_lp import Lb

def cons(z, L, H):
    t = z[-1]; qs = z[:-1].reshape(-1, 3)
    out = []
    for q in qs:
        k = (abs(np.cos(q[2])) + abs(np.sin(q[2]))) / 2
        out += [q[0] - k - t, L - k - q[0] - t, q[1] - k - t, H - k - q[1] - t]
    for i in range(5):
        for j in range(i + 1, 5):
            out.append(gap(qs[i], qs[j]) - t)
    return np.array(out)

def solve(args):
    L, H, seed = args
    rng = np.random.default_rng(seed)
    xs = np.sort(rng.uniform(0.5, L - 0.5, 5))
    z = np.concatenate([np.c_[xs, rng.uniform(0.5, H - 0.5, 5), rng.uniform(0, np.pi / 2, 5)].ravel(), [-0.5]])
    r = minimize(lambda z: -z[-1], z, method='SLSQP', constraints=[dict(type='ineq', fun=lambda z: cons(z, L, H))],
                 options=dict(maxiter=800, ftol=1e-14))
    c = cons(r.x, L, H)
    return (r.x[-1] if c.min() > -1e-8 else -9), r.x

if __name__ == '__main__':
    H = float(sys.argv[1]); n = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    L = Lb(H)
    with Pool(4) as pool:
        res = pool.map(solve, [(L, H, s) for s in range(n)])
    good = [x for t, x in res if t > -1e-7]
    print(f"H={H} L={L:.6f}: {len(good)}/{n} starts reach t >= -1e-7; best t = {max(t for t, _ in res):.3e}")
    # canonical form: sort by x, fold angle to [0, 90), mirror so that the first square is left-lower
    rows = []
    for x in good:
        q = x[:-1].reshape(5, 3).copy(); q[:, 2] = np.degrees(q[:, 2]) % 90
        q = q[np.argsort(q[:, 0])]
        rows.append(q)
    for q in rows[:12]:
        print("  " + "  ".join(f"({a:.3f},{b:.3f},{c:5.1f})" for a, b, c in q))
    if rows:
        A = np.array([q.ravel() for q in rows])
        print("spread (max - min) per coordinate over solutions:")
        print("  " + " ".join(f"{v:.3f}" for v in A.max(0) - A.min(0)))
