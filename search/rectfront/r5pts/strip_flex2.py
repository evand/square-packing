"""Rigidity of the zigzag at the strip boundary: re-optimise (max slack t) from perturbed zigzags; if the
t ~ 0 solutions vary beyond solver noise, the boundary packing is flexible.  Usage: python3 strip_flex2.py H"""
import sys, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from strip_flex import cons
from strip_lp import Lb
R2 = np.sqrt(2); r = R2 / 2

def zig(H):
    L = Lb(H); xb = 2 - H + R2
    return L, np.array([[.5, .5, 0], [xb, H - r, np.pi / 4], [L / 2, r, np.pi / 4], [L - xb, H - r, np.pi / 4],
                        [L - .5, .5, 0]])

def solve(args):
    H, seed, eps = args
    L, q = zig(H)
    rng = np.random.default_rng(seed)
    z = np.concatenate([(q + eps * rng.standard_normal(q.shape)).ravel(), [-0.1]])
    res = minimize(lambda z: -z[-1], z, method='SLSQP', constraints=[dict(type='ineq', fun=lambda z: cons(z, L, H))],
                   options=dict(maxiter=1000, ftol=1e-15))
    c = cons(res.x, L, H)
    return (res.x[-1] if c.min() > -1e-9 else -9), res.x[:-1].reshape(5, 3)

if __name__ == '__main__':
    H = float(sys.argv[1])
    L, q0 = zig(H)
    print(f"zigzag at H={H}: slack {cons(np.r_[q0.ravel(), 0], L, H).min():.2e}")
    for eps in (0.01, 0.05):
        with Pool(4) as pool:
            res = pool.map(solve, [(H, s, eps) for s in range(64)])
        good = [q for t, q in res if t > -1e-8]
        ts = [t for t, _ in res]
        print(f"eps={eps}: {len(good)}/64 with t > -1e-8 (max t {max(ts):.2e})")
        if good:
            A = np.array([np.c_[g[:, :2], np.degrees(g[:, 2]) % 90] for g in good])
            d = A - np.c_[q0[:, :2], np.degrees(q0[:, 2]) % 90][None]
            print("  max |deviation from zigzag| per square (x, y, deg):")
            for i in range(5): print("   ", np.round(np.abs(d[:, i, :]).max(0), 5))
            for g in good[:4]:
                print("   ", " ".join(f"({a:.4f},{b:.4f},{np.degrees(c) % 90:6.2f})" for a, b, c in g))
