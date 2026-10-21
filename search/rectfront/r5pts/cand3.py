"""Family: rows y = 1, b-1: weight be at x = c0, a - c0; 1 - be at x = 1, a - 1; interior ga at (xi, yi) + D2 mirrors.
Total 4 + 4 ga."""
import sys, math, numpy as np, ptcover as pc
from multiprocessing import Pool
def cover(a, b, c0, be, xi, yi, ga):
    P = []; W = []
    for x, wt in ((c0, be), (1.0, 1 - be)):
        for X in (x, a - x):
            for y in (1.0, b - 1):
                P.append((X, y)); W.append(wt)
    for X in (xi, a - xi):
        for y in (yi, b - yi):
            P.append((X, y)); W.append(ga)
    return np.array(P), np.array(W)
def line_b(a): return 3.0 if a <= 1 + math.sqrt(2) + 1e-12 else 4 + math.sqrt(2) - a
def check(a, b, prm, pool, nth=1800):
    P, W = cover(a, b, *prm)
    res = pc.min_capture(P, W, a, b, nth=nth, pool=pool)
    v = min(res, key=lambda t: t[0]); return v
if __name__ == '__main__':
    prm = eval(sys.argv[1]); As = [float(s) for s in sys.argv[2:]]
    with Pool(4) as pool:
        for a in As:
            b = line_b(a); v = check(a, b, prm, pool)
            print(f"prm={prm} a={a:.4f} b={b:.4f} total={4+4*prm[-1]:.4f} min_cap={v[0]:.4f} at {np.round(v[1][:2],3)} th={np.degrees(v[1][2]):.2f}", flush=True)
