"""10-point cover C: rows y = 1, b-1: 2/3 at x = c0, a - c0; 1/3 at x = 1, a - 1; 1/3 at (1, b/2), (a-1, b/2).
Total 14/3."""
import sys, math, numpy as np, ptcover as pc
from multiprocessing import Pool
def cover(a, b, c0, wr=2/3, wm=1/3, ym=None):
    ym = b / 2 if ym is None else ym
    P, W = [], []
    for X, wt in ((c0, wr), (a - c0, wr), (1.0, 1 - wr), (a - 1, 1 - wr)):
        for y in (1.0, b - 1): P.append((X, y)); W.append(wt)
    for X in (1.0, a - 1): P.append((X, ym)); W.append(wm)
    return np.array(P), np.array(W)
def line_b(a): return 3.0 if a <= 1 + math.sqrt(2) + 1e-12 else 4 + math.sqrt(2) - a
def mc(args):
    a, c0, th = args
    b = line_b(a); P, W = cover(a, b, c0)
    v, pose = pc.min_capture_angle(P, W, a, b, th)
    return v, pose
def check(a, c0, nth, pool):
    res = pool.map(mc, [(a, c0, th) for th in np.linspace(0, np.pi/2, nth, endpoint=False)], chunksize=50)
    return min(res, key=lambda t: t[0])
if __name__ == '__main__':
    c0s = [eval(x) for x in sys.argv[1].split(',')]
    As = [eval(x) for x in sys.argv[2].split(',')]
    nth = int(sys.argv[3])
    with Pool(8) as pool:
        for c0 in c0s:
            out = []
            for a in As:
                v, pose = check(a, c0, nth, pool); out.append(v)
                if v < 1 - 1e-9: print(f"  c0={c0:.5f} a={a:.5f} fails {v:.3f} at x={pose[0]:.4f} y={pose[1]:.4f} th={math.degrees(pose[2]):.2f}")
            print(f"c0={c0:.5f} min over a: {min(out):.4f}", flush=True)
