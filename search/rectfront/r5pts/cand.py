import sys, numpy as np, ptcover as pc
from multiprocessing import Pool
def cover(a, b, xr=(0.6, 0.9), xi=0.9, yi=1.2, wi=1/6):
    P = []; W = []
    for x in list(xr) + [1.0]:
        for X in (x, a - x):
            for y in (1.0, b - 1):
                P.append((X, y)); W.append(1/3)
    for X in (xi, a - xi):
        for y in (yi, b - yi):
            P.append((X, y)); W.append(wi)
    return np.array(P), np.array(W)
if __name__ == '__main__':
    a = float(sys.argv[1])
    with Pool(4) as pool:
        for b in map(float, sys.argv[2:]):
            P, W = cover(a, b)
            res = pc.min_capture(P, W, a, b, nth=1800, pool=pool)
            v = min(res, key=lambda t: t[0])
            print(f"a={a} b={b} total={W.sum():.4f} min_capture={v[0]:.4f} at x={v[1][0]:.4f} y={v[1][1]:.4f} th={np.degrees(v[1][2]):.2f}" if v[1] else v, flush=True)
