import sys, numpy as np, ptcover as pc
from multiprocessing import Pool
def cover(a, b, c0=0.8, al=0.6, be=0.4, ga=0.1, yi=1.2):
    P = []; W = []
    for x, wt in ((c0, be), (1.0, al)):
        for X in (x, a - x):
            for y in (1.0, b - 1):
                P.append((X, y)); W.append(wt)
    for X in (1.0, a - 1):
        for y in (yi, b - yi):
            P.append((X, y)); W.append(ga)
    return np.array(P), np.array(W)
if __name__ == '__main__':
    kw = eval(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].startswith('{') else {}
    pairs = [tuple(map(float, s.split(','))) for s in sys.argv[2 if kw or (len(sys.argv)>1 and sys.argv[1]=='{}') else 1:]]
    with Pool(4) as pool:
        for a, b in pairs:
            P, W = cover(a, b, **kw)
            res = pc.min_capture(P, W, a, b, nth=1800, pool=pool)
            v = min(res, key=lambda t: t[0])
            print(f"a={a:.4f} b={b:.4f} total={W.sum():.4f} min_capture={v[0]:.4f} at {np.round(v[1][:2],4)} th={np.degrees(v[1][2]):.2f}", flush=True)
