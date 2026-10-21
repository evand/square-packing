"""Check candidate parameterised covers for R5's strip part on the boundary curve (ptcover.min_capture).
Usage: python3 strip_cand.py NAME nth H1 H2 ...   (L = curve value unless given as H,L)"""
import sys, numpy as np
from multiprocessing import Pool
import ptcover as pc
from strip_lp import Lb

def mirror(pts, L, H):
    """D2 orbit of quadrant points (x <= L/2, y <= H/2) with weights; duplicates merged."""
    out = {}
    for (x, y), w in pts:
        for X in {round(x, 12), round(L - x, 12)}:
            for Y in {round(y, 12), round(H - y, 12)}:
                out[(X, Y)] = w
    P = np.array(list(out.keys())); W = np.array(list(out.values()))
    return P, W

def cover(name, L, H):
    if name == 'mid4':          # (i, H/2), i = 1..4, weight 1
        return np.array([(i, H / 2) for i in range(1, 5)], float), np.ones(4)
    if name == 'two_rows':      # (i, y0), (i, H - y0), weight 1/2
        y0 = float(name_args.get('y0', 0.6))
        return mirror([((1, y0), .5), ((2, y0), .5)], L, H)
    raise ValueError(name)

name_args = {}

if __name__ == '__main__':
    name = sys.argv[1]; nth = int(sys.argv[2])
    for kv in [a for a in sys.argv[3:] if '=' in a]:
        k, v = kv.split('='); name_args[k] = v
    Hs = [a for a in sys.argv[3:] if '=' not in a]
    with Pool(4) as pool:
        for s in Hs:
            if ',' in s: H, L = map(float, s.split(','))
            else: H = float(s); L = Lb(H)
            P, W = cover(name, L, H)
            res = pc.min_capture(P, W, L, H, nth=nth, pool=pool)
            v, pose = min(res, key=lambda t: t[0])
            print(f"{name} H={H:.5f} L={L:.5f} total={W.sum():.4f} min_capture={v:.6f} at {tuple(np.round(pose, 4))}",
                  flush=True)
