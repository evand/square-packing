#!/usr/bin/env python3
"""Rung-2 diagnostic: evaluate a weighted closed cover of [0,4]^2 at the 'fan' poses around an
axis-parallel pose whose edges lie on grid lines, i.e. (cx, cy + delta, eps) for a range of
ratios delta/eps, and at the wall/corner poses.  Closed semantics: a point on the boundary
counts (tolerance 1e-12 in the inclusive direction).

usage: python3 search/zeromargin_fan.py runs/closed4_best.txt
"""
import sys, math
import numpy as np

def read_cert(path):
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    a = np.array(list(map(int, tok[5:5 + 3 * n])), dtype=float).reshape(n, 3)
    return sn / sd, a[:, 0] / D, a[:, 1] / D, a[:, 2] / W

def captured(X, Y, Wt, cx, cy, th, tol=1e-12):
    c, s = math.cos(th), math.sin(th)
    dx, dy = X - cx, Y - cy
    x = dx * c + dy * s; y = -dx * s + dy * c
    m = (np.abs(x) <= 0.5 + tol) & (np.abs(y) <= 0.5 + tol)
    return Wt[m].sum()

def main():
    s, X, Y, Wt = read_cert(sys.argv[1])
    print(f"container {s}, {len(X)} points, total {Wt.sum():.6f}")
    print("\n== axis-parallel poses with edges on grid lines (exact), captured weight ==")
    for cx in (0.5, 1.0, 1.5, 2.0):
        for cy in (0.5, 1.0, 1.5, 2.0):
            print(f"  centre ({cx},{cy}) angle 0: {captured(X, Y, Wt, cx, cy, 0.0):.4f}")
    print("\n== fan around (1.5, 1.5, 0): pose (1.5, 1.5+delta, eps), rows delta, cols eps ==")
    epss = [0.0, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2]
    print("  delta \\ eps  " + " ".join(f"{e:>8.0e}" for e in epss))
    for d in [0.0, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]:
        row = [captured(X, Y, Wt, 1.5, 1.5 + d, e) for e in epss]
        print(f"  {d:8.0e}     " + " ".join(f"{v:8.4f}" for v in row))
    print("\n== same fan around the wall square (1.5, 0.5 + delta, eps) [centre must be >= w/2] ==")
    print("  delta \\ eps  " + " ".join(f"{e:>8.0e}" for e in epss))
    for d in [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 3e-2]:
        row = []
        for e in epss:
            w2 = (math.cos(e) + math.sin(e)) / 2
            cy = max(0.5 + d, w2)
            row.append(captured(X, Y, Wt, 1.5, cy, e))
        print(f"  {d:8.0e}     " + " ".join(f"{v:8.4f}" for v in row))
    print("\n== corner square (w/2 + d, w/2 + d, eps) ==")
    for e in epss:
        w2 = (math.cos(e) + math.sin(e)) / 2
        print(f"  eps {e:.0e}: " + " ".join(f"d={d:.0e}:{captured(X, Y, Wt, w2 + d, w2 + d, e):.4f}" for d in (0.0, 1e-4, 1e-2)))
    # worst pose over a dense fan scan around every axis-parallel grid pose in the fundamental domain
    print("\n== minimum over fans (delta, eps in log grids, both signs of delta in x and y) around each grid pose ==")
    worst = (9, None)
    for cx in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
        for cy in (0.5, 1.0, 1.5, 2.0):
            for e in [10 ** k for k in np.linspace(-6, -1.3, 24)]:
                w2 = (math.cos(e) + math.sin(e)) / 2
                for d in [0] + [10 ** k for k in np.linspace(-6, -1.3, 24)]:
                    for sx in (-1, 0, 1):
                        for sy in (-1, 0, 1):
                            px, py = cx + sx * d, cy + sy * d
                            px = min(max(px, w2), 4 - w2); py = min(max(py, w2), 4 - w2)
                            v = captured(X, Y, Wt, px, py, e)
                            if v < worst[0]: worst = (v, (px, py, e))
    print(f"  worst: {worst[0]:.4f} at centre ({worst[1][0]:.6f}, {worst[1][1]:.6f}) angle {math.degrees(worst[1][2]):.4f} deg")

if __name__ == '__main__':
    main()
