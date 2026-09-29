#!/usr/bin/env python3
"""qx2_tools.py -- float diagnostics for qx2_lp solutions (task quadrant-exact-w2).  HEURISTIC."""
import sys, os, math, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import quadrant_lp as Q
import qx2_lp as X


def load(path):
    d = np.load(path, allow_pickle=True)
    args = argparse.Namespace(**json.loads(str(d['args'])))
    m = X.build(args)
    return m, d['x'], args, d


def true_mass(m, x, P):
    """capture without the LP margin"""
    um = m.use_margin; m.use_margin = False
    try: return m.contrib(P, xval=x)
    finally: m.use_margin = um


def growth_scan(m, x, thetas=(1e-5, 1e-4, 1e-3, 1e-2, 3e-2), h=0.01):
    """(f(c,th) - 1)/th on a grid of centres (plus wall-resting and corner-resting poses) for small th, both signs."""
    ex, ey = Q._ext(m.G)
    xs = np.arange(0.5, ex + 1e-9, h); ys = np.arange(0.5, ey + 1e-9, h)
    X_, Y_ = np.meshgrid(xs, ys, indexing='ij'); X_ = X_.ravel(); Y_ = Y_.ravel()
    k = Y_ <= X_ + 1e-9; X_ = X_[k]; Y_ = Y_[k]
    out = []
    for th in thetas:
        for sg in (1, -1):
            t = (sg * th) % (math.pi / 2)
            hw = Q.halfwidth(np.array([t]))[0]
            P = np.c_[np.maximum(X_, hw + 1e-12), np.maximum(Y_, hw + 1e-12), np.full(len(X_), t)]
            v = true_mass(m, x, P)
            um = m.use_margin; m.use_margin = False; ar = m.lebesgue(P); m.use_margin = um
            out_u = np.clip(1 - ar, 0, 1)
            g = np.where(out_u > 1e-3, (v - 1) / (th * np.maximum(out_u, 1e-3)), np.inf)
            j = int(np.argmin(g))
            out.append((th, sg, float(g[j]), P[j].tolist(), int((g < 0).sum())))
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('sol'); ap.add_argument('--h', type=float, default=0.01)
    a = ap.parse_args()
    m, x, args, d = load(a.sol)
    print(a.sol, 'D', float(d['D']))
    for th, sg, g, p, nneg in growth_scan(m, x, h=a.h):
        print(f"th={sg*th:+.0e}: min (f-1)/(th*(1-areaU)) = {g:+.4f} at ({p[0]:.4f},{p[1]:.4f},{math.degrees(p[2]):.5f}deg)  #neg={nneg}")
