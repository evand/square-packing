#!/usr/bin/env python3
"""unavoid13 calibration [heuristic]: the fractional cover value at container side m, estimated
by a float cutting-plane LP -- columns = the (dominance-reduced) arrangement vertices of a family
file, rows = grid poses added when violated -- and its HONEST cost (LP total / min captured
weight over a finer scan).  Also reports the LP relaxation of the hitting-set IP over the family
alone (a float lower bound on nu_f(m), since its dual is a packing measure on F that is feasible
at every arrangement vertex, hence everywhere).

    python3 search/unavoid13_calib.py FAMILY.txt [--dth 2] [--pitch 0.04] [--threads 2]
"""
import sys, os, math, time, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L


def captured(P, W, cx, cy, th):
    c, s = math.cos(th), math.sin(th)
    dx = P[None, :, 0] - cx[:, None]; dy = P[None, :, 1] - cy[:, None]
    xp = dx * c + dy * s; yp = -dx * s + dy * c
    inside = (np.abs(xp) <= 0.5 + 1e-12) & (np.abs(yp) <= 0.5 + 1e-12)
    return inside @ W, inside


def grid_poses(m, dth, pitch):
    out = []
    for th in np.arange(0.0, math.pi / 2, dth):
        w = L.w_of_theta(th); lo, hi = w / 2, m - w / 2
        g = np.arange(lo, hi + 1e-12, pitch)
        if g[-1] < hi - 1e-12: g = np.append(g, hi)
        CX, CY = np.meshgrid(g, g, indexing='ij')
        out.append((th, CX.ravel(), CY.ravel()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family')
    ap.add_argument('--dth', type=float, default=2.0)
    ap.add_argument('--pitch', type=float, default=0.04)
    ap.add_argument('--fine-dth', type=float, default=0.25)
    ap.add_argument('--fine-pitch', type=float, default=0.005)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--rounds', type=int, default=40)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    m, fam = L.read_family(a.family)
    mf = float(m)
    sq, V, B, rep = L.build_candidates(m, fam)
    K = len(rep)
    # LP relaxation over F alone
    r = L.solve_hitting_set(B, threads=a.threads)
    print(f"hitting-set LP relaxation over F ({len(fam)} squares, {K} candidates): {r['lp']:.6f}; IP h(F) = {r['obj']}")
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix, vstack
    rows = [csr_matrix(B.T.astype(float))]          # F's own rows
    G = grid_poses(mf, math.radians(a.dth), a.pitch)
    P = rep
    t0 = time.time()
    val = None
    for it in range(a.rounds):
        A = vstack(rows)
        res = linprog(np.ones(K), A_ub=-A, b_ub=-np.ones(A.shape[0]), bounds=(0, None), method='highs',
                      options=dict(disp=False))
        W = res.x; val = res.fun
        # scan the grid for violated poses
        worst = []; nviol = 0; minw = 9
        for (th, cx, cy) in G:
            cw, inside = captured(P, W, cx, cy, th)
            minw = min(minw, cw.min())
            bad = np.nonzero(cw < 1 - 1e-7)[0]
            nviol += len(bad)
            if len(bad):
                order = bad[np.argsort(cw[bad])][:60]
                worst.append(inside[order])
        print(f"  round {it}: LP = {val:.6f}, rows {A.shape[0]}, grid min {minw:.6f}, violated {nviol}, {time.time()-t0:.0f}s", flush=True)
        if nviol == 0: break
        rows.append(csr_matrix(np.concatenate(worst, 0).astype(float)))
    # honest cost on a finer scan
    Gf = grid_poses(mf, math.radians(a.fine_dth), a.fine_pitch)
    minw = 9; arg = None
    for (th, cx, cy) in Gf:
        cw, _ = captured(P, W, cx, cy, th)
        k = int(cw.argmin())
        if cw[k] < minw: minw = float(cw[k]); arg = (float(cx[k]), float(cy[k]), math.degrees(th))
    honest = val / minw
    print(f"dense-grid cover LP value {val:.6f}; fine-scan min captured {minw:.6f} at {arg}; honest cost {honest:.6f}")
    print(f"support: {int((W > 1e-9).sum())} points with positive weight")
    out = a.out or os.path.splitext(a.family)[0] + '_calib.json'
    with open(out, 'w') as f:
        json.dump(dict(family=a.family, m=str(m), lp_relax_F=r['lp'], h_F=r['obj'], grid_lp=val,
                       fine_min=minw, fine_arg=arg, honest=honest, dth=a.dth, pitch=a.pitch,
                       fine_dth=a.fine_dth, fine_pitch=a.fine_pitch, ncand=K), f, indent=1)
    print("written", out)


if __name__ == '__main__':
    main()
