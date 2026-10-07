#!/usr/bin/env python3
"""clique_lp.py -- the packing LP over the poses of a support file, with clique cuts added lazily.

Value without cuts = the pure fractional packing mass over these poses (>= the certified L(t));
with cuts = the clique-strengthened value over the SAME poses (an upper bound on the true
strengthened value over these poses since the point constraints are a grid, and no statement about
other poses).  Shows how much the Helly gap is worth on the extremal measure.

    python3 search/clique_lp.py runs/dual_PA2_support.txt [--pitch 0.01] [--rounds 30]
"""
import sys, math, argparse, time
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
sys.path.insert(0, __import__('os').path.dirname(__file__))
import clique_check as CC


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('path'); ap.add_argument('--pitch', type=float, default=0.01)
    ap.add_argument('--rounds', type=int, default=30); ap.add_argument('--time', type=float, default=300); ap.add_argument('--per-round', type=int, default=30); ap.add_argument('--branch', action='store_true', help='input is a branch.py dual dump'); ap.add_argument('--kmass', type=float, default=None, help='with --branch: fix the mass of box poses to this (the leaf k)')
    a = ap.parse_args()
    if a.branch:
        t, poses = CC.read_branch(a.path)
        flags = [int(l.split()[4]) > 0 for l in open(a.path) if not l.startswith('#')]
    else:
        t, poses = CC.read_support(a.path); flags = None
    n = len(poses)
    # expanded images with their parent index
    img = []  # (cx, cy, th, h, parent)
    for i, (cx, cy, th, h, mu) in enumerate(poses):
        for g in range(8):
            x, y, ang = cx, cy, th
            if g & 1: x, y, ang = t - x, y, -ang
            if g & 2: x, y, ang = x, t - y, -ang
            if g & 4: x, y, ang = y, x, math.pi / 2 - ang
            img.append((x, y, ang % (math.pi / 2), h, i))
    IM = np.array([[p[0], p[1], p[2], p[3]] for p in img]); par = np.array([p[4] for p in img])
    # grid points in the fundamental domain 0 <= x <= y <= t/2
    g = np.arange(0, t / 2 + 1e-9, a.pitch); X, Y = np.meshgrid(g, g, indexing='ij'); pts = np.c_[X.ravel(), Y.ravel()]; pts = pts[pts[:, 0] <= pts[:, 1] + 1e-12]
    # coverage matrix: C[p, i] = (#images of pose i containing p) / 8
    rows = []; cols = []; vals = []
    ct = np.cos(IM[:, 2]); st = np.sin(IM[:, 2])
    for j in range(0, len(pts), 2000):
        Pp = pts[j:j + 2000]
        u = (Pp[:, 0:1] - IM[:, 0]) * ct + (Pp[:, 1:2] - IM[:, 1]) * st
        v = -(Pp[:, 0:1] - IM[:, 0]) * st + (Pp[:, 1:2] - IM[:, 1]) * ct
        ins = np.maximum(np.abs(u), np.abs(v)) <= IM[:, 3] + 1e-9
        r, c = np.nonzero(ins)
        rows.append(r + j); cols.append(par[c]); vals.append(np.full(len(r), 1 / 8))
    C = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(len(pts), n)); C.sum_duplicates()
    print(f"{a.path}: t={t} poses={n} (x8 images) grid points={len(pts)} pitch={a.pitch}")
    cuts = []  # each: array of parent indices with multiplicity (coefficients /8)
    for rd in range(a.rounds):
        A = C
        if cuts:
            K = sp.csr_matrix((np.concatenate([np.full(len(c), 1 / 8) for c in cuts]), (np.concatenate([np.full(len(c), k) for k, c in enumerate(cuts)]), np.concatenate(cuts))), shape=(len(cuts), n)); K.sum_duplicates()
            A = sp.vstack([C, K]).tocsr()
        Aeq = beq = None
        if flags is not None and a.kmass is not None:
            Aeq = np.array([[1.0 if f else 0.0 for f in flags]]); beq = [a.kmass]     # corner mass = k (leaf)
        res = linprog(c=-np.ones(n), A_ub=A, b_ub=np.ones(A.shape[0]), A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * n, method='highs')
        if res.status != 0: print('LP failed:', res.message); break
        y = res.x; mass = y.sum()
        # expanded measure for the clique finder
        P = [(IM[k, 0], IM[k, 1], IM[k, 2], IM[k, 3], y[par[k]] / 8) for k in range(len(IM))]
        keep = [k for k in range(len(P)) if P[k][4] > 1e-9]
        Pk = [P[k] for k in keep]; w = np.array([p[4] for p in Pk])
        adj = CC.overlap_graph(Pk)
        bw, bc, nodes, dt = CC.max_weight_clique(adj, w, time_limit=a.time)
        print(f"round {rd}: mass={mass:.6f} cuts={len(cuts)} max clique mass={bw:.6f} size={len(bc)} ({dt:.0f}s)", flush=True)
        if bw <= 1 + 1e-7: print("no violated clique: done"); break
        # one cut per clique; in parent coordinates every dihedral image of the clique is the same cut
        cuts.append(par[np.array([keep[i] for i in bc])])
        # more cuts this round: the best clique through each of the heaviest poses (restricted B&B)
        t1 = time.time(); added = 1
        for v in np.argsort(-w)[:a.per_round]:
            if time.time() - t1 > a.time: break
            nb = np.nonzero(adj[v])[0]
            if nb.size == 0: continue
            sub = adj[np.ix_(nb, nb)]
            bw2, bc2, _, _ = CC.max_weight_clique(sub, w[nb], time_limit=20)
            if bw2 + w[v] > 1 + 1e-7:
                cuts.append(par[np.array([keep[v]] + [keep[nb[i]] for i in bc2])]); added += 1
        print(f"         +{added} cuts ({time.time()-t1:.0f}s)", flush=True)
    print(f"final mass with {len(cuts)} clique cuts: {mass:.6f}  (pure LP over these poses: see round 0)")
    heavy = np.argsort(-y)[:12]
    for i in heavy:
        cx, cy, th, h, mu = poses[i]; print(f"   ({cx:.4f}, {cy:.4f}) theta={math.degrees(th):6.2f} y={y[i]:.4f}")


if __name__ == '__main__':
    main()
