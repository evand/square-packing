#!/usr/bin/env python3
"""clique_check.py -- does a certified fractional packing use non-Helly cliques?

The cover LP is the fractional *point*-clique cover of the pose overlap graph: the poses through a
point pairwise overlap, so sum y <= 1 over them.  Rotated squares are not a Helly family, so there
are cliques (pairwise interior-overlapping sets of poses) with no common point, and each is a
constraint  sum_{S in K} y_S <= 1  the point-cover LP does not have.  This script takes a fractional
packing (a support file of packing_dual.py, or a dual dump of branch.py), expands the D4 symmetry,
builds the interior-overlap graph and finds the maximum-mass clique.  Mass > 1 on some clique means
the pure method's ceiling is the ceiling of the Helly sub-relaxation; a certified measure has mass
<= 1 on every point-clique, so such a clique is necessarily non-Helly.

    python3 search/clique_check.py runs/dual_PA2_support.txt
    python3 search/clique_check.py runs/branch_TAG_dual.txt --branch
"""
import sys, math, argparse, time
import numpy as np
from scipy.optimize import linprog


def read_support(path):
    t = None; poses = []
    for line in open(path):
        if line.startswith('# t='): t = float(line.split('=')[1].split()[0])
        q = line.split()
        if q and q[0] == 'pose':
            cx, cy, thd, mu = map(float, q[1:5]); poses.append((cx, cy, math.radians(thd), 0.5, mu))
    return t, poses


def read_branch(path):
    t = None; poses = []
    for line in open(path):
        if line.startswith('#'):
            for tok in line.split():
                if tok.startswith('s='):
                    a, b = tok[2:].split('/'); t = int(a) / int(b)
            continue
        cx, cy, th, h, fl, y = line.split(); poses.append((float(cx), float(cy), float(th), float(h), float(y)))
    return t, poses


def expand_d4(poses, t):
    """mass mu/8 on each of the 8 dihedral images; coincident images accumulate"""
    acc = {}
    for cx, cy, th, h, mu in poses:
        for g in range(8):
            x, y, a = cx, cy, th
            if g & 1: x, y, a = t - x, y, -a
            if g & 2: x, y, a = x, t - y, -a
            if g & 4: x, y, a = y, x, math.pi / 2 - a
            a = a % (math.pi / 2)
            if a > math.pi / 2 - 1e-9: a = 0.0
            key = (round(x, 7), round(y, 7), round(a, 7), round(h, 6))
            acc[key] = acc.get(key, 0.0) + mu / 8
    return [(k[0], k[1], k[2], k[3], m) for k, m in acc.items()]


def overlap_graph(P, eps=1e-9):
    """adjacency: interiors intersect (separating axis test on the 4 edge normals, strict)"""
    n = len(P); C = np.array([[p[0], p[1]] for p in P]); TH = np.array([p[2] for p in P]); H = np.array([p[3] for p in P])
    adj = np.zeros((n, n), dtype=bool)
    # quick reject: centre distance >= h_i*sqrt2 + h_j*sqrt2
    D = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=2)
    cand = D < (H[:, None] + H[None, :]) * math.sqrt(2) + 1e-9
    np.fill_diagonal(cand, False)
    for i in range(n):
        js = np.nonzero(cand[i])[0]
        if js.size == 0: continue
        ok = np.ones(js.size, dtype=bool)
        for src in (i, None):
            angs = np.array([TH[i], TH[i] + math.pi / 2]) if src is not None else None
            for k in range(2):
                if src is not None:
                    ax = np.array([math.cos(angs[k]), math.sin(angs[k])])[None, :].repeat(js.size, 0)
                else:
                    a = TH[js] + k * math.pi / 2; ax = np.c_[np.cos(a), np.sin(a)]
                # projection half-extent of a square with half side h at angle th onto unit axis u:
                #   h*(|cos(th-phi)| + |sin(th-phi)|), phi = angle of u
                phi = np.arctan2(ax[:, 1], ax[:, 0])
                ei = H[i] * (np.abs(np.cos(TH[i] - phi)) + np.abs(np.sin(TH[i] - phi)))
                ej = H[js] * (np.abs(np.cos(TH[js] - phi)) + np.abs(np.sin(TH[js] - phi)))
                ci = C[i, 0] * ax[:, 0] + C[i, 1] * ax[:, 1]; cj = C[js, 0] * ax[:, 0] + C[js, 1] * ax[:, 1]
                ok &= (np.abs(ci - cj) < ei + ej - eps)
        adj[i, js[ok]] = True
    adj = adj & adj.T
    return adj


def max_weight_clique(adj, w, time_limit=600):
    """branch and bound with a greedy-colouring weight bound (Kumlander-style)"""
    n = len(w); nbr = [np.nonzero(adj[i])[0] for i in range(n)]
    nbrset = [set(v.tolist()) for v in nbr]
    best = [0.0, []]; t0 = time.time(); nodes = [0]
    def colour_bound(cands):
        # greedy colouring; bound = sum over classes of max weight
        classes = []; cw = []
        for v in sorted(cands, key=lambda v: -w[v]):
            for ci, cl in enumerate(classes):
                if not (nbrset[v] & cl):
                    cl.add(v); break
            else:
                classes.append({v}); cw.append(w[v])
        return sum(cw)
    def expand(cands, cur, curw):
        nodes[0] += 1
        if time.time() - t0 > time_limit: return
        if not cands:
            if curw > best[0]: best[0] = curw; best[1] = list(cur)
            return
        if curw + colour_bound(cands) <= best[0] + 1e-12: return
        cl = sorted(cands, key=lambda v: -w[v])
        while cl:
            v = cl.pop(0)
            rest = [u for u in cl if u in nbrset[v]]
            expand(rest, cur + [v], curw + w[v])
            if curw + colour_bound(cl) <= best[0] + 1e-12: return
    order = sorted(range(n), key=lambda v: -w[v])
    expand(order, [], 0.0)
    return best[0], best[1], nodes[0], time.time() - t0


def common_point(P, idx):
    """is there a point in all the closed squares idx?  LP feasibility on 4 half-planes per square"""
    A = []; b = []
    for i in idx:
        cx, cy, th, h, _ = P[i]
        for k in range(4):
            a = th + k * math.pi / 2; u = (math.cos(a), math.sin(a))
            A.append([u[0], u[1]]); b.append(cx * u[0] + cy * u[1] + h)
    res = linprog(c=[0, 0], A_ub=A, b_ub=b, bounds=[(None, None)] * 2, method='highs')
    return res.status == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path'); ap.add_argument('--branch', action='store_true'); ap.add_argument('--time', type=float, default=900)
    ap.add_argument('--top', type=int, default=5, help='also report the best clique through each of the top-k heaviest poses')
    a = ap.parse_args()
    t, poses = (read_branch if a.branch else read_support)(a.path)
    P = expand_d4(poses, t)
    w = np.array([p[4] for p in P])
    print(f"{a.path}: t={t} poses={len(poses)} -> {len(P)} after D4, mass={w.sum():.6f}")
    adj = overlap_graph(P)
    deg = adj.sum(1)
    print(f"overlap graph: edges={int(adj.sum() // 2)} mean degree={deg.mean():.1f} max degree={deg.max()}")
    bw, bc, nodes, dt = max_weight_clique(adj, w, time_limit=a.time)
    helly = common_point(P, bc)
    print(f"max-mass clique: mass={bw:.6f} size={len(bc)} nodes={nodes} time={dt:.0f}s  common point: {helly}")
    for i in sorted(bc, key=lambda i: -w[i]):
        cx, cy, th, h, m = P[i]; print(f"   ({cx:.4f}, {cy:.4f}) theta={math.degrees(th):6.2f} h={h:.4f} mass={m:.5f}")
    # non-Helly triples with large mass: brute force over edges of heavy vertices
    heavy = [i for i in np.argsort(-w)[:200]]
    best3 = (0.0, None)
    for ii, i in enumerate(heavy):
        for j in heavy[ii + 1:]:
            if not adj[i, j]: continue
            for k in heavy:
                if k <= j or not (adj[i, k] and adj[j, k]): continue
                m3 = w[i] + w[j] + w[k]
                if m3 > best3[0] and not common_point(P, [i, j, k]): best3 = (m3, (i, j, k))
    if best3[1]:
        print(f"heaviest non-Helly triple among the 200 heaviest poses: mass={best3[0]:.6f}")
        for i in best3[1]:
            cx, cy, th, h, m = P[i]; print(f"   ({cx:.4f}, {cy:.4f}) theta={math.degrees(th):6.2f} mass={m:.5f}")
    else:
        print("no non-Helly triple among the 200 heaviest poses")


if __name__ == '__main__':
    main()
