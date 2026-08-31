#!/usr/bin/env python3
"""Separate anchor cliques on a dual (fractional packing) dump of search/branch.py.

    python3 search/anchorsep.py runs/branch_TAG_dual.txt [--pitch 0.01] [--top 200] [--fracs ...]

The dual of the cover LP is a fractional packing `y` on the rows.  A clique column `K` is worth
adding exactly when `ybar(K) > 1` (its reduced cost is `cost*(1 - ybar(K))`), which is the same
thing as: the packing violates the clique constraint `mu(K) <= 1`.  The interesting case is
`ybar(K) > 1 >= ybar(P_p)`, i.e. a violation the coverage row at `p` does NOT already see: that is
the non-Helly mass the family is supposed to catch (search/CLIQUE_CONTINUUM.md).

Scans a grid of wall points `p` (Lemma 1: nowhere else can a clique beat the coverage row at `p`),
keeps the heaviest, and tries the Lemma-2 anchors for several `eps` fractions.  Floats: this is a
heuristic that only chooses candidates; the certificate is exact and the verifier checks it.
"""
import sys, os, argparse, math
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anchorclique as AC


def load_dual(path):
    rows = []; ys = []; s = None
    for l in open(path):
        if l.startswith('#'):
            for t in l.split():
                if t.startswith('s='):
                    a, b = t[2:].split('/'); s = float(a) / float(b)
            continue
        q = l.split()
        if len(q) < 6: continue
        rows.append((float(q[0]), float(q[1]), float(q[2]), float(q[3]))); ys.append(float(q[5]))
    return np.array(rows), np.array(ys), s


def cover_grid(R, y, G):
    """coverage sum_r y_r [g in square r] for every point g of G (n,2)"""
    ct = np.cos(R[:, 2]); st = np.sin(R[:, 2]); h = R[:, 3]
    u0 = R[:, 0] * ct + R[:, 1] * st; u1 = -R[:, 0] * st + R[:, 1] * ct
    out = np.zeros(len(G))
    for i in range(0, len(R), 512):
        c, s_, hh, a0, a1, yy = ct[i:i+512], st[i:i+512], h[i:i+512], u0[i:i+512], u1[i:i+512], y[i:i+512]
        q0 = G[:, 0:1] * c + G[:, 1:2] * s_; q1 = -G[:, 0:1] * s_ + G[:, 1:2] * c
        out += (np.maximum(np.abs(q0 - a0), np.abs(q1 - a1)) <= hh + 1e-12) @ yy
    return out


def separate(R, y, s, D, pitch=0.01, top=300, fracs=(0.95, 0.8, 0.6, 0.4, 0.2),
             slack=0.10, pad=0.002, band=1.0, sym=True):
    """Scan the wall bands of the container for anchor cliques the packing `y` violates.
    Returns [(ybar(K), ybar(P_p), (X, Y, wall, eps, rho), clique), ...] sorted by ybar(K).
    `ybar` is the D4 average when `sym` (an orbit column of cost = number of images), which is what
    an orbit column's reduced cost uses; a symmetric model then needs only one wall band scanned."""
    if len(R) == 0: return []
    n = int(round(s / pitch)) + 1
    ax = np.linspace(0, s, n)
    G = []
    for w in range(1 if sym else 4):
        bnd = ax[(ax > 0.02) & (ax < band)]
        for v in bnd:
            for u in ax:
                G.append((v, u) if w == 0 else (s - v, u) if w == 1 else (u, v) if w == 2 else (u, s - v))
    G = np.array(G)
    gs = [lambda p: p, lambda p: np.c_[s - p[:, 0], p[:, 1]], lambda p: np.c_[p[:, 0], s - p[:, 1]],
          lambda p: np.c_[s - p[:, 0], s - p[:, 1]], lambda p: np.c_[p[:, 1], p[:, 0]],
          lambda p: np.c_[s - p[:, 1], p[:, 0]], lambda p: np.c_[p[:, 1], s - p[:, 0]],
          lambda p: np.c_[s - p[:, 1], s - p[:, 0]]]
    cov = (sum(cover_grid(R, y, g(G)) for g in gs) / 8.0) if sym else cover_grid(R, y, G)
    rot = AC._rot(R); sfr = F(int(round(s * D)), D)
    out = []; seen = set()
    for i in np.argsort(-cov)[:top]:
        X, Y = int(round(G[i][0] * D)), int(round(G[i][1] * D))
        for frac in fracs:
            for (wall, en, rn, d) in AC.cand_params(sfr, D, X, Y, frac=frac, slack=slack, pad=pad):
                if (X, Y, wall, en, rn) in seen: continue
                seen.add((X, Y, wall, en, rn))
                cl = AC.kpa(sfr, D, X, Y, wall, en, rn)
                if cl is None: continue
                imgs = AC.images(sfr, cl) if sym else [cl]
                pcl = ((('P', F(X, D), F(Y, D)),), ((0, ()),))
                pim = AC.images(sfr, pcl) if sym else [pcl]
                m = AC.coeff(imgs, R, rot); mp = AC.coeff(pim, R, rot)
                out.append((float(m @ y) / len(imgs), float(mp @ y) / len(pim), (X, Y, wall, en, rn), cl))
    out.sort(key=lambda t: -t[0])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dual'); ap.add_argument('--pitch', type=float, default=0.01)
    ap.add_argument('--top', type=int, default=300); ap.add_argument('--show', type=int, default=25)
    ap.add_argument('--fracs', default='0.95,0.8,0.6,0.4,0.2')
    ap.add_argument('--slack', type=float, default=0.10); ap.add_argument('--pad', type=float, default=0.002)
    ap.add_argument('--band', type=float, default=1.0, help='only wall distances below this')
    ap.add_argument('--out', default=None, help='write the best cliques as a branch clique checkpoint')
    ap.add_argument('--D', type=int, default=784000)
    ap.add_argument('--asym', action='store_true', help='the model is not D4-symmetric: price single cliques, scan all four walls')
    a = ap.parse_args()
    R, y, s = load_dual(a.dual)
    print(f"{len(R)} rows, mass {y.sum():.6f}, s = {s}")
    D = a.D
    best = separate(R, y, s, D, pitch=a.pitch, top=a.top, fracs=[float(t) for t in a.fracs.split(',')],
                    slack=a.slack, pad=a.pad, band=a.band, sym=not a.asym)
    print(f"{'ybar(K)':>9} {'ybar(P_p)':>10} {'extra':>9}  {'p':>17} wall {'eps':>7} {'rho':>7}")
    for (mk, mp, (X, Y, wall, en, rn), cl) in best[:a.show]:
        print(f"{mk:9.5f} {mp:10.5f} {mk-mp:9.5f}  ({X/D:7.4f},{Y/D:7.4f}) {wall:4d} {en/D:7.4f} {rn/D:7.4f}")
    if a.out:
        ncl = 0
        with open(a.out, 'w') as f:
            f.write(f"{int(round(s*D))} {D} {0 if a.asym else 1}\n")
            for (mk, mp, (X, Y, wall, en, rn), cl) in best:
                if mk <= 1.0: break
                f.write(f"{X} {Y} {wall} {en} {rn}\n"); ncl += 1
        print(f"wrote {ncl} cliques with ybar(K) > 1 to {a.out}")


if __name__ == '__main__':
    main()
