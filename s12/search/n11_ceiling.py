#!/usr/bin/env python3
"""Find the smallest container at which a D4-symmetrised fractional packing of mass 11 exists.

A `packing_dual.py`-style search is overkill once the shape of the optimum is known: at
`t = 3.85` the exact measure of `search/N11_ANATOMY.md` has only FOUR poses (the corner square,
a wall square, a tilted square and a near-axis square), so the whole object is 8 numbers.  This
script optimises those numbers directly:

  inner (exact-in-spirit, floats): given poses, the arrangement of their 8k dihedral images is
        enumerated (square corners + pairwise edge intersections), and the masses are chosen by
        the LP   max sum mu  s.t.  cov(v) <= 1 at every arrangement vertex v,  mu >= 0
        -- the same LP `dual_exact.py build` polishes with, but in floats and on one pose set;
  outer: Nelder-Mead / random restarts over the pose parameters, maximising the LP mass.

Nothing here is rigorous.  Its output is a candidate pose set; the statement is made by
`dual_exact.py build --src <candidate>` + `dual_exact.py check --full`, which redo the whole
thing in exact rational arithmetic.

    python3 search/n11_ceiling.py --t 191/50 --restarts 40 --out runs/n11_cand_3.82.txt
    python3 search/n11_ceiling.py --bisect 3.80 3.85 --steps 8
"""
import argparse
import math
import os
import sys
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import linprog, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(os.path.dirname(HERE), 'runs')
TOL = 1e-9


def images(cx, cy, th, t):
    return [(cx, cy, th), (t - cx, cy, -th), (cx, t - cy, -th), (t - cx, t - cy, th),
            (cy, cx, -th), (t - cy, cx, th), (cy, t - cx, th), (t - cy, t - cx, -th)]


def corners_of(cx, cy, th):
    c, s = math.cos(th), math.sin(th)
    return [(cx + c * a - s * b, cy + s * a + c * b) for a, b in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]


def admissible(cx, cy, th, t):
    w = (abs(math.cos(th)) + abs(math.sin(th))) / 2
    return w - 1e-12 <= cx <= t - w + 1e-12 and w - 1e-12 <= cy <= t - w + 1e-12


def arrangement(imgs):
    """vertices of the arrangement: square corners + pairwise edge intersections (floats)"""
    V = []
    segs = []
    for cx, cy, th in imgs:
        c = corners_of(cx, cy, th)
        V += c
        segs += [(c[i], c[(i + 1) % 4]) for i in range(4)]
    P1 = np.array([s[0] for s in segs]); D1 = np.array([s[1] for s in segs]) - P1
    n = len(segs)
    for i in range(n):
        p, d = P1[i], D1[i]
        det = d[0] * D1[:, 1] - d[1] * D1[:, 0]
        ok = np.abs(det) > 1e-9
        w = P1 - p
        tt = np.where(ok, (w[:, 0] * D1[:, 1] - w[:, 1] * D1[:, 0]) / np.where(ok, det, 1), -1)
        uu = np.where(ok, (w[:, 0] * d[1] - w[:, 1] * d[0]) / np.where(ok, det, 1), -1)
        sel = ok & (tt >= -1e-12) & (tt <= 1 + 1e-12) & (uu >= -1e-12) & (uu <= 1 + 1e-12)
        sel[:i + 1] = False
        if sel.any():
            V += list(p + tt[sel, None] * d)
    return np.array(V)


def inside_mask(V, imgs):
    """A[v, k] = 1 if vertex v is in image k (closed, with tolerance)"""
    A = np.zeros((len(V), len(imgs)))
    for k, (cx, cy, th) in enumerate(imgs):
        c, s = math.cos(th), math.sin(th)
        dx = V[:, 0] - cx; dy = V[:, 1] - cy
        A[:, k] = (np.maximum(np.abs(dx * c + dy * s), np.abs(-dx * s + dy * c)) <= .5 + 1e-8)
    return A


def lp_mass(poses, t):
    """max total mass of a D4-symmetrised measure on `poses` with coverage <= 1 everywhere"""
    for p in poses:
        if not admissible(p[0], p[1], p[2], t):
            return -1.0, None
    imgs = []; own = []
    for k, (cx, cy, th) in enumerate(poses):
        for im in images(cx, cy, th, t):
            imgs.append(im); own.append(k)
    V = arrangement(imgs)
    V = V[(V[:, 0] >= -1e-9) & (V[:, 1] >= -1e-9) & (V[:, 0] <= t + 1e-9) & (V[:, 1] <= t + 1e-9)]
    A = inside_mask(V, imgs)
    own = np.array(own)
    B = np.zeros((len(V), len(poses)))
    for k in range(len(poses)):
        B[:, k] = A[:, own == k].sum(axis=1) / 8.0
    keep = B.max(axis=1) > 0
    B = B[keep]
    res = linprog(-np.ones(len(poses)), A_ub=B, b_ub=np.ones(len(B)), bounds=(0, None), method='highs')
    if res.status != 0:
        return -1.0, None
    return -res.fun, np.maximum(res.x, 0.0)


FAMILY = ['corner (1/2,1/2,0)', 'wall (1/2,y,0)', 'tilted (x,y,th)', 'near-axis (x,y,th)']


def unpack(z, t):
    """z = [y2, x3, y3, th3, x4, y4, th4]  ->  the four poses"""
    return [(0.5, 0.5, 0.0), (0.5, z[0], 0.0), (z[1], z[2], z[3]), (z[4], z[5], z[6])]


def objective(z, t):
    m, _ = lp_mass(unpack(z, t), t)
    return -m


def search(t, restarts, seed, z0=None, log=print):
    rng = np.random.default_rng(seed)
    best = (-1.0, None)
    starts = []
    if z0 is not None:
        starts.append(np.array(z0))
    for _ in range(restarts):
        starts.append(np.array([rng.uniform(1.0, t / 2),
                                rng.uniform(0.9, t / 2), rng.uniform(t / 2, t - 0.9), rng.uniform(0.2, 0.8),
                                rng.uniform(0.9, t / 2), rng.uniform(0.9, t - 0.9), rng.uniform(0.0, 0.35)]))
    for i, z in enumerate(starts):
        r = minimize(objective, z, args=(t,), method='Nelder-Mead',
                     options=dict(maxiter=1500, xatol=1e-7, fatol=1e-9))
        if -r.fun > best[0]:
            best = (-r.fun, r.x)
            log(f"  restart {i}: mass {-r.fun:.9f}  z={np.round(r.x, 6).tolist()}")
    return best


def write_candidate(path, t, poses, mu):
    with open(path, 'w') as f:
        f.write(f"# t={t} candidate pose set (masses are re-chosen by dual_exact.py build)\n")
        f.write("# D4-symmetrised measure: pose cx cy theta_deg mu\n")
        for (cx, cy, th), m in zip(poses, mu):
            if m <= 1e-9:
                continue
            f.write(f"pose {cx:.13f} {cy:.13f} {math.degrees(th):.11f} {m:.12f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--t', default=None, help='single container side (exact fraction)')
    ap.add_argument('--bisect', nargs=2, type=float, default=None, metavar=('LO', 'HI'))
    ap.add_argument('--steps', type=int, default=8)
    ap.add_argument('--restarts', type=int, default=40)
    ap.add_argument('--seed', type=int, default=20260913)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    z_seed = [1.5002, 1.30193, 2.449725, math.radians(33.69997), 1.565, 1.625, math.radians(5.0)]

    def run(tf, z0):
        m, z = search(tf, a.restarts, a.seed, z0=z0)
        poses = unpack(z, tf)
        mass, mu = lp_mass(poses, tf)
        print(f"t = {tf:.6f}:  best LP mass = {mass:.9f}   {'>= 11' if mass >= 11 - 1e-9 else '< 11'}")
        for (cx, cy, th), mm in zip(poses, mu):
            print(f"    pose {cx:.6f} {cy:.6f} {math.degrees(th):9.5f}  mu {mm:.9f}")
        return mass, z, poses, mu

    if a.bisect:
        lo, hi = a.bisect
        z0 = z_seed
        # rescale the seed to the container being tried
        for k in range(a.steps):
            mid = (lo + hi) / 2
            zz = list(z0)
            mass, z, poses, mu = run(mid, zz)
            if mass >= 11 - 1e-9:
                hi = mid; z0 = list(z)
            else:
                lo = mid
            print(f"  bracket: [{lo:.6f}, {hi:.6f}]")
        print(f"FINAL bracket for the mass-11 threshold of this family: [{lo:.6f}, {hi:.6f}]")
    else:
        tf = float(Fr(a.t))
        mass, z, poses, mu = run(tf, z_seed)
        if a.out:
            write_candidate(a.out, tf, poses, mu)
            print(f"wrote {a.out}")


if __name__ == '__main__':
    main()
