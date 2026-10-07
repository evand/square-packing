#!/usr/bin/env python3
"""How many squares of a packing can be centred in a level-2 region?  (`alpha(R)`)

For a region R (an axis-parallel box of centres) inside [0,t]^2 this maximises the
*disjointness margin* of a pair (or triple) of admissible poses with centres in R:

    margin(P1, P2) = max over the four edge normals n of ( |<c2-c1, n>| - (1 + w(delta))/2 )

where delta is the angle difference and w(a) = |cos a| + |sin a|.  Two closed unit squares have
disjoint interiors iff margin >= 0 (separating-axis theorem for rectangles: along a normal n of
square i the half-widths are 1/2 and w(delta)/2).  So

    max margin < 0  =>  alpha(R) = 1   (no two squares of a packing are centred in R)

and a witness with margin >= 0 exhibits two.  Admissibility: cx, cy in [w/2, t - w/2].

Search: a coarse grid over (cx1, cy1, th1, cx2, cy2, th2) followed by Nelder-Mead polish from the
best starts; heuristic, one-sided (it can only find pairs, never prove their absence -- but a
clearly negative maximum with a large search is strong evidence, and the note records the
analytic argument that settles the cases that matter).

Usage:  python3 search/level2_capacity.py --t 3.98 --box 1 1.99 0 1        # a wall slot
        python3 search/level2_capacity.py --t 3.98 --box 1 2.98 0 1 --k 3  # a full wall strip
"""
import argparse, itertools, math

import numpy as np
from scipy.optimize import minimize


def wid(a):
    return abs(math.cos(a)) + abs(math.sin(a))


def margin_pair(p, q):
    """max over the 4 normals of |<d,n>| - (1 + w(delta))/2 ; >= 0 iff interiors are disjoint"""
    d = np.array([q[0] - p[0], q[1] - p[1]])
    dl = q[2] - p[2]
    half = (1.0 + wid(dl)) / 2.0
    best = -1e9
    for th in (p[2], q[2]):
        for n in (np.array([math.cos(th), math.sin(th)]), np.array([-math.sin(th), math.cos(th)])):
            best = max(best, abs(float(d @ n)) - half)
    return best


def clip(p, t, box):
    """project a pose into the admissible part of the box"""
    cx, cy, th = p
    w2 = wid(th) / 2
    lox = max(box[0], w2); hix = min(box[1], t - w2)
    loy = max(box[2], w2); hiy = min(box[3], t - w2)
    if lox > hix or loy > hiy:
        return None
    return (min(max(cx, lox), hix), min(max(cy, loy), hiy), th)


def boxes_for(boxes, k):
    return [boxes[min(i, len(boxes) - 1)] for i in range(k)]


def score(x, t, boxes, k):
    bs = boxes_for(boxes, k)
    poses = []
    for i in range(k):
        p = clip((x[3 * i], x[3 * i + 1], x[3 * i + 2]), t, bs[i])
        if p is None:
            return 1e3
        poses.append(p)
    return -min(margin_pair(poses[i], poses[j]) for i, j in itertools.combinations(range(k), 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--t', type=float, required=True)
    ap.add_argument('--box', type=float, nargs=4, action='append', required=True,
                    metavar=('X0', 'X1', 'Y0', 'Y1'),
                    help='repeat for a mixed question ("2 in this box and 1 in that one"): '
                         'the i-th square is confined to the i-th box, the last box is reused')
    ap.add_argument('--k', type=int, default=2, help='how many squares to try to place')
    ap.add_argument('--grid', type=int, default=7, help='grid points per centre coordinate')
    ap.add_argument('--angles', type=int, default=13, help='angle samples in [0, 90)')
    ap.add_argument('--starts', type=int, default=400, help='local polishes from the best grid starts')
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    t, boxes, k = a.t, a.box, a.k
    bs = boxes_for(boxes, k)

    ths = np.linspace(0, math.pi / 2, a.angles, endpoint=False)
    single = []
    for b in bs:
        xs = np.linspace(b[0], b[1], a.grid)
        ys = np.linspace(b[2], b[3], a.grid)
        pool = [clip((x, y, th), t, b) for x in xs for y in ys for th in ths]
        single.append([p for p in pool if p is not None])
    rng = np.random.default_rng(a.seed)

    # coarse: random k-subsets of the grid, keep the best starts
    cand = []
    for _ in range(200000):
        ps = [single[i][rng.integers(0, len(single[i]))] for i in range(k)]
        v = -min(margin_pair(ps[i], ps[j]) for i, j in itertools.combinations(range(k), 2))
        cand.append((v, np.array([c for p in ps for c in p])))
    cand.sort(key=lambda z: z[0])

    best = (1e9, None)
    for v, x0 in cand[:a.starts]:
        r = minimize(score, x0, args=(t, boxes, k), method='Nelder-Mead',
                     options={'maxiter': 4000, 'xatol': 1e-9, 'fatol': 1e-12})
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    v, x = best
    poses = [clip((x[3 * i], x[3 * i + 1], x[3 * i + 2]), t, bs[i]) for i in range(k)]
    print(f't = {t}  k = {k}  boxes = ' + '; '.join(f'[{b[0]}, {b[1]}] x [{b[2]}, {b[3]}]' for b in bs))
    print(f'  best min pairwise margin = {-v:+.6f}   ({"FEASIBLE: k squares fit" if -v >= 0 else "no k squares (heuristic)"})')
    for p in poses:
        print(f'    pose ({p[0]:.6f}, {p[1]:.6f}, {math.degrees(p[2]) % 90:.3f} deg)')
    for i, j in itertools.combinations(range(k), 2):
        print(f'    margin({i},{j}) = {margin_pair(poses[i], poses[j]):+.6f}')


if __name__ == '__main__':
    main()
