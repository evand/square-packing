#!/usr/bin/env python3
"""can a SOUND, positive-volume clique rule credit the fan cliques?  (task `helly-allmeet`)

Doc: `search/ALLMEET.md`.  Companion to `search/honestcost.py` / `search/HONEST.md`, which showed
that every clique row carrying dual on the best pure `t = 4` instances is non-Helly with an empty
core (so the verifier's sound `clique_of_cores` rule credits it nothing) while the pricer's rule
("credit `z_K` to any pose that closed-meets every member") is not sound (the credited set is not
a clique).  This module asks the question in between:

    for a clique row `K`, is there a positive-volume BOX of poses `B_K` such that `K u B_K` is
    pairwise closed-meeting?  If so, "credit `z_K` iff `S in K u B_K`" IS a sound clique rule with
    positive-volume membership, and we can measure the honest cost of the dual under it.

Semantics (`search/ZEROMARGIN.md` §1, unchanged): `t = 4`; closed unit squares; closed containment
and closed meeting (a shared boundary point counts); a pose is `(cx, cy, u)` with `u = tan(th/2)`
and `th in [0, 90)` (a unit square is quarter-turn invariant); admissible iff the closed square
lies in `[0,4]^2`.  A BOX is `[x0,x1] x [y0,y1] x [u0,u1]`; positive volume means all three sides
strictly positive.

Two closed squares meet iff none of their 4 distinct edge normals separates them
(`leaf_ceiling.sq_meets_sq`).  Writing the SIGNED margin on an axis as
`1/2 + w/2 - |<c_j - c_i, e>|` with `w = |cos(th_j - th_i)| + |sin(th_j - th_i)|`, the squares meet
iff all four margins are `>= 0`.

WHAT IS EXACT HERE
------------------
`box_meets_square` and `box_pairwise_meets` and `box_admissible` are exact (`Fraction` /
integer): every `vol(B_K) > 0` reported by `boxes` is backed by running them on every member.
The search that PROPOSES a box (the max-min-margin LP sweep, the growth bisection) is float, and
so is the honest-cost measurement of §`honest`, exactly as in `HONEST.md`.

Commands
--------
    boxes  TAG   per dual-carrying clique row: the max-min meet margin over all admissible poses
                 (float LP sweep), the inscribed box grown from the best pose and CERTIFIED
                 EXACTLY, and the row's grazing anatomy   -> runs/am_TAG_boxes.json
    honest TAG   the honest cost of the dual under the sound rule `credit z_K iff S in K u B_K`
                 (honestcost's scan machinery with that rule)  -> runs/am_TAG_honest.json
    report TAG   the tables of `search/ALLMEET.md` from the two json files
    helly        is there a ROBUST non-Helly triple at `t = 4`?  (three unit squares in `[0,4]^2`,
                 pairwise meeting with strictly positive margin, empty triple intersection, each
                 thickened to a positive-volume pose box)  -> runs/am_helly.json
"""
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import leaf_ceiling as lc          # noqa: E402
import honestcost as hc            # noqa: E402

RUNS = lc.RUNS
T = Fr(4)
HALF = Fr(1, 2)


# ====================================================================== exact pose bookkeeping
def read_exact_poses(tag, base=RUNS, t=T, r=Fr(1), delta=1e-6):
    """the exact poses `(p, q, cx, cy)` of a `cl_TAG_poses.txt` checkpoint, in the index order
    `cliquelever.Poses` assigns them (the order the clique member lists use).

    This replays `Poses.add`'s boundary snap and admissibility clamp rather than building a whole
    `Lever` (which, on `E2Pg`, means re-reading a 75 MB polygon file).  `boxes` asserts that the
    float shadow of the result equals the `pf` array `honestcost.dump_dual` wrote, so the
    correspondence is checked, not assumed."""
    tt, sym, poses = lc.read_measure(os.path.join(base, f'cl_{tag}_poses.txt'))
    assert tt == t and sym == 1, (tt, sym)

    def bsnap(x):
        for b in (r, t / 2, t - r):
            if abs(float(x - b)) <= delta:
                return b
        return x

    P, key = [], {}
    for (p, q, cx, cy, m) in poses:
        cx, cy = bsnap(Fr(cx)), bsnap(Fr(cy))
        a, b, rr = lc.cos_sin(p, q)
        w2 = lc.half_width(a, b, rr)
        cx = min(max(cx, w2), t - w2)
        cy = min(max(cy, w2), t - w2)
        k = (p, q, cx, cy)
        if k not in key:
            key[k] = len(P)
            P.append(k)
    return P


def fold_pose(pose, t=T):
    """the SAME closed square, re-parametrised so that `u = tan(th/2)` lies in `[0, 1)`.

    A unit square is invariant under a quarter turn, so `th -> th +- 90` leaves the square alone
    as a point set while moving `u` by `u -> (1+u)/(1-u)` / `u -> (u-1)/(1+u)`.  Both checkpoints
    carry poses with `th < 0` (`84` of `SUPEPGF`'s 956, `1,998` of `E2Pg`'s 14,611, down to
    `u = -0.414`), and against those `w_lo_exact`'s branch analysis -- which needs
    `|th - th_k| < 90` -- is still SOUND (`cos D + sin D <= |cos D| + |sin D|` for every `D`) but
    badly loose: on E2Pg row 31 it bounds `w >= 0.427` where the truth is `1.348`, and the box
    test then refuses boxes that plainly meet the member.  Folding first removes that entirely.

    Asserted, per pose, by comparing the exact corner sets of the two squares."""
    p, q, cx, cy = pose
    if q < 0:
        p, q = -p, -q
    while p < 0:                       # th -> th + 90
        p, q = p + q, q - p
    while p >= q:                      # th -> th - 90
        p, q = p - q, q + p
    g = math.gcd(abs(p), abs(q))
    if g > 1:
        p, q = p // g, q // g
    return (p, q, cx, cy)


def same_square(a, b, t=T):
    """exact: do two pose records describe the same closed square (as a point set)?"""
    sa, sb = lc.make_square(a[2], a[3], a[0], a[1], t), lc.make_square(b[2], b[3], b[0], b[1], t)
    ca = sorted((Fr(X, sa[7]), Fr(Y, sa[7])) for X, Y in sa[6])
    cb = sorted((Fr(X, sb[7]), Fr(Y, sb[7])) for X, Y in sb[6])
    return ca == cb


def pose_cs(p, q):
    """exact `(cos th, sin th)` of the pose with `u = tan(th/2) = p/q`"""
    a, b, r = lc.cos_sin(p, q)
    return Fr(a, r), Fr(b, r)


def cs_of_u(u):
    """exact `(cos th, sin th)` from `u = tan(th/2)` (a Fraction)"""
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d


# ====================================================================== the exact box tests
def _max_quad(a0, a1, a2, u0, u1):
    """exact max of `a0 + a1 u + a2 u^2` over `[u0, u1]` (Fractions); `zeromargin._max_quad`."""
    v = a0 + a1 * u0 + a2 * u0 * u0
    v1 = a0 + a1 * u1 + a2 * u1 * u1
    if v1 > v:
        v = v1
    if a2 < 0:
        uv = -a1 / (2 * a2)
        if u0 < uv < u1:
            vv = a0 + a1 * uv + a2 * uv * uv
            if vv > v:
                v = vv
    return v


def w_lo_exact(uk, u0, u1):
    """an exact rational LOWER bound on `w(th, th_k) = |cos(th - th_k)| + |sin(th - th_k)|` over
    the poses of a `u`-interval `[u0, u1]`, where the member's angle has `tan(th_k/2) = uk`.

    Both `th` and `th_k` lie in `[0, 90)`, so `D = th - th_k in (-90, 90)` and `cos D > 0`:
    `w = cos D + sin D` when `D >= 0` and `cos D - sin D` when `D <= 0`.  Each of those is
    `sqrt2 cos(D -+ 45)` with the argument in `[-45, 45]`, where `cos` is CONCAVE, so its minimum
    over a `D`-interval is attained at an endpoint -- and `D` is monotone in `u`.  Hence:

    * `u0 >= uk`: `w >= min(w(u0), w(u1))` with `w(u) = cos D + sin D`, exact;
    * `u1 <= uk`: the same with `cos D - sin D`;
    * the interval straddles `uk`: `w >= 1` (`|cos| + |sin| >= 1` always, with equality exactly at
      `D = 0`, which the interval contains -- so this is the exact minimum there too).

    A lower bound is what the meet test needs: `w` enters the margin as `+w/2`, so replacing it by
    a lower bound only makes "the squares meet" harder to certify.  Sound, never optimistic."""
    if u0 <= uk <= u1:
        return Fr(1)
    ck, sk = cs_of_u(uk)
    best = None
    for u in (u0, u1):
        c, s = cs_of_u(u)
        cD = c * ck + s * sk
        sD = s * ck - c * sk
        v = cD + sD if u0 >= uk else cD - sD
        if best is None or v < best:
            best = v
    return best


def w_hi_exact(uk, u0, u1):
    """an exact rational UPPER bound on `w = |cos D| + |sin D|`, `D = th - th_k`, over the poses of
    a `u`-interval.  `D in (-90, 90)`, so `w = cos D + |sin D|`, whose only interior maxima are at
    `D = +-45` (value `sqrt 2`).  Whether the interval contains `D = +-45` is decided exactly by
    the sign of `cos D -+ sin D` at the two endpoints (that difference is monotone in `D` on the
    relevant half), and whether it contains `D = 0` by the sign of `sin D`; the maximum is then the
    largest endpoint value, with `1` thrown in when `D = 0` is straddled and the rational
    `1.41422 > sqrt 2` when `+-45` is."""
    ck, sk = cs_of_u(uk)
    vals, cD, sD = [], [], []
    for u in (u0, u1):
        c, s = cs_of_u(u)
        a, b = c * ck + s * sk, s * ck - c * sk
        cD.append(a)
        sD.append(b)
        vals.append(a + (b if b >= 0 else -b))
    if (cD[0] - sD[0] > 0) != (cD[1] - sD[1] > 0):     # D = +45 inside
        return Fr(141422, 100000)
    if (cD[0] + sD[0] > 0) != (cD[1] + sD[1] > 0):     # D = -45 inside
        return Fr(141422, 100000)
    if (sD[0] >= 0) != (sD[1] >= 0):                   # D = 0 inside
        vals.append(Fr(1))
    return max(vals)


def boxes_meet(b, b2):
    """EXACT sufficient criterion that EVERY pose of box `b` closed-meets EVERY pose of box `b2`:
    the largest possible centre separation between the two centre rectangles is at most 1.

    Same argument as `box_pairwise_meets`: every separating-axis margin is at least
    `1 - ||c' - c||_2` because `w >= 1`.  The extreme centre pair is a pair of rectangle corners."""
    d2 = max((x - x2) ** 2 + (y - y2) ** 2
             for x in (b[0], b[1]) for y in (b[2], b[3])
             for x2 in (b2[0], b2[1]) for y2 in (b2[2], b2[3]))
    return d2 <= 1


def triple_empty_certificate(boxes, evecs, angles):
    """EXACT Farkas certificate that NO point lies in all three squares, for EVERY choice of one
    pose from each of the three boxes.

    For a fixed unit direction `e` and a pose `(c', th')` of box `k`, every point `x` of the square
    obeys `<x, e> <= <c', e> + w'/2` with `w' = |cos(th' - phi_e)| + |sin(th' - phi_e)|`, so with
    `C_k = max over the centre rectangle of <c', e_k>  +  w_hi/2` (both exact, the first a max over
    four corners, the second `w_hi_exact`) the square is contained in `{<x, e_k> <= C_k}`.

    Three half planes have empty intersection as soon as there are `lam_k >= 0`, not all zero, with
    `sum lam_k e_k = 0` and `sum lam_k C_k < 0`: a common point `x` would give
    `0 = <x, sum lam_k e_k> = sum lam_k <x, e_k> <= sum lam_k C_k < 0`.  The canonical `lam` for
    three planar vectors is `lam = (e1 x e2, e2 x e0, e0 x e1)` (the identity
    `(e1 x e2) e0 + (e2 x e0) e1 + (e0 x e1) e2 = 0`); it is used with whichever overall sign makes
    every component non-negative.  All rational.

    `evecs[k]` is the exact rational unit outward normal `e_k` to use for square `k` (one of its
    own four edge normals, so the nominal `w` is exactly 1) and `angles[k] = (p, q)` its
    tan-half-angle, which `w_hi_exact` needs (`w` has period 90 in the angle, so the square's own
    parameter serves for any of its four normals)."""
    e = list(evecs)

    def cross(a, b):
        return a[0] * b[1] - a[1] * b[0]

    lam = [cross(e[1], e[2]), cross(e[2], e[0]), cross(e[0], e[1])]
    if all(v <= 0 for v in lam):
        lam = [-v for v in lam]
    if not all(v >= 0 for v in lam) or all(v == 0 for v in lam):
        return False, None
    tot = Fr(0)
    Cs = []
    for k, (bx, (p, q)) in enumerate(zip(boxes, angles)):
        x0, x1, y0, y1, u0, u1 = bx
        ex, ey = e[k]
        C = max(cx * ex + cy * ey for cx in (x0, x1) for cy in (y0, y1))
        C += w_hi_exact(Fr(p, q), u0, u1) / 2
        Cs.append(C)
        tot += lam[k] * C
    return tot < 0, (lam, Cs, tot)


def box_meets_square(box, mem, wlo=None):
    """EXACT: does EVERY pose in the box `B = [x0,x1] x [y0,y1] x [u0,u1]` closed-meet the fixed
    square `mem = (p, q, ax, ay)`?

    The four separating-axis conditions, with `d = (ax - cx, ay - cy)` and `K = 1/2 + w_lo/2`
    (`w_lo` an exact lower bound on `w` over the box, `w_lo_exact`):

        axis 1  (the member's own normal `(C, S)`):  |d.x C + d.y S|      <= K
        axis 2  (the member's other normal):         |-d.x S + d.y C|     <= K
        axis 3  (the pose's own normal):             |d.x cos th + d.y sin th| <= K
        axis 4  (the pose's other normal):           |-d.x sin th + d.y cos th| <= K

    REQUIRES `0 <= u0 <= u1 <= 1`, i.e. `th in [0, 90]` -- the convention everywhere in this repo
    (a unit square is quarter-turn invariant, so that range is all of pose space).  Outside it
    `w_lo_exact`'s branch analysis, which needs `|th - th_k| < 90` so that `cos` of the difference
    is non-negative, is not valid.  `certify_box` enforces the range.

    Axes 1-2 are `|affine in (cx, cy)| <= const`, so the maximum over the centre rectangle is at
    one of its four CORNERS.  Axes 3-4, multiplied through by `1 + u^2 > 0`, become

        +- ( d.x (1 - u^2) + 2 d.y u )  -  K (1 + u^2)  <= 0        (axis 3)
        +- ( -2 d.x u + d.y (1 - u^2) ) -  K (1 + u^2)  <= 0        (axis 4)

    i.e. polynomials AFFINE in `(cx, cy)` and QUADRATIC in `u`, exactly the class of
    `zeromargin.py`'s CHAIN violation polynomials (`search/ZEROMARGIN.md`, `search/RUNG2.md` §6).
    An affine function of `(cx, cy)` is maximised over a rectangle at a corner, so the maximum over
    the whole box is `max over the 4 corners of (exact quadratic max in u)` -- no Bernstein bound,
    no subdivision, no slack.  Everything is `Fraction`.

    Returns True iff all eight polynomials are `<= 0` on the box."""
    x0, x1, y0, y1, u0, u1 = box
    p, q, ax, ay = mem
    C, S = pose_cs(p, q)
    if wlo is None:
        wlo = w_lo_exact(Fr(p, q), u0, u1)
    K = HALF + wlo / 2
    for cx in (x0, x1):
        for cy in (y0, y1):
            dx, dy = ax - cx, ay - cy
            # axes 1, 2: constant in u
            if abs(dx * C + dy * S) > K:
                return False
            if abs(-dx * S + dy * C) > K:
                return False
            # axis 3: +-(dx (1-u^2) + 2 dy u) - K (1+u^2) <= 0
            if _max_quad(dx - K, 2 * dy, -dx - K, u0, u1) > 0:
                return False
            if _max_quad(-dx - K, -2 * dy, dx - K, u0, u1) > 0:
                return False
            # axis 4: +-(-2 dx u + dy (1-u^2)) - K (1+u^2) <= 0
            if _max_quad(dy - K, -2 * dx, -dy - K, u0, u1) > 0:
                return False
            if _max_quad(-dy - K, 2 * dx, dy - K, u0, u1) > 0:
                return False
    return True


def box_pairwise_meets(box):
    """EXACT: do any two poses of the box closed-meet each other?

    For poses `S, S'` in the box with centres `c, c'` and angles `th, th'`, each of the four
    separating-axis margins is `1/2 + w/2 - |<c' - c, e>|` with `|e| = 1` and
    `w = |cos(th' - th)| + |sin(th' - th)| >= 1`, so every margin is at least
    `1 - ||c' - c||_2`.  The worst pair is at OPPOSITE CORNERS of the centre rectangle, where
    `||c' - c||_2 = sqrt((x1-x0)^2 + (y1-y0)^2)`.  Criterion, exact in rationals:

        (x1 - x0)^2 + (y1 - y0)^2 <= 1.

    (The `u`-extent is irrelevant: rotating one square only ever raises `w`.)  This is sufficient,
    not necessary -- a long thin box can also be pairwise meeting -- but it is what the grown
    boxes need, and it is free."""
    x0, x1, y0, y1, u0, u1 = box
    return (x1 - x0) ** 2 + (y1 - y0) ** 2 <= 1


def box_admissible(box, t=T):
    """EXACT: is every pose of the box admissible, i.e. does its closed square lie in `[0,t]^2`?

    `cx, cy in [w(th)/2, t - w(th)/2]` with `w(th) = cos th + sin th` on `[0, 90]`.  `w` is
    concave there, so its MAXIMUM over a `u`-interval is at an endpoint unless the interval
    contains `u = tan(45/2) = sqrt2 - 1`, where `w = sqrt2`.  `sqrt2` is irrational, so the
    straddling case is bounded above by the rational `1.41422` (`> sqrt 2`), which is sound."""
    x0, x1, y0, y1, u0, u1 = box
    U45 = Fr(41421356, 100000000)          # < sqrt2 - 1 < ...357/10^8
    if u0 <= U45 + Fr(1, 100000000) and u1 >= U45:
        # sqrt2 = 1.414213562373...; this rational exceeds it by 4.4e-9, so a pose box may be
        # refused only within 2.2e-9 of the true wall (with 1.41422 it was 3.2e-6, and that alone
        # cost 14 of SUPEPGF's 317 rows their box)
        whi = Fr(141421357, 100000000)
    else:
        whi = max(sum(cs_of_u(u)) for u in (u0, u1))
    return x0 >= whi / 2 and x1 <= t - whi / 2 and y0 >= whi / 2 and y1 <= t - whi / 2


def certify_box(box, members, verbose=False):
    """run the three exact tests on a proposed box: every pose admissible, every pose meets every
    member, and every two poses of the box meet.  Returns (ok, reason)."""
    if not (box[1] > box[0] and box[3] > box[2] and box[5] > box[4]):
        return False, 'degenerate'
    if not (0 <= box[4] and box[5] <= 1):
        return False, 'u outside [0, 1] (th outside [0, 90], where w_lo_exact is not valid)'
    if not box_admissible(box):
        return False, 'inadmissible'
    if not box_pairwise_meets(box):
        return False, 'box not pairwise meeting'
    for k, mem in enumerate(members):
        if not box_meets_square(box, mem):
            return False, f'member {k} not met'
    return True, 'ok'


# ====================================================================== float margins
def margins(poses, M):
    """float signed meet margins of each pose against each member: shape (len(poses), len(M)).

    `M` is the float member array `(ax, ay, cos, sin)`.  The margin on an axis is
    `1/2 + w/2 - |<c_j - c_i, e>|`; the four axes are the two squares' own normals; the squares
    meet iff the least of the four is `>= 0`.  Same arithmetic as
    `honestcost.pair_meets` / `cliquelever.Poses.meet_margin`."""
    P = np.asarray(poses, dtype=float).reshape(-1, 3)
    ci = np.cos(P[:, 2])[:, None]
    si = np.sin(P[:, 2])[:, None]
    cj = M[None, :, 2]
    sj = M[None, :, 3]
    ddx = M[None, :, 0] - P[:, 0][:, None]
    ddy = M[None, :, 1] - P[:, 1][:, None]
    cr = ci * cj + si * sj
    sr = ci * sj - si * cj
    w = 0.5 * (np.abs(cr) + np.abs(sr))
    m1 = 0.5 + w - np.abs(ddx * ci + ddy * si)
    m2 = 0.5 + w - np.abs(-ddx * si + ddy * ci)
    m3 = 0.5 + w - np.abs(ddx * cj + ddy * sj)
    m4 = 0.5 + w - np.abs(-ddx * sj + ddy * cj)
    return np.minimum(np.minimum(m1, m2), np.minimum(m3, m4))


def maxmin_at_theta(M, th, x0=0.0, x1=4.0, y0=0.0, y1=4.0, sub=None, maxit=12):
    """the exact-for-this-theta answer to "how robustly can one square at angle `th` meet every
    member?": maximise `t` over centres `(cx, cy)` subject to every one of the 8 half-planes of
    every member holding with slack `t`.

    At fixed `th` each margin is `1/2 + w_k/2 - |<(a_k) - c, e>|` with `w_k` a constant, so the
    constraint `margin >= t` is a pair of LINEAR inequalities in `(cx, cy, t)` -- the credited
    region is an intersection of octagons and the max-min is a 3-variable LP.  Solved by cutting
    planes (start from a subset of members, add the violated ones) because `E2Pg`'s rows have up
    to 1,806 members.  Float (HiGHS)."""
    c, s = math.cos(th), math.sin(th)
    cj, sj = M[:, 2], M[:, 3]
    cr = c * cj + s * sj
    sr = c * sj - s * cj
    w = 0.5 * (np.abs(cr) + np.abs(sr))
    rhs0 = 0.5 + w                                   # |<a - c, e>| <= rhs0 - t
    # the 4 axes as unit vectors, per member: (c,s), (-s,c), (cj,sj), (-sj,cj)
    AX = np.stack([np.full(len(M), c), np.full(len(M), -s), cj, -sj], axis=1)
    AY = np.stack([np.full(len(M), s), np.full(len(M), c), sj, cj], axis=1)
    PR = M[:, 0:1] * AX + M[:, 1:2] * AY             # <a_k, e> per (member, axis)
    n = len(M)
    idx = np.arange(n) if sub is None else np.asarray(sub)
    # the centre must stay ADMISSIBLE at this angle: cx, cy in [w(th)/2, t - w(th)/2]
    w2 = 0.5 * (abs(c) + abs(s))
    bounds = [(max(x0, w2), min(x1, 4.0 - w2)), (max(y0, w2), min(y1, 4.0 - w2)), (None, None)]
    for _ in range(maxit):
        rows, rhs = [], []
        for k in idx:
            for j in range(4):
                # <c, e> >= <a,e> - rhs0 + t  ->  -ax cx - ay cy + t <= rhs0 - PR
                rows.append([-AX[k, j], -AY[k, j], 1.0])
                rhs.append(rhs0[k] - PR[k, j])
                #  <c, e> <= <a,e> + rhs0 - t  ->  ax cx + ay cy + t <= rhs0 + PR
                rows.append([AX[k, j], AY[k, j], 1.0])
                rhs.append(rhs0[k] + PR[k, j])
        r = linprog(c=[0.0, 0.0, -1.0], A_ub=np.array(rows), b_ub=np.array(rhs),
                    bounds=bounds, method='highs')
        if not r.success:
            return -1e9, None
        cx, cy, t = r.x
        viol = np.minimum(rhs0 - np.abs(PR - (cx * AX + cy * AY)).max(axis=1) - t, 0.0)
        bad = np.nonzero(viol < -1e-12)[0]
        if not len(bad):
            return float(t), (cx, cy, th)
        idx = np.union1d(idx, bad[np.argsort(viol[bad])[:200]])
    return float(t), (cx, cy, th)


def best_allmeet_pose(M, nth=181, refine=3, sub0=64, log=None):
    """maximise over ALL admissible poses the least meet margin against the members: a sweep of
    `maxmin_at_theta` over `th in [0, 90)` plus a local bisection refinement in `th`.

    If the answer is `<= 0`, no pose meets every member with room to spare, `A(K)` has empty
    interior, and no positive-volume box exists.  If it is `> 0`, a whole neighbourhood of the
    maximiser meets every member and a positive-volume box does exist.  Float."""
    n = len(M)
    sub = None if n <= sub0 else np.argsort(-M[:, 0] - M[:, 1])[:sub0]
    ths = np.linspace(0.0, math.pi / 2, nth, endpoint=False)
    vals, poses = [], []
    for th in ths:
        v, p = maxmin_at_theta(M, th, sub=sub)
        vals.append(v)
        poses.append(p)
    vals = np.array(vals)
    i = int(np.argmax(vals))
    best, bp = float(vals[i]), poses[i]
    h = (ths[1] - ths[0]) if nth > 1 else 0.0
    for _ in range(refine):
        h *= 0.5
        for th in (bp[2] - h, bp[2] + h):
            if not (0.0 <= th < math.pi / 2):
                continue
            v, p = maxmin_at_theta(M, th, sub=sub)
            if v > best:
                best, bp = v, p
    return best, bp


# ====================================================================== growing a certified box
def _maxquad_v(a0, a1, a2, u0, u1):
    """vectorised float max of `a0 + a1 u + a2 u^2` over `[u0, u1]` (`_max_quad`, elementwise)"""
    v = np.maximum(a0 + a1 * u0 + a2 * u0 * u0, a0 + a1 * u1 + a2 * u1 * u1)
    with np.errstate(divide='ignore', invalid='ignore'):
        uv = np.where(a2 < 0, -a1 / (2 * a2), u0 - 1.0)
    inr = (uv > u0) & (uv < u1)
    return np.where(inr, a0 + a1 * uv + a2 * uv * uv, v)


class FloatBoxTest:
    """a fast FLOAT surrogate for `box_meets_square` over a whole member list at once, used to
    steer the search inside `grow_box`.  It evaluates exactly the same eight polynomials at the
    same four centre-rectangle corners, in doubles and vectorised over the members.

    It decides nothing.  `grow_box` exactly certifies the box it finally returns (shrinking it
    until the exact certificate holds, should a double have rounded the wrong way), and
    `report --verify` re-runs the exact certificate on every box in the json."""

    def __init__(self, members, tol=1e-12):
        self.ax = np.array([float(m[2]) for m in members])
        self.ay = np.array([float(m[3]) for m in members])
        self.uk = np.array([float(Fr(m[0], m[1])) for m in members])
        cs = np.array([[float(v) for v in pose_cs(m[0], m[1])] for m in members])
        self.C, self.S = cs[:, 0], cs[:, 1]
        self.tol = tol

    def wlo(self, u0, u1):
        """the float image of `w_lo_exact`, vectorised over the members"""
        uk = self.uk
        ck, sk = (1 - uk ** 2) / (1 + uk ** 2), 2 * uk / (1 + uk ** 2)
        w = None
        for u in (u0, u1):
            c, s = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
            cD, sD = c * ck + s * sk, s * ck - c * sk
            v = np.where(u0 >= uk, cD + sD, cD - sD)
            w = v if w is None else np.minimum(w, v)
        return np.where((u0 <= uk) & (uk <= u1), 1.0, w)

    def ok(self, box):
        x0, x1, y0, y1, u0, u1 = (float(v) for v in box)
        if not (0.0 <= u0 and u1 <= 1.0 and x1 > x0 and y1 > y0 and u1 > u0):
            return False
        K = 0.5 + self.wlo(u0, u1) / 2
        t = -self.tol
        for cx in (x0, x1):
            for cy in (y0, y1):
                dx, dy = self.ax - cx, self.ay - cy
                if (np.abs(dx * self.C + dy * self.S) - K > t).any():
                    return False
                if (np.abs(-dx * self.S + dy * self.C) - K > t).any():
                    return False
                for (a0, a1, a2) in ((dx - K, 2 * dy, -dx - K), (-dx - K, -2 * dy, dx - K),
                                     (dy - K, -2 * dx, -dy - K), (-dy - K, 2 * dx, dy - K)):
                    if (_maxquad_v(a0, a1, a2, u0, u1) > t).any():
                        return False
        return True


def frac_near(x, den):
    return Fr(int(round(x * den)), den)


def grow_box(centre, members, den=10 ** 9, seed=1e-7, rounds=4, cap=(4, 4, 1)):
    """a maximal-ish certified box around `centre = (cx, cy, u)`: seed a tiny box that passes
    `certify_box`, then push each of its SIX faces outwards as far as the exact certificate
    allows (double until it fails, then bisect).

    Six independent faces rather than one symmetric half-side, because the best all-meeting pose
    of a row is often pressed against the admissibility wall (`cx = w(th)/2`), where no symmetric
    box is admissible but a one-sided one is.  The seed is likewise tried in all four `(+-, +-)`
    centre directions.

    Float search, EXACT verdict: every candidate box is tested by `certify_box`, and the box
    returned is the last one that passed."""
    cx, cy, cu = (frac_near(v, den) for v in centre)
    e = frac_near(seed, den)
    if e <= 0:
        e = Fr(1, den)

    fast = FloatBoxTest(members)

    def exact(b):
        return certify_box(b, members)[0]

    def ok(b):
        # the search test: the exact structural checks (cheap) plus the float surrogate
        if not (b[1] > b[0] and b[3] > b[2] and b[5] > b[4] and 0 <= b[4] and b[5] <= 1):
            return False
        return box_admissible(b) and box_pairwise_meets(b) and fast.ok(b)

    # (a) a BALANCED seed: the largest `centre +- s (1, 1, 1)` that certifies, by bisection.  It
    # has to come first -- a greedy face sweep from a degenerate seed spends all of a row's meet
    # slack on whichever face it happens to push first, and returns a sliver.
    def sym(v):
        u0, u1 = cu - v, cu + v
        return (cx - v, cx + v, cy - v, cy + v, max(u0, Fr(0)), min(u1, Fr(1)))

    box = None
    lo, hi = Fr(0), frac_near(0.5, den)
    if ok(sym(hi)):
        lo = hi
    else:
        for _ in range(40):
            m = frac_near(float(lo + hi) / 2, den)
            if m <= lo or m >= hi:
                break
            if ok(sym(m)):
                lo = m
            else:
                hi = m
    if lo > 0:
        box = list(sym(lo))
    else:
        # (b) the centre is on an admissibility wall (or the row is that tight): a one-sided seed
        for sx in (1, -1):
            for sy in (1, -1):
                x0, x1 = (cx, cx + e) if sx > 0 else (cx - e, cx)
                y0, y1 = (cy, cy + e) if sy > 0 else (cy - e, cy)
                u0, u1 = (cu, cu + e) if cu + e <= 1 else (cu - e, cu)
                u0, u1 = max(u0, Fr(0)), min(u1, Fr(1))
                if u1 <= u0:
                    continue
                if ok((x0, x1, y0, y1, u0, u1)):
                    box = [x0, x1, y0, y1, u0, u1]
                    break
            if box:
                break
    if box is None:
        return None
    # face k of coordinate c: 0 -> lower (decrease), 1 -> upper (increase)
    faces = [(0, 0, Fr(0)), (1, 0, Fr(cap[0])), (2, 0, Fr(0)), (3, 0, Fr(cap[1])),
             (4, 0, Fr(0)), (5, 0, Fr(cap[2]))]
    for _ in range(rounds):
        moved = False
        for (k, _z, lim) in faces:
            sgn = -1 if k % 2 == 0 else 1
            step = frac_near(1e-4, den)
            base = box[k]
            grown = base
            for _ in range(28):                      # double outwards while the certificate holds
                cand = list(box)
                nxt = grown + sgn * step
                if (sgn < 0 and nxt < lim) or (sgn > 0 and nxt > lim):
                    nxt = lim
                if nxt == grown:
                    break
                cand[k] = nxt
                if ok(tuple(cand)):
                    grown = nxt
                    if nxt == lim:
                        break
                    step *= 2
                else:
                    break
            hiq = grown + sgn * step                 # bisect the last failed step
            if (sgn < 0 and hiq < lim) or (sgn > 0 and hiq > lim):
                hiq = lim
            loq = grown
            for _ in range(26):
                m = frac_near(float(loq + hiq) / 2, den)
                if m == loq or m == hiq:
                    break
                cand = list(box)
                cand[k] = m
                if ok(tuple(cand)):
                    loq = m
                else:
                    hiq = m
            if loq != base:
                box[k] = loq
                moved = True
        if not moved:
            break
    # the surrogate steered the search; the EXACT certificate decides.  If a double rounded the
    # wrong way anywhere, shrink the box towards its centre until the exact test passes.
    box = tuple(box)
    mid = [(box[0] + box[1]) / 2, (box[2] + box[3]) / 2, (box[4] + box[5]) / 2]
    for _ in range(60):
        if exact(box):
            return box
        box = tuple(mid[k // 2] + (box[k] - mid[k // 2]) * Fr(4, 5) for k in range(6))
        if box[1] <= box[0] or box[3] <= box[2] or box[5] <= box[4]:
            return None
    return None


def box_vol(box):
    return float((box[1] - box[0]) * (box[3] - box[2]) * (box[5] - box[4]))


def box_min_margin(box, M, nsamp=64, rng=None):
    """float least meet margin over a random sample of the box's poses against every member (a
    sanity companion to the exact certificate, and the number the tables quote)"""
    rng = rng or np.random.default_rng(20260912)
    x0, x1, y0, y1, u0, u1 = (float(v) for v in box)
    corners = np.array([[x, y, 2 * math.atan(u)] for x in (x0, x1) for y in (y0, y1)
                        for u in (u0, u1)])
    S = np.column_stack([rng.uniform(x0, x1, nsamp), rng.uniform(y0, y1, nsamp),
                         2 * np.arctan(rng.uniform(u0, u1, nsamp))])
    return float(margins(np.vstack([corners, S]), M).min())


# ====================================================================== per-row anatomy
def core_deficit(M, sub=None):
    """the max-margin LP over the member half-planes: the largest `t` with a point `x` at distance
    `>= t` inside EVERY member square.  `t > 0` means the row has a fat common core (Helly),
    `t = 0` a single common point, `t < 0` an empty core with that deficit.  This re-derives
    `HONEST.md` §0 item 2's `-0.0007 .. -0.042` independently of `honestcost.clique_core`'s
    polygon clipping.  Float."""
    idx = np.arange(len(M)) if sub is None else np.asarray(sub)
    for _ in range(12):
        rows, rhs = [], []
        for k in idx:
            ax, ay, c, s = M[k]
            for (ex, ey) in ((c, s), (-s, c)):
                rows.append([ex, ey, 1.0])
                rhs.append(0.5 + ax * ex + ay * ey)
                rows.append([-ex, -ey, 1.0])
                rhs.append(0.5 - ax * ex - ay * ey)
        r = linprog(c=[0.0, 0.0, -1.0], A_ub=np.array(rows), b_ub=np.array(rhs),
                    bounds=[(0, 4), (0, 4), (None, None)], method='highs')
        if not r.success:
            return float('nan')
        x, y, t = r.x
        dx, dy = M[:, 0] - x, M[:, 1] - y
        u = np.abs(dx * M[:, 2] + dy * M[:, 3])
        v = np.abs(-dx * M[:, 3] + dy * M[:, 2])
        slack = 0.5 - np.maximum(u, v)
        bad = np.nonzero(slack < t - 1e-12)[0]
        if not len(bad):
            return float(t)
        idx = np.union1d(idx, bad[np.argsort(slack[bad])[:200]])
    return float(t)


def grazing_anatomy(M, tol=1e-6, cap=600, rng=None):
    """the pairwise meet-margin structure of the row: how many member pairs GRAZE (`|margin| <
    tol`), and is the grazing graph a FAN (a few members touched by very many)?  On rows with more
    than `cap` members a random subsample of `cap` is used and the counts are reported as rates."""
    n = len(M)
    sel = np.arange(n)
    if n > cap:
        rng = rng or np.random.default_rng(20260912)
        sel = np.sort(rng.choice(n, cap, replace=False))
    Ms = M[sel]
    G = margins(np.column_stack([Ms[:, 0], Ms[:, 1], np.arctan2(Ms[:, 3], Ms[:, 2])]), Ms)
    G = np.minimum(G, G.T)
    np.fill_diagonal(G, 1.0)
    graze = np.abs(G) < tol
    m = len(sel)
    npair = m * (m - 1) // 2
    ngraze = int(graze.sum() // 2)
    deg = graze.sum(axis=1)
    o = np.argsort(-deg)
    return dict(n_members=n, n_sampled=m, n_pairs=npair, n_grazing=ngraze,
                graze_frac=float(ngraze / npair) if npair else 0.0,
                min_pair_margin=float(G[np.triu_indices(m, 1)].min()) if m > 1 else 0.0,
                max_degree=int(deg[o[0]]) if m else 0,
                top_degrees=[int(deg[i]) for i in o[:6]],
                deg_mean=float(deg.mean()) if m else 0.0,
                n_isolated=int((deg == 0).sum()),
                # how tight the row really is: the pairwise margin never reaches 0 on these rows,
                # so a single `graze` tolerance says nothing -- quote the distribution instead
                pair_quantiles=[float(v) for v in
                                np.quantile(G[np.triu_indices(m, 1)], [0, .01, .25, .5])]
                if m > 1 else [],
                n_below=[int((G[np.triu_indices(m, 1)] < t).sum())
                         for t in (1e-6, 1e-4, 1e-3, 1e-2)] if m > 1 else [])


# ====================================================================== cover with the box rule
class BoxCover(hc.Cover):
    """`honestcost.Cover` with the SOUND clique rule of this task: `S` is credited `z_K` iff
    `S in K u B_K`, i.e. `S` is (exactly) one of the row's member poses, or `S` lies in the
    certified box `B_K`.  `K u B_K` is pairwise meeting by construction (`certify_box`), so
    crediting this way IS a valid weighting -- unlike `honestcost`'s `meet`."""

    def attach_boxes(self, boxes, rowsel=None, members=True):
        """`boxes[j]` = the float box of the `j`-th row of THIS cover object (or None).  `rowsel`,
        when given, is the index in the parent cover of each of this object's rows -- it is what
        keeps the box list in step with `honestcost.one_row_cover`'s slicing of `cz`/`cmem`.
        `members` switches on the `S in K` half of the rule, which needs a dictionary lookup per
        pose and is only worth paying for on the pose set `P` itself (elsewhere an exact
        coincidence with a member pose is a measure-zero event that never happens)."""
        self.boxes = list(boxes)
        self.rowsel = list(range(len(self.boxes))) if rowsel is None else list(rowsel)
        self.members = members
        self.memz = {}
        if members:
            for j in range(self.n_cliques):
                for i in self.cmem[j]:
                    k = tuple(self.F[i])
                    self.memz[k] = self.memz.get(k, 0.0) + float(self.cz[j])
        return self

    def subrow(self, j):
        """the one-row restriction of this cover, boxes included (`one_row_cover` alone would
        leave `self.boxes[0]` pointing at the wrong row)"""
        sub = BoxCover.__new__(BoxCover)      # NOT hc.one_row_cover: that builds a plain Cover,
        sub.__dict__.update(self.__dict__)    # whose clique_capture would silently fall back to
        sub.cz = self.cz[j:j + 1]             # the `meet` rule and report the wrong rows
        sub.cmem = self.cmem[j:j + 1]
        sub.ccore = self.ccore[j:j + 1]
        sub.n_cliques = 1
        sub.pz = np.zeros(0)
        sub.ppts = []
        sub.n_pgons = 0
        sub.boxes = self.boxes[j:j + 1]
        sub.rowsel = self.rowsel[j:j + 1]
        sub.members = False
        sub.memz = {}
        return sub

    def clique_capture(self, poses, rule='boxK'):
        if rule != 'boxK':
            return super().clique_capture(poses, rule)
        P = np.ascontiguousarray(np.asarray(poses, dtype=float).reshape(-1, 3))
        out = np.zeros(len(P))
        if not len(P) or not self.n_cliques:
            return out
        u = np.tan(0.5 * np.mod(P[:, 2], math.pi / 2))
        for j in range(self.n_cliques):
            b = self.boxes[j]
            if b is None:
                continue
            inb = ((P[:, 0] >= b[0]) & (P[:, 0] <= b[1]) & (P[:, 1] >= b[2])
                   & (P[:, 1] <= b[3]) & (u >= b[4]) & (u <= b[5]))
            out[inb] += self.cz[j]
        if self.members and self.memz:
            for i in range(len(P)):
                k = (P[i, 0], P[i, 1], math.cos(P[i, 2]), math.sin(P[i, 2]))
                z = self.memz.get(k)
                if z is not None:
                    out[i] += z
        return out


def anatomy_box(cov, pose):
    """`honestcost.anatomy` for a `BoxCover` (it needs `subrow`, not `one_row_cover`)"""
    P = np.array(pose, dtype=float).reshape(1, 3)
    tot, pt, cq, pg = cov.capture(P, 'boxK', 1, parts=True)
    rows = [(int(j), float(cov.cz[j])) for j in range(cov.n_cliques)
            if cov.subrow(j).clique_capture(P, 'boxK')[0] > 0]
    inp = []
    for i in range(cov.n_atoms):
        dx, dy = cov.ax[i] - P[0, 0], cov.ay[i] - P[0, 1]
        co, si = math.cos(P[0, 2]), math.sin(P[0, 2])
        if max(abs(dx * co + dy * si), abs(-dx * si + dy * co)) <= 0.5 + 1e-8:
            inp.append((float(cov.ax[i]), float(cov.ay[i]), float(cov.aw[i])))
    inp.sort(key=lambda r: -r[2])
    return dict(pose=[float(v) for v in P[0]], theta_deg=math.degrees(P[0, 2]) % 90.0,
                capture=float(tot[0]), points=float(pt[0]), cliques=float(cq[0]),
                pgons=float(pg[0]), n_points=len(inp), top_points=inp[:8],
                n_cliques=len(rows), clique_rows=rows[:10])


# ====================================================================== commands
def _logger(t0):
    def log(m):
        print(f'[{time.time() - t0:7.0f}s] {m}', flush=True)
    return log


def cmd_boxes(a):
    t0 = time.time()
    log = _logger(t0)
    d = np.load(os.path.join(RUNS, f'hc_{a.TAG}_dual.npz'), allow_pickle=True)
    cov = hc.Cover(d, a.geo_tol, log=log)
    log(f'# allmeet boxes {a.TAG}')
    log(cov.describe())
    PE = read_exact_poses(a.TAG, a.base)
    FE = np.array([[float(cx), float(cy), *[float(v) for v in pose_cs(p, q)]]
                   for (p, q, cx, cy) in PE])
    assert FE.shape == cov.F.shape and np.abs(FE - cov.F).max() == 0.0, \
        'exact pose reconstruction does not match the dual dump'
    log(f'   {len(PE)} exact poses reconstructed and matched to the dual dump bit for bit')
    nfold = sum(1 for q in PE if Fr(q[0], q[1]) < 0)
    PE = [fold_pose(q) for q in PE]
    log(f'   {nfold} of them re-parametrised to u in [0,1) by a quarter turn (same square)')

    order = np.argsort(-cov.cz)
    rows = []
    nbox = 0
    massbox = 0.0
    for rank, ci in enumerate(order[:a.nrows]):
        ci = int(ci)
        mem = cov.cmem[ci]
        M = cov.F[mem]
        z = float(cov.cz[ci])
        best, bp = best_allmeet_pose(M, nth=a.nth, sub0=a.sub0)
        rec = dict(row=ci, rank=rank, z=z, size=int(len(mem)), maxmin_margin=best,
                   best_pose=None if bp is None else [float(v) for v in bp])
        box = None
        if best > a.mtol and bp is not None:
            # tightest member first: `certify_box` exits on the first failure, and the bisection
            # spends most of its calls on boxes that fail
            om = np.argsort(margins(np.array([bp]), M)[0])
            members = [PE[int(mem[k])] for k in om]
            box = grow_box((bp[0], bp[1], math.tan(0.5 * bp[2])), members, den=a.den)
            if box is None:
                # the maximiser sits ON the admissibility wall `c = w(th)/2` (a 45-degree pose has
                # it at the irrational `sqrt2 / 2`), where no box centred there is admissible.
                # Nudge the centre into the strict interior by `pad` and retry once.
                lo, hi = 0.70710679 + a.pad, 3.29289321 - a.pad
                box = grow_box((min(max(bp[0], lo), hi), min(max(bp[1], lo), hi),
                                math.tan(0.5 * bp[2])), members, den=a.den)
            if box is not None:
                ok, why = certify_box(box, members)
                assert ok, (ci, why)
                rec['box'] = [[str(v) for v in box], [float(v) for v in box]]
                rec['vol'] = box_vol(box)
                rec['box_min_margin'] = box_min_margin(box, M)
                rec['box_sides'] = [float(box[1] - box[0]), float(box[3] - box[2]),
                                    float(box[5] - box[4])]
                nbox += 1
                massbox += z
        if box is None:
            rec['vol'] = 0.0
        if rank < a.nanat:
            rec['core_deficit'] = core_deficit(M, sub=None if len(M) <= a.sub0
                                               else np.arange(min(len(M), a.sub0)))
            rec['graze'] = grazing_anatomy(M, tol=a.gtol, cap=a.cap)
        rows.append(rec)
        if rank % 10 == 0 or box is not None:
            log(f'   row {ci:4d} rank {rank:3d} z={z:.6f} |K|={len(mem):5d} '
                f'maxmin={best:+.3e} vol={rec["vol"]:.3e}')
    out = dict(tag=a.TAG, theta=cov.theta, lp=cov.lp, clique_total=cov.clique_total,
               n_rows=len(rows), n_with_box=nbox, mass_with_box=massbox,
               mass_without_box=float(sum(r['z'] for r in rows)) - massbox, rows=rows)
    log(f'   {nbox} of {len(rows)} rows have vol(B_K) > 0, carrying {massbox:.6f} of '
        f'{cov.clique_total:.6f} clique dual')
    p = os.path.join(RUNS, f'am_{a.TAG}_boxes.json')
    json.dump(out, open(p, 'w'), indent=1)
    log(f'   wrote {p}')


def cmd_honest(a):
    t0 = time.time()
    log = _logger(t0)
    d = np.load(os.path.join(RUNS, f'hc_{a.TAG}_dual.npz'), allow_pickle=True)
    B = json.load(open(os.path.join(RUNS, f'am_{a.TAG}_boxes.json')))
    cov = BoxCover(d, a.geo_tol, log=log)
    boxes = [None] * cov.n_cliques
    for r in B['rows']:
        if r.get('box'):
            boxes[r['row']] = [float(v) for v in r['box'][1]]
    cov.attach_boxes(boxes, members=True)
    log(f'# allmeet honest {a.TAG}  (sound rule: credit z_K iff S in K u B_K)')
    log(cov.describe())
    log(f'   {sum(b is not None for b in boxes)} rows carry a certified positive-volume box, '
        f'carrying {sum(float(cov.cz[j]) for j in range(cov.n_cliques) if boxes[j]):.6f} of '
        f'{cov.clique_total:.6f} clique dual')

    def honest_of(v):
        return cov.theta / v if v > 0 else float('inf')

    out = dict(tag=a.TAG, rule='boxK', theta=cov.theta, lp=cov.lp,
               n_boxes=sum(b is not None for b in boxes))
    PF = d['pf']
    PP = np.column_stack([PF[:, 0], PF[:, 1], np.arctan2(PF[:, 3], PF[:, 2])])
    VP = cov.capture(PP, 'boxK', a.threads)
    out['poses_of_P'] = dict(n=len(PP), min=float(VP.min()), honest=honest_of(float(VP.min())),
                             under_one=int((VP < 1 - 1e-9).sum()))
    log(f'   SELFTEST: min capture over the {len(PP)} poses of P = {VP.min():.9f} '
        f'({int((VP < 1 - 1e-9).sum())} of {len(PP)} under 1 - 1e-9)')
    cov.members = False        # off the pose set an exact member coincidence never happens
    allP, allV = [PP], [VP]
    L = hc.lattice(a.pitch, a.dth)
    VL = cov.capture(L, 'boxK', a.threads)
    out['lattice'] = dict(n=len(L), min=float(VL.min()), honest=honest_of(float(VL.min())))
    log(f'   lattice: {len(L)} poses, min capture {VL.min():.9f} -> honest '
        f'{honest_of(float(VL.min())):.6f}')
    o = np.argsort(VL)[:a.nseed]
    allP.append(L[o])
    allV.append(VL[o])
    if a.families and os.path.exists(a.families):
        FP = np.array([[float(v) for v in q[:3]] for q in
                       (line.split() for line in open(a.families))
                       if len(q) >= 3 and not q[0].startswith('#')])
        FP = hc.admissible_clamp(FP)
        FP[:, 2] = np.mod(FP[:, 2], math.pi / 2)
        FP = hc.admissible_clamp(FP)
        VF = cov.capture(FP, 'boxK', a.threads)
        out['families'] = dict(n=len(FP), min=float(VF.min()), honest=honest_of(float(VF.min())))
        log(f'   families: {len(FP)} poses, min capture {VF.min():.9f}')
        o = np.argsort(VF)[:a.nseed]
        allP.append(FP[o])
        allV.append(VF[o])
    seeds = np.concatenate(allP)
    sv = np.concatenate(allV)
    o = np.argsort(sv)[:a.ndescend]
    DP, DV = hc.pattern_search(cov, seeds[o], 'boxK', a.threads, log=log)
    out['descent'] = dict(n=len(DP), min=float(DV.min()), honest=honest_of(float(DV.min())))
    FP2, FV2 = hc.eps_refine(cov, DP[np.argsort(DV)[:a.nfine]], 'boxK', a.threads, log=log)
    out['knife_edge'] = dict(n=len(FP2), min=float(FV2.min()), honest=honest_of(float(FV2.min())))
    DP = np.concatenate([DP, FP2])
    DV = np.concatenate([DV, FV2])
    i = int(np.argmin(DV))
    mn = float(DV[i])
    out['min_capture'] = mn
    out['honest'] = honest_of(mn)
    out['minimiser'] = anatomy_box(cov, DP[i])
    log(f'   MIN capture {mn:.9f} at {DP[i].tolist()} -> HONEST = {honest_of(mn):.6f}')

    # ---- how much of this is the PLACEMENT of the boxes?  Re-grow every box that can be re-grown
    # AROUND the minimiser (the placement most favourable to the cover at that one pose) and
    # re-descend.  Any single pose can be credited by every row whose A(K) it lies in, so one
    # round of this is an optimistic bound on what re-placing the boxes can buy; the new minimum
    # moves elsewhere, which is the point.
    out['refit'] = []
    PE = [fold_pose(q) for q in read_exact_poses(a.TAG, a.base)] if a.refit else None
    for rnd in range(a.refit):
        S = DP[int(np.argmin(DV))]
        mg = margins(np.array([S]), cov.F)[0]
        nre = 0
        for j in range(cov.n_cliques):
            mem = cov.cmem[j]
            if mg[mem].min() <= 1e-9:
                continue
            om = np.argsort(mg[mem])
            b = grow_box((S[0], S[1], math.tan(0.5 * (S[2] % (math.pi / 2)))),
                         [PE[int(mem[k])] for k in om], den=10 ** 9)
            if b is not None:
                cov.boxes[j] = [float(v) for v in b]
                nre += 1
        log(f'   refit round {rnd}: {nre} boxes re-grown around {S.tolist()}')
        RP, RV = hc.pattern_search(cov, np.concatenate([DP[np.argsort(DV)[:a.ndescend]],
                                                        L[np.argsort(VL)[:a.nseed]]]),
                                   'boxK', a.threads, log=log)
        RP2, RV2 = hc.eps_refine(cov, RP[np.argsort(RV)[:a.nfine]], 'boxK', a.threads, log=log)
        DP = np.concatenate([RP, RP2])
        DV = np.concatenate([RV, RV2])
        k = int(np.argmin(DV))
        out['refit'].append(dict(round=rnd, n_refit=nre, min=float(DV[k]),
                                 honest=honest_of(float(DV[k])),
                                 pose=[float(v) for v in DP[k]]))
        log(f'   refit round {rnd}: MIN capture {DV[k]:.9f} -> HONEST '
            f'{honest_of(float(DV[k])):.6f}')
    p = os.path.join(RUNS, f'am_{a.TAG}_honest.json')
    json.dump(out, open(p, 'w'), indent=1)
    log(f'   wrote {p}')


# ====================================================================== the Helly question
def _sq_from_float(cx, cy, th_deg, den=10 ** 6, udev=10 ** 6):
    """an exact pose `(p, q, cx, cy)` near a float one"""
    u = math.tan(math.radians(th_deg) / 2)
    f = Fr(int(round(u * udev)), udev)
    return (f.numerator, f.denominator, Fr(int(round(cx * den)), den),
            Fr(int(round(cy * den)), den))


def cmd_helly(a):
    """Is there a ROBUST non-Helly triple at `t = 4`?  I.e. a POSITIVE-MEASURE family of closed
    unit squares in `[0,4]^2`, pairwise closed-meeting, with no point common to all of them --
    which is exactly the negation of the Helly-type statement `TODO.md` puts up for proof.

    Construction.  Let `n_k`, `k = 0, 1, 2`, be three unit vectors at `120` degrees, so
    `n_0 + n_1 + n_2 = 0`.  For `d > 0` the half planes `H_k = { x : <x, n_k> <= -d }` have empty
    triple intersection (a common point would give `0 = sum_k <x, n_k> <= -3d < 0`) while any two
    meet in a full wedge whose apex is within `O(d)` of the origin.  Take `S_k` = the unit square
    with one EDGE on the line `<x, n_k> = -d`, lying on the `H_k` side.  Then
    `S_0 n S_1 n S_2 = {}`, and consecutive squares overlap in a region of positive area because
    their wedge apices are all near the origin.  Translate by `(2, 2)` into `[0,4]^2`.

    What is certified, all in `Fraction`s and all on the THICKENED family
    `F = B_0 u B_1 u B_2` (three positive-volume pose boxes), not merely on the three nominal
    squares:

    1. every pose of every `B_k` is admissible                       (`box_admissible`);
    2. every two poses of `F` closed-meet: within a box by
       `box_pairwise_meets`, across boxes by `boxes_meet`;
    3. no point lies in all three squares, for EVERY choice of one
       pose from each box                                            (`triple_empty_certificate`).

    `F` then has positive measure in pose space, is pairwise closed-meeting, and has no common
    point.  (`sq_meets_sq` and an exact 12-half-plane clip of the nominal triple are run as
    independent cross-checks of 2 and 3.)"""
    t0 = time.time()
    log = _logger(t0)
    log('# allmeet helly: a robust non-Helly triple of unit squares in [0,4]^2?')
    dd = a.d
    poses, evecs = [], []
    for k in range(3):
        ang = 120.0 * k
        nx, ny = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        poses.append(_sq_from_float(2.0 - (dd + 0.5) * nx, 2.0 - (dd + 0.5) * ny,
                                    ang % 90.0, den=a.den, udev=a.den))
    for k, (p, q, cx, cy) in enumerate(poses):
        C, S = pose_cs(p, q)
        # the exact edge normal of square k closest to n_k = (cos 120k, sin 120k)
        nx, ny = math.cos(math.radians(120.0 * k)), math.sin(math.radians(120.0 * k))
        evecs.append(max(((C, S), (-C, -S), (-S, C), (S, -C)),
                         key=lambda e: float(e[0]) * nx + float(e[1]) * ny))
    sqs = [lc.make_square(cx, cy, p, q, T) for (p, q, cx, cy) in poses]
    M = np.array([[float(cx), float(cy), *[float(v) for v in pose_cs(p, q)]]
                  for (p, q, cx, cy) in poses])
    log(f'   d = {dd}; poses (cx, cy, theta deg): ' +
        ', '.join(f'({float(c[2]):.6f}, {float(c[3]):.6f}, '
                  f'{math.degrees(2 * math.atan(float(Fr(c[0], c[1])))):.4f})' for c in poses))

    # ---- the nominal triple, as a cross-check
    pm = []
    for i in range(3):
        for j in range(i + 1, 3):
            ok = lc.sq_meets_sq(sqs[i], sqs[j])
            mg = float(margins(np.array([[M[i, 0], M[i, 1],
                                          math.atan2(M[i, 3], M[i, 2])]]), M[j:j + 1])[0, 0])
            pm.append(dict(pair=[i, j], meets_exact=bool(ok), margin=mg))
            log(f'   nominal pair {i}{j}: exact sq_meets_sq = {ok}, float margin = {mg:+.6f}')
    HP = []
    for (p, q, cx, cy) in poses:
        C, S = pose_cs(p, q)
        for (ex, ey) in ((C, S), (-C, -S), (-S, C), (S, -C)):
            HP.append((ex, ey, HALF + cx * ex + cy * ey))
    Pg = [(Fr(0), Fr(0)), (T, Fr(0)), (T, T), (Fr(0), T)]
    for (ex, ey, c) in HP:
        Pg = hc.poly_clip(Pg, ex, ey, c, Fr(0))
        if not Pg:
            break
    nominal_empty = not Pg
    log(f'   nominal triple intersection empty (exact clip of 12 half planes): {nominal_empty}')

    # ---- the thickened family
    h = Fr(int(round(a.half * a.den)), a.den)
    hu = Fr(int(round(a.halfu * a.den)), a.den)
    def urange(p, q):
        """the box's `u`-interval, shifted to stay inside `[0, 1]` -- outside it the branch
        analysis of `w_lo_exact` / `w_hi_exact` (which needs `|th - th_k| < 90`) is not valid"""
        u = Fr(p, q)
        u0, u1 = u - hu, u + hu
        if u0 < 0:
            u0, u1 = Fr(0), 2 * hu
        if u1 > 1:
            u0, u1 = 1 - 2 * hu, Fr(1)
        return u0, u1

    boxes = [(cx - h, cx + h, cy - h, cy + h, *urange(p, q)) for (p, q, cx, cy) in poses]
    assert all(0 <= b[4] and b[5] <= 1 for b in boxes)
    adm = [box_admissible(b) for b in boxes]
    inside = [box_pairwise_meets(b) for b in boxes]
    across = [[i, j, boxes_meet(boxes[i], boxes[j])] for i in range(3) for j in range(i + 1, 3)]
    empty, cert = triple_empty_certificate(boxes, evecs, [(p, q) for (p, q, _, _) in poses])
    vols = [box_vol(b) for b in boxes]
    log(f'   box half-sides {float(h)} (centre) / {float(hu)} (u); volumes ' +
        ', '.join(f'{v:.3e}' for v in vols))
    log(f'   every pose admissible: {adm}')
    log(f'   pairwise meeting inside each box: {inside}; across boxes: '
        f'{[c[2] for c in across]}')
    log(f'   triple intersection empty for EVERY choice of poses (Farkas, exact): {empty}' +
        (f'   [sum lam_k C_k = {float(cert[2]):+.6f} < 0]' if cert else ''))
    ok = (all(adm) and all(inside) and all(c[2] for c in across) and empty
          and all(v > 0 for v in vols))
    log(f'   ROBUST NON-HELLY FAMILY AT t = 4: {"YES" if ok else "NO"}')
    out = dict(d=dd, poses=[[str(v) for v in p] for p in poses],
               poses_float=[[float(p[2]), float(p[3]),
                             math.degrees(2 * math.atan(float(Fr(p[0], p[1]))))] for p in poses],
               nominal_pairs=pm, nominal_triple_empty=bool(nominal_empty),
               boxes=[[str(v) for v in b] for b in boxes], box_vols=vols,
               admissible=adm, pairwise_in_box=inside,
               pairwise_across=[[c[0], c[1], bool(c[2])] for c in across],
               triple_empty_all=bool(empty),
               farkas=None if not cert else dict(lam=[str(v) for v in cert[0]],
                                                 C=[str(v) for v in cert[1]],
                                                 total=str(cert[2]), total_float=float(cert[2])),
               robust=bool(ok))
    pth = os.path.join(RUNS, 'am_helly.json')
    json.dump(out, open(pth, 'w'), indent=1)
    log(f'   wrote {pth}')


def pose_volume_sample(n, t=T, rng=None):
    """`n` poses drawn uniformly from the admissible set `{(cx, cy, u) : the closed square lies in
    [0,t]^2, u in [0,1]}`, together with that set's exact-enough volume
    `int_0^1 (t - w(u))^2 du` (trapezoid, float) -- the denominator that turns a credited fraction
    into a credited volume."""
    rng = rng or np.random.default_rng(20260912)
    tf = float(t)
    u = rng.uniform(0.0, 1.0, n)
    th = 2 * np.arctan(u)
    w = np.cos(th) + np.sin(th)
    lo, hi = w / 2, tf - w / 2
    P = np.column_stack([rng.uniform(0, 1, n) * (hi - lo) + lo,
                         rng.uniform(0, 1, n) * (hi - lo) + lo, th])
    ug = np.linspace(0.0, 1.0, 200001)
    tg = 2 * np.arctan(ug)
    vol = float(np.trapezoid((tf - np.cos(tg) - np.sin(tg)) ** 2, ug))
    # the sample is uniform on u but NOT on (x, y) area, so weight each draw by its slice area
    return P, vol, (hi - lo) ** 2


def cmd_report(a):
    B = json.load(open(os.path.join(RUNS, f'am_{a.TAG}_boxes.json')))
    hp = os.path.join(RUNS, f'am_{a.TAG}_honest.json')
    H = json.load(open(hp)) if os.path.exists(hp) else None
    R = B['rows']
    zt = B['clique_total']
    zs = sum(r['z'] for r in R)
    wb = [r for r in R if r['vol'] > 0]
    print(f"## {B['tag']}: Theta = {B['theta']:.9f}, clique dual {zt:.6f}, "
          f"{B['n_rows']} rows examined ({zs:.6f} of it)")
    print(f"   rows with vol(B_K) > 0 : {len(wb)} / {len(R)}, dual {sum(r['z'] for r in wb):.6f} "
          f"= {100 * sum(r['z'] for r in wb) / zt:.2f} % of the clique dual")
    print(f"   rows with vol(B_K) = 0 : {len(R) - len(wb)}, dual {zs - sum(r['z'] for r in wb):.6f}"
          f" = {100 * (zs - sum(r['z'] for r in wb)) / zt:.2f} %")
    if a.measure:
        d = np.load(os.path.join(RUNS, f'hc_{a.TAG}_dual.npz'), allow_pickle=True)
        cov = hc.Cover(d, 1e-9, log=lambda m: None)
        S, vol, wgt = pose_volume_sample(a.measure)
        print(f"   admissible pose volume (cx, cy, u) = {vol:.6f}; estimating vol(A(K)) from "
              f"{a.measure} uniform poses")
        for r in R[:a.top]:
            sub = hc.one_row_cover(cov, r['row'])
            cred = sub.clique_capture(S, 'meet') > 0
            vA = float((cred * wgt).sum() / wgt.sum() * vol)
            r['volA'] = vA
        print()
    if a.graze:
        d3 = np.load(os.path.join(RUNS, f"hc_{B['tag']}_dual.npz"), allow_pickle=True)
        cov3 = hc.Cover(d3, 1e-9, log=lambda m: None)
        print('   pairwise meet margins inside the top rows (there is no grazing at 1e-6: these '
              'rows are ROBUSTLY pairwise overlapping and non-Helly)')
        print('| rank | row | \\|K\\| | min pair margin | 1 % | median | pairs < 1e-6 / 1e-4 / '
              '1e-3 / 1e-2 |')
        print('|---|---|---|---|---|---|---|')
        for r in R[:a.top]:
            g = grazing_anatomy(cov3.F[cov3.cmem[r['row']]], cap=a.cap)
            r['graze'] = g
            q = g['pair_quantiles']
            print(f"| {r['rank']} | {r['row']} | {r['size']} | {q[0]:+.6f} | {q[1]:+.6f} | "
                  f"{q[3]:+.4f} | " + ' / '.join(str(v) for v in g['n_below']) + f" of "
                  f"{g['n_pairs']} |")
        print()
    if a.verify:
        PE = [fold_pose(q) for q in read_exact_poses(B['tag'], a.base)]
        d2 = np.load(os.path.join(RUNS, f"hc_{B['tag']}_dual.npz"), allow_pickle=True)
        cov2 = hc.Cover(d2, 1e-9, log=lambda m: None)
        nv, nbad, vmax = 0, 0, 0.0
        for r in R:
            if not r.get('box'):
                continue
            box = tuple(Fr(v) for v in r['box'][0])
            ok, why = certify_box(box, [PE[int(k)] for k in cov2.cmem[r['row']]])
            nv += 1
            vmax = max(vmax, box_vol(box))
            if not ok:
                nbad += 1
                print(f"   !! row {r['row']}: {why}")
        print(f"   EXACT RE-VERIFICATION: {nv} boxes re-certified from the json "
              f"({nbad} failures); largest vol(B_K) = {vmax:.6f}")
    mm = np.array([r['maxmin_margin'] for r in R])
    print(f"   max-min meet margin over A(K): max {mm.max():+.3e}, median {np.median(mm):+.3e}, "
          f"min {mm.min():+.3e}; {int((mm > 1e-12).sum())} rows with a strictly positive one")
    print()
    print('| rank | row | z_K | \\|K\\| | core deficit | grazing pairs | max graze degree | '
          'max-min margin | vol(B_K) | min margin on B_K |'
          + (' vol(A(K)) |' if a.measure else ''))
    print('|---|---|---|---|---|---|---|---|---|---|' + ('---|' if a.measure else ''))
    for r in R[:a.top]:
        g = r.get('graze', {})
        print(f"| {r['rank']} | {r['row']} | {r['z']:.6f} | {r['size']} | "
              f"{r.get('core_deficit', float('nan')):+.6f} | "
              f"{g.get('n_grazing', 0)} / {g.get('n_pairs', 0)} "
              f"({100 * g.get('graze_frac', 0):.1f} %) | "
              f"{g.get('max_degree', 0)} of {g.get('n_sampled', 0)} | "
              f"{r['maxmin_margin']:+.4f} | "
              f"{r['vol']:.3e} | "
              f"{r.get('box_min_margin', float('nan')):+.3e} |"
              + (f" {r['volA']:.3e} |" if a.measure else ''))
    if H:
        print()
        print(f"### honest cost under `credit z_K iff S in K u B_K` (sound)")
        for k in ('poses_of_P', 'lattice', 'families', 'descent', 'knife_edge'):
            if k in H:
                print(f"   {k:12s} n={H[k]['n']:>9d}  min capture {H[k]['min']:.9f}  "
                      f"honest {H[k]['honest']:.6f}")
        print(f"   MIN capture {H['min_capture']:.9f}  ->  HONEST = {H['honest']:.6f}")
        m = H['minimiser']
        print(f"   at ({m['pose'][0]:.6f}, {m['pose'][1]:.6f}, {m['theta_deg']:.4f} deg): "
              f"points {m['points']:.6f} ({m['n_points']}) + cliques {m['cliques']:.6f} "
              f"({m['n_cliques']}) + polygons {m['pgons']:.6f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    c = sub.add_parser('boxes')
    c.add_argument('TAG')
    c.add_argument('--base', default=RUNS)
    c.add_argument('--geo-tol', type=float, default=1e-9)
    c.add_argument('--nrows', type=int, default=10 ** 9)
    c.add_argument('--nth', type=int, default=181)
    c.add_argument('--sub0', type=int, default=64)
    c.add_argument('--mtol', type=float, default=1e-12,
                   help='a max-min margin above this is taken as "A(K) has interior"')
    c.add_argument('--den', type=int, default=10 ** 9)
    c.add_argument('--pad', type=float, default=1e-6,
                   help='retry distance from the admissibility wall when the first grow fails')
    c.add_argument('--gtol', type=float, default=1e-6)
    c.add_argument('--cap', type=int, default=600)
    c.add_argument('--nanat', type=int, default=20)

    c = sub.add_parser('honest')
    c.add_argument('TAG')
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--geo-tol', type=float, default=1e-9)
    c.add_argument('--pitch', type=float, default=0.04)
    c.add_argument('--dth', type=float, default=2.5)
    c.add_argument('--families', default=os.path.join(RUNS, 'hc_families.txt'))
    c.add_argument('--nseed', type=int, default=400)
    c.add_argument('--ndescend', type=int, default=400)
    c.add_argument('--nfine', type=int, default=100)
    c.add_argument('--base', default=RUNS)
    c.add_argument('--refit', type=int, default=0,
                   help='rounds of re-growing every box around the current minimiser')

    c = sub.add_parser('report')
    c.add_argument('TAG')
    c.add_argument('--top', type=int, default=10)
    c.add_argument('--measure', type=int, default=0,
                   help='estimate vol(A(K)) for the listed rows from this many uniform poses')
    c.add_argument('--base', default=RUNS)
    c.add_argument('--verify', action='store_true', default=True,
                   help='re-run the exact certificate on every box in the json (default on)')
    c.add_argument('--graze', action='store_true', default=True,
                   help='recompute the pairwise-margin distribution of the top rows (default on)')
    c.add_argument('--no-graze', dest='graze', action='store_false')
    c.add_argument('--cap', type=int, default=400)
    c.add_argument('--no-verify', dest='verify', action='store_false')

    c = sub.add_parser('helly')
    c.add_argument('--d', type=float, default=0.02)
    c.add_argument('--den', type=int, default=10 ** 6)
    c.add_argument('--half', type=float, default=0.02, help='box half-side in cx, cy')
    c.add_argument('--halfu', type=float, default=0.02, help='box half-side in u')

    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    {'boxes': cmd_boxes, 'honest': cmd_honest, 'report': cmd_report, 'helly': cmd_helly}[a.cmd](a)


if __name__ == '__main__':
    main()
