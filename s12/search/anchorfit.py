#!/usr/bin/env python3
"""anchorfit.py -- fit a CERTIFIABLE anchor-graph family to an exact window cut.

`search/locality.py check` produces, for every window with `z_j > 1`, the exact optimal dual `pi`
of the window subproblem: the local facet of the packing polytope that the measure violates, as a
weight per support pose.  That object is exact but it is not certifiable -- a verifier cannot check
"this list of 33 poses has independence number 2".  The certifiable form is the ANCHOR GRAPH.

**The family.**  Let `H` be a graph, and to each vertex `v` attach an ANCHOR `A_v` -- a point or a
closed segment of the container -- with

    X_v = { S admissible : A_v subset S }          (S convex: a segment is inside iff both ends are)

and let `pi_v >= 0` be weights.  Assign every square to the LEAST `v` with `A_v subset S`, so the
pieces are pairwise disjoint by a rule a verifier can apply.  Then for every packing `P`:

  * two squares of `P` cannot both lie in `X_v` (both contain `A_v`, so they meet);
  * if `uv` is an edge of `H`, certified by `A_u cap A_v != {}`, a square of `P` in `X_u` and one in
    `X_v` both contain a common point, so they meet -- and being in disjoint pieces they are
    distinct, which is impossible.

So `S -> (its piece)` injects `P cap (union X_v)` into an INDEPENDENT set of `H`, each square
paying its own piece's weight, and

    sum_v pi_v mu(X_v)  <=  alpha_pi(H) := max { sum_{v in I} pi_v : I independent in H }

is valid for every packing of disjoint closed unit squares -- with no geometry beyond "a convex set
contains a segment iff it contains both endpoints" and the exact intersection test on `H`'s edges.
`H = C_k`, `pi = 1` is `search/rankfamily.py`'s odd polygon (`alpha_pi = (k-1)/2`); `H = K join C_k`
with `pi = 1` on the clique and `1/2` on the cycle is the odd wheel, `alpha_pi = 1` for `k = 5`.

**What this tool does.**  Given a cut, it searches for anchors whose pieces cover the cut's heavy
poses, builds `H` from the EXACT anchor intersections (so `H` is whatever the geometry gives -- a
non-edge is never assumed), computes `alpha_pi(H)` exactly by brute force over the independent sets
of a graph with at most 8 vertices, extends the pieces over the WHOLE loaded pose set (maximal rows,
as `rankfamily.row_members` does), and reports how much of the cut's violation the certified family
captures:

    captured = (sum_v pi_v mu(X_v) - alpha_pi(H)) / (pi . mu - 1)

`1.0` means the anchor family is as strong as the exact local facet there; `0` means the cut has no
description in this family.

    python3 search/anchorfit.py fit TAG --check runs/loc_QA2_check.json \\
        --poses runs/cl_A0101L2_poses.txt --measure runs/cl_A0101L2_exact.txt
    python3 search/anchorfit.py selftest
"""
import argparse
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import leaf_ceiling as lc                                            # noqa: E402
import rankfamily as rf                                              # noqa: E402
from cliquelever import Poses                                        # noqa: E402

RUNS = os.path.join(REPO, 'runs')
DM = 10 ** 9


# ====================================================================== anchors
def seg_meets_seg(a, b):
    """exact: do two closed anchors meet?  An anchor is (P, Q) with P, Q reduced integer triples
    (X, Y, D); P == Q is a point anchor.  Segment/segment intersection over the rationals."""
    (p0, p1), (q0, q1) = a, b
    if p0 == p1 and q0 == q1:
        return p0 == q0
    e1 = (p0[0] * p1[2], p0[1] * p1[2], p1[0] * p0[2], p1[1] * p0[2], p0[2] * p1[2])
    e2 = (q0[0] * q1[2], q0[1] * q1[2], q1[0] * q0[2], q1[1] * q0[2], q0[2] * q1[2])
    if p0 == p1:
        return _point_on_seg(p0, e2)
    if q0 == q1:
        return _point_on_seg(q0, e1)
    if lc.seg_intersect(e1, e2) is not None:
        return True
    # collinear overlap: `seg_intersect` returns None when the determinant vanishes
    return (_point_on_seg(p0, e2) or _point_on_seg(p1, e2)
            or _point_on_seg(q0, e1) or _point_on_seg(q1, e1))


def _point_on_seg(P, e):
    """exact: is the point P = (X, Y, D) on the closed segment e = (x1, y1, x2, y2, De)?"""
    X, Y, D = P
    x1, y1, x2, y2, De = e
    ax, ay = X * De, Y * De
    cross = (x2 - x1) * (ay - y1 * D) - (y2 - y1) * (ax - x1 * D)
    if cross != 0:
        return False
    dot = (ax - x1 * D) * (x2 - x1) + (ay - y1 * D) * (y2 - y1)
    if dot < 0:
        return False
    L = (x2 - x1) ** 2 + (y2 - y1) ** 2
    return dot <= L * D


def piece_members(ps, anchor):
    """boolean over ALL loaded poses: the square contains the whole anchor (both endpoints).
    Float test with the exact integer test inside the `rankfamily.FTOL` band."""
    P, Q = anchor
    m = rf.contains_vec(ps, *P)
    if Q != P:
        m = m & rf.contains_vec(ps, *Q)
    return m


def alpha_pi(H, pi):
    """exact weighted independence number of a small graph: max sum of pi over independent sets"""
    n = len(pi)
    best = 0.0
    for r in range(n + 1):
        for I in combinations(range(n), r):
            if all((u, v) not in H and (v, u) not in H for u, v in combinations(I, 2)):
                best = max(best, sum(pi[v] for v in I))
    return best


def build_H(anchors):
    """the edge set certified by EXACT anchor intersection (a non-edge is never assumed: if two
    anchors do meet, the edge is there, which can only make alpha_pi smaller)"""
    E = set()
    for u, v in combinations(range(len(anchors)), 2):
        if seg_meets_seg(anchors[u], anchors[v]):
            E.add((u, v))
    return E


def evaluate(ps, anchors, pi, mu):
    """assign every loaded pose to the LEAST piece containing it, then value the family.
    Returns (value, alpha, H, sizes, masses)."""
    mem = [piece_members(ps, A) for A in anchors]
    taken = np.zeros(ps.n, dtype=bool)
    sizes, masses = [], []
    for k in range(len(anchors)):
        own = mem[k] & ~taken
        taken |= own
        sizes.append(int(own.sum()))
        masses.append(float(mu[own].sum()))
    H = build_H(anchors)
    a = alpha_pi(H, pi)
    val = sum(pi[k] * masses[k] for k in range(len(anchors)))
    return val, a, H, sizes, masses


# ====================================================================== candidate anchors
def cut_candidates(ps, poses, w, ncand=90):
    """candidate anchor POINTS for a cut: the arrangement vertices of the cut's support, kept by
    the weight of the poses containing them (`rankfamily._cand_masks`'s rule), plus the pose
    centres.  Returns (M, pts) with `M[j, i]` = "cut pose i contains candidate j"."""
    squares = [ps.sq[i][:8] + (max(int(round(w[k] * 10 ** 9)), 1),) for k, i in enumerate(poses)]
    V = lc.enumerate_vertices(squares, ps.t, verbose=False, procs=1)
    Inc = lc.incidences(squares, V, verbose=False, procs=1)
    ww = np.array([s[8] for s in squares], dtype=np.int64)
    return rf._cand_masks(squares, V, Inc, ww, ncand)


def fit_cycle(ps, poses, w, mu, M, pts, k, rng, restarts=60, deadline=None):
    """`H = C_k`, `pi = 1`: the odd polygon.  `rankfamily._climb` maximises the cut-weighted mass
    of the union of the k pieces; the anchors are the k sides of the polygon."""
    wi = np.maximum(np.round(np.asarray(w) * 10 ** 9).astype(np.int64), 0)
    dl = deadline if deadline is not None else time.time() + 60
    found = rf._climb(M, wi, k, 0, rng, restarts, [], dl)
    out = []
    for _v, idx in found[:6]:
        anchors = [(pts[idx[i]], pts[idx[(i + 1) % k]]) for i in range(k)]
        pi = [1.0] * k
        out.append((f'C{k}', anchors, pi) + evaluate(ps, anchors, pi, mu))
    return out


def fit_wheel(ps, poses, w, mu, M, pts, k, rng, restarts=60, deadline=None):
    """`H = K1 join C_k` (the odd wheel) and its clique-hub generalisation: a POINT anchor for the
    hub, the k sides of a polygon for the rim, `pi = 1` on the hub and `1/2` on the rim.  The hub
    is the candidate point of heaviest cut mass among the poses the cycle does not already take."""
    wi = np.maximum(np.round(np.asarray(w) * 10 ** 9).astype(np.int64), 0)
    dl = deadline if deadline is not None else time.time() + 60
    found = rf._climb(M, wi, k, 0, rng, restarts, [], dl)
    out = []
    for _v, idx in found[:4]:
        rim = [(pts[idx[i]], pts[idx[(i + 1) % k]]) for i in range(k)]
        used = np.zeros(len(poses), dtype=bool)
        for i in range(k):
            used |= M[idx[i]] & M[idx[(i + 1) % k]]
        rest = wi * (~used)
        hub_scores = M @ rest
        for h in np.argsort(-hub_scores)[:3]:
            hub = (pts[int(h)], pts[int(h)])
            anchors = [hub] + rim
            pi = [1.0] + [0.5] * k
            out.append((f'K1+C{k}', anchors, pi) + evaluate(ps, anchors, pi, mu))
    return out


# ====================================================================== driver
def fit_cut(ps, mu, cut, a, log, tag=''):
    """every family in the menu, best first"""
    poses = [int(v) for v in cut['poses']]
    pi_c = [float(v) for v in cut['pi']]
    # the climb maximises the CONTRIBUTION a pose makes to the cut, pi_i * mu_i, not pi_i alone:
    # an anchor is worth having only where the measure actually sits
    w = [pi_c[k] * float(mu[i]) for k, i in enumerate(poses)]
    if len(poses) < 3 or max(w) <= 0:
        return None
    t0 = time.time()
    M, pts = cut_candidates(ps, poses, w, a.ncand)
    rng = random.Random(12345)
    cand = []
    dl = time.time() + a.fit_time
    for k in a.cycles:
        cand += fit_cycle(ps, poses, w, mu, M, pts, k, rng, a.restarts, dl)
        cand += fit_wheel(ps, poses, w, mu, M, pts, k, rng, a.restarts, dl)
    viol = cut['z'] - 1.0
    best = None
    rows = []
    for (name, anchors, pi, val, al, H, sizes, masses) in cand:
        exc = val - al
        frac = exc / viol if viol > 1e-12 else 0.0
        rec = dict(family=name, value=val, alpha=al, excess=exc, captured=frac,
                   edges=len(H), sizes=sizes, masses=[round(m, 6) for m in masses],
                   anchors=[[list(A[0]), list(A[1])] for A in anchors], pi=pi)
        rows.append(rec)
        if best is None or exc > best['excess']:
            best = rec
    rows.sort(key=lambda r: -r['excess'])
    # keep the best row of EVERY family as well as the overall top few, so the table can say what
    # each family was worth and not only which one won
    perfam = {}
    for r in rows:
        if r['family'] not in perfam:
            perfam[r['family']] = r
    log(f'   {tag}cut z = {cut["z"]:.9f} ({len(poses)} poses, pi values '
        f'{sorted(set(round(v, 4) for v in pi_c), reverse=True)}), violation {viol:.6f}')
    for r in list(perfam.values())[:a.show] if a.per_family else rows[:a.show]:
        log(f'     {r["family"]:8s} value {r["value"]:.6f} vs alpha_pi {r["alpha"]:.2f} '
            f'-> excess {r["excess"]:+.6f}, captured {100 * r["captured"]:.1f} %, '
            f'|E(H)| = {r["edges"]}, pieces {r["sizes"]}')
    return dict(z=cut['z'], violation=viol, poses=len(poses),
                pi_values=sorted(set(round(v, 6) for v in pi_c), reverse=True),
                best=best, by_family=perfam, all=rows[:a.show], secs=time.time() - t0)


def cmd_fit(a):
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'af_{a.TAG}.log'), 'a')

    def log(m):
        print(m, flush=True)
        logf.write(str(m) + '\n')
        logf.flush()

    ps = Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    ps.add_exact_file(a.poses)
    if os.path.abspath(a.measure) != os.path.abspath(a.poses):
        ps.add_exact_file(a.measure)
    ps.finish()
    t, sym, poses = lc.read_measure(a.measure)
    mu = np.zeros(ps.n)
    for (p_, q_, cx, cy, m) in poses:
        i = ps.key[(p_, q_, ps.bsnap(cx), ps.bsnap(cy))]
        mu[i] += float(m)
    log(f'# anchorfit {a.TAG}: {ps.n} poses, measure mass {mu.sum():.9f}')
    D = json.load(open(a.check))
    out = {}
    for key, blk in sorted(D.items()):
        for c in blk['cert']:
            an = c.get('anatomy')
            if not an or an['z'] <= 1 + 1e-9:
                continue
            r = fit_cut(ps, mu, an, a, log, tag=f'{key} window {c["window"]} ')
            if r:
                r['window'] = c['window']
                r['geo'] = c['geo']
                out.setdefault(key, []).append(r)
            if a.limit and sum(len(v) for v in out.values()) >= a.limit:
                break
    json.dump(out, open(os.path.join(RUNS, f'af_{a.TAG}.json'), 'w'), indent=1)
    log('# summary: cut -> fitted H -> captured fraction')
    for key, rs in sorted(out.items()):
        for r in rs:
            b = r['best']
            log(f'   {key} w{r["window"]:<3d} z = {r["z"]:.6f} viol {r["violation"]:.6f} '
                f'pi {r["pi_values"]} -> {b["family"]:8s} excess {b["excess"]:+.6f} '
                f'captured {100 * b["captured"]:.1f} %')


# ====================================================================== selftest
def selftest():
    ok = [0, 0]

    def chk(name, c, extra=''):
        ok[1] += 1
        ok[0] += bool(c)
        print(f'  {"ok  " if c else "FAIL"} {name} {extra}')

    P = lambda x, y, d=1: (x, y, d)                                   # noqa: E731
    chk('point anchors meet iff equal', seg_meets_seg((P(1, 1), P(1, 1)), (P(1, 1), P(1, 1)))
        and not seg_meets_seg((P(1, 1), P(1, 1)), (P(1, 2), P(1, 2))))
    chk('crossing segments meet', seg_meets_seg((P(0, 0), P(2, 2)), (P(0, 2), P(2, 0))))
    chk('disjoint segments do not', not seg_meets_seg((P(0, 0), P(1, 0)), (P(2, 0), P(3, 0))))
    chk('touching segments meet', seg_meets_seg((P(0, 0), P(1, 0)), (P(1, 0), P(2, 0))))
    chk('collinear overlap meets', seg_meets_seg((P(0, 0), P(2, 0)), (P(1, 0), P(3, 0))))
    chk('point on segment', seg_meets_seg((P(1, 0), P(1, 0)), (P(0, 0), P(2, 0))))
    # alpha_pi: C5 with pi = 1 is 2; the 5-wheel with hub 1 and rim 1/2 is 1
    C5 = {(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)}
    chk('alpha_pi(C5, 1) = 2', alpha_pi(C5, [1.0] * 5) == 2.0)
    W5 = set(C5) | {(0, v) for v in range(1, 6)}
    W5 = {(0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
          (1, 2), (2, 3), (3, 4), (4, 5), (1, 5)}
    chk('alpha_pi(W5, hub 1 rim 1/2) = 1', alpha_pi(W5, [1.0] + [0.5] * 5) == 1.0)
    C7 = {(i, (i + 1) % 7) for i in range(7)}
    chk('alpha_pi(C7, 1) = 3', alpha_pi(C7, [1.0] * 7) == 3.0)
    W7 = {(0, v) for v in range(1, 8)} | {(i, i % 7 + 1) for i in range(1, 8)}
    chk('alpha_pi(W7, hub 1 rim 1/2) = 3/2', alpha_pi(W7, [1.0] + [0.5] * 7) == 1.5)
    print(f'  {ok[0]}/{ok[1]}')
    return 0 if ok[0] == ok[1] else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('fit')
    c.add_argument('TAG')
    c.add_argument('--check', required=True, help='a loc_*_check.json with cut anatomies')
    c.add_argument('--poses', required=True)
    c.add_argument('--measure', required=True)
    c.add_argument('--t', default='4')
    c.add_argument('--r', default='1')
    c.add_argument('--bnd-delta', type=float, default=1e-6)
    c.add_argument('--cycles', type=int, nargs='+', default=[5, 7])
    c.add_argument('--ncand', type=int, default=90)
    c.add_argument('--restarts', type=int, default=60)
    c.add_argument('--fit-time', type=float, default=90.0)
    c.add_argument('--show', type=int, default=4)
    c.add_argument('--per-family', action='store_true',
                   help='log the best row of each family instead of the overall best few')
    c.add_argument('--limit', type=int, default=0)
    sub.add_parser('selftest')
    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    if a.cmd == 'fit':
        cmd_fit(a)
    else:
        sys.exit(selftest())


if __name__ == '__main__':
    main()
