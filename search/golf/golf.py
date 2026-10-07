#!/usr/bin/env python3
"""Golf transforms for D4-invariant mixed covers (mixed v1, points + axis-parallel segments).

Every transform returns an exactly D4-invariant cover (checked), all integers, masses rounded UP; nothing
here is trusted: the output is a candidate to be checked by zmx2 / zm_mixed.py.

  info   COVER                         counts, total, weight quantiles
  scale  COVER OUT (--factor F | --total T)   multiply every mass by F (or by T/total), rounding up
  merge  COVER OUT --h H [--wall W]     cluster the points by D4 orbit: orbit representatives in the fundamental
                                        triangle 0 <= x <= y <= s/2 are bucketed on a grid of pitch H (anchored at 0,
                                        so grid lines are bucket boundaries); each bucket becomes ONE orbit at the
                                        mass-weighted centroid (rounded to the 1/D grid) carrying the bucket's total
                                        orbit mass (8 images of M/8 each, coincident images summed)
                                        --wall W: leave points within W of the container walls unmerged
  coarsen COVER OUT --k K               merge K consecutive segment pieces on each line (blocks aligned at 0;
                                        K*len must divide the line length's symmetric structure): one segment
                                        spanning the block's present pieces, mass = the block's total
  drop   COVER OUT --min-w W            delete points of mass < W/W_den (whole orbits)
"""
import argparse
import os
import sys
from collections import defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import mixed_cover as MC  # noqa: E402


def side(cv):
    S = cv['s_num'] * cv['D']
    assert S % cv['s_den'] == 0
    return S // cv['s_den']


def images(S, X, Y):
    out = []
    for (a, b) in ((X, Y), (Y, X)):
        for (c, d) in ((a, b), (S - a, b), (a, S - b), (S - a, S - b)):
            out.append((c, d))
    return out


def canon(S, X, Y):
    a, b = min(X, S - X), min(Y, S - Y)
    return (min(a, b), max(a, b))


def seg_norm(e):
    X0, Y0, X1, Y1, w = e
    if (X1, Y1) < (X0, Y0):
        X0, Y0, X1, Y1 = X1, Y1, X0, Y0
    return (X0, Y0, X1, Y1, w)


def seg_images(S, e):
    X0, Y0, X1, Y1, _ = e
    out = []
    for (p, q) in zip(images(S, X0, Y0), images(S, X1, Y1)):
        a, b = (p, q) if p <= q else (q, p)
        out.append(a + b)
    return out


def d4_ok(cv):
    S = side(cv)
    P = defaultdict(int)
    for X, Y, w in cv['points']:
        P[(X, Y)] += w
    for (X, Y), w in P.items():
        for im in images(S, X, Y):
            if P.get(im) != w:
                return False
    G = defaultdict(int)
    for e in cv['segments']:
        e = seg_norm(e)
        G[e[:4]] += e[4]
    for k, w in G.items():
        for im in seg_images(S, k + (0,)):
            if G.get(im) != w:
                return False
    return True


def finish(cv, out, comment):
    # aggregate coincident points / segments, drop zeros, sort
    P = defaultdict(int)
    for X, Y, w in cv['points']:
        P[(X, Y)] += w
    cv['points'] = sorted((X, Y, w) for (X, Y), w in P.items() if w > 0)
    G = defaultdict(int)
    for e in cv['segments']:
        e = seg_norm(e)
        G[e[:4]] += e[4]
    cv['segments'] = sorted(k + (w,) for k, w in G.items() if w > 0)
    MC.validate(cv)
    assert d4_ok(cv), "not D4-invariant"
    MC.write(out, cv, comment)
    t = MC.total(cv)
    print(f"{out}: points {len(cv['points'])} segments {len(cv['segments'])} total {t} = {float(t):.9f}")
    return cv


def ceil_div(a, b):
    return -((-a) // b)


def cmd_info(a):
    cv = MC.load(a.cover)
    ws = sorted((w for _, _, w in cv['points']), reverse=True)
    W = cv['W']
    tp = sum(ws) / W
    tsg = sum(e[4] for e in cv['segments']) / W
    print(f"points {len(ws)} (mass {tp:.6f}), segments {len(cv['segments'])} (mass {tsg:.6f}), "
          f"total {float(MC.total(cv)):.9f}, D4 {d4_ok(cv)}")
    for q in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9):
        print(f"  weight quantile {q}: {ws[int(q * (len(ws) - 1))] / W:.3e}")
    for thr in (1e-3, 1e-4, 1e-5):
        n = sum(1 for w in ws if w >= thr * W)
        print(f"  points >= {thr:g}: {n}, mass {sum(w for w in ws if w >= thr * W) / W:.5f}")


def cmd_scale(a):
    cv = MC.load(a.cover)
    if a.total is not None:
        f = Fraction(a.total) / MC.total(cv)
    else:
        f = Fraction(a.factor)
    f = f.limit_denominator(10 ** 9) if a.total is None else Fraction(int(f * 10 ** 9), 10 ** 9)
    cv['points'] = [(x, y, ceil_div(w * f.numerator, f.denominator)) for x, y, w in cv['points']]
    cv['segments'] = [(p, q, r, s_, ceil_div(w * f.numerator, f.denominator)) for p, q, r, s_, w in cv['segments']]
    finish(cv, a.out, f"scaled by {f} = {float(f):.9f} from {os.path.basename(a.cover)} (search/golf/golf.py)")


def cmd_merge(a):
    cv = MC.load(a.cover)
    S = side(cv)
    D = cv['D']
    h = int(round(float(Fraction(a.h)) * D))
    wall = int(round(float(Fraction(a.wall)) * D)) if a.wall else 0
    # orbit representatives with their orbit mass
    orb = defaultdict(int)
    for X, Y, w in cv['points']:
        orb[canon(S, X, Y)] += w          # sums w over the orbit: orbit mass
    keep, buckets = {}, defaultdict(list)
    for (x, y), m in orb.items():
        if wall and x < wall:              # x = min coordinate distance to a wall
            keep[(x, y)] = m
            continue
        buckets[(x // h, y // h)].append((x, y, m))
    newpts = defaultdict(int)

    def put(x, y, M):
        ims = images(S, x, y)
        each = ceil_div(M, 8)
        for im in ims:
            newpts[im] += each
    for (x, y), m in keep.items():
        # re-emit the orbit exactly as it was (orbit mass m over its distinct images)
        ims = set(images(S, x, y))
        for im in ims:
            newpts[im] += m // len(ims)
            assert m % len(ims) == 0
    for key, lst in buckets.items():
        M = sum(m for _, _, m in lst)
        cx = sum(x * m for x, _, m in lst) / M
        cy = sum(y * m for _, y, m in lst) / M
        put(int(round(cx)), int(round(cy)), M)
    cv['points'] = [(X, Y, w) for (X, Y), w in newpts.items()]
    finish(cv, a.out, f"points merged by D4 orbit on a {a.h} grid (wall band {a.wall or 0} kept) from "
           f"{os.path.basename(a.cover)} (search/golf/golf.py)")


def cmd_coarsen(a):
    cv = MC.load(a.cover)
    S = side(cv)
    lines = defaultdict(list)
    for e in cv['segments']:
        X0, Y0, X1, Y1, w = seg_norm(e)
        if X0 == X1:
            lines[(0, X0)].append((Y0, Y1, w))
        else:
            assert Y0 == Y1
            lines[(1, Y0)].append((X0, X1, w))
    lens = {q[1] - q[0] for v in lines.values() for q in v}
    assert len(lens) == 1, lens
    L = lens.pop() * a.k
    assert S % L == 0 or True
    out = []
    for (d, K), pcs in lines.items():
        blk = defaultdict(list)
        for lo, hi, w in pcs:
            # blocks symmetric under t -> S - t: anchor at 0 if S/L is an integer, else at the centre
            if S % L == 0:
                b = lo // L
            else:
                c = S // 2
                b = (lo - c) // L if lo >= c else -((c - hi) // L) - 1
            blk[b].append((lo, hi, w))
        for b, lst in blk.items():
            lo = min(q[0] for q in lst)
            hi = max(q[1] for q in lst)
            w = sum(q[2] for q in lst)
            out.append((K, lo, K, hi, w) if d == 0 else (lo, K, hi, K, w))
    cv['segments'] = out
    finish(cv, a.out, f"segment pieces merged {a.k} to 1 from {os.path.basename(a.cover)} (search/golf/golf.py)")


def cmd_drop(a):
    cv = MC.load(a.cover)
    thr = Fraction(a.min_w) * cv['W']
    cv['points'] = [(X, Y, w) for X, Y, w in cv['points'] if w >= thr]
    finish(cv, a.out, f"points of mass < {a.min_w} dropped from {os.path.basename(a.cover)} (search/golf/golf.py)")


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('info'); p.add_argument('cover')
    p = sp.add_parser('scale'); p.add_argument('cover'); p.add_argument('out')
    p.add_argument('--factor'); p.add_argument('--total')
    p = sp.add_parser('merge'); p.add_argument('cover'); p.add_argument('out'); p.add_argument('--h', required=True)
    p.add_argument('--wall', default=None)
    p = sp.add_parser('coarsen'); p.add_argument('cover'); p.add_argument('out'); p.add_argument('--k', type=int, required=True)
    p = sp.add_parser('drop'); p.add_argument('cover'); p.add_argument('out'); p.add_argument('--min-w', required=True)
    a = ap.parse_args()
    dict(info=cmd_info, scale=cmd_scale, merge=cmd_merge, coarsen=cmd_coarsen, drop=cmd_drop)[a.cmd](a)


if __name__ == '__main__':
    main()
