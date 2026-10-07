#!/usr/bin/env python3
"""t3_chord_geom: exact line-transversal (chord) geometry for unit squares, and the
leaf-linear rows it supplies.  Nothing in search/ is modified; s6skel is imported.

A *chord row* is a linear inequality in the centres that is valid only on a leaf.  The
three kinds used here:

  (C1) endpoint rows       alpha_i(L) = x_i + A_i,  beta_i(L) = x_i + B_i   (line L)
       with A_i, B_i piecewise-linear in the centre; on a leaf fixing which of the two
       slab constraints is active, alpha_i and beta_i are exactly linear in (x_i, y_i).
  (C2) chord-length rows   beta_i - alpha_i = h_i >= H  on a leaf "|f_i| <= d_0".
  (C3) chord-chain rows    for k squares meeting one line, in order,
       sum_r h_{i_r} <= T - (k+1) delta.

Commands
    selftest        all the checks of notes/t3-chord.md section 1
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6skel                                                        # noqa: E402


# --------------------------------------------------------------------------- basic chord
def sq_vertices(x, y, th):
    C, S = math.cos(th), math.sin(th)
    out = []
    for (a, b) in ((.5, .5), (-.5, .5), (-.5, -.5), (.5, -.5)):
        out.append((x + a * C - b * S, y + a * S + b * C))
    return out


def chord(x, y, th, a, axis):
    """exact chord of the unit square (x, y, th) on the line {axis = a}.

    axis = 'v': the vertical line x = a; the chord is an interval in y.
    axis = 'h': the horizontal line y = a; the chord is an interval in x.

    -> (alpha, beta, branch_lo, branch_hi) with alpha > beta meaning 'empty'.
    branch_* in {0, 1} says which of the two slab constraints is active (0 = the
    constraint from the square's own first edge normal, 1 = from the second); the
    endpoints are linear in (x, y) on the leaf that fixes the two branches.
    """
    C, S = math.cos(th), math.sin(th)
    if axis == 'v':
        e = a - x                       # offset along the line's normal
        base = y
        # cond1 |eC + n S| <= 1/2  -> n in [(-1/2 - eC)/S, (1/2 - eC)/S]
        # cond2 |n C - e S| <= 1/2 -> n in [(eS - 1/2)/C, (eS + 1/2)/C]
        p1, q1, p2, q2 = C, S, S, C
    else:
        e = a - y
        base = x
        # cond1 |x'C + eS| <= 1/2  -> x' in [(-1/2 - eS)/C, (1/2 - eS)/C]
        # cond2 |eC - x'S| <= 1/2  -> x' in [(eC - 1/2)/S, (eC + 1/2)/S]
        p1, q1, p2, q2 = S, C, C, S
    los, his = [], []
    # branch 0 : ( -+1/2 - e*p1 ) / q1     (empty constraint when q1 == 0)
    if abs(q1) > 1e-15:
        lo = (-0.5 - e * p1) / q1
        hi = (0.5 - e * p1) / q1
        los.append(min(lo, hi)); his.append(max(lo, hi))
    else:
        los.append(-math.inf); his.append(math.inf)
    if abs(q2) > 1e-15:
        lo = (e * p2 - 0.5) / q2
        hi = (e * p2 + 0.5) / q2
        los.append(min(lo, hi)); his.append(max(lo, hi))
    else:
        los.append(-math.inf); his.append(math.inf)
    bl = 0 if los[0] >= los[1] else 1
    bh = 0 if his[0] <= his[1] else 1
    return base + los[bl], base + his[bh], bl, bh


def chord_len_exact(x, y, th, a, axis):
    al, be, _, _ = chord(x, y, th, a, axis)
    return max(0.0, be - al)


def chord_mid(x, y, th, a, axis):
    al, be, _, _ = chord(x, y, th, a, axis)
    return 0.5 * (al + be)


# --------------------------------------------------------------------------- the band
def band_full(th):
    """(C - S)/2 for th in [0, 45 deg]: the half-width of the set of distances d for which
    the chord is at its maximum 1/max(|C|,|S|) (= sec(tilt)).  Beyond it the chord decays
    linearly with slope 1/(|C||S|)."""
    C, S = abs(math.cos(th)), abs(math.sin(th))
    return abs(C - S) / 2.0


def band_unit(th):
    """D(th) = p - |C||S|: the half-width of the set of d for which the chord is >= 1."""
    return s6skel.chord_D(th)


def chord_of_d(th, d):
    return s6skel.chord_len(th, d)


def _poly_chord(x, y, th, a, axis):
    """endpoints of the intersection of the square's boundary polygon with the line,
    by walking the four edges.  -> (lo, hi) or None if the line misses the square."""
    V = sq_vertices(x, y, th)
    k = 0 if axis == 'v' else 1
    o = 1 - k
    hits = []
    for i in range(4):
        p, q = V[i], V[(i + 1) % 4]
        if (p[k] - a) * (q[k] - a) <= 0 and abs(p[k] - q[k]) > 1e-14:
            s_ = (a - p[k]) / (q[k] - p[k])
            if -1e-12 <= s_ <= 1 + 1e-12:
                hits.append(p[o] + s_ * (q[o] - p[o]))
        elif abs(p[k] - a) < 1e-14:
            hits.append(p[o])
    if not hits:
        return None
    return (min(hits), max(hits))


# --------------------------------------------------------------------------- selftest
def _rand_in_square(rng, x, y, th, npts):
    C, S = math.cos(th), math.sin(th)
    uv = rng.uniform(-0.5, 0.5, size=(npts, 2))
    return np.stack([x + uv[:, 0] * C - uv[:, 1] * S,
                     y + uv[:, 0] * S + uv[:, 1] * C], axis=1)


def selftest():
    rng = np.random.default_rng(20260922)
    ok = [0]
    bad = [0]

    def chk(name, cond, extra=''):
        (ok if cond else bad).__getitem__(0)
        if cond:
            ok[0] += 1
            print('  ok   %-56s %s' % (name, extra))
        else:
            bad[0] += 1
            print('  FAIL %-56s %s' % (name, extra))

    # 1. chord() agrees with s6skel.chord_len on the distance to the line ---------------
    worst = 0.0
    for _ in range(4000):
        th = rng.uniform(0, math.pi / 2)
        x, y = rng.uniform(0.6, 2.4, 2)
        for axis in ('v', 'h'):
            a = rng.uniform(0, 3)
            d = abs(a - x) if axis == 'v' else abs(a - y)
            ref = s6skel.chord_len(th, d)
            got = chord_len_exact(x, y, th, a, axis)
            worst = max(worst, abs(ref - got))
    chk('chord length == s6skel.chord_len(th, d)', worst < 1e-9, 'max err %.2e' % worst)

    # 2. the chord really is the polygon-line intersection ----------------------------
    worst = 0.0
    for _ in range(4000):
        th = rng.uniform(0, math.pi / 2)
        x, y = rng.uniform(0.6, 2.4, 2)
        axis = 'v' if rng.random() < .5 else 'h'
        a = rng.uniform(x - .75, x + .75) if axis == 'v' else rng.uniform(y - .75, y + .75)
        al, be, _, _ = chord(x, y, th, a, axis)
        ref = _poly_chord(x, y, th, a, axis)
        if ref is None:
            worst = max(worst, (be - al) if be > al else 0.0)
        else:
            worst = max(worst, abs(al - ref[0]), abs(be - ref[1]))
    chk('chord == exact polygon-line intersection', worst < 1e-10, 'max err %.2e' % worst)

    # 3. THE BRIEF'S CLAIM: is the chord centred at the centre? ------------------------
    worstmid = 0.0
    arg = None
    for _ in range(20000):
        th = rng.uniform(0, math.pi / 2)
        x, y = 1.5, 1.5
        axis = 'v' if rng.random() < .5 else 'h'
        a = rng.uniform(-0.71, 0.71) + (x if axis == 'v' else y)
        al, be, _, _ = chord(x, y, th, a, axis)
        if be - al <= 1e-9:
            continue
        c0 = y if axis == 'v' else x
        off = abs(0.5 * (al + be) - c0)
        if off > worstmid:
            worstmid, arg = off, (math.degrees(th), a - (x if axis == 'v' else y))
    chk('brief premise REFUTED: chord midpoint != centre coordinate', worstmid > 1e-3,
        'max offset %.6f at th=%.2f deg, e=%.4f' % (worstmid, arg[0], arg[1]))

    # 4. ... and the *correct* statement: midpoints, not centres, are what separate -----
    #    counterexample search for  |y_i - y_j| >= (h_i + h_j)/2  at a genuine packing
    ce = _midpoint_counterexample()
    chk('brief premise REFUTED: |dy| >= (h_i+h_j)/2 is FALSE', ce is not None,
        'counterexample: %s' % (ce,) if ce else 'none found')

    # 5. leaf-linearity: on a fixed branch pair the endpoints are affine in the centre --
    worst = 0.0
    for _ in range(3000):
        th = rng.uniform(0.05, math.pi / 2 - 0.05)
        x, y = rng.uniform(0.6, 2.4, 2)
        axis = 'v' if rng.random() < .5 else 'h'
        a = rng.uniform(x - .4, x + .4) if axis == 'v' else rng.uniform(y - .4, y + .4)
        al0, be0, bl0, bh0 = chord(x, y, th, a, axis)
        if be0 - al0 < 1e-3:
            continue
        for _ in range(4):
            dx, dy = rng.normal(0, 3e-4, 2)
            al1, be1, bl1, bh1 = chord(x + dx, y + dy, th, a, axis)
            if (bl1, bh1) != (bl0, bh0):
                continue
            # predicted from the exact affine forms
            pa = _affine_endpoints(th, a, axis, bl0, bh0)
            p0 = pa[0][0] * x + pa[0][1] * y + pa[0][2]
            p1 = pa[0][0] * (x + dx) + pa[0][1] * (y + dy) + pa[0][2]
            q0 = pa[1][0] * x + pa[1][1] * y + pa[1][2]
            q1 = pa[1][0] * (x + dx) + pa[1][1] * (y + dy) + pa[1][2]
            worst = max(worst, abs(p0 - al0), abs(p1 - al1), abs(q0 - be0), abs(q1 - be1))
    chk('endpoints affine in (x, y) on a fixed branch leaf', worst < 1e-10,
        'max err %.2e' % worst)

    # 6. the band identities ------------------------------------------------------------
    w1 = w2 = 0.0
    for k in range(2001):
        th = k * (math.pi / 4) / 2000
        w1 = max(w1, abs(chord_of_d(th, band_full(th)) - 1.0 / max(math.cos(th), math.sin(th))))
        w2 = max(w2, abs(chord_of_d(th, band_unit(th)) - 1.0))
    chk('chord((C-S)/2) = sec(tilt)', w1 < 1e-9, 'max err %.2e' % w1)
    chk('chord(D(th))   = 1', w2 < 1e-9, 'max err %.2e' % w2)

    print('\n  %d ok, %d FAIL' % (ok[0], bad[0]))
    return bad[0]


def _affine_endpoints(th, a, axis, bl, bh):
    """(alpha, beta) as affine forms (cx, cy, const) in the centre, on the leaf (bl, bh)."""
    C, S = math.cos(th), math.sin(th)
    if axis == 'v':
        # e = a - x;  alpha = y + lo_branch(e),  beta = y + hi_branch(e)
        #   branch 0 :  eta in [(-1/2 - eC)/S, (1/2 - eC)/S]
        #   branch 1 :  eta in [( eS - 1/2)/C, ( eS + 1/2)/C]
        lo = (C / S, 1.0, (-0.5 - a * C) / S) if bl == 0 else (-S / C, 1.0, (a * S - 0.5) / C)
        hi = (C / S, 1.0, (0.5 - a * C) / S) if bh == 0 else (-S / C, 1.0, (a * S + 0.5) / C)
        return lo, hi
    else:
        if bl == 0:
            lo = (1.0, S / C, (-0.5 - a * S) / C)
        else:
            lo = (1.0, -C / S, (a * C - 0.5) / S)
        if bh == 0:
            hi = (1.0, S / C, (0.5 - a * S) / C)
        else:
            hi = (1.0, -C / S, (a * C + 0.5) / S)
        return lo, hi


def _midpoint_counterexample():
    """two disjoint unit squares both meeting one line with |dy| < (h_i + h_j)/2."""
    rng = np.random.default_rng(7)
    best = None
    for _ in range(200000):
        th1 = rng.uniform(0, math.pi / 2)
        th2 = rng.uniform(0, math.pi / 2)
        x1, y1 = 1.5, 1.5
        x2 = x1 + rng.uniform(-1.4, 1.4)
        y2 = y1 + rng.uniform(-1.4, 1.4)
        z = [0.0, x1, y1, th1, x2, y2, th2]
        if s6skel.pair_gap(z, 0, 1) < 0:
            continue
        a = rng.uniform(min(x1, x2) - .7, max(x1, x2) + .7)
        h1 = chord_len_exact(x1, y1, th1, a, 'v')
        h2 = chord_len_exact(x2, y2, th2, a, 'v')
        if h1 <= 1e-6 or h2 <= 1e-6:
            continue
        v = abs(y2 - y1) - 0.5 * (h1 + h2)
        if best is None or v < best[0]:
            best = (v, dict(th1=round(math.degrees(th1), 3), th2=round(math.degrees(th2), 3),
                            c1=(round(x1, 4), round(y1, 4)), c2=(round(x2, 4), round(y2, 4)),
                            line=round(a, 4), h1=round(h1, 4), h2=round(h2, 4),
                            dy=round(abs(y2 - y1), 4), violation=round(v, 5)))
    if best is not None and best[0] < -1e-6:
        return best[1]
    return None


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'selftest'
    if cmd == 'selftest':
        sys.exit(1 if selftest() else 0)
    print(__doc__)
