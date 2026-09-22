#!/usr/bin/env python3
"""t3_chord_oblique: the chord row on lines of ARBITRARY direction.

General chord row (valid at delta >= 0).  Let L be any line and K a set of k squares
whose interiors meet L.  The open chords L n int(Q_i) are pairwise disjoint open
subintervals of L n [delta, T-delta]^2, and consecutive ones are separated by at least
delta.  Hence

    sum_{i in K} h_i(L)  <=  len( L n [delta, T-delta]^2 )  -  (k-1) delta
                         =   len( L n [0,T]^2 ) - kappa(L) delta - (k-1) delta,

kappa(L) > 0 the rate at which the container chord shrinks.  A certificate therefore
exists on L iff  sum_{i in K} h_i(L)  >=  len( L n [0,T]^2 ).  For a vertical or
horizontal L, kappa = 2 and this is the (C3) row of t3_chord_scan.

This script maximises  D(L) = sum_i h_i(L) - len(L n [0,T]^2)  over all directions and
offsets at a given configuration.  D(L) > 0 is impossible for a packing; D(L) = 0 with
k >= 2 is a tight certificate (delta <= 0); D(L) < 0 everywhere means no chord
certificate of any direction exists.

Commands
    pinwheel      the 3x3 grid minus a permutation: the delta = 0 packing at which the
                  axis-parallel chord clause fails
    packings      the same functional over the recorded margin-0 packings (slow)
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t3_chord_geom as G                                            # noqa: E402

T = 3.0


def poly_chord_len(poly, nu, s, tol=1e-12):
    """length of {z : nu.z = s} n poly, poly a convex polygon (list of vertices).
    Also returns whether the line meets the interior."""
    vals = [nu[0] * v[0] + nu[1] * v[1] - s for v in poly]
    pts = []
    m = len(poly)
    for i in range(m):
        a, b = vals[i], vals[(i + 1) % m]
        if abs(a) <= tol:
            pts.append(poly[i])
        if (a > tol and b < -tol) or (a < -tol and b > tol):
            u = a / (a - b)
            p, q = poly[i], poly[(i + 1) % m]
            pts.append((p[0] + u * (q[0] - p[0]), p[1] + u * (q[1] - p[1])))
    if len(pts) < 2:
        return 0.0, False
    d = max(math.hypot(p[0] - q[0], p[1] - q[1]) for p in pts for q in pts)
    interior = any(v > tol for v in vals) and any(v < -tol for v in vals)
    return d, interior


def sq_poly(x, y, th):
    return G.sq_vertices(x, y, th)


BOX = [(0.0, 0.0), (T, 0.0), (T, T), (0.0, T)]


def scan(z, n, nphi=721, ns=3000, verbose=False):
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    polys = [sq_poly(X[i], Y[i], TH[i]) for i in range(n)]
    best = (-1e9, None)
    for a in range(nphi):
        phi = a * math.pi / nphi
        nu = (math.cos(phi), math.sin(phi))
        lo = min(nu[0] * v[0] + nu[1] * v[1] for v in BOX)
        hi = max(nu[0] * v[0] + nu[1] * v[1] for v in BOX)
        for b in range(1, ns):
            s = lo + (hi - lo) * b / ns
            lb, _ = poly_chord_len(BOX, nu, s)
            if lb <= 1e-9:
                continue
            tot = 0.0
            k = 0
            for i in range(n):
                h, inter = poly_chord_len(polys[i], nu, s)
                if inter and h > 1e-12:
                    tot += h
                    k += 1
            if k >= 2:
                d = tot - lb
                if d > best[0]:
                    best = (d, (phi, s, k, tot, lb))
    return best


def cmd_pinwheel(argv):
    for name, cells in (('pinwheel (minus identity)', [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]),
                        ('minus transposition (01)', [(0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1)]),
                        ('minus 3-cycle', [(0, 0), (0, 1), (1, 1), (1, 2), (2, 0), (2, 2)]),
                        ('two columns of three', [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)])):
        z = [0.0]
        for (cx, cy) in cells:
            z += [cx + 0.5, cy + 0.5, 0.0]
        b = scan(z, 6)
        print('%-26s  max_L [sum h - len(L n box)] = %+.9f   (phi = %.2f deg, k = %d, '
              'sum h = %.6f, len = %.6f)'
              % (name, b[0], math.degrees(b[1][0]), b[1][2], b[1][3], b[1][4]))


def cmd_packings(argv):
    import t3_chord_scan as SC
    rows = SC.load(['runs/t3_exist_row1.jsonl', 'runs/t3_exist_hunt1.jsonl',
                    'runs/t3_exist_scope1.jsonl', 'runs/t3_chain_scan1.jsonl'])
    pk = [r for r in rows if r.get('delta', -1) >= -1e-9]
    lim = int(argv[0]) if argv else 40
    worst = (1e9, None)
    for r in pk[:lim]:
        b = scan(r['z'], 6, nphi=181, ns=800)
        if b[0] < worst[0]:
            worst = (b[0], r.get('fam'))
    print('over %d packings: min over packings of max_L D = %+.9f  (%s)'
          % (min(lim, len(pk)), worst[0], worst[1]))


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'pinwheel'
    globals()['cmd_' + cmd](sys.argv[2:])


if __name__ == '__main__':
    main()
