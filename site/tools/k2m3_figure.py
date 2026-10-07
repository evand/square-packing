#!/usr/bin/env python3
"""Inline SVG of the s(k^2-3) fixed-profile family, for docs/k2m3.html.

Reads the family file (search/qx2_data/L4_k02_family.txt: profile pi, corner module nu, R = w = 2,
Lebesgue from a = 9/5), builds the box measure mu_k exactly as QUADRANT_EXACT.md section 2 and
notes/lean-bentz-reduction.md describe it (nu at the four corners, pi on [R, k-R] along each wall with
the closed end x = k-R a phase-0 cross-section, Lebesgue on [a, k-a]^2), and prints an SVG drawing.

Self-check: for k = 7 the segments must equal those of L4_k02_box7.txt (same multiset of unit
1/5-grid segments and masses), and the total must be k^2 - 4D.  The script refuses to print otherwise.

Usage: python3 site/tools/k2m3_figure.py [k] > /tmp/fig.svg   (default k = 9)
"""
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FAM = ROOT / 'search/qx2_data/L4_k02_family.txt'
BOX7 = ROOT / 'search/qx2_data/L4_k02_box7.txt'
R, A = 2, F(9, 5)
D = F(423621306389, 500000000000)


def read_family():
    prof, corner = [], []
    for line in FAM.read_text().splitlines():
        t = line.split('#')[0].split()
        if not t:
            continue
        kind, vals = t[0], [F(x) for x in t[1:]]
        (prof if kind in 'hv' else corner).append((kind, *vals))
    return prof, corner


def box_segments(k):
    """Segments of mu_k as (x0, y0, x1, y1, mass), coordinates as Fractions."""
    prof, corner = read_family()
    base = []                                   # bottom-left corner module + bottom band, before symmetry
    for kind, a, b, c, m in corner:
        base.append((b, a, c, a, m) if kind == 'H' else (a, b, a, c, m))
    for kind, a, b, c, m in prof:
        for n in range(R, k - R):              # half-open periods [n, n+1)
            base.append((n + b, a, n + c, a, m) if kind == 'h' else (n + a, b, n + a, c, m))
        if kind == 'v' and a == 0:              # closed end x = k - R (phase 0)
            base.append((k - R, b, k - R, c, m))
    segs = []
    # corners: all four images of nu; bands: bottom band and its images under the four
    # symmetries taking the bottom wall to each wall (nu is diagonal symmetric, pi mirror symmetric).
    ncorner = len(corner)
    maps_corner = [lambda x, y: (x, y), lambda x, y: (k - x, y), lambda x, y: (x, k - y), lambda x, y: (k - x, k - y)]
    maps_band = [lambda x, y: (x, y), lambda x, y: (x, k - y), lambda x, y: (y, x), lambda x, y: (k - y, x)]
    for i, (x0, y0, x1, y1, m) in enumerate(base):
        for g in (maps_corner if i < ncorner else maps_band):
            p, q = g(x0, y0), g(x1, y1)
            segs.append((*p, *q, m))
    return segs


def key(s):
    x0, y0, x1, y1, m = s
    a, b = sorted([(x0, y0), (x1, y1)])
    return (a, b, m)


def check7():
    segs = box_segments(7)
    t = BOX7.read_text().split('\n')
    toks = ' '.join(l.split('#')[0] for l in t).split()
    i = 2
    def nx():
        nonlocal i
        i += 1
        return toks[i - 1]
    sn, sd = int(nx()), int(nx()); dd = int(nx()); W = int(nx())
    assert (sn, sd, dd, W) == (7, 1, 5, 10**12), (sn, sd, dd, W)
    assert int(nx()) == 0                       # no points
    ns = int(nx())
    filesegs = [tuple(F(int(nx()), dd) for _ in range(4)) + (F(int(nx()), W),) for _ in range(ns)]
    assert Counter(map(key, segs)) == Counter(map(key, filesegs)), 'family -> box 7 mismatch'
    return len(segs)


def total(k):
    return sum(s[4] for s in box_segments(k)) + (k - 2 * A) ** 2


def svg(k):
    n7 = check7()
    for kk in (6, 7, 8, 9, 12):
        assert total(kk) == kk * kk - 4 * D, kk
    segs = box_segments(k)
    px = 60                                     # user units per unit length
    pad = 14
    Z = k * px + 2 * pad
    X = lambda x: pad + float(x) * px
    Y = lambda y: pad + (k - float(y)) * px
    dens = [float(m) / float(abs(x1 - x0) + abs(y1 - y0)) for x0, y0, x1, y1, m in segs]
    dmax = max(dens)
    out = [f'<svg viewBox="0 0 {Z:g} {Z:g}" class="fam" role="img" aria-label="The measure for k = {k}: '
           f'mass on short grid segments near the walls and corners, even area density inside.">']
    out.append(f'<rect x="{X(A):g}" y="{Y(k - A):g}" width="{float(k - 2 * A) * px:g}" height="{float(k - 2 * A) * px:g}" class="leb"/>')
    out.append('<g class="seg">')
    for (x0, y0, x1, y1, m), d in sorted(zip(segs, dens), key=lambda t: t[1]):
        if m == 0:
            continue
        w = max(0.6, 7.0 * d / dmax)
        out.append(f'<path d="M{X(x0):g} {Y(y0):g}L{X(x1):g} {Y(y1):g}" stroke-width="{w:.2f}"/>')
    out.append('</g>')
    out.append(f'<rect x="{pad}" y="{pad}" width="{k * px}" height="{k * px}" class="frame"/>')
    for (cx, cy) in ((0, 0), (k - R, 0), (0, k - R), (k - R, k - R)):
        out.append(f'<rect x="{X(cx):g}" y="{Y(cy + R):g}" width="{R * px}" height="{R * px}" class="zone"/>')
    out.append('</svg>')
    sys.stderr.write(f'k = {k}: {len(segs)} segments, max density {dmax:.4f}; box 7 check OK ({n7} segments); '
                     f'totals k^2 - 4D OK for k = 6, 7, 8, 9, 12\n')
    return '\n'.join(out)


if __name__ == '__main__':
    print(svg(int(sys.argv[1]) if len(sys.argv) > 1 else 9))
