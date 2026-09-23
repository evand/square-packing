#!/usr/bin/env python3
"""bandcut_scan_svg.py -- read a published packing from kingbird.myphotos.cc and turn it into
an explicit list of (centre_x, centre_y, angle) for our own margin instrument.

The pages at <https://kingbird.myphotos.cc/packing/square-N.svg> are the current-records
successor to Friedman's Packing Center.  Each file gives the packing exactly: blocks of
axis-parallel squares as filled rectilinear regions on the integer grid, and every other
square as `<use href="#one" transform="translate(cx cy) rotate(a) translate(-.5 -.5)"/>`.

    python3 search/bandcut_scan_svg.py square-132.svg --T 12

prints the square count, the bounding box, `delta` at the file's own side `s`, and `delta`
after the dilation of `S6_SKELETON.md` S0 to the integer container `[0,T]^2`.  A positive
number there is a feasible point for `delta*_T`, i.e. `s(T^2-T) < T`.
"""
import argparse
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bandcut_scan import float_value                                  # noqa: E402

SVG = '{http://www.w3.org/2000/svg}'


def mat_mul(A, B):
    (a, b, c, d, e, f), (p, q, r, s, t, u) = A, B
    return (a * p + c * q, b * p + d * q, a * r + c * s, b * r + d * s,
            a * t + c * u + e, b * t + d * u + f)


def parse_transform(txt):
    M = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    if not txt:
        return M
    for name, args in re.findall(r'(translate|scale|rotate|matrix)\s*\(([^)]*)\)', txt):
        v = [float(x) for x in re.split(r'[,\s]+', args.strip()) if x]
        if name == 'translate':
            T = (1, 0, 0, 1, v[0], v[1] if len(v) > 1 else 0.0)
        elif name == 'scale':
            sx = v[0]; sy = v[1] if len(v) > 1 else v[0]
            T = (sx, 0, 0, sy, 0, 0)
        elif name == 'rotate':
            a = math.radians(v[0]); c, s = math.cos(a), math.sin(a)
            T = (c, s, -s, c, 0, 0)
        else:
            T = tuple(v)
        M = mat_mul(M, T)
    return M


def apply(M, x, y):
    a, b, c, d, e, f = M
    return (a * x + c * y + e, b * x + d * y + f)


def path_polygon(d):
    """rectilinear path `M x,y H.. V.. H..` -> closed polygon (list of points)."""
    toks = re.findall(r'([MHVLZmhvlz])\s*([-0-9.eE,\s]*)', d)
    pts, cur = [], (0.0, 0.0)
    for cmd, arg in toks:
        nums = [float(x) for x in re.split(r'[,\s]+', arg.strip()) if x]
        if cmd in 'Mm':
            cur = (nums[0], nums[1]); pts.append(cur)
            for k in range(2, len(nums), 2):
                cur = (nums[k], nums[k + 1]); pts.append(cur)
        elif cmd in 'Hh':
            for x in nums:
                cur = (x, cur[1]); pts.append(cur)
        elif cmd in 'Vv':
            for y in nums:
                cur = (cur[0], y); pts.append(cur)
        elif cmd in 'Ll':
            for k in range(0, len(nums), 2):
                cur = (nums[k], nums[k + 1]); pts.append(cur)
    return pts


def cells_in(poly):
    xs = [p[0] for p in poly]; ys = [p[1] for p in poly]
    out = []
    for ix in range(int(math.floor(min(xs))), int(math.ceil(max(xs)))):
        for iy in range(int(math.floor(min(ys))), int(math.ceil(max(ys)))):
            cx, cy = ix + 0.5, iy + 0.5
            inside, n = False, len(poly)
            for k in range(n):
                x1, y1 = poly[k]; x2, y2 = poly[(k + 1) % n]
                if (y1 > cy) != (y2 > cy):
                    xin = x1 + (cy - y1) * (x2 - x1) / (y2 - y1)
                    if xin > cx:
                        inside = not inside
            if inside:
                out.append((cx, cy))
    return out


def walk(el, M, defs, out):
    M = mat_mul(M, parse_transform(el.get('transform')))
    tag = el.tag.replace(SVG, '')
    if tag == 'use':
        href = el.get('{http://www.w3.org/1999/xlink}href') or el.get('href') or ''
        tgt = href.lstrip('#')
        if tgt == 'one':
            M2 = mat_mul(M, (1, 0, 0, 1, 0.5, 0.5))     # #one is 1x1 at origin, then translate(-.5,-.5)
            cx, cy = apply(mat_mul(M, (1, 0, 0, 1, 0, 0)), 0.0, 0.0)
            # the element's own transform already ends with translate(-.5 -.5): centre is
            # the image of the local point (0.5, 0.5)
            cx, cy = apply(M, 0.5, 0.5)
            ang = math.atan2(M[1], M[0])
            out.append((cx, cy, ang))
        elif tgt in defs and tgt != 'outer':
            walk(defs[tgt], M, defs, out)
        return
    if tag == 'rect':
        w = float(el.get('width', 0) or 0); h = float(el.get('height', 0) or 0)
        if abs(w - 1.0) < 1e-9 and abs(h - 1.0) < 1e-9:
            x0 = float(el.get('x', 0) or 0); y0 = float(el.get('y', 0) or 0)
            cx, cy = apply(M, x0 + 0.5, y0 + 0.5)
            out.append((cx, cy, math.atan2(M[1], M[0])))
        return
    if tag == 'path':
        st = el.get('style', '')
        if 'stroke:none' in st:
            poly = path_polygon(el.get('d', ''))
            for (cx, cy) in cells_in(poly):
                X, Y = apply(M, cx, cy)
                out.append((X, Y, math.atan2(M[1], M[0])))
        return
    for ch in el:
        walk(ch, M, defs, out)


def load(path):
    src = open(path).read()
    ents = dict(re.findall(r'<!ENTITY\s+(\w+)\s+"([^"]+)"', src))
    body = src[src.index('<svg'):]
    for k, v in ents.items():
        body = body.replace('&%s;' % k, v)
    root = ET.fromstring(body)
    defs = {}
    for el in root.iter():
        if el.get('id'):
            defs[el.get('id')] = el
    out = []
    for ch in root:
        if ch.tag.replace(SVG, '') in ('defs', 'script'):
            continue
        walk(ch, (1, 0, 0, 1, 0, 0), defs, out)
    s = float(ents.get('s', '0'))
    return s, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('--T', type=float, default=None)
    ap.add_argument('--map', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    s, sq = load(a.file)
    # the containing group is translate(-hs,-hs) for the tilted squares and 0 for the blocks
    X = np.array([p[0] for p in sq]); Y = np.array([p[1] for p in sq])
    TH = np.array([p[2] for p in sq])
    if X.min() < 0:                       # mixed frames: shift everything to [0,s]^2
        pass
    lo = min(X.min(), Y.min())
    if lo < -1e-9:
        X = X + s / 2.0; Y = Y + s / 2.0
    p = 0.5 * (np.abs(np.cos(TH)) + np.abs(np.sin(TH)))
    print('file       : %s' % a.file)
    print('side s     : %.20f' % s)
    print('squares    : %d' % len(sq))
    print('bbox x     : %.12f .. %.12f' % ((X - p).min(), (X + p).max()))
    print('bbox y     : %.12f .. %.12f' % ((Y - p).min(), (Y + p).max()))
    v, kp, wall = float_value(X, Y, TH, s)
    print('delta @ s  : %.6e   (worst pair %s, wall %.3e)' % (v, kp, wall))
    if a.T:
        lam = a.T / s
        v, kp, wall = float_value(X * lam, Y * lam, TH, a.T)
        print('delta @ T=%g: %.9e  (worst pair %s, wall %.3e)' % (a.T, v, kp, wall))
    deg = (np.degrees(TH) % 90.0)
    deg = np.where(deg > 89.999, deg - 90.0, deg)
    import collections
    cnt = collections.Counter(np.round(deg, 2))
    print('angle census (tilt mod 90, deg):')
    for k, c in sorted(cnt.items()):
        print('   %8.3f : %3d' % (k, c))
    if a.map:
        T = int(round(s))
        grid = [['.'] * T for _ in range(T)]
        for (x, y, th) in zip(X, Y, TH):
            ix, iy = int(math.floor(x * T / s)), int(math.floor(y * T / s))
            d = (math.degrees(th) % 90.0)
            ch = '#' if (d < 3 or d > 87) else ('/' if d < 45 else '\\')
            if 0 <= ix < T and 0 <= iy < T:
                grid[iy][ix] = ch
        print('occupancy map (row 0 at the bottom):')
        for r in range(T - 1, -1, -1):
            print('  %2d %s' % (r, ''.join(grid[r])))
    if a.out:
        import json
        json.dump({'s': s, 'sq': [[float(x), float(y), float(t)]
                                  for x, y, t in zip(X, Y, TH)]}, open(a.out, 'w'))
        print('wrote %s' % a.out)


if __name__ == '__main__':
    main()
