#!/usr/bin/env python3
"""Parse one of David Ellsworth's packing SVGs into a list of unit squares.

Independent implementation (not derived from svgDisp.js). Conventions of the SVG dialect:
  * the <rect id="outer"> gives the container side s (it may be positioned at x=-s/2 when centred);
  * unit squares are <rect width=1 height=1>, larger integer rects (blocks), or filled rectilinear
    <path>s (polyominoes) - all in a local frame reached through nested <g transform> and <use>;
  * anything with fill:none is decoration (edge lines) and is skipped;
  * XML entities in the DOCTYPE carry the exact constants (30+ digits).
Output: container side, and squares as (cx, cy, theta_deg) with theta in [0, 90).
All arithmetic in mpmath at 50 digits.
"""
import re, sys, json
import xml.etree.ElementTree as ET
from mpmath import mp, mpf, cos, sin, atan2, pi, radians, degrees, floor, nint
mp.dps = 50

NS = {'svg': 'http://www.w3.org/2000/svg', 'xlink': 'http://www.w3.org/1999/xlink'}
XLINK = '{http://www.w3.org/1999/xlink}href'

def tag(el):
    return el.tag.split('}')[-1]

# ---------- affine transforms: (a b c d e f) maps (x,y) -> (a x + c y + e, b x + d y + f) ----------
def num(x):
    x = x.strip()
    x = re.sub(r'^([-+]?)\.', r'\g<1>0.', x)
    return mpf(x)

def T_id(): return (mpf(1), mpf(0), mpf(0), mpf(1), mpf(0), mpf(0))
def T_mul(A, B):  # A after B  (apply B first, then A)
    a, b, c, d, e, f = A; g, h, i, j, k, l = B
    return (a*g + c*h, b*g + d*h, a*i + c*j, b*i + d*j, a*k + c*l + e, b*k + d*l + f)
def T_apply(T, x, y):
    a, b, c, d, e, f = T
    return (a*x + c*y + e, b*x + d*y + f)

NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
def parse_transform(s):
    T = T_id()
    for name, args in re.findall(r'(\w+)\s*\(([^)]*)\)', s or ''):
        v = [num(x) for x in re.findall(NUM, args)]
        if name == 'translate':
            tx = v[0]; ty = v[1] if len(v) > 1 else mpf(0)
            M = (mpf(1), mpf(0), mpf(0), mpf(1), tx, ty)
        elif name == 'scale':
            sx = v[0]; sy = v[1] if len(v) > 1 else sx
            M = (sx, mpf(0), mpf(0), sy, mpf(0), mpf(0))
        elif name == 'rotate':
            th = radians(v[0]); c, s_ = cos(th), sin(th)
            M = (c, s_, -s_, c, mpf(0), mpf(0))
            if len(v) == 3:
                cx, cy = v[1], v[2]
                M = T_mul((mpf(1), mpf(0), mpf(0), mpf(1), cx, cy), T_mul(M, (mpf(1), mpf(0), mpf(0), mpf(1), -cx, -cy)))
        elif name == 'matrix':
            M = tuple(v)
        elif name == 'skewX' or name == 'skewY':
            raise ValueError('skew not supported')
        else:
            raise ValueError('unknown transform ' + name)
        T = T_mul(T, M)
    return T

# ---------- paths -> rectilinear polygons (list of vertex lists) ----------
def parse_path(d):
    toks = re.findall(r'[MmLlHhVvZz]|' + NUM, d)
    subpaths, cur = [], []
    x = y = sx = sy = mpf(0); cmd = None; i = 0
    while i < len(toks):
        t = toks[i]
        if re.match(r'[A-Za-z]', t):
            cmd = t; i += 1
            if cmd in 'Zz':
                if cur: subpaths.append(cur); cur = []
                x, y = sx, sy
            continue
        if cmd is None: raise ValueError('path without command')
        if cmd in 'Mm':
            nx, ny = num(toks[i]), num(toks[i+1]); i += 2
            if cmd == 'm': nx += x; ny += y
            if cur: subpaths.append(cur)
            cur = [(nx, ny)]; x, y = nx, ny; sx, sy = x, y
            cmd = 'L' if cmd == 'M' else 'l'   # subsequent pairs are implicit lineto
        elif cmd in 'Ll':
            nx, ny = num(toks[i]), num(toks[i+1]); i += 2
            if cmd == 'l': nx += x; ny += y
            x, y = nx, ny; cur.append((x, y))
        elif cmd in 'Hh':
            nx = num(toks[i]); i += 1
            if cmd == 'h': nx += x
            x = nx; cur.append((x, y))
        elif cmd in 'Vv':
            ny = num(toks[i]); i += 1
            if cmd == 'v': ny += y
            y = ny; cur.append((x, y))
        else:
            raise ValueError('unsupported path command ' + cmd)
    if cur: subpaths.append(cur)
    return subpaths

def winding(polys, px, py):
    w = 0
    for poly in polys:
        n = len(poly)
        for k in range(n):
            x1, y1 = poly[k]; x2, y2 = poly[(k+1) % n]
            if y1 <= py < y2 or y2 <= py < y1:
                xi = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
                if xi > px: w += 1 if y2 > y1 else -1
    return w

def drop_collinear(poly):
    out = list(poly)
    changed = True
    while changed and len(out) > 3:
        changed = False
        for k in range(len(out)):
            (xa, ya), (xb, yb), (xc, yc) = out[k-1], out[k], out[(k+1) % len(out)]
            if (abs(xa-xb) < mpf('1e-30') and abs(xb-xc) < mpf('1e-30')) or (abs(ya-yb) < mpf('1e-30') and abs(yb-yc) < mpf('1e-30')) or (abs(xa-xb) < mpf('1e-30') and abs(ya-yb) < mpf('1e-30')):
                del out[k]; changed = True; break
    return out

EPS = mpf('1e-30')

def rect_minus(r, q):
    """r minus q for axis-aligned rects (x0,y0,x1,y1); returns list of rects and the removed area."""
    x0, y0, x1, y1 = r; a0, b0, a1, b1 = q
    ix0, iy0, ix1, iy1 = max(x0, a0), max(y0, b0), min(x1, a1), min(y1, b1)
    if ix0 >= ix1 - EPS or iy0 >= iy1 - EPS:
        return [r], mpf(0)
    out = []
    if iy0 > y0 + EPS: out.append((x0, y0, x1, iy0))
    if iy1 < y1 - EPS: out.append((x0, iy1, x1, y1))
    if ix0 > x0 + EPS: out.append((x0, iy0, ix0, iy1))
    if ix1 < x1 - EPS: out.append((ix1, iy0, x1, iy1))
    return out, (ix1 - ix0) * (iy1 - iy0)

def path_cells(d, fill_rule='nonzero'):
    """Decompose the filled region of a rectilinear path (a union of closed unit squares, not
    necessarily on one lattice) into unit squares; returns their lower-left corners."""
    polys = [drop_collinear(p) for p in parse_path(d)]
    pts = [p for poly in polys for p in poly]
    xs = sorted(set(p[0] for p in pts)); ys = sorted(set(p[1] for p in pts))
    # arrangement cells that are inside
    rects = []
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            cx, cy = (xs[i] + xs[i+1]) / 2, (ys[j] + ys[j+1]) / 2
            w = winding(polys, cx, cy)
            if (w != 0) if fill_rule == 'nonzero' else (w % 2 == 1):
                rects.append((xs[i], ys[j], xs[i+1], ys[j+1]))
    cells = []
    while rects:
        x, y = min((r[1], r[0]) for r in rects)[::-1]
        q = (x, y, x + 1, y + 1)
        new, removed = [], mpf(0)
        for r in rects:
            parts, a = rect_minus(r, q); new.extend(parts); removed += a
        if abs(removed - 1) > mpf('1e-25'):
            raise ValueError('region is not a union of unit squares (removed %s at %s,%s): %s' % (removed, x, y, d[:60]))
        rects = new; cells.append((x, y))
        if len(cells) > 100000: raise ValueError('runaway')
    return cells

# ---------- walk the document ----------
def fill_none(el):
    st = el.get('style') or ''
    return re.search(r'fill\s*:\s*none', st) is not None

def parse(path):
    raw = open(path, encoding='utf-8').read()
    root = ET.fromstring(raw)          # expat expands the internal DTD entities
    # duplicate ids do occur; getElementById (and hence the browser) keeps the FIRST in document order
    ids = {}
    for el in root.iter():
        i = el.get('id')
        if i is not None and i not in ids: ids[i] = el
    outer = ids.get('outer')
    if outer is not None:
        s = num(outer.get('width'))
        ox = num(outer.get('x') or '0'); oy = num(outer.get('y') or '0')
    else:
        vb = [num(v) for v in re.findall(NUM, root.get('viewBox') or '')]
        if len(vb) != 4 or abs(vb[2] - vb[3]) > mpf('1e-30'): raise ValueError('no #outer and no square viewBox')
        ox, oy, s = vb[0], vb[1], vb[2]
    squares = []   # (cx, cy, theta) in outer-frame coords with origin at outer's corner

    trace = []
    def emit(T, x, y):
        cx, cy = T_apply(T, x + mpf('0.5'), y + mpf('0.5'))
        a, b, c, d, e, f = T
        th = atan2(b, a)                                   # rotation of the local x axis
        squares.append((cx - ox, cy - oy, th)); trace.append('/'.join(chain))

    chain = []
    def walk(el, T, inherited_none, depth=0):
        if depth > 50: raise ValueError('recursion')
        chain.append(el.get('id') or tag(el))
        try: _walk(el, T, inherited_none, depth)
        finally: chain.pop()
    def _walk(el, T, inherited_none, depth):
        t = tag(el)
        if t in ('defs', 'script', 'title', 'desc', 'metadata', 'style'): return
        if el.get('id') == 'outer': return
        none = inherited_none or fill_none(el)
        T = T_mul(T, parse_transform(el.get('transform')))
        if t == 'use':
            ux, uy = num(el.get('x') or '0'), num(el.get('y') or '0')
            if ux or uy: T = T_mul(T, (mpf(1), mpf(0), mpf(0), mpf(1), ux, uy))
            ref = (el.get(XLINK) or el.get('href') or '')
            if not ref.startswith('#'): return
            target = ids.get(ref[1:])
            if target is None: raise ValueError('missing ref ' + ref)
            if target.get('id') == 'outer': return
            walk(target, T, none, depth + 1)
            return
        if none and t in ('rect', 'path'): return
        if t == 'rect':
            w, h = num(el.get('width')), num(el.get('height'))
            x, y = num(el.get('x') or '0'), num(el.get('y') or '0')
            W, H = int(nint(w)), int(nint(h))
            if abs(w - W) > mpf('1e-30') or abs(h - H) > mpf('1e-30') or W < 1 or H < 1:
                raise ValueError('non-integer rect %sx%s' % (w, h))
            for ix in range(W):
                for iy in range(H):
                    emit(T, x + ix, y + iy)
        elif t == 'path':
            fr = 'evenodd' if re.search(r'fill-rule\s*:\s*evenodd', el.get('style') or '') or el.get('fill-rule') == 'evenodd' else 'nonzero'
            for (x, y) in path_cells(el.get('d'), fr):
                emit(T, x, y)
        elif t in ('g', 'svg', 'a', 'switch'):
            for ch in el:
                walk(ch, T, none, depth + 1)
        elif t in ('circle', 'line', 'text', 'polyline', 'polygon', 'ellipse', 'image'):
            return
        else:
            sys.stderr.write('note: skipping <%s>\n' % t)

    walk(root, T_id(), False)
    # convert to a y-up frame (Ellsworth's SVGs flip y once at the top; the maths frame is y-up)
    out = []
    for (cx, cy, th) in squares:
        thd = (-degrees(th)) % 90
        if thd > 90 - mpf('1e-25'): thd -= 90
        if thd < 0: thd = mpf(0)
        out.append((cx, s - cy, thd))
    parse.trace = trace
    return s, out

def fmt(x, digits=30):
    return mp.nstr(x, digits, strip_zeros=True)

if __name__ == '__main__':
    for f in sys.argv[1:]:
        s, sq = parse(f)
        print(json.dumps({'file': f, 's': fmt(s), 'n': len(sq),
                          'squares': [[fmt(a), fmt(b), fmt(c)] for a, b, c in sq]}))
