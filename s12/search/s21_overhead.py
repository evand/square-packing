#!/usr/bin/env python3
"""S21_OVERHEAD: region decomposition of the m=4 covers and the nu_f measure (analysis only).

  python3 search/s21_overhead.py regions            # per-cell / per-class tables
Reads (never writes) certificates/rung2/s13_closed_cover_4.txt, runs/closed4_best.txt,
search/cover4_exact_support.txt.  Exact Fraction sums for the covers; the measure is exact too.
"""
import sys
from fractions import Fraction as F

def read_cover(path):
    toks = open(path).read().split()
    it = iter(toks)
    a, b = int(next(it)), int(next(it)); s = F(a, b)
    D = int(next(it)); W = int(next(it)); n = int(next(it))
    pts = []
    for _ in range(n):
        x, y, w = int(next(it)), int(next(it)), int(next(it))
        pts.append((F(x, D), F(y, D), F(w, W)))
    return s, pts

def read_measure(path):
    t = None; poses = []
    for line in open(path):
        if line.startswith('# t ='):
            t = F(line.split('=')[1].split()[0])
        if line.startswith('pose'):
            _, p, q, cx, cy, m = line.split()
            poses.append((int(p), int(q), F(cx), F(cy), F(m)))
    imgs = []
    for p, q, cx, cy, m in poses:
        for (x, y) in [(cx, cy), (t-cx, cy), (cx, t-cy), (t-cx, t-cy), (cy, cx), (t-cy, cx), (cy, t-cx), (t-cy, t-cx)]:
            imgs.append((x, y, m / 8, p, q))
    return t, imgs

def cells_of(x, m):
    """list of (cell index, share) along one axis; a point on an interior integer line is split."""
    if x.denominator == 1 and 0 < x < m:
        return [(int(x) - 1, F(1, 2)), (int(x), F(1, 2))]
    i = min(int(x), m - 1)
    return [(i, F(1))]

def klass(i, j, m):
    e = (i in (0, m-1)) + (j in (0, m-1))
    return {2: 'corner', 1: 'edge', 0: 'interior'}[e]

def per_cell(pts, m):
    g = {}
    for x, y, w in pts:
        for i, a in cells_of(x, m):
            for j, b in cells_of(y, m):
                g[(i, j)] = g.get((i, j), 0) + w * a * b
    return g

def linetype(x, y, m):
    on = lambda z: z.denominator == 1
    wall = lambda z: z == 0 or z == m
    kx = 'W' if wall(x) else ('L' if on(x) else '-')
    ky = 'W' if wall(y) else ('L' if on(y) else '-')
    k = ''.join(sorted(kx + ky))
    return {'--': 'off-line', '-L': 'on 1 interior line', '-W': 'on wall', 'LL': 'grid vertex',
            'LW': 'wall∩line', 'WW': 'container corner'}[k]

def regions():
    m = 4
    covers = [('certified (x1.05)', 'certificates/rung2/s13_closed_cover_4.txt'),
              ('closed4_best (LP 12.4175)', 'runs/closed4_best.txt')]
    for name, path in covers:
        s, pts = read_cover(path)
        tot = sum(w for _, _, w in pts)
        g = per_cell(pts, m)
        print(f'== {name}: {len(pts)} pts, total {float(tot):.6f}')
        for j in reversed(range(m)):
            print('   ' + '  '.join(f'{float(g.get((i, j), 0)):.4f}' for i in range(m)))
        cls = {}
        for (i, j), v in g.items():
            cls.setdefault(klass(i, j, m), []).append(v)
        for k in ('corner', 'edge', 'interior'):
            v = cls[k]; print(f'   {k:9s} n={len(v):2d} sum={float(sum(v)):.4f} mean={float(sum(v)/len(v)):.4f}')
        lt = {}
        for x, y, w in pts:
            k = linetype(x, y, m); lt[k] = lt.get(k, 0) + w
        print('   by line type: ' + ', '.join(f'{k} {float(v):.4f}' for k, v in sorted(lt.items(), key=lambda t: -t[1])))
        # distance band of each point to nearest wall
        bands = {}
        for x, y, w in pts:
            d = min(x, y, m - x, m - y)
            b = '<1' if d < 1 else ('=1' if d == 1 else ('(1,2)' if d < 2 else '=2'))
            bands[b] = bands.get(b, 0) + w
        print('   by wall distance d of point: ' + ', '.join(f'd{k}: {float(v):.4f}' for k, v in bands.items()))
    t, imgs = read_measure('search/cover4_exact_support.txt')
    mass = sum(mm for _, _, mm, _, _ in imgs)
    g = {}; ang = {}
    for x, y, mm, p, q in imgs:
        i = min(int(x), m - 1); j = min(int(y), m - 1)
        g[(i, j)] = g.get((i, j), 0) + mm
        k = klass(i, j, m) + ('/axis' if p == 0 else '/tilted')
        ang[k] = ang.get(k, 0) + mm
    print(f'== nu_f measure (pose centres), mass {float(mass):.6f}')
    for j in reversed(range(m)):
        print('   ' + '  '.join(f'{float(g.get((i, j), 0)):.4f}' for i in range(m)))
    cls = {}
    for (i, j), v in g.items():
        cls.setdefault(klass(i, j, m), []).append(v)
    for k in ('corner', 'edge', 'interior'):
        v = cls[k]; print(f'   {k:9s} n={len(v):2d} sum={float(sum(v)):.4f} mean={float(sum(v)/len(v)):.4f}')
    print('   ' + ', '.join(f'{k} {float(v):.4f}' for k, v in sorted(ang.items())))

if __name__ == '__main__':
    {'regions': regions}[sys.argv[1]]()
