#!/usr/bin/env python3
"""Parse every fetched SVG, validate (count, containment, no overlaps), write data/packings.json."""
import os, sys, json, re, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_svg import parse, fmt
from mpmath import mpf, cos, sin, radians, mp
D = os.path.dirname(os.path.abspath(__file__)) + '/../data/'

def corners(cx, cy, th):
    c, s = cos(radians(th)), sin(radians(th)); h = mpf('0.5')
    return [(cx + c*dx - s*dy, cy + s*dx + c*dy) for dx, dy in ((-h,-h),(h,-h),(h,h),(-h,h))]

def sat_gap(A, B):
    """Separating-axis distance between two convex polygons (>0 apart, <0 overlapping)."""
    best = mpf('-inf')
    for P in (A, B):
        for k in range(len(P)):
            x1, y1 = P[k]; x2, y2 = P[(k+1) % len(P)]
            nx, ny = y2 - y1, x1 - x2
            L = (nx*nx + ny*ny) ** mpf('0.5'); nx, ny = nx/L, ny/L
            pa = [nx*x + ny*y for x, y in A]; pb = [nx*x + ny*y for x, y in B]
            gap = max(min(pb) - max(pa), min(pa) - max(pb))
            if gap > best: best = gap
    return best

# Some of Ellsworth's packings are numerical (simulated annealing) rather than analytic, and their
# constants are given rounded to ~10-16 digits, so they miss exactness by ~1e-11.  Violations below
# TOL_ERR are reported as warnings only; anything above is a real error.
TOL_ERR = mpf('1e-9')     # above this: an error
TOL_WARN = mpf('1e-12')   # above this (but below TOL_ERR): a 'warn:' note

def validate(s, sq):
    errs = []
    polys = [corners(*q) for q in sq]
    worst_out = mpf(0)
    for i, P in enumerate(polys):
        for x, y in P:
            v = max(-x, -y, x - s, y - s)
            if v > worst_out: worst_out = v
            if v > TOL_ERR:
                errs.append('square %d outside container by %s' % (i, fmt(v, 5))); break
    worst_ovl = mpf(0)
    for i in range(len(polys)):
        for j in range(i+1, len(polys)):
            ci, cj = sq[i], sq[j]
            if (ci[0]-cj[0])**2 + (ci[1]-cj[1])**2 > 2.01: continue   # farther than sqrt2 apart: cannot overlap
            g = sat_gap(polys[i], polys[j])
            if -g > worst_ovl: worst_ovl = -g
            if -g > TOL_ERR: errs.append('overlap %d,%d depth %s' % (i, j, fmt(-g, 5)))
    if TOL_WARN < worst_out <= TOL_ERR:
        errs.append('warn: max containment violation %s (numeric, not analytic)' % fmt(worst_out, 5))
    if TOL_WARN < worst_ovl <= TOL_ERR:
        errs.append('warn: max overlap depth %s (numeric, not analytic)' % fmt(worst_ovl, 5))
    return errs

def expected_n(name):
    """The number of unit squares the file should contain.  Usually the number in the filename, but
    the "*_no_min_rot" variants pack one square fewer into the s(n) container (Ellsworth's
    Göbel-strips page labels square-369_no_min_rot.svg as "368.", 586 as "585.", etc.)."""
    m = re.match(r'square-(\d+)', name)
    if not m: return None
    n = int(m.group(1))
    return n - 1 if '_no_min_rot' in name else n

if __name__ == '__main__':
    # optional argv: only re-parse the named files (basenames or paths) and merge into packings.json
    only = [os.path.basename(a) for a in sys.argv[1:]]
    files = sorted(glob.glob(D + 'ellsworth/svg/*.svg'))
    if only:
        files = [f for f in files if os.path.basename(f) in only]
        missing = set(only) - set(os.path.basename(f) for f in files)
        if missing: sys.exit('no such svg: %s' % ', '.join(sorted(missing)))
    out, bad = {}, []
    for f in files:
        name = os.path.basename(f)
        n_expected = expected_n(name)
        try:
            s, sq = parse(f)
        except Exception as e:
            bad.append((name, 'parse: %s' % e)); continue
        errs = validate(s, sq)
        if n_expected is not None and len(sq) != n_expected:
            # too many is as much a bug as too few (e.g. a <use> resolved to the wrong duplicate id)
            errs.append('count %d != %d' % (len(sq), n_expected))
        hard = [e for e in errs if not e.startswith('warn:')]
        if hard: bad.append((name, '; '.join(hard[:3])))
        out[name] = {'s': fmt(s), 'n': len(sq), 'n_file': n_expected,
                     'squares': [[fmt(a), fmt(b), fmt(c)] for a, b, c in sq], 'errors': errs}
    if only and os.path.exists(D + 'packings.json'):
        merged = json.load(open(D + 'packings.json')); merged.update(out)
    else:
        merged = out
    json.dump(merged, open(D + 'packings.json.tmp', 'w'))
    os.replace(D + 'packings.json.tmp', D + 'packings.json')   # atomic: export.py may be reading it
    print('parsed %d files, %d with problems' % (len(out), len(bad)))
    for name, e in bad: print('  ', name, '-', e[:160])
