#!/usr/bin/env python3
"""Structure of the best-known packings mirrored on the site (WISHLIST section P).

Reads the site's best-known table (site/www/data/index.json) and drawings (site/www/data/p/) and prints:
  1. per k: non-integer records, how many use <= 3 rotation groups, lowest fractional part using >= 5 groups;
  2. the margin table k - s(k^2 - c) near the k^2 - c frontier ("=k": best known is the grid; "?": n not on the site,
     i.e. a trivial fill-in with s = k or below the k-interval);
  3. fixed-fractional-part families: fractional parts shared by >= 4 records, with the n that carry them.
Rotation group = distinct angle mod 90 deg (rounded to 1e-5 deg; a mirror pair theta, 90 - theta counts as two).
Data are best-known packings (Ellsworth's table as mirrored by the site build), not proved optima.
"""
import glob, json, math, os, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PDIR = os.path.join(HERE, '..', '..', 'site', 'www', 'data', 'p')


def groups_of(path):
    angs = set()
    for _x, _y, a in json.load(open(path))['squares']:
        a %= 90
        if a > 90 - 1e-7:
            a = 0.0
        angs.add(round(a, 5))
    return len(angs)


def load():
    """Best-known table: n -> (s, rotation groups).  Main-table entries (index.json) give s; their drawing
    (pictured_n, the svg's json) gives the groups.  Drawn n <= 324 missing from the main table are the grids
    (s = sqrt n) and the trivial k^2 - c fill-ins (s = k).  Larger drawn n that are not in the main table are
    sub-page drawings (Goebel strips/squares), not records: returned separately."""
    idx = json.load(open(os.path.join(PDIR, '..', 'index.json')))['records']
    rec, sub = {}, {}
    for v in idx.values():
        js = os.path.join(PDIR, v['svg'].replace('.svg', '.json'))
        if not os.path.exists(js):
            js = os.path.join(PDIR, f"square-{v['pictured_n']}.json")
        rec[v['n']] = (float(v['s_dec']), groups_of(js) if os.path.exists(js) else None)
    for f in glob.glob(os.path.join(PDIR, 'square-*.json')):
        m = re.fullmatch(r'square-(\d+)\.json', os.path.basename(f))
        if not m or int(m.group(1)) in rec:
            continue
        d = json.load(open(f))
        (rec if d['n'] <= 324 else sub)[d['n']] = (float(d['s']), groups_of(f))
    return rec, sub


def main():
    rec, sub = load()
    print(f'{len(rec)} best-known entries ({len(sub)} sub-page drawings above n = 324 set aside: {sorted(sub)})\n')
    print('1. rotation groups per k')
    for k in range(4, 18):
        r = [(n, s - k, g) for n, (s, g) in rec.items() if k + 1e-9 < s < k + 1 - 1e-12]
        messy = [x for x in r if x[2] is not None and x[2] >= 5]
        lo = min(messy, key=lambda x: x[1]) if messy else None
        print(f'  k={k:2d}  non-int {len(r):2d}  <=3 groups {sum(g is not None and g <= 3 for *_, g in r):2d}  '
              f'lowest messy: {"-" if lo is None else f"n={lo[0]} frac={lo[1]:.3f} groups={lo[2]}"}')
    print('\n2. margin k - s(k^2 - c)')
    cs = range(8, 17)
    print('  k   ' + ' '.join(f'c={c:<4d}' for c in cs))
    for k in range(8, 19):
        row = []
        for c in cs:
            n = k * k - c
            if n not in rec:
                row.append('  ?   ')
            elif rec[n][0] < k - 1e-9:
                row.append(f'{k - rec[n][0]:.4f}')
            else:
                row.append('  =k  ')
        print(f'  {k:2d}  ' + ' '.join(row))
    print('\n3. fractional parts shared by >= 4 records')
    fam = defaultdict(list)
    for n, (s, g) in sorted(rec.items()):
        f = s - math.floor(s + 1e-9)
        if f > 1e-9:
            fam[round(f, 8)].append(n)
    for f, ns in sorted(fam.items(), key=lambda x: -len(x[1])):
        if len(ns) >= 4:
            print(f'  frac {f:.8f}: {len(ns):2d} records, n = {ns}')


if __name__ == '__main__':
    main()
