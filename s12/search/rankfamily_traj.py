#!/usr/bin/env python3
"""Trajectory of a `cliquelever.py --master --lattice-every` run: LP value against the loaded
column count, with the lattice injections marked.

    python3 search/rankfamily_traj.py runs/E2A.out [runs/E2B.out ...] [--every N] [--tail N]

Each row is one loop iteration: the LP value (an upper bound on the QSTAB + rank value of the
loaded column set at every iteration, whatever the master did -- `CLMASTER.md` §1), the support
size, the loaded column count, the clique and polygon row counts, `M`, the max clique, and the
best violated polygon that iteration separated.  A `>>>` line is a lattice injection: how many
poses came in, the pricing gap, and how many new clique and polygon memberships the existing
rows picked up (the polygon count is `rankfamily.regrow` re-deriving every row over the enlarged
pose set).
"""
import argparse
import os
import re
import sys

ITHEAD = re.compile(r'\[(\w+)\.(\d+)\] LP=([\d.]+) M=([\d.]+) kmax=([\d.]+)')
F = {
    'rows': re.compile(r' rows=(\d+)'),
    'cq': re.compile(r' cq=(\d+)'),
    'sup': re.compile(r' sup=(\d+)'),
    'pg': re.compile(r' pg=(\d+)\(\+(\d+)\)'),
    'master': re.compile(r' master=(\d+)/(\d+)'),
    'pbest': re.compile(r' pbest=k(\d+):([\d.]+)/([\d.]+)\(\+([\d.]+)\)/(\d+)'),
    'secs': re.compile(r'(\d+)s\)$'),
}
LAT = re.compile(r'\[(\w+)\.(\d+)\] lattice: \+(\d+) poses -> (\d+) poses, (\d+) columns '
                 r'\(gap ([-+\d.]+), (\d+) new clique memberships'
                 r'(?:, (\d+) new polygon memberships)?')
RES = re.compile(r'^RESULT (.*)$')
FIN = re.compile(r'FINAL EXACT: mass = \S+ = ([\d.]+); M = \S+ = ([\d.]+)\s+(\w+); '
                 r'max clique = \S+ = ([\d.]+).*?(\w+); regions (\w+)')
PGF = re.compile(r'FINAL EXACT: (\d+) rank-family rows: (\w+) \(worst violation ([-+\d.e]+), '
                 r'(\d+) tight, all re-derived from their exact anchors: (\w+)\)')


def parse_it(s):
    m = ITHEAD.search(s)
    if not m:
        return None
    d = dict(tag=m.group(1), it=int(m.group(2)), LP=float(m.group(3)),
             M=float(m.group(4)), kmax=float(m.group(5)))
    for k, rx in F.items():
        g = rx.search(s)
        d[k] = g.groups() if g else None
    return d


def rows_of(path):
    out = []
    for line in open(path, errors='replace'):
        s = line.rstrip('\n')
        m = LAT.search(s)
        if m:
            out.append(('lat', m.groups()))
            continue
        d = parse_it(s.strip())
        if d:
            out.append(('it', d))
            continue
        m = PGF.search(s)
        if m:
            out.append(('pgfin', m.groups()))
            continue
        m = FIN.search(s)
        if m:
            out.append(('fin', m.groups()))
            continue
        m = RES.match(s)
        if m:
            out.append(('res', m.groups()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--every', type=int, default=1, help='print every Nth iteration')
    ap.add_argument('--tail', type=int, default=0, help='only the last N iterations')
    args = ap.parse_args()
    for path in args.files:
        rec = rows_of(path)
        its = [r for k, r in rec if k == 'it']
        if not its:
            print(f"\n=== {os.path.basename(path)}: no iterations yet")
            continue
        lo = its[-args.tail]['it'] if args.tail and len(its) > args.tail else its[0]['it']
        print(f"\n=== {os.path.basename(path)}: {len(its)} iterations, "
              f"{sum(1 for k, _ in rec if k == 'lat')} lattice injections, "
              f"LP {its[0]['LP']:.6f} -> {its[-1]['LP']:.6f}")
        print(f"{'it':>5} {'LP':>12} {'sup':>6} {'cols':>7} {'cq':>6} {'pg':>5} "
              f"{'M':>10} {'kmax':>9} {'best polygon':>24} {'secs':>7}")
        for k, r in rec:
            if k == 'it':
                if r['it'] < lo:
                    continue
                if r['it'] % args.every and r is not its[-1]:
                    continue
                b = r['pbest']
                best = (f"k{b[0]}:{float(b[1]):.4f}>{b[2]}(+{float(b[3]):.4f})/{b[4]}"
                        if b else '')
                print(f"{r['it']:>5} {r['LP']:>12.6f} "
                      f"{r['sup'][0] if r['sup'] else '':>6} "
                      f"{r['master'][1] if r['master'] else '':>7} "
                      f"{r['cq'][0] if r['cq'] else '':>6} "
                      f"{r['pg'][0] if r['pg'] else '-':>5} "
                      f"{r['M']:>10.6f} {r['kmax']:>9.6f} {best:>24} "
                      f"{r['secs'][0] if r['secs'] else '':>7}")
            elif k == 'lat':
                (tag, it, npos, poses, cols, gap, nmem, pmem) = r
                if int(it) < lo:
                    continue
                print(f"  >>> lattice after it {it}: +{npos} poses -> {poses} poses / {cols} "
                      f"columns, pricing gap {float(gap):+.6f}, +{nmem} clique / "
                      f"+{pmem or 0} polygon memberships")
            elif k == 'pgfin':
                n, ok, worst, tight, red = r
                print(f"  FINAL: {n} polygon rows {ok} (worst violation {worst}, {tight} tight, "
                      f"all re-derived: {red})")
            elif k == 'fin':
                mass, M, Mok, Kq, Kok, rok = r
                print(f"  FINAL EXACT: mass {mass}  M {M} {Mok}  max clique {Kq} {Kok}  "
                      f"regions {rok}")
            elif k == 'res':
                print(f"  RESULT {r[0]}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
