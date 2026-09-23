#!/usr/bin/env python3
"""t3_exist_scope: how many squares of a REAL packing of six unit squares in [0,3]^2 can be
strongly tilted -- i.e. what the existence clause (E3) actually has to cover.

Task tasks/t3-existence/README.md.  Nothing in search/ is modified; t3_chain is imported.

For each (k, tau): k angles drawn uniformly from [tau, 90 deg - tau] (tilt >= tau), the other
6 - k set to exactly 0 (the most favourable choice for fitting), and delta*(theta) measured by
the multistart LP ascent of t3_chain.delta_star.  delta* is a FEASIBLE value, so
"delta* = 0 was reached" is a reliable positive statement: a packing with k squares of tilt
>= tau exists.  "max delta* < 0 over the samples" is [measured], not a proof.

usage:
    python3 search/t3_exist_scope.py scan --out runs/t3_exist_scope1.jsonl [--reps 30]
    python3 search/t3_exist_scope.py report runs/t3_exist_scope1.jsonl
"""
import argparse
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3_chain as TC                                             # noqa: E402

N, T = 6, 3
HALF = math.pi / 2


def tiltd(t):
    r = abs(t) % HALF
    return math.degrees(min(r, HALF - r))


def job(arg):
    k, tau, rep, seed = arg
    rng = np.random.default_rng(seed)
    th = np.zeros(N)
    lo, hi = math.radians(tau), HALF - math.radians(tau)
    if hi < lo:
        lo = hi = math.radians(45.0)
    idx = rng.permutation(N)[:k]
    th[idx] = rng.uniform(lo, hi, k)
    v, z, _ = TC.delta_star(list(th), T, rng, limit=250, extra=350)
    rec = dict(k=k, tau=tau, rep=rep, theta=[float(x) for x in th], delta=float(v),
               tilts=[tiltd(x) for x in th])
    if v > -1e-9:
        rec['z'] = [float(x) for x in z]
        A, b, tags = TC.all_rows(z, N, T, v)
        bh = TC.best_h(z, N, T, v, modes=('C', 'CW', 'H'), rows=(A, b, tags))
        rec['C'] = bh['C']['bound']
        rec['CW'] = bh['CW']['bound']
        rec['H'] = bh['H']['bound']
        rec['nchains'] = bh['nchains']
        rec['Cpath'] = bh['C']['path']
        rec['Caxis'] = bh['C']['axis']
    return rec


def cmd_scan(a):
    args = []
    s = a.seed
    for k in (2, 3, 4, 5, 6):
        for tau in (10, 15, 20, 25, 30, 35, 40, 45):
            for rep in range(a.reps):
                s += 1
                args.append((k, tau, rep, s))
    with Pool(a.nproc) as p, open(a.out, 'w') as f:
        for i, rec in enumerate(p.imap_unordered(job, args, chunksize=1)):
            f.write(json.dumps(rec) + '\n')
            f.flush()
            if rec['delta'] > -1e-9:
                print(f"[{i+1}/{len(args)}] PACKING k={rec['k']} tau={rec['tau']} "
                      f"tilts={[round(x, 1) for x in rec['tilts']]} "
                      f"C={rec.get('C', float('nan')):+.5f} H={rec.get('H', float('nan')):+.5f}",
                      flush=True)


def cmd_report(path):
    rows = [json.loads(l) for l in open(path)]
    cells = {}
    for r in rows:
        cells.setdefault((r['k'], r['tau']), []).append(r)
    taus = sorted({t for (_, t) in cells})
    print("max delta* (and #packings found / #samples) by (k tilted squares, tilt floor tau):\n")
    print("  k \\ tau " + "".join(f"{t:>16}" for t in taus))
    for k in sorted({kk for (kk, _) in cells}):
        line = f"  {k:>7} "
        for t in taus:
            rs = cells.get((k, t), [])
            if not rs:
                line += f"{'-':>16}"
                continue
            mx = max(r['delta'] for r in rs)
            nz = sum(1 for r in rs if r['delta'] > -1e-9)
            line += f"{mx:>+10.2e}{nz:>3}/{len(rs):<3}"
        print(line)
    print("\ncertificates at the packings found (delta = 0):")
    ok = [r for r in rows if r['delta'] > -1e-9]
    print(f"  packings: {len(ok)}   C<=0: {sum(1 for r in ok if r.get('C', 1) <= 1e-9)}"
          f"   CW<=0: {sum(1 for r in ok if r.get('CW', 1) <= 1e-9)}"
          f"   H<=0: {sum(1 for r in ok if r.get('H', 1) <= 1e-9)}"
          f"   no chain of three: {sum(1 for r in ok if r.get('nchains', 0) == 0)}")
    ok.sort(key=lambda r: -sorted(r['tilts'], reverse=True)[2])
    print("\n  the most-tilted packings found (by 3rd largest tilt):")
    for r in ok[:10]:
        print(f"    tilts={[round(x, 2) for x in sorted(r['tilts'], reverse=True)]} "
              f"C={r.get('C', float('nan')):+.5f} H={r.get('H', float('nan')):+.5f} "
              f"nchains={r.get('nchains')}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--out', default='runs/t3_exist_scope1.jsonl')
    s.add_argument('--reps', type=int, default=24)
    s.add_argument('--nproc', type=int, default=2)
    s.add_argument('--seed', type=int, default=1000)
    r = sub.add_parser('report')
    r.add_argument('path')
    a = ap.parse_args()
    if a.cmd == 'scan':
        cmd_scan(a)
    else:
        cmd_report(a.path)


if __name__ == '__main__':
    main()
