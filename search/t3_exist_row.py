#!/usr/bin/env python3
"""t3_exist_row: stress-test of the candidate existence clause (E3-row) over MANY margin-0
packings at each dangerous angle vector.

Task tasks/t3-existence/README.md.  Nothing in search/ is modified; t3_chain is imported.

(E3-row)  Every packing of six unit squares in [0,3]^2 contains squares s1, s2, s3 and an axis ax
          such that both consecutive pairs are separated in the EXACT ax direction (a pair row
          with normal exactly e_x or e_y, which is available when one endpoint of the link is
          exactly axis-parallel).  Lemma Z'' then gives delta <= 0.

At a theta with delta*(theta) = 0 the optimal set is usually positive-dimensional, so a single
optimum per theta (which is all runs/t3_chain_scan1.jsonl and runs/t3_exist_scope1.jsonl have)
is a weak test.  Here every distinct optimum the multistart produces is tested.

usage:
    python3 search/t3_exist_row.py test [--reps 20] [--out runs/t3_exist_row1.txt]
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3_chain as TC                                             # noqa: E402
import t3_exist_zero as EZ                                        # noqa: E402

N, T = 6, 3
D = math.pi / 180.0


def families():
    """angle vectors known (or likely) to admit a packing, hardest first."""
    F = []
    for a in (45, 44, 40, 35, 30, 25, 20, 15, 10, 5, 2, 1, 0.5):
        F.append((f'two at {a} deg', [0, 0, 0, 0, a * D, a * D]))
        F.append((f'two at +-{a} deg', [0, 0, 0, 0, a * D, -a * D]))
    for a in (45, 30, 20, 10, 5, 1):
        F.append((f'one at {a} deg', [0, 0, 0, 0, 0, a * D]))
        F.append((f'three at {a} deg', [0, 0, 0, a * D, a * D, a * D]))
    F.append(('all axis-parallel', [0] * 6))
    F.append(('mixed 43/6.7/0.65', [0, 0, 0, 43.07 * D, 6.70 * D, 0.65 * D]))
    F.append(('mixed 30.9/33.4/0.3', [0, 0, 0, 0.31 * D, 30.91 * D, 33.43 * D]))
    F.append(('mixed 45/30', [0, 0, 0, 0, 45 * D, 30 * D]))
    F.append(('mixed 45/20', [0, 0, 0, 0, 45 * D, 20 * D]))
    F.append(('mixed 45/5', [0, 0, 0, 0, 45 * D, 5 * D]))
    return F


def distinct(zs, tol=1e-3):
    out = []
    for z in zs:
        c = np.array(list(z[1::3]) + list(z[2::3]))
        if all(np.abs(c - o).max() > tol for o in out):
            out.append(c)
            yield z


def test_one(theta, rng, reps):
    v, z0, tops = TC.delta_star(list(theta), T, rng, limit=300, extra=400, keep=reps)
    recs = []
    for z in distinct(tops):
        lvl = float(z[0])
        if lvl < -1e-9:
            continue
        A, b, tags = TC.all_rows(z, N, T, lvl)
        TH = np.array(z[3::3])
        zc = EZ.zchain(z, tags, TH, lvl)
        bh = TC.best_h(z, N, T, lvl, modes=('C', 'H'), rows=(A, b, tags))
        F, _ = TC.dual_bound(A, b, range(len(b)))
        recs.append(dict(delta=lvl, zchain=(None if zc is None else zc[0]),
                         zpath=(None if zc is None else [zc[1], zc[2]]),
                         C=bh['C']['bound'], H=bh['H']['bound'], F=F,
                         nchains=bh['nchains'],
                         nzero=int(sum(1 for t in theta if EZ.tilt(t) < 1e-12)),
                         z=[float(x) for x in z]))
    return v, recs


def hunt_families(rng, n):
    """angle vectors with only 2 or 3 exact zeros -- the thinnest part of the scope."""
    F = []
    D3 = [(43.07, 6.70, 0.65), (0.31, 30.91, 33.43), (29.60, 2.48, 10.19), (1.10, 40.67, 3.31),
          (45, 20, 2), (45, 10, 1), (40, 30, 0.5), (35, 5, 0.2), (44, 2, 0.1), (20, 10, 0.5)]
    for tri in D3:
        F.append((f'3 zeros {tri}', [0, 0, 0] + [x * D for x in tri]))
    for _ in range(n):
        tri = sorted(rng.uniform(0.1, 45, 3), reverse=True)
        F.append((f'3 zeros rand {tuple(round(x,2) for x in tri)}',
                  [0, 0, 0] + [x * D for x in tri]))
    for _ in range(n):
        q = sorted(rng.uniform(0.1, 45, 4), reverse=True)
        F.append((f'2 zeros rand {tuple(round(x,2) for x in q)}',
                  [0, 0] + [x * D for x in q]))
    return F


def cmd_hunt(a):
    rng = np.random.default_rng(a.seed)
    tot = bad = 0
    minap = 6
    with open(a.out, 'w') as f:
        print(f"{'family':>42} {'delta*':>11} {'#pack':>6} {'row':>8} {'#axis-par':>10}")
        for (name, th) in hunt_families(rng, a.nrand):
            v, recs = test_one(th, rng, a.reps)
            ok = sum(1 for r in recs if r['zchain'] is not None)
            tot += len(recs)
            bad += len(recs) - ok
            aps = [int(sum(1 for t in r['z'][3::3] if EZ.tilt(t) < 1e-12)) for r in recs]
            if aps:
                minap = min(minap, min(aps))
            for r in recs:
                r['fam'] = name
                f.write(json.dumps(r) + '\n')
            f.flush()
            print(f"{name:>42} {v:+11.3e} {len(recs):>6} {ok:>3}/{len(recs):<4} "
                  f"{(min(aps) if aps else -1):>10}", flush=True)
    print(f"\npackings tested: {tot};  (E3-row) fails at {bad};  "
          f"fewest exactly-axis-parallel squares seen: {minap}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    t = sub.add_parser('test')
    t.add_argument('--reps', type=int, default=20)
    t.add_argument('--seed', type=int, default=7)
    t.add_argument('--out', default='runs/t3_exist_row1.jsonl')
    h = sub.add_parser('hunt')
    h.add_argument('--reps', type=int, default=40)
    h.add_argument('--nrand', type=int, default=12)
    h.add_argument('--seed', type=int, default=11)
    h.add_argument('--out', default='runs/t3_exist_hunt1.jsonl')
    a = ap.parse_args()
    if a.cmd == 'hunt':
        cmd_hunt(a)
        return
    rng = np.random.default_rng(a.seed)
    tot = bad = 0
    with open(a.out, 'w') as f:
        print(f"{'family':>22} {'delta*':>11} {'#packings':>10} {'row exists':>11} "
              f"{'max row val':>12} {'max C':>10}")
        for (name, th) in families():
            v, recs = test_one(th, rng, a.reps)
            ok = sum(1 for r in recs if r['zchain'] is not None)
            tot += len(recs)
            bad += len(recs) - ok
            for r in recs:
                r['fam'] = name
                f.write(json.dumps(r) + '\n')
            f.flush()
            mv = max((r['zchain'] for r in recs if r['zchain'] is not None), default=float('nan'))
            mc = max((r['C'] for r in recs), default=float('nan'))
            print(f"{name:>22} {v:+11.3e} {len(recs):>10} {ok:>6}/{len(recs):<4} "
                  f"{mv:12.2e} {mc:10.2e}", flush=True)
    print(f"\nmargin-0 packings tested: {tot};  (E3-row) fails at {bad}")


if __name__ == '__main__':
    main()
