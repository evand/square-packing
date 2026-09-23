#!/usr/bin/env python3
"""Lift D4-symmetrised pose sets from a smaller container to a larger one, and merge them.

A pose `(cx, cy, theta)` of a closed unit square inside `[0,t]^2` is also a pose inside `[0,t']^2`
for any `t' >= t` after the *centred* translation `(cx, cy) -> (cx + (t'-t)/2, cy + (t'-t)/2)`,
which commutes with the dihedral group of the container (it maps the centre to the centre), so a
D4-symmetrised measure lifts to a D4-symmetrised measure of the same mass and the same maximum
coverage.  This lets the pose sets found at `3.82` be reused at `3.85` and `3.87` — only the poses
matter, since `dual_exact.py build` re-chooses all the masses.

    python3 search/n11_lift.py --t 77/20 --out runs/n11_union_385.txt runs/dual_N385_support.txt \
                                                                     runs/dual_N382_support.txt
"""
import argparse
import math
import os
from fractions import Fraction as Fr


def read(path):
    t = None
    out = []
    for line in open(path):
        s = line.strip()
        if s.startswith('#'):
            if t is None and 't=' in s:
                t = float(Fr(s.split('t=')[1].split()[0]))
            elif t is None and 't =' in s:
                t = float(Fr(s.split('t =')[1].split()[0]))
            continue
        q = s.split('#')[0].split()
        if not q or q[0] != 'pose':
            continue
        if len(q) == 5:
            out.append((float(q[1]), float(q[2]), float(q[3]), float(q[4])))
        elif len(q) == 6:
            p, qq = int(q[1]), int(q[2])
            out.append((float(Fr(q[3])), float(Fr(q[4])), math.degrees(2.0 * math.atan2(p, qq)), float(Fr(q[5]))))
    if t is None:
        raise SystemExit(f'{path}: no container side in the header')
    return t, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--t', required=True, help='target container side (exact fraction)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--eps', type=float, default=2e-4, help='dedupe tolerance on (cx, cy, theta_deg)')
    a = ap.parse_args()
    T = float(Fr(a.t))
    seen = {}
    for f in a.files:
        t, poses = read(f)
        if t > T + 1e-12:
            raise SystemExit(f'{f}: container {t} > target {T}; a packing does not shrink')
        d = (T - t) / 2.0
        n = 0
        for cx, cy, th, mu in poses:
            k = (round((cx + d) / a.eps), round((cy + d) / a.eps), round(th / a.eps))
            if k in seen:
                continue
            seen[k] = (cx + d, cy + d, th, mu)
            n += 1
        print(f"  {f}: t={t} -> {T} (shift {d:+.6f}), {len(poses)} poses, {n} new")
    with open(a.out, 'w') as fh:
        fh.write(f"# t={T} lifted union of {' '.join(os.path.basename(x) for x in a.files)}\n")
        fh.write("# D4-symmetrised measure: pose cx cy theta_deg mu  (masses are re-chosen downstream)\n")
        for cx, cy, th, mu in seen.values():
            fh.write(f"pose {cx:.13f} {cy:.13f} {th:.11f} {max(mu, 1e-9):.12f}\n")
    print(f"wrote {a.out}: {len(seen)} poses at t={T}")


if __name__ == '__main__':
    main()
