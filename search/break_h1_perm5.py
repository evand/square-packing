#!/usr/bin/env python3
"""break_h1_perm5: the T = 5, n = 20, k = 15 top-of-hole chain-free family WITH a permutation
hole set -- the true analogue of the T = 4 `eps^3` cell.

`BREAK_H1.md` §2.1: with the merge-type level map, `k = (T-1)^2` at `T = 5` cannot have a
permutation hole set, but `k = (T-1)^2 - 1 = 15` can, and this file writes out the unique
combinatorial shape that does:

  merged cols {a, a+1}, merged rows {b, b+1};  C, R the three non-merged cols / rows.
  central level cell     4 tiling cells, a pinwheel at tilt eps: 1 cut-set rep, 3 far
  (a/a+1, r) for r in R  1 cut rep + 1 removed; two rows' removed are HOLES in DIFFERENT
                         columns, the third row r* has its removed cell FAR
  (c, b/b+1) for c in C  1 cut rep + 1 removed; two cols' removed are HOLES in DIFFERENT rows,
                         the third col c1 has its removed cell FAR
  (c, r), c in C, r in R the level cell (c1, r*) is DROPPED -- that tiling cell is the fifth
                         HOLE; the other eight carry axis-parallel cut-set squares
  => cut 15, far 5, holes 5, one hole per tiling row and per tiling column.

Imports only; nothing in search/ is modified."""
import itertools
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import break_h1 as B                                             # noqa: E402

T, n, K = 5.0, 20, 15


def families():
    Ti = int(round(T))
    out = []
    for a, bq in itertools.product(range(Ti - 1), repeat=2):
        mcol, mrow = (a, a + 1), (bq, bq + 1)
        C = [c for c in range(Ti) if c not in mcol]
        R = [r for r in range(Ti) if r not in mrow]
        Lx = {c: (c if c <= a else (a if c == a + 1 else c - 1)) for c in range(Ti)}
        My = {r: (r if r <= bq else (bq if r == bq + 1 else r - 1)) for r in range(Ti)}
        for pin_rep in range(4):
            pin = [(mcol[0], mrow[0]), (mcol[1], mrow[0]),
                   (mcol[0], mrow[1]), (mcol[1], mrow[1])]
            for rstar in R:
                Rh = [r for r in R if r != rstar]
                for assign_r in itertools.permutations(mcol):       # hole column per hole row
                    for c1 in C:
                        Ch = [c for c in C if c != c1]
                        for assign_c in itertools.permutations(mrow):
                            for far_rcol in mcol:                   # which col of r* is far
                                for far_crow in mrow:               # which row of c1 is far
                                    cut, far, holes = [], [], []
                                    cut.append(pin[pin_rep])
                                    far += [p for q, p in enumerate(pin) if q != pin_rep]
                                    for i, r in enumerate(Rh):
                                        hc = assign_r[i]
                                        holes.append((hc, r))
                                        cut.append(((mcol[1] if hc == mcol[0] else mcol[0]), r))
                                    far.append((far_rcol, rstar))
                                    cut.append(((mcol[1] if far_rcol == mcol[0]
                                                 else mcol[0]), rstar))
                                    for i, c in enumerate(Ch):
                                        hr = assign_c[i]
                                        holes.append((c, hr))
                                        cut.append((c, (mrow[1] if hr == mrow[0]
                                                        else mrow[0])))
                                    far.append((c1, far_crow))
                                    cut.append((c1, (mrow[1] if far_crow == mrow[0]
                                                     else mrow[0])))
                                    holes.append((c1, rstar))
                                    for c in C:
                                        for r in R:
                                            if (c, r) != (c1, rstar):
                                                cut.append((c, r))
                                    if len(cut) != K or len(far) != 5 or len(holes) != 5:
                                        continue
                                    cells = cut + far
                                    if len(set(cells + holes)) != 25:
                                        continue
                                    if (len({h[0] for h in holes}) != 5
                                            or len({h[1] for h in holes}) != 5):
                                        continue
                                    pat = [(Lx[c], My[r]) for (c, r) in cut]
                                    if len(set(pat)) != K:
                                        continue
                                    out.append((cut, far, [list(p) for p in pat], set(pin),
                                                f'perm5 m{a}{bq} p{pin_rep} r*{rstar} c1{c1}'))
    return out


def main():
    fams = families()
    seen, ded = set(), []
    for f in fams:
        key = (tuple(sorted(f[0])), tuple(sorted(f[1])))
        if key in seen:
            continue
        seen.add(key); ded.append(f)
    print(f'# {len(fams)} shapes, {len(ded)} distinct cell/far assignments', flush=True)
    rng = np.random.default_rng(97)
    nsamp = int(os.environ.get('NSAMP', 250))
    if len(ded) > nsamp:
        ded = [ded[q] for q in rng.permutation(len(ded))[:nsamp]]
        print(f'# sampled down to {len(ded)}', flush=True)
    out = open('runs/break_h1_t5_perm.jsonl', 'w')
    carry = None
    for epsdeg in (10.0, 5.0, 2.0, 1.0, 0.5, 0.25):
        eps = math.radians(epsdeg)
        best = None
        nhit = 0
        for (cut, far, pat, pin, lab) in ded:
            cells = cut + far
            for sgn in (1.0, -1.0):
                th = [(sgn * eps if c in pin or c in set(far) else 0.0) for c in cells]
                X = np.array([c + 0.5 for (c, r) in cells], dtype=float)
                Y = np.array([r + 0.5 for (c, r) in cells], dtype=float)
                seeds = [(X, Y)]
                for _ in range(8):
                    seeds.append((X + rng.normal(0, 0.05, n), Y + rng.normal(0, 0.05, n)))
                if carry is not None:
                    seeds.append((carry[0].copy(), carry[1].copy()))
                v, z, hit = B.delta_star_at(th, T, rng, seeds=seeds, limit=0, jitters=0,
                                            extra=0, pattern=[tuple(p) for p in pat], k=K)
                if z is None or not B.chainfree_ok(z, n, T, v, range(K)):
                    continue
                nhit += 1
                if best is None or v > best[0]:
                    best = (v, z, lab, pat)
        if best is None:
            print(f'  eps={epsdeg:g}: nothing feasible', flush=True)
            continue
        v, z, lab, pat = best
        carry = (z[1::3].copy(), z[2::3].copy())
        print(f'  eps={epsdeg:6g}  delta* = {v:+.8e}   /e {-v/eps:.5f}  /e^2 {-v/eps**2:.5f}  '
              f'/e^3 {-v/eps**3:.4f}  /e^4 {-v/eps**4:.3f}   feas {nhit}   {lab}', flush=True)
        out.write(json.dumps(dict(lab=lab, k=K, eps=epsdeg, delta=float(v), cfree=True,
                                  pattern=pat, z=[float(q) for q in z])) + '\n')
        out.flush()
    out.close()


if __name__ == '__main__':
    main()
