#!/usr/bin/env python3
"""Angle homotopy: rotate chosen squares from their angles in A to their angles in B in steps; at each step those angles are
frozen and everything else (positions incl. theirs, other angles, side) re-minimised by slp2.  Prints s(t)."""
import math, sys
from contacts import load, align
from slp2 import slp2, repair

A, B = sys.argv[1], sys.argv[2]
ids = [int(x) for x in sys.argv[3].split(',')]
steps = int(sys.argv[4]) if len(sys.argv) > 4 else 10
sA, SA = load(A); sB, SB = load(B); SB = align(SB, SA)
s, sq = sA, list(SA)
print(f't=0.00 s={s:.9f}')
for k in range(1, steps + 1):
    t = k / steps
    sq = [(x, y, (SA[i][2] * (1 - t) + SB[i][2] * t) if i in ids else a) for i, (x, y, a) in enumerate(sq)]
    s, sq = repair(s, sq)
    fixed = [3 * i + 2 for i in ids]
    s, sq, ds0 = slp2(s, sq, R=1e-3, rmin=1e-9, fixed=fixed)
    print(f't={t:.2f} s={s:.9f} angles={[round(math.degrees(sq[i][2]) % 90, 4) for i in ids]} jammed={ds0 is not None and ds0 > -1e-12}', flush=True)
# release
s, sq, ds0 = slp2(s, sq, R=1e-3, rmin=1e-9)
print(f'released: s={s:.9f}')
