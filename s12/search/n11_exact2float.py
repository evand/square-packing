#!/usr/bin/env python3
"""n11_exact2float.py SRC DST T -- an exact support (pose p q cx cy mass) as a float warm-start file for clique_ceiling.py --warm (pose cx cy theta_deg mass, with a "# t=" header)."""
import sys, math
from fractions import Fraction as Fr
src, dst, t = sys.argv[1], sys.argv[2], sys.argv[3]
out = [f"# t={float(Fr(t))} converted from {src} (pose cx cy theta_deg mass)"]
for line in open(src):
    q = line.split()
    if q and q[0] == 'pose':
        p, qq = int(q[1]), int(q[2]); cx, cy = float(Fr(q[3])), float(Fr(q[4]))
        th = math.degrees(2*math.atan2(p, qq))
        out.append(f"pose {cx:.9f} {cy:.9f} {th:.6f} {float(Fr(q[5])):.6f}")
open(dst,'w').write('\n'.join(out)+'\n'); print('\n'.join(out))
