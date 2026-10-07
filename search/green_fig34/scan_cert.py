#!/usr/bin/env python3
"""Independent float cross-check of a unit-weight certificate: smax (largest point-free square at
any angle inside the container, smax.py) must be < 1.  Shares no code with zeromargin.py.

    python3 scan_cert.py CERT [nth]
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from smax import smax

tok = open(sys.argv[1]).read().split()
sn, sd, D, W, n = map(int, tok[:5])
P = [(int(tok[5 + 3 * i]) / D, int(tok[6 + 3 * i]) / D) for i in range(n)]
t = sn / sd
nth = int(sys.argv[2]) if len(sys.argv) > 2 else 720
v, th, _ = smax(P, t, nth=nth, nloc=6)
print(f"{sys.argv[1]}: side {t:.7f}, {n} points, smax = {v:.9f} at {math.degrees(th):.3f} deg "
      f"-> {'OK (< 1)' if v < 1 else 'FAILS'}")
