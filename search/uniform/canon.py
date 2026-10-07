#!/usr/bin/env python3
"""Canonicalise a certificate in place: container side in lowest terms, D reduced by the gcd of D
and all coordinates, points sorted.  Pure integer arithmetic; the asserted statement is unchanged.

    python3 search/uniform/canon.py CERT [CERT ...]
"""
import math, sys
from fractions import Fraction
from pathlib import Path

for path in sys.argv[1:]:
    t = Path(path).read_text().split()
    sn, sd, D, W, m = map(int, t[:5])
    rows = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    s = Fraction(sn, sd)
    g = math.gcd(D, *(v for r in rows for v in r[:2]))
    D //= g; rows = sorted((x // g, y // g, w) for x, y, w in rows)
    assert (s * D).denominator == 1
    with open(path, "w") as f:
        f.write(f"{s.numerator} {s.denominator}\n{D}\n{W}\n{len(rows)}\n")
        for x, y, w in rows: f.write(f"{x} {y} {w}\n")
    print(f"{path}: s={s} D={D} m={m} k={W}")
