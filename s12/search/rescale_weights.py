#!/usr/bin/env python3
"""Uniformly rescale a certificate's weights so that a verifier minimum a (a fraction) becomes >= 1.
w'_i = ceil(w_i * Wout / (Win * a)) over Wout: every placement captured >= a (old weights) now captures >= 1.
usage: python3 search/rescale_weights.py IN OUT A_NUM A_DEN [Wout]
(A_NUM/A_DEN: the exact minimum reported by `verify` for IN; zero-weight points are dropped)"""
import sys
from fractions import Fraction
inp, out, an, ad = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
Wout = int(sys.argv[5]) if len(sys.argv) > 5 else 10**7
a = Fraction(an, ad)
tok = open(inp).read().split()
sn, sd, D, Win, m = map(int, tok[:5]); rest = list(map(int, tok[5:]))
assert len(rest) == 3 * m
lines = []; tot = 0
for i in range(m):
    X, Y, w = rest[3*i:3*i+3]
    if w == 0: continue
    q = Fraction(w * Wout, Win) / a
    wn = -(-q.numerator // q.denominator)
    lines.append(f"{X} {Y} {wn}"); tot += wn
with open(out, 'w') as f:
    f.write(f"{sn} {sd}\n{D}\n{Wout}\n{len(lines)}\n" + "\n".join(lines) + "\n")
print(f"points {len(lines)} total {Fraction(tot, Wout)} = {tot / Wout:.7f}")
