#!/usr/bin/env python3
"""Upper-bound certificate from an f64 packing (no KKT solve): spread the centres by lam = 1 + eps about the origin
(opens a gap ~ eps * distance at every contact; squares keep their size), round centres and t = tan(theta/2) to
rationals, side S' = lam * S + 1e-12.  The result is checked exactly by ../exact/verify_cert.py; if valid, s(n) <= S'.

  mkcert.py in.txt out.cert [--eps 1e-6]
"""
import argparse, math
from fractions import Fraction as Fr
import mcmin

ap = argparse.ArgumentParser()
ap.add_argument('inp'); ap.add_argument('out')
ap.add_argument('--eps', type=float, default=1e-6)
a = ap.parse_args()
s, sq = mcmin.load_deg(a.inp)
lam = Fr(1) + Fr(a.eps).limit_denominator(10 ** 12)
S = (lam * Fr(s)).limit_denominator(10 ** 15) + Fr(1, 10 ** 12)
with open(a.out, 'w') as f:
    f.write(f'# upper-bound certificate from {a.inp} (centres x {lam}); check with verify_cert.py\n')
    f.write(f'{len(sq)} {S}\n')
    for x, y, d in sq:
        t = Fr(math.tan(math.radians(d % 90) / 2)).limit_denominator(10 ** 15)
        f.write(f'{(lam * Fr(x)).limit_denominator(10 ** 15)} {(lam * Fr(y)).limit_denominator(10 ** 15)} {t}\n')
print(f'S\' = {float(S)!r}  ({S})')
