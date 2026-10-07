#!/usr/bin/env python3
"""Write unit-weight certificates (certificates/FORMAT.md) for the Figure-34 sets.

    python3 mkcert.py 17 TNUM TDEN OUT [--set NAME]
    python3 mkcert.py 19 TNUM TDEN OUT [--u U]

The exact (irrational) set at its nominal side t is scaled by lambda = t'/t to the
rational side t' = TNUM/TDEN and every coordinate is rounded to the grid 1/(TDEN*k).  Rounding is
not assumed sound: the certificate is checked as it stands by the exact checkers.
"""
import sys, os, math, argparse
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family import green17, T17
from family19 import friedman19

R2 = math.sqrt(2)
R65 = math.sqrt(65)

SETS17 = {
    # Green's set as derived from the bound (GREEN_FIG34.md sec. 2): s = 1, b = sqrt2 - 1/2,
    # V = (2 sqrt2 + 12)/17, a = sqrt2 - V/2, e = (8 sqrt2 - 3)/17.
    'green': lambda: (T17, green17(T17, (16 * R2 - 6) / 17, R2 - 0.5, 1.0, (2 * R2 + 12) / 17)),
    # the family optimum (GREEN_FIG34.md sec. 3): s = 8/sqrt65, e = 4/sqrt65, V = 7/sqrt65,
    # t = 2 sqrt2 + 13/sqrt65
    'hstar': lambda: (2 * R2 + 13 / R65, green17(2 * R2 + 13 / R65, R2 - 3.5 / R65, R2 - 4 / R65,
                                                  8 / R65, 7 / R65)),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int)
    ap.add_argument('tnum', type=int); ap.add_argument('tden', type=int)
    ap.add_argument('out')
    ap.add_argument('--set', default='hstar', choices=sorted(SETS17))
    ap.add_argument('--u', type=float, default=None)
    ap.add_argument('-k', type=int, default=1, help='coordinate grid 1/(TDEN*k)')
    ap.add_argument('--params', type=str, default=None,
                    help='17 only: a,b,s,v at side T17 (floats), overrides --set')
    a = ap.parse_args()
    tp = F(a.tnum, a.tden)
    if a.n == 17:
        if a.params:
            t = T17
            pa, pb, ps, pv = map(float, a.params.split(','))
            P = green17(t, pa, pb, ps, pv)
        else:
            t, P = SETS17[a.set]()
    else:
        t = 6 * R2 - 4
        d = 2 * R2 - 2
        u = a.u if a.u is not None else (d + math.sqrt(1 - d * d) / 2) / 2
        P = friedman19(t, u)
    lam = float(tp) / t
    D = a.tden * a.k
    rows = []
    for (x, y) in P:
        X = round(x * lam * D); Y = round(y * lam * D)
        rows.append((X, Y))
    with open(a.out, 'w') as f:
        f.write(f"{a.tnum} {a.tden}\n{D}\n1\n{len(rows)}\n")
        for X, Y in rows:
            f.write(f"{X} {Y} 1\n")
    print(f"wrote {a.out}: {len(rows)} points, side {a.tnum}/{a.tden} = {float(tp):.9f}, "
          f"lambda = {lam:.9f} (t = {t:.9f})")


if __name__ == '__main__':
    main()
