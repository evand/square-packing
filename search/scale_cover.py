#!/usr/bin/env python3
"""Scale every weight of a certificate-format cover by an exact rational factor lambda = p/q > 1.

Rationale (task H, coordinator follow-up 2026-08-30): if a closed cover has true minimum captured
weight mu > 0 at every admissible pose, then multiplying every point's weight by 1/mu gives a
valid closed cover of total (original total)/mu -- captured weight scales linearly with the
weights, so scaling by any lambda with lambda*mu >= 1 (lambda >= 1/mu) is again a valid cover.
Since the true global minimum mu is only known approximately (from a stress test / exact spot
checks), scaling by a safety factor above the reciprocal of the WORST FOUND captured weight gives
a candidate cover that the exact zeromargin.py checker can then confirm or refute exhaustively.

Exact scaling: weight_i/W becomes (weight_i * p)/(W * q) -- multiply every integer weight
numerator by p and the weight denominator by q; no rounding, no floats.

Usage: python3 search/scale_cover.py IN.txt P Q OUT.txt
   (scales every weight by the exact fraction P/Q, e.g. 51 50 for x1.02, 103 100 for x1.03)
"""
import sys

def main():
    path, p, q, outpath = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    tok = open(path).read().split()
    sn, sd, D, W, n = map(int, tok[:5])
    lines = []
    tot_before = 0
    tot_after = 0
    for i in range(n):
        X, Y, w = map(int, tok[5 + 3 * i: 8 + 3 * i])
        tot_before += w
        w2 = w * p
        lines.append((X, Y, w2))
        tot_after += w2
    Wn = W * q
    with open(outpath, 'w') as f:
        f.write(f"{sn} {sd}\n{D}\n{Wn}\n{n}\n")
        for X, Y, w2 in lines:
            f.write(f"{X} {Y} {w2}\n")
    print(f"scaled {path} by {p}/{q}: total {tot_before}/{W} = {tot_before/W:.6f}  ->  "
          f"{tot_after}/{Wn} = {tot_after/Wn:.6f}")
    print(f"written to {outpath}")

if __name__ == '__main__':
    main()
