#!/usr/bin/env python3
"""Scale a finished certificate up to its critical container size.

This is the step that took the bound from 3.92 to 3.931795, and it is pure
arithmetic -- no LP, no new points.

The idea
--------
A certificate is a set of integer points (X_i, Y_i) with denominator D, asserting a
covering property for the container [0, s]^2 where the unit square has side 1.
Scaling the whole picture (points and container together) by lambda is, in this file
format, just::

    D  ->  D / lambda            s  ->  lambda * s

so the product ``K = s * D`` is invariant, and the integer coordinates never change.
K is the container side measured in grid units; D is how many grid units make one
unit square.  Shrinking D therefore blows the container up relative to the unit
square, which is exactly the direction that strengthens the claim -- and eventually
breaks it, when some placement's captured weight falls below 1.

So the certificate does not prove one bound; it proves a *family* of bounds, and the
best of them is at the critical D.  We binary-search for it.

Because the weights are untouched, the total-weight side of the argument (sum w < n)
is unaffected by scaling; only the covering property can fail.

Usage
-----
    python3 search/scale_to_critical.py certificates/s12_lower_3.931795.txt

    # wider bracket, coarser angle grid for speed while bracketing:
    python3 search/scale_to_critical.py CERT --n 12 --N 2000 --lo 390000 --hi 400000

Reproducing the published number
--------------------------------
Starting from the raw LP output at container 3.92 and searching down from D = 400000
recovers D = 398800, i.e. s = 1568000/398800 = 3920/997 = 3.931795386 -- the
certificate shipped in ``certificates/s12_lower_3.931795.txt`` (which is the same
point set written at D = 1994, every coordinate having been divisible by 200).
"""

import argparse
import math
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_VERIFIER = REPO / "verify" / "target" / "release" / "verify"


def read_cert(path):
    """Return (s: Fraction, D: int, W: int, rows: list[(x, y, w)])."""
    toks = Path(path).read_text().split()
    it = iter(toks)
    s_num = int(next(it))
    s_den = int(next(it))
    D = int(next(it))
    W = int(next(it))
    m = int(next(it))
    rows = [(int(next(it)), int(next(it)), int(next(it))) for _ in range(m)]
    return Fraction(s_num, s_den), D, W, rows


def write_cert(path, s: Fraction, D, W, rows):
    with open(path, "w") as f:
        f.write(f"{s.numerator} {s.denominator}\n{D}\n{W}\n{len(rows)}\n")
        for x, y, w in rows:
            f.write(f"{x} {y} {w}\n")


def verifies(verifier, cert_path, n, N, threads):
    """True iff the verifier accepts. Note 'NOT VERIFIED' contains 'VERIFIED'."""
    out = subprocess.run(
        [str(verifier), str(cert_path), str(n), str(N), str(threads), "0"],
        capture_output=True, text=True,
    )
    if out.returncode != 0:
        raise RuntimeError(f"verifier failed:\n{out.stdout}\n{out.stderr}")
    if "NOT VERIFIED" in out.stdout:
        return False
    if "VERIFIED" in out.stdout:
        return True
    raise RuntimeError(f"could not parse verifier output:\n{out.stdout}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cert")
    ap.add_argument("--n", type=int, default=12, help="number of squares to rule out")
    ap.add_argument("--N", type=int, default=2000, help="angle parameter for the verifier")
    ap.add_argument("--threads", type=int, default=0, help="0 = os.cpu_count()")
    ap.add_argument("--verifier", default=str(DEFAULT_VERIFIER))
    ap.add_argument("--lo", type=int, default=None,
                    help="smallest D to consider (largest container). default: 0.97 * D")
    ap.add_argument("--hi", type=int, default=None,
                    help="largest D to consider (smallest container). default: cert's own D")
    args = ap.parse_args()

    import os
    threads = args.threads or os.cpu_count() or 4

    verifier = Path(args.verifier)
    if not verifier.exists():
        sys.exit(f"verifier not built: {verifier}\nrun:  (cd verify && cargo build --release)")

    s0, D0, W, rows = read_cert(args.cert)
    K = s0 * D0
    if K.denominator != 1:
        sys.exit(f"s*D = {K} is not an integer; this certificate is not on a scalable grid")
    K = K.numerator

    print(f"certificate : {args.cert}")
    print(f"  {len(rows)} points, s = {s0} = {float(s0):.9f}, D = {D0}")
    print(f"  invariant K = s*D = {K} (container side in grid units)")
    print(f"  total weight = {sum(r[2] for r in rows)}/{W} = {sum(r[2] for r in rows)/W:.7f}")
    print()

    # hi = smallest container we are confident about (the certificate's own D).
    # lo = an ambitious D (bigger container) that we expect to fail.
    hi = args.hi if args.hi is not None else D0
    lo = args.lo if args.lo is not None else int(D0 * 0.97)
    if lo >= hi:
        sys.exit("need lo < hi")

    def side(D):
        return Fraction(K, D)

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "scaled.txt"

        def ok(D):
            write_cert(tmp, side(D), D, W, rows)
            r = verifies(verifier, tmp, args.n, args.N, threads)
            print(f"  D = {D:>9}   s = {float(side(D)):.9f}   "
                  f"{'VERIFIED' if r else 'rejected'}")
            return r

        print(f"bracketing in D ∈ [{lo}, {hi}]  (smaller D = larger container = stronger)")
        if not ok(hi):
            sys.exit("the certificate does not even verify at its own D -- nothing to scale")
        if ok(lo):
            print("\nlo end still verifies; widen the bracket with a smaller --lo")
            best = lo
        else:
            # invariant: ok(hi) true, ok(lo) false. Find smallest verifying D.
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if ok(mid):
                    hi = mid
                else:
                    lo = mid
            best = hi

        s_best = side(best)
        print()
        print(f"critical D = {best}")
        print(f"best container side = {s_best} = {float(s_best):.9f}")

        # Reduce: if every coordinate shares a factor with D, write it smaller.
        g = math.gcd(best, *(x for r in rows for x in r[:2]))
        if g > 1:
            print(f"all coordinates and D share a factor {g}; reduced form is "
                  f"D = {best // g}, s = {s_best}")

        print()
        print(f"==> s({args.n}) >= {s_best} = {float(s_best):.9f}")
        if s_best > s0:
            gain = float(s_best) - float(s0)
            print(f"    (an improvement of {gain:.9f} over the input certificate)")


if __name__ == "__main__":
    main()
