#!/usr/bin/env python3
"""Scale a uniform certificate to its critical container and write the scaled file.

Same arithmetic as search/scale_to_critical.py (integer coordinates fixed, D shrinks, s = K/D
grows), but (1) the coordinates are first multiplied by M so that D >= DMIN (finer steps in s),
(2) the scaled certificate is WRITTEN to OUT (reduced by the gcd of D and all coordinates), and
(3) the result is re-checked with the verifier at N=2000 and N=8000.

    python3 search/uniform/critical.py CERT OUT [--dmin 200000] [--N 2000] [--threads 32] [--n 12]
"""
import argparse, math, os, subprocess, sys, tempfile
from fractions import Fraction
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
VERIFY = REPO / "verify" / "target" / "release" / "verify"


def read_cert(path):
    t = Path(path).read_text().split()
    sn, sd, D, W, m = map(int, t[:5])
    rows = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    return Fraction(sn, sd), D, W, rows


def write_cert(path, s, D, W, rows):
    with open(path, "w") as f:
        f.write(f"{s.numerator} {s.denominator}\n{D}\n{W}\n{len(rows)}\n")
        for x, y, w in rows: f.write(f"{x} {y} {w}\n")


def verifies(cert, n, N, threads):
    r = subprocess.run([str(VERIFY), str(cert), str(n), str(N), str(threads), "0"], capture_output=True, text=True)
    return r.returncode == 0 and "VERIFIED" in r.stdout and "NOT VERIFIED" not in r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cert"); ap.add_argument("out")
    ap.add_argument("--dmin", type=int, default=200000); ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 8); ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--lo", type=float, default=0.93, help="lower end of the D bracket as a fraction of D")
    a = ap.parse_args()
    s0, D0, W, rows = read_cert(a.cert)
    M = max(1, -(-a.dmin // D0))
    D0 *= M; rows = [(x * M, y * M, w) for x, y, w in rows]
    K = s0 * D0; assert K.denominator == 1; K = K.numerator
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "c.txt"
        def ok(D):
            write_cert(tmp, Fraction(K, D), D, W, rows)
            return verifies(tmp, a.n, a.N, a.threads)
        hi = D0; lo = int(D0 * a.lo)
        if not ok(hi): sys.exit("does not verify at its own D")
        if ok(lo): sys.exit("lo end verifies; widen --lo")
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ok(mid): hi = mid
            else: lo = mid
        best = hi
    g = math.gcd(best, *(v for r in rows for v in r[:2]))
    D = best // g; rows = [(x // g, y // g, w) for x, y, w in rows]
    s = Fraction(K, best)
    write_cert(a.out, s, D, W, rows)
    v2 = verifies(a.out, a.n, 2000, a.threads); v8 = verifies(a.out, a.n, 8000, a.threads)
    print(f"CRITICAL {a.cert}: m={len(rows)} k={W} s0={float(s0):.6f} -> s={s} = {float(s):.9f} D={D}  verify N=2000:{v2} N=8000:{v8}  -> {a.out}")


if __name__ == "__main__":
    main()
