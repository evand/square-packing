#!/bin/sh
# Re-verify every certificate in this repo, from scratch.
set -e
cd "$(dirname "$0")"
echo "=== building verifier ==="
( cd verify && cargo build --release )
V=verify/target/release/verify
echo
echo "=== main certificate: s(12) >= 3920/997 ==="
$V certificates/s12_lower_3.931795.txt 12 6000 "$(nproc)" 0
echo
echo "=== same certificate, finer angle grid (N=12000) ==="
$V certificates/s12_lower_3.931795.txt 12 12000 "$(nproc)" 0 | tail -3
echo
echo "=== 56-point certificate: s(12) >= 19/5 ==="
$V certificates/s12_56points_3.8.txt 12 2000 "$(nproc)" 0
echo
echo "=== independent exact re-check (Python, exact rationals, EVERY angle bin) ==="
# Exhaustive: all 829 bins of the N=2000 net (the same net the Rust run above used, so the
# per-bin minima are directly comparable).  56 points -> well under a minute even on 2 cores;
# the 788-point certificate at N=6000 takes ~1 min on 32 cores, too slow for CI, so it is
# left to `python3 xcheck.py certificates/s12_lower_3.931795.txt 6000 --n 12` by hand.
python3 xcheck.py certificates/s12_56points_3.8.txt 2000 --all --n 12
echo
echo "=== rejection tests (a verifier that never says no is worthless) ==="
./tests/rejection_tests.sh
