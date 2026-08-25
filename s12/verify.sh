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
echo "=== independent exact re-check (Python Fractions, sampled bins) ==="
python3 xcheck.py certificates/s12_56points_3.8.txt 2000 5
echo
echo "=== rejection tests (a verifier that never says no is worthless) ==="
./tests/rejection_tests.sh
echo
echo "=== points.json companions: json -> txt reproduces the shipped .txt byte for byte ==="
python3 search/export_points.py --roundtrip certificates/s12_lower_3.931795.txt certificates/s12_lower_3.931795.json
python3 search/export_points.py --roundtrip certificates/s12_56points_3.8.txt certificates/s12_56points_3.8.json
