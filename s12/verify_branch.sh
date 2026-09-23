#!/bin/sh
# Re-verify the branch (corner-occupancy) certificates of search/BRANCH.md.  Not part of
# verify.sh / CI: the three files have 10-12k points each and the arrangement sweep is quadratic
# in the point count, so this takes about an hour on 32 cores (N = 2000; N = 4000 doubles it).
# They prove nothing about s(12) by themselves -- the all-corners leaf k = 4 is open at 3.98 --
# they are shipped as worked examples of the format: each one refutes every packing of 12 unit
# squares with exactly k squares centred in the corner boxes [0,1]^2 of a container of side < 3.98.
set -e
cd "$(dirname "$0")"
( cd verify && cargo build --release )
V=verify/target/release/verify
N=${N:-2000}
chk() { out=$("$@") || { echo "$out"; echo "verifier failed: $*"; exit 1; }
        echo "$out" | grep -E '^(BRANCH|==>|min |VERIFIED|NOT VERIFIED)' | cut -c1-200
        echo "$out" | grep -q '^VERIFIED:' || { echo "REJECTED: $*"; exit 1; }; }
for k in 0 1 2; do
  echo "=== corner leaf k=$k at s = 199/50 = 3.98 (N=$N) ==="
  chk $V certificates/branch/s12_t3.98_corner_k$k.txt 12 $N "$(nproc)" 0
done
echo "all branch certificates verified (leaves k=0,1,2 at 3.98; k=3 needs per-box multipliers and k=4 is open, see search/BRANCH.md)"
