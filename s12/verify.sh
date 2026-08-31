#!/bin/sh
# Re-verify every certificate in this repo, from scratch.
set -e
cd "$(dirname "$0")"
echo "=== building verifier ==="
( cd verify && cargo build --release )
V=verify/target/release/verify
# `verify` exits 0 on NOT VERIFIED (it is a verdict, not an error), so under `set -e` the
# verdict must be checked explicitly or CI would stay green on a rejected certificate.
chk() { out=$("$@") || { echo "$out"; echo "verifier failed: $*"; exit 1; }
        echo "$out" | grep -E '^(==>|    i\.e\.|min |VERIFIED|NOT VERIFIED)'
        echo "$out" | grep -q '^VERIFIED:' || { echo "REJECTED: $*"; exit 1; }; }
echo
echo "=== main certificate: s(12) >= 15680/3951 = 3.968616 ==="
chk $V certificates/s12_lower_3.9686.txt 12 6000 "$(nproc)" 0
echo
echo "=== same certificate, finer angle grid (N=12000) ==="
chk $V certificates/s12_lower_3.9686.txt 12 12000 "$(nproc)" 0
echo
echo "=== 764-point certificate: s(12) >= 980/247 = 3.967611 ==="
chk $V certificates/s12_lower_3.9676.txt 12 6000 "$(nproc)" 0
echo
echo "=== the earlier bound 3920/997: 224-point and original 788-point certificates ==="
chk $V certificates/s12_lower_3.931795_sparse.txt 12 6000 "$(nproc)" 0
chk $V certificates/s12_lower_3.931795.txt 12 6000 "$(nproc)" 0
echo
echo "=== 56-point certificate: s(12) >= 19/5 ==="
chk $V certificates/s12_56points_3.8.txt 12 2000 "$(nproc)" 0
echo
echo "=== 81-point uniform certificate: s(12) >= 35/9 (every unit square contains 7 of 81 points) ==="
chk $V certificates/s12_uniform_7of81_3.888.txt 12 2000 "$(nproc)" 0
python3 xcheck.py certificates/s12_uniform_7of81_3.888.txt 2000 --all --n 12
echo
echo "=== the other uniform certificates (k of m points, see search/uniform/UNIFORM.md) ==="
for c in certificates/s12_uniform_*.txt; do echo "$c"; chk $V "$c" 12 2000 "$(nproc)" 0 | grep -E "^(VERIFIED|NOT)"; done
echo
echo "=== n = 11: s(11) >= 3040/797 = 3.814304 (680 points, total weight 10.8146708 < 11) ==="
chk $V certificates/s11_lower_3.8143.txt 11 6000 "$(nproc)" 0
chk $V certificates/s11_lower_3.8143.txt 11 12000 "$(nproc)" 0
echo
echo "=== box-clique demonstration certificate (points + 8 clique orbits; FORMAT.md 'Clique certificates'; net N=2000 is part of the file) ==="
chk $V certificates/s12_boxclique_demo_3.9318_N2000.txt 12 2000 "$(nproc)" 0
python3 xcheck.py certificates/s12_boxclique_demo_3.9318_N2000.txt 2000 --all --n 12
echo
echo "=== anchor-clique demonstration certificate (223 points + 3 zero-weight anchor atoms + one"
echo "    K(p,A); the clique carries 0.1198 of the total and IS load-bearing; FORMAT.md 'Anchor"
echo "    cliques'.  Unlike a box clique it does not refer to the angle net, so any N works) ==="
chk $V certificates/s12_anchorclique_demo_3.9318.txt 12 6000 "$(nproc)" 0
chk $V certificates/s12_anchorclique_demo_3.9318.txt 12 12000 "$(nproc)" 0
echo
echo "=== independent exact re-check (Python, exact rationals, EVERY angle bin) ==="
# Exhaustive: all 829 bins of the N=2000 net (the same net the Rust run above used, so the
# per-bin minima are directly comparable).  56 points -> well under a minute even on 2 cores;
# the weighted certificates at N=6000 take 1-5 min on 32 cores, too slow for CI, so they are
# left to `python3 xcheck.py certificates/s12_lower_3.9686.txt 6000 --all --n 12` by hand.
python3 xcheck.py certificates/s12_56points_3.8.txt 2000 --all --n 12
echo
echo "=== rejection tests (a verifier that never says no is worthless) ==="
./tests/rejection_tests.sh
echo
echo "=== points.json companions: json -> txt reproduces the shipped .txt byte for byte ==="
for c in certificates/s12_lower_*.txt; do python3 search/export_points.py --roundtrip "$c" "${c%.txt}.json"; done
python3 search/export_points.py --roundtrip certificates/s12_56points_3.8.txt certificates/s12_56points_3.8.json
python3 search/export_points.py --roundtrip certificates/s11_lower_3.8143.txt certificates/s11_lower_3.8143.json
for c in certificates/s12_uniform_*.txt; do python3 search/export_points.py --roundtrip "$c" "${c%.txt}.json"; done
