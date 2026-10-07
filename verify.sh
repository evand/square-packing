#!/bin/sh
# Re-verify the certificates in this repo, from scratch.
#
#   ./verify.sh --spot   spot tier (CI on every push): build both checkers, hashes, the s(12)/s(11)
#                        point certificates, the shipped run records of every bundle re-summarised;
#                        no fresh sweeps, no rejection suites.  ~1 min here, a few on a CI runner.
#   ./verify.sh          fast tier (run locally after touching a checker or a certificate; CI only
#                        by hand): adds the bundles' fresh zmx2 sweeps, the unavoid13 checks and all
#                        rejection suites.  ~4 min on 16 cores, ~18 min on a 4-core runner.
#   ./verify.sh --full   adds the slow sweeps: the N = 12000 re-runs, the box-clique Python
#                        re-check, and the full-domain zero-margin sweep for s(13) = 4
#                        (~35 min at 4 threads here; 3 h 40 min on a GitHub runner).
set -e
cd "$(dirname "$0")"
FULL=0; SPOT=
case "${1:-}" in --full) FULL=1 ;; --spot) SPOT=1 ;; --fast|"") ;; *) echo "usage: $0 [--spot|--fast|--full]"; exit 2 ;; esac
export SPOT   # the bundle scripts skip their fresh sweeps when it is set
full() { [ "$FULL" = 1 ]; }
spot() { [ -n "$SPOT" ]; }
echo "=== building verifiers ==="
( cd verify && cargo build --release )
( cd verify2 && cargo build --release )
V=verify/target/release/verify
Z=verify2/target/release/zmcheck
# Exit status of `verify` and `zmcheck` (since 2026-10-07; jlevy/squares#238): 0 VERIFIED,
# 1 NOT VERIFIED, 2 ERROR (bad input), 3 internal error, 4 partial run (no verdict), 101 panic.
# chk still requires the VERIFIED: line too, so a verdict needs both the line and exit 0.
chk() { rc=0; out=$("$@") || rc=$?
        [ $rc -le 1 ] || [ $rc -eq 4 ] || { echo "$out"; echo "verifier failed (exit $rc): $*"; exit 1; }
        echo "$out" | grep -E '^(==>|    i\.e\.|min |VERIFIED|NOT VERIFIED)'
        [ $rc -eq 0 ] && echo "$out" | grep -q '^VERIFIED:' || { echo "REJECTED (exit $rc): $*"; exit 1; }; }
# Same for `zmcheck`: `NOT VERIFIED` (exit 1) and `PARTIAL SWEEP` (exit 4: a restricted sweep
# is a refusal to give a verdict) both fail the script here.
chkzm() { rc=0; out=$("$@") || rc=$?
          [ $rc -le 1 ] || [ $rc -eq 4 ] || { echo "$out"; echo "zmcheck failed (exit $rc): $*"; exit 1; }
          echo "$out" | grep -E '^(total weight|done in |  leaves:|VERIFIED|NOT VERIFIED|PARTIAL SWEEP)'
          [ $rc -eq 0 ] && echo "$out" | grep -q '^VERIFIED:' || { echo "REJECTED (exit $rc): $*"; exit 1; }; }
echo
echo "=== main certificate: s(12) >= 15680/3951 = 3.968616 ==="
chk $V certificates/s12_lower_3.9686.txt 12 6000 "$(nproc)" 0
echo
if full; then echo "=== same certificate, finer angle grid (N=12000) ==="
chk $V certificates/s12_lower_3.9686.txt 12 12000 "$(nproc)" 0; fi
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
full && chk $V certificates/s11_lower_3.8143.txt 11 12000 "$(nproc)" 0
echo
echo "=== box-clique demonstration certificate (points + 8 clique orbits; FORMAT.md 'Clique certificates'; net N=2000 is part of the file) ==="
chk $V certificates/s12_boxclique_demo_3.9318_N2000.txt 12 2000 "$(nproc)" 0
full && python3 xcheck.py certificates/s12_boxclique_demo_3.9318_N2000.txt 2000 --all --n 12
echo
echo "=== anchor-clique demonstration certificate (223 points + 3 zero-weight anchor atoms + one"
echo "    K(p,A); the clique carries 0.1198 of the total and IS load-bearing; FORMAT.md 'Anchor"
echo "    cliques'.  Unlike a box clique it does not refer to the angle net, so any N works) ==="
chk $V certificates/s12_anchorclique_demo_3.9318.txt 12 6000 "$(nproc)" 0
full && chk $V certificates/s12_anchorclique_demo_3.9318.txt 12 12000 "$(nproc)" 0
echo
echo "=== independent exact re-check (Python, exact rationals, EVERY angle bin) ==="
# Exhaustive: all 829 bins of the N=2000 net (the same net the Rust run above used, so the
# per-bin minima are directly comparable).  56 points -> well under a minute even on 2 cores;
# the weighted certificates at N=6000 take 1-5 min on 32 cores, too slow for CI, so they are
# left to `python3 xcheck.py certificates/s12_lower_3.9686.txt 6000 --all --n 12` by hand.
python3 xcheck.py certificates/s12_56points_3.8.txt 2000 --all --n 12
echo
echo "=== s(13) = 4, case-free: the weighted closed cover of [0,4]^2, total weight"
echo "    2591194431/200000000 = 12.955972155 < 13, checked at margin zero over the FULL"
echo "    pose domain (cx, cy in [0,4], u = tan(th/2) in [0,1]) by exhaustive subdivision."
echo "    README 's(13) = 4 without case analysis'; search/RUNG2_XCHECK.md. ~17 min on 8 threads. ==="
if full; then chkzm $Z cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --threads "$(nproc)"
else echo "    (full-domain sweep skipped in the fast tier; the rung-2 rejection tests below still sweep"
     echo "    the certificate's tightest band.  Run ./verify.sh --full for the whole domain.)"; fi
# Slow path, not run here: the independent Python checker (exact `fractions.Fraction`, its own
# subdivision and its own primitive set), over its symmetry-reduced domain.  1 h 42 min on 8
# processes -- too slow for this script and for CI, so it is left to be run by hand, the way the
# N=6000 `xcheck.py` runs above are:
#   python3 search/zeromargin.py cert certificates/rung2/s13_closed_cover_4.txt \
#           --depth 18 --nproc 8 --disj --chain-from 0 --dump runs/leaves.txt
#   -> done in 6138s: boxes 16872, max depth 13
#      leaves: ADM 2867  CORE 0  P1 0  MIX 0  CHAIN 5320  TRI 0  EMPTY 3449  UNCERTIFIED 0
#      VERIFIED
echo
echo "=== s(32) = 6: weighted closed cover of [0,6]^2, total 3171350535386/10^11 = 31.713505 < 32,"
echo "    certified over the D4 fundamental region by zeromargin.py (certificates/s32/README.md) ==="
# The full re-sweep is ~2.8 CPU-h per cover; CI skips it (S32_SEPARATE=1): run it locally.
if full && [ -z "${S32_SEPARATE:-}" ]; then certificates/s32/verify.sh --full; else certificates/s32/verify.sh; fi
echo
echo "=== s(21) = 5: mixed cover of [0,5]^2 (7,536 points + 1,872 grid-line segments), total"
echo "    522368729933/25000000000 = 20.894749 < 21 (certificates/s21/README.md).  Fast tier: shipped records of"
echo "    both checkers re-summarised, and a fresh zmx2 run over the whole pose space with no symmetry assumed ==="
# The zm_mixed.py re-sweep is ~20 CPU-h; CI skips it (S21_SEPARATE=1): run it locally.
if full && [ -z "${S21_SEPARATE:-}" ]; then certificates/s21/verify.sh --full; else certificates/s21/verify.sh; fi
echo
echo "=== s(45) = 7: mixed cover of [0,7]^2 (19,989 points + 3,912 grid-line segments), total"
echo "    2238676387/50000000 = 44.773528 < 45 (certificates/s45/README.md; no Lean).  Fast tier: shipped records of"
echo "    both checkers re-summarised, and a fresh zmx2 run over the whole pose space with no symmetry assumed ==="
# The zm_mixed.py re-sweep is ~16.5 CPU-h; CI skips it (S45_SEPARATE=1): run it locally.
if full && [ -z "${S45_SEPARATE:-}" ]; then certificates/s45/verify.sh --full; else certificates/s45/verify.sh; fi
echo
echo "=== s(60) = 8: mixed cover of [0,8]^2 (23,744 points + 5,216 grid-line segments), total"
echo "    748233441/12500000 = 59.858675 < 60 (certificates/s60/README.md; no Lean).  Fast tier: shipped records of"
echo "    both checkers re-summarised, and a fresh zmx2 run over the whole pose space with no symmetry assumed ==="
# The zm_mixed.py re-sweep is ~19.6 CPU-h; CI skips it (S60_SEPARATE=1): run it locally.
if full && [ -z "${S60_SEPARATE:-}" ]; then certificates/s60/verify.sh --full; else certificates/s60/verify.sh; fi
echo
if ! spot; then
echo "=== rung-2 rejection tests (23 checks: mutations, invalid historical covers, malformed input) ==="
./tests/rung2/rejection_tests.sh
echo
echo "=== unavoidable point sets (notes/unavoid13.md, certificates/unavoid13/README.md): p(3) = 7, and no"
echo "    half-turn-symmetric 13-point set for [0,4]^2.  Lower bounds: exact branch-and-bound certificates"
echo "    checked by search/unavoid13_exactcheck.py (Fractions; own vertices, incidence, domination, tree) ==="
echo "    p(3) <= 7: the rational Kearney-Shiu set ks7_rational_3.txt is unavoidable (zero-margin checker + SEG)"
out=$(python3 search/unavoid13_check.py cert certificates/unavoid13/ks7_rational_3.txt --tri --seg --depth 16) \
  || { echo "$out"; echo "checker failed: unavoid13_check.py"; exit 1; }
echo "$out" | grep -E '^(done in |  leaves:|VERIFIED|NOT VERIFIED)'
echo "$out" | grep -qx 'VERIFIED' || { echo "REJECTED: ks7_rational_3.txt"; exit 1; }
echo "    p(3) >= 7: no 6 points meet the 87 squares of unavoid3_lower7_family.txt"
chk python3 search/unavoid13_exactcheck.py certificates/unavoid13/unavoid3_lower7_family.txt \
    certificates/unavoid13/unavoid3_lower7_bb.txt --threads 4
echo "    no half-turn-symmetric set of <= 13 points meets the 441 squares of unavoid4_C2_family.txt (~25 s)"
chk python3 search/unavoid13_exactcheck.py certificates/unavoid13/unavoid4_C2_family.txt \
    certificates/unavoid13/unavoid4_C2_bb.txt --threads 4
echo
echo "=== unavoid13 rejection tests (22 checks: mutated families, tampered duals, dropped subtrees) ==="
./tests/unavoid13/rejection_tests.sh
echo
echo "=== rejection tests (a verifier that never says no is worthless) ==="
./tests/rejection_tests.sh
echo
fi   # ! spot
echo "=== points.json companions: json -> txt reproduces the shipped .txt byte for byte ==="
for c in certificates/s12_lower_*.txt; do python3 search/export_points.py --roundtrip "$c" "${c%.txt}.json"; done
python3 search/export_points.py --roundtrip certificates/s12_56points_3.8.txt certificates/s12_56points_3.8.json
python3 search/export_points.py --roundtrip certificates/s11_lower_3.8143.txt certificates/s11_lower_3.8143.json
for c in certificates/s12_uniform_*.txt; do python3 search/export_points.py --roundtrip "$c" "${c%.txt}.json"; done
echo
if full; then echo "=== ALL CHECKS PASSED (full tier) ==="; elif spot; then echo "=== ALL CHECKS PASSED (spot tier; ./verify.sh adds the fresh sweeps and rejection suites) ==="; else echo "=== ALL CHECKS PASSED (fast tier; ./verify.sh --full adds the slow sweeps) ==="; fi
