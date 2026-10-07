#!/bin/sh
# Rejection tests for search/unavoid13_exactcheck.py (the exact hitting-set certificates in
# certificates/unavoid13/).  Every input below is a mutation of a shipped family or certificate,
# built in a temporary directory; each must produce the stated verdict.
#   VERIFIED  exit 0 and a "VERIFIED:" line
#   REJECT    exit 0 and "NOT VERIFIED"
#   ERROR     non-zero exit, an "ERROR:" line, and no verdict word
# The C2 mutations use --float-filter (same verdict as the all-Fraction default, ~10x faster).
cd "$(dirname "$0")/../.."
X="python3 search/unavoid13_exactcheck.py"
D=certificates/unavoid13
F3=$D/unavoid3_lower7_family.txt;  C3=$D/unavoid3_lower7_bb.txt
F4=$D/unavoid4_C2_family.txt;      C4=$D/unavoid4_C2_bb.txt
T=$(mktemp -d /tmp/unavoid13rej.XXXXXX) || exit 1
trap 'rm -rf "$T"' EXIT
LAST=$T/last.out
pass=0; fail=0
check() {  # name, expect, family, cert, [extra args]
  name=$1; want=$2; shift 2
  $X "$@" --threads 4 >"$LAST" 2>&1; rc=$?
  if grep -q '^NOT VERIFIED' "$LAST" && [ $rc -eq 0 ]; then got=REJECT
  elif grep -q '^VERIFIED:' "$LAST" && [ $rc -eq 0 ]; then got=VERIFIED
  elif [ $rc -ne 0 ] && grep -q '^ERROR:' "$LAST" && ! grep -q 'VERIFIED' "$LAST"; then got=ERROR
  else got="UNKNOWN(rc=$rc)"; fi
  if [ "$got" = "$want" ]; then pass=$((pass+1)); printf '  ok    %-60s -> %s\n' "$name" "$got"
    grep -E '^(NOT VERIFIED|ERROR)' "$LAST" | head -1 | cut -c1-110 | sed 's/^/          /'
  else fail=$((fail+1)); printf '  FAIL  %-60s -> %s (wanted %s)\n' "$name" "$got" "$want"
    tail -3 "$LAST" | sed 's/^/        | /'; fi
}
# line numbers of the first pose / first tree records, for the mutations
first_pose=$(grep -n '^pose' $F3 | head -1 | cut -d: -f1)
treeline=$(grep -n '^tree$' $C3 | cut -d: -f1)
endline=$(grep -n '^end$' $C3 | cut -d: -f1)

echo "-- baseline (T = 3)"
check "unmodified T = 3 certificate"                      VERIFIED $F3 $C3

echo "-- mutated families"
# square 0 (the corner tile [0,1]^2) replaced by a copy of square 1: 6 points now suffice
awk -v L=$first_pose 'NR==L{getline nxt; print nxt; print nxt; next}{print}' $F3 >$T/f1.txt
check "tile [0,1]^2 replaced by a copy of another square"  REJECT $T/f1.txt $C3
awk -v L=$first_pose 'NR!=L' $F3 >$T/f2.txt
check "tile [0,1]^2 deleted (rows shift)"                 REJECT $T/f2.txt $C3
sed "${first_pose}s|.*|pose 0 1/3 1/2|" $F3 >$T/f3.txt
check "a square pushed outside the container"             REJECT $T/f3.txt $C3
sed '1s|m = 3|m = 4|' $F3 >$T/f4.txt
check "family relabelled as m = 4"                        REJECT $T/f4.txt $C3

echo "-- tampered duals"
awk -v a=$treeline '/^D /&&NR>a&&!done{$2=$2*2; done=1}{print}' $C3 >$T/c1.txt
check "first leaf's duals halved (den doubled)"           REJECT $F3 $T/c1.txt
awk -v a=$treeline '/^D /&&NR>a&&!done{split($3,p,":"); $3=p[1]":-"p[2]; done=1}{print}' $C3 >$T/c2.txt
check "a negative dual"                                   REJECT $F3 $T/c2.txt
awk -v a=$treeline '/^D /&&NR>a&&!done{$0="D 1"; done=1}{print}' $C3 >$T/c3.txt
check "first leaf's duals zeroed"                         REJECT $F3 $T/c3.txt

echo "-- tree structure"
awk -v e=$endline 'NR!=e-1' $C3 >$T/t1.txt
check "last leaf dropped (x_j = 0 subtree missing)"       REJECT $F3 $T/t1.txt
awk -v a=$treeline 'NR>a && /^D / && !done {done=1; next}{print}' $C3 >$T/t2.txt
check "first leaf dropped (subtrees misaligned)"          REJECT $F3 $T/t2.txt
awk -v e=$endline 'NR==e{print "D 1 0:7"}{print}' $C3 >$T/t3.txt
check "an extra leaf after the root is closed"            REJECT $F3 $T/t3.txt
j1=$(awk -v a=$treeline 'NR>a && /^B /{print $2; exit}' $C3)
awk -v a=$treeline -v j=$j1 'NR>a && /^B /{n++; if(n==2){$2=j}}{print}' $C3 >$T/t4.txt
check "second branching on an already-fixed candidate"    REJECT $F3 $T/t4.txt
awk -v a=$treeline -v e=$endline 'NR<=a{print} NR==e{print}' $C3 | sed 's/^tree$/tree\nD 1 0:1/' >$T/t5.txt
check "whole tree replaced by one weak leaf"              REJECT $F3 $T/t5.txt

echo "-- candidates and claim"
lastcand=$(grep -n '^cand' $C3 | tail -1 | cut -d: -f1)
awk -v L=$lastcand 'NR!=L' $C3 >$T/k1.txt
check "last candidate deleted (a vertex is undominated)"  REJECT $F3 $T/k1.txt
sed 's/^k 6$/k 7/' $C3 >$T/k2.txt
check "claim strengthened to k = 7 (7 points do exist)"   REJECT $F3 $T/k2.txt
sed 's/^k 6$/k 5/' $C3 >$T/k3.txt
check "claim weakened to k = 5 (still true)"              VERIFIED $F3 $T/k3.txt

echo "-- C2 certificate"
sed 's/^k 13$/k 14/' $C4 >$T/s1.txt
check "C2 claim at k = 14 (Friedman's 14 is half-turn symmetric)" REJECT $F4 $T/s1.txt --float-filter
sed 's/^claim C2$/claim plain/' $C4 >$T/s2.txt
check "C2 certificate read as an asymmetric claim"        REJECT $F4 $T/s2.txt --float-filter
c1=$(grep -n '^cand' $C4 | head -1 | cut -d: -f1)
sed "${c1}s|.*|cand 2 2001/1000|" $C4 >$T/s3.txt
check "centre candidate moved off the centre"             REJECT $F4 $T/s3.txt --float-filter

echo "-- malformed input"
head -n $((endline-3)) $C3 >$T/e1.txt
check "truncated certificate (no end record)"             ERROR $F3 $T/e1.txt
sed 's/^claim plain$/claim D4/' $C3 >$T/e2.txt
check "unknown claim"                                     ERROR $F3 $T/e2.txt
check "missing family file"                               ERROR $T/nope.txt $C3

printf '\n%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
