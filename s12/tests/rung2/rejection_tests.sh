#!/bin/sh
# Rejection tests for verify2/ (`zmcheck`), the independent zero-margin checker.
#
# A checker that always says VERIFIED is worthless.  Every file below is a mutation of the shipped
# rung-2 certificate `certificates/rung2/s13_closed_cover_4.txt`, or an invalid cover from the
# development history; each must produce the stated verdict.
#
# Verdicts:  VERIFIED  exit 0 and a "VERIFIED:" line
#            REJECT    exit 0 and "NOT VERIFIED"   (a clean, reasoned refusal to certify)
#            VIOLATION exit 1 and "*** VIOLATION"  (`pose` mode: an exact captured weight < 1,
#                                                   i.e. a *disproof*, not just a refusal)
#            ERROR     non-zero exit, an "ERROR:" message, and NO verdict word in the output
#                      (malformed / illegal input: the checker refuses to give a verdict)
#            PANIC     a Rust panic or abnormal exit.  Never acceptable, even on garbage input:
#                      a panic is not a rejection.  Any PANIC fails the whole script.
#
# The full-domain sweep of the shipped certificate takes ~20 min on 8 threads, so the mutation
# tests here are run over a *restricted* centre band, `cx in [0.5,0.6]`, `cy in [1.0,2.0]`
# (`--xlo/--xhi/--ylo/--yhi`), which contains the pose `(1/2, 3/2, 0)` that `RUNG2.md` sec 2
# identifies as the worst monotone-witness pose of the container -- the place where the
# certificate has the least room.  A restricted sweep never prints VERIFIED (it prints
# `PARTIAL SWEEP` and then `NOT VERIFIED`), so for those runs the test compares the
# **uncertified-box count**: `0` for the baseline, `> 0` for a mutation that the checker refuses.
# The whole-file verdicts (VERIFIED / REJECT / ERROR / VIOLATION) are exercised by the `pose`
# mode and the malformed inputs, which are instant.
#
#   sh tests/rung2/rejection_tests.sh          # ~4 min
#   ZM=path/to/zmcheck sh tests/rung2/rejection_tests.sh
cd "$(dirname "$0")/../.." || exit 1
Z=${ZM:-verify2/target/release/zmcheck}
[ -x "$Z" ] || { echo "build it first:  (cd verify2 && cargo build --release)"; exit 1; }
C=certificates/rung2/s13_closed_cover_4.txt
X103=${X103:-tests/rung2/inputs/closed4_best_x103.txt}
T=$(mktemp -d /tmp/rung2rej.XXXXXX) || exit 1
trap 'rm -rf "$T"' EXIT
LAST=$T/last.out
BAND="--xlo 0.5 --xhi 0.6 --ylo 1.0 --yhi 2.0"
pass=0; fail=0; panics=0
ok()  { pass=$((pass+1)); printf '  ok    %-54s -> %s\n' "$1" "$2"; }
bad() { fail=$((fail+1)); printf '  FAIL  %-54s -> %s (wanted %s)\n' "$1" "$2" "$3"; }

verdict() {  # sets $got from $LAST and $rc
  if grep -q 'panicked at' "$LAST" || [ "$rc" -ge 100 ]; then got=PANIC; panics=$((panics+1))
  elif grep -q '\*\*\* VIOLATION' "$LAST" && [ "$rc" -eq 1 ]; then got=VIOLATION
  elif grep -q '^NOT VERIFIED' "$LAST"; then got=REJECT
  elif grep -q '^VERIFIED:' "$LAST" && [ "$rc" -eq 0 ]; then got=VERIFIED
  elif [ "$rc" -ne 0 ] && ! grep -q 'VERIFIED' "$LAST"; then got=ERROR
  else got="UNKNOWN(rc=$rc)"; fi
}
check() {  # name, expect, args...
  name=$1; want=$2; shift 2
  $Z "$@" >"$LAST" 2>&1; rc=$?
  verdict
  if [ "$got" = "$want" ]; then ok "$name" "$got"; else bad "$name" "$got" "$want"
    [ "$got" = PANIC ] && sed 's/^/        | /' "$LAST" | head -3; fi
}
band() {  # name, want (zero|nonzero), file
  name=$1; want=$2; f=$3
  # shellcheck disable=SC2086
  $Z cert "$f" --threads 4 --depth 18 $BAND >"$LAST" 2>&1; rc=$?
  if grep -q 'panicked at' "$LAST"; then panics=$((panics+1)); bad "$name" PANIC "$want"; return; fi
  n=$(sed -n 's/^NOT VERIFIED: \([0-9]*\) uncertified.*/\1/p' "$LAST")
  [ -n "$n" ] || { bad "$name" "no verdict line" "$want"; return; }
  cen=$(sed -n 's/^  leaves: *//p' "$LAST" | tr -s ' ')
  case "$want:$n" in
    zero:0)    ok "$name" "0 uncertified  [$cen]" ;;
    nonzero:0) bad "$name" "0 uncertified  [$cen]" "> 0 uncertified" ;;
    zero:*)    bad "$name" "$n uncertified" "0 uncertified" ;;
    *)         ok "$name" "$n uncertified  [$cen]" ;;
  esac
}

echo "-- exact single-pose cross-checks against search/RUNG2.md (no boxes, no subdivision)"
$Z pose "$C" --x 7/2 --y 7/2 --u 0 >"$LAST" 2>&1; rc=$?; verdict
[ "$got" = VERIFIED ] || [ "$rc" -eq 0 ] && ok "shipped cover at (7/2,7/2,0) captures >= 1" \
  "$(sed -n 's/.*EXACT captured weight = //p' "$LAST")" || bad "shipped cover at (7/2,7/2,0)" "$got" OK

echo "-- baseline: the shipped certificate over the hard band must be fully certified"
band "unmodified certificate, cx in [0.5,0.6] cy in [1,2]" zero "$C"

echo "-- the five mutations of the brief"
# (a) one point's weight reduced by 1 %.  Line 273 is the point (1.000, 1.560), weight
#     8722560/10^9: it lies on the heavy grid line x = 1 next to the worst pose (1/2, 3/2, 0).
#     NOTE: this file is STILL A VALID COVER -- the cover's minimum captured weight is 1.03138
#     (RUNG2.md sec 4.7) and the change is 8.7e-5 -- so a checker that "refused" it would be
#     refusing a true statement.  The certification survives too, which is the honest outcome.
awk 'NR==273{print $1, $2, int($3*0.99); next}{print}' "$C" >"$T/a.txt"
band "(a) one point's weight reduced by 1 % [still valid]" zero "$T/a.txt"
# (a') the same mutation with teeth: every weight reduced by 5 % (total 12.308 < 13).  Now the
#     cover is false, and at the corner pose (1/2,1/2,0) -- where it captures 1.0500063, its
#     tightest axis-parallel pose -- the exact capture drops below 1 and is exhibited.
awk 'NR<=4{print} NR>4{print $1, $2, int($3*19/20)}' "$C" >"$T/a2.txt"
check "(a') every weight reduced by 5 %: violating pose" VIOLATION \
      pose "$T/a2.txt" --x 1/2 --y 1/2 --u 0
grep -q '997505942/1000000000' "$LAST" \
  && ok "(a')  ...exact capture 0.997505942 < 1" seen \
  || bad "(a')  exact capture" "$(sed -n 's/.*EXACT captured weight = //p' "$LAST")" "0.997505942"
# (b) the same point deleted -- again within the cover's 3.1 % capture margin, so still valid;
#     the certification survives, but visibly harder: the band needs ~50 % more boxes.
awk 'NR==4{print $1-1; next} NR!=273{print}' "$C" >"$T/b.txt"
band "(b) one point deleted [still valid]" zero "$T/b.txt"
# (b') all 3621 weights halved: total 6.478, which cannot cover anything.
awk 'NR<=4{print} NR>4{print $1, $2, int($3/2)}' "$C" >"$T/b2.txt"
check "(b') all weights halved: violating pose" VIOLATION pose "$T/b2.txt" --x 1/2 --y 1/2 --u 0
# (c) the whole point set scaled by 0.995 (exactly: X -> 995 X, D -> 10^6), the container
#     unchanged -- i.e. the cover is asked to cover [0, 4/0.995]^2 = [0, 4.0201]^2.  This one is
#     spectacularly false: all of the shipped cover's corner-tile weight sits on the lines x = 3
#     and y = 3 (no point has both coordinates > 3.015), so after the shrink the corner pose
#     (7/2, 7/2, 0) captures nothing at all.
awk 'NR==2{print 1000000; next} NR>4{print $1*995, $2*995, $3; next}{print}' "$C" >"$T/c.txt"
check "(c) set scaled by 0.995: violating pose (7/2,7/2,0)" VIOLATION \
      pose "$T/c.txt" --x 7/2 --y 7/2 --u 0
grep -q '= 0/1000000000' "$LAST" \
  && ok "(c)  ...captures exactly nothing there" seen \
  || bad "(c)  exact capture" "$(sed -n 's/.*EXACT captured weight = //p' "$LAST")" "0"
# (d) one point moved by 0.01 -- inside the capture margin again, so still a valid cover.
awk 'NR==273{print $1+10, $2, $3; next}{print}' "$C" >"$T/d.txt"
band "(d) one point moved by 0.01 [still valid]" zero "$T/d.txt"
# (d') the whole set translated by 0.01: the strip x < 0.01 loses its cover and the corner pose
#      goes under.
awk 'NR>4{print $1+10, $2+10, $3; next}{print}' "$C" >"$T/d2.txt"
check "(d') whole set moved by 0.01: violating pose" VIOLATION pose "$T/d2.txt" --x 1/2 --y 1/2 --u 0
# (e) closed4_best_x103.txt: invalid at (3/2, 1461/2000, u = 1/40000) by 0.9703 (RUNG2.md sec 4.3)
if [ -r "$X103" ]; then
  check "(e) closed4_best_x103: violating pose exhibited" VIOLATION \
        pose "$X103" --x 3/2 --y 1461/2000 --u 1/40000
  grep -q '970282351/1000000000' "$LAST" \
    && ok "(e)  ...and the exact weight is RUNG2.md's 0.970282351" "seen" \
    || bad "(e)  exact weight" "$(sed -n 's/.*EXACT captured weight = //p' "$LAST")" \
           "970282351/1000000000"
  check "(e)  ...and closed4_best.txt likewise, 0.9420217" VIOLATION \
        pose "${X103%_x103.txt}.txt" --x 3/2 --y 1461/2000 --u 1/40000
else
  bad "(e) closed4_best_x103.txt input" "missing at $X103" "readable"
fi

echo "-- malformed input: refuse, do not judge (exit 2, ERROR:, no verdict word)"
awk 'NR==5{print $1, $2, -1; next}{print}' "$C" >"$T/e1.txt"
check "negative weight"                    ERROR cert "$T/e1.txt" $BAND
awk 'NR==5{print 4001, $2, $3; next}{print}' "$C" >"$T/e2.txt"
check "point outside the container"        ERROR cert "$T/e2.txt" $BAND
head -20 "$C" >"$T/e3.txt"
check "truncated point list"               ERROR cert "$T/e3.txt" $BAND
awk 'NR==1{print 4, 3; next}{print}' "$C" >"$T/e4.txt"
check "s_den does not divide s_num*D"      ERROR cert "$T/e4.txt" $BAND
{ cat "$C"; echo "1 2"; } >"$T/e5.txt"
check "trailing data after the last point" ERROR cert "$T/e5.txt" $BAND
awk 'NR==3{print 0; next}{print}' "$C" >"$T/e6.txt"
check "W = 0"                              ERROR cert "$T/e6.txt" $BAND
printf 'x y\n1000\n1\n0\n' >"$T/e7.txt"
check "non-integer header"                 ERROR cert "$T/e7.txt" $BAND
check "missing file"                       ERROR cert "$T/nope.txt" $BAND
check "u outside [0,1] in pose mode"       ERROR pose "$C" --x 1/2 --y 1/2 --u 3/2

printf '\n%d passed, %d failed, %d panics\n' "$pass" "$fail" "$panics"
[ "$fail" -eq 0 ] && [ "$panics" -eq 0 ]
