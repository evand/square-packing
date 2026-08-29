#!/bin/sh
# A verifier that always says VERIFIED is worthless.  Every file below is a mutation of a
# shipped certificate, built in a temporary directory; each must produce the stated verdict.
#
# Verdicts:  VERIFIED  exit 0 and a "VERIFIED:" line
#            REJECT    exit 0 and "NOT VERIFIED"   (a clean, reasoned rejection)
#            ERROR     non-zero exit, an ERROR: message, and NO verdict word in the output
#                      (malformed / illegal input: the verifier refuses to give a verdict)
#            PANIC     a Rust panic or abnormal exit.  Never acceptable, even on garbage input:
#                      a panic is not a rejection.  Any PANIC fails the whole script.
cd "$(dirname "$0")/.."
V=${VERIFY:-verify/target/release/verify}   # override with VERIFY=... to test another binary
C=certificates/s12_56points_3.8.txt
T=$(mktemp -d /tmp/rejtest.XXXXXX) || exit 1
trap 'rm -rf "$T"' EXIT
LAST=$T/last.out
pass=0; fail=0; panics=0
ok()  { pass=$((pass+1)); printf '  ok    %-52s -> %s\n' "$1" "$2"; }
bad() { fail=$((fail+1)); printf '  FAIL  %-52s -> %s (wanted %s)\n' "$1" "$2" "$3"; }
run() {  # file, n, [topk, witness-file]   -> sets $got, output in $LAST
  if [ -n "$3" ]; then $V "$1" "$2" 2000 4 "$3" "$4" >"$LAST" 2>&1; else $V "$1" "$2" 2000 4 0 >"$LAST" 2>&1; fi
  rc=$?
  if grep -q 'panicked at' "$LAST" || [ "$rc" -ge 100 ]; then got=PANIC; panics=$((panics+1))
  elif grep -q '^NOT VERIFIED' "$LAST"; then got=REJECT
  elif grep -q '^VERIFIED:' "$LAST" && [ "$rc" -eq 0 ]; then got=VERIFIED
  elif [ "$rc" -ne 0 ] && ! grep -q 'VERIFIED' "$LAST"; then got=ERROR
  else got="UNKNOWN(rc=$rc)"; fi
}
check() {  # name, expect(VERIFIED|REJECT|ERROR), file, n, [topk, witness-file]
  run "$3" "$4" "$5" "$6"
  if [ "$got" = "$2" ]; then ok "$1" "$got"; else bad "$1" "$got" "$2"
    [ "$got" = PANIC ] && sed 's/^/        | /' "$LAST" | head -3; fi
}
expect_out() {  # name, grep pattern that must appear in the LAST output
  if grep -q -- "$2" "$LAST"; then ok "$1" "seen"; else bad "$1" "not seen" "'$2'"; fi
}
angles() {  # name, witness file, above|below : every violated placement must be on that side of 45 deg
  if python3 - "$2" "$3" <<'EOF'
import sys, math
th = [float(l.split()[1]) for l in open(sys.argv[1]) if l.strip()]
if not th: print("      no violated placements written"); sys.exit(1)
print("      violated placements at %.2f..%.2f deg (%d)" % (math.degrees(min(th)), math.degrees(max(th)), len(th)))
side = (lambda t: t > math.pi/4) if sys.argv[2] == "above" else (lambda t: t < math.pi/4)
sys.exit(0 if all(side(t) for t in th) else 1)
EOF
  then ok "$1" "all $3 45 deg"; else bad "$1" "not all $3 45 deg" "all $3 45 deg"; fi
}

echo "-- baseline"
check "unmodified certificate"                    VERIFIED "$C" 12
# 1. drop one point  -> covering must break
awk 'NR==4{print $1-1; next} NR<=4{print} NR>5{print}' "$C" > $T/t1.txt
check "one point deleted"                         REJECT  $T/t1.txt 12
# 2. move one point.  NB small displacements may legitimately still verify: this
#    certificate has ~0.76% of slack (it survives scale-up to 1520/397), so the test
#    uses a displacement well beyond that margin.
awk 'NR<=4{print} NR==5{print $1+100, $2, $3} NR>5{print}' "$C" > $T/t2.txt
check "one point displaced by 100/400 = 0.25"     REJECT  $T/t2.txt 12
# 3. claim a larger container with the same points
awk 'NR==1{print 39,10; next}{print}' "$C" > $T/t3.txt
check "container inflated to 3.9"                 REJECT  $T/t3.txt 12
# 4. weights too heavy to prove anything for n=12
awk 'NR==3{print 4; next}{print}' "$C" > $T/t4.txt
check "weights inflated (total 14 > 12)"          REJECT  $T/t4.txt 12
# 5. same certificate cannot prove the STRONGER claim for n=11 (11.2 > 11)
check "same certificate claimed for n=11"         REJECT  "$C" 11
# 6. it should still prove the weaker claim for n=13
check "same certificate claimed for n=13"         VERIFIED "$C" 13

echo "-- symmetry reduction: non-D4-symmetric sets must be checked over all of [0,90) deg"
# 7. one extra weighted point off every symmetry axis: still a valid certificate (adding
#    weight never breaks covering; total 11.4 < 12) but no longer D4-symmetric, so the
#    verifier must take the full-range path instead of the [0,45] reduction.
awk 'NR==4{print $1+1; next}{print} END{print 1000, 700, 1}' "$C" > $T/t7.txt
check "extra point breaks D4 symmetry (valid)"    VERIFIED $T/t7.txt 12
expect_out "  ...symmetry check reports false"        'D4-symmetric atom set: false'
expect_out "  ...angles cover \[0,90) deg"           'angles cover \[0,90) deg'
expect_out "  ...all N bins enumerated (k=0..2000)"   'angles: k=0\.\.2000 '
# 8. A set that is fully covered for every angle in [0,45] and fails ONLY above 45 deg.
#    A verifier that wrongly applied the symmetry reduction to it would say VERIFIED.
#    Construction (all coordinates are the 56-point set scaled x10 with denominator 3972,
#    i.e. container 15200/3972 = 3800/993, just past the critical scale):
#      base:  D4-symmetric; the verifier's complete witness list shows it fails in exactly
#             four tiny families (centre extents < 0.05) at 39.80-39.85 deg, the four
#             rotations of one family near centre (2.60,1.44).  By the D4 symmetry that
#             the verifier checks exactly, the full failure set is these four families plus
#             their four mirror images at 90-39.8 = 50.2 deg (e.g. near (2.38,1.22)).
#      full:  base + the D4 orbit of r = (11636,7030)/3972 = (2.93,1.77): r lies within 0.46
#             of every centre of the (2.60,1.44) family, hence inside every one of its
#             squares (inscribed radius 0.5), and its images do the same for the other seven
#             families -> VERIFIED, i.e. the eight points repair everything.
#      A:     full minus the ONE image (8170,3564) = (2.06,0.90) that serves the 50.2-deg
#             family near (2.38,1.22).  Every family at angles <= 45 deg keeps its own repair
#             point, so A is covered on all of [0,45]; the 50.2-deg family near (2.38,1.22)
#             is left with only its four original points, and the remaining repair points
#             are all > 0.73 from its centres (> circumradius 0.707): its squares miss them.
#             Hence A fails only above 45 deg -> must be REJECTED, with every violated
#             placement above 45 deg.
#      A':    control: full minus r itself -> fails only at 39.8 deg (below 45).
#    The eight repair points weigh 8/5, so these four claims use n = 13 (total 12.8 < 13);
#    the covering check is identical for every n.
python3 - "$C" "$T" <<'EOF'
import sys
t = open(sys.argv[1]).read().split(); T = sys.argv[2]
pts = [(10*int(t[5+3*i]), 10*int(t[6+3*i])) for i in range(56)]
SD, r, rmirror = 15200, (11636, 7030), (8170, 3564)
orbit = sorted({q for p in [r, (SD-r[0], r[1]), (r[0], SD-r[1]), (SD-r[0], SD-r[1])] for q in (p, (p[1], p[0]))})
assert len(orbit) == 8 and rmirror in orbit
for name, extra in [("base", []), ("full", orbit), ("A", [p for p in orbit if p != rmirror]), ("Actl", [p for p in orbit if p != r])]:
    with open(f"{T}/t8_{name}.txt", "w") as f:
        f.write(f"3800 993\n3972\n5\n{56+len(extra)}\n")
        for x, y in pts + extra: f.write(f"{x} {y} 1\n")
EOF
check "symmetric base at D=3972 (fails at 39.8 deg)"  REJECT   $T/t8_base.txt 13 100000 $T/t8_base.sep
angles "  ...base violations lie below 45 deg"         $T/t8_base.sep below
check "base + 8 repair points (D4 orbit)"              VERIFIED $T/t8_full.txt 13
expect_out "  ...symmetric again: [0,45] path"          'angles cover \[0,45\] deg'
check "7 repair points: fails ONLY above 45 deg"       REJECT   $T/t8_A.txt 13 100000 $T/t8_A.sep
expect_out "  ...taken over [0,90)"                     'angles cover \[0,90) deg'
angles "  ...every violation above 45 deg"              $T/t8_A.sep above
check "control: other 7 points, fails below 45 deg"    REJECT   $T/t8_Actl.txt 13 100000 $T/t8_Actl.sep
angles "  ...every violation below 45 deg"              $T/t8_Actl.sep below

echo "-- input validation: illegal certificates get an ERROR, never a verdict, never a panic"
# 9. a point outside the container.  FORMAT.md defines the container as the closed square
#    [0,s]^2 and the points as (X/D, Y/D) in it; a point outside it can never be captured, so
#    it only pads the total weight and signals a broken file (wrong D, wrong s).  Refused.
awk 'NR==4{print $1+1; next}{print} END{print 5000, 5000, 1}' "$C" > $T/t9.txt
check "point outside the container"               ERROR   $T/t9.txt 12
awk 'NR==4{print $1+1; next}{print} END{print -1, 700, 1}' "$C" > $T/t9b.txt
check "point with a negative coordinate"          ERROR   $T/t9b.txt 12
# 10. weights.  A zero weight is harmless (the rest still covers -> VERIFIED; it also breaks
#     the symmetry, so this exercises the full-range path on a valid set).  A negative weight
#     is illegal: the reduction (Lean: hypothesis hw of packing_le_weight) needs w >= 0.
#     The last case is the exploit that a missing sign check allows: a negative weight where
#     no square can see it lowers the total below n without touching the covering, and would
#     "prove" s(11) >= 3.8 from a certificate of total weight 11.2.  It must not verify.
awk 'NR==4{print $1+1; next}{print} END{print 1000, 700, 0}' "$C" > $T/t10a.txt
check "zero-weight extra point (harmless)"        VERIFIED $T/t10a.txt 12
expect_out "  ...checked over [0,90) (symmetry broken)" 'angles cover \[0,90) deg'
awk 'NR<=4{print} NR==5{print $1, $2, -1} NR>5{print}' "$C" > $T/t10b.txt
check "negative weight on a point"                ERROR   $T/t10b.txt 12
awk 'NR==4{print $1+1; next}{print} END{print 5000, 5000, -2}' "$C" > $T/t10c.txt
check "negative weight outside container, n=11"   ERROR   $T/t10c.txt 11
# 11. malformed headers and files
printf '19 5\n400\n5\n' > $T/t11a.txt
check "header only, no point count"               ERROR   $T/t11a.txt 12
awk 'NR==1{print "19.0", 5; next}{print}' "$C" > $T/t11b.txt
check "non-integer field (19.0)"                  ERROR   $T/t11b.txt 12
awk 'NR==1{print 19, 5, 400; next}{print}' "$C" > $T/t11c.txt
check "wrong field count on line 1"               ERROR   $T/t11c.txt 12
awk 'NR==5{print $1, $2; next}{print}' "$C" > $T/t11d.txt
check "point line with 2 fields instead of 3"     ERROR   $T/t11d.txt 12
awk 'NR==4{print $1+1; next}{print}' "$C" > $T/t11e.txt
check "m = 57 but only 56 points"                 ERROR   $T/t11e.txt 12
awk 'NR==4{print $1-1; next}{print}' "$C" > $T/t11f.txt
check "m = 55 but 56 points (trailing data)"      ERROR   $T/t11f.txt 12
: > $T/t11g.txt
check "empty file"                                ERROR   $T/t11g.txt 12
check "missing file"                              ERROR   $T/does_not_exist.txt 12
$V >"$LAST" 2>&1; rc=$?
if grep -q 'panicked at' "$LAST" || [ "$rc" -ge 100 ]; then panics=$((panics+1)); bad "no arguments" PANIC ERROR
elif [ "$rc" -ne 0 ] && ! grep -q VERIFIED "$LAST"; then ok "no arguments" ERROR; else bad "no arguments" "rc=$rc" ERROR; fi
# 12. denominators.  The container side in atom units, s*D, must be an integer (s_den | s_num*D):
#     every legitimate rescaling keeps s*D fixed (FORMAT.md, "Scaling"), and the exact
#     integer arithmetic relies on it.  D=399 with s=19/5 gives s*D = 7581/5: refused.
#     Doubling D together with all coordinates is the same set of points: still verifies.
awk 'NR==2{print 399; next}{print}' "$C" > $T/t12a.txt
check "D=399 does not divide 19*D by 5"           ERROR   $T/t12a.txt 12
awk 'NR==2{print 800; next} NR>4{print 2*$1, 2*$2, $3; next}{print}' "$C" > $T/t12b.txt
check "D and all coordinates doubled (same set)"  VERIFIED $T/t12b.txt 12
# 13. no points at all: nothing is covered
printf '19 5\n400\n5\n0\n' > $T/t13.txt
check "empty point list"                          REJECT  $T/t13.txt 12

echo "-- multiset semantics: a point listed twice carries the sum of its weights"
# 14. FORMAT.md sums "the total weight of the certificate points lying in Q": two lines with the
#     same coordinates are two atoms whose weights add.  All weights doubled (W=10) and the
#     first point split into two lines of weight 1 is the same weighting -> VERIFIED.  Drop
#     one of the two lines and that point only weighs 1/10: squares through it drop to 9/10.
awk 'NR==3{print 10; next} NR==4{print $1+1; next} NR==5{print $1, $2, 1; print $1, $2, 1; next} NR>5{print $1, $2, 2; next} {print}' "$C" > $T/t14a.txt
check "duplicated point, weight split 1+1 of 10"  VERIFIED $T/t14a.txt 12
awk 'NR==3{print 10; next} NR==5{print $1, $2, 1; next} NR>5{print $1, $2, 2; next} {print}' "$C" > $T/t14b.txt
check "same file with the duplicate line removed" REJECT  $T/t14b.txt 12

echo "-- centre box of an angle bin above 45 deg must use the bin's minimum width"
# 15. For a bin [theta_k, theta_k+1] the admissible-centre box must be taken for the SMALLEST
#     bounding-box width w = cos+sin over the bin.  w increases up to 45 deg and decreases
#     after it, so above 45 deg the minimum is at theta_k+1, not theta_k.  Using theta_k there
#     silently drops the centres closest to the container edge for the whole bin.  On the
#     inflated-container mutant (non-symmetric, so bins above 45 deg are actually run) bin
#     k=1530 (74.83 deg) has a placement in that omitted strip -- a sigma-square in the
#     top-right corner at about (3.287, 3.287) -- covering only 1/5, while every centre of
#     the too-small box covers >= 2/5.  A verifier with the wrong box reports 2/5 for the bin.
#     (The overall verdict is REJECT either way: the neighbouring bin k=1531 sees the same
#     corner.  Consecutive bins overlap, which is why this bug never flipped a verdict, but
#     it does break the stated invariant "every admissible centre of every bin is checked".)
check "inflated container, witness mode"          REJECT  $T/t3.txt 12 20 $T/t15.sep
if python3 - $T/t15.sep <<'EOF'
import sys, math
best = {}
for l in open(sys.argv[1]):
    v, th, cx, cy = map(float, l.split()[:4]); k = round(math.tan(th/2)*2000)
    best[k] = min(best.get(k, 9), v)
print("      per-bin minimum near 74.8 deg:", "  ".join("k=%d:%.1f" % (k, best[k]) for k in (1529, 1530, 1531) if k in best))
sys.exit(0 if best.get(1530, 9) <= 0.2 + 1e-9 else 1)
EOF
then ok "  ...bin k=1530 reaches its edge placement (1/5)" "1/5"; else bad "  ...bin k=1530 reaches its edge placement (1/5)" "min > 1/5" "1/5"; fi

echo "-- branch certificates: the region / lambda / k trailer (FORMAT.md, 'Branch certificates')"
# 16. A trailer with lambda = 0 changes nothing: same verdict as the plain file.
( cat "$C"; printf 'region corner 6 5\nlambda 0\nk 4\n' ) > $T/t16a.txt
check "trailer with lambda=0, k=4"                VERIFIED $T/t16a.txt 12
expect_out "  ...reported as a branch certificate"  "^VERIFIED: (branch k=4)"
#     r = 13/10: the corner box could hold two centres ((2r-1)^2 = 2.56 >= 2): refused.
( cat "$C"; printf 'region corner 13 10\nlambda 0\nk 4\n' ) > $T/t16b.txt
check "r=13/10 (two squares could share a box)"    ERROR   $T/t16b.txt 12
( cat "$C"; printf 'region corner 6 5\nlambda 0\nk 5\n' ) > $T/t16c.txt
check "k=5 (more squares than boxes)"             ERROR   $T/t16c.txt 12
( cat "$C"; printf 'region square 6 5\nlambda 0\nk 4\n' ) > $T/t16d.txt
check "unknown region kind"                       ERROR   $T/t16d.txt 12
( cat "$C"; printf 'region corner 6 5\nlambda 0\nk 4\n7\n' ) > $T/t16e.txt
check "data after the trailer"                    ERROR   $T/t16e.txt 12
( cat "$C"; printf 'region corner 6 5\nk 4\nlambda 0\n' ) > $T/t16f.txt
check "trailer keywords out of order"             ERROR   $T/t16f.txt 12
#     lambda = +1 (W=5): squares centred in a corner box must capture 2 -- they do not.
( cat "$C"; printf 'region corner 6 5\nlambda 5\nk 4\n' ) > $T/t16g.txt
check "lambda=+1, k=4: corners must capture 2"    REJECT  $T/t16g.txt 12 8 $T/t16g.sep
if python3 - $T/t16g.sep <<'EOF2'
import sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip()]
s = 19/5; r = 6/5
ok = rows and all(len(q) == 5 for q in rows) and all(q[4] in ('1', '2', '3', '4') for q in rows)   # flag = index of the box met
nin = sum((float(q[2]) <= r + 1e-6 or float(q[2]) >= s - r - 1e-6) and (float(q[3]) <= r + 1e-6 or float(q[3]) >= s - r - 1e-6) for q in rows)
print("      %d witnesses, all flagged with a box: %s, centred in a corner box: %d (cells straddling the box edge fall back to their centroid)" % (len(rows), ok, nin))
sys.exit(0 if ok and nin >= 0.95 * len(rows) else 1)
EOF2
then ok "  ...witnesses flagged 1, >=95% centred in the boxes" "yes"; else bad "  ...witnesses flagged 1, >=95% centred in the boxes" "no" "yes"; fi
#     lambda = -1, k = 4: the covering holds (corners may capture 0) but W - lambda*k = 15.2 >= 12.
( cat "$C"; printf 'region corner 6 5\nlambda -5\nk 4\n' ) > $T/t16h.txt
check "lambda=-1, k=4: weight bound fails"        REJECT  $T/t16h.txt 12
expect_out "  ...for the stated reason"            "^WEIGHT NOT"
#     lambda = -1, k = 0: corner poses are free, everything else still covered: verifies as branch k=0.
( cat "$C"; printf 'region corner 6 5\nlambda -5\nk 0\n' ) > $T/t16i.txt
check "lambda=-1, k=0 (corner free)"              VERIFIED $T/t16i.txt 12
expect_out "  ...reported as branch k=0"           "^VERIFIED: (branch k=0)"
#     the same with the points inside the corner boxes deleted: squares centred OUTSIDE the
#     boxes still reach into the corners, so the plain covering breaks -> REJECT, with flag-0
#     witnesses only (flag-1 violations cannot occur at lambda = -1).
awk 'NR<=4{print; next} { if (($1<=480 || $1>=1040) && ($2<=480 || $2>=1040)) next; print }' "$C" > $T/t16j.txt
n=$(($(wc -l < $T/t16j.txt) - 4)); awk -v n=$n 'NR==4{print n; next}{print}' $T/t16j.txt > $T/t16k.txt
( cat $T/t16k.txt; printf 'region corner 6 5\nlambda -5\nk 0\n' ) > $T/t16l.txt
check "k=0 with the corner points deleted"        REJECT  $T/t16l.txt 12 8 $T/t16l.sep
if python3 - $T/t16l.sep <<'EOF2'
import sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip()]
ok = rows and all(len(q) == 5 and q[4] == '0' for q in rows)
print("      %d witnesses, all flagged 0: %s" % (len(rows), ok)); sys.exit(0 if ok else 1)
EOF2
then ok "  ...witnesses all flagged 0" "yes"; else bad "  ...witnesses all flagged 0" "no" "yes"; fi

echo "-- per-box trailer: lambda L1 L2 L3 L4 / k K1 K2 K3 K4"
# 17. Four multipliers, one per corner box (box 1 = [0,r]^2, 2 = bottom-right, 3 = top-left, 4 = top-right).
( cat "$C"; printf 'region corner 6 5\nlambda 0 0 0 0\nk 1 1 1 1\n' ) > $T/t17a.txt
check "per-box, all lambda=0, pattern 1111"       VERIFIED $T/t17a.txt 12
expect_out "  ...reported with the pattern"         "^VERIFIED: (branch k=1111)"
#     unequal lambdas break the D4 symmetry of the claim: the full angle range must be swept
( cat "$C"; printf 'region corner 6 5\nlambda -5 -5 0 0\nk 0 0 1 1\n' ) > $T/t17b.txt
check "per-box, boxes 1,2 free, 3,4 threshold 1"  VERIFIED $T/t17b.txt 12
expect_out "  ...swept over [0,90) deg"             "angles cover \[0,90) deg"
#     box 1 must capture 2: rejected, and every witness names box 1
( cat "$C"; printf 'region corner 6 5\nlambda 5 0 0 0\nk 1 0 0 0\n' ) > $T/t17c.txt
check "per-box, lambda_1=+1: box 1 must capture 2" REJECT  $T/t17c.txt 12 8 $T/t17c.sep
if python3 - $T/t17c.sep <<'EOF2'
import sys
rows = [l.split() for l in open(sys.argv[1]) if l.strip()]
ok = rows and all(len(q) == 5 and q[4] == '1' for q in rows)
nin = sum(float(q[2]) <= 1.2 + 1e-6 and float(q[3]) <= 1.2 + 1e-6 for q in rows)
print("      %d witnesses, all flagged box 1: %s, %d centred in box 1" % (len(rows), ok, nin))
sys.exit(0 if ok and nin >= 0.95 * len(rows) else 1)
EOF2
then ok "  ...witnesses all name box 1, >=95% inside it" "yes"; else bad "  ...witnesses all name box 1, >=95% inside it" "no" "yes"; fi
#     weight check uses sum_j lambda_j k_j: 11.2 - (-1)(1+1) = 13.2 >= 12
( cat "$C"; printf 'region corner 6 5\nlambda -5 -5 -5 -5\nk 1 1 0 0\n' ) > $T/t17d.txt
check "per-box, lambda=-1 on two occupied boxes"   REJECT  $T/t17d.txt 12
expect_out "  ...for the stated reason"             "^WEIGHT NOT"
( cat "$C"; printf 'region corner 6 5\nlambda 0 0\nk 1 1\n' ) > $T/t17e.txt
check "two lambdas (neither 1 nor 4)"             ERROR   $T/t17e.txt 12
( cat "$C"; printf 'region corner 6 5\nlambda 0 0 0 0\nk 1 1 2 0\n' ) > $T/t17f.txt
check "per-box k=2"                               ERROR   $T/t17f.txt 12
( cat "$C"; printf 'region corner 6 5\nlambda 0 0 0 0\nk 1 1 1\n' ) > $T/t17g.txt
check "four lambdas, three k"                     ERROR   $T/t17g.txt 12

echo "  ---- $pass passed, $fail failed, $panics panics"
[ "$fail" -eq 0 ] && [ "$panics" -eq 0 ]
