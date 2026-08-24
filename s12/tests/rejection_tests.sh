#!/bin/sh
# A verifier that always says VERIFIED is worthless.  These mutations MUST be rejected.
cd "$(dirname "$0")/.."
V=verify/target/release/verify
C=certificates/s12_56points_3.8.txt
pass=0; fail=0
check() {  # name, expect(VERIFIED|REJECT), file, n
  out=$($V "$3" "$4" 2000 4 0 2>&1)
  if echo "$out" | grep -q "NOT VERIFIED"; then got=REJECT
  elif echo "$out" | grep -q "^VERIFIED:"; then got=VERIFIED
  else got=REJECT; fi
  if [ "$got" = "$2" ]; then pass=$((pass+1)); printf '  ok    %-42s -> %s\n' "$1" "$got"
  else fail=$((fail+1)); printf '  FAIL  %-42s -> %s (wanted %s)\n' "$1" "$got" "$2"; fi
}
mut() { python3 - "$@" ; }
check "unmodified certificate"                    VERIFIED "$C" 12
# 1. drop one point  -> covering must break
awk 'NR==4{print $1-1; next} NR<=4{print} NR>5{print}' "$C" > /tmp/t1.txt
check "one point deleted"                         REJECT  /tmp/t1.txt 12
# 2. move one point.  NB small displacements may legitimately still verify: this
#    certificate has ~0.76% of slack (it survives scale-up to 1520/397), so the test
#    uses a displacement well beyond that margin.
awk 'NR<=4{print} NR==5{print $1+100, $2, $3} NR>5{print}' "$C" > /tmp/t2.txt
check "one point displaced by 100/400 = 0.25"     REJECT  /tmp/t2.txt 12
# 3. claim a larger container with the same points
awk 'NR==1{print 39,10; next}{print}' "$C" > /tmp/t3.txt
check "container inflated to 3.9"                 REJECT  /tmp/t3.txt 12
# 4. weights too heavy to prove anything for n=12
awk 'NR==3{print 4; next}{print}' "$C" > /tmp/t4.txt
check "weights inflated (total 14 > 12)"          REJECT  /tmp/t4.txt 12
# 5. same certificate cannot prove the STRONGER claim for n=11 (11.2 > 11)
check "same certificate claimed for n=11"         REJECT  "$C" 11
# 6. it should still prove the weaker claim for n=13
check "same certificate claimed for n=13"         VERIFIED "$C" 13
echo "  ---- $pass passed, $fail failed"
[ "$fail" -eq 0 ]
