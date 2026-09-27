#!/usr/bin/env bash
# zmx2_tests.sh -- tests of the independent mixed-cover checker zmx2 (search/ZMX2.md sec 8).
# Run from s12/:  taskset -c 0-9 bash search/zmx2_tests.sh [--quick]
# Writes scratch files to $ZT (default runs/zmx2_tests/).  Exit 0 iff every test passes.
set -u
cd "$(dirname "$0")/.."
QUICK=0; [ "${1:-}" = "--quick" ] && QUICK=1
ZT=${ZT:-runs/zmx2_tests}; mkdir -p "$ZT"
TH=${TH:-10}
(cd verify2 && cargo build --release --bin zmx2 2>&1 | grep -E '^(error|warning: unused)' ) && { echo "build failed"; exit 1; }
Z=verify2/target/release/zmx2
T="python3 search/zmx2_tools.py"
C=runs/line-cover_m5_candidate_x1003.txt
S13=certificates/rung2/s13_closed_cover_4.txt
npass=0; nfail=0
ok()   { echo "PASS  $*"; npass=$((npass+1)); }
bad()  { echo "FAIL  $*"; nfail=$((nfail+1)); }
verdict() { grep -E '^(VERIFIED|VERIFIED-D4|NOT VERIFIED|REGION CLEAN|ERROR|INCOMPLETE)' "$1" | head -1 | sed -E 's/[:(].*//; s/ +$//'; }
fminpose() { grep UNCERT "$1" | awk -F'float-min ' '{print $2}' | sort -n | head -1; }

echo "== T0 zmcheck untouched (source and manifest identical to HEAD)"
if git diff --quiet HEAD -- verify2/src/main.rs verify2/Cargo.toml verify2/Cargo.lock; then ok "verify2/src/main.rs, Cargo.toml, Cargo.lock unchanged"; else bad "zmcheck sources changed"; fi

echo "== T1 malformed / unsupported files are refused (ERROR, exit 2)"
mk() { printf "mixed 1\n5 1\n1000\n100\n%b" "$2" > "$ZT/bad_$1.txt"; }
mk negw     "1\n100 100 -1\n0\n0\n"
mk outside  "1\n6000 100 1\n0\n0\n"
mk degen    "0\n1\n1000 10 1000 10 5\n0\n"
mk diag     "0\n1\n1000 10 1010 20 5\n0\n"
mk polygon  "0\n0\n1\n3 1 0 0 1000 0 0 1000\n"
mk trailing "1\n100 100 1\n0\n0\n7\n"
mk short    "2\n100 100 1\n"
for b in negw outside degen diag polygon trailing short; do
  out=$($Z info "$ZT/bad_$b.txt" 2>&1); rc=$?
  if [ $rc = 2 ] && echo "$out" | grep -q '^ERROR'; then ok "refuses $b: $(echo "$out" | head -1 | cut -c1-70)"; else bad "accepted $b (rc $rc)"; fi
done

echo "== T2 D4 invariance is checked exactly"
out=$($Z d4 $C); [ $? = 0 ] && ok "candidate: $out" || bad "candidate D4"
$T perturb $C $ZT/cand_bump.txt --op bump >/dev/null
out=$($Z cert $ZT/cand_bump.txt --d4 --threads 1 2>&1); rc=$?
[ $rc = 2 ] && echo "$out" | grep -q 'not D4-invariant' && ok "one weight +1 refused under --d4" || bad "bumped cover accepted under --d4"

echo "== T3 soundness harness: exact mu at random admissible rational poses >= box bound"
N=400; [ $QUICK = 1 ] && N=80
$T sound $C --n $N --poses 20 --seed 11 > $ZT/sound_c.txt; tail -1 $ZT/sound_c.txt
grep -q ' 0 FAIL' $ZT/sound_c.txt && ok "candidate random boxes" || bad "candidate random boxes"
$T sound $C --n $N --poses 20 --seed 12 --refl > $ZT/sound_cr.txt; tail -1 $ZT/sound_cr.txt
grep -q ' 0 FAIL' $ZT/sound_cr.txt && ok "reflected candidate random boxes" || bad "reflected candidate"
python3 - "$S13" "$ZT/s13m.txt" <<'EOF'
import sys
t = []
for l in open(sys.argv[1]):
    t += l.split('#')[0].split()
s_num, s_den, D, W, m = map(int, t[:5])
with open(sys.argv[2], 'w') as f:
    f.write('mixed 1\n%d %d\n%d\n%d\n%d\n' % (s_num, s_den, D, W, m))
    for i in range(m):
        f.write(' '.join(t[5 + 3 * i:8 + 3 * i]) + '\n')
    f.write('0\n0\n')
EOF
$T sound $ZT/s13m.txt --n $N --poses 20 --seed 13 > $ZT/sound_s13.txt; tail -1 $ZT/sound_s13.txt
grep -q ' 0 FAIL' $ZT/sound_s13.txt && ok "s(13) point cover (line atoms) random boxes" || bad "s(13) random boxes"
$Z cert $C --d4 --threads $TH --tight 0.0005 --dump-tight $ZT/tight_c.txt > /dev/null
$T tight $ZT/tight_c.txt $C --max $((N/2)) --poses 12 --seed 14 > $ZT/tight_c.out; tail -1 $ZT/tight_c.out
grep -q ' 0 FAIL' $ZT/tight_c.out && ok "candidate: tightest certified leaves (bound < 1.0005)" || bad "candidate tight leaves"
$Z cert $ZT/s13m.txt --d4 --threads $TH --tight 0.001 --dump-tight $ZT/tight_s13.txt > /dev/null
$T tight $ZT/tight_s13.txt $ZT/s13m.txt --max $((N/2)) --poses 12 --seed 15 > $ZT/tight_s13.out; tail -1 $ZT/tight_s13.out
grep -q ' 0 FAIL' $ZT/tight_s13.out && ok "s(13): tightest certified leaves (bound < 1.001)" || bad "s(13) tight leaves"

echo "== T4 toy mixed covers: certify at 1.02 x (float min)^-1, reject at 0.98 x"
for spec in "2 lines 1 1.662215" "3 lines 1 1.656854" "4 wave 1 2.011632" "3 wavec 1 2.243322"; do
  set -- $spec; s=$1; k=$2; extra=""
  [ $k = wavec ] && { k=wave; extra="--centre-pts --centre-w 0.2"; }
  $T toy $ZT/toy_${s}_$2.txt --s $s --kind $k --rho $3 $extra > /dev/null
  for tag in ok rej; do
    m=1.02; [ $tag = rej ] && m=0.98
    f=$(python3 -c "from fractions import Fraction as F; print(F($m/$4).limit_denominator(10**6))")
    $T perturb $ZT/toy_${s}_$2.txt $ZT/toy_${s}_$2_$tag.txt --op scale --f $f > /dev/null
    for mode in d4 full; do
      $Z cert $ZT/toy_${s}_$2_$tag.txt --$mode --threads $TH > $ZT/toy_${s}_$2_${tag}_$mode.log 2>&1
      v=$(verdict $ZT/toy_${s}_$2_${tag}_$mode.log)
      if [ $tag = ok ]; then
        [ "$v" = VERIFIED ] || [ "$v" = VERIFIED-D4 ] && ok "toy s=$s $2 x$m $mode: $v" || bad "toy s=$s $2 x$m $mode: $v"
      else
        [ "$v" = "NOT VERIFIED" ] && ok "toy s=$s $2 x$m $mode: refused, float-min $(fminpose $ZT/toy_${s}_$2_${tag}_$mode.log)" || bad "toy s=$s $2 x$m $mode: $v"
      fi
    done
  done
done

echo "== T5 rejection tests on the candidate"
$T perturb $C $ZT/cand_x0994.txt --op scale --f 994/1000 > /dev/null
$Z cert $ZT/cand_x0994.txt --d4 --threads $TH > $ZT/cand_x0994.log 2>&1
v=$(verdict $ZT/cand_x0994.log); p=$(fminpose $ZT/cand_x0994.log)
[ "$v" = "NOT VERIFIED" ] && ok "candidate x0.994: refused; deepest float-min $p (known dip (0.573,1.441,9.13 deg))" || bad "candidate x0.994: $v"
ex=$($T mu $ZT/cand_x0994.txt 573593/1000000 1440624/1000000 805601/10000000)
echo "$ex" | awk '{exit !($5 < 1)}' && ok "  exact mu there: $(echo $ex | awk '{print $5}') < 1" || bad "  exact mu there not < 1: $ex"
$T perturb $C $ZT/cand_x09955.txt --op scale --f 9955/10000 > /dev/null
$Z cert $ZT/cand_x09955.txt --d4 --threads $TH > $ZT/cand_x09955.log 2>&1
v=$(verdict $ZT/cand_x09955.log)
[ "$v" = VERIFIED-D4 ] && ok "candidate x0.9955 (true min ~1.0009): still VERIFIED-D4 (checker loses < 0.1 %)" || bad "candidate x0.9955: $v"
if [ $QUICK = 0 ]; then
  $T perturb $C $ZT/cand_hole.txt --op zero-seg --line x --at 2 --lo 1.40 --hi 1.60 > /dev/null
  $Z cert $ZT/cand_hole.txt --full --threads $TH > $ZT/cand_hole.log 2>&1
  v=$(verdict $ZT/cand_hole.log); p=$(fminpose $ZT/cand_hole.log)
  [ "$v" = "NOT VERIFIED" ] && ok "candidate, x=2 zeroed on y in [1.4,1.6] (--full): refused; deepest float-min $p (germ (1.5,1.5), theta -> 0+)" || bad "hole: $v"
  e0=$($T mu $ZT/cand_hole.txt 3/2 3/2 0 | awk '{print $5}'); e1=$($T mu $ZT/cand_hole.txt 1500024/1000000 1500016/1000000 264/10000000 | awk '{print $5}')
  python3 -c "import sys; sys.exit(not ($e0 >= 1 > $e1))" && ok "  exact mu: $e0 at theta=0 (edges count), $e1 at theta=0.003 deg" || bad "  exact mu $e0 / $e1"
fi

echo "== T6 point covers: agreement with zmcheck on the shipped s(13) cover"
$Z cert $S13 --full --threads $TH --xlo 12 --xhi 12 > $ZT/s13col_zmx2.log 2>&1
v=$(verdict $ZT/s13col_zmx2.log)
[ "$v" = "REGION CLEAN" ] && ok "zmx2 --full column c_x in [1.2,1.3] (320 roots): clean; $(grep '^done' $ZT/s13col_zmx2.log | sed 's/.*boxes/boxes/' | cut -d, -f1)" || bad "zmx2 s13 column: $v"
if [ -x verify2/target/release/zmcheck ]; then
  verify2/target/release/zmcheck cert $S13 --threads $TH --xlo 1.2 --xhi 1.2 > $ZT/s13col_zmcheck.log 2>&1
  grep -q 'UNCERTIFIED 0' $ZT/s13col_zmcheck.log && ok "zmcheck same column (320 roots): 0 uncertified; $(grep 'boxes' $ZT/s13col_zmcheck.log | head -1)" || bad "zmcheck s13 column"
fi
$Z cert $S13 --full --threads $TH --xlo 12 --xhi 12 --no-atoms > $ZT/s13col_noatoms.log 2>&1
v=$(verdict $ZT/s13col_noatoms.log)
[ "$v" = "NOT VERIFIED" ] && ok "zmx2 --no-atoms (points only, no pair lemma for line points): refused at one-cut germs, float-min $(fminpose $ZT/s13col_noatoms.log)" || bad "no-atoms column: $v"
$Z cert $S13 --d4 --threads $TH > $ZT/s13_d4.log 2>&1
v=$(verdict $ZT/s13_d4.log)
[ "$v" = VERIFIED-D4 ] && ok "zmx2 --d4 on s(13) (whole region): VERIFIED-D4; $(grep '^done' $ZT/s13_d4.log | sed 's/.*boxes/boxes/' | cut -d, -f1)" || bad "s13 d4: $v"
$T perturb $ZT/s13m.txt $ZT/s13_x0975.txt --op scale --f 975/1000 > /dev/null
$Z cert $ZT/s13_x0975.txt --d4 --threads $TH > $ZT/s13_x0975.log 2>&1
v=$(verdict $ZT/s13_x0975.log)
[ "$v" = "NOT VERIFIED" ] && ok "s(13) x0.975 (min ~0.995): refused, float-min $(fminpose $ZT/s13_x0975.log)" || bad "s13 x0.975: $v"

echo "== T7 the candidate"
$Z cert $C --d4 --threads $TH > $ZT/cand_d4.log 2>&1
v=$(verdict $ZT/cand_d4.log)
[ "$v" = VERIFIED-D4 ] && ok "candidate --d4: $v; $(grep '^done' $ZT/cand_d4.log | sed 's/.*boxes/boxes/' | cut -d, -f1)" || bad "candidate d4: $v"

echo
echo "zmx2_tests: $npass passed, $nfail failed"
[ $nfail = 0 ]
