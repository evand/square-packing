#!/usr/bin/env bash
# zmx2_area_tests.sh -- tests of zmx2's area densities and theta = 0 mode (search/ZMX2_AREA.md sec 10).
# Run from s12/:  taskset -c 10-11 bash search/zmx2_area_tests.sh [--quick]     (TH threads, default 2)
# Also run at the end of search/zmx2_tests.sh.  Scratch files in $ZT (default runs/zmx2_area_tests/).
# Exit 0 iff every test passes.
set -u
cd "$(dirname "$0")/.."
QUICK=0; [ "${1:-}" = "--quick" ] && QUICK=1
ZT=${ZT:-runs/zmx2_area_tests}; mkdir -p "$ZT"
TH=${TH:-2}
(cd verify2 && cargo build --release --bin zmx2 2>&1 | grep -E '^(error|warning: unused)' ) && { echo "build failed"; exit 1; }
Z=verify2/target/release/zmx2
T="python3 search/zmx2_tools.py"
K7=certificates/k2m3/L4_k02_box7.txt
C21=certificates/s21/s21_mixed_cover_5.txt
npass=0; nfail=0
ok()   { echo "PASS  $*"; npass=$((npass+1)); }
bad()  { echo "FAIL  $*"; nfail=$((nfail+1)); }
verdict() { grep -E '^(VERIFIED|VERIFIED-D4|NOT VERIFIED|REGION CLEAN|ERROR|INCOMPLETE)' "$1" | head -1 | sed -E 's/[:(].*//; s/ +$//'; }
census() { grep -E '^done' "$1" | sed 's/.*roots/roots/; s/, max depth.*//'; }

echo "== A1 parser: axis-parallel rectangles (ccw) accepted, everything else refused"
mk() { printf "mixed 1\n5 1\n1000\n100\n%b" "$2" > "$ZT/area_$1.txt"; }
mk rect_ok    "0\n0\n1\n4 100 0 0 1000 0 1000 1000 0 1000\n"
mk rect_rot   "0\n0\n1\n4 100 1000 1000 0 1000 0 0 1000 0\n"
mk rect_cw    "0\n0\n1\n4 100 0 0 0 1000 1000 1000 1000 0\n"
mk rect_degen "0\n0\n1\n4 100 0 0 1000 0 1000 0 0 0\n"
mk rect_tri   "0\n0\n1\n3 100 0 0 1000 0 0 1000\n"
mk rect_quad  "0\n0\n1\n4 100 0 0 1000 0 1000 1000 0 900\n"
mk rect_big   "0\n0\n1\n4 1125899906842624 0 0 1000 0 1000 1000 0 1000\n"
mk rect_out   "0\n0\n1\n4 100 0 0 6000 0 6000 1000 0 1000\n"
mk rect_negw  "0\n0\n1\n4 -1 0 0 1000 0 1000 1000 0 1000\n"
for b in rect_ok rect_rot; do
  out=$($Z info "$ZT/area_$b.txt" 2>&1); rc=$?
  [ $rc = 0 ] && echo "$out" | grep -q 'area density' && ok "accepts $b" || bad "refused $b: $out"
done
for b in rect_cw rect_degen rect_tri rect_quad rect_big rect_out rect_negw; do
  out=$($Z info "$ZT/area_$b.txt" 2>&1); rc=$?
  if [ $rc = 2 ] && echo "$out" | grep -q '^ERROR'; then ok "refuses $b: $(echo "$out" | head -1 | cut -c1-70)"; else bad "accepted $b (rc $rc)"; fi
done

echo "== A2 D4 invariance with a rectangle (exact)"
out=$($Z d4 $K7); [ $? = 0 ] && ok "k2m3 box: $out" || bad "k2m3 box D4"
$T perturb $K7 $ZT/k7_bump.txt --op seg-delta --at 9,13,9,12 --dw 1 > /dev/null
out=$($Z d4 $ZT/k7_bump.txt 2>&1); rc=$?
[ $rc = 2 ] && echo "$out" | grep -q 'not D4-invariant' && ok "one segment of the box +1: refused under --d4" || bad "bumped box accepted"
printf "mixed 1\n7 1\n5\n100\n0\n0\n1\n4 100 9 9 26 9 26 25 9 25\n" > $ZT/area_asym.txt
out=$($Z d4 $ZT/area_asym.txt 2>&1); rc=$?
[ $rc = 2 ] && echo "$out" | grep -q 'rectangles' && ok "non-symmetric rectangle refused under --d4: $(echo "$out" | cut -c1-80)" || bad "asymmetric rectangle accepted"

echo "== A3 toy: pure Lebesgue on the container (mu = 1 at every pose: zero margin everywhere)"
for s in 2 3; do
  $T toy $ZT/leb$s.txt --s $s --kind lebesgue > /dev/null
  for mode in d4 full; do
    $Z cert $ZT/leb$s.txt --$mode --threads $TH > $ZT/leb${s}_$mode.log 2>&1
    v=$(verdict $ZT/leb${s}_$mode.log)
    [ "$v" = VERIFIED ] || [ "$v" = VERIFIED-D4 ] && ok "lebesgue s=$s cert --$mode: $v ($(census $ZT/leb${s}_$mode.log))" || bad "lebesgue s=$s --$mode: $v"
    $Z cert0 $ZT/leb$s.txt --$mode > $ZT/leb${s}_0$mode.log 2>&1
    v=$(grep -E '^(VERIFIED|NOT)' $ZT/leb${s}_0$mode.log | head -1 | cut -c1-24)
    echo "$v" | grep -q '^VERIFIED' && ok "lebesgue s=$s cert0 --$mode: $v" || bad "lebesgue s=$s cert0 --$mode: $v"
  done
  # density lowered by 1e-6 (W = 10^9, mass s^2 W): every pose has mu = 1 - 1e-6
  $T perturb $ZT/leb$s.txt $ZT/leb${s}m.txt --op rect-delta --dw -$((s*s*1000)) > /dev/null
  $Z cert $ZT/leb${s}m.txt --full --threads $TH --uncert-cap 5 > $ZT/leb${s}m.log 2>&1
  v=$(verdict $ZT/leb${s}m.log)
  r=$($T uncert $ZT/leb${s}m.log $ZT/leb${s}m.txt --max 5)
  [ "$v" = "NOT VERIFIED" ] && echo "$r" | awk '{exit !($7 < 1)}' && ok "lebesgue s=$s, rho = 1 - 1e-6: refused; $r" || bad "lebesgue s=$s - 1e-6: $v / $r"
  $Z cert0 $ZT/leb${s}m.txt --full --uncert-cap 5 > $ZT/leb${s}m0.log 2>&1
  grep -q '^NOT VERIFIED' $ZT/leb${s}m0.log && ok "  cert0: refused ($(grep '^worst' $ZT/leb${s}m0.log | cut -c1-90))" || bad "  cert0 accepted rho - 1e-6"
done

echo "== A4 the k2m3 box at theta = 0 (cert0, Lemma Z0): D4 region and whole square"
for mode in d4 full; do
  $Z cert0 $K7 --$mode > $ZT/k7_0$mode.log 2>&1
  v=$(grep -E '^(VERIFIED|NOT)' $ZT/k7_0$mode.log | head -1 | cut -c1-30)
  echo "$v" | grep -q '^VERIFIED' && ok "k2m3 box cert0 --$mode: $v; $(grep '^done0' $ZT/k7_0$mode.log)" || bad "k2m3 box cert0 --$mode: $v"
done

echo "== A5 rejection: the box with a piece lowered (theta = 0 and theta > 0)"
# area density lowered by 1e-6 (W = 10^12, mass 11.56): squares inside the area square have mu = 1 - 1e-6
$T perturb $K7 $ZT/k7_rho.txt --op rect-delta --dw -11560000 > /dev/null
$Z cert0 $ZT/k7_rho.txt --d4 --uncert-cap 3 > $ZT/k7_rho0.log 2>&1
grep -q '^NOT VERIFIED' $ZT/k7_rho0.log && ok "rho - 1e-6, cert0: refused ($(grep '^worst' $ZT/k7_rho0.log | cut -c1-100))" || bad "rho - 1e-6 accepted by cert0"
$Z cert $ZT/k7_rho.txt --d4 --threads $TH --umin 10 --xlo 30 --xhi 31 --ylo 30 --yhi 31 --uncert-cap 3 > $ZT/k7_rho.log 2>&1
v=$(verdict $ZT/k7_rho.log); r=$($T uncert $ZT/k7_rho.log $ZT/k7_rho.txt --max 6)
[ "$v" = "NOT VERIFIED" ] && echo "$r" | awk '{exit !($7 < 1)}' && ok "rho - 1e-6, cert (centres [3,3.2]^2, all angles): refused; $r" || bad "rho - 1e-6 cert: $v / $r"
# a segment piece in the wall band lowered by 1e-4: Q = [0,1] x [2.5,3.5] has mu = 1 - 1e-4 at theta = 0
$T perturb $K7 $ZT/k7_v4.txt --op seg-delta --at 4,14,4,15 --dw -100000000 > /dev/null
$Z cert0 $ZT/k7_v4.txt --full --uncert-cap 3 > $ZT/k7_v40.log 2>&1
ex=$($T mu $ZT/k7_v4.txt 1/2 3 0 | awk '{print $5}')
grep -q '^NOT VERIFIED' $ZT/k7_v40.log && python3 -c "import sys; sys.exit(not ($ex < 1))" && ok "wall-band piece x=0.8, y in [2.8,3.0] - 1e-4: cert0 refused at $(grep -m1 '^UNCERT0' $ZT/k7_v40.log | cut -c9-80); exact mu(1/2, 3, 0) = $ex" || bad "v4 piece: $(tail -1 $ZT/k7_v40.log) $ex"
# a piece inside Q = [1,2] x [2.8,3.8] lowered by 1e-4: the germ (1.5 + th/2, 3.3 + 0.28 th, th) drops below 1 for
# th < ~0.03 deg (exact mu 0.99994 at th = 0.014 deg); refused at theta = 0 (c_x -> 1.5+) and at theta >= 0.014 deg
$T perturb $K7 $ZT/k7_v8.txt --op seg-delta --at 8,16,8,17 --dw -100000000 > /dev/null
$Z cert0 $ZT/k7_v8.txt --full --uncert-cap 3 > $ZT/k7_v80.log 2>&1
grep -q '^NOT VERIFIED' $ZT/k7_v80.log && ok "piece x=1.6, y in [3.2,3.4] - 1e-4: cert0 refused ($(grep '^worst' $ZT/k7_v80.log | cut -c1-90))" || bad "v8 piece accepted by cert0"
$Z cert $ZT/k7_v8.txt --full --threads $TH --umin 10 --xlo 15 --xhi 15 --ylo 32 --yhi 33 --bins 0-0 --uncert-cap 3 > $ZT/k7_v8.log 2>&1
ex=$($T mu $ZT/k7_v8.txt 12289/8192 2162733/655360 1/8192 | awk '{print $5}')
[ "$(verdict $ZT/k7_v8.log)" = "NOT VERIFIED" ] && python3 -c "import sys; sys.exit(not ($ex < 1))" && ok "  cert --umin 10 (theta >= 0.014 deg) near (1.5, 3.3): refused; exact mu at (1.5 + th/2, 3.3 + 0.28 th, th = 0.014 deg) = $ex" || bad "  v8 cert: $(verdict $ZT/k7_v8.log) $ex"

echo "== A6 soundness harness: exact mu at random admissible rational poses >= the box bound"
N=300; [ $QUICK = 1 ] && N=60
$T sound0 $K7 --n $N --poses 15 --seed 21 > $ZT/k7_sound0.txt; r=$(tail -1 $ZT/k7_sound0.txt)
echo "$r" | grep -q ' 0 FAIL' && ok "theta = 0 (Lemma Z0), k2m3 box: $r" || bad "sound0: $r"
$T sound0 $ZT/k7_v8.txt --n $((N/2)) --poses 15 --seed 22 > $ZT/k7v8_sound0.txt; r=$(tail -1 $ZT/k7v8_sound0.txt)
echo "$r" | grep -q ' 0 FAIL' && ok "theta = 0, weakened box: $r" || bad "sound0 weakened: $r"
$T sound $K7 --n $N --poses 15 --seed 23 --area > $ZT/k7_sound_a.txt; r=$(tail -1 $ZT/k7_sound_a.txt)
echo "$r" | grep -q ' 0 FAIL' && ok "boxes at the area square's sides/corners/inside (Lemmas 1-3, K), k2m3 box: $r" || bad "sound --area: $r"
$T sound $K7 --n $N --poses 15 --seed 24 --area --refl > $ZT/k7_sound_ar.txt; r=$(tail -1 $ZT/k7_sound_ar.txt)
echo "$r" | grep -q ' 0 FAIL' && ok "same, reflected cover: $r" || bad "sound --area --refl: $r"
$T sound $K7 --n $N --poses 15 --seed 25 > $ZT/k7_sound.txt; r=$(tail -1 $ZT/k7_sound.txt)
echo "$r" | grep -q ' 0 FAIL' && ok "generic random boxes (germs, walls), k2m3 box: $r" || bad "sound: $r"
$T sound $ZT/k7_rho.txt --n $((N/2)) --poses 15 --seed 26 --area > $ZT/k7rho_sound.txt; r=$(tail -1 $ZT/k7rho_sound.txt)
echo "$r" | grep -q ' 0 FAIL' && ok "weakened box (rho - 1e-6), area boxes: $r" || bad "sound rho: $r"
$Z cert $K7 --d4 --threads $TH --umin 10 --xlo 16 --xhi 19 --ylo 16 --yhi 30 --tight 0.0005 --dump-tight $ZT/k7_tight.txt > $ZT/k7_tightrun.log 2>&1
$T tight $ZT/k7_tight.txt $K7 --max $((N/2)) --poses 12 --seed 27 > $ZT/k7_tight.out; r=$(tail -1 $ZT/k7_tight.out)
echo "$r" | grep -q ' 0 FAIL' && ok "tightest certified leaves (bound < 1.0005) around the area square's left side: $r" || bad "tight: $r"

echo "== A8 --first-order (Lemma U, ZMX2_AREA.md sec 11): zero-limit germs at theta -> 0+"
# the wall band and the column c_x = 1.5 (single germs): certified down to theta = 0
$Z cert $K7 --d4 --first-order --threads $TH --xlo 4 --xhi 5 --ylo 26 --yhi 30 --bins 0-0 > $ZT/k7_fo_wall.log 2>&1
v=$(verdict $ZT/k7_fo_wall.log)
[ "$v" = "REGION CLEAN" ] && ok "wall band c_x in [0.4,0.6], c_y in [2.6,3.1], theta in [0, 14.25 deg]: clean ($(census $ZT/k7_fo_wall.log))" || bad "wall band first-order: $v"
$Z cert $K7 --d4 --first-order --threads $TH --xlo 15 --xhi 15 --ylo 26 --yhi 27 --bins 0-0 > $ZT/k7_fo_col.log 2>&1
v=$(verdict $ZT/k7_fo_col.log)
[ "$v" = "REGION CLEAN" ] && ok "column c_x in [1.5,1.6], c_y in [2.6,2.8]: clean ($(census $ZT/k7_fo_col.log))" || bad "column first-order: $v"
$Z cert $K7 --d4 --threads $TH --xlo 15 --xhi 15 --ylo 26 --yhi 27 --bins 0-0 --uncert-cap 3 > $ZT/k7_nofo_col.log 2>&1
[ "$(verdict $ZT/k7_nofo_col.log)" = "NOT VERIFIED" ] && ok "  the same cells without --first-order: refused (the obstruction of sec 11)" || bad "  column without first-order not refused"
# rejection: the wall-band first-order gain removed (h-lines y = 2.6..3.4 zeroed on x in [1, 1.2]): mu = 1 at
# theta = 0 on the wall, mu < 1 for small theta > 0 against the wall
python3 - $K7 $ZT/k7_gain.txt <<'EOF'
import sys
sys.path.insert(0, 'search')
import zmx2_tools as T
cv = T.load(sys.argv[1])
cv['segments'] = [(a, b, c, d, 0 if (b == d and b in (13, 14, 15, 16, 17) and {a, c} == {5, 6}) else w)
                  for (a, b, c, d, w) in cv['segments']]
T.write(sys.argv[2], cv, comment='h-lines y = 2.6..3.4 zeroed on x in [1, 1.2]')
EOF
$Z cert $ZT/k7_gain.txt --full --first-order --threads $TH --xlo 5 --xhi 5 --ylo 29 --yhi 30 --bins 0-0 --uncert-cap 3 > $ZT/k7_gain.log 2>&1
ex=$($T mu $ZT/k7_gain.txt 50001/100000 3 1/100000 | awk '{print $5}')
[ "$(verdict $ZT/k7_gain.log)" = "NOT VERIFIED" ] && python3 -c "import sys; sys.exit(not ($ex < 1))" && ok "wall gain removed: refused with --first-order; exact mu at (0.50001, 3, u = 1e-5) = $ex" || bad "gain removed: $(verdict $ZT/k7_gain.log) $ex"
$Z cert $ZT/k7_v8.txt --full --first-order --threads $TH --xlo 15 --xhi 15 --ylo 32 --yhi 33 --bins 0-0 --uncert-cap 3 > $ZT/k7_v8fo.log 2>&1
[ "$(verdict $ZT/k7_v8fo.log)" = "NOT VERIFIED" ] && ok "column piece x=1.6 lowered by 1e-4: refused with --first-order" || bad "v8 first-order: accepted"
# harness: exact mu at random poses with tiny angles (germ variables) of boxes near the germs
NH=$([ $QUICK = 1 ] && echo 40 || echo 150)
sd=80
for near in 0.5,3.0,0 0.5,2.9,0 1.5,3.0,0 1.5,3.1,0 1.5,2.3,0 1.5,1.5,0 0.5,0.5,0 2.48,2.75,0; do
  sd=$((sd+1))
  r=$($T sound $K7 --n $NH --poses 15 --seed $sd --near $near --zflags first-order --small-u | tail -1)
  echo "$r" | grep -q ' 0 FAIL' && ok "first-order harness near $near: $r" || bad "first-order harness near $near: $r"
done
r=$($T sound $ZT/k7_gain.txt --n $NH --poses 15 --seed 91 --near 0.5,3.0,0 --zflags first-order --small-u | tail -1)
echo "$r" | grep -q ' 0 FAIL' && ok "first-order harness, gain-removed cover, near the wall: $r" || bad "harness gain: $r"
r=$($T sound $K7 --n $N --poses 15 --seed 92 --area --zflags first-order --small-u | tail -1)
echo "$r" | grep -q ' 0 FAIL' && ok "first-order harness, area boxes: $r" || bad "harness area fo: $r"
r=$($T sound $K7 --n $N --poses 15 --seed 93 --refl --zflags first-order --small-u | tail -1)
echo "$r" | grep -q ' 0 FAIL' && ok "first-order harness, generic boxes, reflected cover: $r" || bad "harness refl fo: $r"

echo "== A9 double germs (Lemma V and the per-cell Lemma U, ZMX2_AREA.md sec 11.5-11.7)"
# the cells of sec 11.4 that the old first-order bound left: certified down to theta = 0
$Z cert $K7 --d4 --first-order --threads $TH --xlo 15 --xhi 15 --ylo 15 --yhi 34 --bins 0-0 > $ZT/k7_dg_col.log 2>&1
v=$(verdict $ZT/k7_dg_col.log)
[ "$v" = "REGION CLEAN" ] && ok "double-germ column c_x in [1.5,1.6], c_y in [1.5,3.5], theta in [0, 14.25 deg]: clean ($(census $ZT/k7_dg_col.log))" || bad "double-germ column: $v"
$Z cert $K7 --d4 --first-order --threads $TH --xlo 23 --xhi 23 --ylo 23 --yhi 23 --bins 0-0 > $ZT/k7_dg_23.log 2>&1
v=$(verdict $ZT/k7_dg_23.log)
[ "$v" = "REGION CLEAN" ] && ok "R-corner germ cell (2.3, 2.3): clean ($(census $ZT/k7_dg_23.log))" || bad "R-corner cell: $v"
$Z cert $K7 --full --first-order --threads $TH --xlo 15 --xhi 15 --ylo 54 --yhi 54 --bins 0-0 > $ZT/k7_dg_55.log 2>&1
v=$(verdict $ZT/k7_dg_55.log)
[ "$v" = "REGION CLEAN" ] && ok "  its D4 image (1.5, 5.5), both passes of --full: clean ($(census $ZT/k7_dg_55.log))" || bad "(1.5, 5.5) full: $v"
# harness: exact mu at tiny-angle poses of boxes at the double-germ scale (sides 1/(10 2^a), a = 12..18)
NF=$([ $QUICK = 1 ] && echo 60 || echo 300)
sd=400
for near in 1.5,1.5,0 1.5,2.3,0 1.5,2.5,0 1.5,2.9,0 1.5,3.1,0 1.5,3.5,0 2.3,2.3,0 1.5,5.5,0; do
  for refl in "" "--refl"; do
    sd=$((sd+1))
    r=$($T sound $K7 --n $NF --poses 15 --seed $sd --near $near --r 0.00003 --ru 0.00003 --fine --zflags first-order --small-u $refl | tail -1)
    echo "$r" | grep -q ' 0 FAIL' && ok "fine harness near $near $refl: $r" || bad "fine harness near $near $refl: $r"
  done
done
# rejection 1 (pair cut, Lemma V): mass moved so that mu >= 1 at theta = 0 and for theta >= theta_min near the germ,
# but the pair y = 2.6 / 3.6 cut in the middle misses the added mass: mu < 1 for 0 < theta < ~0.003 deg at (1.5, 3.1)
python3 - $K7 $ZT/k7_dg.txt <<'EOF'
import sys
sys.path.insert(0, 'search')
import zmx2_tools as T
cv = T.load(sys.argv[1])
A, dl = 20000000, 10000000   # masses 2e-5 and 1e-5 (W = 1e12)
segs = []
for (a, b, c, d, w) in cv['segments']:
    if b == d and b == 18 and {a, c} == {5, 6}: w += A        # y = 3.6 on x in [1.0, 1.2]
    if b == d and b == 13 and {a, c} == {8, 9}: w += A        # y = 2.6 on x in [1.6, 1.8]
    if a == c and a == 8 and {b, d} == {16, 17}: w -= dl      # x = 1.6 on y in [3.2, 3.4]
    segs.append((a, b, c, d, w))
cv['segments'] = segs
T.write(sys.argv[2], cv, comment='double-germ rejection: +2e-5 on y=3.6 x[1,1.2] and y=2.6 x[1.6,1.8], -1e-5 on x=1.6 y[3.2,3.4]')
EOF
$Z cert0 $ZT/k7_dg.txt --full > $ZT/k7_dg_cert0.log 2>&1
$Z cert $ZT/k7_dg.txt --full --umin 10 --threads $TH --xlo 14 --xhi 16 --ylo 29 --yhi 32 > $ZT/k7_dg_umin.log 2>&1
$Z cert $ZT/k7_dg.txt --full --first-order --threads $TH --xlo 15 --xhi 15 --ylo 30 --yhi 31 --bins 0-0 --uncert-cap 3 > $ZT/k7_dg_fo.log 2>&1
ex=$($T mu $ZT/k7_dg.txt 15000020011/10000000000 31000004/10000000 1/500000 | awk '{print $5}')
[ "$(verdict $ZT/k7_dg_cert0.log)" = "VERIFIED" ] && [ "$(verdict $ZT/k7_dg_umin.log)" = "REGION CLEAN" ] && ok "pair-cut rejection cover: theta = 0 VERIFIED, theta >= 0.014 deg clean around the germ" || bad "pair-cut cover: $(verdict $ZT/k7_dg_cert0.log) / $(verdict $ZT/k7_dg_umin.log)"
[ "$(verdict $ZT/k7_dg_fo.log)" = "NOT VERIFIED" ] && grep -q '^UNCERT .*x\[1.5000000,1.5000122\] y\[3.1000000' $ZT/k7_dg_fo.log && python3 -c "import sys; sys.exit(not ($ex < 1))" && ok "  refused by --first-order at (1.5, 3.1); exact mu at (1.5000020011, 3.1000004, u = 2e-6) = $ex" || bad "pair-cut cover first-order: $(verdict $ZT/k7_dg_fo.log) $ex"
r=$($T sound $ZT/k7_dg.txt --n $NF --poses 15 --seed 431 --near 1.5,3.1,0 --r 0.00003 --ru 0.00003 --fine --zflags first-order --small-u | tail -1)
echo "$r" | grep -q ' 0 FAIL' && ok "  fine harness on it near (1.5, 3.1): $r" || bad "  harness pair-cut cover: $r"
# rejection 2 (R corner, Lemma K (d)): y = 1.8 zeroed on x in [1.8, 2]: a square inside R touching its left side,
# rotated so that its lowest vertex pokes below y = 1.8 near x = 1.8, loses area with no line mass
python3 - $K7 $ZT/k7_rc.txt <<'EOF'
import sys
sys.path.insert(0, 'search')
import zmx2_tools as T
cv = T.load(sys.argv[1])
cv['segments'] = [(a, b, c, d, 0 if (b == d and b == 9 and {a, c} == {9, 10}) else w) for (a, b, c, d, w) in cv['segments']]
T.write(sys.argv[2], cv, comment='R-corner rejection: y = 1.8 zeroed on x in [1.8, 2.0]')
EOF
$Z cert0 $ZT/k7_rc.txt --full > $ZT/k7_rc_cert0.log 2>&1
$Z cert $ZT/k7_rc.txt --full --first-order --threads $TH --xlo 23 --xhi 23 --ylo 23 --yhi 23 --bins 0-0 --uncert-cap 3 > $ZT/k7_rc_fo.log 2>&1
ex=$($T mu $ZT/k7_rc.txt 23001/10000 230008/100000 1/10000 | awk '{print $5}')
[ "$(verdict $ZT/k7_rc_cert0.log)" = "VERIFIED" ] && [ "$(verdict $ZT/k7_rc_fo.log)" = "NOT VERIFIED" ] && python3 -c "import sys; sys.exit(not ($ex < 1))" && ok "R-corner rejection cover: theta = 0 VERIFIED, refused by --first-order at (2.3, 2.3); exact mu at (2.3001, 2.30008, u = 1e-4) = $ex" || bad "R-corner cover: $(verdict $ZT/k7_rc_cert0.log) / $(verdict $ZT/k7_rc_fo.log) $ex"
r=$($T sound $ZT/k7_rc.txt --n $NF --poses 15 --seed 432 --near 2.3,2.3,0 --r 0.00003 --ru 0.00003 --fine --zflags first-order --small-u | tail -1)
echo "$r" | grep -q ' 0 FAIL' && ok "  fine harness on it near (2.3, 2.3): $r" || bad "  harness R-corner cover: $r"

echo "== A7 covers without area densities: behaviour unchanged"
$Z cert $C21 --d4 --threads $TH --log $ZT/s21_d4.log > /dev/null 2>&1
python3 search/zmx2_census_cmp.py $ZT/s21_d4.log certificates/s21/zmx2_d4/roots.log > $ZT/s21_cmp.txt && ok "s(21) --d4 census = shipped zmx2_d4/roots.log root for root ($(cat $ZT/s21_cmp.txt | sed 's/.*: //'))" || bad "s(21) census differs: $(cat $ZT/s21_cmp.txt)"
[ "$(head -1 $ZT/s21_d4.log)" = "$(head -1 certificates/s21/zmx2_d4/roots.log)" ] && ok "  log header unchanged" || bad "  log header changed"

echo
echo "zmx2_area_tests: $npass passed, $nfail failed"
[ $nfail = 0 ]
