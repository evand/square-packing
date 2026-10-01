#!/usr/bin/env bash
# zmx2_curve.sh -- zmx2 --d4 box count (the cheap Lean-cost proxy, GOLF_PILOT.md sec. 3) of the s(21) cover scaled to
# several totals.  Cores 6-11.  Output: search/golf/data/zmx2_curve.txt
cd "$(dirname "$0")/../.."
C=certificates/s21/s21_mixed_cover_5.txt
O=search/golf/zbis; mkdir -p $O
out=search/golf/data/zmx2_curve.txt; : > $out
for T in 20.80 20.82 20.85 20.8947 20.92 20.95 20.97 20.99 20.999; do
  python3 search/golf/golf.py scale $C $O/curve_$T.txt --total $T > /dev/null
  R=$(taskset -c 6-11 verify2/target/release/zmx2 cert $O/curve_$T.txt --d4 --threads 6 --uncert-cap 50 2>&1 | grep -E "^(done in|VERIFIED|NOT VERIFIED)" | tr '\n' ' ')
  echo "total $T : $R" | tee -a $out
  rm -f $O/curve_$T.txt
done
for F in search/golf/covers/cand_A_scale.txt search/golf/covers/cand_B_drop.txt; do
  R=$(taskset -c 6-11 verify2/target/release/zmx2 cert $F --d4 --threads 6 2>&1 | grep -E "^(done in|VERIFIED|NOT VERIFIED)" | tr '\n' ' ')
  echo "$F : $R" | tee -a $out
done
