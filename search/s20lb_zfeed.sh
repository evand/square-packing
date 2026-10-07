#!/usr/bin/env bash
# s20lb_zfeed.sh TAG [LAST_DONE_ROUND] -- for each new round cover runs/lc_TAG_r{N}.txt of a running line_cover.py loop: zmx2 --d4 --sym-atoms
# on the cover scaled by 1.005 and 1.012, rows from the UNCERT boxes (s20lb_zrows.py), handed to the loop via
# runs/lc_TAG_inject.txt (read once per phase-B round).  search/S20_LB.md.  Stops when the loop's log says DONE/limit.
set -u
cd "$(dirname "$0")/.."
T=$1; Z=verify2/target/release/zmx2; D=runs/s20lb_zmx/feed_$T; mkdir -p $D
done_r=${2:--1}
while true; do
  grep -qE "DONE|round limit|time limit|STABLE|LP failed" runs/lc_$T.log 2>/dev/null && break
  last=$(ls runs/lc_${T}_r*.txt 2>/dev/null | sed 's/.*_r\([0-9]*\)\.txt/\1/' | sort -n | tail -1)
  if [[ -n "$last" && $last -gt $done_r ]]; then
    sleep 5; C=runs/lc_${T}_r$last.txt; logs=""
    for F in 1.005 1.012; do
      O=$D/r${last}_x$F.txt; L=$D/r${last}_x$F.log; rm -f $L
      python3 search/line_cover.py scale $C $O --factor $F > /dev/null
      R=$($Z cert $O --d4 --sym-atoms --threads ${ZTHREADS:-6} --uncert-cap 10 --log $L 2>&1 | grep -E "^(VERIFIED|NOT VERIFIED)" | head -1)
      echo "$(date +%T) r$last x$F: $R"; logs="$logs $L"; rm -f $O
    done
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 python3 search/s20lb_zrows.py $C $D/rows_r$last.txt $logs --thr 1.0 --n 400
    cat $D/rows_r$last.txt >> runs/lc_${T}_inject.txt
    done_r=$last
  fi
  sleep 30
done
echo "$(date +%T) loop finished; feed stops"
