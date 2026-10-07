#!/usr/bin/env bash
# zm_sample.sh TAG COVER -- zm_mixed.py --d4 --cert-mode (the shipped s(21) settings) on a sample of 512 roots:
# the tile-germ block [1.4,1.6]^2 (256 roots), the two hottest tilted cells (theta ~ 30-50 deg) [2.4,2.5]x[1.2,1.3]
# and [1.35,1.45]x[1.2,1.3], the tilted-dip cells [0.5,0.6]x[1.4,1.5] (9.1 deg) and [2.3,2.4]x[1.6,1.7] (31.9 deg),
# 64 roots each.  Cores 6-11, 6 processes.  Records: search/golf/zm/TAG_REGION.jsonl.
set -u
cd "$(dirname "$0")/../.."
TAG=$1; C=$2
O=search/golf/zm; mkdir -p $O
run() { # name xlo xhi ylo yhi
  taskset -c ${ZC:-6-11} python3 search/zm_mixed.py cert $C --d4 --cert-mode --disj --chain-from 0 --depth 24 --pitch 1/20 \
      --ubins 16 --nproc ${NP:-6} --progress 100000 --cx-lo $2 --cx-hi $3 --cy-lo $4 --cy-hi $5 \
      --resume $O/${TAG}_$1.jsonl > $O/${TAG}_$1.out 2>&1
  echo "$1: $(grep -E '^(done in|  leaves|VERIFIED|NOT VERIFIED)' $O/${TAG}_$1.out | tr '\n' ' ')"
}
run germ 7/5 8/5 7/5 8/5
run tilt2412 12/5 5/2 6/5 13/10
run tilt1412 27/20 29/20 6/5 13/10
run dip9 1/2 3/5 7/5 3/2
run dip32 23/10 12/5 8/5 17/10
