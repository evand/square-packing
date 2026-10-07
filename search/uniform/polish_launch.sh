#!/bin/sh
# polish_launch.sh TIME THREADS SEED CERT [CERT ...]  -> runs/pol_<stem>_s<seed>.txt, log runs/plog_<stem>_s<seed>.txt
cd "$(dirname "$0")/../.." || exit 1
T=$1; TH=$2; SEED=$3; shift 3
for c in "$@"; do
  stem=$(basename "$c" .txt)
  nohup python3 search/uniform/polish.py "$c" "runs/pol_${stem}_s${SEED}.txt" --time "$T" --threads "$TH" --seed "$SEED" $EXTRA > "runs/plog_${stem}_s${SEED}.txt" 2>&1 &
done
echo launched
