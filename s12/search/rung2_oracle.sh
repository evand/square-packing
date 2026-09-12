#!/usr/bin/env bash
# rung2_oracle.sh -- one turn of the separation loop of search/RUNG2.md sec 10:
#   scale the cover LP's current weights up by P/Q, run the EXACT checker on the result at a
#   shallow depth, and write its uncertified-box poses as new hard rows for the LP.
#
# Where the cover is valid the checker reports 0 uncertified at depth 12 (RUNG2.md sec 4.3), so a
# box still uncertified at depth 8-10 is either a real violation or a near-violation -- in both
# cases exactly the row the LP is missing.  The oracle file is in the `cx cy theta_rad` format that
# `closed4.py --seed-rows` reads.
#
# usage: search/rung2_oracle.sh CHECKPOINT TAG [P Q] [DEPTH] [NPROC] [CORES]
set -u
CK=${1:?checkpoint (e.g. runs/closed4_r5_last.txt)}
TAG=${2:?tag}
P=${3:-103}; Q=${4:-100}
DEPTH=${5:-8}; NP=${6:-8}; CORES=${7:-0-7}
SC="runs/${TAG}_scaled.txt"
python3 search/scale_cover.py "$CK" "$P" "$Q" "$SC" || exit 1
taskset -c "$CORES" python3 search/zeromargin.py cert "$SC" --depth "$DEPTH" --nproc "$NP" \
        --disj --chain-from 0 --oracle "runs/${TAG}_oracle.txt" > "runs/${TAG}_check.log" 2>&1
grep -E "^done in|leaves:|VERIFIED|oracle rows" "runs/${TAG}_check.log"
