#!/usr/bin/env bash
# make_covers.sh -- the golf variants of the s(21) cover (unscaled transforms; scaling is a separate step).
# Run from s12/.  Outputs in search/golf/covers/ (gitignored? no: small enough, but regenerable; see GOLF_PILOT.md).
set -eu
cd "$(dirname "$0")/../.."
G=search/golf/golf.py
C=certificates/s21/s21_mixed_cover_5.txt
O=search/golf/covers
mkdir -p $O
T="taskset -c ${CORES:-11}"
$T python3 $G coarsen $C $O/v_seg2.txt --k 2
$T python3 $G coarsen $C $O/v_seg5.txt --k 5
$T python3 $G merge $C $O/v_merge05.txt --h 1/20
$T python3 $G merge $C $O/v_merge10.txt --h 1/10
$T python3 $G merge $C $O/v_merge10w.txt --h 1/10 --wall 3/10
$T python3 $G drop $C $O/v_drop1e4.txt --min-w 1/10000
$T python3 $G drop $C $O/v_drop1e3.txt --min-w 1/1000
