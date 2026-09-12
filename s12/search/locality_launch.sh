#!/bin/bash
# locality_launch.sh -- the LOC_D table.  One detached run per (pose set, D, shape), at most
# --jobs at a time, 2 LP threads each (tasks/locality-ceiling budget).
#
#   bash search/locality_launch.sh A   # the leaf 01010101   (cl_A0101L2_poses.txt)
#   bash search/locality_launch.sh B   # the corner k = 4    (cl_B40KL2_poses.txt)
set -u
cd "$(dirname "$0")/.."
SET=${1:-A}
TIME=${2:-2400}
SHAPES=${3:-box}
DS=${4:-"1.5 2 2.5 3 3.5"}
case "$SET" in
  A) POSES=runs/cl_A0101L2_poses.txt; PINS="--corners 1111 --patterns 01010101 --chord" ;;
  B) POSES=runs/cl_B40KL2_poses.txt; PINS="--corners 1111 --patterns ........ --chord" ;;
  *) echo "unknown pose set $SET"; exit 1 ;;
esac
for SH in $SHAPES; do
  for D in $DS; do
    TAG="${SET}$(echo "$D" | tr -d .)${SH:0:1}"
    rm -f "runs/loc_${TAG}.log"
    setsid nohup python3 search/locality.py run "$TAG" --poses "$POSES" $PINS \
        --D "$D" --shape "$SH" --threads 2 --procs 2 --iters 400 --time "$TIME" \
        > "runs/loc_${TAG}.out" 2>&1 &
    echo "launched $TAG (D=$D $SH) pid $!"
    sleep 2
  done
done
