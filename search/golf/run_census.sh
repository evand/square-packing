#!/usr/bin/env bash
# run_census.sh TAG COVER -- Lean-tree census of the stratified sample for one cover (cores 6-11, 6 processes):
# the 18 sampled zm roots of the heavy cells (timeout $RT s each), then the 51 sampled light cells.
set -u
cd "$(dirname "$0")/../.."
TAG=$1; C=$2
D=search/golf/data
taskset -c 6-11 python3 search/golf/lean_census.py $C --n 21 --roots-file $D/sample_roots.txt --nproc 6 \
    --timeout ${RT:-1200} --out $D/lc_${TAG}_roots.jsonl 2>> $D/lc_${TAG}.log
taskset -c 6-11 python3 search/golf/lean_census.py $C --n 21 --cells-file $D/sample_cells_light.txt --nproc 6 \
    --timeout ${CT:-1200} --out $D/lc_${TAG}_cells.jsonl 2>> $D/lc_${TAG}.log
echo "DONE $TAG" >> $D/lc_${TAG}.log
