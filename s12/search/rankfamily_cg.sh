#!/bin/sh
# tasks/rank-diag round 5 -- drive the column generation on the 0.04 / 2.5 deg lattice until the
# best lattice reduced cost is small, so the reported number is CG-converged and not just
# "wherever the clock stopped".
#
#   sh search/rankfamily_cg.sh P|B        resume E2P / E2B from its checkpoint (poses + coverage
#                                         rows + clique rows + polygon rows, the last rebuilt from
#                                         their exact anchors) with --lattice-every 5 still on.
#
# Watch `gap +X` in the lattice lines: that field IS rc.max() over the priced lattice.  Stop when
# it is below the target, then certify by the support route (search/rankfamily_e1.sh style: take
# the run's support measure, converge QSTAB + polygons on it, certify).
cd "$(dirname "$0")/.." || exit 1
CL="python3 search/cliquelever.py run"
PENT="--pent 5,7,9 --pent-restarts 2500 --pent-cands 400 --pent-want 10 --pent-time 300"
BASE="--threads 2 --procs 2 --lp-tlim 0 --master --lattice-every 5 --cg-want 400 \
      --clique-time 60 --cq-want 60 --cq-top 40 --cq-age 0 --row-pitch 0 --iters 100000 --ckpt 5"

case "$1" in
  P) SRC=E2P; TAG=E2Pg; BR="--corners .... --patterns ........" ;;
  B) SRC=E2B; TAG=E2Bg; BR="--corners 1111 --patterns ........ --chord" ;;
  *) echo "usage: $0 {P|B}"; exit 1 ;;
esac

R=runs
setsid nohup $CL "$TAG" --exact $R/cl_${SRC}_poses.txt --resume-rows $R/cl_${SRC}_rows.txt \
    --resume-cliques $R/cl_${SRC}_cliques.txt --resume-pgons $R/cl_${SRC}_pgons.json \
    $BR $BASE $PENT --price 0 --time "${T:-21600}" > "$R/$TAG.out" 2>&1 < /dev/null &
echo "$TAG (from $SRC) pid $!"
