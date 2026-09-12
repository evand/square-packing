#!/bin/sh
# tasks/rank-diag round 4 -- finish-and-certify, then the E1 fine-lattice test.
#
#   sh search/rankfamily_e1.sh finish A|B|P      resume a stopped E2 run's checkpoint
#                                                (poses + coverage rows + clique rows + polygon
#                                                rows, the last rebuilt from their exact anchors)
#                                                with NO lattice pricing, and run to convergence
#                                                so the final measure certifies.
#   sh search/rankfamily_e1.sh fine   A|B|P      the same, plus ONE pricing stage on the FINE
#                                                lattice (pitch 0.02, dtheta 1.25 deg,
#                                                1,350,226 candidates), then converge and certify.
#
# The number E1 asks for is the RISE from the finish value to the fine value: ~0 means the
# 0.04 / 2.5 deg lattice was already resolving the continuum optimum, a rise of 0.05+ means it
# was not.
cd "$(dirname "$0")/.." || exit 1
CL="python3 search/cliquelever.py run"
PENT="--pent 5,7,9 --pent-restarts 2500 --pent-cands 400 --pent-want 10 --pent-time 300"
BASE="--threads 2 --procs 2 --lp-tlim 0 --master --lattice-every 0 --clique-time 60 \
      --cq-want 60 --cq-top 40 --cq-age 0 --row-pitch 0 --iters 100000 --ckpt 5"

case "$2" in
  A) SRC=E2A; BR="--corners 1111 --patterns 01010101 --chord" ;;
  B) SRC=E2B; BR="--corners 1111 --patterns ........ --chord" ;;
  P) SRC=E2P; BR="--corners .... --patterns ........" ;;
  *) echo "usage: $0 {finish|fine} {A|B|P}"; exit 1 ;;
esac

case "$1" in
  finish) TAG=${SRC}f;  FROM=$SRC;      EXTRA="--price 0" ;;
  fine)   TAG=${SRC}n;  FROM=${SRC}f;   EXTRA="--price 1 --price-pitch 0.02 --price-dth 1.25 \
                                               --cg-want 400" ;;
  *) echo "usage: $0 {finish|fine} {A|B|P}"; exit 1 ;;
esac

R=runs
setsid nohup $CL "$TAG" --exact $R/cl_${FROM}_poses.txt --resume-rows $R/cl_${FROM}_rows.txt \
    --resume-cliques $R/cl_${FROM}_cliques.txt --resume-pgons $R/cl_${FROM}_pgons.json \
    $BR $BASE $PENT $EXTRA --time "${T:-3600}" > "$R/$TAG.out" 2>&1 < /dev/null &
echo "$TAG (from $FROM) pid $!"
