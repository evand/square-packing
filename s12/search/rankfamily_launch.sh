#!/bin/sh
# tasks/rank-diag round 2 -- every run quoted in search/RANKDIAG.md sections 10-11.
# Usage:  sh search/rankfamily_launch.sh {support|price|loop|all}
# Inputs are the read-only checkpoints of the task brief, copied into runs/ (gitignored):
#   /home/evand/math/square-packing/s12/runs/inputs-2026-09-11/
# 2 LP threads per run; the loop runs are detached and take hours, the rest minutes.
cd "$(dirname "$0")/.." || exit 1
CL="python3 search/cliquelever.py run"
launch() { tag=$1; shift; setsid nohup $CL "$tag" "$@" > "runs/$tag.out" 2>&1 < /dev/null & echo "$tag pid $!"; }

# the rank family: odd polygons of 5, 7 and 9 anchors.  --pent-restarts is the setting that
# matters (300 per k converges the corner support at 11.519638, 2500 at 11.470839).
PENT="--pent 5,7,9 --pent-restarts 2500 --pent-cands 400 --pent-want 10 --pent-time 300"

# ---- the SUPPORT experiment: the same LP on the pose set the converged measure itself uses
SUP="--chord --row-pitch 0 --threads 2 --procs 2 --lp-tlim 10 --clique-time 30 --cq-want 40 \
     --cq-top 30 --cq-age 0 --iters 400 --time ${T:-3000}"
support() {
  launch SUPA0 --exact runs/cl_A0101L2_exact.txt --corners 1111 --patterns 01010101 $SUP
  launch SUPA  --exact runs/cl_A0101L2_exact.txt --corners 1111 --patterns 01010101 $SUP $PENT
  launch SUPB0 --exact runs/cl_B40KL2_exact.txt --corners 1111 --patterns ........ $SUP
  launch SUPB  --exact runs/cl_B40KL2_exact.txt --corners 1111 --patterns ........ $SUP $PENT
}

# ---- one PRICING stage on top of the leaf support's converged 11.000000
price() {
  launch SUPAP --exact runs/cl_A0101L2_exact.txt --corners 1111 --patterns 01010101 \
        $SUP $PENT --price 1
}

# ---- the LOOP experiment: the full recorded pose sets, resumed as runs/launch_2026-09-11.sh
#      does, on the restricted master (search/CLMASTER.md 4)
LOOP="--threads 2 --procs 2 --lp-tlim 0 --master --clique-time 60 --cq-want 60 --cq-top 40 \
      --cq-age 0 --row-pitch 0 --chord --time ${T:-9000}"
loop() {
  launch PGA --exact runs/cl_A0101L2_poses.txt --resume-rows runs/cl_A0101L2_rows.txt \
        --resume-cliques runs/cl_A0101L2_cliques.txt --corners 1111 --patterns 01010101 \
        $LOOP $PENT --price 1
  launch PGB --exact runs/cl_B40KL2_poses.txt --resume-rows runs/cl_B40KL2_rows.txt \
        --resume-cliques runs/cl_B40KL2_cliques.txt --corners 1111 --patterns ........ \
        $LOOP $PENT --price 0
}

case "${1:-all}" in
  support) support ;;
  price)   price ;;
  loop)    loop ;;
  all)     support; price; loop ;;
  *) echo "usage: $0 {support|price|loop|all}"; exit 1 ;;
esac
