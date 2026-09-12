#!/bin/sh
# tasks/rank-diag round 3 -- the E2 measurement: QSTAB + the odd-polygon rank family on the
# restricted master WITH continuous lattice pricing (search/CLMASTER.md 2.3), so the pose set
# keeps growing while the rank rows keep cutting.  Three runs, 2 LP threads each.
#
#   E2A  the leaf 01010101              (resume cl_A0101L2_*)
#   E2B  the corner leaf k = 4          (resume cl_B40KL2_*)
#   E2P  NO BRANCH AT ALL               (tl_L2PURE_*, no corners, no patterns, no chord)
#
# E2P is the one that matters: the clique-only twin of it (`PUREM`, runs/launch_master.sh) is
# climbing toward 12 under lattice pricing.  If pure + polygons stays clearly below 12, the
# no-tree proof shape is back on the table.
#
# The lattice pricer charges polygon-row duals (rankfamily.pgon_cost) and every injection
# re-derives every polygon row over the enlarged pose set (rankfamily.regrow), so the rows stay
# maximal and the reduced costs stay honest.
cd "$(dirname "$0")/.." || exit 1
CL="python3 search/cliquelever.py run"
COMMON="--threads 2 --procs 2 --lp-tlim 0 --master --lattice-every 5 --cg-want 400 \
        --clique-time 60 --cq-want 60 --cq-top 40 --cq-age 0 --iters 100000 \
        --time ${T:-25200} --ckpt 5"
# --pent-restarts high enough that the separator is not the limit (300 per k converges the corner
# support at 11.519638, 2500 at 11.470839; see RANKDIAG.md 10.1)
PENT="--pent 5,7,9 --pent-restarts 2500 --pent-cands 400 --pent-want 10 --pent-time 300"
launch() { tag=$1; shift; setsid nohup $CL "$tag" "$@" $COMMON $PENT > "runs/$tag.out" 2>&1 < /dev/null & echo "$tag pid $!"; }

launch E2A --exact runs/cl_A0101L2_poses.txt --resume-rows runs/cl_A0101L2_rows.txt \
           --resume-cliques runs/cl_A0101L2_cliques.txt --corners 1111 --patterns 01010101 \
           --chord --row-pitch 0
launch E2B --exact runs/cl_B40KL2_poses.txt --resume-rows runs/cl_B40KL2_rows.txt \
           --resume-cliques runs/cl_B40KL2_cliques.txt --corners 1111 --patterns ........ \
           --chord --row-pitch 0
launch E2P --load-poses runs/tl_L2PURE_poses.txt --load-rows runs/tl_L2PURE_dual.txt \
           --corners .... --patterns ........
