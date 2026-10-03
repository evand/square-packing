#!/usr/bin/env bash
# qx2_xloop.sh -- one exact-check pass of the qx2 "zmx2 as oracle" loop (2026-10-02).
# usage: qx2_xloop.sh RUNDIR K [CORES]
#   RUNDIR = runs/qx2_<tag> holding sol.npz (a converged qx2_lp / qx2_lp_cached run); K = box size (2R + 3).
# Rows: see the filter in the python block (theta range, admissibility, no images).
# Steps: project (exact, support fixed) -> axis (Lemma Z, exact) -> family + box cover -> germ scan (exact)
#        -> zmx2 d4 + cert0 + cert --d4 --first-order -> rows file of the uncertified boxes' float-min poses
#        (plus diagonal images) for the next warm LP.  Prints a one-line verdict per step.
set -euo pipefail
RD=${1:?rundir}; K=${2:?box size}; CORES=${3:-8-15}
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
NC=$(( $(echo "$CORES" | cut -d- -f2) - $(echo "$CORES" | cut -d- -f1) + 1 ))
Z=../verify2/target/release/zmx2
J=$RD/sol_exact.json; B=$RD/sol_exact_box$K.txt
taskset -c "$CORES" python3 qx2_exact.py project $RD/sol.npz > $RD/project.out 2>&1
grep -E "rank|exact:" $RD/project.out
taskset -c "$CORES" python3 qx2_exact.py axis $J | tail -1
python3 qx2_exact.py family $J | tail -1
python3 qx2_exact.py cover $J --k $K | tail -1
HI="$K/2"   # R + 3/2 = K/2 with R = (K-3)/2
taskset -c "$CORES" python3 qx2_germscan.py $B --hi $HI --nproc $NC > $RD/germscan.out 2>&1; grep "germ scan" $RD/germscan.out
$Z d4 $B
taskset -c "$CORES" $Z cert0 $B --d4 | tail -1
taskset -c "$CORES" $Z cert $B --d4 --first-order --threads $NC > $RD/zmx2_fo.log 2>&1 || true
tail -1 $RD/zmx2_fo.log; grep '^done' $RD/zmx2_fo.log
python3 - "$RD" <<'EOF'
import re, math, sys, numpy as np
rd = sys.argv[1]; P = []; nb = 0; fm = []
for l in open(rd + '/zmx2_fo.log'):
    if not l.startswith('UNCERT'): continue
    m = re.search(r'float-min ([\d.]+) at \(([\d.]+),([\d.]+),u ([\d.]+)', l)
    if not m: continue
    f, x, y, u = map(float, m.groups()); th = 2 * math.atan(u); fm.append(f)
    # zmx2 prints poses to ~1e-6: useless at near-tangent tiny-theta poses (covered by the exact germ / axis rows anyway).
    # Keep theta in [0.05, 89.95] deg and LP-admissible poses; no diagonal images (the quadrant model is diagonal-symmetric,
    # and theta = 90 deg is outside its domain [0, 90)).
    if not (math.radians(0.05) <= th <= math.radians(89.95)): continue
    hw = (math.cos(th) + math.sin(th)) / 2
    if x < hw or y < hw: continue
    P.append((x, y, th))
if P:
    P = np.unique(np.round(np.array(P), 12), axis=0); np.savez(rd + '/zmx2_rows.npz', rows=P)
    print(f"rows file: {len(P)} poses; uncert float-min: min {min(fm):.6f}, below 1: {sum(v < 1 for v in fm)}")
else:
    print("rows file: none (no uncertified boxes)")
EOF
