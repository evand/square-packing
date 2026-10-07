#!/usr/bin/env bash
# run_audit.sh -- reproduce the zmx2 audit tests (search/ZMX2_AUDIT.md).  Run from the repo root:
#   CORES=10,11 bash search/zmx2_audit/run_audit.sh
# Needs verify2/target/release/zmx2 (or ZMX2=...).  Scratch files in $ZT (default runs/zmx2_audit).
set -u
cd "$(dirname "$0")/../.."
CORES=${CORES:-10,11}; export ZT=${ZT:-runs/zmx2_audit}; mkdir -p $ZT
Z=${ZMX2:-verify2/target/release/zmx2}; export ZMX2=$Z
A=search/zmx2_audit; C=runs/line-cover_m5_candidate_x1003.txt
run() { taskset -c $CORES "$@"; }
python3 $A/mkcovers.py $C
echo "== F1 integer wrap: s_num = 5 + 2^125 (expected today: VERIFIED-D4 for s ~ 4e37, i.e. the bug)"
$Z info $ZT/wrap_s.txt | head -1
run $Z cert $ZT/wrap_s.txt --d4 --threads 2 | tail -1
echo "== F1 integer wrap: weights 2^127-1, 2^127-1, 2 at (0,0) (expected today: total printed unchanged)"
$Z info $ZT/wrap_w.txt | tail -1
echo "== A1 candidate x 0.9943 (--d4): must be refused at (0.5736, 1.4406, 9.21 deg)"
run $Z cert $ZT/x9943.txt --d4 --threads 2 --uncert-cap 50 > $ZT/x9943.log; tail -1 $ZT/x9943.log
grep UNCERT $ZT/x9943.log | awk -F'float-min ' '{print $2}' | sort -n | head -1
echo "== A2 germ hole g1 (--full, c in [2.4,2.6]x[2.3,2.7], bin 0): must be refused at (2.5, 2.5, theta -> 0+)"
run $Z cert $ZT/g1.txt --full --threads 2 --xlo 24 --xhi 25 --ylo 23 --yhi 26 --bins 0-0 --uncert-cap 200 > $ZT/g1.log; tail -1 $ZT/g1.log
grep UNCERT $ZT/g1.log | awk -F'float-min ' '{print $2}' | sort -n | head -1
echo "== A3 D4-symmetric germ hole g2a (--d4): must be refused"
run $Z cert $ZT/g2a.txt --d4 --threads 2 --uncert-cap 100 | tail -1
echo "== A4 wall-germ hole w1 (--full, c in [0.5,0.6]x[2.4,2.6], bin 0): must be refused at (0.5, 2.5, theta -> 0+)"
run $Z cert $ZT/w1.txt --full --threads 2 --xlo 5 --xhi 5 --ylo 24 --yhi 25 --bins 0-0 --uncert-cap 100000 > $ZT/w1.log; tail -1 $ZT/w1.log
grep UNCERT $ZT/w1.log | awk -F'float-min ' '{print $2}' | sort -n | head -1
echo "== A5 nested-box probes at the violating germ / wall / dip poses (bound <= exact mu)"
(cd $A && run python3 probe.py ../../$ZT/g1.txt 250000084/100000000,5/2,1/1000000 5/2,5/2,0)
(cd $A && run python3 probe.py ../../$ZT/g2a.txt 249999914/100000000,5/2,1/1000000 5/2,5/2,1/1000000 5/2,5/2,0 5/2,249999/100000,1/100000000)
(cd $A && run python3 probe.py ../../$ZT/w1.txt 1000001/2000000,5/2,1/2000000 1/2,5/2,0 1/2,501/200,0)
(cd $A && run python3 probe.py ../../$ZT/x9943.txt 573593/1000000,1440624/1000000,805601/10000000)
echo "== A6 random adversarial covers, nested-box probes (seeds 10-29: 800 covers x 60 poses)"
for sd in $(seq 10 29); do run python3 $A/campaign.py $sd 40 60 | grep -E 'FAIL:|campaign'; done
echo "== A7 differential cert: random covers scaled to exact mu = 1 - 1e-5 at a known pose; must be refused"
run python3 $A/diffcert.py 2 150 1/100000 | tail -1
run python3 $A/diffcert.py 3 40 1/100000 | tail -1
echo "== A8 tightest certified leaves of the candidate (--full, bound < 1.001), independent exact mu"
run $Z cert $C --full --threads 2 --tight 0.001 --dump-tight $ZT/tight_full.txt > /dev/null
run python3 $A/tightchk.py $C $ZT/tight_full.txt 300
