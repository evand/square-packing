#!/bin/sh
# stage 4 (spare time): try to lift k = 6, 8 with the k = 7 ring structure as the seed.
cd "$(dirname "$0")/../.." || exit 1
B7=certificates/s12_uniform_7of81_3.888.txt
I=search/uniform/ilp_uniform.py
for k in 6 8; do
  for s in 3.85 3.82; do
    mm=$((12 * k - 1))
    nohup python3 $I $k $s --fine 0.0125 --near $B7 --radius 0.06 --feas --mmax $mm --tl 600 --rounds 40 --threads 4 --tag near7_k${k}_s${s} > runs/log_near7_k${k}_s${s}.txt 2>&1 &
  done
done
sed "3s/.*/8/" $B7 > runs/from7b_k8.txt
sed "3s/.*/6/" $B7 > runs/from7b_k6.txt
EXTRA="--mmax 93 --threads 8" sh search/uniform/polish_launch.sh 1200 8 21 runs/from7b_k8.txt
EXTRA="--mmax 69 --threads 8" sh search/uniform/polish_launch.sh 1200 8 21 runs/from7b_k6.txt
