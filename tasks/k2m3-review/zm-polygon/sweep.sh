#!/usr/bin/env bash
# $1 core, $2 seed: random straddle roots, 120 s each, total budget 1500 s
core=$1; seed=$2; out=sweep_$core.out; : > $out; start=$(date +%s)
python3 -c "
import random; r=random.Random($seed)
for _ in range(200):
    i=r.randrange(11,25); j=r.randrange(11,35); k=r.randrange(0,8)
    print(f'{i}/10 {j}/10 {k}')" | while read x y k; do
  [ $(( $(date +%s) - start )) -gt 1500 ] && break
  timeout 120 taskset -c $core python3 test_leaves.py $x $y $k >> $out 2>&1 || echo "timeout/err $x $y $k" >> $out
done
echo DONE >> $out
