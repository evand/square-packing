#!/bin/sh
# phase 2 of the t = 4 screen: one process per table row, 1 thread each (<= 8 cores in total),
# each in chunks of 600 s so no single process outlives the 15-minute cap.
M=/home/evand/math/square-packing/s12/runs
W="--warm $M/dual_PC1_support.txt --warm $M/cqx_PURE99_support.txt"
C="--cliques --cq-wall --cq-interior --cq-max 800 --cq-age 4 --row-age 4 --threads 1 --rowloops 10 --stages 4"
N=${N:-3}
S=${S:-600}

for P in 01010101 01010110 01011010 01100110; do
  nohup sh search/t4screen_loop.sh 4.0 L$P $N $S --corners 1111 --patterns $P --price-pattern $P \
      $W $C > runs/L$P.out 2>&1 &
done
nohup sh search/t4screen_loop.sh 4.0 K40 $N $S --corners 1111 --patterns "........" $W $C > runs/K40.out 2>&1 &
nohup sh search/t4screen_loop.sh 4.0 U40 $N $S --corners .... --patterns "........" $W $C > runs/U40.out 2>&1 &
nohup sh search/t4screen_loop.sh 3.98 V98 $N $S --corners 1111 --patterns "........" \
    --warm $M/branch_J16i_dual.txt --cq-load $M/branch_J16i_cliques.txt \
    --cliques --cq-interior --cq-max 800 --cq-age 4 --row-age 4 --threads 1 --rowloops 10 --stages 4 > runs/V98.out 2>&1 &
echo "launched 7"
