#!/bin/sh
# stage 2: polish the stage-1 results with fill-up to the maximal point count and kick moves
cd "$(dirname "$0")/../.." || exit 1
L=search/uniform/polish_launch.sh
T=${T:-1800}
# D4-symmetric, m = 12k-3
EXTRA="--mmax 21" sh $L $T 2 2 runs/pol_uni_k2_s3.65_f0.05_s1.txt
EXTRA="--mmax 21" sh $L $T 2 3 runs/pol_uni_k2_s3.65_f0.05_s1.txt
EXTRA="--mmax 33" sh $L $T 2 2 runs/pol_uni_k3_s3.75_f0.05_s1.txt
EXTRA="--mmax 33" sh $L $T 2 3 runs/pol_uni_k3_s3.75_f0.05_s1.txt
EXTRA="--mmax 45" sh $L $T 2 2 runs/pol_uni_k4_s3.70_f0.05_s1.txt
EXTRA="--mmax 45" sh $L $T 2 3 runs/pol_uni_k4_s3.70_f0.05_s1.txt
EXTRA="--mmax 57" sh $L $T 2 3 runs/pol_s12_56points_3.8_s1.txt
EXTRA="--mmax 57" sh $L $T 2 4 runs/pol_s12_56points_3.8_s1.txt
EXTRA="--mmax 57" sh $L $T 2 3 runs/pol_uni_t5b_s2.txt
EXTRA="--mmax 69" sh $L $T 2 2 runs/pol_uni_k6_s3.75_f0.05_s1.txt
EXTRA="--mmax 69" sh $L $T 2 3 runs/pol_uni_k6_s3.75_f0.05_s1.txt
EXTRA="--mmax 81" sh $L $T 2 2 runs/pol_uni_k7_s3.85_f0.05_s1.txt
EXTRA="--mmax 81" sh $L $T 2 3 runs/pol_uni_k7_s3.85_f0.05_s1.txt
EXTRA="--mmax 93" sh $L $T 2 2 runs/pol_k8fill_test.txt
EXTRA="--mmax 93" sh $L $T 2 3 runs/pol_k8fill_test.txt
# no symmetry, m = 12k-1
EXTRA="--sym none --mmax 11" sh $L $T 2 2 runs/pol_k1snap370_s1.txt
EXTRA="--sym none --mmax 11" sh $L $T 2 3 runs/pol_k1snap370_s1.txt
EXTRA="--sym none --mmax 11" sh $L $T 2 2 runs/pol_k1snap375_s1.txt
EXTRA="--sym none --mmax 23" sh $L $T 2 5 runs/pol_uni_k2_s3.65_f0.05_s1.txt
EXTRA="--sym none --mmax 35" sh $L $T 2 5 runs/pol_uni_k3_s3.75_f0.05_s1.txt
EXTRA="--sym none --mmax 59" sh $L $T 2 5 runs/pol_s12_56points_3.8_s1.txt
EXTRA="--sym none --mmax 83" sh $L $T 2 5 runs/pol_uni_k7_s3.85_f0.05_s1.txt
EXTRA="--sym none --mmax 95" sh $L $T 2 5 runs/pol_k8fill_test.txt
