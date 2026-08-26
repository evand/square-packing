#!/bin/sh
# stage 3: transfer structures between k (change the weight denominator k in the file, then
# fill or drop points to 12k-3 / 12k-1) and re-polish the best of each k with new seeds.
cd "$(dirname "$0")/../.." || exit 1
L=search/uniform/polish_launch.sh
T=${T:-1500}
B7=runs/pol_pol_uni_k7_s3.85_f0.05_s1_s3.txt     # k=7  m=81  3.8879
B5=runs/pol_pol_s12_56points_3.8_s1_s3.txt        # k=5  m=57  3.8520
B3=runs/pol_pol_uni_k3_s3.75_f0.05_s1_s5.txt      # k=3  m=35  3.8081 (none)
B3d=runs/pol_uni_k3_s3.75_f0.05_s1.txt            # k=3  m=33  3.8081 (d4)
B2=runs/pol_pol_uni_k2_s3.65_f0.05_s1_s3.txt      # k=2  m=21
B1=runs/pol_pol_k1snap370_s1_s3.txt               # k=1  m=11
mk() { sed "3s/.*/$3/" "$1" > "$2"; }   # copy with the weight denominator replaced by k
mk $B7 runs/from7_k8.txt 8
mk $B7 runs/from7_k6.txt 6
mk $B7 runs/from7_k5.txt 5
mk $B5 runs/from5_k4.txt 4
mk $B5 runs/from5_k6.txt 6
mk $B3d runs/from3_k4.txt 4
mk $B5 runs/from5_k3.txt 3
mk $B3d runs/from3_k2.txt 2
mk $B2 runs/from2_k1.txt 1
EXTRA="--mmax 81" sh $L $T 2 11 $B7
EXTRA="--mmax 81" sh $L $T 2 12 $B7
EXTRA="--mmax 93" sh $L $T 2 11 runs/from7_k8.txt
EXTRA="--mmax 93" sh $L $T 2 12 runs/from7_k8.txt
EXTRA="--mmax 69" sh $L $T 2 11 runs/from7_k6.txt
EXTRA="--mmax 69" sh $L $T 2 12 runs/from7_k6.txt
EXTRA="--mmax 69" sh $L $T 2 11 runs/from5_k6.txt
EXTRA="--mmax 57" sh $L $T 2 11 runs/from7_k5.txt
EXTRA="--mmax 57" sh $L $T 2 11 $B5
EXTRA="--mmax 57" sh $L $T 2 12 $B5
EXTRA="--mmax 45" sh $L $T 2 11 runs/from5_k4.txt
EXTRA="--mmax 45" sh $L $T 2 12 runs/from5_k4.txt
EXTRA="--mmax 45" sh $L $T 2 11 runs/from3_k4.txt
EXTRA="--mmax 33" sh $L $T 2 11 runs/from5_k3.txt
EXTRA="--mmax 33" sh $L $T 2 11 $B3d
EXTRA="--mmax 21" sh $L $T 2 11 runs/from3_k2.txt
EXTRA="--mmax 21" sh $L $T 2 11 $B2
EXTRA="--sym none --mmax 11" sh $L $T 2 11 $B1
EXTRA="--sym none --mmax 11" sh $L $T 2 12 $B1
EXTRA="--sym none --mmax 11" sh $L $T 2 11 runs/from2_k1.txt
