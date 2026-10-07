#!/bin/sh
# install and check the selected uniform certificates (see UNIFORM.md)
cd "$(dirname "$0")/../.." || exit 1
F=search/uniform/finalize.sh
sh $F runs/pol_pol_pol_k1snap370_s1_s3_s11.txt      s12_uniform_1of11_3.775
sh $F runs/pol_pol_uni_k2_s3.65_f0.05_s1_s5.txt     s12_uniform_2of23_3.767
sh $F runs/pol_pol_uni_k3_s3.75_f0.05_s1_s2.txt     s12_uniform_3of33_3.808
sh $F runs/pol_pol_uni_k4_s3.70_f0.05_s1_s3.txt     s12_uniform_4of45_3.793
sh $F runs/pol_s12_56points_3.8_s1.txt              s12_uniform_5of56_3.852
sh $F runs/pol_pol_s12_56points_3.8_s1_s3.txt       s12_uniform_5of57_3.852
sh $F runs/pol_pol_uni_k6_s3.75_f0.05_s1_s2.txt     s12_uniform_6of69_3.815
sh $F runs/pol_pol_uni_k7_s3.85_f0.05_s1_s3.txt     s12_uniform_7of81_3.888
sh $F runs/pol_uni_near7_k8_s3.82_s31.txt          s12_uniform_8of93_3.848
