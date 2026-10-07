#!/bin/sh
# post-hoc rigorous lower bounds from every final certificate of the batch (run from repo root)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for c in ${1:-runs/nuf_e*_best.txt runs/cert_c*.txt}; do
  tag=$(basename "$c" .txt)
  nohup python3 search/nu_f.py lower "$c" "$tag" 0.01 3000 > "runs/stdout_lower_$tag.txt" 2>&1 &
done
echo launched
