#!/bin/sh
# Launch the nu_f bracket batch used for search/CEILING.md.  Run from the repo root.
#   9 exact-mode runs (search/nu_f.py, 3 procs each)  +  4 cell-mode runs (search/lp_search.py, 1 proc each)
# Each run is time-boxed to ~35 min; results land in runs/.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p runs
TL=${TL:-2100}
MODE=${MODE:-all}
[ "$MODE" = cell ] || for s in 3920/997 3.94 3.945 3.95 3.955 3.96 3.97 3.98 4; do
  tag=$(echo "$s" | tr '/' '_')
  nohup python3 search/nu_f.py run "$s" "$TL" "e$tag" > "runs/stdout_e$tag.txt" 2>&1 &
done
[ "$MODE" = exact ] || for s in 3.93 3.94 3.95 3.96; do
  nohup python3 -c "
import sys; sys.path.insert(0,'search'); import lp_search as L
s=$s; best=L.run(s,fine=0.005,eta=0.005,dt=0.005,tlimit=$TL-100,tag='c$s')
val,x,m=best; tot=L.write_cert(m,x,s,'runs/cert_c$s.txt'); print(f'RESULT s={s} TOTAL={val:.5f} exported={tot:.5f}')
" > "runs/stdout_c$s.txt" 2>&1 &
done
echo launched
