#!/usr/bin/env bash
# run one sos_probe_sdp.py job single-threaded; its JSON line goes to runs/sos_probe/res/<unique>.json
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1
mkdir -p runs/sos_probe/res
ulimit -v ${SOS_VMEM:-16000000}   # default 16 GB address-space cap per job
f=runs/sos_probe/res/$(date +%s%N)_$$.json
runs/sos_venv/bin/python search/sos_probe_sdp.py "$@" 2>/dev/null | grep '^{' > "$f"
python3 -c "import json,sys
for l in open('$f'):
    d=json.loads(l); print(d['cfg'],d['sub'],'cls='+d['classes'],'B' if d['bound'] else '','orth='+d['orthant'],'ord',d['order'],'T',d['T'],'prod',d['prod'],d['pdeg'],'csp',d['csp'],d['csp0'],d['status'],'lam=%s'%d['lam'],'mono',d['n_mono'],'gram',d['gram'][:3],d['n_gram'],'scal',d['n_scal'],'t',d['t_build'],d['t_solve'])"
