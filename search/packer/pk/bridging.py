"""Do recombination children land between their parents?  For battery trials with mates and a polished sub-k child:
d(child, parent), d(child, nearest mate), parent-mate distance; binned by mean mate distance.

  bridging.py RUN [RUN...] [--max 300]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, collections, json, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pk.store import Store
from pk.packing import Packing
from pk.crossing import match_distance

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'battery')


def job(r):
    st = Store()
    c = Packing(r['s'], np.array(r['fin']))
    dp = match_distance(c, st.get(r['pid']))
    dm = min(match_distance(c, st.get(m)) for m in r['mates'])
    return r['mate_dist'], dp, dm


ap = argparse.ArgumentParser(); ap.add_argument('runs', nargs='+'); ap.add_argument('--max', type=int, default=300)
a = ap.parse_args()
for run in a.runs:
    R = [json.loads(l) for l in open(f'{ROOT}/{run}/trials.jsonl')]
    R = [r for r in R if 'fin' in r and r.get('mates')][:a.max]
    with ProcessPoolExecutor(16) as ex:
        out = list(ex.map(job, R, chunksize=4))
    B = collections.defaultdict(list)
    for md, dp, dm in out:
        b = '<0.02' if md < 0.02 else '<0.05' if md < 0.05 else '<0.1' if md < 0.1 else '>=0.1'
        B[b].append((md, dp, dm))
    print(f'{run}: {len(out)} sub-k children with mates')
    print('   mate dist   N   median d(child,parent)  median d(child,nearest mate)  between (both < 0.8 x mate dist)')
    for b in ['<0.02', '<0.05', '<0.1', '>=0.1']:
        v = B[b]
        if v:
            between = sum(dp < 0.8 * md and dm < 0.8 * md for md, dp, dm in v)
            print(f'   {b:8s} {len(v):4d}   {np.median([x[1] for x in v]):.4f}                  {np.median([x[2] for x in v]):.4f}'
                  f'                        {between}/{len(v)}')
