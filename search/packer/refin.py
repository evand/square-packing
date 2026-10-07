#!/usr/bin/env python3
"""Re-finish packings with slp2 in parallel: refin.py out_dir file1 file2 ..."""
import sys, os, multiprocessing as mp
from rigid import load
from slp2 import slp2
from slp import save

def job(p):
    s0, sq = load(p)
    s, sq2, ds0 = slp2(s0, sq, R=1e-3, rmin=1e-9, budget=600)
    out = os.path.join(sys.argv[1], os.path.basename(os.path.dirname(p)) + '_' + os.path.basename(p))
    save(out, s, sq2)
    return p, s0, s, ds0

if __name__ == '__main__':
    from jobpool import run_jobs
    if True:
        for p, s0, s, ds0 in sorted((r for r in run_jobs(job, sys.argv[2:], procs=14, timeout=1200) if r), key=lambda r: r[2]):
            print(f'{s0:.7f} -> {s:.7f}  jammed={ds0 > -1e-12}  {p}', flush=True)
