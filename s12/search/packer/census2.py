#!/usr/bin/env python3
"""Basin census with the LP squeeze: reference + Gaussian sigma (positions; angles 20*sigma deg) -> slp2 -> side, jam status.
  census2.py ref.txt --sigmas 1e-4,1e-3,1e-2,3e-2 --trials 8 --procs 14 --out prefix
"""
import argparse, math, random, multiprocessing as mp
from rigid import load
from slp2 import slp2
from slp import save


def job(args):
    ref, sig, t, out = args
    s0, C = load(ref)
    rng = random.Random(1000 * t + int(sig * 1e6))
    P = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), a + math.radians(rng.gauss(0, 20 * sig))) for x, y, a in C]
    s, sq, ds0 = slp2(s0 + 4 * sig, P, R=1e-3, rmin=1e-9, budget=600)
    save(f'{out}_{sig:g}_{t}.txt', s, sq)
    return sig, t, s, ds0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('ref'); ap.add_argument('--sigmas', default='1e-4,1e-3,1e-2,3e-2')
    ap.add_argument('--trials', type=int, default=8); ap.add_argument('--procs', type=int, default=14)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    jobs = [(a.ref, float(sg), t, a.out) for sg in a.sigmas.split(',') for t in range(a.trials)]
    from jobpool import run_jobs
    res = [r for r in run_jobs(job, jobs, procs=a.procs, timeout=1200) if r is not None]
    for sg in a.sigmas.split(','):
        r = sorted(x[2] for x in res if x[0] == float(sg))
        jam = sum(1 for x in res if x[0] == float(sg) and x[3] > -1e-12)
        print(f'sigma={sg}: sides {[round(v, 7) for v in r]}  jammed {jam}/{len(r)}', flush=True)
