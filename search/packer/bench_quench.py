#!/usr/bin/env python3
"""Head-to-head quench benchmark (10-07): fq (smooth separating-line ALM, Rust) vs mcmin.quench (soft squeeze + slp2).

Same kicked state into both (register packing ../exact/batch/inputs/n-N.txt, Gaussian kick sigma, angles 20 sigma deg).
Also fq followed by slp2 (is fq's output already a jammed local minimum, and how much does slp2 still find?).
Per trial: side, wall time, grid obstruction (layout.full_lines, incl. staggered chains).

  bench_quench.py --n 110 270 --sigmas 0.01 0.03 0.1 --trials 6 --procs 15 --out runs/bq0
"""
import argparse, json, math, os, random, subprocess, tempfile, time
from jobpool import run_jobs
import mcmin

HERE = os.path.dirname(os.path.abspath(__file__))
FQ = os.path.join(HERE, 'target/release/fq')
BATCH = os.path.join(HERE, '../exact/batch/inputs')


def fq(path_in, path_out, extra=()):
    t = time.time()
    r = subprocess.run([FQ, 'quench', '--in', path_in, '--out', path_out, *extra], capture_output=True, text=True, timeout=600)
    d = json.loads(r.stdout)
    d['wall'] = time.time() - t
    return d


def slp2_polish(path, budget):
    from rigid import load as load_rad
    from slp2 import slp2
    s, sq = load_rad(path)
    t = time.time()
    s3, sq3, _ = slp2(s, sq, R=1e-3, rmin=1e-8, budget=budget, info={})
    return s3, time.time() - t


def job(j):
    n, sig, t, methods, budget = j
    from layout import full_lines
    s0, sq0 = mcmin.load_deg(f'{BATCH}/n-{n}.txt')
    rng = random.Random(n * 100003 + int(sig * 1e5) * 101 + t)
    sq = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), a + rng.gauss(0, 20 * sig)) for x, y, a in sq0]
    k = math.ceil(math.sqrt(n) - 1e-12)
    tmp = tempfile.mkdtemp()
    kin = os.path.join(tmp, 'kick.txt')
    mcmin.write_deg(kin, s0, sq)
    out = dict(n=n, sig=sig, t=t)
    if 'fq' in methods:
        fo = os.path.join(tmp, 'fq.txt')
        d = fq(kin, fo)
        out['fq'] = dict(s=d['s'], sec=d['wall'], evals=d['evals'], lines=bool(full_lines(fo, k)))
        if 'fq+slp2' in methods:
            s3, dt = slp2_polish(fo, budget)
            out['fq+slp2'] = dict(s=s3, sec=d['wall'] + dt)
    if 'slp2' in methods:
        mcmin.SLP_BUDGET = budget
        t0 = time.time()
        r = mcmin.quench(s0, sq, k, tmp)
        out['slp2'] = dict(s=r[0], sec=time.time() - t0, lines=bool(r[2]['lines'])) if r else dict(s=None)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+', default=[110])
    ap.add_argument('--sigmas', type=float, nargs='+', default=[0.01, 0.03, 0.1])
    ap.add_argument('--trials', type=int, default=6)
    ap.add_argument('--methods', nargs='+', default=['fq', 'fq+slp2', 'slp2'])
    ap.add_argument('--budget', type=float, default=180)
    ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--out', default='runs/bq0')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = [(n, sg, t, a.methods, a.budget) for n in a.n for sg in a.sigmas for t in range(a.trials)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=1800)
    res = [r for r in res if r]
    json.dump(res, open(f'{a.out}/bq.json', 'w'), indent=0)
    for n in a.n:
        for sg in a.sigmas:
            R = [r for r in res if r['n'] == n and r['sig'] == sg]
            print(f'n {n} sigma {sg}')
            for r in R:
                row = '  '.join(f"{m} {r[m]['s']:.10f} ({r[m]['sec']:.1f}s)" for m in a.methods if m in r and r[m].get('s'))
                print('   ', row)
