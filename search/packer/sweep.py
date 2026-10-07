#!/usr/bin/env python3
"""Neighbour sweep over register records (2026-10-06): cheap probes for a better packing at each open n.

Probes per n (quench = mcmin.quench: loosen 2 %, soft squeeze, slp2):
  rm    record(n+1) minus one random square
  kick  record(n) + Gaussian sigma (0.01, 0.03; angles 20 sigma deg)
Flags any jammed, line-free result below the register value of n (exact/batch/candidates_all.json, register as of 10-05:
re-check jlevy's register before believing anything).

  sweep.py --lo 30 --hi 200 --rm 8 --kick 4 --out runs/sw1 --procs 4
"""
import argparse, json, math, os, random, tempfile
from jobpool import run_jobs
import mcmin

SIGMAS = [0.01, 0.03]
BATCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../exact/batch')


def job(j):
    n, kind, seed, out = j
    rng = random.Random(seed)
    k = math.ceil(math.sqrt(n) - 1e-12)
    if kind == 'rm':
        s, sq = mcmin.load_deg(f'{BATCH}/inputs/n-{n + 1}.txt')
        i = rng.randrange(len(sq)); sq = sq[:i] + sq[i + 1:]; desc = f'rm {i}'
    else:
        sig = SIGMAS[seed % len(SIGMAS)]
        s, sq = mcmin.load_deg(f'{BATCH}/inputs/n-{n}.txt')
        sq = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), t + rng.gauss(0, 20 * sig)) for x, y, t in sq]; desc = f'kick {sig}'
    r = mcmin.quench(s, sq, k, tempfile.mkdtemp())
    if r is None:
        return dict(n=n, kind=kind, seed=seed, desc=desc, status='quench failed')
    s2, sq2, info = r
    fn = f'{out}/n{n}_{kind}_{seed}.txt'
    mcmin.write_deg(fn, s2, sq2)
    return dict(n=n, kind=kind, seed=seed, desc=desc, s=s2, lines=info['lines'], jam=info['jammed'], path=fn, status='ok')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int, default=30); ap.add_argument('--hi', type=int, default=200)
    ap.add_argument('--rm', type=int, default=8); ap.add_argument('--kick', type=int, default=4)
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--out', required=True); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--sigmas', default='0.01,0.03')
    a = ap.parse_args()
    SIGMAS[:] = [float(x) for x in a.sigmas.split(',')]
    os.makedirs(a.out, exist_ok=True)
    reg = {d['n']: d for d in json.load(open(f'{BATCH}/candidates_all.json'))}
    ns = [n for n in range(a.lo, a.hi + 1) if reg[n]['status'] == 'open' and abs(float(reg[n]['register']) - round(float(reg[n]['register']))) > 1e-9
          and os.path.exists(f'{BATCH}/inputs/n-{n}.txt') and os.path.exists(f'{BATCH}/inputs/n-{n + 1}.txt')]
    jobs = [(n, 'rm', a.seed * 10 ** 6 + n * 100 + t, a.out) for n in ns for t in range(a.rm)] + \
           [(n, 'kick', a.seed * 10 ** 6 + n * 100 + 50 + t, a.out) for n in ns for t in range(a.kick)] if a.kick else []
    jobs.sort(key=lambda j: -j[0])                  # longest (largest n) first
    print(f'{len(ns)} n values, {len(jobs)} jobs', flush=True)
    res = run_jobs(job, jobs, procs=a.procs, timeout=900, on_timeout=lambda j: dict(n=j[0], kind=j[1], seed=j[2], status='killed'),
                   on_error=lambda j, e: dict(n=j[0], kind=j[1], seed=j[2], status=f'failed {e}'))
    res = [r or dict(n=j[0], kind=j[1], seed=j[2], status='failed') for r, j in zip(res, jobs)]
    json.dump(res, open(f'{a.out}/sweep.json', 'w'), indent=0)
    print('n  register  best(rm)  best(kick)  delta', flush=True)
    for n in ns:
        R = float(reg[n]['register'])
        b = {kd: min((r['s'] for r in res if r['n'] == n and r.get('status') == 'ok' and r['kind'] == kd and not r['lines']),
                     default=float('inf')) for kd in ('rm', 'kick')}
        m = min(b.values())
        flag = '  <<< BELOW REGISTER' if m < R - 1e-10 else ''
        print(f'{n} {R:.10f} {b["rm"]:.10f} {b["kick"]:.10f} {m - R:+.2e}{flag}', flush=True)
