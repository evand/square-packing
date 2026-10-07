#!/usr/bin/env python3
"""Deep census on a few n (10-06, the deep leg of the "T": wide = sweep.py / rediscover.py over n < 324).

For each (n, start): Gaussian kicks (sigma list, angles 20 sigma deg) -> mcmin.quench (slp2 budget raised) -> every output saved,
with side, grid obstruction, jam, timeout, and matching distance to the start and to the best known (rediscover.distance).
Record hunting and sampler characterisation from the same data: per (n, start, sigma) best side, distinct minima (f64 side to
1e-9), f1/N, distance distribution.

  deepcensus.py --n 266 270 272 --trials 25 --procs 15 --out runs/dc1
"""
import argparse, collections, json, math, os, random, tempfile
from jobpool import run_jobs
import mcmin, rediscover

SIGMAS = (0.003, 0.01, 0.02, 0.03, 0.05)
HERE = os.path.dirname(os.path.abspath(__file__))
OURS = {266: 'candidates/sw2_n266.exact.txt', 270: 'candidates/sw2_n270.exact.txt', 272: 'candidates/rd1_n272.exact.txt'}


def starts(n):
    out = {}
    f = f'{rediscover.SITE}/square-{n}.json'
    reg_s, reg = mcmin.load_deg(f'{rediscover.BATCH}/inputs/n-{n}.txt')
    if os.path.exists(f):
        s, sq = rediscover.load_site(n)
        if s > reg_s + 1e-9:
            out['mirror'] = (s, sq)
    out['register'] = (reg_s, reg)
    if n in OURS:
        out['ours'] = mcmin.load_deg(os.path.join(HERE, OURS[n]))
    return out


def job(j):
    n, label, sig, t, out, budget = j
    mcmin.SLP_BUDGET = budget
    S = starts(n)
    s0, sq0 = S[label]
    best_s, best = min(S.values(), key=lambda v: v[0])
    rng = random.Random(n * 1000003 + hash(label) % 9973 * 101 + int(sig * 1e5) * 7 + t)
    k = math.ceil(math.sqrt(n) - 1e-12)
    sq = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), a + rng.gauss(0, 20 * sig)) for x, y, a in sq0]
    r = mcmin.quench(s0, sq, k, tempfile.mkdtemp())
    if r is None:
        return dict(n=n, start=label, sig=sig, t=t, status='quench failed')
    s, sq2, info = r
    fn = f'{out}/n{n}_{label}_{sig:g}_{t}.txt'
    mcmin.write_deg(fn, s, sq2)
    return dict(n=n, start=label, sig=sig, t=t, s=s, lines=info['lines'], jam=info['jammed'], timeout=info['timeout'],
                d_start=rediscover.distance(sq0, s0, sq2, s), d_best=rediscover.distance(best, best_s, sq2, s), path=fn,
                status='ok')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+', required=True); ap.add_argument('--trials', type=int, default=25)
    ap.add_argument('--procs', type=int, default=15); ap.add_argument('--budget', type=float, default=180)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = [(n, lab, sg, t, a.out, a.budget) for n in a.n for lab in starts(n) for sg in SIGMAS for t in range(a.trials)]
    print(f'{len(jobs)} jobs', flush=True)
    res = run_jobs(job, jobs, procs=a.procs, timeout=4 * a.budget + 300,
                   on_timeout=lambda j: dict(n=j[0], start=j[1], sig=j[2], t=j[3], status='killed'),
                   on_error=lambda j, e: dict(n=j[0], start=j[1], sig=j[2], t=j[3], status=f'failed {e}'))
    res = [r or dict(n=j[0], start=j[1], sig=j[2], t=j[3], status='failed') for r, j in zip(res, jobs)]
    json.dump(res, open(f'{a.out}/dc.json', 'w'), indent=0)
    for n in a.n:
        S = starts(n)
        ref = min(v[0] for v in S.values())
        print(f'n = {n}: best known {ref:.10f}')
        for lab in S:
            for sg in SIGMAS:
                R = [r for r in res if r['n'] == n and r['start'] == lab and r['sig'] == sg and r['status'] == 'ok' and not r['lines']]
                if not R:
                    continue
                c = collections.Counter(round(r['s'], 9) for r in R)
                f1 = sum(1 for v in c.values() if v == 1)
                b = min(r['s'] for r in R)
                print(f'  {lab:8s} σ={sg:<5g} {len(R):3d} ok  distinct {len(c):3d}  f1/N {f1 / len(R):.2f}  best {b:.10f} ({b - ref:+.2e})'
                      f'  below best-known {sum(r["s"] < ref - 1e-9 for r in R)}  med d_start {sorted(r["d_start"] for r in R)[len(R) // 2]:.3f}'
                      f'  timeouts {sum(r["timeout"] for r in R)}', flush=True)
