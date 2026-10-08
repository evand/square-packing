#!/usr/bin/env python3
"""Disk -> square shape annealing as a start constructor (10-08): anneal (Rust, `anneal run`) then the normal fq finish
(ALM + polish, loosen 1.02) and the grid-obstruction test; compared with controls (same MC with squares throughout;
random placement + fq).  Novelty vs the s(110) references (runs/known110_all.json certified sides, runs/ref110_ex10.json
archive-level sides; 2e-9).

  anneal_test.py --n 110 --variants mid early late sq rand --count 48 --procs 8 --out runs/an1
  anneal_test.py --report runs/an1
"""
import argparse, bisect, collections, json, math, os, statistics as st, subprocess, sys, tempfile, time
from jobpool import run_jobs
import mcmin

HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.join(HERE, 'target-dev/release/anneal')
FQ = os.path.join(HERE, 'target/release/fq')

# variant -> anneal args (schedule over progress t: bP = bp0 (bp1/bp0)^t; r = 0.5 (1 - clamp((t-t0)/(t1-t0)))^p)
V = {
    'mid':   ['--sweeps', '20000', '--t0', '0.2', '--t1', '0.8'],
    'early': ['--sweeps', '20000', '--t0', '0.0', '--t1', '0.4'],
    'late':  ['--sweeps', '20000', '--t0', '0.6', '--t1', '0.95'],
    'slow':  ['--sweeps', '100000', '--t0', '0.2', '--t1', '0.8'],
    'fast':  ['--sweeps', '4000', '--t0', '0.2', '--t1', '0.8'],
    's300':  ['--sweeps', '300000', '--t0', '0.2', '--t1', '0.8'],
    's1m':   ['--sweeps', '1000000', '--t0', '0.2', '--t1', '0.8'],
    's300e': ['--sweeps', '300000', '--t0', '0.0', '--t1', '0.5'],
    's300l': ['--sweeps', '300000', '--t0', '0.5', '--t1', '0.95'],
    's300p': ['--sweeps', '300000', '--t0', '0.2', '--t1', '0.8', '--bp1', '30000'],
    'sq':    ['--sweeps', '20000', '--t0', '0', '--t1', '0'],              # control: squares throughout, same MC
    'rand':  ['--sweeps', '0', '--t0', '0', '--t1', '0'],                  # control: random placement only
}


def job(j):
    n, var, seed, out, extra = j
    from layout import full_lines
    a = f'{out}/{var}_{seed:03d}.txt'; b = f'{out}/{var}_{seed:03d}.q.txt'
    t0 = time.time()
    r1 = subprocess.run([AN, 'run', '--n', str(n), '--seed', str(seed), '--out', a, *V[var], *extra], capture_output=True, text=True, timeout=3600)
    da = json.loads(r1.stdout)
    t1 = time.time()
    r2 = subprocess.run([FQ, 'quench', '--in', a, '--out', b], capture_output=True, text=True, timeout=3600)
    d = json.loads(r2.stdout)
    t2 = time.time()
    ok = d['min_gap'] >= 0 and d['min_wall'] >= 0
    k = math.ceil(math.sqrt(n) - 1e-12)
    s, sq = mcmin.load_deg(b)
    import explore
    return dict(var=var, seed=seed, s_an=da['s'], an=da, s=d['s'], ok=ok, lines=int(full_lines(b, k)), roles=explore.roles(sq),
                t_an=t1 - t0, t_fq=t2 - t1, path=b)


def report(out):
    R = [json.loads(l) for l in open(f'{out}/res.jsonl')]
    meta = json.load(open(f'{out}/meta.json'))
    n = meta['n']; k = math.ceil(math.sqrt(n) - 1e-12)
    refs = {}
    for name in ('known110_all', 'ref110_ex10'):
        p = f'{HERE}/runs/{name}.json'
        if n == 110 and os.path.exists(p):
            refs[name] = sorted(json.load(open(p)))
    def novel(s, ref):
        i = bisect.bisect_left(ref, s - 2e-9)
        return not (i < len(ref) and abs(ref[i] - s) < 2e-9)
    by = collections.defaultdict(list)
    for r in R:
        by[r['var']].append(r)
    print(f'== {out}: n = {n} (k = {k}); per variant: count, sides p10/50/90, best, frac < k, frac grid (full lines / chain), '
          f'distinct sub-k (2e-9), novel vs refs, CPU s per result (anneal + fq)')
    for v, L in by.items():
        L = [r for r in L if r['ok']]
        ss = sorted(r['s'] for r in L)
        sub = [r for r in L if r['s'] < k and not r['lines']]
        grid = sum(r['lines'] > 0 for r in L)
        ds = sorted(set(round(r['s'], 9) for r in sub))
        nov = {nm: sum(novel(s, ref) for s in ds) for nm, ref in refs.items()}
        q = lambda p: ss[min(len(ss) - 1, int(p * len(ss)))]
        roles = collections.Counter(tuple(r['roles']) for r in sub).most_common(3)
        print(f'  {v:6} {len(L):3d}  {q(.1):.5f}/{q(.5):.5f}/{q(.9):.5f}  best {ss[0]:.9f}  sub-k {len(sub) / len(L):.2f}  '
              f'grid {grid / len(L):.2f}  distinct sub-k {len(ds)}  novel {nov}  CPU {st.mean(r["t_an"] + r["t_fq"] for r in L):.1f} s '
              f'(an {st.mean(r["t_an"] for r in L):.1f})  top roles {roles}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=110); ap.add_argument('--variants', nargs='+', default=['mid', 'sq', 'rand'])
    ap.add_argument('--count', type=int, default=48); ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--out'); ap.add_argument('--report'); ap.add_argument('--seed0', type=int, default=1)
    ap.add_argument('--extra', default='', help='extra anneal args for every variant (e.g. "--bp1 1e4")')
    a = ap.parse_args()
    if a.report:
        report(a.report); raise SystemExit
    os.makedirs(a.out, exist_ok=True)
    json.dump(dict(n=a.n, variants={v: V[v] for v in a.variants}, extra=a.extra, argv=sys.argv), open(f'{a.out}/meta.json', 'w'))
    jobs = [(a.n, v, a.seed0 + i, a.out, a.extra.split()) for i in range(a.count) for v in a.variants]
    res = run_jobs(job, jobs, procs=a.procs, timeout=7200)
    with open(f'{a.out}/res.jsonl', 'a') as f:
        for r in res:
            if isinstance(r, dict):
                f.write(json.dumps(r) + '\n')
    report(a.out)
