#!/usr/bin/env python3
"""Search benchmark (10-07): best side vs time from fixed seeds, R replicates per seed, for comparing search variants.

Seeds (bench/seeds/, made once by `--make-seeds`, committed): neighbour removals in SQUISH's lineage and two records:
  s108   = record(110) - 2 (best of 4 random removals, fq-quenched)      SQUISH 10.909940073445
  s179   = Couzo's 180 - 1 (best of 4)                                   SQUISH 13.883795490512
  r126   = register 126 (10-05), fq polish only (--no-alm)             SQUISH 11.773303606607
  r237   = register 237 (10-05), fq polish only (--no-alm)             SQUISH 15.903676235191
Variants are named functions in VARIANTS; each replicate is one process with a wall budget (= CPU budget).

  bench_search.py --make-seeds
  bench_search.py --variant base --reps 5 --budget 480 --procs 15 --out runs/bs_base
  bench_search.py --summary runs/bs_base [runs/bs_other ...]
"""
import argparse, json, os, random, statistics as st, tempfile
from jobpool import run_jobs
import mcmin, hop

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = os.path.join(HERE, 'bench/seeds')
TARGET = {'s108': 10.909940073445, 's179': 13.883795490512, 'r126': 11.773303606607, 'r237': 15.903676235191}
CHECK = (60, 120, 240, 480)


def make_seeds():
    os.makedirs(SEEDS, exist_ok=True)
    hop.EXTRA[:] = ['--loosen', '1.0']
    tmp = tempfile.mkdtemp()
    for name, n, k in (('s108', 108, 2), ('s179', 179, 1)):
        s, sq = hop.seed_remove(n, k, 4, random.Random(12345 + n), tmp)
        mcmin.write_deg(f'{SEEDS}/{name}.txt', s, sq); print(name, s)
    for name, n in (('r126', 126), ('r237', 237)):
        s0, sq0 = mcmin.load_deg(f'{hop.BATCH}/n-{n}.txt')
        s, sq, _ = hop.quench(s0, sq0, tmp, extra=('--no-alm',))
        mcmin.write_deg(f'{SEEDS}/{name}.txt', s, sq); print(name, s, 'register', s0)


def v_base(s, sq, budget, rng, tmp, log):
    hop.EXTRA[:] = ['--loosen', '1.0']
    return hop.search(s, sq, budget, rng, tmp, log)


def v_l102(s, sq, budget, rng, tmp, log):
    hop.EXTRA[:] = ['--loosen', '1.02']
    return hop.search(s, sq, budget, rng, tmp, log)


def v_da(s, sq, budget, rng, tmp, log):
    hop.EXTRA[:] = ['--loosen', '1.0']
    return hop.search_da(s, sq, budget, rng, tmp, log)


def v_da102(s, sq, budget, rng, tmp, log):
    hop.EXTRA[:] = ['--loosen', '1.02']
    return hop.search_da(s, sq, budget, rng, tmp, log)


VARIANTS = {'base': v_base, 'l102': v_l102, 'da': v_da, 'da102': v_da102}


def job(j):
    variant, seed, rep, budget, out = j
    s, sq = mcmin.load_deg(f'{SEEDS}/{seed}.txt')
    rng = random.Random(hash((seed, rep)) & 0xffffffff)
    tmp = tempfile.mkdtemp()
    with open(f'{out}/{seed}_r{rep}.jsonl', 'w') as log:
        log.write(json.dumps(dict(seed=seed, rep=rep, start=s, variant=variant)) + '\n')
        b = VARIANTS[variant](s, sq, budget, rng, tmp, log)
    mcmin.write_deg(f'{out}/{seed}_r{rep}_best.txt', *b)
    return dict(seed=seed, rep=rep, best=b[0])


def curves(out):
    res = {}
    for f in sorted(os.listdir(out)):
        if not f.endswith('.jsonl'):
            continue
        L = [json.loads(l) for l in open(os.path.join(out, f))]
        if not L:
            continue
        h, P = L[0], L[1:]
        at = [min([h['start']] + [p['best'] for p in P if p['t'] <= c]) for c in CHECK]
        res.setdefault(h['seed'], []).append(dict(start=h['start'], at=at, props=len(P)))
    return res


def summary(outs):
    for out in outs:
        print(f'== {out}   (median over replicates of best - start; [min, max]; frac = median gain / SQUISH gap)')
        for seed, R in sorted(curves(out).items()):
            st0 = R[0]['start']; gap = st0 - TARGET[seed]
            cols = []
            for i, c in enumerate(CHECK):
                g = [st0 - r['at'][i] for r in R]
                cols.append(f'{c:>4}s {st.median(g):8.1e} [{min(g):.0e},{max(g):.0e}]')
            print(f'  {seed} start {st0:.6f} gap {gap:.1e} | ' + ' | '.join(cols) +
                  f' | frac {st.median([st0 - r["at"][-1] for r in R]) / gap:.2f} | props {st.median([r["props"] for r in R]):.0f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--make-seeds', action='store_true'); ap.add_argument('--variant', default='base')
    ap.add_argument('--seeds', nargs='+', default=list(TARGET)); ap.add_argument('--reps', type=int, default=5)
    ap.add_argument('--budget', type=float, default=480); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--out'); ap.add_argument('--summary', nargs='+')
    a = ap.parse_args()
    if a.make_seeds:
        make_seeds(); raise SystemExit
    if a.summary:
        summary(a.summary); raise SystemExit
    os.makedirs(a.out, exist_ok=True)
    run_jobs(job, [(a.variant, sd, r, a.budget, a.out) for sd in a.seeds for r in range(a.reps)], procs=a.procs,
             timeout=a.budget + 1800)
    summary([a.out])
