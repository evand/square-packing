#!/usr/bin/env python3
"""Ground-truth check of the explorer's --kick-retry (10-11).  n = 110 census points with a known destination:
  notmin   stored explorer outputs that exactsolve called not jammed and descend.py later carried to a certified minimum
           (runs/census/n110/descend.jsonl: pid -> S of the minimum reached)
  control  certified minimum points (descend outputs, store rows source 'descend:census'): kick-retry should not move them
Each point: explore.kick_retry with --tries kicked re-quenches (side after each try recorded).  Outcome per point: side after
t tries vs the destination S (reached = within 1e-8), vs every known minimum (table `minimum`), measured CPU.

  kr_eval.py --out runs/kr_eval --per 60 --tries 4 --procs 4
"""
import argparse, json, math, os, random, sys, tempfile, time, bisect
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def job(j):
    pid, kind, S_dest, tries, sigma, seed = j
    import explore
    from pk.store import Store
    from pk.manifest import cpu_seconds
    p = Store().get(pid)
    sq = [tuple(map(float, r)) for r in p.sq]
    c0 = cpu_seconds()
    with tempfile.TemporaryDirectory() as tmp:
        s, _, info = explore.kick_retry(float(p.s), sq, random.Random(seed), tmp, math.ceil(math.sqrt(p.n) - 1e-12), tries, sigma)
    return dict(pid=pid, kind=kind, s_in=float(p.s), S_dest=S_dest, s_out=s, traj=info.get('traj', []), raw=info.get('raw', []),
                drops=info['drops'], cpu=round(cpu_seconds() - c0, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='runs/kr_eval'); ap.add_argument('--per', type=int, default=60)
    ap.add_argument('--tries', type=int, default=4); ap.add_argument('--sigma', type=float, default=1e-4)
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--seed', type=int, default=1011)
    ap.add_argument('--report')
    a = ap.parse_args()
    if a.report:
        return report(a.report)
    import scratch; scratch.scope('kr_eval')
    from pk import manifest
    from pk.store import Store
    man = manifest.start(a.out, 'kr_eval', seed=a.seed, args=vars(a), procs=a.procs)
    rng = random.Random(a.seed)
    D = [json.loads(l) for l in open(f'{HERE}/runs/census/n110/descend.jsonl')]
    D = [d for d in D if d.get('cls') == 'certified' and d.get('S') and d.get('pid')]
    notmin = rng.sample(D, min(a.per, len(D)))
    st = Store()
    ctrl_ids = [r['id'] for r in st.query(n=110, source='descend:census', status='certified', limit=None)]
    ctrl = rng.sample(ctrl_ids, min(a.per, len(ctrl_ids)))
    S_of = {r['pid_out']: r['S'] for r in D if r.get('pid_out')}
    jobs = [(d['pid'], 'notmin', d['S'], a.tries, a.sigma, rng.randrange(1 << 30)) for d in notmin]
    jobs += [(pid, 'control', S_of.get(pid), a.tries, a.sigma, rng.randrange(1 << 30)) for pid in ctrl]
    rng.shuffle(jobs)
    with open(f'{a.out}/trials.jsonl', 'w') as f, ProcessPoolExecutor(a.procs) as ex:
        for r in ex.map(job, jobs):
            f.write(json.dumps(r) + '\n'); f.flush()
    man.finish(summary=dict(points=len(jobs)))
    report(a.out)


def report(out):
    from pk.store import Store
    R = [json.loads(l) for l in open(f'{out}/trials.jsonl')]
    M = sorted(float(m['S'][:20]) for m in Store().minima(110))
    def on_min(s, tol=1e-8):
        i = bisect.bisect_left(M, s - tol)
        return i < len(M) and abs(M[i] - s) < tol
    import statistics as st
    for kind in ('notmin', 'control'):
        X = [r for r in R if r['kind'] == kind]
        if not X: continue
        T = max(len(r['traj']) for r in X)
        print(f'{kind}: {len(X)} points, CPU median {st.median(r["cpu"] for r in X):.1f}s (mean {st.mean(r["cpu"] for r in X):.1f}s for {T} tries)')
        print(f'  input on a known minimum (1e-8): {sum(on_min(r["s_in"]) for r in X)}')
        for t in range(T + 1):
            s_t = [(r['traj'][t - 1] if t else r['s_in']) for r in X]
            dest = [r['S_dest'] and abs(s - float(r['S_dest'][:20])) < 1e-8 for r, s in zip(X, s_t)]
            below = [r['S_dest'] and s < float(r['S_dest'][:20]) - 1e-8 for r, s in zip(X, s_t)]
            onm = [on_min(s) for s in s_t]
            print(f'  after {t} tries: at descend destination {sum(map(bool, dest))}, below it {sum(map(bool, below))}, '
                  f'on any known minimum {sum(onm)}, moved from input {sum(abs(s - r["s_in"]) > 2e-9 for r, s in zip(X, s_t))}')
        dr = sorted(r['s_in'] - r['s_out'] for r in X)
        print(f'  side drop: median {dr[len(dr)//2]:.2e}, max {dr[-1]:.2e}; min {dr[0]:.2e}')


if __name__ == '__main__':
    main()
