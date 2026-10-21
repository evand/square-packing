#!/usr/bin/env python3
"""Junction moves A/B (10-21; Evan: "Have we tried local melt, or delete and refill, at the junction region?").

Paired single proposals from fixed parents, through the explorer's own screen + polish (explore._screen_polish):
  melt     explore.melt_side, random centre (the explorer's melt; control)
  jmelt    the same draws (radius, rounding, seed), centre = one of the 8 squares nearest the two-band junction
  jrefill  explore.mv_jrefill: drop the 3-10 squares nearest the (jittered) junction, re-add as many at pose holes
Same seed for the three arms of a trial.  Parents: store ids (--ids / --ids-file), one line each.
Per trial: status, polished side, lines, ds vs parent, matching distance to the parent (other basin or not).

  junction_ab.py run --ids-file F --seeds 10 --procs 5 --out runs/jab/j1 [--cpu-h 6]
  junction_ab.py report runs/jab/j1
"""
import argparse, collections, json, math, os, random, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ARMS = ('melt', 'jmelt', 'jrefill')


def job(j):
    path, seed, arm = j
    import explore, mcmin
    from pk.manifest import cpu_seconds
    c0 = cpu_seconds()
    s, sq = mcmin.load_deg(path)
    n = len(sq); k = math.ceil(math.sqrt(n) - 1e-12)
    rng = random.Random(seed)
    if arm == 'jrefill':
        s2, prop, desc = explore.mv_jrefill(s, sq, rng, explore.A)
    else:
        s2, prop, desc = explore.melt_side(s, sq, rng, explore.A, junction=(arm == 'jmelt'))
    smax = max(k + 0.05, s + 0.05)
    with tempfile.TemporaryDirectory() as tmp:
        r = explore._screen_polish(s2, prop, arm, desc, rng, k, [], 0.0, smax, float('inf'), (), tmp, None)
    out = dict(path=path, seed=seed, arm=arm, n=n, s_par=s, st=r['status'], desc=desc, s_screen=r.get('s_screen'),
               s=r.get('s'), lines=r.get('lines'))
    if r.get('sq') is not None and r.get('s') is not None and not r.get('lines'):
        from pk.crossing import match_distance
        from pk.packing import Packing
        import numpy as np
        out['ds'] = r['s'] - s
        out['self'] = abs(r['s'] - s) < 2e-9
        out['d'] = match_distance(Packing(r['s'], np.array(r['sq'])), Packing(s, np.array(sq)))
        if r['s'] < s - 2e-9 or not out['self']:
            os.makedirs(os.path.join(A_OUT, 'kids'), exist_ok=True)
            f = os.path.join(A_OUT, 'kids', f'{os.path.basename(path)[:-4]}_{arm}_{seed}.txt')
            mcmin.write_deg(f, r['s'], r['sq']); out['kid'] = f
    out['cpu'] = round(cpu_seconds() - c0, 2)
    return out


A_OUT = None


def _init(out):
    global A_OUT
    A_OUT = out


def run(a):
    from pk.store import Store
    from pk import manifest
    os.makedirs(f'{a.out}/parents', exist_ok=True)
    ids = list(a.ids or []) + ([l.split()[0] for l in open(a.ids_file) if l.strip()] if a.ids_file else [])
    st = Store(); P = []
    for i in ids:
        p = f'{a.out}/parents/p{int(i)}.txt'; st.get(int(i)).write(p); P.append(p)
    rng = random.Random(a.seed)
    trials = [[(p, sd, arm) for arm in ARMS] for p in P for sd in [rng.randrange(1 << 30) for _ in range(a.seeds)]]
    rng.shuffle(trials)                                     # triples stay together: partial runs stay balanced
    flat = [t for tr in trials for t in tr]
    man = manifest.start(a.out, 'junction_ab', seed=a.seed, args=vars(a), procs=a.procs)
    budget = a.cpu_h * 3600 if a.cpu_h else float('inf'); spent = 0.0
    log = open(f'{a.out}/trials.jsonl', 'a'); it = iter(flat); pend = set(); stop = None
    with ProcessPoolExecutor(a.procs, initializer=_init, initargs=(a.out,)) as ex:
        while True:
            while len(pend) < a.procs and spent < budget:
                j = next(it, None)
                if j is None: break
                pend.add(ex.submit(job, j))
            if not pend: break
            done, pend = wait(pend, return_when=FIRST_COMPLETED)
            for f in done:
                try:
                    r = f.result(); spent += r.get('cpu', 0)
                except Exception as e:
                    r = dict(error=repr(e)[:300])
                log.write(json.dumps(r) + '\n'); log.flush()
    if spent >= budget:
        stop = 'cpu'
    man.finish(summary=dict(stop=stop or 'done', cpu_trials_s=round(spent, 1)))
    report(a.out)


def report(out):
    R = [json.loads(l) for l in open(f'{out}/trials.jsonl')]
    err = sum('error' in r for r in R); R = [r for r in R if 'arm' in r]
    by = collections.defaultdict(dict)
    for r in R: by[(r['path'], r['seed'])][r['arm']] = r
    full = [v for v in by.values() if len(v) == len(ARMS)]
    print(f'{out}: {len(R)} trials ({err} errors), {len(full)} complete triples')
    groups = collections.defaultdict(list)
    for v in full: groups[v['melt']['n']].append(v)
    for n in sorted(groups):
        T = groups[n]
        print(f'\nn = {n}: {len(T)} triples, {len(set(v["melt"]["path"] for v in T))} parents')
        for arm in ARMS:
            X = [v[arm] for v in T]
            ok = [x for x in X if 'ds' in x]
            away = [x for x in ok if not x['self']]
            below = [x for x in ok if x['ds'] < -2e-9]
            keys = {round(x['s'], 9) for x in away}
            best = min((x['ds'] for x in ok), default=float('nan'))
            cpu = sum(x.get('cpu', 0) for x in X)
            print(f'  {arm:8s} line-free {len(ok):4d}/{len(X)}  other basin {len(away):4d}  below parent {len(below):4d}  '
                  f'distinct new sides {len(keys):4d}  best ds {best:+.2e}  cpu {cpu / 3600:.2f} h')
        for a_, b_ in (('jmelt', 'melt'), ('jrefill', 'melt')):
            w = sum(1 for v in T if v[a_].get('ds', 9) < v[b_].get('ds', 9) - 2e-9)
            l_ = sum(1 for v in T if v[b_].get('ds', 9) < v[a_].get('ds', 9) - 2e-9)
            print(f'  paired {a_} vs {b_}: {a_} lower side {w}, {b_} lower {l_}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run'); r.add_argument('--ids', nargs='*'); r.add_argument('--ids-file')
    r.add_argument('--seeds', type=int, default=10); r.add_argument('--procs', type=int, default=4)
    r.add_argument('--seed', type=int, default=2110); r.add_argument('--cpu-h', type=float, default=None)
    r.add_argument('--out', required=True)
    p = sub.add_parser('report'); p.add_argument('out')
    a = ap.parse_args()
    if a.cmd == 'report':
        report(a.out)
    else:
        import scratch
        scratch.scope('junction_ab')
        run(a)
