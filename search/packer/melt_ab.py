#!/usr/bin/env python3
"""Paired A/B of the explorer melt move before / after the 10-11 side fix (PACKER.md 10-11).

Same parent, same seed -> identical melt (centre, radius, rmax, anneal seed).  Arms differ only in the side the melted
coordinates are screened at: 'old' = parent's side (the 10-08..10-11 explorer bug), 'new' = the anneal's returned side.
Both go through explore._screen_polish with the explorer defaults (same = 0, smax = k + 0.05, no polish margin).

Per proposal: status, screened / polished side, side change vs parent, grid lines, parent-basin return (|ds| < 2e-9),
matching distance child -> parent (pk.crossing.match_distance, as move_spread.py), anneal box growth.

  melt_ab.py --out runs/melt_ab --procs 4 --minutes 12 [--per-n 20]
  melt_ab.py --report runs/melt_ab
"""
import argparse, glob, json, math, os, random, sys, tempfile, time, collections
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

POOLS = {110: 'runs/lineage/L2/*/archive.jsonl',            # sub-11 explorer entries (hunt-recipe archives, 10-10)
         130: 'runs/be_melt/base/s130_r*/archive.jsonl'}    # the be_melt s130 task, base arm (no melt in their lineage)


def parents(n, per, rng):
    k = math.ceil(math.sqrt(n) - 1e-12)
    E = []
    for f in sorted(glob.glob(os.path.join(HERE, POOLS[n]))):
        d = os.path.dirname(f)
        for l in open(f):
            e = json.loads(l)
            if e['s'] < k:
                E.append((e['s'], os.path.join(d, os.path.basename(e['path']))))
    E = sorted(set(E))
    step = max(1, len(E) // per)                             # spread over the side range
    pick = E[::step][:per]
    return [p for _, p in pick]


def job(j):
    path, seed, arm = j
    import explore, mcmin
    s, sq = mcmin.load_deg(path)
    n = len(sq); k = math.ceil(math.sqrt(n) - 1e-12)
    rng = random.Random(seed)
    s2, prop, desc = explore.melt_side(s, sq, rng, explore.A)
    s_screen_at = s2 if arm == 'new' else s
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        r = explore._screen_polish(s_screen_at, prop, 'melt', desc, rng, k, [], 0.0, k + 0.05, float('inf'), (), tmp, None)
    out = dict(path=path, seed=seed, arm=arm, n=n, s_par=s, s_melt=s2, grow=s2 / s, st=r['status'], desc=desc,
               s_screen=r.get('s_screen'), s=r.get('s'), lines=r.get('lines'), sec=round(time.time() - t0, 2))
    if r.get('sq') is not None and r['status'] == 'new?' and not r.get('lines'):
        from pk.crossing import match_distance
        from pk.packing import Packing
        import numpy as np
        out['d'] = match_distance(Packing(r['s'], np.array(r['sq'])), Packing(s, np.array(sq)))
        out['ds'] = r['s'] - s
        out['self'] = abs(r['s'] - s) < 2e-9
    return out


def run(a):
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    tasks = []
    for n in (110, 130):
        for p in parents(n, a.per_n, rng):
            for _ in range(a.seeds):
                sd = rng.randrange(1 << 30)
                tasks.append([(p, sd, 'old'), (p, sd, 'new')])
    rng.shuffle(tasks)                                       # pairs stay together; partial runs stay balanced
    flat = [t for pair in tasks for t in pair]
    json.dump(dict(argv=sys.argv, pools=POOLS, tasks=len(flat)), open(f'{a.out}/meta.json', 'w'), indent=1)
    from pk import manifest
    man = manifest.start(a.out, 'melt_ab', seed=a.seed, args=vars(a), procs=a.procs, extra=dict(pools=POOLS))
    t_end = time.time() + 60 * a.minutes
    log = open(f'{a.out}/trials.jsonl', 'a')
    it = iter(flat); pend = set()
    with ProcessPoolExecutor(a.procs) as ex:
        while True:
            while len(pend) < a.procs and time.time() < t_end:
                j = next(it, None)
                if j is None: break
                pend.add(ex.submit(job, j))
            if not pend: break
            done, pend = wait(pend, return_when=FIRST_COMPLETED)
            for f in done:
                try:
                    log.write(json.dumps(f.result()) + '\n'); log.flush()
                except Exception as e:
                    log.write(json.dumps(dict(error=repr(e)[:200])) + '\n'); log.flush()
    man.finish()
    report(a.out)


def report(out):
    R = [json.loads(l) for l in open(f'{out}/trials.jsonl')]
    R = [r for r in R if 'arm' in r]
    key = lambda r: (r['path'], r['seed'])
    by = collections.defaultdict(dict)
    for r in R: by[key(r)][r['arm']] = r
    pairs = [v for v in by.values() if len(v) == 2]
    print(f'{out}: {len(R)} trials, {len(pairs)} complete pairs')
    import statistics as st
    for n in (110, 130):
        P = [p for p in pairs if p['old']['n'] == n]
        if not P: continue
        g = [p['new']['grow'] for p in P]
        print(f'\nn = {n}: {len(P)} pairs, {len(set(p["old"]["path"] for p in P))} parents; anneal box growth median {st.median(g):.4f} '
              f'(q10 {sorted(g)[len(g)//10]:.4f}, q90 {sorted(g)[9*len(g)//10]:.4f})')
        for arm in ('old', 'new'):
            X = [p[arm] for p in P]
            c = collections.Counter(x['st'] for x in X)
            ok = [x for x in X if 'd' in x]
            ret = sum(x['self'] for x in ok); away = [x for x in ok if not x['self']]
            below = sum(x['ds'] < -2e-9 for x in ok)
            d = sorted(x['d'] for x in away)
            ds = sorted(x['ds'] for x in ok)
            sec = st.mean(x['sec'] for x in X)
            print(f'  {arm}: status {dict(c)}; kept non-grid {len(ok)}/{len(X)}, parent-basin return {ret}, other basin {len(away)} '
                  f'(below parent {below}); d median {d[len(d)//2] if d else float("nan"):.4f} q90 {d[9*len(d)//10] if d else float("nan"):.4f}; '
                  f'ds median {ds[len(ds)//2] if ds else float("nan"):+.2e}; screen+polish {sec:.1f} s')
        same = sum(1 for p in P if p['old'].get('s') is not None and p['new'].get('s') is not None
                   and abs(p['old']['s'] - p['new']['s']) < 2e-9)
        both = sum(1 for p in P if 'd' in p['old'] and 'd' in p['new'])
        print(f'  paired: both kept {both}; same final side in both arms {same}/{len(P)}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='runs/melt_ab'); ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--minutes', type=float, default=12); ap.add_argument('--per-n', type=int, default=20)
    ap.add_argument('--seeds', type=int, default=8); ap.add_argument('--seed', type=int, default=1011)
    ap.add_argument('--report')
    a = ap.parse_args()
    if a.report:
        report(a.report)
    else:
        import scratch
        scratch.scope('melt_ab')
        run(a)
