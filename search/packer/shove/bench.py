#!/usr/bin/env python3
"""Corner-shove bench (10-16).  One task = (source packing at n as of cutoff T, corner whose block is shoved);
each task runs at m = 2k, 2k + 1, 2k + 2 added squares (k = floor s: the L-step adds 2k + 1 at side s + 1).

  python3 shove/bench.py build --cutoff reg1005 --out runs/shove/b1005 [--nmin 100 --nmax 300 --min-block 4]
  python3 shove/bench.py run --bench runs/shove/b1005 --fill greedy --out runs/shove/b1005/greedy
         [--d 1.0] [--limit N] [--tasks-every K] [--procs 4] [--seed 0] [--cpu-h 4]
  python3 shove/bench.py score runs/shove/b1005/greedy [more run dirs ...]

--fill NAME loads shove/fill_NAME.py (see harness.py for the interface).  Built in: 'lstep' (Friedman's L, no shove:
the reference), 'uprights' (shove + exact upright refill of the whole box around the tilted squares: transplant.py's
operator), 'none' (shove, no added squares: plumbing check).

Per evaluated candidate (results.jsonl): n, corner, block size, m, n' (its count), s' (quenched side), grid flag,
vsL = s' - (s + 1) (< 0: beats the L-step from this source), vsT = s' - frontier_at(n', T) (< 0: below the record
as of the cutoff), vsNow = s' - frontier_now(n'), rediscovered = vsNow <= 1e-9 < vsT (a later record re-found).
Every s' < frontier_now - 1e-9 is printed as a LIVE candidate (run candidate.py; frontier_now is the event-bench
snapshot, so live-check it).
"""
import argparse, glob, importlib, json, math, os, random, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK); sys.path.insert(0, HERE)
import harness as H


def read_pk(path):
    with open(path) as f:
        n, s = f.readline().split()[:2]
        sq = np.array([[float(v.replace('np.float64(', '').rstrip(')')) for v in l.split()[:3]] for l in f if l.strip()])
    from fractions import Fraction
    return float(Fraction(s)), sq[:int(n)]


def build(a):
    from pk import eventbench as EB
    os.makedirs(a.out, exist_ok=True)
    EB.starts(a.cutoff, os.path.join(a.out, 'starts'))
    t = EB.cutoff(a.cutoff)
    tasks = []
    for f in sorted(glob.glob(os.path.join(a.out, 'starts', 'n-*.txt')), key=lambda p: int(p.split('-')[-1][:-4])):
        n = int(f.split('-')[-1][:-4])
        if not (a.nmin <= n <= a.nmax):
            continue
        s, sq = read_pk(f)
        if abs(s - round(s)) < 1e-9:
            continue                                       # grid packings: nothing to shove
        k = int(math.floor(s + 1e-9))
        B = H.corner_blocks(s, sq)
        for c, blk in B.items():
            if len(blk) < a.min_block:
                continue
            tg = {}
            for m in (2 * k, 2 * k + 1, 2 * k + 2):
                n2 = n + m
                tg[str(m)] = dict(n2=n2, fT=EB.frontier_at(n2, t), fNow=EB.frontier_now(n2))
            tasks.append(dict(n=n, s=s, k=k, corner=c, block=len(blk), src=os.path.relpath(f, a.out), targets=tg))
    json.dump(dict(cutoff=a.cutoff, t=t.isoformat(), tasks=tasks), open(os.path.join(a.out, 'tasks.json'), 'w'), indent=0)
    print(f'{len(tasks)} tasks over {len({x["n"] for x in tasks})} n')


def get_fill(name):
    if name in ('lstep', 'uprights', 'none'):
        return None
    return importlib.import_module(f'fill_{name}').fill


def one(args):
    bench, task, fill_name, d, seed, max_cands = args
    T = json.load(open(os.path.join(bench, 'tasks.json')))
    from pk import eventbench as EB
    s, sq = read_pk(os.path.join(bench, task['src']))
    rng = random.Random(seed)
    import resource
    def _cpu():
        a_, b_ = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
        return a_.ru_utime + a_.ru_stime + b_.ru_utime + b_.ru_stime
    t0 = _cpu()
    rows = []
    cands = []                                             # (m, array, S)
    if fill_name == 'lstep':
        S, X = H.lstep(s, sq); cands.append((2 * task['k'] + 1, X, S))
    else:
        if task.get('ishove'):                             # 10-17 IS1: island shove (shove/island.py)
            import island
            S, sq2, meta = island.from_task(s, sq, task, d)
        else:
            S, sq2, meta = H.shove(s, sq, task['corner'], d)
        if fill_name == 'none':
            cands.append((0, sq2, S))
        elif fill_name == 'uprights':
            import mondrian
            ax = (np.minimum(sq2[:, 2] % 90, 90 - sq2[:, 2] % 90) < 1.0)
            F = mondrian.exact_fill([tuple(q) for q in sq2[~ax]], S)
            X = np.asarray(F, float)
            cands.append((len(X) - len(sq2), X, S))
        else:
            fill = get_fill(fill_name)
            for m in task.get('ms') or (2 * task['k'], 2 * task['k'] + 1, 2 * task['k'] + 2):   # 'ms': per-task counts
                try:
                    out = fill(S, sq2.copy(), dict(meta), m, rng)
                except Exception as ex:
                    rows.append(dict(n=task['n'], corner=task['corner'], m=m, error=repr(ex)[:200]))
                    continue
                for X in out[:max_cands]:
                    cands.append((m, np.asarray(X, float), S))
    for ci, (m, X, S) in enumerate(cands):
        try:
            e = H.evaluate(S, X)
        except Exception as ex:                            # 10-17: one bad candidate must not kill the bench
            rows.append(dict(n=task['n'], corner=task['corner'], m=m, error='evaluate: ' + repr(ex)[:200]))
            continue
        n2 = len(X)
        tg = task['targets'].get(str(n2 - task['n']))
        fT = tg['fT'] if tg else EB.frontier_at(n2, EB.ts(T['t']))
        fN = tg['fNow'] if tg else EB.frontier_now(n2)
        r = dict(n=task['n'], s0=s, corner=task['corner'], block=task['block'], m=m, n2=n2, S=S, s2=e['s'], grid=e['grid'],
                 vsL=None if e['s'] is None else e['s'] - (s + 1), vsT=None if e['s'] is None else e['s'] - fT,
                 vsNow=None if e['s'] is None else e['s'] - fN)
        r.update(task.get('tags') or {})                   # 10-17: per-task features (parent gap etc.) for parent-value fits
        r['rediscovered'] = bool(e['s'] is not None and r['vsNow'] <= 1e-9 < r['vsT'])
        # 10-17 (agent A): a real win over the L-step = more squares than it at a side below s + 1, or two more at s + 1
        k0 = task['k']
        r['beatL'] = bool(e['s'] is not None and not e['grid'] and
                          ((n2 - task['n'] >= 2 * k0 + 1 and r['vsL'] < -1e-9) or (n2 - task['n'] >= 2 * k0 + 2 and r['vsL'] <= 1e-9)))
        r['live'] = bool(e['s'] is not None and r['vsNow'] < -1e-9)
        if e['s'] is not None and (r['live'] or r['beatL'] or r['rediscovered'] or r['vsT'] < 1e-3):
            r['file'] = f"{task['n']}_{task['corner']}_{n2}_{seed}_c{ci}" + (f"_t{task['id']}" if 'id' in task else '') + '.txt'
            r['_sq'] = e['sq']
        rows.append(r)
    for r in rows:
        r['cpu'] = round((_cpu() - t0) / max(1, len(rows)), 2)
    return rows


def run(a):
    from concurrent.futures import ProcessPoolExecutor
    import mcmin
    from pk import manifest
    os.makedirs(a.out, exist_ok=True); os.makedirs(os.path.join(a.out, 'hits'), exist_ok=True)
    man = manifest.start(a.out, f'shove_bench_{a.fill}', args=vars(a), procs=a.procs)
    T = json.load(open(os.path.join(a.bench, 'tasks.json')))['tasks'][::a.tasks_every]
    if a.limit:
        T = T[:a.limit]
    frm = a.from_task or 0                                 # 10-17: shards keep each task's seed (global index)
    T = T[frm:a.to_task]
    print(f'{len(T)} tasks, fill {a.fill}, d {a.d}', flush=True)
    t0 = time.time(); cpu = 0.0; nrow = 0
    with open(os.path.join(a.out, 'results.jsonl'), 'a') as out, ProcessPoolExecutor(a.procs) as ex:
        jobs = [(a.bench, t, a.fill, a.d, a.seed * 100003 + frm + i, a.max_cands) for i, t in enumerate(T)]
        for rows in ex.map(one, jobs):
            for r in rows:
                sqv = r.pop('_sq', None)
                if sqv is not None:
                    mcmin.write_deg(os.path.join(a.out, 'hits', r['file']), r['s2'], [tuple(map(float, q)) for q in sqv])
                out.write(json.dumps(r) + '\n'); nrow += 1
                cpu += r.get('cpu', 0)
                if r.get('live'):
                    print(f"** LIVE? n'={r['n2']} s'={r['s2']:.10f} vsNow {r['vsNow']:+.2e} ({r['file']}): live-check + candidate.py", flush=True)
            out.flush()
            if a.cpu_h and cpu > a.cpu_h * 3600:
                print('cpu budget reached', flush=True); break
    man.finish(summary=dict(rows=nrow))
    score_dirs([a.out])


def score_dirs(dirs):
    for d in dirs:
        R = [json.loads(l) for l in open(os.path.join(d, 'results.jsonl'))]
        ok = [r for r in R if r.get('s2') is not None]
        err = [r for r in R if 'error' in r]
        ng = [r for r in ok if not r['grid']]
        cpu = sum(r.get('cpu', 0) for r in R) / 3600
        bL = [r for r in ng if r.get('beatL', r['vsL'] < -1e-9 and r['n2'] - r['n'] > 2 * int(r['s0']))]
        bT = [r for r in ng if r['vsT'] < -1e-9]
        red = [r for r in ng if r['rediscovered']]
        live = [r for r in ng if r['live']]
        print(f"{d}: {len(R)} rows ({len(err)} fill errors), {len(ok)} feasible, {len(ng)} non-grid; CPU {cpu:.2f} h")
        print(f"   beat L-step from its source (more squares below s + 1, or 2 more at s + 1): {len(bL)} on {len({(r['n'], r['corner']) for r in bL})} tasks   below record as of cutoff: {len(bT)}   rediscovered later "
              f"records: {len(red)}   below event-bench now: {len(live)}")
        if ng:
            v = sorted(r['vsL'] for r in ng)
            print(f"   vsL quantiles (s' - (s + 1)): min {v[0]:+.3e}  10% {v[len(v) // 10]:+.3e}  median {v[len(v) // 2]:+.3e}")
            for r in sorted(ng, key=lambda r: r['vsT'])[:8]:
                print(f"     n {r['n']} {r['corner']} blk {r['block']} m {r['m']} -> n' {r['n2']} s' {r['s2']:.10f}  vsL {r['vsL']:+.2e}"
                      f"  vsT {r['vsT']:+.2e}  vsNow {r['vsNow']:+.2e}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('build'); b.add_argument('--cutoff', default='reg1005'); b.add_argument('--out', required=True)
    b.add_argument('--nmin', type=int, default=100); b.add_argument('--nmax', type=int, default=300)
    b.add_argument('--min-block', type=int, default=4)
    r = sp.add_parser('run'); r.add_argument('--bench', required=True); r.add_argument('--fill', required=True)
    r.add_argument('--out', required=True); r.add_argument('--d', type=float, default=1.0)
    r.add_argument('--limit', type=int); r.add_argument('--tasks-every', type=int, default=1)
    r.add_argument('--procs', type=int, default=4); r.add_argument('--seed', type=int, default=0)
    r.add_argument('--max-cands', type=int, default=3)
    r.add_argument('--from-task', type=int); r.add_argument('--to-task', type=int); r.add_argument('--cpu-h', type=float, default=4.0)
    s = sp.add_parser('score'); s.add_argument('dirs', nargs='+')
    a = ap.parse_args()
    {'build': build, 'run': run, 'score': lambda a: score_dirs(a.dirs)}[a.cmd](a)
