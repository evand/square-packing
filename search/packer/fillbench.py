#!/usr/bin/env python3
"""Refill benchmark for pk/fill.py (10-20).  A task = a good packing with some squares deleted; a fill config puts
squares back (placer + parameters, optional relax), the CS1 evaluate quenches the result (shove/harness.evaluate:
fq SCREEN at loosen 1.0 / 1.02 / 1.05 from the overlapping state, then POLISH, grid check).  The answer is known (the
source packing), so this measures fills, separately from any construction that feeds them.

Task kinds (per source packing = the store's best at n):
  uprights   delete every upright square, keep the tilted ones (exact upright fill restores these exactly at
             s + 1e-4 for 110 / 154 / 211 / 241 / 308: the easy case for exact, a hard one for random placers)
  disc:R     delete every square within R of a centre on the band boundary (a tilted square with an upright
             neighbour), so the hole mixes tilted and upright squares

Configs: "placer[:k=v,...][+schedule]", e.g. bestofK:K=64  leastbad:K=16+round.15c  exact  greedy:tau=0.03
  bestofK:K=16,stop=full  (stop=full: m=None, the placer stops when full, capped at m + 2).
Per row: count reached (it may differ from m: the result is scored at its own count), s' after evaluate, grid flag,
recovered (count = n and s' <= s0 + 1e-9), d_ref = s' - reference side at that count (the source at n, else the
store frontier), pen before the quench, CPU.

  python3 fillbench.py build --ns 110,154,211,241,308 --out runs/fillbench/b1 [--discs 1.5,2.5 --per 2]
  python3 fillbench.py run --bench runs/fillbench/b1 --configs "bestofK:K=64;exact" --seeds 2 --out runs/fillbench/b1/r1 [--procs 2 --cpu-h 1]
  python3 fillbench.py score runs/fillbench/b1/r1
"""
from __future__ import annotations
import argparse, json, math, os, random, resource, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, 'shove'))


def _tilted(sq, eps=1e-6):
    a = np.asarray(sq)[:, 2] % 90.0
    return (a > eps) & (a < 90.0 - eps)


def build(a):
    from pk.store import Store
    st = Store()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    tasks = []
    for n in [int(x) for x in a.ns.split(',')]:
        pid = st.query(n=n, limit=1)[0]['id']; p = st.get(pid)
        sq = np.asarray(p.sq, float); s0 = float(p.s)
        fr = st.frontier(n).get(n)
        src = dict(n=n, s0=s0, pid=int(pid), frontier=float(fr['s']) if fr else None)
        np.save(os.path.join(a.out, f'src_{n}.npy'), sq)
        T = _tilted(sq)
        tasks.append(dict(src, kind='uprights', delete=np.nonzero(~T)[0].tolist(), region=None))
        # band boundary: tilted squares with an upright neighbour within 1.6
        D = np.hypot(sq[:, None, 0] - sq[None, :, 0], sq[:, None, 1] - sq[None, :, 1])
        bnd = [i for i in np.nonzero(T)[0] if np.any((D[i] < 1.6) & ~T)]
        for R in [float(x) for x in a.discs.split(',')]:
            for _ in range(a.per if bnd else 0):
                c = sq[rng.choice(bnd), :2]
                dele = np.nonzero(np.hypot(sq[:, 0] - c[0], sq[:, 1] - c[1]) < R)[0].tolist()
                reg = [(float(c[0] - R), float(c[1] - R), float(c[0] + R), float(c[1] + R))]
                tasks.append(dict(src, kind=f'disc:{R:g}', centre=c.tolist(), delete=dele, region=reg))
    for i, t in enumerate(tasks):
        t['id'] = i; t['m'] = len(t['delete'])
    json.dump(dict(tasks=tasks, built=time.strftime('%Y-%m-%d %H:%M')), open(os.path.join(a.out, 'tasks.json'), 'w'), indent=0)
    for t in tasks:
        print(f"task {t['id']:3d}: n {t['n']} {t['kind']:9s} delete {t['m']:3d}  s0 {t['s0']:.9f}")


def parse_config(c):
    relax = None
    if '+' in c:
        c, relax = c.split('+', 1)
    name, _, kv = c.partition(':')
    P = {}
    for item in filter(None, kv.split(',')):
        k, v = item.split('=')
        P[k] = v if not v.replace('.', '', 1).replace('-', '', 1).isdigit() else (float(v) if '.' in v else int(v))
    return name, P, relax


def _cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def one(job):
    bench, task, config, seed = job
    import harness as H
    from pk import fill as F
    from pk.store import Store
    t0 = _cpu()
    sq = np.load(os.path.join(bench, f"src_{task['n']}.npy"))
    keep = np.setdiff1d(np.arange(len(sq)), task['delete'])
    S = task['s0'] + 1e-6
    name, P, relax = parse_config(config)
    stop = P.pop('stop', 'count')
    m = task['m']
    if stop == 'full':
        P.setdefault('m_max', m + 2)
    row = dict(task=task['id'], n=task['n'], kind=task['kind'], m=m, config=config, seed=seed)
    try:
        X, info = F.fill(S, sq[keep], random.Random(seed * 1000003 + task['id']), placer=name,
                         m=None if stop == 'full' else m, region=task['region'], relax_with=relax, **P)
        row.update(count=info['count'], added=info['added'], pen=info['pen'], s_relax=info['s'], params=info['params'])
        e = H.evaluate(S if info['s'] is None else info['s'], X)
        row.update(s2=e['s'], grid=e['grid'])
        if e['s'] is not None:
            c = info['count']
            ref = task['s0'] if c == task['n'] else (Store().frontier(c).get(c) or {}).get('s')
            row.update(ref=ref, d_ref=None if ref is None else e['s'] - ref,
                       recovered=bool(c == task['n'] and e['s'] <= task['s0'] + 1e-9))
    except Exception as ex:                                 # one failed fill never kills the run
        row['error'] = f'{type(ex).__name__}: {ex}'[:300]
    row['cpu'] = round(_cpu() - t0, 2)
    return row


def run(a):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from pk import manifest
    import scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'fillbench', args=vars(a), procs=a.procs)
    scratch.scope('fillbench')
    T = json.load(open(os.path.join(a.bench, 'tasks.json')))['tasks']
    if a.tasks:
        keep = {int(x) for x in a.tasks.split(',')}
        T = [t for t in T if t['id'] in keep]
    jobs = [(a.bench, t, c, sd) for sd in range(a.seeds) for t in T for c in a.configs.split(';')]
    print(f'{len(jobs)} jobs ({len(T)} tasks x {len(a.configs.split(";"))} configs x {a.seeds} seeds), {a.procs} procs', flush=True)
    cpu, nrow = 0.0, 0
    with open(os.path.join(a.out, 'results.jsonl'), 'a') as f, ProcessPoolExecutor(a.procs) as ex:
        futs = [ex.submit(one, j) for j in jobs]
        for fu in as_completed(futs):
            r = fu.result(); f.write(json.dumps(r) + '\n'); f.flush(); nrow += 1; cpu += r['cpu']
            if a.cpu_h and cpu > a.cpu_h * 3600:
                print('cpu budget reached', flush=True)
                ex.shutdown(wait=True, cancel_futures=True)      # CS1 harness note: don't keep running queued work
                break
    man.finish(summary=dict(rows=nrow, cpu_s=round(cpu, 1)))
    score([a.out])


def score(dirs):
    for d in dirs:
        R = [json.loads(l) for l in open(os.path.join(d, 'results.jsonl'))]
        print(f'{d}: {len(R)} rows, {sum("error" in r for r in R)} errors, CPU {sum(r["cpu"] for r in R) / 3600:.2f} h')
        by = {}
        for r in R:
            by.setdefault((r['config'], r['kind'].split(':')[0]), []).append(r)
        print(f"  {'config':32s} {'kind':9s} rows  count-n (min/med/max)  feasible  non-grid  recovered  d_ref<1e-3  med d_ref   med pen  CPU s")
        for (c, k), rs in sorted(by.items()):
            ok = [r for r in rs if r.get('s2') is not None]
            dc = sorted(r['count'] - r['n'] for r in rs if 'count' in r)
            dr = sorted(r['d_ref'] for r in ok if r.get('d_ref') is not None)
            pen = sorted(r['pen'] for r in rs if 'pen' in r)
            med = lambda v: v[len(v) // 2] if v else float('nan')
            print(f"  {c:32s} {k:9s} {len(rs):4d}  {(dc[0] if dc else 0):+d}/{med(dc):+d}/{(dc[-1] if dc else 0):+d}"
                  f"           {len(ok):4d}     {sum(not r['grid'] for r in ok):4d}      {sum(bool(r.get('recovered')) for r in ok):4d}"
                  f"       {sum(x < 1e-3 for x in dr):4d}   {med(dr):+.2e}  {med(pen):7.3f}  {med(sorted(r['cpu'] for r in rs)):6.1f}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('build'); b.add_argument('--ns', required=True); b.add_argument('--out', required=True)
    b.add_argument('--discs', default='1.5,2.5'); b.add_argument('--per', type=int, default=2); b.add_argument('--seed', type=int, default=0)
    r = sp.add_parser('run'); r.add_argument('--bench', required=True); r.add_argument('--out', required=True)
    r.add_argument('--configs', required=True, help='";"-separated config strings (see the module docstring)')
    r.add_argument('--seeds', type=int, default=1); r.add_argument('--tasks', help='comma list of task ids')
    r.add_argument('--procs', type=int, default=2); r.add_argument('--cpu-h', type=float, default=1.0)
    s = sp.add_parser('score'); s.add_argument('dirs', nargs='+')
    a = ap.parse_args()
    {'build': build, 'run': run, 'score': lambda a: score(a.dirs)}[a.cmd](a)
