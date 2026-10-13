#!/usr/bin/env python3
"""Soft-to-rigid positive control (PLAN-2026-10-12 item 1; 10-12).

Question: does fq's rounded-square pre-stage (--round-r0 ...) reproduce Ostashev's #575 band merges from the same seeds?
Seeds (store ids) and targets (#575):  241: Couzo 15.98813244 (id 347) -> 15.987134249051;
                                       211: de Winter 14.99796070 (id 94) -> 14.993694186584.
Success: any output below the seed side by > 1e-8 (secondary: below the #575 value).

  python3 round_ctl.py run --out runs/round_ctl/r1 --procs 14 [--ns 241,211] [--seeds 3]
  python3 round_ctl.py report runs/round_ctl/r1
"""
from __future__ import annotations
import argparse, itertools, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FQ = os.path.join(HERE, 'target', 'release', 'fq')
SEEDS = {241: (347, 15.98813243952499, 15.987134249051), 211: (94, 14.997960704967715, 14.993694186584)}


def export(ref, path):
    from fractions import Fraction
    subprocess.run([sys.executable, os.path.join(HERE, 'pk.py'), 'get', str(ref), '-o', path + '.tmp'], check=True,
                   capture_output=True)
    L = open(path + '.tmp').read().split('\n')
    n, s = L[0].split()
    L[0] = f'{n} {float(Fraction(s))!r}'
    open(path, 'w').write('\n'.join(L))
    os.remove(path + '.tmp')


def grid2(seeds):
    """Ostashev's README wording (10-11): soften slightly, hold the box at its side while re-hardening, release + compress
    (= our final quench).  No compression stage; larger kicks during the soft phase stand in for GPU noise / replicas."""
    for r0, steps, kick, seed in itertools.product([0.05, 0.1, 0.2], [10, 40], [0.03, 0.05, 0.1], range(seeds)):
        yield dict(r0=r0, comp=False, hold=1.0, steps=steps, p=1.0, kick=kick, seed=seed)


def grid(seeds, name='g1'):
    if name == 'g2':
        yield from grid2(seeds); return
    for r0, comp, hold, steps, pw, kick, seed in itertools.product(
            [0.03, 0.06, 0.1, 0.15, 0.2, 0.3], [True, False], [1.0, 0.999, 0.997], [10, 40], [1.0, 2.0],
            [0.0, 0.01], range(seeds)):
        if kick == 0.0 and seed > 0:
            continue
        if not comp and hold == 1.0 and kick == 0.0:
            continue          # nothing drives the configuration
        yield dict(r0=r0, comp=comp, hold=hold, steps=steps, p=pw, kick=kick, seed=seed)


def job(t):
    inp, outp, c = t
    a = [FQ, 'quench', '--in', inp, '--out', outp, '--loosen', '1.0', '--round-r0', str(c['r0']),
         '--hold-f', str(c['hold']), '--round-steps', str(c['steps']), '--round-p', str(c['p']),
         '--round-kick', str(c['kick']), '--seed', str(c['seed'] + 1)]
    if c['comp']:
        a.append('--round-compress')
    t0 = time.process_time()
    r = subprocess.run(a, capture_output=True, text=True, timeout=600)
    try:
        d = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return dict(c, err=r.stderr[-300:])
    return dict(c, s=d['s'], round=d.get('round'), min_gap=d['min_gap'], out=outp, sec=d['sec'])


def run(a):
    from pk import manifest
    import scratch
    os.makedirs(a.out, exist_ok=True)
    m = manifest.start(a.out, 'round_ctl', args=vars(a), procs=a.procs)
    ns = [int(x) for x in a.ns.split(',')]
    tasks = []
    for n in ns:
        sid, s0, tgt = SEEDS[n]
        inp = os.path.join(a.out, f'seed_{n}.txt')
        export(sid, inp)
        os.makedirs(os.path.join(a.out, f'n{n}'), exist_ok=True)
        for k, c in enumerate(grid(a.seeds, a.grid)):
            c = dict(c, n=n)
            tasks.append((inp, os.path.join(a.out, f'n{n}', f'{k:04d}.txt'), c))
    scratch.scope('round_ctl')
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'results.jsonl'), 'w') as f:
        for res in ex.map(job, tasks, chunksize=1):
            f.write(json.dumps(res) + '\n'); f.flush()
    m.finish(summary=dict(tasks=len(tasks)))


def report(a):
    R = [json.loads(l) for l in open(os.path.join(a.dir, 'results.jsonl'))]
    for n, (sid, s0, tgt) in SEEDS.items():
        rs = [r for r in R if r.get('n') == n and 's' in r]
        if not rs:
            continue
        rs.sort(key=lambda r: r['s'])
        below = [r for r in rs if r['s'] < s0 - 1e-8]
        print(f'n={n} seed {s0:.10f} target {tgt:.10f}: {len(rs)} runs, {len(below)} below seed, '
              f'{sum(r["s"] < tgt - 1e-9 for r in rs)} below target; best {rs[0]["s"]:.10f}')
        for r in rs[:8]:
            print(f'   {r["s"]:.10f} r0 {r["r0"]} comp {r["comp"]} hold {r["hold"]} steps {r["steps"]} p {r["p"]} '
                  f'kick {r["kick"]} seed {r["seed"]}  {r["out"]}')
        # by r0 / comp
        for key in ('r0', 'comp', 'hold'):
            vals = sorted({r[key] for r in rs}, key=str)
            print(f'   by {key}: ' + ', '.join(f'{v}: {sum(r["s"] < s0 - 1e-8 for r in rs if r[key] == v)}/'
                                             f'{sum(1 for r in rs if r[key] == v)}' for v in vals))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--out', required=True); r.add_argument('--procs', type=int, default=14)
    r.add_argument('--ns', default='241,211'); r.add_argument('--seeds', type=int, default=3)
    r.add_argument('--grid', default='g1')
    q = sp.add_parser('report'); q.add_argument('dir')
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else report(a)
