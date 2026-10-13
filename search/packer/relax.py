#!/usr/bin/env python3
"""Relax compositions (or any packings) with fq's soft / shrink operators (10-13; Ostashev's README: "softening every
square of a composition can also help ... For 304, the composer's two-piece layout (17.939518) relaxed to 17.929428").

Inputs: compose2 near.jsonl rows (n, S_min, sq: a layout at its smallest fitting side) or packing files.  Each input x
operator x seed -> fq quench; the result side vs the input's frontier target.

  python3 relax.py run --near runs/relax/comp304.jsonl --target 17.929565946849 --out runs/relax/r304 [--ops ...]
  python3 relax.py run --ids 347 --out runs/relax/bm241 [--rounds R]        (store packings, e.g. two-band records)
"""
from __future__ import annotations
import argparse, itertools, json, math, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FQ = os.path.join(HERE, 'target', 'release', 'fq')
OPS = {
    'polish': ['--loosen', '1.0', '--no-alm'],                                   # rigid baseline
    'shrink': ['--loosen', '1.0', '--shrink'],                                   # fixed-side overlap relax + shrink hops
    'kshrink': ['--loosen', '1.0', '--kick', '0.02', '--shrink'],
    'soft05': ['--loosen', '1.0', '--round-r0', '0.05', '--hold-f', '1.0', '--round-steps', '40', '--round-kick', '0.03'],
    'soft15': ['--loosen', '1.0', '--round-r0', '0.15', '--hold-f', '1.0', '--round-steps', '40', '--round-kick', '0.03'],
    'soft15c': ['--loosen', '1.0', '--round-r0', '0.15', '--round-compress', '--round-steps', '40', '--round-kick', '0.02'],
    'loose': ['--loosen', '1.02', '--kick', '0.03'],                             # jitter + ALM from a dilated box
}


def job(t):
    inp, outp, op, seed = t
    r = subprocess.run([FQ, 'quench', '--in', inp, '--out', outp, '--seed', str(seed)] + OPS[op],
                       capture_output=True, text=True, timeout=1800)
    try:
        d = json.loads(r.stdout.strip().splitlines()[-1])
        return dict(op=op, seed=seed, s=d['s'], ok=d['min_gap'] > -1e-9 and d['min_wall'] > -1e-9, out=outp, sec=d.get('sec'))
    except Exception:
        return dict(op=op, seed=seed, err=r.stderr[-200:])


def run(a):
    from pk import manifest
    import scratch, mondrian
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'relax', args=vars(a), procs=a.procs)
    scratch.scope('relax')
    from pk import eventbench as E
    rows, per = [], {}
    near = [json.loads(l) for l in open(a.near)] if a.near else []
    if a.ids:                                       # store packings as inputs (band-merge relax of published packings)
        from pk.store import Store
        st = Store()
        for i in a.ids:
            p = st.get(i)
            near.append(dict(n=p.n, S_min=float(p.s), sq=np.asarray(p.sq).tolist(), need=0.0, A=f'store:{i}', B='-'))
    for r in sorted(near, key=lambda r: r['need']):
        if not (a.nmin <= r['n'] <= a.nmax and a.min_need < r['need'] < a.max_need) or \
                (a.ns and r['n'] not in {int(x) for x in a.ns.split(',')}):
            continue
        L = per.setdefault(r['n'], [])
        if len(L) < a.per_n and all((r['A'], r['B']) != x for x in L):
            L.append((r['A'], r['B'])); rows.append(r)
    print(f'{len(rows)} inputs over {len(per)} n', flush=True)
    meta, cur = {}, {}
    for i, r in enumerate(rows):
        X = mondrian.trim(np.array(r['sq'], float), r['n'])
        inp = os.path.join(a.out, f'in{i:03d}.txt')
        with open(inp, 'w') as f:
            f.write(f"{r['n']} {float(r['S_min'])!r}\n" + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in X))
        meta[i] = dict(n=r['n'], S0=r['S_min'], A=r.get('A'), B=r.get('B'),
                       target=a.target if a.target else E.frontier_now(r['n']))
        cur[i] = (inp, float(r['S_min']))
    # rounds: monotonic basin hopping with soft kicks; each input restarts from its best result so far
    tasks = []
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'results.jsonl'), 'w') as f:
        for rd in range(a.rounds):
            T = [(cur[i][0], os.path.join(a.out, f'out{i:03d}_r{rd}_{op}_{seed}.txt'), op, seed + 1000 * rd, i)
                 for i in cur for op, seed in itertools.product(a.ops.split(','), range(1, a.seeds + 1))
                 if not (op == 'polish' and seed > 1)]
            tasks += T
            for t, r in zip(T, ex.map(job, [t[:4] for t in T])):
                r.update(meta[t[4]], i=t[4], round=rd, s_in=cur[t[4]][1]); f.write(json.dumps(r) + '\n'); f.flush()
                if r.get('ok') and r['s'] < cur[t[4]][1] - 1e-12:
                    cur[t[4]] = (r['out'], r['s'])
            if a.rounds > 1:
                print(f'round {rd}: ' + ' '.join(f"{meta[i]['n']}:{cur[i][1] - meta[i]['target']:+.2e}" for i in cur),
                      flush=True)
    man.finish(summary=dict(tasks=len(tasks)))
    report(a)


def report(a):
    R = [json.loads(l) for l in open(os.path.join(a.out, 'results.jsonl'))]
    R = [dict(r, target=r.get('target', a.target)) for r in R if r.get('ok')]
    print(f'{len(R)} ok runs')
    for r in sorted((r for r in R if r['s'] < r['target'] - 1e-9), key=lambda r: r['n']):
        print(f"  BELOW n={r['n']} s {r['s']:.12f} target {r['target']:.12f} gain {r['target'] - r['s']:.2e} {r['op']} {r['out']}")
    for op in a.ops.split(','):
        rs = [r for r in R if r['op'] == op]
        if not rs:
            continue
        d = [r['S0'] - r['s'] for r in rs]
        b = min(rs, key=lambda r: r['s'])
        print(f'{op:8s} {len(rs):4d} runs  best {b["s"]:.10f} (in {b["i"]})  drop: max {max(d):.2e} median {np.median(d):.2e}'
              f'  >5e-3: {sum(x > 5e-3 for x in d)}  below target: {sum(r["s"] < r["target"] - 1e-9 for r in rs)}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--near')
    r.add_argument('--ids', type=int, nargs='*', help='store packing ids as inputs (instead of / besides --near)'); r.add_argument('--out', required=True)
    r.add_argument('--target', type=float, default=0.0, help='default: frontier_now(n) per input')
    r.add_argument('--ops', default=','.join(OPS))
    r.add_argument('--nmin', type=int, default=0); r.add_argument('--nmax', type=int, default=10 ** 6)
    r.add_argument('--min-need', type=float, default=-1.0); r.add_argument('--max-need', type=float, default=1.0)
    r.add_argument('--per-n', type=int, default=10 ** 6); r.add_argument('--rounds', type=int, default=1)
    r.add_argument('--ns', help='comma list: only these n')
    r.add_argument('--seeds', type=int, default=4); r.add_argument('--procs', type=int, default=14)
    q = sp.add_parser('report'); q.add_argument('--out', required=True); q.add_argument('--target', type=float, default=0.0)
    q.add_argument('--ops', default=','.join(OPS))
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else report(a)
