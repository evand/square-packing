#!/usr/bin/env python3
"""Soft-to-rigid re-hardening driver (10-14; engine `src/harden.rs`, Nakajima's protocol as Ostashev #575 uses it).

Each input (store packing or file) x transform x parameter cell x seed -> `harden reharden` -> fq polish
(`quench --loosen 1.0 --no-alm`) -> side vs the target (default: frontier_now(n') for the output count n').

  python3 harden.py run --ids 347 --out runs/harden/bm241 --grid "r0=0.02,0.05,0.1;hold=8000,30000" --seeds 8
  python3 harden.py run --ids 4711 --xform insert:1 --hold-f ... --out ...      (n -> n + 1: the border / growth test)
  python3 harden.py run --ids 4711 --xform remove:all --out ...                 (n + 1 -> n, every square: plateaus)
  python3 harden.py report --out DIR

Transforms: `none`; `insert:K` (K squares at uniform random positions / angles in the box, per seed); `remove:all`
(each square in turn) or `remove:R` (R random single removals, per seed).  Grid keys are harden flags without the dashes
(r0, pre, hold, release, settle, hold_f, noise, mu, lr, p); a cell is one combination.
"""
from __future__ import annotations
import argparse, itertools, json, math, os, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
HARDEN = os.path.join(HERE, 'target', 'release', 'harden')
FQ = os.path.join(HERE, 'target', 'release', 'fq')


def parse_grid(g):
    keys, vals = [], []
    for part in filter(None, (g or '').split(';')):
        k, v = part.split('=')
        keys.append(k.strip()); vals.append([x.strip() for x in v.split(',')])
    return [dict(zip(keys, c)) for c in itertools.product(*vals)] if keys else [{}]


def write_cfg(path, n, s, sq):
    with open(path, 'w') as f:
        f.write(f'{n} {float(s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in sq))


def make_input(sq, s, xform, seed, k):
    """(n', side, squares) after the transform; k = index for remove:all."""
    sq = np.asarray(sq, float)
    rng = np.random.default_rng(seed * 7919 + k)
    kind, _, arg = xform.partition(':')
    if kind == 'none':
        return len(sq), s, sq
    if kind == 'insert':
        m = int(arg or 1)
        add = np.column_stack([rng.uniform(0.5, s - 0.5, m), rng.uniform(0.5, s - 0.5, m), rng.uniform(0, 90, m)])
        out = np.vstack([sq, add])
        return len(out), s, out
    if kind == 'remove':
        i = k if arg == 'all' else int(rng.integers(len(sq)))
        out = np.delete(sq, i, axis=0)
        return len(out), s, out
    raise ValueError(xform)


def load_out(path):
    L = open(path).read().split('\n')
    n, s = L[0].split()[:2]
    sq = np.array([[float(v) for v in l.split()[:3]] for l in L[1:1 + int(n)]])
    return f'file:{path}', int(n), float(s), sq


def job(t):
    inp, raw, pol, cell, seed = t
    cmd = [HARDEN, 'reharden', '--in', inp, '--out', raw, '--seed', str(seed)]
    for k, v in cell.items():
        cmd += ['--' + k.replace('_', '-'), v]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    try:
        h = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        return dict(err='harden: ' + r.stderr[-300:])
    q = subprocess.run([FQ, 'quench', '--in', raw, '--out', pol, '--seed', str(seed), '--loosen', '1.0', '--no-alm'],
                       capture_output=True, text=True, timeout=3600)
    try:
        d = json.loads(q.stdout.strip().splitlines()[-1])
    except Exception:
        return dict(s_raw=h['s'], h_sec=h['sec'], err='fq: ' + q.stderr[-300:])
    return dict(s_raw=h['s'], h_sec=h['sec'], s=d['s'], ok=d['min_gap'] > -1e-9 and d['min_wall'] > -1e-9,
                q_sec=d.get('sec'), out=pol)


def run(a):
    from pk import manifest
    import scratch
    from pk.store import Store
    from pk import eventbench as E
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'harden', args=vars(a), procs=a.procs)
    scratch.scope('harden')
    st = Store()
    srcs = []
    for i in a.ids or []:
        p = st.get(i)
        srcs.append((f'store:{i}', p.n, float(p.s), np.asarray(p.sq, float)))
    for fpath in a.files or []:
        L = open(fpath).read().split('\n')
        n, s = L[0].split()[:2]
        s = float(eval(s)) if '/' in s else float(s)
        sq = np.array([[float(eval(v)) if '/' in v else float(v) for v in l.split()[:3]] for l in L[1:1 + int(n)]])
        srcs.append((f'file:{fpath}', int(n), s, sq))
    cells = parse_grid(a.grid)
    targets, best, ntask = {}, {}, 0
    # round 0: sources x transform; later rounds: the --keep best distinct outputs per source (transform none)
    cur = {si: [(lab, n, s, sq, a.xform)] for si, (lab, n, s, sq) in enumerate(srcs)}
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'results.jsonl'), 'a') as f:
        for rd in range(a.rounds):
            tasks, meta = [], []
            for si, ins in cur.items():
                for ii, (lab, n, s, sq, xf) in enumerate(ins):
                    ks = range(n) if xf == 'remove:all' else range(1)
                    for k in ks:
                        for seed in range(1, a.seeds + 1):
                            n2, s2, sq2 = make_input(sq, s, xf, seed, k)
                            inp = os.path.join(a.out, f'in{si}_r{rd}_{ii}_{k}_{seed}.txt')
                            write_cfg(inp, n2, s2, sq2)
                            if n2 not in targets:
                                targets[n2] = a.target if a.target else E.frontier_now(n2)
                            for ci, cell in enumerate(cells):
                                tag = f'{si}_r{rd}_{ii}_{k}_{seed}_c{ci}'
                                tasks.append((inp, os.path.join(a.out, f'raw{tag}.txt'), os.path.join(a.out, f'pol{tag}.txt'),
                                              cell, seed + 1000 * k + 100000 * rd))
                                meta.append(dict(src=lab, si=si, round=rd, n_in=n, s_in=s, n=n2, k=k, seed=seed, cell=ci,
                                                 **{'p_' + x: y for x, y in cell.items()}))
            print(f'round {rd}: {len(cur)} sources, {len(cells)} cells, {len(tasks)} runs; targets ' +
                  ' '.join(f'{n}:{t:.10f}' for n, t in sorted(targets.items())), flush=True)
            ntask += len(tasks)
            outs = {}
            for m, r in zip(meta, ex.map(job, tasks, chunksize=1)):
                r.update(m, target=targets[m['n']])
                f.write(json.dumps(r) + '\n'); f.flush()
                if not r.get('ok'):
                    continue
                outs.setdefault(m['si'], []).append((r['s'], r['out'], m['n']))
                if r['s'] < best.get(m['n'], 1e9):
                    best[m['n']] = r['s']
                    gap = r['s'] - targets[m['n']]
                    print(f"  new best n={m['n']} s={r['s']:.10f} ({gap:+.2e} vs target) round {rd} cell {m['cell']} seed {m['seed']}"
                          + ('  ** BELOW TARGET **' if gap < -1e-8 else ''), flush=True)
            # next round's inputs: previous inputs and outputs pooled, keep the best distinct (side 1e-7 apart)
            nxt = {}
            for si in cur:
                pool = [(s, None, (lab, n, s, sq)) for (lab, n, s, sq, xf) in cur[si] if xf == 'none']
                pool += [(s_, o, None) for s_, o, n_ in outs.get(si, [])]
                pool.sort(key=lambda t: t[0])
                keep = []
                for s_, o, src in pool:
                    if all(abs(s_ - k[0]) > 1e-7 for k in keep):
                        keep.append((s_, o, src))
                    if len(keep) == a.keep:
                        break
                L = []
                for s_, o, src in keep:
                    if src is not None:
                        L.append(src + ('none',))
                    else:
                        lab, n_, s2, sq2 = load_out(o)
                        L.append((lab, n_, s2, sq2, 'none'))
                if L:
                    nxt[si] = L
            cur = nxt
            if a.rounds > 1:
                print(f'round {rd} done: ' + ' '.join(f'{n}:{v - targets[n]:+.2e}' for n, v in sorted(best.items())), flush=True)
    man.finish(summary=dict(runs=ntask, best={str(k): v for k, v in best.items()}))
    report(a)


def report(a):
    R = [json.loads(l) for l in open(os.path.join(a.out, 'results.jsonl'))]
    bad = [r for r in R if not r.get('ok')]
    R = [r for r in R if r.get('ok')]
    print(f'{len(R)} ok runs, {len(bad)} failed')
    for r in sorted((r for r in R if r['s'] < r['target'] - 1e-8), key=lambda r: r['s']):
        print(f"  BELOW n={r['n']} s {r['s']:.12f} target {r['target']:.12f} gain {r['target'] - r['s']:.2e} cell {r['cell']} {r['out']}")
    cells = sorted({r['cell'] for r in R})
    for c in cells:
        rs = [r for r in R if r['cell'] == c]
        par = {k[2:]: v for k, v in rs[0].items() if k.startswith('p_')}
        g = np.array([r['s'] - r['s_in'] for r in rs if r['n'] == r['n_in']] or [np.nan])
        b = min(rs, key=lambda r: r['s'] - r['target'])
        print(f"cell {c:3d} {par}  runs {len(rs):4d}  best vs target {b['s'] - b['target']:+.2e} (n={b['n']})"
              f"  vs input: min {np.nanmin(g):+.2e} median {np.nanmedian(g):+.2e}"
              f"  sec {np.median([r['h_sec'] + (r['q_sec'] or 0) for r in rs]):.1f}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--ids', type=int, nargs='*'); r.add_argument('--files', nargs='*')
    r.add_argument('--out', required=True)
    r.add_argument('--xform', default='none')
    r.add_argument('--grid', default='')
    r.add_argument('--seeds', type=int, default=4)
    r.add_argument('--target', type=float, default=0.0, help='default: frontier_now(n) of the output count')
    r.add_argument('--procs', type=int, default=4)
    r.add_argument('--rounds', type=int, default=1, help='monotone rounds: each source restarts from its --keep best')
    r.add_argument('--keep', type=int, default=4)
    q = sp.add_parser('report'); q.add_argument('--out', required=True)
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else report(a)
