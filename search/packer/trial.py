#!/usr/bin/env python3
"""Trials (10-22): grade an idea (a chain head, pk/stages.py) on a frozen suite, after the suite's standard tail, against
the suite's baseline head given the same tail.  The raw head output is logged but never the verdict (CLAUDE.md rule:
a generator's output is a distribution of starts for refinement).

  trial.py run NAME ARM --head SPEC --suite SUITE [--starts 8] [--procs 4]
  trial.py run NAME baseline --suite SUITE                  the suite's baseline head (run it once per trial)
  trial.py report NAME

NAME must have a "### Pre-registration: ... (NAME" block in PACKER.md or a plugin NOTEBOOK (hypothesis, metric, rule)
before it runs.
Output: runs/trials/NAME/ARM.jsonl, one row per start: n, seed, head / tail / chain specs, the stage log, raw side
(after the head), final side, best feasible side along the chain, measured CPU.  The report compares each arm with the
baseline per n: fraction of starts whose best side beats the baseline's median, best side, CPU, raw side quantiles.

TODO (10-22, TODO.md): sequential stopping and rate posteriors (benchkit.race / battery statistics); more suites.
"""
from __future__ import annotations
import argparse, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, 'runs', 'trials')

# The standard tail: release, side-free rounded relax, shallow anneal, polish.  Change it only with a new suite name
# (trials on one suite must stay comparable).  PROVISIONAL (10-22 smoke test, n = 548): one perturb with no acceptance
# only worsens a good start (composer 23.89211 -> anneal 24.83 -> polish 23.97), so `best` is the raw side.  Replace it
# by an accept-loop tail (basin hopping: repeat(perturb | relax) with Metropolis acceptance) before grading anything.
TAIL = 'constrain.release | relax.round.15c | perturb.anneal(sweeps=1000) | relax.polish'

# Suites: frozen (n list, tail, baseline head, reference sides).  None are public yet; plugins add theirs as TRIAL_SUITES
# (pk/stages.py PK_PLUGINS; a suite with tail None gets TAIL).  Fields: ns, tail, baseline (a chain head),
# ref = (json with per-n 'ours' / 'pub' sides, description).
SUITES = {}


def _plugin_suites():
    from pk import stages as S
    for m in S.PLUGINS:
        for k, v in getattr(m, 'TRIAL_SUITES', {}).items():
            SUITES[k] = dict(v, tail=v.get('tail') or TAIL)


def ref_sides(suite):
    path, _ = SUITES[suite]['ref']
    c = json.load(open(os.path.join(HERE, path)))
    return {int(k): min(c['ours'][k], c['pub'][k]) for k in c['ours']}


def preregistered(name):
    """a '### Pre-registration ... (NAME' heading in PACKER.md or a plugin's notebook (plugin NOTEBOOK = path)"""
    from pk import stages as S
    paths = [os.path.join(HERE, 'PACKER.md')] + [m.NOTEBOOK for m in S.PLUGINS if hasattr(m, 'NOTEBOOK')]
    txt = '\n'.join(open(p).read() for p in paths)
    return any(l.startswith('### Pre-registration') and f'({name}' in l for l in txt.splitlines())


def one(args):
    n, seed, head, tail = args
    from pk import stages as S
    t0 = S._cpu()
    st, log, err = None, [], None
    raw = None
    try:
        st, log = S.run_chain(head, dict(n=n), seed)
        raw = float(st.p.s) if st.p.n == n and st.p.max_pen() <= 1e-9 else None
        st, log2 = S.run_chain(tail, dict(n=n), seed + 7919, st=st)
        log += log2
    except Exception as e:
        err = repr(e)[:300]
    feas = [r['s'] for r in log if r['count'] == n and r['frozen'] == 0 and r['kind'] == 'relax' and not r['info'].get('failed')]
    final = float(st.p.s) if st is not None and st.p.n == n and st.p.max_pen() <= 1e-9 else None
    best = min([x for x in feas + [raw, final] if x is not None], default=None)
    return dict(n=n, seed=seed, head=head, tail=tail, raw=raw, final=final, best=best, cpu=round(S._cpu() - t0, 2),
                err=err, log=log), (st.p if st is not None and best is not None and final == best else None)


def cmd_run(a):
    if not preregistered(a.name):
        sys.exit(f'{a.name}: no "### Pre-registration: ... ({a.name}" block in PACKER.md; write hypothesis, metric and rule first')
    su = SUITES[a.suite]
    head = su['baseline'] if a.arm == 'baseline' else a.head
    if not head:
        sys.exit('--head SPEC needed (or arm "baseline")')
    from pk import stages as S
    S.parse(head); S.parse(su['tail'])                     # fail fast on a bad spec
    from pk import manifest
    d = os.path.join(OUT, a.name); os.makedirs(d, exist_ok=True)
    man = manifest.start(d, f'trial:{a.name}:{a.arm}', args=vars(a), procs=a.procs, fname=f'manifest_{a.arm}.json')
    starts = 1 if a.arm == 'baseline' and 'gen.file' in head and 'perturb' not in head else a.starts
    tasks = [(n, 1000 * n + k + a.seed0, head, su['tail']) for k in range(starts) for n in su['ns']]
    pdir = os.path.join(d, a.arm); os.makedirs(pdir, exist_ok=True)
    rows = 0
    with open(os.path.join(d, f'{a.arm}.jsonl'), 'a') as fh, ProcessPoolExecutor(a.procs) as ex:
        for rec, p in ex.map(one, tasks):
            rec['suite'] = a.suite; rec['arm'] = a.arm
            fh.write(json.dumps(rec) + '\n'); fh.flush(); rows += 1
            if p is not None:
                p.write(os.path.join(pdir, f'n{rec["n"]}_s{rec["seed"]}.txt'))
    man.finish(dict(rows=rows))
    print(f'{rows} starts -> {d}/{a.arm}.jsonl')


def cmd_report(a):
    import numpy as np
    d = os.path.join(OUT, a.name)
    arms = {}
    for f in sorted(os.listdir(d)):
        if f.endswith('.jsonl'):
            arms[f[:-6]] = [json.loads(l) for l in open(os.path.join(d, f))]
    if 'baseline' not in arms:
        sys.exit('run the baseline arm first')
    suite = arms['baseline'][0]['suite']; ref = ref_sides(suite)
    base = {}
    for r in arms['baseline']:
        if r['best'] is not None:
            base.setdefault(r['n'], []).append(r['best'])
    bmed = {n: float(np.median(v)) for n, v in base.items()}
    print(f'suite {suite}; tail: {arms["baseline"][0]["tail"]}')
    print(f'{"arm":14s} {"n":>4s} {"starts":>6s} {"ok":>3s} {"beat base":>9s} {"best - base":>12s} {"best - ref":>11s} '
          f'{"raw med - base":>14s} {"CPU s":>7s} errors')
    for arm, R in arms.items():
        for n in sorted({r['n'] for r in R}):
            rr = [r for r in R if r['n'] == n]
            b = [r['best'] for r in rr if r['best'] is not None]
            raw = [r['raw'] for r in rr if r['raw'] is not None]
            bm = bmed.get(n, math.inf)
            beat = sum(x < bm - 1e-9 for x in b)
            fmt = lambda x: f'{x:+.5f}' if x is not None and math.isfinite(x) else '-'
            print(f'{arm:14s} {n:4d} {len(rr):6d} {len(b):3d} {beat:4d}/{len(rr):<4d} {fmt(min(b) - bm if b else None):>12s} '
                  f'{fmt(min(b) - ref[n] if b else None):>11s} {fmt(float(np.median(raw)) - bm if raw else None):>14s} '
                  f'{sum(r["cpu"] for r in rr):7.0f} {sum(r["err"] is not None for r in rr)}')


def main():
    _plugin_suites()
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run'); r.add_argument('name'); r.add_argument('arm'); r.add_argument('--head')
    r.add_argument('--suite', required=True, choices=sorted(SUITES)); r.add_argument('--starts', type=int, default=8)
    r.add_argument('--procs', type=int, default=4); r.add_argument('--seed0', type=int, default=0)
    p = sub.add_parser('report'); p.add_argument('name')
    a = ap.parse_args()
    {'run': cmd_run, 'report': cmd_report}[a.cmd](a)


if __name__ == '__main__':
    main()
