#!/usr/bin/env python3
"""Basin-entry benchmark at n = 110 (10-13; Evan: "withhold one whole basin and the adjacent region"; PACKER "E0").

The deep 110 minima are two basins (A: the record + 10, B: 39 from 10.9968870; 0.049 apart).  Holdout
refsets/h110_basinA withholds A + every minimum within 0.03 of it (78 of 939); the best visible is B's 10.9968869868.
A run "enters" when it reaches any side below that: the analog of finding a new record by finding a new basin.

Arms (starts; everything else = census v4's recipe, 1 proc per replicate, all replicates concurrent):
  rand     40 random visible minima (refsets/h110_basinA/starts)
  above    40 random n = 110 store packings with 11 < s < 11.05 (grid crossing: the chain must first get below 11)
Report: per replicate, time (and measured CPU) to the first side below the best visible, the proposal kind that made
it, the best side reached; per proposal kind over all replicates: proposals, CPU, sub-11 rate, new-minimum rate.

  python3 bench110.py run --arms rand,above --reps 8 --minutes 60 --out runs/bench110/e0
  python3 bench110.py report runs/bench110/e0
"""
from __future__ import annotations
import argparse, collections, glob, json, os, random, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HOLD = 'h110_basinA'
BEST_VISIBLE = 10.9968869868
RECIPE = ['--n', '110', '--adapt', '--reward', 'legacy', '--sparse-share', '0.5']


def above_starts(d, k=40, seed=0, lo=11 + 1e-7, hi=11.05):
    from pk.store import Store
    st = Store()
    ids = [r[0] for r in st.db.execute("SELECT id FROM packing WHERE n=110 AND s > ? AND s < ? AND status != 'infeasible'",
                                       (lo, hi))]
    os.makedirs(d, exist_ok=True)
    for i in random.Random(seed).sample(ids, k):
        p = st.get(i)
        with open(f'{d}/p{i}.txt', 'w') as f:
            f.write(f'110 {float(p.s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in p.sq))
    return sorted(glob.glob(f'{d}/*.txt'))


def run(a):
    from pk import refset, manifest, fingerprint
    import scratch
    scratch.scope('bench110')
    os.makedirs(a.out, exist_ok=True)
    arms = a.arms.split(',')
    man = manifest.start(a.out, 'bench110', args=vars(a), procs=len(arms) * a.reps)
    hd = f'{refset.ROOT}/{HOLD}'
    bank = os.path.join(a.out, 'fp_bank.npz')
    if not os.path.exists(bank):
        fingerprint.bank(HOLD, bank)
    starts = {'rand': sorted(glob.glob(f'{hd}/starts/*.txt'))}
    if 'above' in arms:
        starts['above'] = above_starts(os.path.join(a.out, 'starts_above'))
    P = []
    for arm in arms:
        for r in range(a.reps):
            out = os.path.join(a.out, f'{arm}_r{r}')
            os.makedirs(out, exist_ok=True)
            cmd = [sys.executable, os.path.join(HERE, 'explore.py'), *RECIPE, *a.extra.split(), '--fp-bank', bank,
                   '--starts', *starts[arm], '--minutes', str(a.minutes), '--procs', '1',
                   '--seed', str(7000 + 10 * r + arms.index(arm)), '--novel-ref', f'{hd}/novel_ref.json', '--out', out]
            P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
    for p in P:
        p.wait()
    man.finish(summary=dict(runs=len(P)))
    report(argparse.Namespace(dir=a.out))


def report(a):
    rows, kinds = [], collections.defaultdict(lambda: collections.Counter())
    for d in sorted(glob.glob(os.path.join(a.dir, '*_r*'))):
        pp = os.path.join(d, 'proposals.jsonl')
        if not os.path.exists(pp):
            continue
        cpu = 0.0
        entry, best, sides = None, 99.0, set()
        for l in open(pp):
            p = json.loads(l)
            cpu += p.get('cpu') or 0.0
            k = p.get('kind', '?')
            kinds[k]['n'] += 1; kinds[k]['cpu'] += p.get('cpu') or 0.0
            s = p.get('s')
            if s is None or s >= 11 - 1e-7 or p.get('st') in ('screen-grid', 'discard', 'screen-discard'):
                continue
            kinds[k]['sub11'] += 1
            if p.get('st') == 'new?':
                kinds[k]['new'] += 1
            best = min(best, s)
            if entry is None and s < BEST_VISIBLE - 1e-9:
                entry = dict(t=p.get('t'), cpu=cpu, kind=k, s=s, desc=p.get('desc'))
                kinds[k]['entries'] += 1
        rows.append(dict(run=os.path.basename(d), cpu_h=cpu / 3600, entry=entry, best=best))
    print(f'{"run":10s} {"CPU-h":>6s}  entry (t s / CPU s / kind / side)                      best side')
    for r in rows:
        e = r['entry']
        es = f"{e['t']:6.0f} {e['cpu']:7.0f} {e['kind']:10s} {e['s']:.10f}" if e else '   -'
        print(f"{r['run']:10s} {r['cpu_h']:6.2f}  {es:52s} {r['best']:.10f}")
    by = collections.defaultdict(list)
    for r in rows:
        by[r['run'].rsplit('_r', 1)[0]].append(r)
    for arm, R in by.items():
        ent = [r['entry']['cpu'] for r in R if r['entry']]
        print(f'{arm}: entered {len(ent)} / {len(R)}; CPU s to entry {sorted(round(x) for x in ent)}; '
              f'entries per CPU-h {len(ent) / max(sum(r["cpu_h"] for r in R), 1e-9):.2f}')
    print('\nproposal kinds (all runs): proposals, CPU h, sub-11 per 1000, new per 1000, new per CPU-h, entries')
    for k, c in sorted(kinds.items(), key=lambda kv: -kv[1]['n']):
        print(f"  {k:12s} {c['n']:7d} {c['cpu'] / 3600:6.2f} {c['sub11'] / c['n'] * 1000:7.1f} {c['new'] / c['n'] * 1000:7.1f}"
              f" {c['new'] / max(c['cpu'] / 3600, 1e-9):8.1f} {c['entries']:3d}")
    json.dump(dict(rows=rows, kinds={k: dict(v) for k, v in kinds.items()}), open(os.path.join(a.dir, 'report.json'), 'w'),
              indent=1, default=str)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--arms', default='rand,above'); r.add_argument('--reps', type=int, default=8)
    r.add_argument('--minutes', type=float, default=60); r.add_argument('--out', required=True)
    r.add_argument('--extra', default='', help='extra explore.py arguments for every replicate (ablations)')
    q = sp.add_parser('report'); q.add_argument('dir')
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else report(a)
