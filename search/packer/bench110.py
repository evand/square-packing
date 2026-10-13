#!/usr/bin/env python3
"""Basin-entry benchmark at n = 110 (10-13; PACKER "E0", "R1"): the n = 110 specifics for pk/benchkit races.

The deep 110 minima are two basins (A: the record + 10, B: 39 from 10.9968870; 0.049 apart).  Holdout
refsets/h110_basinA withholds A + every minimum within fingerprint distance 0.03 of it (78 of 939); the best visible is
B's 10.9968869868.  Event = the first explorer proposal with a sub-11 side below that (entering a withheld basin: the
analog of finding a new record by finding a new basin).

Search = explore.py with census v4's recipe, 1 proc per replicate.  Start pools (a fresh random sample per replicate):
  rand    40 random visible minima of the holdout (representative packings)
  above   40 random n = 110 store packings with 11 < s < 11.05 (grid crossing)
Arms = name:pool[:flags], flags = explore.py flags separated by ';', each value after its first '=':
  control:rand   nomirror:rand:--weights=mirror=0   finish3:rand:--finish=3   both:rand:--weights=mirror=0;--finish=3

  python3 bench110.py race --arms control:rand nomirror:rand:--weights=mirror=0 --out runs/bench110/r1 \
                           [--budget 16 --cap-min 30 --check 2 --margin 1.5 --slots 16]
  python3 bench110.py report runs/bench110/r1        (benchkit summary + entering proposal kinds)
  python3 bench110.py e0report runs/bench110/e0      (the fixed-hour E0 layout)
"""
from __future__ import annotations
import argparse, collections, glob, json, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HOLD = 'h110_basinA'
BEST_VISIBLE = 10.9968869868
RECIPE = ['--n', '110', '--adapt', '--reward', 'legacy', '--sparse-share', '0.5']
NSTART = 40


# ---------------------------------------------------------------- start pools

def _write(path, p):
    with open(path, 'w') as f:
        f.write(f'110 {float(p.s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in p.sq))


class Pools:
    """Candidate start packings per pool (store ids), sampled per replicate seed."""
    def __init__(self):
        from pk import refset
        from pk.store import Store
        self.st = Store()
        H = json.load(open(f'{refset.ROOT}/{HOLD}/holdout.json'))
        R = refset.load(H['ref'])
        vis = set(H['visible'])
        self.ids = {'rand': [m['rep'] for m in R['minima'] if m['key'] in vis and m['rep']]}
        self.ids['above'] = [r[0] for r in self.st.db.execute(
            "SELECT id FROM packing WHERE n=110 AND s > ? AND s < 11.05 AND status != 'infeasible'", (11 + 1e-7,))]

    def starts(self, pool):
        def fn(seed, d):
            os.makedirs(os.path.join(d, 'starts'), exist_ok=True)
            out = []
            for i in random.Random(seed).sample(self.ids[pool], NSTART):
                f = os.path.join(d, 'starts', f'p{i}.txt'); _write(f, self.st.get(i)); out.append(f)
            return out
        return fn


# ---------------------------------------------------------------- event + launch

_OFF = {}


def entry_event(d):
    """First proposal in d/proposals.jsonl with a sub-11 side below the best visible (read incrementally)."""
    p = os.path.join(d, 'proposals.jsonl')
    if not os.path.exists(p):
        return None
    with open(p) as f:
        f.seek(_OFF.get(d, 0))
        for l in f:
            if not l.endswith('\n'):
                break
            _OFF[d] = _OFF.get(d, 0) + len(l)
            try:
                q = json.loads(l)
            except ValueError:
                continue
            s = q.get('s')
            if s is not None and s < BEST_VISIBLE - 1e-9 and q.get('st') not in ('screen-grid', 'discard', 'screen-discard'):
                return dict(s=s, kind=q.get('kind'), desc=q.get('desc'), t=q.get('t'))
    return None


def launcher(bank, novel_ref):
    def fn(arm, d, seed, starts):
        return [sys.executable, os.path.join(HERE, 'explore.py'), *RECIPE, *arm.args, '--fp-bank', bank,
                '--starts', *starts, '--minutes', '600', '--procs', '1', '--seed', str(seed), '--novel-ref', novel_ref,
                '--out', d]
    return fn


def parse_arm(spec, pools):
    """name:pool[:args]  args = explore flags joined by ';' with '=' for values: --weights=mirror=0;--finish=3"""
    from pk.benchkit import Arm
    parts = spec.split(':', 2)
    name, pool = parts[0], parts[1]
    args = []
    for tok in (parts[2].split(';') if len(parts) > 2 and parts[2] else []):
        flag, _, val = tok.partition('=')
        args += [flag] + ([val] if val else [])
    return Arm(name, args, pools.starts(pool))


# ---------------------------------------------------------------- commands

def race_cmd(a):
    from pk import benchkit, refset, manifest, fingerprint
    import scratch
    scratch.scope('bench110')
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'bench110-race', args=vars(a), procs=a.slots)
    bank = os.path.join(a.out, 'fp_bank.npz')
    if not os.path.exists(bank):
        fingerprint.bank(HOLD, bank)
    pools = Pools()
    arms = [parse_arm(s, pools) for s in a.arms]
    json.dump([dict(name=x.name, args=x.args) for x in arms], open(os.path.join(a.out, 'arms.json'), 'w'), indent=1)
    S = benchkit.race(arms, launcher(bank, f'{refset.ROOT}/{HOLD}/novel_ref.json'), entry_event, a.out, slots=a.slots,
                      cap_cpu_s=a.cap_min * 60, budget_cpu_h=a.budget, check_cpu_h=a.check, margin=a.margin,
                      control=arms[0].name)
    kinds(a.out)
    man.finish(summary=S)


def kinds(out):
    from pk import benchkit
    R = benchkit.read_results(out)
    c = collections.defaultdict(collections.Counter)
    for r in R:
        if r['event']:
            c[r['arm']][r['ev'].get('kind')] += 1
    for arm, k in c.items():
        print(f'entering kinds, {arm}: ' + ', '.join(f'{x} {v}' for x, v in k.most_common()))


def report_cmd(a):
    from pk import benchkit
    arms = json.load(open(os.path.join(a.dir, 'arms.json')))
    benchkit.summary(benchkit.read_results(a.dir), arms[0]['name'], margin=a.margin)
    kinds(a.dir)


def e0report(a):
    """E0 (fixed 60-min replicates, 10-13): CPU to entry per replicate from proposals.jsonl."""
    for d in sorted(glob.glob(os.path.join(a.dir, '*_r*'))):
        cpu, ev = 0.0, None
        for l in open(os.path.join(d, 'proposals.jsonl')):
            q = json.loads(l); cpu += q.get('cpu') or 0.0
            s = q.get('s')
            if ev is None and s is not None and s < BEST_VISIBLE - 1e-9 and q.get('st') not in ('screen-grid', 'discard', 'screen-discard'):
                ev = (round(cpu), q.get('kind'), s)
        print(os.path.basename(d), ev)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('race'); r.add_argument('--arms', nargs='+', required=True); r.add_argument('--out', required=True)
    r.add_argument('--budget', type=float, default=16, help='CPU-h, all arms'); r.add_argument('--cap-min', type=float, default=30)
    r.add_argument('--check', type=float, default=2, help='CPU-h per arm between sequential checks')
    r.add_argument('--margin', type=float, default=1.5); r.add_argument('--slots', type=int, default=16)
    q = sp.add_parser('report'); q.add_argument('dir'); q.add_argument('--margin', type=float, default=1.5)
    e = sp.add_parser('e0report'); e.add_argument('dir')
    a = ap.parse_args()
    {'race': race_cmd, 'report': report_cmd, 'e0report': e0report}[a.cmd](a)
