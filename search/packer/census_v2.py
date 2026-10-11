#!/usr/bin/env python3
"""Census v2 at n = 110 (pre-registration: PACKER.md "Census v2"; claims CLAIMS.md; Evan 10-11).

Holdout: refsets/h110_rand50 (random half of census110_1011 withheld: 370 minima incl. the record).  Arms (explorer census
recipe = insert110's: --adapt --reward legacy, global novelty vs the holdout's visible minima, no elite share):
  control   40 random visible starts (the holdout's starts/)
  melt      control + --melt 2 (fixed melt, C10)
  finish    control + --finish 3 (inline kick + polish-only on archive-new sub-k results, C11)
  far       control recipe, starts = the 40 visible minima farthest (matching distance) from the best visible minimum (C05)
Each arm x --reps replicates, 1 proc each, all concurrent (same load), --minutes each.

  census_v2.py run [--minutes 60 --reps 2]     launch (wait for all)
  census_v2.py certify [--procs 8]               import the runs into the store, certify every new sub-11 side
  census_v2.py report                            per run / arm: refset.score vs the frozen holdout; rate ratios vs control
"""
import argparse, glob, json, math, os, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HOLD = 'h110_rand50'
OUT = f'{HERE}/runs/census_v2'
RECIPE = ['--n', '110', '--adapt', '--reward', 'legacy']
ARMS = {'control': [], 'melt': ['--melt', '2'], 'finish': ['--finish', '3'], 'far': []}


def far_starts(k=40):
    from pk import refset
    from pk.store import Store
    from pk.crossing import match_distance
    H = json.load(open(f'{refset.ROOT}/{HOLD}/holdout.json'))
    R = refset.load(H['ref'])
    vis = [m for m in R['minima'] if m['key'] in set(H['visible']) and m['rep']]
    st = Store()
    best = min(vis, key=lambda m: float(m['S'][:20]))
    pb = st.get(best['rep'])
    D = sorted(((match_distance(st.get(m['rep']), pb), m) for m in vis), key=lambda t: -t[0])[:k]
    d = f'{OUT}/far_starts'; os.makedirs(d, exist_ok=True)
    for dist, m in D:
        st.get(m['rep']).write(f'{d}/m{float(m["S"][:20]):.10f}_d{dist:.3f}.txt')
    json.dump(dict(best=best['key'], picked=[(m['key'], dd) for dd, m in D]), open(f'{d}/far.json', 'w'), indent=0)
    print(f'far starts: {len(D)}, matching distance {D[-1][0]:.3f} .. {D[0][0]:.3f} from the best visible ({float(best["S"][:20]):.10f})')
    return sorted(glob.glob(f'{d}/*.txt'))


def run(a):
    from pk import refset, manifest
    import scratch; scratch.scope('census_v2')
    os.makedirs(OUT, exist_ok=True)
    man = manifest.start(OUT, 'census_v2', args=vars(a), procs=len(ARMS) * a.reps, extra=dict(holdout=HOLD))
    hd = f'{refset.ROOT}/{HOLD}'
    starts = sorted(glob.glob(f'{hd}/starts/*.txt'))
    far = far_starts()
    P = []
    for arm, extra in ARMS.items():
        for r in range(a.reps):
            out = f'{OUT}/{arm}_r{r}'
            st = far if arm == 'far' else starts
            cmd = [sys.executable, f'{HERE}/explore.py', *RECIPE, '--starts', *st, '--minutes', str(a.minutes), '--procs', '1',
                   '--seed', str(1000 + 10 * r + list(ARMS).index(arm)), '--novel-ref', f'{hd}/novel_ref.json', '--out', out, *extra]
            os.makedirs(out, exist_ok=True)
            P.append((arm, r, subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE)))
            print(f'{arm} r{r}: {len(st)} starts {" ".join(extra)}', flush=True)
    for arm, r, p in P:
        p.wait()
        print(f'{arm} r{r} exit {p.returncode}', flush=True)
    man.finish()


def certify(a):
    runs = sorted(d for d in glob.glob(f'{OUT}/*_r*') if os.path.exists(f'{d}/archive.jsonl'))
    subprocess.run([sys.executable, f'{HERE}/pk.py', 'import-run', *runs], cwd=HERE, check=True)
    subprocess.run([sys.executable, f'{HERE}/census_certify.py', 'certify', '--n', '110', '--procs', str(a.procs)], cwd=HERE, check=True)


def report(a):
    from pk import refset
    rows = {}
    for d in sorted(glob.glob(f'{OUT}/*_r*')):
        if not os.path.exists(f'{d}/archive.jsonl'):
            continue
        arm = os.path.basename(d).rsplit('_r', 1)[0]
        o = refset.score(HOLD, [d])
        P = [json.loads(l) for l in open(f'{d}/proposals.jsonl')]
        o['proposals'] = len(P)
        rows.setdefault(arm, []).append(o)
    json.dump(rows, open(f'{OUT}/report.json', 'w'), indent=1)
    from scipy.stats import beta
    def ratio_ci(x1, t1, x0, t0):
        # conditional binomial CI for a Poisson rate ratio (x1 events in t1 vs x0 in t0)
        n = x1 + x0
        if n == 0: return float('nan'), float('nan')
        lo = beta.ppf(0.025, x1, x0 + 1) if x1 > 0 else 0.0
        hi = beta.ppf(0.975, x1 + 1, x0) if x0 > 0 else float('inf')
        f = lambda p: p / (1 - p) * t0 / t1 if p < 1 else float('inf')
        return f(lo), f(hi)
    c = rows.get('control', [])
    cw, cn, ch = sum(o['withheld_rediscovered'] for o in c), sum(o['new'] for o in c), sum(o['cpu_h'] for o in c)
    print(f'\n{"arm":8} {"reps":>4} {"CPU-h":>6} {"withheld":>9} {"/CPU-h":>7} {"ratio (95% CI) vs control":>27} {"new":>4} {"not-min sides":>13}  per-rep withheld')
    for arm, R in rows.items():
        w, nw, h = sum(o['withheld_rediscovered'] for o in R), sum(o['new'] for o in R), sum(o['cpu_h'] for o in R)
        lo, hi = ratio_ci(w, h, cw, ch) if arm != 'control' and ch else (float('nan'), float('nan'))
        rr = (w / h) / (cw / ch) if arm != 'control' and ch and cw else float('nan')
        print(f'{arm:8} {len(R):4d} {h:6.2f} {w:9d} {w / h if h else 0:7.1f} {rr:9.2f} ({lo:5.2f}, {hi:5.2f})      {nw:4d} '
              f'{sum(o["notmin_sides"] for o in R):13d}  {[o["withheld_rediscovered"] for o in R]}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('run'); p.add_argument('--minutes', type=float, default=60); p.add_argument('--reps', type=int, default=2)
    p = sp.add_parser('certify'); p.add_argument('--procs', type=int, default=8)
    sp.add_parser('report')
    a = ap.parse_args()
    {'run': run, 'certify': certify, 'report': report}[a.cmd](a)
