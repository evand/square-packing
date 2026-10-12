#!/usr/bin/env python3
"""Census v3 at n = 110: near-record holdout, far vs random starts x distance-rule parents (pre-registration: PACKER.md
"Census v3"; Evan 10-11).

Holdout: refsets/h110_near_b (snapshot census110_1011b, 780 minima; the 62 within 1e-3 of the record withheld, record
included).  Explorer census recipe as census v2 (--adapt --reward legacy, novelty vs the visible minima, 1 proc per run).
2 x 2 arms:
  far          starts = the 40 visible minima farthest (matching distance) from the best visible minimum
  far_sparse   far + --sparse-share 0.5 (Archive.pick_sparse: parents weighted by fingerprint distance to the nearest
               known configuration, bank = the visible minima; Evan's distance rule)
  rand         starts = 40 random visible minima (the holdout's starts/)
  rand_sparse  rand + --sparse-share 0.5
Each arm x --reps, all concurrent, --minutes each.

  census_v3.py run [--minutes 180 --reps 2]
  census_v3.py certify [--procs 8]
  census_v3.py report
"""
import argparse, glob, json, math, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
HOLD = 'h110_near_b'
OUT = f'{HERE}/runs/census_v3'
BANK = f'{OUT}/fp_bank.npz'
RECIPE = ['--n', '110', '--adapt', '--reward', 'legacy']
SPARSE = ['--sparse-share', '0.5', '--fp-bank', BANK]
ARMS = {'far': ('far', []), 'far_sparse': ('far', SPARSE), 'rand': ('rand', []), 'rand_sparse': ('rand', SPARSE)}


def far_starts(k=40):
    from pk import refset
    from pk.store import Store
    from pk.crossing import match_distance
    d = f'{OUT}/far_starts'
    if os.path.exists(f'{d}/far.json'):
        return sorted(glob.glob(f'{d}/*.txt'))
    H = json.load(open(f'{refset.ROOT}/{HOLD}/holdout.json'))
    R = refset.load(H['ref'])
    vis = [m for m in R['minima'] if m['key'] in set(H['visible']) and m['rep']]
    st = Store()
    best = min(vis, key=lambda m: float(m['S'][:20]))
    pb = st.get(best['rep'])
    from jobpool import run_jobs
    D = sorted(zip(run_jobs(_md, [(m['rep'], best['rep']) for m in vis], procs=8, timeout=60), vis),
               key=lambda t: -(t[0] or 0))[:k]
    os.makedirs(d, exist_ok=True)
    for dist, m in D:
        st.get(m['rep']).write(f'{d}/m{float(m["S"][:20]):.10f}_d{dist:.3f}.txt')
    json.dump(dict(best=best['key'], best_S=best['S'], picked=[(m['key'], dd) for dd, m in D]), open(f'{d}/far.json', 'w'), indent=0)
    print(f'far starts: {len(D)}, matching distance {D[-1][0]:.3f} .. {D[0][0]:.3f} from the best visible ({float(best["S"][:20]):.10f})')
    return sorted(glob.glob(f'{d}/*.txt'))


def _md(arg):
    from pk.store import Store
    from pk.crossing import match_distance
    st = Store()
    return match_distance(st.get(arg[0]), st.get(arg[1]))


def run(a):
    from pk import refset, manifest, fingerprint
    import scratch; scratch.scope('census_v3')
    os.makedirs(OUT, exist_ok=True)
    if not os.path.exists(BANK):
        fingerprint.bank(HOLD, BANK)
    man = manifest.start(OUT, 'census_v3', args=vars(a), procs=len(ARMS) * a.reps, extra=dict(holdout=HOLD))
    hd = f'{refset.ROOT}/{HOLD}'
    S = dict(rand=sorted(glob.glob(f'{hd}/starts/*.txt')), far=far_starts())
    P = []
    for arm, (sk, extra) in ARMS.items():
        for r in range(a.reps):
            out = f'{OUT}/{arm}_r{r}'
            cmd = [sys.executable, f'{HERE}/explore.py', *RECIPE, '--starts', *S[sk], '--minutes', str(a.minutes), '--procs', '1',
                   '--seed', str(2000 + 10 * r + list(ARMS).index(arm)), '--novel-ref', f'{hd}/novel_ref.json', '--out', out, *extra]
            os.makedirs(out, exist_ok=True)
            P.append((arm, r, subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE)))
            print(f'{arm} r{r}: {len(S[sk])} {sk} starts {" ".join(extra)}', flush=True)
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
    rec = None
    rows = {}
    for d in sorted(glob.glob(f'{OUT}/*_r*')):
        if not os.path.exists(f'{d}/archive.jsonl'):
            continue
        arm = os.path.basename(d).rsplit('_r', 1)[0]
        o = refset.score(HOLD, [d])
        H = json.load(open(f'{refset.ROOT}/{HOLD}/holdout.json')); R = refset.load(H['ref'])
        rec = min(float(m['S'][:20]) for m in R['minima'])
        E = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
        t = [e['t'] for e in E if abs(e['s'] - rec) < 1e-9 and e.get('parent', 0) >= 0]
        o['record_t'] = min(t) if t else None
        o['best'] = min(e['s'] for e in E)
        nn = sorted(e['nnd'] for e in E if e['s'] < 11 and e.get('nnd') is not None and math.isfinite(e['nnd']))
        o['nnd_median'] = nn[len(nn) // 2] if nn else None
        rows.setdefault(arm, []).append(o)
    json.dump(rows, open(f'{OUT}/report.json', 'w'), indent=1)
    from scipy.stats import beta
    def ratio_ci(x1, t1, x0, t0):
        n = x1 + x0
        if n == 0: return float('nan'), float('nan')
        lo = beta.ppf(0.025, x1, x0 + 1) if x1 > 0 else 0.0
        hi = beta.ppf(0.975, x1 + 1, x0) if x0 > 0 else float('inf')
        f = lambda p: p / (1 - p) * t0 / t1 if p < 1 else float('inf')
        return f(lo), f(hi)
    print(f'\n{"arm":12} {"CPU-h":>6} {"withheld":>8} {"/CPU-h":>7} {"vs base (95% CI)":>22} {"record t (s)":>16} {"best":>14} {"new":>4}  per-rep withheld')
    for arm, R in rows.items():
        base = rows.get(arm.replace('_sparse', '')) if arm.endswith('_sparse') else None
        w, h = sum(o['withheld_rediscovered'] for o in R), sum(o['cpu_h'] for o in R)
        if base:
            bw, bh = sum(o['withheld_rediscovered'] for o in base), sum(o['cpu_h'] for o in base)
            lo, hi = ratio_ci(w, h, bw, bh); rr = (w / h) / (bw / bh) if bw else float('nan')
            vs = f'{rr:5.2f} ({lo:4.2f}, {hi:4.2f})'
        else:
            vs = '—'
        print(f'{arm:12} {h:6.2f} {w:8d} {w / h if h else 0:7.1f} {vs:>22} {str([o["record_t"] for o in R]):>16} '
              f'{min(o["best"] for o in R):14.10f} {sum(o["new"] for o in R):4d}  {[o["withheld_rediscovered"] for o in R]}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('run'); p.add_argument('--minutes', type=float, default=180); p.add_argument('--reps', type=int, default=2)
    p = sp.add_parser('certify'); p.add_argument('--procs', type=int, default=8)
    sp.add_parser('report')
    a = ap.parse_args()
    {'run': run, 'certify': certify, 'report': report}[a.cmd](a)
