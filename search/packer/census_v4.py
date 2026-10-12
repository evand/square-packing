#!/usr/bin/env python3
"""Census v4 at n = 110: region-count holdouts ("can you find a family you have never seen"; pre-registration: PACKER.md
"Census v4"; Evan 10-11).

Holdouts (refset.holdout --rule groups on census110_1011b; groups = explore.roles_j counts, L / B mirror pairs merged):
  J0     every group with no ~45 deg junction squares (186 minima, record included): recovery needs a structural change
  grand  random groups to >= 24 % of minima (217 minima, record visible): the matched unbiased comparison
  recnb  the record's group + every group within 1 in each count (156 minima, record included): second batch
Recipe: census v3's plus the distance rule (--sparse-share 0.5, bank = the holdout's visible minima); starts = 40 random
visible minima of each holdout; 1 proc per run, all runs of a batch concurrent.

  census_v4.py run --batch a [--minutes 240 --reps 4]     batch a = J0 + grand;  batch b = recnb
  census_v4.py certify [--procs 8]
  census_v4.py report
"""
import argparse, glob, json, math, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = f'{HERE}/runs/census_v4'
HOLDS = {'J0': 'h110_J0', 'grand': 'h110_grand', 'recnb': 'h110_recnb'}
BATCH = {'a': ['J0', 'grand'], 'b': ['recnb']}
RECIPE = ['--n', '110', '--adapt', '--reward', 'legacy', '--sparse-share', '0.5']


def run(a):
    from pk import refset, manifest, fingerprint
    import scratch; scratch.scope(f'census_v4{a.batch}')
    os.makedirs(OUT, exist_ok=True)
    man = manifest.start(f'{OUT}/batch_{a.batch}', f'census_v4{a.batch}', args=vars(a), procs=len(BATCH[a.batch]) * a.reps)
    P = []
    for arm in BATCH[a.batch]:
        hd = f'{refset.ROOT}/{HOLDS[arm]}'
        bank = f'{OUT}/fp_bank_{arm}.npz'
        if not os.path.exists(bank):
            fingerprint.bank(HOLDS[arm], bank)
        starts = sorted(glob.glob(f'{hd}/starts/*.txt'))
        for r in range(a.reps):
            out = f'{OUT}/{arm}_r{r}'
            cmd = [sys.executable, f'{HERE}/explore.py', *RECIPE, '--fp-bank', bank, '--starts', *starts, '--minutes', str(a.minutes),
                   '--procs', '1', '--seed', str(4000 + 10 * r + list(HOLDS).index(arm)), '--novel-ref', f'{hd}/novel_ref.json', '--out', out]
            os.makedirs(out, exist_ok=True)
            P.append((arm, r, subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE)))
            print(f'{arm} r{r}: {len(starts)} starts, holdout {HOLDS[arm]}', flush=True)
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
        mf = f'{d}/manifest.json'
        if not os.path.exists(f'{d}/archive.jsonl') or not os.path.exists(mf) or json.load(open(mf)).get('status') != 'ok':
            continue                         # unfinished runs (no measured CPU yet)
        arm = os.path.basename(d).rsplit('_r', 1)[0]
        H = json.load(open(f'{refset.ROOT}/{HOLDS[arm]}/holdout.json')); R = refset.load(H['ref'])
        o = refset.score(HOLDS[arm], [d])
        rec = min(float(m['S'][:20]) for m in R['minima'])
        E = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
        t = [e['t'] for e in E if abs(e['s'] - rec) < 1e-9 and e.get('parent', 0) >= 0]
        o['record_t'] = min(t) if t else None
        o['best'] = min(e['s'] for e in E)
        o['frac'] = o['withheld_rediscovered'] / o['withheld_total']
        # withheld groups entered: first time a withheld-group minimum was re-found (any of them)
        o['first_withheld_t'] = min(o['withheld_first_t']) if o['withheld_first_t'] else None
        rows.setdefault(arm, []).append(o)
    json.dump(rows, open(f'{OUT}/report.json', 'w'), indent=1)
    print(f'\n{"arm":6} {"CPU-h":>6} {"withheld":>9} {"frac":>6} {"frac/CPU-h":>10} {"first (s)":>22} {"record t (s)":>28} {"new":>4}  per-rep re-found')
    for arm, R in rows.items():
        w, h, tot = sum(o['withheld_rediscovered'] for o in R), sum(o['cpu_h'] for o in R), R[0]['withheld_total']
        print(f'{arm:6} {h:6.2f} {w:5d}/{tot:<3d} {w / len(R) / tot:6.2f} {w / tot / h:10.3f} {str([o["first_withheld_t"] for o in R]):>22} '
              f'{str([o["record_t"] for o in R]):>28} {sum(o["new"] for o in R):4d}  {[o["withheld_rediscovered"] for o in R]}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('run'); p.add_argument('--batch', choices=list(BATCH), required=True)
    p.add_argument('--minutes', type=float, default=240); p.add_argument('--reps', type=int, default=4)
    p = sp.add_parser('certify'); p.add_argument('--procs', type=int, default=8)
    sp.add_parser('report')
    a = ap.parse_args()
    {'run': run, 'certify': certify, 'report': report}[a.cmd](a)
