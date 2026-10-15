#!/usr/bin/env python3
"""n = 110 holdout test of the parent value model (pre-registration: PACKER.md "H110-PM", 10-14).

Holdout refsets/h110_rand50 (census110_1011; 390 visible / 370 withheld minima, record withheld; 40 visible starts).
Arms (explore.py directly, current explorer config as in the TT3 chain rounds: --adapt --reward value --elite-share 0.3,
polish gate 2e-3, novelty vs the visible minima, frontier_s = best visible side):
  G    gate only (old parent rule)
  PF   G + --pm-share 1 (parent value model)
--reps replicates per arm x --procs procs, all concurrent, --minutes each.

  python3 analysis/h110_pm.py run [--reps 4 --procs 2 --minutes 10]
  python3 analysis/h110_pm.py score [--json runs/h110_pm/score.json]
Score: certified withheld minima re-found (side key -> store side_minima, as pk/refset.score), per measured CPU-h, by
stratum: side gap to the record (<1e-3 / 1e-3..3e-3 / >=3e-3) and matching distance to the record (< / >= 0.05).
"""
import argparse, glob, json, math, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
HOLD = 'h110_rand50'
OUT = os.path.join(PK, 'runs', 'h110_pm')
ARMS = {'G': ['--polish-gate', '2e-3'], 'PF': ['--polish-gate', '2e-3', '--pm-share', '1']}
MD_FAR = 0.05


def holdout():
    from pk import refset
    H = json.load(open(f'{refset.ROOT}/{HOLD}/holdout.json')); R = refset.load(H['ref'])
    assert R['sha256'] == H['ref_sha256']
    return H, R


def run(a):
    from pk import refset, manifest
    H, R = holdout()
    S = {m['key']: float(m['S'][:20]) for m in R['minima']}
    bv = min(S[k] for k in H['visible'])
    hd = f'{refset.ROOT}/{HOLD}'
    starts = sorted(glob.glob(f'{hd}/starts/*.txt'))
    os.makedirs(OUT, exist_ok=True)
    man = manifest.start(OUT, 'h110_pm', args=vars(a), procs=len(ARMS) * a.reps * a.procs, extra=dict(holdout=HOLD, best_visible=bv))
    P = []
    for r in range(a.reps):
        for i, (arm, extra) in enumerate(ARMS.items()):
            out = f'{OUT}/{arm}_r{r}'; os.makedirs(out, exist_ok=True)
            cmd = [sys.executable, f'{PK}/explore.py', '--n', '110', '--adapt', '--reward', 'value', '--elite-share', '0.3',
                   '--frontier-s', repr(bv), '--novel-ref', f'{hd}/novel_ref.json', '--starts', *starts,
                   '--minutes', str(a.minutes), '--procs', str(a.procs), '--seed', str(2000 + 10 * r + i), '--out', out, *extra]
            P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=PK))
    for p in P:
        p.wait()
    man.finish()


def md_to_record(R):
    """matching distance of every withheld / visible minimum's representative to the record's (cached)"""
    cache = f'{OUT}/md_record.json'
    if os.path.exists(cache):
        return json.load(open(cache))
    from pk.store import Store
    from pk import cluster as CL
    st = Store()
    rec = min(R['minima'], key=lambda m: float(m['S'][:20]))
    pr = st.get(rec['rep'])
    D = {}
    for m in R['minima']:
        if not m['rep']:
            continue
        p = st.get(m['rep'])
        D[m['key']] = float(CL.md(pr.sq, pr.s, p.sq, p.s))
    os.makedirs(OUT, exist_ok=True)
    json.dump(D, open(cache, 'w'))
    return D


def found_in(d, SM, k, tmax=None):
    f = set()
    for l in open(f'{d}/archive.jsonl'):
        e = json.loads(l)
        if e['s'] >= k or e.get('parent', 0) < 0 or e.get('unpolished') or (tmax is not None and e.get('t', 0) > tmax):
            continue
        from pk import identity as ID
        sk = ID.side_key(e['s'])
        if sk in SM and SM[sk][1] == ID.CERTIFIED:
            f.add(SM[sk][0])
    return f


def score(a):
    from pk.store import Store
    H, R = holdout()
    S = {m['key']: float(m['S'][:20]) for m in R['minima']}
    rec = min(S.values())
    D = md_to_record(R)
    wh = set(H['withheld'])
    st = Store(); SM = st.side_minima(110)
    strata = {'all': lambda k: True,
              'gap<1e-3': lambda k: S[k] - rec < 1e-3, 'gap 1e-3..3e-3': lambda k: 1e-3 <= S[k] - rec < 3e-3,
              'gap>=3e-3': lambda k: S[k] - rec >= 3e-3,
              f'md<{MD_FAR}': lambda k: D.get(k, 9) < MD_FAR, f'md>={MD_FAR}': lambda k: D.get(k, 9) >= MD_FAR}
    tot = {name: sum(1 for k in wh if f(k)) for name, f in strata.items()}
    res = {}
    for arm in ARMS:
        reps = []
        for d in sorted(glob.glob(f'{OUT}/{arm}_r*')):
            cpu = json.load(open(f'{d}/manifest.json'))['cpu_s'] / 3600
            f = found_in(d, SM, H['k']) & wh
            reps.append(dict(run=d, cpu_h=cpu, found=sorted(f), rec=any(S[x] == rec for x in f),
                             **{name: sum(1 for x in f if g(x)) for name, g in strata.items()}))
        res[arm] = reps
    print(f'withheld per stratum: {tot}')
    out = dict(totals=tot, arms=res, ratios={})
    for name in strata:
        line = []
        for arm in ARMS:
            rates = [r[name] / r['cpu_h'] for r in res[arm]]
            pooled = sum(r[name] for r in res[arm]) / sum(r['cpu_h'] for r in res[arm])
            line.append(f'{arm} {pooled:6.1f}/CPU-h (reps ' + ' '.join(f'{x:.0f}' for x in rates) + ')')
        g = [math.log(max(r[name], 0.5) / r['cpu_h']) for r in res['G']]; p = [math.log(max(r[name], 0.5) / r['cpu_h']) for r in res['PF']]
        mg, mp = sum(g) / len(g), sum(p) / len(p)
        vg = sum((x - mg) ** 2 for x in g) / (len(g) - 1); vp = sum((x - mp) ** 2 for x in p) / (len(p) - 1)
        se = math.sqrt(vg / len(g) + vp / len(p)) or 1e-9
        dfw = se ** 4 / ((vg / len(g)) ** 2 / (len(g) - 1) + (vp / len(p)) ** 2 / (len(p) - 1) + 1e-300)
        from scipy.stats import t as T
        q = T.ppf(0.975, max(dfw, 1))
        rr = math.exp(mp - mg)
        out['ratios'][name] = dict(ratio=rr, lo=math.exp(mp - mg - q * se), hi=math.exp(mp - mg + q * se))
        print(f'{name:16s} ({tot[name]:3d} withheld)  ' + '   '.join(line) +
              f'   PF/G {rr:.2f} [{out["ratios"][name]["lo"]:.2f}, {out["ratios"][name]["hi"]:.2f}]')
    for tm in (150, 300):                          # saturation check: withheld found by wall time tm (s)
        c = {arm: sum(len(found_in(r['run'], SM, H['k'], tm) & wh) for r in res[arm]) for arm in ARMS}
        out[f'by_{tm}s'] = c
        print(f'withheld found by {tm} s (summed over reps): ' + '  '.join(f'{k} {v}' for k, v in c.items()) +
              f'  PF/G {c["PF"] / max(c["G"], 1):.2f}')
    for arm in ARMS:
        print(f'{arm}: record re-found in {sum(r["rec"] for r in res[arm])} / {len(res[arm])} reps; union of withheld found '
              f'{len(set().union(*[set(r["found"]) for r in res[arm]]))}')
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['run', 'score'])
    ap.add_argument('--reps', type=int, default=4); ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--minutes', type=float, default=10); ap.add_argument('--json')
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
