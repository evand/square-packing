#!/usr/bin/env python3
"""Record batch: register witness -> exactsolve; if not certified: slp2 polish -> exactsolve.

  run.py [n ...]        (default: every n in candidates.json)  -> work/n-N/, results.json, table on stdout
Each stage is a subprocess with a timeout; trials run through jobpool (one process per n).
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import json, subprocess, sys, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
EXD = os.path.dirname(HERE)                                   # search/exact
PK = os.path.join(os.path.dirname(EXD), 'packer')               # search/packer (slp2, jobpool; local)
sys.path.insert(0, PK)
from jobpool import run_jobs

PY = sys.executable


def sh(args, timeout, log):
    t = time.time()
    try:
        r = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        out = r.stdout + r.stderr
    except subprocess.TimeoutExpired:
        out = f'TIMEOUT after {timeout}s'
    with open(log, 'a') as f:
        f.write(f'$ {" ".join(args)}  ({time.time() - t:.0f}s)\n{out}\n')
    return out


def exact(inp, outdir, log, timeout):
    sh([PY, os.path.join(EXD, 'exactsolve.py'), inp, '--out', outdir, '-q'], timeout, log)
    p = os.path.join(outdir, os.path.basename(inp)[:-4] + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def job(args):
    n, reg = args
    d = os.path.join(HERE, 'work', f'n-{n}')
    os.makedirs(d, exist_ok=True)
    log = os.path.join(d, 'log.txt')
    w = os.path.join(d, 'witness.txt')
    sh([PY, os.path.join(HERE, 'register.py'), 'witness', str(n), w], 60, log)
    tmo = 1800 if n < 200 else 3600
    stages = []
    r = exact(w, d, log, tmo)
    stages.append(('witness', r and r.get('status'), r and r.get('S_exact')))
    if not (r and r.get('cert_valid')):
        p = os.path.join(d, 'polished.txt')
        sh([PY, os.path.join(PK, 'slp2.py'), w, p, '--R', '1e-4'], tmo, log)
        if os.path.exists(p):
            r = exact(p, d, log, tmo)
            stages.append(('slp2', r and r.get('status'), r and r.get('S_exact')))
    res = dict(n=n, register=reg, stages=stages)
    if r and r.get('cert_valid'):
        S = r['S_exact']
        res.update(S_exact=S, S_cert=r['S_cert_decimal'], lam=r.get('lambda_A_maxmin'),
                   second=(r.get('second_order') or {}).get('status'), vv=r.get('vv_milp_dS'), free=r.get('free_squares'),
                   eq_res=r.get('equilibrium_residual_exact'), forced=r.get('forced_pairs'),
                   kkt=bool(r.get('lambda_A_maxmin') is not None and (r.get('equilibrium_residual_exact') or 0) < 1e-30),
                   delta_vs_register=float(Fraction(r['S_cert']) - Fraction(reg)))
    return res


if __name__ == '__main__':
    args = sys.argv[1:]
    cf = 'candidates_all.json'
    if args and args[0].startswith('--cands='):
        cf = args.pop(0).split('=', 1)[1]
    C = json.load(open(os.path.join(HERE, cf)))
    want = {int(a) for a in args}
    done = set()
    if True:                                         # skip n already certified
        for c in C:
            d = os.path.join(HERE, 'work', f"n-{c['n']}")
            if any(os.path.exists(os.path.join(d, x + '.json')) and json.load(open(os.path.join(d, x + '.json'))).get('cert_valid')
                   for x in ('witness', 'polished')):
                done.add(c['n'])
    jobs = [(c['n'], c['register']) for c in C if (not want or c['n'] in want) and c['n'] not in done]
    jobs.sort(key=lambda j: -j[0])                    # big n first (longest)
    R = run_jobs(job, jobs, procs=14, timeout=3 * 3600, on_timeout=lambda j: dict(n=j[0], register=j[1], stages='KILLED'))
    R = sorted((r for r in R if r), key=lambda r: r['n'])
    old = {}
    p = os.path.join(HERE, 'results.json')
    if os.path.exists(p):
        old = {r['n']: r for r in json.load(open(p))}
    old.update({r['n']: r for r in R})
    json.dump(sorted(old.values(), key=lambda r: r['n']), open(p, 'w'), indent=1)
    for r in R:
        if 'S_exact' in r:
            print(f"{r['n']:4d} reg {r['register']:<22} cert {r['S_cert'][:24]}  Δ {r['delta_vs_register']:+.2e}  "
                  f"λmin {(r['lam'] if r['lam'] is not None else float('nan')):.1e}  {'KKT' if r['kkt'] else 'bound only'}  {(r['second'] or '')[:30]}  via {r['stages'][-1][0]}", flush=True)
        else:
            print(f"{r['n']:4d} reg {r['register']:<22} NOT CERTIFIED  {r['stages']}", flush=True)
