#!/usr/bin/env python3
"""FP1 bench (10-15): fillfast vs gen.mis_fill2 on real fill instances.

Corpus: (a) the final EV2 gen-1 programs re-run in the sandbox on n = 234-240 (seed 1000 + n), every candidate's tilted
squares as one instance; (b) the tilted squares of store-best packings at a spread of n (records: the fill should
give n).  Per instance: reference count + CPU (mis_fill2), fillfast exact / greedy / tolerant (tau) counts + CPU +
kernel stats, and an independence check of every returned set.  Downstream (--downstream): for instances whose
reference count reaches n, trim + fq quench + polish for the reference and each mode, compare sides.

  python3 analysis/fp1_bench.py --out runs/fp1 [--max 200 --procs 12 --downstream]
"""
import argparse, glob, json, math, os, subprocess, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
PROGS = {'P007': 'prog_10.py', 'P008': 'prog_8.py', 'P009': 'prog_8.py', 'P010': 'prog_8.py'}
TAUS = (0.005, 0.01, 0.02)


def corpus_ev2(maxc):
    out = []
    for tag, pf in PROGS.items():
        prog = f'{PK}/runs/evolve2/gen1/{tag}/{pf}'
        for n in range(234, 241):
            with tempfile.TemporaryDirectory() as td:
                r = subprocess.run(['unshare', '-rn', sys.executable, '-I', f'{HERE}/evolve_runner.py', prog, str(n),
                                    str(1000 + n), str(maxc), '1'], cwd=td, capture_output=True, text=True, timeout=120)
            try:
                res = json.loads(r.stdout.strip().splitlines()[-1])
            except Exception:
                continue
            for s, sq in res.get('cands', []):
                a = np.asarray(sq, float); t = np.minimum(a[:, 2] % 90, 90 - a[:, 2] % 90) > 1
                out.append(dict(src=f'{tag}:n{n}', n=n, s=float(s), T=a[t].tolist()))
    return out


def corpus_records(ns):
    import sqlite3
    from pk.store import Store
    from pk import eventbench as EB
    db = sqlite3.connect(f"file:{PK}/runs/store.sqlite?mode=ro", uri=True); st = Store(); out = []
    for n in ns:
        T0 = EB.frontier_at(n, None)
        r = db.execute("SELECT id FROM packing WHERE n=? AND abs(s-?)<2e-9 ORDER BY added LIMIT 1", (n, T0)).fetchone()
        if not r:
            continue
        p = st.get(r[0]); sq = np.asarray(p.sq, float); a = sq[:, 2] % 90
        T = sq[np.minimum(a, 90 - a) > 1]
        if len(T) and float(p.s) < math.ceil(math.sqrt(n) - 1e-12):
            out.append(dict(src=f'record:n{n}', n=n, s=float(p.s), T=T.tolist()))
    return out


def independent(U):
    U = np.asarray(U)
    if len(U) < 2:
        return True
    d = np.abs(U[:, None, :2] - U[None, :, :2])
    bad = (d[..., 0] < 1 - 1e-6) & (d[..., 1] < 1 - 1e-6)
    np.fill_diagonal(bad, False)
    return not bad.any()


def job(inst):
    import gen, fillfast as F
    s, T = inst['s'], np.asarray(inst['T'], float).reshape(-1, 3)
    out = dict(src=inst['src'], n=inst['n'], s=s, nT=len(T))
    t = time.process_time(); L = gen.mis_fill2(s, [tuple(q) for q in T], time_limit=30); out['ref'] = len(L) + len(T)
    out['t_ref'] = time.process_time() - t
    for mode, tau in [('exact', 0.0), ('greedy', 0.0)] + [('greedy', x) for x in TAUS]:
        key = mode if tau == 0 else f'tol{tau}'
        t = time.process_time(); R = F.fill(s, T, mode, tau=tau); dt = time.process_time() - t
        out[key] = len(R); out[f't_{key}'] = dt
        if tau == 0:
            out[f'indep_{key}'] = independent(R[len(T):])
        if key == 'exact':
            out.update(cand=F.LAST['cand'], kernel=F.LAST['kernel'], comps=F.LAST['comps'], maxcomp=F.LAST['maxcomp'])
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--max', type=int, default=200)
    ap.add_argument('--procs', type=int, default=12); ap.add_argument('--maxc', type=int, default=12)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from pk import manifest
    man = manifest.start(a.out, 'fp1_bench', args=vars(a), procs=a.procs)
    C = corpus_ev2(a.maxc)
    rng = np.random.default_rng(0)
    if len(C) > a.max:
        C = [C[i] for i in sorted(rng.choice(len(C), a.max, replace=False))]
    C += corpus_records([60, 80, 100, 110, 127, 131, 150, 156, 178, 182, 200, 207, 211, 236, 241, 259, 270, 292, 303])
    json.dump(C, open(f'{a.out}/corpus.json', 'w'))
    with ProcessPoolExecutor(a.procs) as ex:
        R = list(ex.map(job, C, chunksize=2))
    with open(f'{a.out}/results.jsonl', 'w') as f:
        for r in R:
            f.write(json.dumps(r) + '\n')
    man.finish()
    report(R)


def report(R):
    N = len(R)
    print(f'{N} instances ({sum(1 for r in R if r["src"].startswith("record"))} records)')
    for key in ['exact', 'greedy'] + [f'tol{x}' for x in TAUS]:
        d = np.array([r[key] - r['ref'] for r in R])
        sp = np.array([r['t_ref'] / max(r[f't_{key}'], 1e-6) for r in R])
        print(f'  {key:9s}: count - ref: =0 {np.mean(d == 0):.3f}  <0 {np.mean(d < 0):.3f} (min {d.min()})  >0 {np.mean(d > 0):.3f} '
              f'(max {d.max()}); CPU ref {sum(r["t_ref"] for r in R):.0f} s vs {sum(r[f"t_{key}"] for r in R):.0f} s; '
              f'median speedup {np.median(sp):.1f}x')
    print(f'  independence: exact {all(r["indep_exact"] for r in R)}  greedy {all(r["indep_greedy"] for r in R)}')
    k = np.array([r['kernel'] / max(r['cand'], 1) for r in R]); mc = np.array([r['maxcomp'] for r in R])
    print(f'  kernel share: median {np.median(k):.3f}, p90 {np.quantile(k, 0.9):.3f}, max {k.max():.3f}; '
          f'kernel = 0 on {np.mean(k == 0):.2f}; largest component: median {np.median(mc):.0f}, max {mc.max()}')


if __name__ == '__main__':
    main()
