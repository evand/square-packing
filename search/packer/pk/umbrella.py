"""Free energy along Q4 = mean cos(4 theta) for hard squares at fixed pressure (umbrella sampling + WHAM).

  umbrella.py run NAME --n 110 --P 30 100 300 1000 [--q0 -0.9:1.0:0.1] [--kappa 4] [--sweeps 40000] [--seeds 2]
  umbrella.py report NAME

Each window: anneal sched from scratch (disks at low pressure -> squares by t = 0.2, pressure ramps to P by t = 0.3,
then held), harmonic bias kappa n (Q4 - q0)^2 / 2 throughout, samples (s, Q4) from t >= 0.5.  WHAM per pressure -> F(Q4)
in kT (min = 0) and the mean side per Q4 bin.  Two seeds per window: if their window means disagree by more than the
window width the run is not equilibrated (flagged).
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, collections, json, math, subprocess, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pk.anneal import ANNEAL

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'umbrella')


def window(job):
    n, P, q0, kappa, sweeps, seed = job
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [ANNEAL, 'sched', '--n', str(n), '--seed', str(seed), '--sweeps', str(sweeps), '--out', f'{tmp}/o.txt',
               '--bp', f'0:1,0.3:{P},1:{P}', '--r', '0:0.5,0.2:0,1:0', '--bias', 'harm', '--bias-lam', str(kappa),
               '--bias-q0', str(q0), '--log-every', '5', '--log-from', '0.5']
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        info = json.loads(r.stdout.strip().splitlines()[-1])
        L = np.loadtxt(f'{tmp}/o.txt.log')
    return dict(n=n, P=P, q0=q0, kappa=kappa, seed=seed, s=L[:, 1].tolist(), q=L[:, 2].tolist(), sec=info['sec'])


def wham(W, n, bins, iters=5000, tol=1e-10):
    """W: list of dict(q (samples), q0, kappa).  Returns (centers, F (kT, min 0), counts)."""
    edges = bins
    c = 0.5 * (edges[1:] + edges[:-1])
    H = np.array([np.histogram(w['q'], edges)[0] for w in W], float)        # (K, B)
    N = H.sum(1)
    U = np.array([0.5 * w['kappa'] * n * (c - w['q0']) ** 2 for w in W])     # bias (K, B)
    f = np.zeros(len(W))
    tot = H.sum(0)
    for _ in range(iters):
        den = (N[:, None] * np.exp(f[:, None] - U)).sum(0)
        p = np.where(tot > 0, tot / np.maximum(den, 1e-300), 0)
        fn = -np.log(np.maximum((p[None] * np.exp(-U)).sum(1), 1e-300))
        fn -= fn[0]
        if np.max(np.abs(fn - f)) < tol:
            f = fn; break
        f = fn
    with np.errstate(divide='ignore'):
        F = -np.log(p)
    F -= np.nanmin(F[np.isfinite(F)])
    return c, F, tot


def run(a):
    d = os.path.join(ROOT, a.name); os.makedirs(d, exist_ok=True)
    lo, hi, st_ = map(float, a.q0.split(':'))
    qs = np.round(np.arange(lo, hi + 1e-9, st_), 4)
    jobs = [(a.n, P, float(q0), a.kappa, a.sweeps, 1000 * s + i) for P in a.P for i, q0 in enumerate(qs) for s in range(a.seeds)]
    print(f'{a.name}: {len(jobs)} windows', flush=True)
    t0 = time.time()
    with open(f'{d}/windows.jsonl', 'a') as f, ProcessPoolExecutor(a.procs) as ex:
        for r in ex.map(window, jobs):
            f.write(json.dumps(r) + '\n')
    print(f'done in {time.time() - t0:.0f}s', flush=True)
    report(a)


def report(a):
    d = os.path.join(ROOT, a.name)
    R = [json.loads(l) for l in open(f'{d}/windows.jsonl')]
    byP = collections.defaultdict(list)
    for r in R:
        byP[r['P']].append(r)
    for P, W in sorted(byP.items()):
        n = W[0]['n']
        k = math.ceil(math.sqrt(n) - 1e-12)
        # equilibration check: per window, seeds' means
        g = collections.defaultdict(list)
        for w in W:
            g[w['q0']].append(np.mean(w['q']))
        bad = [q0 for q0, m in g.items() if len(m) > 1 and max(m) - min(m) > 0.1]
        c, F, tot = wham(W, n, np.linspace(-1, 1, 101))
        S = collections.defaultdict(list)
        for w in W:
            idx = np.clip(np.digitize(w['q'], np.linspace(-1, 1, 101)) - 1, 0, 99)
            for i, s in zip(idx, w['s']):
                S[i].append(s)
        print(f'\nP = {P}: {len(W)} windows; seeds disagree (> 0.1) at q0 = {sorted(bad)}')
        print('   Q4      F(kT)   <s>      samples')
        for i in range(0, 100, 4):
            if tot[i] > 0:
                print(f'   {c[i]:+.2f}  {F[i]:7.2f}  {np.mean(S[i]):8.4f}  {int(tot[i])}')
        mn = c[np.nanargmin(np.where(np.isfinite(F), F, np.nan))]
        print(f'   F minimum at Q4 = {mn:+.2f}; F(-0.2) - F(min) = '
              f'{F[np.argmin(abs(c + 0.2))]:.1f} kT (the 110 record has Q4 ~ {a.ref_q4 if hasattr(a, "ref_q4") else "?"})')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('run'); p.add_argument('name'); p.add_argument('--n', type=int, default=110)
    p.add_argument('--P', type=float, nargs='+', default=[30, 100, 300, 1000]); p.add_argument('--q0', default='-0.9:1.0:0.1')
    p.add_argument('--kappa', type=float, default=4.0); p.add_argument('--sweeps', type=int, default=40000)
    p.add_argument('--seeds', type=int, default=2); p.add_argument('--procs', type=int, default=16)
    p = sp.add_parser('report'); p.add_argument('name')
    a = ap.parse_args()
    {'run': run, 'report': report}[a.cmd](a)


if __name__ == '__main__':
    main()
