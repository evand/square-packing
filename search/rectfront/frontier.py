#!/usr/bin/env python3
"""Numerical frontier of R_n = {(w, h): n unit squares fit in w x h}, small n (2026-10-20).

For each n and dx >= 0, minimise h with the container [0, h + dx] x [0, h] (fq rectangle mode), from random
starts (rectsearch.random_start: random / tilt-alphabet / diagonal band) plus kicks of the best.  Compare with the
grid frontier G_n(w) = min over column counts c <= w of ceil(n / c) (axis-parallel rows).  Output JSONL per (n, dx).

  python3 frontier.py --n 2-8 --dx 0:3:0.05 --starts 40 --procs 4 --out runs/rectfront/f1
"""
import argparse, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
PK = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'packer')
sys.path.insert(0, PK)


def grid_h(n, dx):
    """Least h with an axis-parallel grid of r = h rows, ceil(n/r) columns in width h + dx (axis-parallel squares
    in w x h hold floor(w) floor(h), so this is also the axis-parallel optimum on the line w = h + dx)."""
    return min(max(r, math.ceil(n / r) - dx) for r in range(1, n + 1))


def start(n, dx, rng):
    """rectsearch.random_start with h >= 1.05 (it needs room for centres in both directions)."""
    from rectsearch import random_start, h_guess
    loose = 1.6
    while h_guess(n, dx, loose) < 1.05:
        loose *= 1.2
    return random_start(n, dx, rng, loose)


def job(a):
    n, dx, starts, kicks, seed = a
    from rectsearch import random_start
    from pk.engine import FQ, QOpt
    from pk.packing import Packing
    rng = np.random.default_rng(seed)
    fq = FQ(); Q = QOpt(loosen=1.0, rect_dx=float(dx), flip_top=16)
    best = (math.inf, None); allh = []
    for k in range(starts + kicks):
        if k < starts or best[1] is None:
            h0, sq = start(n, dx, rng)
        else:
            sig = [0.02, 0.05, 0.15][k % 3]
            sq = best[1].copy(); sq[:, :2] += rng.normal(0, sig, (n, 2)); sq[:, 2] += rng.normal(0, 20 * sig, n)
            h0 = best[0] * (1 + sig / 5)
        r = fq.quench(Packing(float(h0), sq), Q, timeout=120)
        if not r.ok:
            continue
        allh.append(r.s)
        if r.s < best[0] - 1e-12:
            best = (r.s, np.array(r.p.sq, float))
    fq.close()
    h = best[0]
    return dict(n=n, dx=dx, w=h + dx, h=h, grid=grid_h(n, dx), nq=len(allh),
                hits=sum(1 for x in allh if x < h + 1e-7), sq=None if best[1] is None else best[1].tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', default='2-8'); ap.add_argument('--dx', default='0:3:0.05')
    ap.add_argument('--starts', type=int, default=40); ap.add_argument('--kicks', type=int, default=40)
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    lo, hi = map(int, a.n.split('-')); d0, d1, dd = map(float, a.dx.split(':'))
    if d1 < 0: d1 = hi - 1
    dxs = [round(d0 + i * dd, 6) for i in range(int(round((d1 - d0) / dd)) + 1)]  # d1 < 0: up to n - 1
    os.makedirs(a.out, exist_ok=True)
    jobs = [(n, dx, a.starts, a.kicks, 1000 * n + i) for n in range(lo, hi + 1) for i, dx in enumerate(dxs)
            if d1 >= 0 or dx <= n - 1 + 1e-9]
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'front.jsonl'), 'a') as f:
        for r in ex.map(job, jobs):
            f.write(json.dumps(r) + '\n'); f.flush()
            tag = '  <-- below grid' if r['h'] < r['grid'] - 1e-6 else ''
            print(f"n={r['n']} dx={r['dx']:.2f} w={r['w']:.6f} h={r['h']:.6f} grid={r['grid']} hits={r['hits']}/{r['nq']}{tag}", flush=True)


if __name__ == '__main__':
    main()
