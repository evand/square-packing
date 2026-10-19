#!/usr/bin/env python3
"""Random starts, free vs diagonal-mirror (Md) symmetric (10-18, stage 0).

free:  n poses uniform in [0.5, S-0.5]^2, angle uniform [0, 90).
sym:   m squares on the diagonal (x = y, angle 0 or 45 with p 1/2), m = n mod 2 + 2j (j random so m ~ the record's own
       on-diagonal count, passed in), the other (n - m) / 2 uniform in the lower triangle (y < x) + their images
       (x, y, th) -> (y, x, -th mod 90).
Both then go through the same fq quench (ALM + SLP polish); the quench itself is unconstrained (stage 0 tests whether
symmetric STARTS alone find the record funnel more often).
"""
import sys, os, math, json, argparse, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np


def free_start(n, S, rng):
    return np.c_[rng.uniform(0.5, S - 0.5, (n, 2)), rng.uniform(0, 90, n)]


def sym_start(n, S, m, rng, diag_mode='run'):
    """diag_mode 'iid': each diagonal square 0 / 45 with p 1/2 (smoke 10-18: worse than free everywhere);
    'run' (default): diagonal squares sorted, one contiguous 45-degree run of random start / length, 0 elsewhere (the
    pattern every Md record with tilts has: Goebel strip across the diagonal).  m = None: drawn from
    ceil(S) - 1 .. ceil(S) + 2 (no record information)."""
    if m is None:
        m = int(rng.integers(math.ceil(S) - 1, math.ceil(S) + 3))
    m = max(0, min(n, m)); m += (n - m) % 2
    d = np.sort(rng.uniform(0.5, S - 0.5, m))
    if diag_mode == 'iid':
        ang = np.where(rng.random(m) < 0.5, 0.0, 45.0)
    else:
        L = int(rng.integers(0, m + 1)); a0 = int(rng.integers(0, m - L + 1))
        ang = np.zeros(m); ang[a0:a0 + L] = 45.0
    diag = np.c_[d, d, ang]
    k = (n - m) // 2
    x = rng.uniform(0.5, S - 0.5, (k, 2)); x.sort(axis=1); x = x[:, ::-1]        # y < x
    th = rng.uniform(0, 90, k)
    half = np.c_[x, th]; img = np.c_[x[:, 1], x[:, 0], (-th) % 90]
    return np.r_[diag, half, img]


def ondiag(sq, tol=1e-6):
    return int((np.abs(sq[:, 0] - sq[:, 1]) < tol).sum())


def job(args):
    from pk.engine import fq_quench, QOpt
    from pk.packing import read
    arm, n, S, m, seed, loosen = args
    rng = np.random.default_rng(seed)
    sq = free_start(n, S, rng) if arm == 'free' else sym_start(n, S, m, rng)
    t0 = time.process_time()
    with tempfile.TemporaryDirectory() as tmp:
        run = fq_quench((S, sq), QOpt(loosen=loosen), tmp=tmp, timeout=600)
        p = read(run.out) if run.out else None
    if p is None:
        return dict(arm=arm, n=n, seed=seed, s=None)
    return dict(arm=arm, n=n, seed=seed, s=float(p.s), pen=float(p.max_pen()), sq=p.sq.tolist(),
                cpu=time.process_time() - t0)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, default=8)
    ap.add_argument('--loosen', type=float, default=1.0); ap.add_argument('--S', type=float, default=None)
    a = ap.parse_args()
    from multiprocessing import Pool
    sc = json.load(open('runs/symsub/symscan.json'))[str(a.n)]
    from pk.store import Store
    rec = Store().get(sc['pid'])
    m = ondiag(rec.sq); S = a.S or rec.s * 1.05
    T = [(arm, a.n, S, m, 1000 + i, a.loosen) for i in range(a.k) for arm in ('free', 'sym')]
    t0 = time.time()
    with Pool(min(15, len(T))) as pool:
        R = pool.map(job, T)
    for arm in ('free', 'sym'):
        ss = sorted(r['s'] for r in R if r['arm'] == arm and r['s'])
        print(arm, 'record', round(rec.s, 6), 'sides', [round(x, 4) for x in ss])
    print(f'{time.time()-t0:.0f}s wall, ondiag {m}')
