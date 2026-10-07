#!/usr/bin/env python3
"""Rediscovery calibration (10-06): from an older record, does our search find the newer one, and how does that depend on
the improvement delta and the distance between the two packings?

Pairs: every n where the live register (../exact/batch/candidates_all.json, 10-05) beats the 09-30 site mirror
(../../site/www/data/p/square-<n>.json; mostly Couzo, Sep 2026), plus our sw2 finds (266, 270: register -> ours).
Distance: min over the 8 symmetries of the container of the optimal matching (Hungarian) cost
sqrt(mean_i (dx^2 + dy^2 + (H sqrt2 * dtheta)^2)), dtheta folded mod 90 deg (corner-displacement metric).
Protocol (frozen): kicks from the old packing, sigma in {0.003, 0.01, 0.03, 0.1} (angles 20 sigma deg), mcmin.quench.
Detected: final side <= new + 1e-7 and no grid obstruction.

  rediscover.py --dist-only            # distances and deltas
  rediscover.py --trials 10 --procs 15 --out runs/rd1
"""
import argparse, collections, json, math, os, random, tempfile
import numpy as np
from scipy.optimize import linear_sum_assignment
from jobpool import run_jobs
import mcmin

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.join(HERE, '../../site/www/data/p')
BATCH = os.path.join(HERE, '../exact/batch')
SIGMAS = (0.003, 0.01, 0.03, 0.1)


def load_site(n):
    d = json.load(open(f'{SITE}/square-{n}.json'))
    return float(d['s']), [(float(x), float(y), float(a)) for x, y, a in d['squares']]


def pairs():
    reg = {d['n']: d for d in json.load(open(f'{BATCH}/candidates_all.json'))}
    out = []
    for n in range(1, 325):
        f = f'{SITE}/square-{n}.json'
        if not os.path.exists(f) or n not in reg or not os.path.exists(f'{BATCH}/inputs/n-{n}.txt'):
            continue
        s_old, old = load_site(n)
        s_new, new = mcmin.load_deg(f'{BATCH}/inputs/n-{n}.txt')
        if s_new < s_old - 1e-9:
            out.append(dict(n=n, label=f'{n}:mirror>register', s_old=s_old, s_new=s_new, old=old, new=new))
    for n, f in ((266, 'candidates/sw2_n266.exact.txt'), (270, 'candidates/sw2_n270.exact.txt')):
        s_old, old = mcmin.load_deg(f'{BATCH}/inputs/n-{n}.txt')
        s_new, new = mcmin.load_deg(os.path.join(HERE, f))
        out.append(dict(n=n, label=f'{n}:register>ours', s_old=s_old, s_new=s_new, old=old, new=new))
    return out


def sym(sq, s, k):
    """k-th symmetry of the square container (rotations by 90 deg, optionally mirrored); angles in degrees mod 90."""
    out = []
    for x, y, a in sq:
        if k >= 4:
            x, a = s - x, -a
        for _ in range(k % 4):
            x, y, a = s - y, x, a + 90
        out.append((x, y, a % 90))
    return out


def distance(old, s_old, new, s_new):
    f = s_new / s_old
    A = np.array([(x * f, y * f, a) for x, y, a in old])
    best = float('inf')
    for k in range(8):
        B = np.array(sym(new, s_new, k))
        dx = A[:, None, 0] - B[None, :, 0]; dy = A[:, None, 1] - B[None, :, 1]
        da = (A[:, None, 2] - B[None, :, 2] + 45) % 90 - 45
        C = dx * dx + dy * dy + (0.5 * math.sqrt(2) * np.radians(da)) ** 2
        r, c = linear_sum_assignment(C)
        best = min(best, math.sqrt(C[r, c].mean()))
    return best


def job(j):
    p, sig, t, out = j
    rng = random.Random(p['n'] * 100003 + int(sig * 1e5) * 101 + t)
    k = math.ceil(math.sqrt(p['n']) - 1e-12)
    sq = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), a + rng.gauss(0, 20 * sig)) for x, y, a in p['old']]
    r = mcmin.quench(p['s_old'], sq, k, tempfile.mkdtemp())
    if r is None:
        return dict(label=p['label'], sig=sig, t=t, status='quench failed')
    s, sq2, info = r
    return dict(label=p['label'], sig=sig, t=t, s=s, lines=info['lines'], status='ok')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--dist-only', action='store_true'); ap.add_argument('--trials', type=int, default=10)
    ap.add_argument('--procs', type=int, default=15); ap.add_argument('--out', default='runs/rd1')
    a = ap.parse_args()
    P = pairs()
    for p in P:
        p['delta'] = p['s_old'] - p['s_new']
        p['dist'] = distance(p['old'], p['s_old'], p['new'], p['s_new'])
    if a.dist_only:
        for p in sorted(P, key=lambda p: p['dist']):
            print(f"{p['label']:22s} delta {p['delta']:.2e}  dist {p['dist']:.4f}")
        raise SystemExit
    os.makedirs(a.out, exist_ok=True)
    jobs = [(p, sg, t, a.out) for p in sorted(P, key=lambda p: -p['n']) for sg in SIGMAS for t in range(a.trials)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=900, on_timeout=lambda j: dict(label=j[0]['label'], sig=j[1], t=j[2], status='killed'),
                   on_error=lambda j, e: dict(label=j[0]['label'], sig=j[1], t=j[2], status=f'failed {e}'))
    res = [r or dict(label=j[0]['label'], sig=j[1], t=j[2], status='failed') for r, j in zip(res, jobs)]
    meta = [{k: v for k, v in p.items() if k not in ('old', 'new')} for p in P]
    json.dump(dict(pairs=meta, trials=res, sigmas=SIGMAS), open(f'{a.out}/rd.json', 'w'), indent=0)
    print('label  delta  dist  | detected per sigma (of trials)  | improved-over-old per sigma')
    for p in sorted(P, key=lambda p: p['dist']):
        row, imp = [], []
        for sg in SIGMAS:
            R = [r for r in res if r['label'] == p['label'] and r['sig'] == sg and r['status'] == 'ok']
            row.append(sum(r['s'] <= p['s_new'] + 1e-7 and not r['lines'] for r in R))
            imp.append(sum(r['s'] < p['s_old'] - 1e-9 and not r['lines'] for r in R))
        print(f"{p['label']:22s} {p['delta']:.2e} {p['dist']:.4f} | {row} | {imp}", flush=True)
