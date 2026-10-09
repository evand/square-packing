#!/usr/bin/env python3
"""Settling trajectories for the s(110) explainer (10-08): slp2 squeeze from the record loosened by 3 % and nudged
(Gaussian sigma in position, 20 sigma degrees in angle); optional random starts (slow: slp2 from a 16 x 16 box).
Every accepted slp2 state is recorded (slp2.repair is wrapped; states that do not lower the side are dropped).

  settle.py --out runs/landscape110/settle --procs 6
"""
import argparse, json, math, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import slp2, morph
from jobpool import run_jobs

REC = os.path.join(os.path.dirname(HERE), 'runs/ex0/b00000.txt')


def traj(s, sq, budget):
    rec = []
    orig = slp2.repair

    def rep(s_, q):
        out = orig(s_, q)
        if math.isfinite(out[0]) and (not rec or out[0] < rec[-1][0] - 1e-12):
            rec.append((out[0], out[1]))
        return out
    slp2.repair = rep
    try:
        s2, sq2, ds0 = slp2.slp2(s, sq, R=1e-2, iters=100000, budget=budget)
    finally:
        slp2.repair = orig
    return s2, sq2, ds0, rec


def job(a):
    kind, seed, budget, sig = a
    rng = random.Random(seed)
    if kind == 'kick':                              # loosen 3 % (gaps everywhere), nudge, squeeze
        s, sq = morph.load(REC)
        f = 1.03
        s = s * f
        sq = [(x * f + rng.gauss(0, sig), y * f + rng.gauss(0, sig), t + math.radians(rng.gauss(0, 20 * sig))) for x, y, t in sq]
    else:                                           # random non-overlapping start in a 16 x 16 box
        s, sq = 16.0, []
        while len(sq) < 110:
            q = (rng.uniform(0.75, 15.25), rng.uniform(0.75, 15.25), rng.uniform(0, math.pi / 2))
            if all((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 > 2.0 for p in sq):
                sq.append(q)
    s2, sq2, ds0, rec = traj(s, sq, budget)
    return dict(kind=kind, seed=seed, sigma=sig, s_final=s2, jammed=bool(ds0 is not None and ds0 > -1e-8), n=len(rec),
                traj=[(round(s_, 6), [(round(x, 4), round(y, 4), round(t, 4)) for x, y, t in q]) for s_, q in rec])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--procs', type=int, default=6)
    ap.add_argument('--kick', type=int, default=6); ap.add_argument('--rand', type=int, default=4)
    ap.add_argument('--budget', type=float, default=900)
    ap.add_argument('--sigmas', type=float, nargs='+', default=[0.01, 0.03, 0.1])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = [('kick', 1000 * i + k, a.budget, sg) for i, sg in enumerate(a.sigmas) for k in range(a.kick)] + \
           [('rand', 100 + k, a.budget, 0) for k in range(a.rand)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=a.budget + 300, on_timeout=lambda j: None, on_error=lambda j, x: dict(err=str(x)[:200]))
    for (kind, seed, _, sg), r in zip(jobs, res):
        if r and 'traj' in r:
            json.dump(r, open(f'{a.out}/{kind}_{seed}.json', 'w'))
            print(kind, seed, sg, f"final {r['s_final']:.7f} jammed {r['jammed']} states {r['n']}", flush=True)
        else:
            print(kind, seed, 'failed', r)
