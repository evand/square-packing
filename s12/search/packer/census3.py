#!/usr/bin/env python3
"""Landscape census: reference + Gaussian sigma -> soft penalty squeeze (loosen 2 %, mu0 1e3) -> slp2 finish.
Records final side, full lines, first-order jam, and a structure fingerprint (tilted count; angle-cluster sizes).
  census3.py ref.txt --k 11 --sigmas ... --trials 30 --out prefix
"""
import argparse, math, os, random, subprocess, tempfile, json, collections
from jobpool import run_jobs
from rigid import load
from slp2 import slp2
from slp import save
from layout import full_lines

JOB_BUDGET = 600          # s: slp2 wall budget per trial (soft); jobpool kills at 2x
PACKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'target/release/packer')


def fingerprint(sq):
    angs = sorted(math.degrees(t) % 90 for _, _, t in sq)
    tilt = [a for a in angs if 0.5 < a < 89.5]
    clusters = []
    for a in tilt:
        if clusters and a - clusters[-1][-1] < 2.0:
            clusters[-1].append(a)
        else:
            clusters.append([a])
    return len(tilt), tuple(sorted((len(c) for c in clusters), reverse=True))


def job(args):
    ref, sig, t, out, k, i = args
    s0, C = load(ref)
    rng = random.Random(7919 * t + int(sig * 1e7) + 1000003 * max(i, 0))       # i = ref index (0: as before)
    f = 1.02
    P = [(x * f + rng.gauss(0, sig), y * f + rng.gauss(0, sig), a + math.radians(rng.gauss(0, 20 * sig))) for x, y, a in C]
    tmp = tempfile.mkdtemp()
    p, q = os.path.join(tmp, 'p.txt'), os.path.join(tmp, 'q.txt')
    save(p, s0 * f + 4 * sig, P)
    subprocess.run([PACKER, 'relax', '--in', p, '--squeeze-pen', '--mu0', '1e3', '--out', q], capture_output=True, timeout=300)
    s1, sq = load(q)
    info = {}
    s, sq2, ds0 = slp2(s1, sq, R=1e-3, rmin=1e-8, budget=JOB_BUDGET, info=info)
    path = f'{out}_{sig:g}_{t}.txt' if i < 0 else f'{out}_r{i}_{sig:g}_{t}.txt'
    save(path, s, sq2)
    return dict(ref=ref, path=path, sig=sig, t=t, soft=s1, s=s, jam=info['jammed'], lines=full_lines(path, k), fp=fingerprint(sq2),
                timeout=info['timeout'], switch_steps=info['switch_steps'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('ref', nargs='+'); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--sigmas', default='1e-3,3e-3,1e-2,3e-2,0.1'); ap.add_argument('--trials', type=int, default=30)
    ap.add_argument('--procs', type=int, default=14); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    # one ref: file names and seeds as before (index -1 -> seed offset 0); several: names carry _r<i>
    refs = [(a.ref[0], -1)] if len(a.ref) == 1 else list(zip(a.ref, range(len(a.ref))))
    jobs = [(rf, float(sg), t, a.out, a.k, i) for rf, i in refs for sg in a.sigmas.split(',') for t in range(a.trials)]
    # killed / failed trials stay in the record (status), so census statistics count every trial
    res = run_jobs(job, jobs, procs=a.procs, timeout=2 * JOB_BUDGET,
                   on_timeout=lambda j: dict(ref=j[0], sig=j[1], t=j[2], status='killed'),
                   on_error=lambda j, e: dict(ref=j[0], sig=j[1], t=j[2], status=f'failed: {e}'))
    res = [r if r is not None else dict(ref=j[0], sig=j[1], t=j[2], status='failed') for r, j in zip(res, jobs)]
    json.dump(dict(ref=a.ref, k=a.k, sigmas=a.sigmas, trials=a.trials, budget=JOB_BUDGET, trials_out=res),
              open(a.out + '.json', 'w'))
    for sg in a.sigmas.split(','):
        R = [r for r in res if r['sig'] == float(sg) and 's' in r]
        nf = sum(1 for r in res if r['sig'] == float(sg) and 's' not in r)
        if not R:
            print(f'sigma={sg}: all {nf} trials failed', flush=True)
            continue
        if nf:
            print(f'sigma={sg}: {nf} trials failed/killed', flush=True)
        sides = collections.Counter(round(r['s'], 6) for r in R)
        fps = collections.Counter(r['fp'][0] for r in R)
        below = sum(r['s'] < a.k and r['lines'] == 0 for r in R)
        print(f"sigma={sg}: below k {below}/{len(R)}; distinct sides {len(sides)} (repeats {sum(v for v in sides.values() if v > 1)}); "
              f"jammed {sum(r['jam'] for r in R)}; with lines {sum(r['lines'] > 0 for r in R)}; best {min(sides):.6f}; "
              f"median {sorted(r['s'] for r in R)[len(R) // 2]:.6f}; tilted counts {dict(sorted(fps.items()))}", flush=True)
