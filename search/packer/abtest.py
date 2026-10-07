#!/usr/bin/env python3
"""A/B test of single moves (10-06): from fixed minima, apply one move of a variant, quench, record where it lands.
Judged by where moves lead (full lines, below k, new vs same minimum, side change), not by acceptance.

  abtest.py --starts a.txt b.txt --k 11 --trials 40 --out runs/ab1 --procs 8
"""
import argparse, collections, json, math, os, random, statistics, tempfile, types
from jobpool import run_jobs
import mcmin

VARIANTS = {
    'clear':        dict(kind='reinsert', reinsert='clear', remove='uniform'),
    'clear_loose':  dict(kind='reinsert', reinsert='clear', remove='loose'),
    'pose0':        dict(kind='reinsert', reinsert='pose', tau=0.0, remove='uniform'),
    'pose0_loose':  dict(kind='reinsert', reinsert='pose', tau=0.0, remove='loose'),
    'poseT':        dict(kind='reinsert', reinsert='pose', tau=0.05, remove='uniform'),
    'poseT_loose':  dict(kind='reinsert', reinsert='pose', tau=0.05, remove='loose'),
    'pose0_l05':    dict(kind='reinsert', reinsert='pose', tau=0.0, remove='uniform', loosen=1.005),
    'pose0_l10':    dict(kind='reinsert', reinsert='pose', tau=0.0, remove='uniform', loosen=1.01),
    'pose0_loose_l05': dict(kind='reinsert', reinsert='pose', tau=0.0, remove='loose', loosen=1.005),
    'kick':         dict(kind='kick'),
    'kicksym':      dict(kind='kicksym'),
}


def job(j):
    start, var, t, k = j
    v = VARIANTS[var]
    a = types.SimpleNamespace(alpha=2.0, reinsert=v.get('reinsert', 'clear'), tau=v.get('tau', 0.0), remove=v.get('remove', 'uniform'), k=k)
    rng = random.Random(hash((start, var, t)) & 0xffffffff)
    s0, sq = mcmin.load_deg(start)
    new, desc = mcmin.MOVES[v['kind']](s0, sq, rng, a)
    r = mcmin.quench(s0, new, k, tempfile.mkdtemp(), loosen=v.get('loosen', 1.02))
    if r is None:
        return dict(start=start, var=var, t=t, desc=desc, status='quench failed')
    s, _, info = r
    return dict(start=start, var=var, t=t, desc=desc, s0=s0, s=s, lines=info['lines'], status='ok')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--starts', nargs='+', required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--trials', type=int, default=40); ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--variants', default=','.join(VARIANTS)); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    jobs = [(st, var, t, a.k) for t in range(a.trials) for st in a.starts for var in a.variants.split(',')]
    res = run_jobs(job, jobs, procs=a.procs, timeout=600, on_timeout=lambda j: dict(start=j[0], var=j[1], t=j[2], status='killed'),
                   on_error=lambda j, e: dict(start=j[0], var=j[1], t=j[2], status=f'failed {e}'))
    res = [r or dict(start=j[0], var=j[1], t=j[2], status='failed') for r, j in zip(res, jobs)]
    json.dump(res, open(a.out + '.json', 'w'), indent=0)
    print('start / variant: n, lines, same, new below k, new no-line, median ds (no-line, new), best')
    for st in a.starts:
        print(st)
        for var in a.variants.split(','):
            R = [r for r in res if r['start'] == st and r['var'] == var and r['status'] == 'ok']
            if not R:
                continue
            same = [r for r in R if abs(r['s'] - r['s0']) < 1e-8 and not r['lines']]
            new = [r for r in R if not r['lines'] and abs(r['s'] - r['s0']) >= 1e-8]
            ds = statistics.median([r['s'] - r['s0'] for r in new]) if new else float('nan')
            print(f'  {var:12s} {len(R):3d}  lines {sum(bool(r["lines"]) for r in R) / len(R):.2f}  same {len(same) / len(R):.2f}  '
                  f'newsub{a.k} {sum(r["s"] < a.k for r in new) / len(R):.2f}  new {len(new) / len(R):.2f}  ds {ds:+.4f}  '
                  f'best {min((r["s"] for r in new), default=float("nan")):.6f}')
