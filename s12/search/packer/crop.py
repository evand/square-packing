#!/usr/bin/env python3
"""Crop constructor: a packing at side ~k+1 -> seeds at side ~k (e.g. s(132) record -> n = 110, s(110) record -> n = 90).

Remove the squares whose centres lie in a horizontal strip [cy, cy+1) and a vertical strip [cx, cx+1); shift the squares
above / right of a strip by -1 (side s -> s - 1); fix the count to n (drop the most-overlapped squares, or insert axis
squares at least-overlap positions); penalty squeeze from the overlapping state at side (s-1)*loosen; reject full lines;
slp2 finish.  Cuts are drawn uniformly over the square (through bands, junction, or axis region); mu0 and loosen vary.

  crop.py src.txt --n 110 --k 11 --trials 300 --out runs/crop1/c [--procs 15] [--seed 1]
Output per trial: cut, removed count, count fix, soft side, final side, jam, full lines, role counts.
"""
import argparse, collections, json, math, os, random, tempfile
from ftmc import load, run_packer
from gen import pen, write
from layout import full_lines
from jobpool import run_jobs
from movegen import role

JOB_BUDGET = 300


def crop(C, s, cx, cy):
    out = []
    for x, y, a in C:
        if cy <= y < cy + 1 or cx <= x < cx + 1:
            continue
        out.append((x - (x >= cx + 1), y - (y >= cy + 1), a))
    return out


def fix_count(C, n, s):
    C = list(C)
    over = lambda q, S: sum(max(0.0, pen(q, o)) for o in S if o is not q)
    dropped = added = 0
    while len(C) > n:
        C.remove(max(C, key=lambda q: over(q, C)))
        dropped += 1
    while len(C) < n:
        grid = [(0.5 + i * (s - 1) / 60, 0.5 + j * (s - 1) / 60, 0.0) for i in range(61) for j in range(61)]
        C.append(min(grid, key=lambda q: over(q, C)))
        added += 1
    return C, dropped, added


def job(args):
    src, n, k, seed, out = args
    from slp2 import slp2
    from rigid import load as load_rad
    from slp import save
    rng = random.Random(seed)
    s0, C0 = load(src)
    cx, cy = rng.uniform(0, s0 - 1), rng.uniform(0, s0 - 1)
    mu0 = rng.choice([1e2, 1e3]); loosen = rng.choice([1.0, 1.01, 1.02])
    C = crop(C0, s0, cx, cy)
    s1 = s0 - 1
    rec = dict(seed=seed, cx=cx, cy=cy, mu0=mu0, loosen=loosen, kept=len(C))
    C, rec['dropped'], rec['added'] = fix_count(C, n, s1)
    tmp = tempfile.mkdtemp()
    a, b = os.path.join(tmp, 'c.txt'), os.path.join(tmp, 'cr.txt')
    s = s1 * loosen
    write(a, s, [(min(max(x, 0.5), s - 0.5), min(max(y, 0.5), s - 0.5), t) for x, y, t in C])
    _, ss = run_packer(['relax', '--in', a, '--squeeze-pen', '--mu0', repr(mu0), '--out', b])
    rec['soft'] = ss
    if ss == float('inf'):
        rec['status'] = 'squeeze failed'
        return rec
    if full_lines(b, k) > 0:
        rec['status'] = 'full line after squeeze'
        return rec
    s2, sq = load_rad(b)
    info = {}
    sf, sq2, _ = slp2(s2, sq, R=1e-3, rmin=1e-8, budget=JOB_BUDGET, info=info)
    fn = f'{out}_{seed}.txt'
    save(fn, sf, sq2)
    c = collections.Counter(role(math.degrees(q[2])) for q in sq2)
    rec.update(status='ok', s=sf, jam=info['jammed'], lines=full_lines(fn, k), path=fn,
               counts=dict(L=c['L'], B=c['B'], A=c['A'] + c['J'], O=c['O']))
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--trials', type=int, default=100); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    jobs = [(a.src, a.n, a.k, a.seed * 100000 + t, a.out) for t in range(a.trials)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=3 * JOB_BUDGET,
                   on_timeout=lambda j: dict(seed=j[3], status='killed'),
                   on_error=lambda j, e: dict(seed=j[3], status=f'failed: {e}'))
    res = [r if r is not None else dict(seed=j[3], status='failed') for r, j in zip(res, jobs)]
    json.dump(dict(src=a.src, n=a.n, k=a.k, trials=res), open(a.out + '.json', 'w'), indent=0)
    ok = [r for r in res if r.get('status') == 'ok' and not r['lines']]
    print(f'{len(res)} trials; statuses {dict(collections.Counter(r["status"][:30] for r in res))}')
    print(f'jammed, no full line: {len(ok)}; below k: {sum(r["s"] < a.k for r in ok)}')
    for r in sorted(ok, key=lambda r: r['s'])[:15]:
        print(f"  {r['s']:.7f}  cut ({r['cx']:.2f}, {r['cy']:.2f})  kept {r['kept']} -{r['dropped']} +{r['added']}  "
              f"mu0 {r['mu0']:g} loosen {r['loosen']}  counts {r['counts']}")
