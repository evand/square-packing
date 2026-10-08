#!/usr/bin/env python3
"""Fast crop seeds (10-08): crop.py's strip cut (one unit row + column strip removed, rest shifted) on bigger records, count
fixed to n (crop.fix_count), fq quench (loosen 1.02/1.05), grid check.  Keeps distinct seeds below --keep.

  cropseed.py --src 129 130 131 132 --n 110 --trials 150 --procs 8 --out seeds110c
"""
import argparse, json, math, os, random, tempfile
from jobpool import run_jobs
import crop, hop, mcmin, explore


def job(j):
    src, t, n, k = j
    from layout import full_lines
    rng = random.Random(src * 10007 + t)
    s, sq = mcmin.load_deg(f'{hop.BATCH}/n-{src}.txt')
    cx, cy = rng.uniform(0, s - 1), rng.uniform(0, s - 1)
    C = crop.crop(sq, s, cx, cy)
    C, dr, ad = crop.fix_count(C, n, s - 1)
    hop.EXTRA[:] = ['--loosen', rng.choice(('1.02', '1.05'))]
    tmp = tempfile.mkdtemp()
    r = hop.quench(s - 1, C, tmp)
    if not r:
        return dict(src=src, t=t, s=None)
    p = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p, r[0], r[1])
    return dict(src=src, t=t, s=r[0], sq=r[1], grid=bool(full_lines(p, k)), dropped=dr, added=ad,
                roles=explore.roles(r[1]), cut=(round(cx, 2), round(cy, 2)))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', type=int, nargs='+'); ap.add_argument('--n', type=int, default=110)
    ap.add_argument('--trials', type=int, default=100); ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--keep', type=float, default=11.05); ap.add_argument('--out', default='seeds110c')
    a = ap.parse_args()
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    os.makedirs(a.out, exist_ok=True)
    R = [r for r in run_jobs(job, [(src, t, a.n, k) for src in a.src for t in range(a.trials)], procs=a.procs, timeout=600) if r and r.get('s')]
    ok = sorted([r for r in R if not r['grid'] and r['s'] < a.keep], key=lambda r: r['s'])
    seen = []
    for r in ok:
        if all(abs(r['s'] - x) > 1e-9 for x in seen):
            seen.append(r['s'])
            mcmin.write_deg(f"{a.out}/c{r['src']}_{r['t']}_{r['s']:.7f}.txt", r['s'], r['sq'])
    json.dump([{kk: v for kk, v in r.items() if kk != 'sq'} for r in R], open(f'{a.out}/crops.json', 'w'))
    import collections
    print(f'{len(R)} crops; grid {sum(r["grid"] for r in R)}; kept {len(seen)} distinct < {a.keep}; sub-{k}: {sum(x < k for x in seen)}')
    for src in a.src:
        L = [r for r in R if r['src'] == src]
        g = [r for r in L if not r['grid']]
        print(f'  src {src}: {len(L)} crops, {len(g)} non-grid, best {min([r["s"] for r in g] or [float("nan")]):.6f}, '
              f'sub-11 {sum(r["s"] < k for r in g)}, roles of sub-11: {dict(collections.Counter(tuple(r["roles"]) for r in g if r["s"] < k))}')
