import json, os
from mpmath import mp, mpf
mp.dps = 150
W = '/home/evand/math/square-packing/public/s12/search/exact/batch/work'
def ourS(n):
    for f in (f'{W}/n-{n}/witness.exact.txt', f'{W}/n-{n}/polished.exact.txt'):
        if os.path.exists(f): return mpf(open(f).readline().split()[1]), f.split('/')[-1]
    return None, None
rows = []
for r in json.load(open('kingbird.json')):
    if 'factors' not in r: continue
    n = r['n']; S, src = ourS(n)
    if S is None: rows.append((n, 'no ours')); continue
    best = None
    for c in r['factors']:
        x = S
        for _ in range(60):
            p = mpf(0); dp = mpf(0)
            for ci in c: dp = dp * x + p; p = p * x + ci
            if dp == 0: break
            x -= p / dp
        d = abs(x - S)
        if best is None or d < best[0]: best = (d, len(c) - 1)
    rows.append((n, best[1], mp.nstr(best[0], 3), src, r['num']))
for x in rows: print(*x)
json.dump([{'n': x[0], 'deg': x[1], 'delta': x[2]} for x in rows if len(x) > 2], open('cmp.json', 'w'))
