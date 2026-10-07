#!/usr/bin/env python3
"""Count-pair beam: pick the best packing per role-count key (L, B, axis, other) from movegen runs, then run movegen from them.

  beam.py runs/mv2/m.json [more.json ...] --top 12 --below 11.01 --trials 40 --out runs/mv3/m [--seed 3]
Keys are recounted from the saved packings with movegen.role (so older runs use the current role windows).
"""
import argparse, collections, json, math, subprocess, sys
from movegen import role
from ftmc import load

ap = argparse.ArgumentParser()
ap.add_argument('runs', nargs='+'); ap.add_argument('--top', type=int, default=12); ap.add_argument('--below', type=float, default=11.01)
ap.add_argument('--trials', type=int, default=40); ap.add_argument('--out', required=True); ap.add_argument('--seed', type=int, default=3)
ap.add_argument('--n', default='110'); ap.add_argument('--k', default='11'); ap.add_argument('--dry', action='store_true'); ap.add_argument('--procs', default='15'); ap.add_argument('--focus', default=None)
a = ap.parse_args()
best = {}
for p in a.runs:
    for r in json.load(open(p))['trials']:
        if r.get('status') != 'ok' or r['lines'] or r['s'] >= a.below:
            continue
        _, C = load(r['path'])
        c = collections.Counter(role(q[2]) for q in C)
        key = (c['L'], c['B'], c['A'] + c['J'], c['O'])
        if key not in best or r['s'] < best[key][0]:
            best[key] = (r['s'], r['path'])
sel = sorted(best.items(), key=lambda kv: kv[1][0])[:a.top]
for key, (s, path) in sel:
    print(f'{key}  {s:.7f}  {path}')
sys.stdout.flush()
if not a.dry:
    subprocess.run([sys.executable, 'movegen.py', *[p for _, (_, p) in sel], '--n', a.n, '--k', a.k, '--trials', str(a.trials),
                    '--out', a.out, '--seed', str(a.seed), '--procs', a.procs] + (['--focus', a.focus] if a.focus else []))
