#!/usr/bin/env python3
"""Per-lineage report of an explore.py run: for each start (root), basins, below-k basins, certified (if exact.json),
best side, role groups of below-k basins, and CPU spent on parents of that lineage.

  lineage.py runs/ex7 [--k 11]
"""
import argparse, collections, json, os
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--k', type=float, default=11); a = ap.parse_args()
E = {e['i']: e for e in map(json.loads, open(f'{a.run}/archive.jsonl'))}
C = json.load(open(f'{a.run}/exact.json')) if os.path.exists(f'{a.run}/exact.json') else {}
def root(i):
    e = E[i]
    while e['parent'] >= 0: e = E[e['parent']]
    return e['i']
P = [json.loads(l) for l in open(f'{a.run}/proposals.jsonl')]
sec = collections.Counter()
for p in P: sec[root(p['par'])] += p['sec']
by = collections.defaultdict(list)
for i in E: by[root(i)].append(E[i])
print(f"{'root':38s} {'s0':>9s} {'CPU-h':>6s} {'basins':>6s} {'<k':>4s} {'cert':>4s} {'best':>12s}  roles of <k")
for r, L in sorted(by.items(), key=lambda t: min(e['s'] for e in t[1])):
    sub = [e for e in L if e['s'] < a.k]
    cert = sum(1 for e in sub if C.get(str(e['i']), {}).get('cls') == 'certified')
    print(f"{E[r]['kind'][6:44]:38s} {E[r]['s']:9.5f} {sec[r] / 3600:6.2f} {len(L):6d} {len(sub):4d} {cert:4d} {min(e['s'] for e in L):12.8f}  "
          f"{dict(collections.Counter(tuple(e['roles']) for e in sub).most_common(4))}")
