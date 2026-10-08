#!/usr/bin/env python3
"""Certify an explore.py archive: every basin below k through exactsolve (census_exact.solve); distinct certified minima
keyed by S_exact (1e-20), compared with the known set (runs/known<n>.json, 1e-9), role groups.  Cached per basin.

  explore_exact.py runs/ex0 [--procs 15] [--k 11]
"""
import argparse, collections, json, os
from jobpool import run_jobs
import census_exact

ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--procs', type=int, default=15)
ap.add_argument('--k', type=float, default=11); ap.add_argument('--n', type=int, default=110)
a = ap.parse_args()
E = [json.loads(l) for l in open(f'{a.run}/archive.jsonl')]
sub = [e for e in E if e['s'] < a.k]
cache = f'{a.run}/exact.json'
C = json.load(open(cache)) if os.path.exists(cache) else {}
todo = [e for e in sub if str(e['i']) not in C]
os.makedirs(f'{a.run}/exact', exist_ok=True)
res = run_jobs(census_exact.solve, [(e['path'], f'{a.run}/exact') for e in todo], procs=a.procs, timeout=900,
               on_timeout=lambda j: dict(cls='unresolved', status='killed'), on_error=lambda j, x: dict(cls='unresolved', status=str(x)[:60]))
for e, r in zip(todo, res):
    C[str(e['i'])] = r or dict(cls='unresolved', status='none')
json.dump(C, open(cache, 'w'), indent=0)
known = json.load(open(f'runs/known{a.n}.json'))
cls = collections.Counter(C[str(e['i'])]['cls'] for e in sub)
cert = {}
for e in sub:
    r = C[str(e['i'])]
    if r['cls'] == 'certified':
        cert.setdefault(r['key'], []).append(e)
new = {k: v for k, v in cert.items() if not any(abs(float(k) - x) < 2e-9 for x in known)}
print(f'{a.run}: {len(sub)} sub-{a.k:g} archive basins: {dict(cls)}')
print(f'  distinct certified minima: {len(cert)}; not in the known {len(known)}: {len(new)}; '
      f'union with known: {len(set(round(float(k), 9) for k in cert) | set(known))}')
print('  certified role groups (L, B, axis, other):', dict(collections.Counter(tuple(v[0]['roles']) for v in cert.values()).most_common()))
print('  new certified, lowest 10:', [f'{float(k):.10f} {tuple(v[0]["roles"])}' for k, v in sorted(new.items())[:10]])
print('  archive entries per certified class (f64 duplicates):', collections.Counter(len(v) for v in cert.values()))
