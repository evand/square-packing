#!/usr/bin/env python3
"""Per-move-kind yield of an explore.py run: proposals, CPU seconds, outcome mix, and new certified sub-k minima (certified
by explore_exact.py, not in a reference set, credited to the kind that first produced the archive entry) per CPU hour.

  arm_yield.py runs/ex2 [--ref runs/known110_all.json] [--k 11]
"""
import argparse, collections, json
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--ref', default='runs/known110_all.json')
ap.add_argument('--k', type=float, default=11); a = ap.parse_args()
ref = json.load(open(a.ref))
import bisect
def known(s):
    i = bisect.bisect_left(ref, s - 2e-9)
    return i < len(ref) and abs(ref[i] - s) < 2e-9
E = {e['i']: e for e in map(json.loads, open(f'{a.run}/archive.jsonl'))}
C = json.load(open(f'{a.run}/exact.json'))
P = [json.loads(l) for l in open(f'{a.run}/proposals.jsonl')]
seen = set(); credit = collections.Counter(); credit_all = collections.Counter()
for p in P:
    if p['new'] is None: continue
    e = E[p['new']]; r = C.get(str(p['new']))
    if e['s'] < a.k and r and r['cls'] == 'certified':
        key = round(float(r['key']), 9)
        if key in seen: continue
        seen.add(key)
        credit_all[p['kind']] += 1
        if not known(float(r['key'])): credit[p['kind']] += 1
by = collections.defaultdict(list)
for p in P: by[p['kind']].append(p)
tot_new = sum(credit.values()); tot_sec = sum(p['sec'] for p in P)
print(f"{a.run}: {len(P)} proposals, {tot_sec / 3600:.2f} CPU-h; distinct certified sub-{a.k:g}: {len(seen)}, new vs ref ({len(ref)}): {tot_new}")
print(f"{'kind':11s} {'props':>6s} {'CPU-h':>6s} {'sec/p':>6s} {'ret%':>5s} {'disc%':>6s} {'cert':>5s} {'new':>5s} {'new/CPU-h':>9s}")
for k in sorted(by, key=lambda k: -credit[k] / max(1e-9, sum(p['sec'] for p in by[k]))):
    L = by[k]; sec = sum(p['sec'] for p in L)
    ret = sum(p['st'] == 'return' for p in L) / len(L); disc = sum(p['st'] in ('screen-grid', 'screen-discard', 'discard') or (p['st'] == 'new?' and p['new'] is None) for p in L) / len(L)
    print(f"{k:11s} {len(L):6d} {sec / 3600:6.2f} {sec / len(L):6.2f} {100 * ret:5.0f} {100 * disc:6.0f} {credit_all[k]:5d} {credit[k]:5d} {credit[k] / (sec / 3600):9.1f}")
