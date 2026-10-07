#!/usr/bin/env python3
"""Square roles per census trial (axis / left 10-45 deg / bottom 45-80 deg mod 90 / other) vs the start, by index.
Caveat: the 45 deg cut splits nothing real at n = 110 but junction squares rotating to ~46 deg read as L>B.  Run from packer/."""
import json, math, collections, sys
sys.path.insert(0, '.')
def load(p):
    L = open(p).read().split('\n')
    n = int(L[0].split()[0])
    return [tuple(float(v) for v in l.split()[:3]) for l in L[1:n+1]]
def role(deg):
    a = deg % 90
    d = min(a, 90 - a)
    if d < 3: return 'A'          # axis (incl. junction squares tilted < 3 deg)
    if 10 < a < 45: return 'L'
    if 45 <= a < 80: return 'B'
    return 'O'
def roles(p): return [role(t) for _, _, t in load(p)]
ref = {'rec': roles('seeds/rec110.txt'), 'cand': roles('candidates/s110_10.996783403.txt')}
print('ref counts', {k: dict(collections.Counter(v)) for k, v in ref.items()})
print('rec vs cand role diffs:', [(i, a, b) for i, (a, b) in enumerate(zip(ref['rec'], ref['cand'])) if a != b])
T = json.load(open('runs/cen7/census_exact.json'))
rows = []
for r in T:
    if 's' not in r: continue
    R = roles(r['path']); R0 = ref[r['run']]
    c = collections.Counter(R)
    moves = collections.Counter(f'{a}>{b}' for a, b in zip(R0, R) if a != b)
    cls = r['key'] if r['sub'] and r['ex']['cls'] == 'certified' else ('sub-other' if r['sub'] else ('line' if r.get('lines') else '>=11'))
    rows.append((cls, r['sig'], (c['L'], c['B'], c['A'], c['O']), tuple(sorted(moves.items()))))
# per certified class: count tuples and moves
print('\n== certified sub-11 classes: (L, B, axis, other) counts and role moves vs start')
by = collections.defaultdict(list)
for cls, sg, cnt, mv in rows: by[cls].append((cnt, mv, sg))
for cls in sorted(k for k in by if k.startswith('10.')):
    v = by[cls]
    print(cls, len(v), dict(collections.Counter(x[0] for x in v)), 'moves:', dict(collections.Counter(x[1] for x in v).most_common(3)))
print('\n== other outcomes by sigma: count tuples (top) and any-move fraction')
for cls in ('sub-other', '>=11', 'line'):
    v = by[cls]
    print(cls, len(v), dict(collections.Counter(x[0] for x in v).most_common(6)))
print('\n== role moves by sigma (all trials): fraction with any move; move types')
for sg in sorted({x[1] for x in rows}):
    R = [x for x in rows if x[1] == sg]
    mv = collections.Counter(k for x in R for k, _ in x[3])
    print(sg, f'{sum(bool(x[3]) for x in R)}/{len(R)} trials with a role change;', dict(mv))
