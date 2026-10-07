#!/usr/bin/env python3
"""Coverage check of run V2 record (exact)."""
import gzip, json, hashlib, subprocess
from fractions import Fraction as F
from collections import Counter
P = '/home/evand/math/square-packing/public/s12/search/'
L = [json.loads(l) for l in gzip.open(P + 'qx2_data/cert/runV2_cdade4b6_full.jsonl.gz', 'rt')]
hdr = L[0]; recs = L[1:]
print('header', hdr['sha256'])
print('input sha == box file:', hdr['sha256']['input'] == hashlib.sha256(open(P + 'qx2_data/L4_k02_box7.txt', 'rb').read()).hexdigest())
code = subprocess.run(['git', '-C', P, 'show', '0e303ca:s12/search/qx2_zm.py'], capture_output=True).stdout
if not code:
    code = subprocess.run(['git', '-C', P, 'show', '0e303ca:search/qx2_zm.py'], capture_output=True).stdout
print('0e303ca qx2_zm.py sha', hashlib.sha256(code).hexdigest(), '== header:', hashlib.sha256(code).hexdigest() == hdr['sha256']['qx2_zm.py'])
for fn in ('zm_mixed.py', 'zeromargin.py', 'mixed_cover.py', 'qx2_zm.py'):
    for pre in ('s12/search/', 'search/'):
        c = subprocess.run(['git', '-C', P, 'show', '0e303ca:' + pre + fn], capture_output=True).stdout
        if c: break
    cur = hashlib.sha256(open(P + fn, 'rb').read()).hexdigest()
    print(f'  {fn}: at 0e303ca {hashlib.sha256(c).hexdigest()[:12]}  working tree {cur[:12]}  header {hdr["sha256"][fn][:12]}')
print(len(recs), 'records; keys', Counter(tuple(sorted(r)) for r in recs))
roots = [tuple(F(v) for v in r['root']) for r in recs]
print('distinct roots', len(set(roots)))
# tiling: roots must be boxes whose union is [0,7/2]^2 x [0,1/2] with disjoint interiors.
vol = sum((x1 - x0) * (y1 - y0) * (u1 - u0) for x0, x1, y0, y1, u0, u1 in roots)
print('total volume', vol, '== 49/4 * 1/2:', vol == F(49, 4) * F(1, 2))
xs = sorted(set([r[0] for r in roots] + [r[1] for r in roots])); ys = sorted(set([r[2] for r in roots] + [r[3] for r in roots]))
us = sorted(set([r[4] for r in roots] + [r[5] for r in roots]))
print('x breaks', xs[0], xs[-1], len(xs), ' y', ys[0], ys[-1], len(ys), ' u', us)
grid = set((xs[i], xs[i+1], ys[j], ys[j+1], us[l], us[l+1]) for i in range(len(xs)-1) for j in range(len(ys)-1) for l in range(len(us)-1))
print('roots == full product grid cells:', set(roots) == grid, len(grid))
# per-root stats
leafkinds = ['LEB', 'CAP', 'EXACT', 'EXACT0', 'EXACT45', 'AXIS', 'SYM', 'PIECE', 'ADM', 'P1', 'MIX', 'SPLIT', 'EMPTY', 'UNCERT']
tot = Counter(); bad = []
for r in recs:
    st = r['st']
    nleaf = sum(st.get(k, 0) for k in leafkinds)
    if st['UNCERT'] != 0 or r['unc']: bad.append(('unc', r['root']))
    if st['boxes'] != 2 * nleaf - 1: bad.append(('count', r['root'], st['boxes'], nleaf))
    if st['maxdepth'] > 18: bad.append(('depth', r['root']))
    for k, v in st.items():
        if k != 'cpu': tot[k] += v
print('bad records:', len(bad), bad[:5])
print('totals', dict(tot))
print('unknown stat keys:', set(tot) - set(leafkinds) - {'UNCERT', 'boxes', 'maxdepth', 'THR', 'LIN', 'CORE', 'CHAIN', 'TRI', 'TPTS'})
# where are SYM / EXACT45 leaves (u-bins)?
sym = Counter(); ex45 = Counter(); ax = Counter(); e0 = Counter()
for r, rt in zip(recs, roots):
    sym[rt[4]] += r['st']['SYM']; ex45[rt[4]] += r['st']['EXACT45']; ax[rt[4]] += r['st']['AXIS']; e0[rt[4]] += r['st']['EXACT0']
print('SYM by u0', dict(sym)); print('EXACT45 by u0', dict(ex45)); print('AXIS by u0', dict(ax)); print('EXACT0 by u0', dict(e0))
