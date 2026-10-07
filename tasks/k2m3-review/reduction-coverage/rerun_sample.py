#!/usr/bin/env python3
"""Re-run a sample of V2 roots with leaf dumps (public code imported read-only, unmodified; clip_bin wrapped to log).
Checks: per-root stats reproduce the record; leaves + clipped-away slabs tile the root exactly (volume + disjointness);
every leaf kind is a certifying kind."""
import sys, gzip, json, random, time
from fractions import Fraction as F
from itertools import combinations
sys.path.insert(0, '/home/evand/math/square-packing/public/s12/search')
import zeromargin as zm, zm_mixed as ZM, mixed_cover as MC, qx2_zm as Q
P = '/home/evand/math/square-packing/public/s12/search/'
recs = [json.loads(l) for l in gzip.open(P + 'qx2_data/cert/runV2_cdade4b6_full.jsonl.gz', 'rt')][1:]
cov = ZM.Cover(MC.load(P + 'qx2_data/L4_k02_box7.txt'))
chk = Q.QXChecker(cov, exact_umax=F(1, 2), exact_from=3, use_exact=True, max_depth=18, use_chain=False, cert_mode=False)
chk.dump = True
LOG = []
_orig = zm.clip_bin
def logged(box, m, steps=28):
    c = _orig(box, m, steps); LOG.append((box, c)); return c
zm.clip_bin = logged
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
cheap = [r for r in recs if float(sys.argv[4] if len(sys.argv) > 4 else 0) <= r['st']['cpu'] < float(sys.argv[2] if len(sys.argv) > 2 else 5)]
pick = random.sample(cheap, int(sys.argv[3]) if len(sys.argv) > 3 else 40)
# make sure special kinds are represented
for k in ('SYM', 'EXACT45', 'AXIS', 'EXACT0', 'CAP', 'LEB', 'EXACT'):
    c = [r for r in cheap if r['st'][k] > 0 and r not in pick]
    if c: pick.append(random.choice(c))
GOOD = {'LEB', 'CAP', 'EXACT', 'EXACT0', 'EXACT45', 'AXIS', 'SYM', 'PIECE', 'ADM', 'P1', 'MIX', 'SPLIT', 'EMPTY'}
vol = lambda b: (b[1] - b[0]) * (b[3] - b[2]) * (b[5] - b[4])
def overlap(a, b):
    return all(max(a[2*i], b[2*i]) < min(a[2*i+1], b[2*i+1]) for i in range(3))
nbad = 0; t0 = time.time()
for r in pick:
    root = tuple(F(v) for v in r['root'])
    LOG.clear()
    st, unc, leaves = chk.run_box(root)
    same = all(st[k] == r['st'][k] for k in st if k != 'cpu')
    pieces = [b for b, k, w in leaves]
    kinds_ok = all(k in GOOD for b, k, w in leaves) and not unc
    # clipped-away slabs: box x (cu1, u1] for every clip that shrank the bin
    slabs = [(b[0], b[1], b[2], b[3], c, b[5]) for b, c in LOG if c < b[5]]
    # each slab: check no admissible pose (w(u) > K for u in (c, u1]) -- w increasing on [0,45deg], bin <= 45deg
    slab_ok = True
    for b, c in LOG:
        if c < b[5]:
            K = 2 * min(b[1], 7 - b[0], b[3], 7 - b[2])
            cc, ss = zm.trig(c)
            slab_ok &= (b[5] ** 2 + 2 * b[5] - 1 <= 0) and (cc + ss >= K) and (c >= b[4])
    allp = pieces + slabs
    inside = all(root[0] <= p[0] and p[1] <= root[1] and root[2] <= p[2] and p[3] <= root[3] and root[4] <= p[4] and p[5] <= root[5] for p in allp)
    vsum = sum(vol(p) for p in allp)
    disj = not any(overlap(a, b) for a, b in combinations([p for p in allp if vol(p) > 0], 2)) if len(allp) < 3000 else None
    # degenerate leaves (u0 == u1) are measure-zero: they are covered as limits; report them
    degen = [(b, k) for b, k, w in leaves if b[4] == b[5]]
    ok = same and kinds_ok and slab_ok and inside and vsum == vol(root) and disj is not False
    nbad += not ok
    print(r['root'], 'leaves', len(leaves), 'slabs', len(slabs), 'degen', len(degen), [k for b, k in degen][:3],
          'stats_same', same, 'kinds_ok', kinds_ok, 'slab_ok', slab_ok, 'tile(vol)', vsum == vol(root), 'disjoint', disj,
          f'{time.time()-t0:.0f}s', flush=True)
print('checked', len(pick), 'roots; bad', nbad)
