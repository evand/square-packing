#!/usr/bin/env python3
"""One-off: refresh the parts of an exported analysis that depend on free_region / NEAR_MAX,
without redoing the expensive rigidity search (a full re-export is ~11 h).

Recomputes, per packing:
  * near_misses          -- the gap ceiling NEAR_MAX changed
  * regions[i]           -- the sampler now also finds slides along a zero-width channel
  * free / mobile / wedged / first_order_only / rigid  -- follow from which regions are free

Everything else (contacts, angle groups, symmetry, the verified motions) is carried over unchanged.
Usage: patch_regions.py [names...]   (default: every exported packing)
"""
import os, sys, json, glob, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis as A
from export import compact, build_index

HERE = os.path.dirname(os.path.abspath(__file__))
D = HERE + '/../data/'; W = HERE + '/../www/data/'

def near_misses(sq, s):
    """The near-miss list of analysis.analyze(), recomputed at the current NEAR_MAX."""
    st = np.column_stack([sq[:, 0], sq[:, 1], np.radians(sq[:, 2])])
    P = A._corners_all(st)
    near = []
    for w, (ax, val) in enumerate(((0, 0.), (0, s), (1, 0.), (1, s))):
        Dw = P[:, :, ax] - val if val == 0 else val - P[:, :, ax]
        G = Dw.min(1)
        for i in np.where((G >= A.CONTACT_TOL) & (G < A.NEAR_MAX))[0]:
            near.append({'i': int(i), 'j': None, 'wall': w, 'gap': float(G[i])})
    pi, pj = A.near_pairs(sq[:, :2], math.sqrt(2) + A.NEAR_MAX)
    if len(pi):
        pg = A._sat_batch(P, A._edge_normals_all(P), pi, pj)[0]
        for m in np.where((pg >= A.CONTACT_TOL) & (pg < A.NEAR_MAX))[0]:
            near.append({'i': int(pi[m]), 'j': int(pj[m]), 'gap': float(pg[m])})
    near.sort(key=lambda m: m['gap'])
    return near[:40], P

def patch(name, packings):
    outp = W + 'p/' + name.replace('.svg', '.json')
    if not os.path.exists(outp): return None
    d = json.load(open(outp)); a = d['analysis']
    sq = np.array([[float(x) for x in q] for q in packings[name]['squares']])
    s = float(packings[name]['s'])
    near, P = near_misses(sq, s)
    free_alone = np.array(a['free_alone'], bool)
    # the LP's first-order directions are not exported; the geometric candidates alone find the slides
    regions = {int(i): A.free_region(int(i), sq, P, s) for i in np.where(free_alone)[0]}
    truly_free = sorted(i for i, r in regions.items() if r.get('free', True))
    verified = np.array(a['mobile'], bool).copy()
    for i in truly_free: verified[i] = True
    # `first_order_only` is (firstorder & ~verified_old); since verified only grows, that is all the
    # firstorder information the new `wedged` can need.
    firstorder = set(a['first_order_only'])
    wedged = [i for i in range(len(sq)) if not verified[i] and (free_alone[i] or i in firstorder)]
    was_free = list(a['free'])
    a.update({'near_misses': near, 'regions': regions, 'free': truly_free, 'wedged': wedged,
              'mobile': verified.tolist(), 'rigid': not bool(verified.any()),
              'first_order_only': sorted(i for i in firstorder if not verified[i])})
    d['analysis'] = compact(a)
    json.dump(d, open(outp, 'w'), separators=(',', ':'))
    return was_free, truly_free, wedged

def main(argv):
    packings = json.load(open(D + 'packings.json'))
    names = argv or sorted(k for k in packings if os.path.exists(W + 'p/' + k.replace('.svg', '.json')))
    t0 = time.time(); changed = 0
    for k, name in enumerate(names):
        res = patch(name, packings)
        if res is None: continue
        was_free, now_free, wedged = res
        newly = sorted(set(now_free) - set(was_free))
        if newly:
            changed += 1
            print('%-24s free += %s  (now %d free, %d wedged)' % (name, newly[:8], len(now_free), len(wedged)), flush=True)
        if k % 50 == 0: print('  ... %d/%d  %.0fs' % (k, len(names), time.time() - t0), flush=True)
    print('%d packings gained free squares; %.0fs' % (changed, time.time() - t0))
    build_index(packings)

if __name__ == '__main__':
    main([a for a in sys.argv[1:] if not a.startswith('--')])
