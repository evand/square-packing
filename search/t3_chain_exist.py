#!/usr/bin/env python3
"""t3_chain_exist: does an H actually EXIST in the sampled optimum, and is it long enough?

For every configuration recorded by `t3_chain.py scan` this asks the geometric question the
H-lemma needs:  is there a wall-to-wall chain of T (at the configuration's own margin) whose low
end carries a transverse chain going down to the far wall and whose high end carries one going
up, with a total of L links -- and is L at least the threshold L*(T, t) of
`search/t3_chain_hform.py`?

Everything here is a MEASUREMENT on feasible points found by a multistart, never a proof.
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                    # noqa: E402
import bandcut_k as BK                                           # noqa: E402
from t3_chain_hform import Lstar                                 # noqa: E402


def longest_leg(g, dr, lvl, node, forbid, up, n):
    """longest path of the transverse DAG that ends at (up=False) or starts at (up=True) `node`,
    avoiding `forbid`.  -> (number of links, the path)."""
    best = [0, [node]]

    def walk(path, cur):
        if len(path) - 1 > best[0]:
            best[0], best[1] = len(path) - 1, list(path)
        for j in range(n):
            if j in path or j in forbid:
                continue
            ok = (g[cur, j] >= lvl and dr[cur, j] > 0) if up else \
                 (g[j, cur] >= lvl and dr[j, cur] > 0)
            if ok:
                walk(path + [j], j)
    walk([node], node)
    return best[0], best[1]


def hexist(z, n, T, lvl, tol=1e-7):
    """-> (best L over all chains, description).  L = a + c, the number of transverse links that
    can be hung on the two ends of some wall-to-wall chain of T, all squares distinct."""
    X, Y = np.array(z[1::3]), np.array(z[2::3])
    gx, dxr, gy, dyr = BK.normal_gaps(np.asarray(z, dtype=float), n)
    lv = lvl - tol
    best = (-1, None)
    for (g, dr, co, gt, drt, ax) in ((gx, dxr, X, gy, dyr, 'x'),
                                     (gy, dyr, Y, gx, dxr, 'y')):
        for ch in BK.all_chains(list(range(n)), g, dr, co, lv, int(round(T))):
            f = set(ch)
            a1, pa1 = longest_leg(gt, drt, lv, ch[0], f, False, n)
            c1, pc1 = longest_leg(gt, drt, lv, ch[-1], f | set(pa1), True, n)
            c2, pc2 = longest_leg(gt, drt, lv, ch[-1], f, True, n)
            a2, pa2 = longest_leg(gt, drt, lv, ch[0], f | set(pc2), False, n)
            for (tot, a, c, pa, pc) in ((a1 + c1, a1, c1, pa1, pc1),
                                        (a2 + c2, a2, c2, pa2, pc2)):
                if tot > best[0]:
                    best = (tot, dict(axis=ax, chain=list(ch), a=a, c=c,
                                      down=list(pa), up=list(pc)))
    return best


def main():
    recs = []
    for f in sys.argv[1:]:
        recs += [json.loads(line) for line in open(f)]
    fams = defaultdict(list)
    print('# L = transverse links hangable on the ends of some wall-to-wall chain of 3;')
    print('# Lneed = ceil(L*(3, maxtilt)) from the closed form.  "ok" = L >= Lneed, i.e. the')
    print('# H-lemma has a long enough H available at that optimum (with a COMMON tilt read as')
    print('# the configuration maximum -- a conservative reading).')
    for r in recs:
        z = np.array(r['z'])
        L, info = hexist(z, 6, 3.0, r['delta'])
        tmax = max(abs(s6skel.norm_tilt(q)) for q in r['theta'])
        need = Lstar(3, tmax)
        r['L'] = L
        r['Lneed'] = need
        r['ok'] = (L >= need - 1e-12)
        r['info'] = info
        fams[r['fam']].append(r)
    print(f'# {"family":13s} {"N":>3s} {"medL":>5s} {"minL":>5s} {"medLneed":>9s} '
          f'{"L>=Lneed":>9s} {"H<=0":>7s} {"L>=3":>6s}')
    for fam in sorted(fams):
        rs = fams[fam]
        Ls = sorted(r['L'] for r in rs)
        nd = sorted(r['Lneed'] for r in rs)
        print(f'  {fam:13s} {len(rs):3d} {Ls[len(Ls)//2]:5d} {Ls[0]:5d} '
              f'{nd[len(nd)//2]:9.3f} {sum(1 for r in rs if r["ok"]):4d}/{len(rs):<4d} '
              f'{sum(1 for r in rs if r["H"] <= 1e-9):3d}/{len(rs):<3d} '
              f'{sum(1 for r in rs if r["L"] >= 3):3d}/{len(rs):<3d}')
    bad = [r for r in recs if not r['ok']]
    print(f'\n# no long-enough H at {len(bad)} / {len(recs)} sampled optima.')
    bad.sort(key=lambda r: (r['L'] - r['Lneed']))
    for r in bad[:20]:
        tl = ' '.join(f'{math.degrees(s6skel.norm_tilt(t)):+6.1f}' for t in r['theta'])
        print(f"  {r['fam']:12s} d*={r['delta']:+.4e} H={r['H']:+.4e} L={r['L']} "
              f"need {r['Lneed']:.3f}  {r['info']}")
        print(f'      tilts[{tl}]')
    agree = sum(1 for r in recs if r['ok'] == (r['H'] <= 1e-9))
    print(f'\n# "a long enough H exists" agrees with "the H-restricted dual is <= 0" on '
          f'{agree}/{len(recs)} samples')


if __name__ == '__main__':
    main()
