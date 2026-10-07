#!/usr/bin/env python3
"""t3_chain_shape: read the duals recorded by t3_chain.py scan and classify their SHAPE.

An "H" is a tree: one wall-to-wall chain of T (the crossbar) plus transverse chains hanging off
it.  The question this script answers is whether the true dual is a tree at all -- the support
graph on the squares, with an edge per pair row in the support, has cyclomatic number
   mu  =  (#support pair rows counted as edges)  -  (#squares touched)  +  (#components) .
mu = 0 is a forest (H-shaped or smaller); mu >= 1 means the certificate contains a CYCLE of
links, which no H can express.
"""
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                     # noqa: E402


def support_graph(supp, wtol=1e-6):
    """-> (edges, verts, ncomp, mu, nwall_lo, nwall_hi, axis mix)"""
    edges, verts = [], set()
    walls = defaultdict(float)
    for w, tg, ax in supp:
        if w <= wtol:
            continue
        if tg[0] in ('lo', 'hi'):
            walls[(tg[0], tg[1])] += w
            verts.add(tg[2])
        else:
            edges.append((tg[1], tg[2], ax, w))
            verts.add(tg[1]); verts.add(tg[2])
    par = {v: v for v in verts}

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for (i, j, ax, w) in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            par[ri] = rj
    ncomp = len({find(v) for v in verts}) if verts else 0
    mu = len(edges) - len(verts) + ncomp
    nx = sum(1 for e in edges if e[2] == 'x')
    ny = len(edges) - nx
    return dict(nedge=len(edges), nvert=len(verts), ncomp=ncomp, mu=mu,
                nx=nx, ny=ny, walls=dict(walls), verts=sorted(verts))


def main():
    recs = []
    for f in sys.argv[1:]:
        recs += [json.loads(line) for line in open(f)]
    fams = defaultdict(list)
    for r in recs:
        r['g'] = support_graph(r['supp'])
        fams[r['fam']].append(r)

    print('# SHAPE of the true dual at each sampled optimum.  mu = cyclomatic number of the')
    print('# support graph (edge = pair row, vertex = square).  mu = 0: a forest, so an H can')
    print('# in principle express it.  mu >= 1: the certificate contains a CYCLE of links.')
    print(f'# {"family":13s} {"N":>3s} {"mu=0":>5s} {"mu=1":>5s} {"mu>=2":>5s} '
          f'{"med edges":>9s} {"med verts":>9s} {"H<=0":>7s} {"mu=0 & H>0":>10s}')
    order = sorted(fams, key=lambda f: (f[:4], f))
    tot = defaultdict(int)
    for fam in order:
        rs = fams[fam]
        m = [r['g']['mu'] for r in rs]
        ne = sorted(r['g']['nedge'] for r in rs)
        nv = sorted(r['g']['nvert'] for r in rs)
        ok = sum(1 for r in rs if r['H'] <= 1e-9)
        f0 = sum(1 for r in rs if r['g']['mu'] == 0 and r['H'] > 1e-9)
        for r in rs:
            tot['mu%d' % min(r['g']['mu'], 2)] += 1
            tot['n'] += 1
            if r['H'] > 1e-9:
                tot['Hfail'] += 1
                tot['Hfail_mu%d' % min(r['g']['mu'], 2)] += 1
        print(f'  {fam:13s} {len(rs):3d} {m.count(0):5d} {m.count(1):5d} '
              f'{sum(1 for x in m if x >= 2):5d} {ne[len(ne)//2]:9d} {nv[len(nv)//2]:9d} '
              f'{ok:3d}/{len(rs):<3d} {f0:10d}')
    print(f"\n# totals: {tot['n']} samples; mu=0 {tot['mu0']}, mu=1 {tot['mu1']}, "
          f"mu>=2 {tot['mu2']}")
    print(f"# H fails {tot['Hfail']}: of those mu=0 {tot.get('Hfail_mu0',0)}, "
          f"mu=1 {tot.get('Hfail_mu1',0)}, mu>=2 {tot.get('Hfail_mu2',0)}")

    print('\n# correlation: max tilt of the configuration vs mu')
    bins = defaultdict(lambda: [0, 0, 0])
    for r in recs:
        mt = max(abs(s6skel.norm_tilt(t)) for t in r['theta'])
        b = int(math.degrees(mt) // 5) * 5
        bins[b][min(r['g']['mu'], 2)] += 1
    print(f'# {"maxtilt":>8s} {"mu=0":>6s} {"mu=1":>6s} {"mu>=2":>6s}')
    for b in sorted(bins):
        print(f'  {b:3d}-{b+5:<4d} ' + ''.join(f'{x:6d}' for x in bins[b]))

    print('\n# the cycle-carrying duals: 12 examples, weights and tags')
    ex = [r for r in recs if r['g']['mu'] >= 1]
    ex.sort(key=lambda r: -(r['H'] - r['delta']))
    for r in ex[:12]:
        tl = ' '.join(f'{math.degrees(s6skel.norm_tilt(t)):+6.1f}' for t in r['theta'])
        print(f"  {r['fam']:12s} d*={r['delta']:+.4e} H={r['H']:+.4e} mu={r['g']['mu']} "
              f"E={r['g']['nedge']} V={r['g']['nvert']} tilts[{tl}]")
        print('      ' + '  '.join(f"{w:.3f}{tg[0]}{tuple(tg[1:3])}{ax}"
                                   for w, tg, ax in r['supp'][:9]))


if __name__ == '__main__':
    main()
