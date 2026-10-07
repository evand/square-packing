#!/usr/bin/env python3
"""Force network of a jammed packing: the support of the jam-LP dual (contact forces balancing the pressure on the side).

LP: min ds s.t. linearised contacts (exact: gap 0), |z| <= 1.  At a first-order jam ds* = 0 and the duals lambda >= 0 satisfy
sum_k lambda_k grad g_k = e_s (equilibrium).  Force-bearing contacts = rows with lambda > thr; labelled by features
(contacts.pair_contact) so graphs can be compared across configurations (labels aligned to a reference).

  stress.py ref.txt other.txt ...   -> force-network sizes, diff vs ref, count of distinct force networks
"""
import math, sys, collections
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
from rigid import pair_rows, _sgns
from contacts import load, align, pair_contact, corners, describe

H = 0.5


def force_network(s, sq, delta=1e-7, thr=1e-9):
    n = len(sq)
    rows, labels = [], []
    for i, (x, y, t) in enumerate(sq):
        ct, st = math.cos(t), math.sin(t)
        w = H * (abs(ct) + abs(st))
        for s1 in _sgns(ct, 1e-6):
            for s2 in _sgns(st, 1e-6):
                dw = H * (-s1 * st + s2 * ct)
                for gap, cx, cy, cs, wall in ((x - w, 1, 0, 0, 'L'), (s - x - w, -1, 0, 1, 'R'), (y - w, 0, 1, 0, 'B'), (s - y - w, 0, -1, 1, 'T')):
                    if gap < delta:
                        c = {3 * i + 2: -dw}
                        if cx: c[3 * i] = cx
                        if cy: c[3 * i + 1] = cy
                        if cs: c[3 * n] = 1.0
                        rows.append(c); labels.append(('W', i, wall))
    for i in range(n):
        for j in range(i + 1, n):
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 > (math.sqrt(2) + delta) ** 2:
                continue
            gap, rr = pair_rows(sq[i], sq[j], 1e-6)
            if gap < delta:
                for r in rr:
                    rows.append({3 * i: r[0], 3 * i + 1: r[1], 3 * i + 2: r[2], 3 * j: r[3], 3 * j + 1: r[4], 3 * j + 2: r[5]})
                    labels.append(('P', i, j))
    A = lil_matrix((len(rows), 3 * n + 1))
    for k, c in enumerate(rows):
        for v, a in c.items():
            A[k, v] = -a
    cvec = np.zeros(3 * n + 1); cvec[-1] = 1.0
    res = linprog(cvec, A_ub=A.tocsr(), b_ub=np.zeros(len(rows)), bounds=[(-1, 1)] * (3 * n + 1), method='highs')
    lam = -res.ineqlin.marginals
    force = collections.defaultdict(float)
    for k, l in enumerate(lam):
        if l > thr:
            force[labels[k]] += l
    return res.fun, force


def features(sq, force, tol=1e-6):
    """Feature-labelled force network."""
    out = set()
    for lab, f in force.items():
        if lab[0] == 'W':
            i, wall = lab[1], lab[2]
            s_ = None
            out.add(('W', i, wall))
        else:
            fc = pair_contact(lab[1], lab[2], sq[lab[1]], sq[lab[2]], tol)
            out |= set(fc) if fc else {('P?', lab[1], lab[2])}
    return frozenset(out)


if __name__ == '__main__':
    s0, R = load(sys.argv[1])
    d0, F0 = force_network(s0, R)
    G0 = features(R, F0)
    tilted = lambda sq, i: min(math.degrees(sq[i][2]) % 90, 90 - math.degrees(sq[i][2]) % 90) > 0.5
    print(f'ref s={s0:.10f} ds*={d0:.1e}: force-bearing contacts {len(F0)} (involving {len({i for k in F0 for i in k[1:] if isinstance(i, int)})} squares), '
          f'{len(G0)} feature incidences')
    nets = collections.Counter()
    for p in sys.argv[2:]:
        s, sq = load(p)
        sq = align(sq, R)
        d, F = force_network(s, sq)
        G = features(sq, F)
        nets[G] += 1
        if len(sys.argv) <= 5:
            print(f'\n{p}: s={s:.10f} ds*={d:.1e}: force-bearing {len(F)}, feature incidences {len(G)}; lost {len(G0 - G)}, gained {len(G - G0)}')
            for e in sorted(G0 - G, key=str):
                print('   - ' + (describe(e) if e[0] in ('CS', 'SS') else str(e)))
            for e in sorted(G - G0, key=str):
                print('   + ' + (describe(e) if e[0] in ('CS', 'SS') else str(e)))
    print(f'\n{len(sys.argv) - 2} configurations, {len(nets)} distinct force networks')
