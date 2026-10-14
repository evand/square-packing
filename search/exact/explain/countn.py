#!/usr/bin/env python3
"""Staged unknown/equation count for record N (regularized view): pick an independent set of contacts by Jacobian
rank at the packing (walls first, then upright pairs, then contacts with tilted squares) and classify them.

  python3 search/exact/explain/countn.py N
"""
import sys, numpy as np, sympy as sp, math, render, contacts as C, fullsys

def classify(n):
    Sf, sq = render.load(n)
    eqs, circ, unk, val = fullsys.build(n)
    rows = []
    for ct in C.contacts(Sf, sq):
        e = fullsys.build(n, only=[ct])[0]
        if not e: continue
        kind = 'wall' if ct[0] == 'wall' else 'pair'
        ids = {ct[1]} if kind == 'wall' else {ct[1], ct[3]}
        tilted = any(abs(math.remainder(sq[i][2], math.pi/2)) > 1e-9 for i in ids)
        rows.append((ct, e[0], kind, tilted))
    rows.sort(key=lambda r: (r[3], r[2] != 'wall'))
    J = lambda es: np.array([[float(sp.diff(e, v).subs(val)) for v in unk] for e in es])
    chosen = []
    for r in rows:
        if np.linalg.matrix_rank(J([x[1] for x in chosen] + [r[1]] + circ), tol=1e-9) > len(chosen) + len(circ):
            chosen.append(r)
    return chosen, len(unk), len(circ)

if __name__ == '__main__':
    n = int(sys.argv[1])
    chosen, nu, nc = classify(n)
    triv = [r for r in chosen if not r[3]]; tilt = [r for r in chosen if r[3]]
    tp = sum(r[2] == 'pair' for r in tilt); tw = len(tilt) - tp
    print(f'n={n}: {nu} unknowns; {len(chosen)} independent contacts')
    print(f'  trivial (upright squares only): {len(triv)} ({sum(r[2] == "wall" for r in triv)} wall, '
          f'{sum(r[2] == "pair" for r in triv)} pair) -> {nu - len(triv)} unknowns left')
    print(f'  with tilted squares: {len(tilt)} ({tp} pair: 2 eqs + 1 slide each; {tw} wall)')
    print(f'  with slides: {2*tp + tw + nc} equations in {nu - len(triv) + tp} unknowns ({nc} c^2+s^2=1)')
