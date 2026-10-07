#!/usr/bin/env python3
"""unavoid13-no step 1: maximal cells of a family, D4-closed, saved; LP relaxation value.

    python3 search/unavoid13no_lp.py --init support,t4L --out runs/unavoid13no_lp/full
"""
import sys, os, time, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13_loop as LOOP, unavoid13no_lib as N


def build_family(m, init):
    parts = []
    for item in init.split(','):
        if item == 'support': parts.append(N.support_poses(m))
        elif item.startswith('support:'):
            sup = L.read_support('search/cover4_exact_support.txt', m); sup.sort(key=lambda t: -t[1])
            parts.append(L.d4_closure(m, [p for p, _ in sup[:int(item.split(':')[1])]]))
        elif item == 't4L': parts.append(L.read_family('runs/unavoid13_t4L/F_round01.txt')[1])
        elif item.startswith('file:'): parts.append(L.read_family(item[5:])[1])
        else: parts.append(LOOP.initial_family(m, item))
    return L.d4_closure(m, N.family_union(m, parts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', default='4'); ap.add_argument('--init', default='tiles,grid2,support,t4L')
    ap.add_argument('--out', required=True); ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--no-lp', action='store_true'); ap.add_argument('--no-d4', action='store_true')
    a = ap.parse_args()
    m = F(a.m); os.makedirs(os.path.dirname(a.out), exist_ok=True)
    Fam = build_family(m, a.init)
    print(f"family {a.init}: {len(Fam)} squares", flush=True)
    sq, V, B = N.maximal_cells(m, Fam, threads=a.threads)
    if not a.no_d4:
        V, B, perms = N.d4_close_cells(m, Fam, V, B)
        np.save(a.out + '_perms.npy', perms)
    L.write_family(a.out + '_family.txt', m, Fam)
    np.save(a.out + '_B.npy', np.packbits(B, axis=1)); np.save(a.out + '_V.npy', V)
    print(f"saved {a.out}_*: {B.shape[0]} cells x {B.shape[1]} rows", flush=True)
    if not a.no_lp:
        t0 = time.time(); lp, y = N.lp_value(B, threads=a.threads)
        print(f"LP over all cells: {lp:.6f} ({time.time()-t0:.0f}s); dual support {int((y > 1e-9).sum())} rows", flush=True)
        np.save(a.out + '_dual.npy', y)
        json.dump(dict(init=a.init, n=len(Fam), K=int(B.shape[0]), lp=lp), open(a.out + '_lp.json', 'w'))


if __name__ == '__main__':
    main()
