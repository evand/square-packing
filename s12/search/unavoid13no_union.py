#!/usr/bin/env python3
"""unavoid13-no: h of the union of several families (with a bulk), by highspy; dumps the instance.
    python3 search/unavoid13no_union.py --out runs/unavoid13no_union/u1 --bulk support:80 FAM1 FAM2 ...
"""
import sys, os, time, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13no_lib as N
from unavoid13no_lp import build_family


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('families', nargs='+'); ap.add_argument('--out', required=True); ap.add_argument('--bulk', default='support:80')
    ap.add_argument('--threads', type=int, default=2); ap.add_argument('--k', type=int, default=13); ap.add_argument('--time', type=float, default=7200)
    ap.add_argument('--exclude-bulk', default='', help="drop squares of this family spec (e.g. an A1 bulk) from the union")
    a = ap.parse_args()
    m = F(4)
    parts = [build_family(m, a.bulk)]
    excl = set(build_family(m, a.exclude_bulk)) if a.exclude_bulk else set()
    for f in a.families:
        parts.append([p for p in L.read_family(f)[1] if p not in excl])
    Fam = L.d4_closure(m, N.family_union(m, parts))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    L.write_family(a.out + '_family.txt', m, Fam, header='union of ' + ' '.join(a.families))
    print(f"union: {len(Fam)} squares", flush=True)
    sq, V, B = N.maximal_cells(m, Fam, threads=a.threads)
    V, B, perms = N.d4_close_cells(m, Fam, V, B)
    np.save(a.out + '_B.npy', np.packbits(B, axis=1)); np.save(a.out + '_V.npy', V)
    lp, _ = N.lp_value(B, threads=a.threads); print(f"LP {lp:.4f}", flush=True)
    j0 = Fam.index(L.pose_key(F(0), m / 2, m / 2))
    cut = N.root_cut_cells(m, V, B, j0, perms=perms)
    r = N.solve_ip(B, a.k, threads=a.threads, time_limit=a.time, extra_rows=[cut], feasibility=False)
    print(f"highspy k={a.k}: {r['status']} feasible={r['feasible']} nodes {r['nodes']} bound {r['bound']} {r['time']:.0f}s", flush=True)
    if r['x'] is not None:
        P = V[np.nonzero(r['x'] > 0.5)[0]]; np.save(a.out + '_P.npy', P)
        print("  set: " + " ".join(f"({x:.3f},{y:.3f})" for x, y in P))
    json.dump(dict(n=len(Fam), K=int(B.shape[0]), lp=lp, status=r['status'], feasible=r['feasible'], time=r['time']), open(a.out + '.json', 'w'))


if __name__ == '__main__':
    main()
