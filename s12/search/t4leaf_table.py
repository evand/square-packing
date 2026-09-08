#!/usr/bin/env python3
"""Read `search/t4leaf.py` output: the stage trajectory, and the anatomy of one measure.

    python3 search/t4leaf_table.py --table runs/tl_A01010101.json
    python3 search/t4leaf_table.py --anatomy runs/tl_B40K_measure_pure_nochord.txt --cliques

`--anatomy` prints the mass, the region split (corner boxes / wall slots / interior), the four
wall-strip masses (the chord quantity), the mass sitting on a region boundary, the support grouped
into geometric clusters, and -- with `--cliques` -- the largest anchor-clique mass the separators
of `search/anchorsep.py` can find on it (`> 1` means the measure is NOT clique-feasible).
"""
import argparse, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from level2_regions import classify                                  # noqa: E402

RNAME = ([f'C{i}' for i in range(4)] + [f'W{j}' for j in range(8)] + ['I'])


def read_measure(path):
    t = 4.0
    out = []
    for line in open(path):
        if line.startswith('#'):
            for tok in line.split():
                if tok.startswith('t='):
                    try:
                        t = float(tok[2:])
                    except ValueError:
                        pass
            continue
        q = line.split()
        if q and q[0] == 'pose':
            out.append((float(q[1]), float(q[2]), math.radians(float(q[3])), float(q[4])))
    return t, out


def table(paths):
    rows = []
    for p in paths:
        d = json.load(open(p))
        tag = d['tag']
        for r in d.get('hist', []):
            if r.get('it') != 'final':
                continue
            rows.append((tag, r))
    hdr = (f"{'tag':<11} {'st':>3} {'pattern':<9} {'LP':>11} {'M':>9} {'kmax':>8} "
           f"{'pure':>11} {'nochord':>11} {'p+nc':>11} {'int':>8} {'strips':>22} "
           f"{'cols':>6} {'poses':>6} {'rows':>6} {'cq':>5} {'bnd':>7} conv")
    print(hdr)
    for tag, r in rows:
        st = r.get('strips') or [0] * 4
        print(f"{tag:<11} {r.get('stage', -1):>3} {r.get('pattern', ''):<9} {r['LP']:>11.6f} "
              f"{r['M']:>9.6f} {r['kmax']:>8.5f} "
              f"{r.get('LP_pure', float('nan')):>11.6f} {r.get('LP_nochord', float('nan')):>11.6f} "
              f"{r.get('LP_pure_nochord', float('nan')):>11.6f} {r['interior']:>8.5f} "
              + '[' + ' '.join(f'{v:.3f}' for v in st) + '] '
              + f"{r['cols']:>6} {r.get('poses', 0):>6} {r['rows']:>6} {r['cq']:>5} "
                f"{r.get('bnd', 0):>7.4f} {int(bool(r.get('converged')))}")


def anatomy(path, r0, do_cliques, top, cluster):
    t, meas = read_measure(path)
    total = sum(m for *_, m in meas)
    print(f'{path}\n  t = {t}  r = {r0}  poses in support = {len(meas)}  mass = {total:.9f}')
    reg = {}
    for (cx, cy, th, m) in meas:
        k, j = classify(cx, cy, t, r0)
        lab = j if k == 'C' else (4 + j if k == 'W' else 12)
        reg.setdefault(lab, []).append((m, cx, cy, math.degrees(th)))
    tot = {k: sum(v[0] for v in vs) for k, vs in reg.items()}
    mc = sum(tot.get(i, 0) for i in range(4))
    mw = sum(tot.get(4 + j, 0) for j in range(8))
    mi = tot.get(12, 0.0)
    print(f'  corners  {mc:9.6f}   per box  ' + ' '.join(f'{tot.get(i,0):.5f}' for i in range(4)))
    print(f'  slots    {mw:9.6f}   per slot ' + ' '.join(f'{tot.get(4+j,0):.5f}' for j in range(8)))
    print(f'  interior {mi:9.6f}')
    P = np.array([[cx, cy, th, m] for (cx, cy, th, m) in meas])
    e = 1e-9
    st = [float(P[P[:, 1] <= r0 + e, 3].sum()), float(P[P[:, 1] >= t - r0 - e, 3].sum()),
          float(P[P[:, 0] <= r0 + e, 3].sum()), float(P[P[:, 0] >= t - r0 - e, 3].sum())]
    print(f'  wall strips (centres within {r0} of a wall; chord bound 3): '
          + ' '.join(f'{v:.6f}' for v in st))
    onb = np.zeros(len(P), bool)
    for b in (r0, t / 2.0, t - r0):
        onb |= (np.abs(P[:, 0] - b) <= 1e-7) | (np.abs(P[:, 1] - b) <= 1e-7)
    print(f'  mass on a region boundary (1e-7): {float(P[onb, 3].sum()):.6f} '
          f'on {int(onb.sum())} poses')
    ang = np.degrees(P[:, 2]) % 90.0
    axp = np.minimum(ang, 90 - ang) <= 1e-6
    print(f'  axis-parallel mass: {float(P[axp, 3].sum()):.6f}   tilted: {float(P[~axp, 3].sum()):.6f}')
    # ---- clusters: single-link in (cx, cy, theta) with the given radius
    idx = np.argsort(-P[:, 3])
    used = np.zeros(len(P), bool)
    cl = []
    for i in idx:
        if used[i]:
            continue
        d = (np.abs(P[:, 0] - P[i, 0]) < cluster) & (np.abs(P[:, 1] - P[i, 1]) < cluster) & \
            (np.abs(((np.degrees(P[:, 2]) - np.degrees(P[i, 2])) + 45) % 90 - 45) < 6.0) & ~used
        used |= d
        w = P[d, 3]
        cl.append((float(w.sum()), float((P[d, 0] * w).sum() / w.sum()),
                   float((P[d, 1] * w).sum() / w.sum()),
                   float(np.degrees(P[i, 2])), int(d.sum())))
    cl.sort(reverse=True)
    print(f'  --- {len(cl)} clusters (radius {cluster}, 6 deg), heaviest {min(top,len(cl))} ---')
    for (m, cx, cy, th, n) in cl[:top]:
        k, j = classify(cx, cy, t, r0)
        lab = j if k == 'C' else (4 + j if k == 'W' else 12)
        print(f'      mu={m:9.6f}  ({cx:8.5f}, {cy:8.5f})  {th:7.2f} deg  {RNAME[lab]:>3}  '
              f'({n} poses)')
    print(f'  cluster mass total {sum(c[0] for c in cl):.6f}')
    if do_cliques:
        import anchorclique as AC, anchorsep as ASEP
        Q = np.column_stack([P[:, 0], P[:, 1], P[:, 2], np.full(len(P), 0.5)])
        w = P[:, 3]
        D = 1000000
        found = ASEP.separate(Q, w, t, D, pitch=0.01, top=400, band=1.0, sym=False)
        found += ASEP.separate_interior(Q, w, t, D, pitch=0.02, top=200, ndir=24, neps=8,
                                        nrho=8, sym=False)
        found.sort(key=lambda z: -z[0])
        print(f'  --- anchor cliques: {len(found)} candidates, best mu(K) = '
              f'{found[0][0]:.6f}' if found else '  --- anchor cliques: none found')
        for (mk, mp, par, c) in found[:8]:
            print(f'      mu(K) = {mk:.6f}   mu(P_p) = {mp:.6f}   params {par}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', nargs='*', default=[])
    ap.add_argument('--anatomy', default=None)
    ap.add_argument('--r', type=float, default=1.0)
    ap.add_argument('--cliques', action='store_true')
    ap.add_argument('--top', type=int, default=30)
    ap.add_argument('--cluster', type=float, default=0.06)
    a = ap.parse_args()
    if a.table:
        table(a.table)
    if a.anatomy:
        anatomy(a.anatomy, a.r, a.cliques, a.top, a.cluster)


if __name__ == '__main__':
    main()
