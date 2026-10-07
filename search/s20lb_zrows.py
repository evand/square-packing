#!/usr/bin/env python3
"""s20lb_zrows.py COVER OUT LOG [LOG ...] -- LP rows from zmx2's uncertified boxes (search/S20_LB.md).

zmx2 leaves a box UNCERTIFIED where no point is certain at every pose of the box; near a pose where several
points sit on the boundary of Q at once (a vertex of the arrangement) the float scans and the LP's tolerant
rows (closed, +TOL) count all of them, while just off the vertex some leave.  This is the "near-tight cells as
exact rows" lever (notes/jlevy-s17-techniques.md sec 2 C) in float form: for every UNCERT box, sample poses in
a few shells around the box centre (|dx|, |dy| <= r, |dtheta| <= r), evaluate the float mass of COVER there,
and write the K lowest distinct poses per box (all below --thr) as `cx cy theta_rad mu` rows for
line_cover.py loop --rows-from.  Heuristic only: it chooses rows, it certifies nothing.
"""
import argparse, math, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import line_cover as LC
import closed4 as C


def boxes(path):
    out = []
    for line in open(path):
        if not line.startswith('UNCERT'): continue
        b = line.split('box ')[1].split()[0].split(',')
        v = [int(t.split('/')[0]) / int(t.split('/')[1]) for t in b]
        out.append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cover'); ap.add_argument('out'); ap.add_argument('logs', nargs='+')
    ap.add_argument('--n', type=int, default=4000, help='samples per box')
    ap.add_argument('--k', type=int, default=6, help='rows kept per box')
    ap.add_argument('--thr', type=float, default=1.002, help='keep only poses with float mu below this')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    cov = LC.cover_from_file(a.cover)[0]; s = cov.s; rng = np.random.default_rng(a.seed)
    B = []
    for p in a.logs: B += boxes(p)
    # dedupe boxes by centre (1e-4)
    seen = set(); Bc = []
    for x0, x1, y0, y1, u0, u1 in B:
        c = ((x0 + x1) / 2, (y0 + y1) / 2, 2 * math.atan((u0 + u1) / 2))
        key = (round(c[0], 3), round(c[1], 3), round(c[2], 3))
        if key in seen: continue
        seen.add(key); Bc.append(c)
    QQ, ID = [], []
    for b, (cx, cy, th) in enumerate(Bc):
        for r in (2e-6, 2e-5, 2e-4, 1e-3):
            n = a.n // 4
            QQ.append(np.c_[cx + rng.uniform(-r, r, n), cy + rng.uniform(-r, r, n), th + rng.uniform(-r, r, n)])
            ID.append(np.full(n, b))
    rows = []
    if QQ:
        Q = np.concatenate(QQ); ID = np.concatenate(ID)
        Q[:, 2] = np.clip(Q[:, 2], 0.0, math.pi / 2 - 1e-12)
        lo = (np.abs(np.cos(Q[:, 2])) + np.abs(np.sin(Q[:, 2]))) / 2; hi = s - lo   # closed4.admissible_box
        ok = (Q[:, 0] >= lo) & (Q[:, 0] <= hi) & (Q[:, 1] >= lo) & (Q[:, 1] <= hi)
        Q = Q[ok]; ID = ID[ok]
        v = cov.mass(Q) if len(Q) else np.zeros(0)
        for b in np.unique(ID):
            idx = np.nonzero(ID == b)[0]; o = idx[np.argsort(v[idx])]; kept = []
            for i in o:
                if v[i] >= a.thr or len(kept) >= a.k: break
                if all(np.abs(Q[i] - Q[j]).max() > 1e-6 for j in kept): kept.append(i)
            rows += [(*Q[i], v[i]) for i in kept]
    with open(a.out, 'w') as f:
        f.write(f"# s20lb_zrows.py {a.cover} {' '.join(a.logs)}: {len(Bc)} boxes, {len(rows)} rows\n")
        for r_ in rows: f.write("%.9f %.9f %.12f %.8f\n" % r_)
    mv = min((r_[3] for r_ in rows), default=float('nan'))
    print(f"{len(Bc)} distinct UNCERT boxes -> {len(rows)} rows (min float mu {mv:.6f}) -> {a.out}")


if __name__ == '__main__':
    main()
