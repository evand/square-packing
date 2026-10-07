#!/usr/bin/env python3
"""Run the n = 11 cover LP (`search/tighten.py reopt`, column generation + cuts) at one container
and DUMP ITS DUAL: the fractional packing measure on placements that blocks the method.

The cover LP of `tighten.py` is

    min  sum_k |orbit_k| x_k   s.t.   sum_k count(r, k) x_k >= 1 + margin   for every row r,

rows r being closed squares (cx, cy, theta) of half-side h_r (the verifier's angle-net shrink
sigma_k/2 <= 1/2) and columns the D4 orbits of points.  Its dual is

    max  (1+margin) sum_r y_r   s.t.   sum_r count(r, k) y_r <= |orbit_k|   for every orbit k,

i.e. a nonnegative measure y on placements whose coverage, averaged over each D4 orbit of points,
is at most 1 -- a *fractional packing* of mass  sum_r y_r = LP value / (1+margin).  When the LP
value is 11, y is the mass-11 fractional packing that makes a weight-< 11 unavoidable set
impossible on this row/column set.  (It is not a packing measure in the sense of
`search/packing_dual.py`: its squares are shrunk by the angle net and its coverage is constrained
only at the model's points.  `packing_dual.py` / `dual_exact.py` produce the certified object;
this one is what the cover LP itself sees.)

Writes, at the tag TAG:
    runs/n11_cd_TAG.log              the tighten log
    runs/n11_cd_TAG_dual.txt         the dual measure, `pose cx cy theta_deg mu` + `# h` (NOT D4-
                                     symmetrised: every row is its own placement; `sym 1`)
    runs/n11_cd_TAG_cover.txt        the primal cover (a certificate file at W = 10^7)
    runs/n11_cd_TAG.json             LP value, mass, max symmetrised dual coverage, timings

    python3 search/n11_coverdual.py runs/n11_seed_union2.txt CD385 --Dp 3948 --colgen 12 --threads 2
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tighten as TG                                                  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(REPO, 'runs')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('seed'); ap.add_argument('tag')
    ap.add_argument('--Dp', type=int, required=True)
    ap.add_argument('--colgen', type=int, default=12)
    ap.add_argument('--cg-want', type=int, default=150)
    ap.add_argument('--N', type=int, default=6000)
    ap.add_argument('--topk', type=int, default=6)
    ap.add_argument('--margin', type=float, default=2e-6)
    ap.add_argument('--n', type=int, default=11)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--max-iters', type=int, default=200)
    a = ap.parse_args()

    TG.NPROC = a.threads                      # the verifier and the float separator honour this
    lf = open(os.path.join(RUNS, f'n11_cd_{a.tag}.log'), 'w')

    def log(msg):
        print(msg, flush=True); lf.write(str(msg) + '\n'); lf.flush()

    t0 = time.time()
    m, w0, s = TG.build_model(a.seed, a.Dp)
    log(f"[{a.tag}] {a.seed}: {len(m.orbits)} orbits / {len(m.P)} atoms, container {s} = {float(s):.7f}")
    lr = TG.lattice_rows(m, w0[m.own], thr=1.05, nproc=a.threads)
    m.add_rows(lr); log(f"[{a.tag}] warm start: {len(m.rows)} lattice rows ({time.time()-t0:.0f}s)")
    x, val, info = TG.tighten(m, a.tag, margin=a.margin, N=a.N, topk=a.topk, log=log,
                              colgen=a.colgen, cg_want=a.cg_want, n=a.n, max_iters=a.max_iters)
    if x is None:
        log(f"[{a.tag}] FAILED"); return
    out = m.solve(a.margin)
    val, x, y = out
    mass = float(y.sum())
    cand = TG.price(m, y, pitch=0.01, want=1)
    dualcov = cand[0][0] if cand else 1.0
    nzy = int((y > 1e-12).sum())
    log(f"[{a.tag}] LP value = {val:.7f}   dual mass = {mass:.7f}   dual support = {nzy} rows / {len(m.rows)}"
        f"   max symmetrised dual coverage on the 0.01 grid = {dualcov:.6f}")

    dpath = os.path.join(RUNS, f'n11_cd_{a.tag}_dual.txt')
    with open(dpath, 'w') as f:
        f.write(f"# t = {s}\n# n = {a.n} cover LP dual (fractional packing on placements); "
                f"LP={val:.9f} mass={mass:.9f} margin={a.margin}\n")
        f.write("# NOT D4-symmetrised: each row is one placement.  half-side h <= 1/2 is the "
                "verifier's angle-net shrink\nsym 1\n")
        f.write("# pose cx cy theta_deg mu    h\n")
        for r in np.nonzero(y > 1e-12)[0]:
            cx, cy, th, h = m.rows[r]
            f.write(f"pose {cx:.12f} {cy:.12f} {math.degrees(th):.11f} {y[r]:.12f}    {h:.12f}\n")
    cpath = os.path.join(RUNS, f'n11_cd_{a.tag}_cover.txt')
    tw, npts = TG.export(m, x, cpath, WD=10 ** 7, up=True)
    log(f"[{a.tag}] wrote {dpath} ({nzy} placements) and {cpath} ({npts} points, total {tw:.7f})")
    json.dump(dict(tag=a.tag, t=str(s), tf=float(s), Dp=a.Dp, n=a.n, lp=float(val), mass=mass,
                   dual_support=nzy, rows=len(m.rows), orbits=len(m.orbits), dualcov=float(dualcov),
                   cover_points=npts, cover_total=tw, seconds=time.time() - t0, info=info),
              open(os.path.join(RUNS, f'n11_cd_{a.tag}.json'), 'w'), indent=1)
    log(f"[{a.tag}] done in {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
