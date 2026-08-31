#!/usr/bin/env python3
"""Rung-2 closing loop (task H, step 2): iteratively resolve the cover LP over a FIXED column
set (no pricing -- the column set of the starting certificate, e.g. runs/closed4_best.txt) while
adding, each round, the actual worst poses found by a fine stress scan + local polish as hard
cutting-plane rows.  This targets the specific gap CLOSED4.md left open (LP converged on its own
sampled row set at 12.4174, but a finer stress scan finds 0.99267 at (3.40,1.42,76.4deg) -- the
row lattice just didn't sample near there): each round adds exactly the rows the previous round's
LP was blind to, so it is a genuine (float) separation loop, not another generic refinement.

The exact zero-margin checker (search/zeromargin.py) is the final word: this script produces a
candidate certificate; run zeromargin.py cert on it to see whether the residual gap is now zero
(and where any remaining uncertified boxes are) -- see search/ZEROMARGIN.md.

Usage:  python3 search/rung2_close.py CERT TAG [--rounds 60] [--nproc 6] [--time 10800]
Output: runs/closed4_TAG.log, runs/closed4_TAG_best.txt/.json (same format as closed4.py run)
"""
import sys, os, math, time, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import multiprocessing as mp
import closed4 as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cert'); ap.add_argument('tag')
    ap.add_argument('--rounds', type=int, default=60)
    ap.add_argument('--nproc', type=int, default=6)
    ap.add_argument('--time', type=float, default=10800)
    ap.add_argument('--stress-pitch', type=float, default=0.003)
    ap.add_argument('--nrand', type=int, default=300)
    ap.add_argument('--nseed', type=int, default=24, help='how many of the worst stress poses to polish/add per round')
    ap.add_argument('--cap', type=float, default=13.0, help='abort (report) if LP value would need to exceed this')
    # NOTE: this loop only reweights the fixed column set from `cert` -- it does not price in new
    # points.  The 2026-08-30 run (search/FAMILY.md sec 2) plateaus at stress_min ~0.994-0.996
    # after 67 rounds without reweighting closing the gap; column generation (new points at the
    # residual violation loci) is the natural next step and was not implemented here.
    a = ap.parse_args()
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/closed4_{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"rung2_close.py {' '.join(sys.argv[1:])}")

    s, P0, w0c = C.load_points(a.cert)
    m = C.Model(s, D=1000, sym=True)
    ncol = sum(m.add_point(x, y, 'cert') for x, y in P0)
    wmap = {(int(round(x * m.D)), int(round(y * m.D))): wt for (x, y), wt in zip(P0, w0c)}
    w0 = np.array([wmap.get((int(round(px * m.D)), int(round(py * m.D))), 0.0) for px, py in m.P])
    log(f"[{a.tag}] s={s} columns from {a.cert}: {ncol} orbits / {len(m.P)} atoms, file total {w0.sum():.6f}")

    thetas = [math.radians(d) for d in C.angle_list(0.5)]
    pool = mp.get_context('fork').Pool(a.nproc)
    gmin, nviol, poses, vals = C.separate(m.P, w0, s, thetas, 0.01, thr=1.0 + 0.02, perang=800, block=0.03, pool=pool)
    m.add_rows(poses)
    log(f"[{a.tag}] initial rows (0.01 lattice, thr 1.02): {len(m.rows)}, lattice min {gmin:.6f}")

    t0 = time.time(); hist = []; best = None
    for r in range(a.rounds):
        out = m.solve()
        if out is None: log("LP failed"); break
        val, x, y, Ax = out; wt = x[m.own]
        vmin, worst, wp = C.stress(m.P, wt, s, pitch=a.stress_pitch, nproc=a.nproc, nrand=a.nrand, log=lambda *_a: None)
        rec = dict(round=r, t=round(time.time() - t0, 1), LP=val, rows=len(m.rows), stress_min=vmin,
                   worst_pose=list(wp))
        hist.append(rec)
        log(f"[{a.tag}] round {r} LP={val:.6f} rows={len(m.rows)} stress_min={vmin:.7f} "
            f"worst=({wp[0]:.5f},{wp[1]:.5f},{wp[2]:.4f}deg) t={time.time()-t0:.0f}s")
        json.dump(dict(s=s, args=vars(a), hist=hist), open(f"runs/closed4_{a.tag}.json", 'w'), indent=1)
        # checkpoint every round (2026-08-30 coordinator note): never lose a plateau's weights to
        # a kill/timeout -- the caller is responsible for then scaling this by 1/stress_min (or a
        # safety factor above it) if it wants a candidate valid cover; this file alone is NOT one
        # (stress_min < 1 means it is a known-incomplete cover as exported).
        C.export(m, x, f"runs/closed4_{a.tag}_last.txt", WD=10 ** 7, up=True)
        if vmin >= 1 - 1e-7:
            log(f"[{a.tag}] CONVERGED: stress min {vmin:.9f} >= 1"); best = (val, r); break
        if val > a.cap:
            log(f"[{a.tag}] LP value {val:.6f} exceeds cap {a.cap}: the fixed column set from {a.cert} "
                f"cannot close this gap under weight {a.cap} -- reporting, not forcing further.")
            best = (val, r); break
        # seeds: the worst-K stress poses (already sorted) + the single polished worst
        seeds = [[p[0], p[1], math.radians(p[2])] for _, p in worst[:a.nseed]]
        seeds.append([wp[0], wp[1], math.radians(wp[2])])
        seeds = np.array(seeds)
        fp, fv, cur, curv = C.polish(m.P, wt, s, seeds, rounds=10, nper=96,
                                      rng=np.random.default_rng(1000 + r), thmax=math.pi / 4,
                                      arad0=math.radians(1.0))
        n1 = m.add_rows(fp); n2 = m.add_rows(seeds); n3 = m.add_rows(cur)
        # also re-scan the lattice at the current weights to catch any other now-violated poses
        gmin, nviol, poses, vals2 = C.separate(m.P, wt, s, thetas, 0.01, thr=1.0 - 1e-7, perang=100, block=0.08, pool=pool)
        n4 = m.add_rows(poses)
        added = n1 + n2 + n3 + n4
        log(f"   +rows {n1}(polish-visited) {n2}(seeds) {n3}(polish-final) {n4}(lattice) = {added}; rows now {len(m.rows)}")
        if added == 0:
            log(f"[{a.tag}] no new rows found but stress_min {vmin:.7f} < 1: stalled"); best = (val, r); break
        if time.time() - t0 > a.time:
            log(f"[{a.tag}] time limit"); best = (val, r); break

    pool.close(); pool.join()
    out = m.solve(); val, x, y, Ax = out; wt = x[m.own]
    log(f"[{a.tag}] final LP over all {len(m.rows)} rows: {val:.6f}")
    tw, npts = C.export(m, x, f"runs/closed4_{a.tag}_best.txt", WD=10 ** 7, up=True)
    bd = C.weight_breakdown(m.P, wt, s, m.D)
    js = dict(s=s, LP=val, exported_total=tw, points=npts, breakdown=bd, rows=len(m.rows),
              hist=hist, points_list=[(float(px), float(py), float(wp_)) for (px, py), wp_ in zip(m.P, wt) if wp_ > 1e-9])
    json.dump(js, open(f"runs/closed4_{a.tag}_best.json", 'w'), indent=1)
    log(f"[{a.tag}] exported runs/closed4_{a.tag}_best.txt: {npts} points, total {tw:.6f}")
    vmin, worst, wp = C.stress(m.P, wt, s, pitch=a.stress_pitch, nproc=a.nproc, nrand=a.nrand, log=log)
    log(f"[{a.tag}] RESULT s={s} LP={val:.6f} exported total={tw:.6f} stress min={vmin:.7f} cost=total/min={tw/vmin:.6f}")


if __name__ == '__main__':
    main()
