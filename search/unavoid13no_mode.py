#!/usr/bin/env python3
"""unavoid13-no: targeted killer harvest for the single 13-set mode the loops converged on
(notes/unavoid13-no.md s7): three columns of three points on y = 1, 2, 3 at x ~ (0.9, 1.4, 3.1) and a
column of four at x ~ 2.25 with y ~ (0.95, 1.5, 2.45, 3.0), plus D4 images.

    python3 search/unavoid13no_mode.py --resume FAMILY.txt --tag M1 [--tries 2000] [--add 6]

Each try: random mode parameters -> max-margin LP polish against ALL of F (unavoid13_lib.polish_positions,
no cells, no IP) -> if the polished set hits F (t* >= -1e-9): float violation search, rationalise the
worst poses, add them and their D4 images to F.  No IP is solved; the output family (F_final.txt) is a
row source for the exact loop.  Heuristic: it can only ADD squares (sound for the lower bound); a
failure to find an F-hitting mode set proves nothing.
"""
import sys, os, math, time, argparse
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L, unavoid13no_lib as N


def mode_set(rng):
    x3 = [0.90, 1.40, 3.10] + rng.normal(0, 0.06, 3)
    x4 = 2.25 + rng.normal(0, 0.08)
    y3 = [1.0, 2.0, 3.0]
    y4 = np.array([0.95, 1.52, 2.45, 3.02]) + rng.normal(0, 0.06, 4)
    P = [(x, y + rng.normal(0, 0.03)) for x in x3 for y in y3] + [(x4 + rng.normal(0, 0.03), y) for y in y4]
    P = np.array(P)
    g = rng.integers(8)
    m = 4.0
    if g >= 4: P[:, 0] = m - P[:, 0]
    for _ in range(g % 4): P = np.stack([m - P[:, 1], P[:, 0]], 1)
    return np.clip(P, 0, m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', required=True); ap.add_argument('--tag', required=True)
    ap.add_argument('--tries', type=int, default=3000); ap.add_argument('--add', type=int, default=6)
    ap.add_argument('--min-viol', type=float, default=0.005); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--time', type=float, default=3600)
    a = ap.parse_args()
    m = F(4); mf = 4.0
    _, Fam = L.read_family(a.resume); Fset = set(Fam)
    out = f"runs/unavoid13no_{a.tag}"; os.makedirs(out, exist_ok=True)
    log = open(f"{out}/round_log.txt", 'a')
    def say(s):
        print(s, flush=True); log.write(s + '\n'); log.flush()
    say(f"# mode attack: |F0| = {len(Fam)}, {time.ctime()}")
    rng = np.random.default_rng(a.seed)
    t0 = time.time(); hits = 0; added_total = 0; sq = L.Squares(Fam)
    for it in range(a.tries):
        if time.time() - t0 > a.time: break
        P0 = mode_set(rng)
        P, tstar = L.polish_positions(P0, m, [sq], rounds=8, trust=0.15, verbose=False)
        if tstar is None or tstar < -1e-9:
            continue
        hits += 1
        viol, gmin = L.find_violations(P, m, dth_deg=0.5, pitch=0.01, top=400, want=a.add, verbose=False)
        viol = [v for i, v in enumerate(viol) if i < 2 or v[0] <= -a.min_viol]
        if not viol:
            Pex = L.snap_points(P, 1000); cert = f"{out}/cand_{it:05d}.txt"; L.write_cert(cert, m, Pex)
            say(f"  try {it}: F-hitting mode set with NO float violation (grid min {gmin:+.2e}) -> {cert} (run the checker!)")
            continue
        added = 0
        for (f, cx, cy, th) in viol:
            key = L.rationalise_pose(m, cx, cy, th, den=1000)
            for kk in L.d4_images(m, *key):
                if kk not in Fset: Fset.add(kk); Fam.append(kk); added += 1
        added_total += added
        sq = L.Squares(Fam)
        say(f"  try {it}: hit F (t* {tstar:+.1e}); worst violation {viol[0][0]:+.3f}; +{added} squares -> |F| = {len(Fam)} ({time.time()-t0:.0f}s)")
        if hits % 10 == 0:
            L.write_family(f"{out}/F_partial.txt", m, Fam, header=f"mode attack partial, {hits} hits")
    L.write_family(f"{out}/F_final.txt", m, Fam, header=f"mode attack: {hits} F-hitting mode sets in {it+1} tries, +{added_total} squares")
    say(f"# done: {hits} hits in {it+1} tries, +{added_total} squares, |F| = {len(Fam)}, {time.ctime()}")


if __name__ == '__main__':
    main()
