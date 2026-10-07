#!/usr/bin/env python3
"""LP-only squeeze (no penalty relax): trust-region sequential LP on  min s  s.t. linearised separation.

Step: LP (rigid.shrink_lp, contacts = pairs/walls with gap < 3R, |dz| <= R)  ->  apply  ->  measure true worst overlap p
(SAT, all near pairs, walls)  ->  repair by scaling configuration and side by (1 + p') about the origin corner, where p' just
clears every overlap  ->  accept if the side decreased, else R /= 3.  Converges to a first-order jammed state.

Contact model (default model='inc', inc.py): corner-on-side incidences as in exact/find_contacts; corner-corner touches are
disjunctions (LP relax-and-round when the LP stalls; exact MILP for the final jam test ds0).  The old model='axis' (rigid.py,
one separating axis per pair) forbids feasible rotations/slides and jams falsely; kept for comparison (--axis).
Second-order correction (soc): after a step, the side is held and the linearised contacts are restored (boxed at ~10x the
overlap) before the global rescaling repair, so curvature overlaps cost O(R^4) instead of O(R^2) in side.
info['jammed'] = ds0 > -1e-8 at R = 1e-3, i.e. dS/|z| > -1e-5 (an f64 point; exact/ certifies).

Wall clock: `budget` seconds (default none) caps the whole call; info['timeout'] is set if it ran out.

  slp2.py in.txt out.txt [--R 1e-2] [--iters 2000] [-v]
"""
import math, sys, time
from rigid import load as _l  # noqa: pins BLAS threads before numpy loads
import numpy as np
from rigid import load, shrink_lp, shrink_switch, pair_rows
import inc

H = 0.5


def worst(s, sq):
    """Largest overlap: max over pairs of -gap (SAT) and over walls of protrusion; <= 0 means feasible."""
    w = -1.0
    n = len(sq)
    for i in range(n):
        x, y, t = sq[i]
        e = H * (abs(math.cos(t)) + abs(math.sin(t)))
        w = max(w, e - x, x + e - s, e - y, y + e - s)
    for i in range(n):
        for j in range(i + 1, n):
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 < 2.0:
                w = max(w, -pair_rows(sq[i], sq[j])[0])
    return w


def repair(s, sq):
    """Scale about the origin until feasible; returns (s', sq')."""
    p = worst(s, sq)
    if p <= 0:
        return s, sq
    f = 1.0 + 1.01 * p           # pair distances ~>= 1 so relative growth p clears overlap p (walls scale with s)
    for _ in range(30):
        s2 = s * f
        sq2 = [(x * f, y * f, t) for x, y, t in sq]
        if worst(s2, sq2) <= 0:
            return s2, sq2
        f = 1.0 + 2.0 * (f - 1.0)
    return float('inf'), sq


def slp2(s, sq, R=1e-2, iters=2000, verbose=False, rmin=1e-10, fixed=(), budget=None, info=None, milp=True, model='inc', soc=True):
    info = {} if info is None else info
    info.update(timeout=False, switch_steps=0, iters=0)
    t0 = time.time()
    s, sq = repair(s, sq)

    def step(ds, z):
        new = [(float(x + z[3 * i]), float(y + z[3 * i + 1]), float(t + z[3 * i + 2])) for i, (x, y, t) in enumerate(sq)]
        return repair(float(s + ds), new)

    it = 0
    for it in range(iters):
        if budget is not None and time.time() - t0 > budget:
            info['timeout'] = True
            break
        via_milp = False
        if model == 'inc':
            ds, z, nc, nd = inc.shrink(s, sq, 3 * R, R, fixed=fixed, switch=False)
            if (ds is None or ds > -1e-4 * R) and milp and nd:
                dm, zm, _, _ = inc.shrink(s, sq, 3 * R, R, fixed=fixed, switch=True)
                if dm is not None and dm < -1e-14 and (ds is None or dm < 2 * ds):
                    ds, z, via_milp = dm, zm, True
        else:
            ds, z, nc = shrink_lp(s, sq, delta=3 * R, R=R, fixed=fixed)
        if model != 'inc' and (ds is None or ds > -1e-4 * R) and milp:
            dm, zm, nd = shrink_switch(s, sq, delta=3 * R, R=R, fixed=fixed)
            if dm is not None and dm < -1e-14 and (ds is None or dm < 2 * ds):
                ds, z, via_milp = dm, zm, True
        if ds is None or ds > -1e-14:
            if R > rmin:
                R /= 3; continue
            break
        if model == 'inc' and soc:
            new = [(float(x + z[3 * i]), float(y + z[3 * i + 1]), float(t + z[3 * i + 2])) for i, (x, y, t) in enumerate(sq)]
            for _ in range(2):
                if worst(s + ds, new) <= 0:
                    break
                c = inc.correct(s + ds, new, R, fixed)
                if c is None:
                    break
                new = c
            s2, new = repair(float(s + ds), new)
        else:
            s2, new = step(ds, z)
        if s2 < s - 1e-15:
            s, sq = s2, new
            info['switch_steps'] += via_milp
            R = min(R * 1.5, 0.05)
            if verbose and (it % 10 == 0 or via_milp):
                print(f'it={it} s={s:.10f} R={R:.1e} contacts={nc}{" SWITCH" if via_milp else ""} t={time.time() - t0:.1f}s',
                      flush=True)
        else:
            R /= 3
            if R < rmin:
                break
    info['iters'] = it + 1
    if model == 'inc':
        # final first-order jam test, exhaustive over corner-corner branches (MILP; ~0.1-2 s at n = 110)
        ds0, _, nd = inc.shrink_milp(s, sq, delta=1e-7, R=1e-3, exact=True, fixed=fixed, time_limit=60.0)
        info['n_disj'] = nd
        if ds0 is None or nd == 0:
            ds0, _, _, _ = inc.shrink(s, sq, 1e-7, 1e-3, fixed=fixed, exact=True, switch=False)
    else:
        ds0, _, nc = shrink_lp(s, sq, delta=1e-7, exact=True, fixed=fixed)
        if milp and ds0 is not None:
            dm, _, nd = shrink_switch(s, sq, delta=1e-7, R=1e-3, exact=True, fixed=fixed)
            info['n_disj'] = nd
            if dm is not None:
                ds0 = min(ds0, dm)
    info['jammed'] = bool(ds0 is not None and ds0 > -1e-5 * 1e-3)    # dS/|z| > -1e-5 (f64 point; exact/ certifies)
    info['time'] = time.time() - t0
    return float(s), [(float(x), float(y), float(t)) for x, y, t in sq], ds0


if __name__ == '__main__':
    from slp import save
    s, sq = load(sys.argv[1])
    R = float(sys.argv[sys.argv.index('--R') + 1]) if '--R' in sys.argv else 1e-2
    t0 = time.time()
    info = {}
    s2, sq2, ds0 = slp2(s, sq, R, verbose='-v' in sys.argv, milp='--no-milp' not in sys.argv, info=info,
                        model='axis' if '--axis' in sys.argv else 'inc', soc='--no-soc' not in sys.argv)
    save(sys.argv[2], s2, sq2)
    print(f'{sys.argv[1]}: {s:.10f} -> {s2:.10f}; exact-contact ds* = {ds0:.2e}; {info}; {time.time() - t0:.1f}s')
