#!/usr/bin/env python3
"""t3_exist_zero: the SCOPE of the existence clause (E3) -- the set of angle vectors that admit a
packing of six unit squares in [0,3]^2 at all -- mapped by walking the zero set of delta*, and the
certificate shapes available at the margin-0 packings it finds.

Task tasks/t3-existence/README.md.  Nothing in search/ is modified (t3_chain, t3_chain_hform,
bandcut_k are imported).

Rationale.  (E3) quantifies over configurations with margin delta >= 0, i.e. over genuine
packings.  A packing at angle vector theta exists iff delta*(theta) >= 0, and delta*(theta) <= 0
is the theorem, so the scope of (E3) is exactly Z = {theta : delta*(theta) = 0}.  Every
certificate failure recorded in notes/t3-chain.md and search/T4_CYCLES.md is at a configuration
with delta < 0, hence OUTSIDE the scope.  This script asks what Z actually contains.

Every configuration it reports is a FEASIBLE point (an explicit packing), so "delta >= -1e-9 was
reached" is a one-sided, reliable statement; "the walk got stuck" is not a proof of anything.

usage:
    python3 search/t3_exist_zero.py walk  --out runs/t3_exist_zero1.jsonl [--nwalk 40] [--seed 0]
    python3 search/t3_exist_zero.py report runs/t3_exist_zero1.jsonl
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3_chain as TC                                             # noqa: E402

N, T = 6, 3
HALF = math.pi / 2


def tilt(th):
    """distance of th to 0 mod 90 deg, in radians, in [0, 45 deg]."""
    r = abs(th) % HALF
    return min(r, HALF - r)


def tilts(theta):
    return sorted((tilt(t) for t in theta), reverse=True)


def dstar(theta, rng, limit=200, extra=300, keep=1):
    v, z, tops = TC.delta_star(list(theta), T, rng, limit=limit, extra=extra, keep=keep)
    return v, z, tops


def walk(theta0, rng, steps=300, tol=1e-9, step0=0.20, verbose=False):
    """greedy walk on the zero set: repeatedly try to raise ONE square's tilt while keeping
    delta* >= -tol.  Returns the last angle vector that still admits a packing, and its z."""
    theta = np.array(theta0, dtype=float)
    v, z, _ = dstar(theta, rng)
    if v < -tol:
        return None
    step = step0
    stuck = 0
    for _ in range(steps):
        i = int(rng.integers(N))
        sgn = 1.0 if rng.random() < 0.5 else -1.0
        cand = theta.copy()
        cand[i] += sgn * step
        if tilt(cand[i]) <= tilt(theta[i]) + 1e-12:
            stuck += 1
        vv, zz, _ = dstar(cand, rng)
        if vv >= -tol:
            theta, v, z = cand, vv, zz
            stuck = 0
        else:
            stuck += 1
            if stuck >= 12:
                step *= 0.5
                stuck = 0
                if step < 1e-3:
                    break
    return dict(theta=[float(t) for t in theta], delta=float(v), z=[float(x) for x in z])


def certify(z, lvl):
    """every restricted-dual bound of interest at the configuration z, at level lvl."""
    A, b, tags = TC.all_rows(z, N, T, lvl)
    TH = np.array(z[3::3])
    out = {}
    bh = TC.best_h(z, N, T, lvl, modes=('C', 'CW', 'H', 'HP'), rows=(A, b, tags))
    for m in ('C', 'CW', 'H', 'HP'):
        out[m] = bh[m]['bound']
        out[m + '_path'] = bh[m]['path']
        out[m + '_axis'] = bh[m]['axis']
    out['nchains'] = bh['nchains']
    out['F'], w = TC.dual_bound(A, b, range(len(b)), want_w=True)
    supp = [(float(w[q]), tags[q]) for q in range(len(b)) if w[q] > 1e-7]
    supp.sort(key=lambda t: -t[0])
    out['supp'] = [[s[0], str(s[1])] for s in supp]
    # axis-parallel-normal chain of three: the Lemma Z support
    out['Zchain'] = zchain(z, tags, TH, lvl)
    return out


def zchain(z, tags, TH, lvl):
    """is there a wall-to-wall chain of three whose two link normals are EXACTLY axis-parallel
    (tilt of the normal-owning square = 0 mod 90 deg)?  -> best Lemma-Z value, or None."""
    best = None
    chs = TC.chains_of(z, N, lvl, T)
    for (ax, path) in chs:
        ok = True
        rows = []
        for s in range(len(path) - 1):
            idx = TC.link_rows(tags, TH, ax, path[s], path[s + 1], None)
            good = [q for q in idx if tilt(TH[tags[q][3]]) < 1e-9]
            if not good:
                ok = False
                break
            rows.append(good[0])
        if not ok:
            continue
        u = np.abs(np.cos(TH)) + np.abs(np.sin(TH))
        m = sum(0.5 + 0.5 * (abs(math.cos(TH[tags[q][2]] - TH[tags[q][1]]))
                             + abs(math.sin(TH[tags[q][2]] - TH[tags[q][1]]))) for q in rows)
        val = (T - u[path[0]] / 2 - u[path[-1]] / 2 - m) / (T + 1)
        if best is None or val < best[0]:
            best = (float(val), ax, list(path))
    return best


def cmd_walk(a):
    rng = np.random.default_rng(a.seed)
    out = open(a.out, 'w')
    starts = []
    # start from the axis-parallel point and from every k-zero stratum
    for k in range(7):
        for _ in range(max(1, a.nwalk // 7)):
            th = np.zeros(N)
            free = rng.permutation(N)[:N - k]
            th[free] = rng.uniform(0, HALF, N - k)
            starts.append(th)
    for si, th0 in enumerate(starts):
        v0, z0, _ = dstar(th0, rng)
        rec = None
        if v0 >= -1e-9:
            rec = walk(th0, rng, steps=a.steps)
        if rec is None:
            # start is not in Z: walk from the axis-parallel point instead
            rec = walk(np.zeros(N), rng, steps=a.steps)
        if rec is None:
            continue
        z = np.array(rec['z'])
        rec['tilts_deg'] = [math.degrees(t) for t in tilts(rec['theta'])]
        rec['ntilted_1deg'] = sum(1 for t in rec['theta'] if tilt(t) > math.radians(1))
        rec['cert'] = certify(z, rec['delta'])
        out.write(json.dumps(rec) + '\n')
        out.flush()
        print(f"[{si+1}/{len(starts)}] delta={rec['delta']:+.3e} "
              f"tilts={[round(x,2) for x in rec['tilts_deg']]} "
              f"C={rec['cert']['C']:+.4f} CW={rec['cert']['CW']:+.4f} H={rec['cert']['H']:+.4f} "
              f"Z={'-' if rec['cert']['Zchain'] is None else round(rec['cert']['Zchain'][0],4)}",
              flush=True)
    out.close()


def cmd_report(path):
    rows = [json.loads(l) for l in open(path)]
    rows = [r for r in rows if r['delta'] >= -1e-9]
    print(f"margin-0 packings found: {len(rows)}")
    print(f"{'#tilted>1deg':>12} {'n':>4} {'max 4th tilt':>13} {'C<=0':>7} {'CW<=0':>7} "
          f"{'H<=0':>7} {'Zchain':>8}")
    for k in range(7):
        rs = [r for r in rows if r['ntilted_1deg'] == k]
        if not rs:
            continue
        m4 = max(r['tilts_deg'][3] for r in rs)
        print(f"{k:>12} {len(rs):>4} {m4:13.3f} "
              f"{sum(1 for r in rs if r['cert']['C'] <= 1e-9):>7} "
              f"{sum(1 for r in rs if r['cert']['CW'] <= 1e-9):>7} "
              f"{sum(1 for r in rs if r['cert']['H'] <= 1e-9):>7} "
              f"{sum(1 for r in rs if r['cert']['Zchain'] is not None):>8}")
    print("\nmost-tilted packings found (by 3rd largest tilt):")
    rows.sort(key=lambda r: -r['tilts_deg'][2])
    for r in rows[:12]:
        print(f"  delta={r['delta']:+.2e} tilts={[round(x,2) for x in r['tilts_deg']]} "
              f"C={r['cert']['C']:+.5f} CW={r['cert']['CW']:+.5f} H={r['cert']['H']:+.5f} "
              f"F={r['cert']['F']:+.5f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    w = sub.add_parser('walk')
    w.add_argument('--out', default='runs/t3_exist_zero1.jsonl')
    w.add_argument('--nwalk', type=int, default=35)
    w.add_argument('--steps', type=int, default=200)
    w.add_argument('--seed', type=int, default=0)
    r = sub.add_parser('report')
    r.add_argument('path')
    a = ap.parse_args()
    if a.cmd == 'walk':
        cmd_walk(a)
    else:
        cmd_report(a.path)


if __name__ == '__main__':
    main()
