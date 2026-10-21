#!/usr/bin/env python3
"""Approximate L addition (10-21; Evan: "It's maybe even 240/241 and approximate L addition moves."): the inexact inverse
of lcut.py.  A donor (n0, s) -> n = n0 + d seeds at side ~s + 1: open a vertical strip [cx, cx + w) and a horizontal one
[cy, cy + w) (squares shift by +w by their centre's side of the cut; w = U(1, wmax)), so the empty cross has room for
about floor(s + w) + floor(s) axis unit squares, then lay the d new squares along the cross's centre lines at uniform
spacing (len / count: < 1 when d exceeds the axis capacity, i.e. overlap the quench must absorb), each strip with its own
tilt (0 or U(+-tmax), + per-square jitter) so the refill is staggered / tilted, not a rebuilt grid row.  fq quench
(loosen 1.02 / 1.05), full-line check.  Exact L addition (Friedman L-step) is the w = 1, axis, cut-at-the-edge case.

Cuts: arm 'scored' = the --cuts cleanest distinct cuts on a --step grid (cleanliness = area of squares on the far side of
a cut line from their centre, i.e. what the shift tears), arm 'random' = uniform cuts (edge cuts included).

  ladd.py --donors 291673 --n 274 --cuts 40 --reps 2 --random 80 --procs 8 --out runs/ladd/a274
Useful pairs (store records 10-21): 241 (15.98712) -> 274 (d 33, 16.987 + squeeze), 240 -> 274 (34), 273 -> 308 (35),
272 -> 307 / 308.
"""
import argparse, json, math, os, random, sys, tempfile, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lcut import samples


def tear(P, cx, cy):
    """per square: area fraction on the other side of the cut lines from its centre"""
    X, Y = P[..., 0], P[..., 1]
    Cx, Cy = X.mean(axis=1, keepdims=True), Y.mean(axis=1, keepdims=True)
    return (((X >= cx) != (Cx >= cx)) | ((Y >= cy) != (Cy >= cy))).mean(axis=1)


def best_cuts(sq, s, step, m, sep):
    P = samples(sq)
    R = sorted((float(tear(P, cx, cy).sum()), float(cx), float(cy))
               for cx in np.arange(0, s + 1e-9, step) for cy in np.arange(0, s + 1e-9, step))
    out = []
    for r in R:
        if all(abs(r[1] - q[1]) + abs(r[2] - q[2]) > sep for q in out):
            out.append(r)
        if len(out) >= m:
            break
    return out


def add_l(sq, s, d, cx, cy, w, rng, tmax):
    """shifted donor + d new squares in the cross; -> (S, squares)"""
    S = s + w
    C = [(x + w * (x >= cx), y + w * (y >= cy), a) for x, y, a in sq]
    # segments of the cross centre lines: vertical x = cx + w/2 over [0, S]; horizontal y = cy + w/2 over [0, cx] and
    # [cx + w, S] (the crossing cell belongs to the vertical one)
    segs = [('v', 0.0, S), ('h', 0.0, cx), ('h', cx + w, S)]
    L = [b - a_ for _, a_, b in segs]
    cnt = [0, 0, 0]
    for _ in range(d):                                   # largest remaining length per square (D'Hondt)
        j = max(range(3), key=lambda j: L[j] / (cnt[j] + 1) if L[j] > 0.5 else -1)
        cnt[j] += 1
    for (o, a0, b0), c in zip(segs, cnt):
        if not c:
            continue
        tilt = 0.0 if rng.random() < 0.3 else rng.uniform(-tmax, tmax)
        h = (b0 - a0) / c
        for i in range(c):
            t = a0 + h * (i + 0.5) + rng.gauss(0, 0.05 * h)
            off = rng.gauss(0, 0.03)
            ang = tilt + rng.gauss(0, 1.0)
            C.append((cx + w / 2 + off, t, ang) if o == 'v' else (t, cy + w / 2 + off, ang))
    return S, C


def job(j):
    import hop, mcmin
    from layout import full_lines
    arm, src, n, k, cx, cy, loosen, seed, wmax, tmax = j
    s, sq = mcmin.load_deg(src)
    rng = random.Random(seed)
    if arm == 'random':
        cx, cy = rng.uniform(0, s), rng.uniform(0, s)
    w = rng.uniform(1.0, wmax)
    S, C = add_l(sq, s, n - len(sq), cx, cy, w, rng, tmax)
    C = [(min(max(x, 0.5), S - 0.5), min(max(y, 0.5), S - 0.5), a) for x, y, a in C]
    hop.EXTRA[:] = ['--loosen', str(loosen)]
    tmp = tempfile.mkdtemp()
    r = hop.quench(S, C, tmp)
    base = dict(arm=arm, src=src, cut=(round(cx, 3), round(cy, 3)), w=round(w, 4), loosen=loosen)
    if not r:
        return dict(base, s=None)
    p = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p, r[0], r[1])
    return dict(base, s=r[0], grid=bool(full_lines(p, k)), sq=r[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donors', nargs='*', default=[], help='store ids')
    ap.add_argument('--donor-list', help='file: one store id per line')
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--cuts', type=int, default=40)
    ap.add_argument('--reps', type=int, default=2); ap.add_argument('--random', type=int, default=80)
    ap.add_argument('--step', type=float, default=0.1); ap.add_argument('--sep', type=float, default=0.3)
    ap.add_argument('--wmax', type=float, default=1.06); ap.add_argument('--tmax', type=float, default=20.0)
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    from pk.store import Store
    from pk import manifest
    from jobpool import run_jobs
    import mcmin, scratch
    scratch.scope('ladd')
    os.makedirs(f'{a.out}/donors', exist_ok=True); os.makedirs(f'{a.out}/seeds', exist_ok=True)
    man = manifest.start(a.out, 'ladd', seed=a.seed, args=vars(a), procs=a.procs)
    st = Store(); rng = random.Random(a.seed)
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    if a.donor_list:
        a.donors += [l.split()[0] for l in open(a.donor_list) if l.strip()]
    J = []
    for dn in a.donors:
        src = f'{a.out}/donors/{dn}.txt'; st.get(int(dn)).write(src)
        s, sq = mcmin.load_deg(src)
        cuts = best_cuts(sq, s, a.step, a.cuts, a.sep)
        print(f'{dn}: n0 {len(sq)} s {s:.6f} add {a.n - len(sq)}; cleanest cuts {[round(c[0], 2) for c in cuts[:5]]} '
              f'at {[(round(c[1], 2), round(c[2], 2)) for c in cuts[:5]]}', flush=True)
        for sc, cx, cy in cuts:
            for q in range(a.reps):
                J.append(('scored', src, a.n, k, cx, cy, (1.02, 1.05)[q % 2], rng.randrange(10 ** 9), a.wmax, a.tmax))
        for t in range(a.random):
            J.append(('random', src, a.n, k, 0, 0, (1.02, 1.05)[t % 2], rng.randrange(10 ** 9), a.wmax, a.tmax))
    t0 = time.time()
    R = [r for r in run_jobs(job, J, procs=a.procs, timeout=1200) if isinstance(r, dict)]
    with open(f'{a.out}/seeds.jsonl', 'w') as fo:
        for i, r in enumerate(R):
            if r.get('s') and not r['grid']:
                p = f'{a.out}/seeds/{r["arm"]}_{i}.txt'; mcmin.write_deg(p, r['s'], r['sq']); r['path'] = p
            r.pop('sq', None); fo.write(json.dumps(r) + '\n')
    for arm in ('scored', 'random'):
        A = [r for r in R if r['arm'] == arm]
        ok = [r for r in A if r.get('s')]
        free = sorted(r['s'] for r in ok if not r['grid'])
        print(f'{arm}: {len(A)} seeds, {len(ok)} quenched, line-free {len(free)}; < {k}: {sum(x < k - 1e-9 for x in free)}, '
              f'< {k}.01: {sum(x < k + .01 for x in free)}, < {k}.03: {sum(x < k + .03 for x in free)}; best '
              f'{[round(x, 6) for x in free[:6]]}', flush=True)
    print(f'done [{time.time() - t0:.0f}s]')
    man.finish()


if __name__ == '__main__':
    main()
