#!/usr/bin/env python3
"""L removal, shove, refill + 1 (10-19; Evan: "273 inexact L removal, shove, and refill with an extra square might also
be worth a few tries").  A one-square insertion spread along a whole unit band instead of one pose hole.

Per donor (n - 1 squares, side s) and cut (cx, cy): the horizontal strip [cy, cy + 1) and the vertical one [cx, cx + 1):
  1. inexact removal: drop the d squares whose centre lies in the strip union (squares that only poke in are kept);
  2. shove: squares whose centre is beyond a strip's midline move outward by delta (x and / or y), box s + delta;
  3. refill: d + 1 axis unit squares, one in the crossing cell, the rest split evenly over the two arms, evenly spaced
     on the strip centre lines, angles per --modes (axis refill always makes a full line of k: smoke 10-19, 24 / 24
     line-blocked), tilt = uniform 10-35 deg, band = the donor's dominant tilt; fq quench (loosen 1.02 or 1.05).
Cuts: arm 'scored' = lcut.best_cuts (cleanest cuts: least kept area in the strips + least removed area outside; at 273
these are edge L's), arm 'uniform' = uniform random (cx, cy).  delta in --deltas.  Line check: layout.full_lines at k.

  lrefill.py --n 274 --donors ID... [--donor-list F] --cuts 10 --random 10 --deltas 0,0.03,0.08 --procs 15 --out DIR
"""
import argparse, json, math, os, random, sys, tempfile, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def refill(sq, s, cx, cy, delta, mode='axis', rng=None):
    keep, d = [], 0
    for x, y, a in sq:
        if cx <= x < cx + 1 or cy <= y < cy + 1:
            d += 1; continue
        keep.append((x + delta * (x >= cx + 0.5), y + delta * (y >= cy + 0.5), a))
    S = s + delta
    w = 1 + delta
    mx, my = cx + w / 2, cy + w / 2
    m = d + 1 - 1                                 # one goes in the crossing cell
    h = m // 2 + (m % 2)
    v = m - h
    if mode == 'band':                             # the donor's dominant tilt (mod 90, 1-degree bins)
        A = [a % 90 for _, _, a in sq if 3 < a % 90 < 87]
        ang = float(np.bincount(np.round(A).astype(int)).argmax()) if A else 45.0
    angle = (lambda: 0.0) if mode == 'axis' else (lambda: ang) if mode == 'band' else (lambda: rng.uniform(10, 35))
    new = [(mx, my, angle())]
    # horizontal arm: x in [0, cx) U [cx + w, S), evenly by length
    for cnt, lo1, hi1, lo2, hi2, horiz in ((h, 0, cx, cx + w, S, True), (v, 0, cy, cy + w, S, False)):
        L1, L2 = max(0.0, hi1 - lo1), max(0.0, hi2 - lo2)
        n1 = int(round(cnt * L1 / (L1 + L2))) if L1 + L2 > 0 else 0
        for c, lo, L in ((n1, lo1, L1), (cnt - n1, lo2, L2)):
            for i in range(c):
                t = lo + (i + 0.5) * L / c
                new.append((t, my, angle()) if horiz else (mx, t, angle()))
    C = [(min(max(x, 0.5), S - 0.5), min(max(y, 0.5), S - 0.5), a) for x, y, a in keep + new]
    return S, C, d


def job(j):
    import hop, mcmin
    from layout import full_lines
    arm, src, n, k, cx, cy, delta, loosen, seed, mode = j
    s, sq = mcmin.load_deg(src)
    S, C, d = refill(sq, s, cx, cy, delta, mode, random.Random(seed))
    if len(C) != n:
        return dict(arm=arm, mode=mode, src=src, cut=(cx, cy), delta=delta, d=d, s=None, err=f'count {len(C)}')
    hop.EXTRA[:] = ['--loosen', str(loosen)]
    tmp = tempfile.mkdtemp()
    r = hop.quench(S, C, tmp)
    if not r:
        return dict(arm=arm, mode=mode, src=src, cut=(cx, cy), delta=delta, d=d, s=None)
    p = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p, r[0], r[1])
    return dict(arm=arm, mode=mode, src=src, cut=(round(cx, 3), round(cy, 3)), delta=delta, loosen=loosen, d=d, s=r[0],
                grid=bool(full_lines(p, k)), sq=r[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--donors', nargs='*', default=[]); ap.add_argument('--donor-list')
    ap.add_argument('--cuts', type=int, default=10); ap.add_argument('--random', type=int, default=10)
    ap.add_argument('--modes', default='tilt,band', help='refill angles: axis (0; makes a full line), tilt (10-35 random), band')
    ap.add_argument('--deltas', default='0,0.03,0.08'); ap.add_argument('--step', type=float, default=0.05)
    ap.add_argument('--procs', type=int, default=15); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    from pk.store import Store
    from pk import manifest
    from jobpool import run_jobs
    import mcmin, scratch, lcut
    scratch.scope('lrefill')
    for d in ('donors', 'seeds'):
        os.makedirs(f'{a.out}/{d}', exist_ok=True)
    man = manifest.start(a.out, 'lrefill', seed=a.seed, args=vars(a), procs=a.procs)
    st = Store(); rng = random.Random(a.seed)
    if a.donor_list:
        a.donors += [l.strip() for l in open(a.donor_list) if l.strip()]
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    D = [float(x) for x in a.deltas.split(',')]
    J = []
    for dn in a.donors:
        src = f'{a.out}/donors/{dn}.txt'; st.get(int(dn)).write(src)
        s, sq = mcmin.load_deg(src)
        assert len(sq) == a.n - 1, (dn, len(sq))
        cuts = [('scored', c[1], c[2]) for c in lcut.best_cuts(sq, s, 2 * k - 1, a.step, a.cuts, 0.15)]
        cuts += [('uniform', rng.uniform(0, s - 1), rng.uniform(0, s - 1)) for _ in range(a.random)]
        for arm, cx, cy in cuts:
            for dl in D:
                for mode in a.modes.split(','):
                    J.append((arm, src, a.n, k, cx, cy, dl, (1.02, 1.05)[rng.randrange(2)], rng.randrange(10 ** 9), mode))
    print(f'{len(a.donors)} donors, {len(J)} quenches', flush=True)
    t0 = time.time()
    R = [r for r in run_jobs(job, J, procs=a.procs, timeout=1200) if isinstance(r, dict)]
    with open(f'{a.out}/seeds.jsonl', 'w') as fo:
        for i, r in enumerate(R):
            if r.get('s') and not r['grid']:
                p = f'{a.out}/seeds/{r["arm"]}_{i}.txt'; mcmin.write_deg(p, r['s'], r['sq']); r['path'] = p
            r.pop('sq', None); fo.write(json.dumps(r) + '\n')
    for arm in ('scored', 'uniform'):
      for mode in a.modes.split(','):
        for dl in D:
            A = [r for r in R if r['arm'] == arm and r['delta'] == dl and r['mode'] == mode]
            ok = [r for r in A if r.get('s')]
            free = sorted(r['s'] for r in ok if not r['grid'])
            print(f'{arm} {mode} delta {dl}: {len(A)} runs, {len(ok)} feasible, line-free {len(free)}; < {k}: '
                  f'{sum(x < k - 1e-9 for x in free)}, < {k}.01: {sum(x < k + .01 for x in free)}, < {k}.02: '
                  f'{sum(x < k + .02 for x in free)}; best {[round(x, 6) for x in free[:4]]}', flush=True)
    print(f'done [{time.time() - t0:.0f}s]')
    man.finish()


if __name__ == '__main__':
    main()
