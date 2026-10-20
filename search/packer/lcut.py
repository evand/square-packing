#!/usr/bin/env python3
"""Inexact L removal (10-18; Evan): a packing at side ~k+1 with n0 squares -> n = n0 - d seeds at side ~k by removing a
horizontal unit strip [cy, cy+1) and a vertical one [cx, cx+1), but keeping squares that only poke gently into the strips
(strict removal of everything centred in the strips, e.g. 211 -> 182, leaves gaps where such squares were).

Per cut: every square's area inside the strip union (7x7 sample); remove the d squares with the most; score = area of
the kept squares inside the strips (the overlap the quench must absorb) + the removed squares' area outside them (the
wasted space).  Enumerate cuts on a --step grid, keep the --cuts best-scoring distinct cuts per donor; shift squares
beyond a strip by -1 (by their centre's side of the strip's midline), fq quench (loosen 1.02 / 1.05), full-line check.
Arm 'random' = cropseed.py's rule (uniform cut, centre-in-strip removal, crop.fix_count) on the same donors, as baseline.

  lcut.py --donors ID... --n 183 --cuts 60 --reps 2 --random 120 --procs 4 --out runs/lcut183
"""
import argparse, collections, json, math, os, random, sys, tempfile, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

_U = (np.arange(7) + 0.5) / 7 - 0.5
_UU = np.array([(u, v) for u in _U for v in _U])


def samples(sq):
    """(N, 49, 2) sample points per square (centre x, y, angle deg)"""
    P = []
    for x, y, a in sq:
        t = math.radians(a); c, s = math.cos(t), math.sin(t)
        P.append(np.c_[x + c * _UU[:, 0] - s * _UU[:, 1], y + s * _UU[:, 0] + c * _UU[:, 1]])
    return np.array(P)


def inside(P, cx, cy):
    """area fraction of each square inside the strip union"""
    X, Y = P[..., 0], P[..., 1]
    return ((X >= cx) & (X < cx + 1) | (Y >= cy) & (Y < cy + 1)).mean(axis=1)


def best_cuts(sq, s, d, step, m, sep):
    P = samples(sq)
    R = []
    for cx in np.arange(0, s - 1 + 1e-9, step):
        for cy in np.arange(0, s - 1 + 1e-9, step):
            f = inside(P, cx, cy)
            o = np.argsort(-f)
            rem, kept = o[:d], o[d:]
            R.append((float(f[kept].sum() + (1 - f[rem]).sum()), float(cx), float(cy), float(f[kept].max())))
    R.sort()
    out = []
    for r in R:
        if all(abs(r[1] - q[1]) + abs(r[2] - q[2]) > sep for q in out):
            out.append(r)
        if len(out) >= m:
            break
    return out


def apply_cut(sq, cx, cy, d):
    P = samples(sq)
    f = inside(P, cx, cy)
    rem = set(np.argsort(-f)[:d].tolist())
    C = []
    for i, (x, y, a) in enumerate(sq):
        if i in rem:
            continue
        C.append((x - (x >= cx + 0.5), y - (y >= cy + 0.5), a))
    return C


def job(j):
    import hop, mcmin
    from layout import full_lines
    arm, src, n, k, cx, cy, loosen, seed = j
    s, sq = mcmin.load_deg(src)
    rng = random.Random(seed)
    if arm == 'random':
        import crop
        cx, cy = rng.uniform(0, s - 1), rng.uniform(0, s - 1)
        C = crop.crop(sq, s, cx, cy)
        C, dr, ad = crop.fix_count(C, n, s - 1)
    else:
        C = apply_cut(sq, cx, cy, len(sq) - n)
    s1 = s - 1
    C = [(min(max(x, 0.5), s1 - 0.5), min(max(y, 0.5), s1 - 0.5), a) for x, y, a in C]
    hop.EXTRA[:] = ['--loosen', str(loosen)]
    tmp = tempfile.mkdtemp()
    r = hop.quench(s1, C, tmp)
    if not r:
        return dict(arm=arm, src=src, cut=(cx, cy), s=None)
    p = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p, r[0], r[1])
    return dict(arm=arm, src=src, cut=(round(cx, 3), round(cy, 3)), loosen=loosen, s=r[0], grid=bool(full_lines(p, k)),
                sq=r[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donors', nargs='*', default=[], help="store ids")
    ap.add_argument('--donor-paths', nargs='*', default=[], help='packing files (e.g. dgen.py donors.json paths)')
    ap.add_argument('--donor-list', help='file: one store id per line (added to --donors)')
    ap.add_argument('--n', type=int, default=183); ap.add_argument('--cuts', type=int, default=60)
    ap.add_argument('--reps', type=int, default=2); ap.add_argument('--random', type=int, default=120)
    ap.add_argument('--step', type=float, default=0.05); ap.add_argument('--sep', type=float, default=0.15)
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    from pk.store import Store
    from pk import manifest
    from jobpool import run_jobs
    import mcmin, scratch
    scratch.scope('lcut')
    os.makedirs(f'{a.out}/donors', exist_ok=True); os.makedirs(f'{a.out}/seeds', exist_ok=True)
    man = manifest.start(a.out, 'lcut', seed=a.seed, args=vars(a), procs=a.procs)
    st = Store(); rng = random.Random(a.seed)
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    J, cutinfo = [], {}
    if a.donor_list:
        a.donors += [l.strip() for l in open(a.donor_list) if l.strip()]
    for dn in a.donors + a.donor_paths:
        if os.path.exists(dn):                        # a file: keep its run / name in the donor label
            src = f'{a.out}/donors/' + dn.strip('/').replace('/', '_')
            import shutil; shutil.copy(dn, src)
        else:
            src = f'{a.out}/donors/{dn.replace(":", "")}.txt'; st.get(int(dn)).write(src)
        s, sq = mcmin.load_deg(src)
        d = len(sq) - a.n
        cuts = best_cuts(sq, s, d, a.step, a.cuts, a.sep)
        cutinfo[dn] = cuts[:10]
        print(f'{dn}: n0 {len(sq)} s {s:.6f} remove {d}; best cut scores {[round(c[0], 3) for c in cuts[:5]]} '
              f'at {[(round(c[1], 2), round(c[2], 2)) for c in cuts[:5]]}', flush=True)
        for sc, cx, cy, mx in cuts:
            for q in range(a.reps):
                J.append(('lcut', src, a.n, k, cx, cy, (1.02, 1.05)[q % 2], rng.randrange(10 ** 9)))
        for t in range(a.random):
            J.append(('random', src, a.n, k, 0, 0, (1.02, 1.05)[t % 2], rng.randrange(10 ** 9)))
    t0 = time.time()
    R = run_jobs(job, J, procs=a.procs, timeout=1200)
    R = [r for r in R if isinstance(r, dict)]
    with open(f'{a.out}/seeds.jsonl', 'w') as fo:
        for i, r in enumerate(R):
            if r.get('s') and not r['grid']:
                p = f'{a.out}/seeds/{r["arm"]}_{i}.txt'; mcmin.write_deg(p, r['s'], r['sq']); r['path'] = p
            r.pop('sq', None); fo.write(json.dumps(r) + '\n')
    json.dump(cutinfo, open(f'{a.out}/cuts.json', 'w'))
    for arm in ('lcut', 'random'):
        A = [r for r in R if r['arm'] == arm]
        ok = [r for r in A if r.get('s')]
        free = sorted(r['s'] for r in ok if not r['grid'])
        print(f'{arm}: {len(A)} seeds, {len(ok)} quenched, line-free {len(free)}; < {k}: {sum(x < k - 1e-9 for x in free)}, '
              f'< {k}.01: {sum(x < k + .01 for x in free)}, < {k}.02: {sum(x < k + .02 for x in free)}; best '
              f'{[round(x, 6) for x in free[:6]]}', flush=True)
    print(f'done [{time.time() - t0:.0f}s]')
    man.finish()


if __name__ == '__main__':
    main()
