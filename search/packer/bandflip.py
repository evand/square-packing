#!/usr/bin/env python3
"""Band flip (10-21; Evan: "trying same-tilt more ... Same-tilt also gives a very different set of junction patterns").
Design surgery on a real packing instead of construction from scratch (bulbgen's idealised bands hold ~250 of 274 at
side < k and grid-lock when compressed, with either upright filler: PACKER 10-21 "multi-band bulbgen").

Find the tilted bands (clusters of >= 15 squares of similar signed tilt), pick one, find its wall (nearest wall to its
centroid), and reflect a strip along that wall across the perpendicular line through the strip's middle: every square
with depth (distance to the wall) < D and wall coordinate in [lo, hi] maps t -> lo + hi - t, angle -> -angle.  The band
keeps its wall and its width but its tilt sign flips, so a mixed pair becomes same-tilt (and the reverse), with the
band's ends swapped: the junction it meets is new.  D = band depth + U(0, 1.5); [lo, hi] = the band's wall extent
+ U(0, 1.5) margins, clipped at the other bands' depth on the perpendicular walls (corner blocks stay).  fq quench
(loosen 1.02 / 1.05), full-line check.  Mode 'both' flips every band (the mirror image, a control: same design).

  bandflip.py --donors ID ... --reps 40 --procs 8 --out runs/bandflip/f274
"""
import argparse, json, math, os, random, sys, tempfile, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def bands(sq, min_size=15):
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial.distance import cdist
    sq = np.asarray(sq, float)
    a = (sq[:, 2] + 45) % 90 - 45
    T = np.where(np.abs(a) > 3)[0]
    if not len(T):
        return []
    P, A = sq[T, :2], a[T]
    adj = (cdist(P, P) < 1.6) & (np.sign(A)[:, None] == np.sign(A)[None, :]) & (np.abs(A[:, None] - A[None, :]) < 8)
    nc, lab = connected_components(adj, directed=False)
    out = [T[lab == j] for j in range(nc) if (lab == j).sum() >= min_size]
    return sorted(out, key=len, reverse=True)


def wall_frame(s, w):
    """(depth, along) functions and the inverse for wall w (0 right, 1 bottom, 2 left, 3 top)"""
    if w == 0: return (lambda x, y: (s - x, y)), (lambda d, t: (s - d, t))
    if w == 1: return (lambda x, y: (y, x)), (lambda d, t: (t, d))
    if w == 2: return (lambda x, y: (x, y)), (lambda d, t: (d, t))
    return (lambda x, y: (s - y, x)), (lambda d, t: (t, s - d))


def flip(s, sq, idx, rng, B):
    sq = [tuple(q) for q in np.asarray(sq, float)]
    C = np.asarray(sq)[idx, :2].mean(0)
    w = int(np.argmin([s - C[0], C[1], C[0], s - C[1]]))
    fwd, inv = wall_frame(s, w)
    DT = np.array([fwd(sq[i][0], sq[i][1]) for i in idx])
    D = DT[:, 0].max() + 0.5 + rng.uniform(0, 1.5)
    lo = DT[:, 1].min() - 0.5 - rng.uniform(0, 1.5); hi = DT[:, 1].max() + 0.5 + rng.uniform(0, 1.5)
    # clip at other bands sitting on the perpendicular walls near either end (their corner blocks stay put)
    for J in B:
        if J is idx: continue
        dt = np.array([fwd(sq[i][0], sq[i][1]) for i in J])
        if dt[:, 0].min() < D:                       # the other band reaches into our strip
            if dt[:, 1].mean() < (lo + hi) / 2: lo = max(lo, dt[:, 1].max() + 0.5)
            else: hi = min(hi, dt[:, 1].min() - 0.5)
    lo, hi = max(lo, 0.0), min(hi, s)
    out, nflip = [], 0
    for x, y, a in sq:
        d, t = fwd(x, y)
        if d < D and lo <= t <= hi:
            x, y = inv(d, lo + hi - t); a = -a; nflip += 1
        out.append((x, y, a % 90))
    return out, dict(wall=w, D=round(D, 2), lo=round(lo, 2), hi=round(hi, 2), nflip=nflip)


def job(j):
    import hop, mcmin
    from layout import full_lines
    src, mode, bi, loosen, seed = j
    s, sq = mcmin.load_deg(src)
    n = len(sq); k = math.ceil(math.sqrt(n) - 1e-12)
    rng = random.Random(seed)
    B = bands(sq)
    if len(B) <= bi:
        return dict(src=src, mode=mode, s=None, err='no band')
    if mode == 'both':
        info = []
        for J in B:
            sq, inf = flip(s, sq, J, rng, B); info.append(inf)
    else:
        sq, info = flip(s, sq, B[bi], rng, B)
    hop.EXTRA[:] = ['--loosen', str(loosen)]
    tmp = tempfile.mkdtemp()
    r = hop.quench(s, sq, tmp)
    base = dict(src=src, mode=mode, band=bi, info=info, loosen=loosen, s_par=s, n=n)
    if not r:
        return dict(base, s=None)
    p = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p, r[0], r[1])
    signs = sorted('+' if np.median([(r[1][i][2] + 45) % 90 - 45 for i in J]) > 0 else '-' for J in bands(r[1]))
    return dict(base, s=r[0], grid=bool(full_lines(p, k)), signs=''.join(signs), sq=r[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--donors', nargs='*', default=[]); ap.add_argument('--donor-list')
    ap.add_argument('--reps', type=int, default=20, help='per donor and band')
    ap.add_argument('--both', type=int, default=4, help="mirror-image controls per donor")
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    from pk.store import Store
    from pk import manifest
    from jobpool import run_jobs
    import mcmin, scratch
    scratch.scope('bandflip')
    os.makedirs(f'{a.out}/donors', exist_ok=True); os.makedirs(f'{a.out}/seeds', exist_ok=True)
    man = manifest.start(a.out, 'bandflip', seed=a.seed, args=vars(a), procs=a.procs)
    st = Store(); rng = random.Random(a.seed)
    if a.donor_list:
        a.donors += [l.split()[0] for l in open(a.donor_list) if l.strip()]
    J = []
    for dn in a.donors:
        src = f'{a.out}/donors/{dn}.txt'; st.get(int(dn)).write(src)
        s, sq = mcmin.load_deg(src)
        B = bands(sq)
        print(f'{dn}: n {len(sq)} s {s:.6f} bands {[len(b) for b in B]}', flush=True)
        for bi in range(min(len(B), 2)):
            for q in range(a.reps):
                J.append((src, 'one', bi, (1.02, 1.05)[q % 2], rng.randrange(10 ** 9)))
        for q in range(a.both):
            J.append((src, 'both', 0, (1.02, 1.05)[q % 2], rng.randrange(10 ** 9)))
    t0 = time.time()
    R = [r for r in run_jobs(job, J, procs=a.procs, timeout=1200) if isinstance(r, dict)]
    with open(f'{a.out}/seeds.jsonl', 'w') as fo:
        for i, r in enumerate(R):
            if r.get('s') and not r['grid']:
                p = f'{a.out}/seeds/{r["mode"]}{r["band"]}_{i}.txt'; mcmin.write_deg(p, r['s'], r['sq']); r['path'] = p
            r.pop('sq', None); fo.write(json.dumps(r) + '\n')
    import collections
    G = collections.defaultdict(list)
    for r in R:
        G[(r.get('n'), r['mode'])].append(r)
    for (n, mode), A in sorted(G.items(), key=lambda kv: (kv[0][0] or 0, kv[0][1])):
        free = sorted((r['s'], r['signs'], r['s_par']) for r in A if r.get('s') and not r['grid'])
        same = [f for f in free if len(f[1]) >= 2 and len(set(f[1])) == 1]
        print(f'n {n} {mode}: {len(A)} seeds, line-free {len(free)}, of which same-tilt {len(same)}; best '
              f'{[(round(f[0], 5), f[1]) for f in free[:4]]}; best same-tilt {[(round(f[0], 5), round(f[0] - f[2], 4)) for f in same[:4]]}',
              flush=True)
    print(f'done [{time.time() - t0:.0f}s]')
    man.finish()


if __name__ == '__main__':
    main()


def design_class(s, sq):
    """Band-layout class (10-21, the design bandit): bands = same-sign tilted components >= 15 squares, each split into
    wall arms (nearest wall, >= 6 squares).  'none' / 'one' (one straight band) / 'bent' (one band with arms >= 10 on
    two walls, 211-style) / 'adj-mixed' 'adj-same' 'opp-mixed' 'opp-same' (two bands: main walls adjacent or
    opposite, tilt signs) / '3+' (three or more bands)."""
    sq = np.asarray(sq, float)
    a = (sq[:, 2] + 45) % 90 - 45
    B = [J for J in bands(sq) if abs(np.median(a[J])) >= 8]   # 10-21: drop barely tilted clusters (a -4 deg
    if not B:                                                    # cluster made a same-tilt pair read '3+')
        return 'none'
    walls = []
    for J in B:
        X = sq[J, :2]
        w = np.c_[s - X[:, 0], X[:, 1], X[:, 0], s - X[:, 1]].argmin(1)
        cnt = np.bincount(w, minlength=4)
        walls.append((cnt, int(cnt.argmax()), np.sign(np.median(a[J]))))
    if len(B) >= 3:
        return '3+'
    if len(B) == 1:
        return 'bent' if (walls[0][0] >= 10).sum() >= 2 else 'one'
    (_, w1, g1), (_, w2, g2) = walls
    rel = 'opp' if (w1 - w2) % 2 == 0 and w1 != w2 else 'adj' if w1 != w2 else 'same-wall'
    return f"{rel}-{'same' if g1 == g2 else 'mixed'}"


def relocate(s, sq, idx, rng, B):
    """Move a band to the opposite wall (10-21; Evan: bands on opposite sides / > 2 bands, let the bandit decide): the
    band's wall strip (depth D, wall range [lo, hi] as in flip) and its 180-degree image on the opposite wall are
    swapped by the rotation (x, y) -> (s - x, s - y), which keeps tilt signs.  R = strip u rot(strip) is rot-invariant,
    so this is a permutation of the squares in R."""
    sq = [tuple(q) for q in np.asarray(sq, float)]
    C = np.asarray(sq)[idx, :2].mean(0)
    w = int(np.argmin([s - C[0], C[1], C[0], s - C[1]]))
    fwd, inv = wall_frame(s, w)
    DT = np.array([fwd(sq[i][0], sq[i][1]) for i in idx])
    D = DT[:, 0].max() + 0.5 + rng.uniform(0, 1.0)
    lo = max(0.0, DT[:, 1].min() - 0.5 - rng.uniform(0, 1.0)); hi = min(s, DT[:, 1].max() + 0.5 + rng.uniform(0, 1.0))
    def inR(x, y):
        d, t = fwd(x, y)
        if d < D and lo <= t <= hi:
            return True
        d, t = fwd(s - x, s - y)
        return d < D and lo <= t <= hi
    out, nmov = [], 0
    for x, y, a in sq:
        if inR(x, y):
            x, y = s - x, s - y; nmov += 1
        out.append((x, y, a % 90))
    return out, dict(wall=w, D=round(D, 2), lo=round(lo, 2), hi=round(hi, 2), nmov=nmov)
