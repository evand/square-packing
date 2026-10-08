#!/usr/bin/env python3
"""Recombination of good basins (10-08, Evan): parents of the same type, squares matched to common IDs, each ID's start
pose taken from one parent.  Standalone comparison of variants before anything goes into the explorer.

Parents: sub-k basins of an explorer archive.  Per trial: a random reference R and the P - 1 archive basins nearest to R
in matching distance (best of the 8 container symmetries, Hungarian on the corner metric of rediscover.py; of a random
sample of candidates), each aligned to R's square IDs and scaled to R's side.  Variants (how each ID picks a parent):
  kick        control: one kick of R (mcmin kick), no recombination
  uniform     independent uniform parent per ID
  bestfit:K   random ID order; K parents sampled per ID, keep the pose with the least overlap with squares already placed
  field       smooth random field per parent (3 random plane waves, wavelength s/4..s), ID takes the argmax parent: patches
  softfield   as field, parent sampled with probability ~ exp(4 * field): patches with mixed borders
  half        random line: one side from R, the other from one other parent (cut-and-splice control)
Pipeline per proposal = the explorer's: fq screen (<= 8 SLP it., no flips, loosen 1.02) -> discard (> k + 0.05) / grid ->
full polish.  Scored against every sub-k side seen before (certified census + explorer archives): new, parent return.

  cross.py --archive runs/ex10 --n 110 --trials 150 --P 5 --variants kick uniform bestfit:3 field softfield half --out runs/cx1
  cross.py --report runs/cx1
"""
import argparse, bisect, collections, glob, json, math, os, random, tempfile, time
import numpy as np
from scipy.optimize import linear_sum_assignment
from jobpool import run_jobs
import mcmin, hop, explore
from rediscover import sym

HERE = os.path.dirname(os.path.abspath(__file__))


def align(R, sR, Q, sQ):
    """Q's squares reordered to R's IDs (best of 8 symmetries), scaled to R's side; returns (rms distance, poses)."""
    f = sR / sQ
    A = np.array(R)
    best = (float('inf'), None)
    for k in range(8):
        B = np.array([(x * f, y * f, a) for x, y, a in sym(Q, sQ, k)])
        dx = A[:, None, 0] - B[None, :, 0]; dy = A[:, None, 1] - B[None, :, 1]
        da = (A[:, None, 2] - B[None, :, 2] + 45) % 90 - 45
        C = dx * dx + dy * dy + (0.5 * math.sqrt(2) * np.radians(da)) ** 2
        r, c = linear_sum_assignment(C)
        d = math.sqrt(C[r, c].mean())
        if d < best[0]:
            P = [None] * len(R)
            for i, j in zip(r, c):
                P[i] = tuple(float(z) for z in B[j])
            best = (d, P)
    return best


def corners(x, y, a):
    t = math.radians(a); c, s = math.cos(t) / 2, math.sin(t) / 2
    return np.array([(x + c * u - s * v, y + s * u + c * v) for u, v in ((1, 1), (-1, 1), (-1, -1), (1, -1))])


def sat_pen(p, q):
    """Penetration depth of two unit squares (0 if separated): min over the 4 face normals of the overlap of projections."""
    if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > 2.0:
        return 0.0
    P, Q = corners(*p), corners(*q)
    best = float('inf')
    for a in (p[2], q[2]):
        t = math.radians(a)
        for u in ((math.cos(t), math.sin(t)), (-math.sin(t), math.cos(t))):
            pp, qq = P @ u, Q @ u
            ov = min(pp.max(), qq.max()) - max(pp.min(), qq.min())
            if ov <= 0:
                return 0.0
            best = min(best, ov)
    return best


def wall_pen(p, s):
    C = corners(*p)
    return max(0.0, -C[:, 0].min(), -C[:, 1].min(), C[:, 0].max() - s, C[:, 1].max() - s)


def recombine(variant, s, parents, rng):
    n = len(parents[0]); P = len(parents)
    if variant == 'uniform':
        return [parents[rng.randrange(P)][i] for i in range(n)]
    if variant.startswith('bestfit'):
        K = int(variant.split(':')[1])
        order = list(range(n)); rng.shuffle(order)
        placed = []                                        # (pose)
        out = [None] * n
        for i in order:
            cands = rng.sample(range(P), min(K, P))
            sc = []
            for pidx in cands:
                q = parents[pidx][i]
                pen = wall_pen(q, s) + sum(sat_pen(q, r) for r in placed)
                sc.append((pen, rng.random(), q))
            q = min(sc)[2]
            out[i] = q; placed.append(q)
        return out
    if variant in ('field', 'softfield'):
        waves = [[(rng.uniform(0, 2 * math.pi), 2 * math.pi / rng.uniform(s / 4, s), rng.uniform(0, 2 * math.pi), rng.uniform(0.5, 1))
                  for _ in range(3)] for _ in range(P)]
        def f(pidx, x, y):
            return sum(w * math.cos(kk * (x * math.cos(th) + y * math.sin(th)) + ph) for th, kk, ph, w in waves[pidx]) / 3
        out = []
        for i in range(n):
            x, y = parents[0][i][0], parents[0][i][1]
            v = [f(pidx, x, y) for pidx in range(P)]
            if variant == 'field':
                pidx = max(range(P), key=lambda q: v[q])
            else:
                w = [math.exp(4 * z) for z in v]; pidx = rng.choices(range(P), weights=w)[0]
            out.append(parents[pidx][i])
        return out
    if variant == 'half':
        th = rng.uniform(0, 2 * math.pi); c = (s / 2 + rng.uniform(-s / 4, s / 4), s / 2 + rng.uniform(-s / 4, s / 4))
        other = rng.randrange(1, P)
        return [parents[0][i] if (parents[0][i][0] - c[0]) * math.cos(th) + (parents[0][i][1] - c[1]) * math.sin(th) < 0
                else parents[other][i] for i in range(n)]
    raise ValueError(variant)


def job(j):
    variant, rpath, others, P, k, smax, seed = j
    from layout import full_lines
    rng = random.Random(seed)
    t0 = time.time()
    sR, R = mcmin.load_deg(rpath)
    rec = dict(variant=variant, ref=rpath)
    if variant == 'kick':
        prop, _ = explore.MOVES['kick'](sR, R, rng, explore.A)
        rec['parents'] = [sR]
    else:
        cand = []
        for p in others:
            sQ, Q = mcmin.load_deg(p)
            d, Qa = align(R, sR, Q, sQ)
            cand.append((d, sQ, Qa))
        cand.sort(key=lambda c: c[0])
        parents = [R] + [c[2] for c in cand[:P - 1]]
        rec['parents'] = [sR] + [c[1] for c in cand[:P - 1]]
        rec['dist'] = [round(c[0], 4) for c in cand[:P - 1]]
        prop = recombine(variant, sR, parents, rng)
    tmp = tempfile.mkdtemp()
    hop.EXTRA[:] = ['--loosen', '1.02']
    r = hop.quench(sR, prop, tmp, extra=('--pit', '8', '--flip-top', '0'))
    if r is None:
        return dict(rec, status='fail', sec=time.time() - t0)
    s1, sq1, _ = r
    rec['s_screen'] = s1
    if s1 > smax + 1e-3:
        return dict(rec, status='discard', sec=time.time() - t0)
    p1 = os.path.join(tmp, 's1.txt'); mcmin.write_deg(p1, s1, sq1)
    if full_lines(p1, k):
        return dict(rec, status='grid', sec=time.time() - t0)
    hop.EXTRA[:] = ['--loosen', '1.0']
    r2 = hop.quench(s1, sq1, tmp, extra=('--no-alm',))
    if r2 is None:
        return dict(rec, status='fail', sec=time.time() - t0)
    s2, sq2, _ = r2
    p2 = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p2, s2, sq2)
    if full_lines(p2, k):
        return dict(rec, status='grid', s=s2, sec=time.time() - t0)
    rec.update(status='ok', s=s2, sec=time.time() - t0)
    if s2 < k:
        rec['sq'] = sq2
    return rec


def seen_sides(n):
    S = set()
    if os.path.exists(f'{HERE}/runs/known{n}_all.json'):
        S |= set(json.load(open(f'{HERE}/runs/known{n}_all.json')))
    for f in glob.glob(f'{HERE}/runs/ex*/archive.jsonl'):
        for l in open(f):
            e = json.loads(l)
            S.add(e['s'])
    return sorted(S)


def near(L, s, tol=2e-9):
    i = bisect.bisect_left(L, s - tol)
    return i < len(L) and abs(L[i] - s) < tol


def report(out, n=110):
    R = [json.loads(l) for l in open(f'{out}/res.jsonl')]
    k = math.ceil(math.sqrt(n) - 1e-12)
    seen = seen_sides(n)
    by = collections.defaultdict(list)
    for r in R:
        by[r['variant']].append(r)
    print(f'== {out}: {len(R)} proposals; seen-before reference {len(seen)} sides')
    print(f'  {"variant":12} {"N":>4} {"grid":>5} {"disc":>5} {"ok":>4} {"sub-k":>5} {"=parent":>7} {"new":>4} {"distinct new":>12} {"new/CPU-h":>9} {"deep new":>8} {"best":>14} {"med s_scr":>10}')
    for v, L in by.items():
        st = collections.Counter(r['status'] for r in L)
        sub = [r for r in L if r['status'] == 'ok' and r['s'] < k]
        par = [r for r in sub if any(abs(r['s'] - p) < 2e-9 for p in r['parents'])]
        new = [r for r in sub if not near(seen, r['s'])]
        dn = set(round(r['s'], 9) for r in new)
        cpu = sum(r['sec'] for r in L) / 3600
        ss = sorted(r.get('s_screen', 99) for r in L)
        print(f'  {v:12} {len(L):4d} {st["grid"]:5d} {st["discard"]:5d} {st["ok"]:4d} {len(sub):5d} {len(par):7d} {len(new):4d} '
              f'{len(dn):12d} {len(dn) / cpu:9.1f} {len(set(round(r["s"], 9) for r in new if r["s"] < k - 0.0025)):8d} {min([r["s"] for r in sub] or [99]):14.10f} {ss[len(ss) // 2]:10.5f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive'); ap.add_argument('--n', type=int, default=110); ap.add_argument('--trials', type=int, default=150)
    ap.add_argument('--P', type=int, default=5); ap.add_argument('--cands', type=int, default=30)
    ap.add_argument('--variants', nargs='+', default=['kick', 'uniform', 'bestfit:3', 'field', 'softfield', 'half'])
    ap.add_argument('--procs', type=int, default=8); ap.add_argument('--out'); ap.add_argument('--report')
    ap.add_argument('--seed', type=int, default=3)
    a = ap.parse_args()
    if a.report:
        report(a.report, a.n); raise SystemExit
    k = math.ceil(math.sqrt(a.n) - 1e-12); smax = k + 0.05
    E = [json.loads(l) for l in open(f'{a.archive}/archive.jsonl')]
    pool = [e['path'] for e in E if e['s'] < k and os.path.exists(e['path'])]
    rng = random.Random(a.seed)
    jobs = []
    for t in range(a.trials):                    # same reference + candidate sample for every variant (paired)
        ref = rng.choice(pool)
        others = rng.sample([p for p in pool if p != ref], min(a.cands, len(pool) - 1))
        seed = rng.randrange(1 << 30)
        for v in a.variants:
            jobs.append((v, ref, others, a.P, k, smax, seed))
    os.makedirs(a.out, exist_ok=True)
    res = run_jobs(job, jobs, procs=a.procs, timeout=900)
    with open(f'{a.out}/res.jsonl', 'w') as f:
        for r in res:
            if isinstance(r, dict):
                f.write(json.dumps(r) + '\n')
    report(a.out, a.n)
