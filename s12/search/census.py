#!/usr/bin/env python3
"""census.py -- the integral margin of EVERY leaf of the full pattern tree on Bentz's P0, t = 4.

Why (2026-09-19).  `search/RANK8.md` measured one leaf: leaf A = `(AB)^4` has margin exactly 0 on
a large plateau.  The full pattern tree (`bentz.py tree`) has 86,403 leaves, 10,945 up to D4: a
leaf is a family of twelve pairwise-disjoint realisable patterns, one labelled square per pattern,
each square confined to its pattern region.  Every degree-1 value on these leaves that matters is
`>= 12`, so the LP says nothing; what decides the shape of a proof is the INTEGRAL margin

    margin(leaf) = sup over the (closed relaxation of the) leaf of  min_{i<j} gap(S_i, S_j),

    < 0 : the leaf is strictly infeasible -- a compact statement, closable by exhaustive search
    = 0 : a plateau leaf -- dies only by touching, needs a strictness (chain) argument
    > 0 : a genuine closed packing of twelve unit squares in [0,4]^2, i.e. s(12) < 4.

This file measures it with `rank8.py`'s instrument (assignment-fixed SLSQP multistart), which is
hard-wired to leaf A only through the module global `rank8.LABELS`; a worker swaps that global for
the family at hand.  **Floating point, multistart, local: a value is a LOWER bound on the leaf's
margin, and a negative value is a measurement, never a certificate** (`RANK8.md` §2 records a false
negative of exactly this kind at 8 starts).  The one exact-in-spirit ingredient is `tiling`: the
combinatorial list of leaves whose closed relaxation contains a 4x4-tiling-minus-4 configuration;
those have margin `>= 0` by an explicit configuration with pairwise disjoint interiors.

Sub-commands
------------
  families : enumerate the leaves up to D4, mark the tiling-compatible ones, write runs/census_families.json
  run      : pass 1 -- multistart on every class, streaming to a JSONL (resumable)
  refine   : pass 2 -- many more starts + kicks on the classes a pass left in a band
  report   : the table
"""
import argparse
import itertools
import json
import math
import os
import sys
import time

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bentz                                                         # noqa: E402
import rank8                                                         # noqa: E402

T = 4.0
PTS = bentz.bentz_points()
PNAMES = [n for (n, _x, _y) in PTS]
PXY = rank8.PXY
PERMS = bentz.d4_perm(PTS)
# angle action of the eight maps in bentz.d4_maps order: e mx my r2 d r1 r3 ad
REFLECT = [False, True, True, False, True, False, False, True]
MAPS = [g for (_n, g) in bentz.d4_maps(T)]


def bits(b):
    return tuple(i for i in range(16) if b >> i & 1)


def fam_name(f):
    return ' '.join('{' + ','.join(PNAMES[i] for i in bits(b)) + '}' for b in f)


def img(b, pm):
    return sum(1 << pm[i] for i in range(16) if b >> i & 1)


def canon_family(f):
    """-> (canonical tuple of bitmasks, index of a D4 map carrying f to it)"""
    best, arg = None, None
    for k, pm in enumerate(PERMS):
        key = tuple(sorted(img(b, pm) for b in f))
        if best is None or key < best:
            best, arg = key, k
    return best, arg


def map_config(poses, k):
    """apply D4 map k to a list of (x, y, th)"""
    out = []
    for (x, y, th) in poses:
        X, Y = MAPS[k](x, y)
        out.append((X, Y, -th if REFLECT[k] else th))
    return out


def float_pattern(x, y, th, tol=0.0):
    c, s = math.cos(th), math.sin(th)
    dx, dy = PXY[:, 0] - x, PXY[:, 1] - y
    u, v = np.abs(c * dx + s * dy), np.abs(-s * dx + c * dy)
    m = np.maximum(u, v)
    return m


# ===================================================================== families
def enumerate_families(pitch, dth):
    ok, _drop, _n = bentz.realizable(pitch, dth)
    pat = sorted(b for b in ok if b)
    reps = []

    def rec(start, used, k, cur):
        if k == 12:
            reps.append(tuple(cur))
            return
        if 16 - bin(used).count('1') < 12 - k:
            return
        for idx in range(start, len(pat)):
            b = pat[idx]
            if b & used:
                continue
            cur.append(b)
            rec(idx + 1, used | b, k + 1, cur)
            cur.pop()

    rec(0, 0, 0, [])
    classes = {}
    for f in reps:
        c, _k = canon_family(f)
        classes[c] = classes.get(c, 0) + 1
    wit = {b: (float(x), float(y), 2 * math.atan2(p, q)) for b, (p, q, x, y) in ok.items() if b}
    return pat, wit, reps, classes


def tiling_compatible(patset):
    """every (tiling-minus-4, assignment) -> canonical class -> list of configs in canonical frame.

    A tile's relaxed pattern may be any set  mandatory <= pi <= mandatory + boundary  (P0 points
    in its open interior are mandatory, those on its boundary optional: containment is closed and
    the relaxed exclusion is `not in the open square`)."""
    tiles = [(i, j) for i in range(4) for j in range(4)]
    man, opt = {}, {}
    for (i, j) in tiles:
        m = float_pattern(i + 0.5, j + 0.5, 0.0)
        man[(i, j)] = sum(1 << p for p in range(16) if m[p] < 0.5 - 1e-9)
        opt[(i, j)] = [p for p in range(16) if abs(m[p] - 0.5) <= 1e-9]
    choices = {}
    for t in tiles:
        ch = []
        for r in range(len(opt[t]) + 1):
            for sub in itertools.combinations(opt[t], r):
                b = man[t] | sum(1 << p for p in sub)
                if b and b in patset:
                    ch.append(b)
        choices[t] = ch
    out = {}
    nconf = 0
    for removed in itertools.combinations(tiles, 4):
        keep = [t for t in tiles if t not in removed]

        def rec(k, used, cur):
            nonlocal nconf
            if k == 12:
                nconf += 1
                f = tuple(cur)
                c, g = canon_family(f)
                poses = map_config([(t[0] + 0.5, t[1] + 0.5, 0.0) for t in keep], g)
                # order the poses to match the canonical tuple
                imgs = [img(b, PERMS[g]) for b in f]
                order = [imgs.index(b) for b in c]
                lst = out.setdefault(c, [])
                if len(lst) < 6:
                    lst.append([poses[o] for o in order])
                return
            for b in choices[keep[k]]:
                if b & used:
                    continue
                cur.append(b)
                rec(k + 1, used | b, cur)
                cur.pop()
        rec(0, 0, [])
    return out, nconf


def cmd_families(a):
    t0 = time.time()
    pat, wit, reps, classes = enumerate_families(a.pitch, a.dth)
    print(f'{len(pat)} realisable non-empty patterns; {len(reps)} leaves, {len(classes)} up to D4 '
          f'({time.time() - t0:.0f} s)')
    tc, nconf = tiling_compatible(set(pat))
    missing = [c for c in tc if c not in classes]
    print(f'{nconf} (tiling-minus-4, assignment) pairs land in {len(tc)} classes '
          f'({len(missing)} of them not in the enumerated tree -- should be 0)')
    from collections import Counter
    K = Counter(sum(bin(b).count('1') - 1 for b in c) for c in classes)
    Kt = Counter(sum(bin(b).count('1') - 1 for b in c) for c in tc)
    print('  classes by K = sum(|pi|-1):           ' + str(dict(sorted(K.items()))))
    print('  tiling-compatible classes by K:       ' + str(dict(sorted(Kt.items()))))
    leafA = canon_family(tuple(sum(1 << p for p in s) for (_n, s) in rank8.LABELS))[0]
    print(f'  leaf A is class #{sorted(classes).index(leafA)}, tiling-compatible: {leafA in tc}')
    data = dict(patterns=pat, witness={str(b): wit[b] for b in pat},
                classes=[dict(id=k, fam=list(c), orbit=classes[c],
                              tiling=[[list(p) for p in cfg] for cfg in tc.get(c, [])])
                         for k, c in enumerate(sorted(classes))],
                leafA=sorted(classes).index(leafA))
    with open(a.out, 'w') as f:
        json.dump(data, f)
    print(f'wrote {a.out}')


# ===================================================================== pools and starts
def sample_pools(patset, nsamp, rng):
    pools = {b: [] for b in patset}
    B = 200000
    w16 = (1 << np.arange(16))
    done = 0
    while done < nsamp:
        b = min(B, nsamp - done)
        done += b
        th = rng.uniform(0, math.pi / 2, b)
        # half the samples nearly axis-parallel: that is where the plateau lives
        th[: b // 2] = rng.normal(0, 0.03, b // 2)
        w = 0.5 * (np.abs(np.cos(th)) + np.abs(np.sin(th)))
        x = rng.uniform(0, 1, b) * (T - 2 * w) + w
        y = rng.uniform(0, 1, b) * (T - 2 * w) + w
        c, s = np.cos(th), np.sin(th)
        dx = PXY[None, :, 0] - x[:, None]
        dy = PXY[None, :, 1] - y[:, None]
        inside = np.maximum(np.abs(c[:, None] * dx + s[:, None] * dy),
                            np.abs(-s[:, None] * dx + c[:, None] * dy)) <= 0.5
        key = (inside * w16[None, :]).sum(1)
        for r in range(b):
            k = int(key[r])
            p = pools.get(k)
            if p is not None and len(p) < 4000:
                p.append((float(x[r]), float(y[r]), float(th[r])))
    return pools


_W = {}


def _winit(pools, wit, opts):
    _W.update(pools=pools, wit=wit, opts=opts)


def _set_family(fam):
    rank8.LABELS = [(f'S{k}', bits(b)) for k, b in enumerate(fam)]
    rank8.LNAME = [n for (n, _s) in rank8.LABELS]


def _starts(cls, rng, nrand, nnudge, sigma, sigma_th):
    fam = cls['fam']
    pools, wit = _W['pools'], _W['wit']
    out = []
    for cfg in cls['tiling']:
        out.append(('tiling', rank8.z_from([tuple(p) for p in cfg])))
    for k in range(nnudge if cls['tiling'] else 0):
        cfg = cls['tiling'][rng.integers(len(cls['tiling']))]
        z = rank8.z_from([tuple(p) for p in cfg])
        sc = (1.0, 0.25, 0.05)[k % 3]
        z[1::3] += rng.normal(0, sigma * sc, 12)
        z[2::3] += rng.normal(0, sigma * sc, 12)
        z[3::3] += rng.normal(0, sigma_th * sc, 12)
        out.append(('nudge', z))
    for k in range(nrand):
        poses = []
        for b in fam:
            p = pools.get(b) or []
            poses.append(p[rng.integers(len(p))] if p else tuple(wit[b]))
        out.append(('rand', rank8.z_from(poses)))
    return out


def _clip(z):
    """push every square back inside the container (a start only)"""
    th = z[3::3]
    w = 0.5 * (np.abs(np.cos(th)) + np.abs(np.sin(th)))
    z[1::3] = np.minimum(np.maximum(z[1::3], w), T - w)
    z[2::3] = np.minimum(np.maximum(z[2::3], w), T - w)
    return z


def _solve_many(cls, starts, o):
    labels = list(range(12))
    res = []
    for (tag, z0) in starts:
        try:
            val, z, _pr = rank8.solve(labels, _clip(z0), eps_excl=0.0, rounds=o['rounds'],
                                      maxiter=o['maxiter'])
        except Exception as e:                       # noqa: BLE001
            val, z = -float('inf'), z0
            tag = f'{tag}!{type(e).__name__}'
        res.append((val, tag, z))
    return res


def _class_job(cls):
    o = _W['opts']
    t0 = time.time()
    _set_family(cls['fam'])
    rng = np.random.default_rng(o['seed'] + 7919 * cls['id'])
    starts = _starts(cls, rng, o['nrand'], o['nnudge'], o['sigma'], o['sigma_th'])
    seeds = [tuple(s) for s in cls.get('seeds', [])]
    for s in seeds:
        starts.append(('seed', np.array(s)))
    res = _solve_many(cls, starts, o)
    # kicks around the incumbents (basin hopping, RANK8.md 2: local optimisation alone lies)
    for _round in range(o['bh']):
        res.sort(key=lambda r: -r[0])
        inc = [r for r in res if r[0] > -float('inf')][:o['bhkeep']]
        kicks = []
        for k in range(o['bhn']):
            if not inc:
                break
            z = inc[rng.integers(len(inc))][2].copy()
            sc = (1.0, 0.3, 0.1)[k % 3]
            z[1::3] += rng.normal(0, o['sigma'] * sc, 12)
            z[2::3] += rng.normal(0, o['sigma'] * sc, 12)
            z[3::3] += rng.normal(0, o['sigma_th'] * sc, 12)
            kicks.append(('kick', z))
        res += _solve_many(cls, kicks, o)
    res.sort(key=lambda r: -r[0])
    fin = [r for r in res if r[0] > -float('inf')]
    best = res[0]
    z = best[2]
    inc_s, exc_s = rank8.pattern_slack(z, list(range(12)))
    return dict(id=cls['id'], best=best[0], tag=best[1], nstart=len(res), nfin=len(fin),
                tilt=rank8.max_tilt_deg(z), exc=exc_s,
                top=[r[0] for r in fin[:5]],
                z=[float(v) for v in z], sec=time.time() - t0)


def _pair_job(arg):
    b1, b2 = arg
    o = _W['opts']
    pools, wit = _W['pools'], _W['wit']
    rank8.LABELS = [('S0', bits(b1)), ('S1', bits(b2))]
    rank8.LNAME = ['S0', 'S1']
    rng = np.random.default_rng(o['seed'] + 65537 * b1 + b2)
    best, bz = -float('inf'), None
    for k in range(o['nrand']):
        poses = []
        for b in (b1, b2):
            pl = pools.get(b) or []
            poses.append(pl[rng.integers(len(pl))] if pl and k else tuple(wit[b]))
        try:
            val, z, _pr = rank8.solve([0, 1], _clip(rank8.z_from(poses)), rounds=o['rounds'],
                                      maxiter=o['maxiter'])
        except Exception:                            # noqa: BLE001
            continue
        if val > best:
            best, bz = val, z
        if best > o['enough']:
            break
    # a pair still <= 0 gets the long treatment: kicks around the incumbent and around both
    # witnesses (tiny pattern regions have tiny pools, and the first pass stalls on them)
    if best <= 1e-9 and bz is not None:
        for k in range(o['nkick']):
            z = bz.copy() if k % 2 else _clip(rank8.z_from([tuple(wit[b1]), tuple(wit[b2])]))
            sc = (0.3, 0.1, 0.03, 0.01)[k % 4]
            z[1::3] += rng.normal(0, sc, 2)
            z[2::3] += rng.normal(0, sc, 2)
            z[3::3] += rng.normal(0, sc, 2)
            try:
                val, z, _pr = rank8.solve([0, 1], _clip(z), rounds=o['rounds'],
                                          maxiter=o['maxiter'])
            except Exception:                        # noqa: BLE001
                continue
            if val > best:
                best, bz = val, z
    return (b1, b2, best, [float(v) for v in bz] if bz is not None else None)


def cmd_pairs(a):
    """the margin of every PAIR of disjoint patterns: two squares, six unknowns, where multistart
    is reliable.  A leaf's margin is at most the smallest of its 66 pair margins, so this table
    (i) is an upper bound for every class and (ii) says when a class value is a stall."""
    import multiprocessing as mp
    data = json.load(open(a.families))
    pat = data['patterns']
    rng = np.random.default_rng(a.seed)
    pools = sample_pools(set(pat), a.nsamp, rng)
    wit = {int(b): v for b, v in data['witness'].items()}
    opts = dict(seed=a.seed, nrand=a.nrand, rounds=8, maxiter=300, enough=a.enough, nkick=a.nkick)
    jobs = [(b1, b2) for i, b1 in enumerate(pat) for b2 in pat[i + 1:] if not b1 & b2]
    print(f'# {len(jobs)} disjoint pattern pairs', flush=True)
    t0 = time.time()
    with mp.Pool(a.nproc, initializer=_winit, initargs=(pools, wit, opts)) as pool:
        res = pool.map(_pair_job, jobs, chunksize=8)
    neg = [r for r in res if r[2] <= 1e-9]
    print(f'{len(res)} pairs in {time.time() - t0:.0f} s; {len(neg)} with margin <= 0:')
    for (b1, b2, v, _z) in sorted(neg, key=lambda r: r[2]):
        print(f'  {v:+.6e}  {fam_name([b1])} | {fam_name([b2])}')
    json.dump([dict(a=b1, b=b2, v=v, z=z) for (b1, b2, v, z) in res], open(a.out, 'w'))
    print(f'wrote {a.out}')


def _load_done(path):
    done = {}
    if os.path.exists(path):
        for line in open(path):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r['id'] not in done or r['best'] > done[r['id']]['best']:
                done[r['id']] = r
    return done


def _run(a, todo, data, outpath):
    import multiprocessing as mp
    rng = np.random.default_rng(a.seed)
    t0 = time.time()
    pools = sample_pools(set(data['patterns']), a.nsamp, rng)
    small = sorted((len(v), b) for b, v in pools.items())[:8]
    print(f'# pools from {a.nsamp} samples in {time.time() - t0:.0f} s; smallest: '
          + ', '.join(f'{fam_name([b])}:{n}' for n, b in small), flush=True)
    wit = {int(b): v for b, v in data['witness'].items()}
    opts = dict(seed=a.seed, nrand=a.nrand, nnudge=a.nnudge, sigma=a.sigma, sigma_th=a.sigma_th,
                rounds=a.rounds, maxiter=a.maxiter, bh=a.bh, bhn=a.bhn, bhkeep=a.bhkeep)
    print(f'# {len(todo)} classes to do, options {opts}', flush=True)
    n = 0
    with mp.Pool(a.nproc, initializer=_winit, initargs=(pools, wit, opts)) as pool, \
            open(outpath, 'a') as f:
        for r in pool.imap_unordered(_class_job, todo, chunksize=1):
            f.write(json.dumps(r) + '\n')
            f.flush()
            n += 1
            if r['best'] > 1e-9:
                print(f'!!! POSITIVE  class {r["id"]}  {r["best"]:+.6e}', flush=True)
            if n % 100 == 0 or n == len(todo):
                el = time.time() - t0
                print(f'  {n}/{len(todo)}  {el:.0f} s  eta {el / n * (len(todo) - n):.0f} s',
                      flush=True)


def cmd_run(a):
    data = json.load(open(a.families))
    done = _load_done(a.out)
    todo = [c for c in data['classes'] if c['id'] not in done]
    if a.only:
        ids = set(int(x) for x in a.only.split(','))
        todo = [c for c in data['classes'] if c['id'] in ids]
    if a.limit:
        todo = todo[:a.limit]
    _run(a, todo, data, a.out)


def cmd_refine(a):
    """pass 2: classes whose best value so far lies in [lo, hi], with many more starts, seeded
    from the incumbent"""
    data = json.load(open(a.families))
    done = _load_done(a.inp)
    for p in a.also or []:
        for k, r in _load_done(p).items():
            if k not in done or r['best'] > done[k]['best']:
                done[k] = r
    todo = []
    for c in data['classes']:
        r = done.get(c['id'])
        if r is None or not (a.lo <= r['best'] <= a.hi):
            continue
        if a.skip_tiling and c['tiling']:
            continue
        c = dict(c)
        c['seeds'] = [r['z']]
        todo.append(c)
    already = _load_done(a.out)
    todo = [c for c in todo if c['id'] not in already]
    _run(a, todo, data, a.out)


def cmd_report(a):
    data = json.load(open(a.families))
    done = {}
    for p in a.FILES:
        for k, r in _load_done(p).items():
            if k not in done or r['best'] > done[k]['best']:
                done[k] = r
    cl = {c['id']: c for c in data['classes']}
    print(f'{len(done)} of {len(cl)} classes measured')
    vals = np.array([r['best'] for r in done.values()])
    z = a.zero
    pos = [r for r in done.values() if r['best'] > z]
    zer = [r for r in done.values() if -z <= r['best'] <= z]
    neg = [r for r in done.values() if r['best'] < -z]
    ntile = sum(1 for c in cl.values() if c['tiling'])
    print(f'tiling-compatible classes (margin >= 0 by an explicit configuration): {ntile}')
    print(f'margin > +{z:g}: {len(pos)}    |margin| <= {z:g}: {len(zer)}    margin < -{z:g}: {len(neg)}')
    zt = sum(1 for r in zer if cl[r['id']]['tiling'])
    print(f'  of the zero classes, {zt} are tiling-compatible and {len(zer) - zt} are not')
    tn = [r for r in neg if cl[r['id']]['tiling']]
    print(f'  tiling-compatible classes the optimiser left NEGATIVE (instrument failures): {len(tn)}')
    edges = [-10, -0.5, -0.2, -0.1, -0.05, -0.02, -0.01, -0.005, -0.002, -0.001, -1e-4, -1e-5, -z]
    print('histogram of the negative margins:')
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = int(((vals >= lo) & (vals < hi)).sum())
        print(f'  [{lo:+.0e}, {hi:+.0e})  {k}')
    for title, rows in (('POSITIVE', sorted(pos, key=lambda r: -r['best'])[:a.n]),
                        ('zero, not tiling-compatible', sorted(
                            [r for r in zer if not cl[r['id']]['tiling']],
                            key=lambda r: -r['tilt'])[:a.n]),
                        ('closest to zero from below', sorted(neg, key=lambda r: -r['best'])[:a.n])):
        print(f'\n{title}:')
        for r in rows:
            c = cl[r['id']]
            print(f'  #{r["id"]:<6} {r["best"]:+.6e}  tilt {r["tilt"]:7.3f}  orbit {c["orbit"]}  '
                  f'tiling {len(c["tiling"]) > 0!s:<5}  {fam_name(c["fam"])}')
    if a.by_k:
        from collections import defaultdict
        t = defaultdict(lambda: [0, 0, 0])
        for r in done.values():
            K = sum(bin(b).count('1') - 1 for b in cl[r['id']]['fam'])
            t[K][0 if r['best'] > z else 1 if r['best'] >= -z else 2] += 1
        print('\nby K = sum(|pi|-1) (= 4 - u):   K: positive / zero / negative')
        for K in sorted(t):
            print(f'  {K}: {t[K][0]} / {t[K][1]} / {t[K][2]}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('families')
    c.add_argument('--pitch', type=float, default=0.01)
    c.add_argument('--dth', type=float, default=0.5)
    c.add_argument('--out', default='runs/census_families.json')
    for name in ('run', 'refine'):
        c = sub.add_parser(name)
        c.add_argument('--families', default='runs/census_families.json')
        c.add_argument('--out', default=f'runs/census_{name}.jsonl')
        c.add_argument('--seed', type=int, default=20260919)
        c.add_argument('--nsamp', type=int, default=6000000)
        c.add_argument('--nrand', type=int, default=24)
        c.add_argument('--nnudge', type=int, default=12)
        c.add_argument('--sigma', type=float, default=0.08)
        c.add_argument('--sigma-th', type=float, default=0.12)
        c.add_argument('--rounds', type=int, default=8)
        c.add_argument('--maxiter', type=int, default=300)
        c.add_argument('--bh', type=int, default=1)
        c.add_argument('--bhn', type=int, default=12)
        c.add_argument('--bhkeep', type=int, default=4)
        c.add_argument('--nproc', type=int, default=12)
        if name == 'run':
            c.add_argument('--only', default=None)
            c.add_argument('--limit', type=int, default=0)
        else:
            c.add_argument('--inp', default='runs/census_run.jsonl')
            c.add_argument('--also', nargs='*')
            c.add_argument('--lo', type=float, default=-0.02)
            c.add_argument('--hi', type=float, default=-1e-7)
            c.add_argument('--skip-tiling', action='store_true')
    c = sub.add_parser('pairs')
    c.add_argument('--families', default='runs/census_families.json')
    c.add_argument('--out', default='runs/census_pairs.json')
    c.add_argument('--seed', type=int, default=20260919)
    c.add_argument('--nsamp', type=int, default=6000000)
    c.add_argument('--nrand', type=int, default=60)
    c.add_argument('--enough', type=float, default=0.05)
    c.add_argument('--nkick', type=int, default=400)
    c.add_argument('--nproc', type=int, default=12)
    c = sub.add_parser('report')
    c.add_argument('--pairs', default='runs/census_pairs.json')
    c.add_argument('FILES', nargs='+')
    c.add_argument('--families', default='runs/census_families.json')
    c.add_argument('--zero', type=float, default=1e-7)
    c.add_argument('-n', type=int, default=25)
    c.add_argument('--by-k', action='store_true')
    a = ap.parse_args()
    dict(families=cmd_families, pairs=cmd_pairs, run=cmd_run, refine=cmd_refine, report=cmd_report)[a.cmd](a)


if __name__ == '__main__':
    main()
