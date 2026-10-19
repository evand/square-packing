#!/usr/bin/env python3
"""Rational-angle block enumerator (10-18; Evan: "a family of exact constructions, compare side lengths, look for records
or near records, then optimize with those as start points").

Design = Pythagorean angle (a, b, c) + block of C rows x L squares at that angle (unsheared at the start) + placement
(fx, fy: fractions of the free range) + loose start side K + f.  Pipeline: block in the box -> fillfast exact upright
fill -> fixed-angle LP polish (rational/fixlp.py) -> refill at the new side (more squares?) -> polish, until the count
stops growing.  Output per design: (N, s*), compared with the frontier at N.

  ratgen.py run --out DIR [--K 10-15] [--triples 3,4,5:20,21,29] [--procs 14]
  ratgen.py report DIR
"""
import sys, os, math, json, argparse, itertools, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import fillfast, fixlp

TRIPLES = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (12, 35, 37), (9, 40, 41), (28, 45, 53),
           (33, 56, 65), (16, 63, 65), (48, 55, 73), (36, 77, 85), (39, 80, 89), (65, 72, 97), (119, 120, 169)]


def block(s, th, C, L, fx, fy, shear=None):
    """C rows of L squares at angle th; row j shifted by shear[j] along the row direction (sheared block)"""
    t = math.radians(th); u = np.array([math.cos(t), math.sin(t)]); v = np.array([-math.sin(t), math.cos(t)])
    sh = list(shear) if shear is not None else [0.0] * C
    P = np.array([(i + sh[j] - (L - 1) / 2) * u + (j - (C - 1) / 2) * v for j in range(C) for i in range(L)])
    e = 0.5 * (abs(u[0]) + abs(u[1]))
    lo = P.min(0) - e; hi = P.max(0) + e
    free = s - (hi - lo)
    if (free < 0).any():
        return None
    off = np.array([fx, fy]) * free - lo
    P = P + off
    return np.c_[P, np.full(len(P), th % 90)]


def nontilted(sq):
    return np.abs(((sq[:, 2] + 45) % 90) - 45) < 1e-6


def ladder(sq, s, steps=12, drop=1.2):
    """removal ladder: remove the most-stressed square (largest LP contact force), polish, refill, polish; record the
    best (s, sq) per n on the way.  -> {n: (s, sq)}"""
    out = {len(sq): (s, sq)}
    s_top = s
    for _ in range(steps):
        stv = getattr(fixlp.polish, 'stress', None)
        if stv is None or len(stv) != len(sq) or stv.max() <= 0:
            break
        sq = np.delete(sq, int(stv.argmax()), 0)
        sq, s = fixlp.polish(sq, s)
        for _r in range(3):
            sq2 = fillfast.fill(s + 1e-9, sq[~nontilted(sq)])
            if len(sq2) <= len(sq):
                break
            sq, s = fixlp.polish(sq2, s)
        if fixlp.max_pen(sq, s) > 1e-9:
            break
        n = len(sq)
        if n not in out or s < out[n][0] - 1e-12:
            out[n] = (s, sq)
        if s < s_top - drop:
            break
    return out


def design(args, steps=0):
    th, C, L, fx, fy, s0 = args[:6]
    shear = args[6] if len(args) > 6 else None
    T = block(s0, th, C, L, fx, fy, shear)
    if T is None:
        return None
    sq = fillfast.fill(s0, T); s = s0
    hist = []
    for rnd in range(6):
        sq, s = fixlp.polish(sq, s)
        hist.append((len(sq), s))
        T2 = sq[~nontilted(sq)]
        sq2 = fillfast.fill(s + 1e-9, T2)
        if len(sq2) <= len(sq):
            break
        sq = sq2
    pen = fixlp.max_pen(sq, s)
    lad = {}
    if steps and pen <= 1e-9:
        lad = {int(k): (float(v[0]), [list(map(float, r)) for r in v[1]]) for k, v in ladder(sq, s, steps).items()}
    return dict(th=th, C=C, L=L, fx=fx, fy=fy, s0=s0, shear=shear, ladder=lad, n=len(sq), s=s, pen=pen, hist=hist,
                sq=[list(map(float, r)) for r in sq])


def tasks(Ks, triples, fracs=(0.3, 0.55, 0.8, 0.97), grid=5):
    for K in Ks:
        for (a, b, c) in triples:
            for th in (math.degrees(math.atan2(a, b)), math.degrees(math.atan2(b, a))):
                for C in range(1, 7):
                    for L in range(C, int(1.45 * (K + 1)) + 1):
                        for f in fracs:
                            s0 = K + f
                            for fx, fy in itertools.product(np.linspace(0, 1, grid), repeat=2):
                                yield (th, C, L, float(fx), float(fy), s0)


def frontier():
    from pk.store import Store
    st = Store()
    out = {n: v[0] for n, v in st.public_best().items()}
    for n, r in st.claims().items():
        out[n] = min(out.get(n, 1e9), r['s'])
    return out


def run(a):
    from multiprocessing import Pool
    os.makedirs(a.out, exist_ok=True)
    Ks = range(int(a.K.split('-')[0]), int(a.K.split('-')[1]) + 1)
    triples = [tuple(map(int, t.split(','))) for t in a.triples.split(':')] if a.triples else TRIPLES
    T = list(tasks(Ks, triples, grid=a.grid))
    print(len(T), 'designs', flush=True)
    best = {}
    t0 = time.time()
    with Pool(a.procs) as pool, open(os.path.join(a.out, 'designs.jsonl'), 'w') as fo:
        for k, r in enumerate(pool.imap_unordered(design, T, chunksize=8)):
            if r is None or r['pen'] > 1e-9:
                continue
            key = r['n']
            if key not in best or r['s'] < best[key]['s'] - 1e-12:
                best[key] = r
                fo.write(json.dumps(r) + '\n'); fo.flush()
            if k % 5000 == 0:
                print(k, f'{time.time()-t0:.0f}s', flush=True)
    json.dump(best, open(os.path.join(a.out, 'best.json'), 'w'))


def rows_block(th, rows):
    """block with per-row (length, shift): rows = [(L_j, d_j)], centred at the origin"""
    t = math.radians(th); u = np.array([math.cos(t), math.sin(t)]); v = np.array([-math.sin(t), math.cos(t)])
    C = len(rows)
    return np.array([(i + d) * u + (j - (C - 1) / 2) * v for j, (L, d) in enumerate(rows) for i in range(L)])


def place_multi(s, th, blocks):
    """blocks = [(rows, fx, fy)]: each block placed at fractions of its free range; later blocks' squares that
    overlap earlier ones are dropped (the fill / LP sort out the rest)"""
    t = math.radians(th); e = 0.5 * (abs(math.cos(t)) + abs(math.sin(t)))
    out = []
    for rows, fx, fy in blocks:
        P = rows_block(th, rows)
        lo = P.min(0) - e; hi = P.max(0) + e
        free = s - (hi - lo)
        if (free < 0).any():
            return None
        P = P + np.array([fx, fy]) * free - lo
        if out:
            Q = np.array(out)[:, :2]
            P = P[(np.abs(P[:, None, :] - Q[None]).max(-1) > 2 * e - 1e-9).all(1)]
        out += [(x, y, th % 90) for x, y in P]
    return np.array(out)


def design_multi(args):
    th, blocks, s0, tri = args
    T = place_multi(s0, th, blocks)
    if T is None or len(T) == 0:
        return None
    sq = fillfast.fill(s0, T); s = s0
    for rnd in range(6):
        sq, s = fixlp.polish(sq, s)
        sq2 = fillfast.fill(s + 1e-9, sq[~nontilted(sq)])
        if len(sq2) <= len(sq):
            break
        sq = sq2
    return dict(th=th, blocks=blocks, s0=s0, tri=list(tri), n=len(sq), s=s, pen=fixlp.max_pen(sq, s),
                C=len(blocks[0][0]), L=max(L for L, d in blocks[0][0]), shear=None, fx=blocks[0][1], fy=blocks[0][2],
                sq=[list(map(float, r)) for r in sq])


TRIPLES2 = [(20, 21, 29), (48, 55, 73), (65, 72, 97), (119, 120, 169), (3, 4, 5), (28, 45, 53)]


def tasks2(Ks, triples, samples, Cs=(2, 3, 4), seed=0):
    """one orientation per triple (the placement grid and shear set are mirror-symmetric, so 90 - th is redundant);
    per (K, triple, C, L): `samples` random (row shears on the 1/4 grid in [-1, 1), placement, start side)"""
    rng = np.random.default_rng(seed)
    grid = np.arange(-4, 4) / 4
    for K in Ks:
        for (a, b, c) in triples:
            th = math.degrees(math.atan2(a, b))
            for C in Cs:
                for L in range(max(C, K - 4), int(1.45 * (K + 1)) + 1):
                    for _ in range(samples):
                        sh = (0.0,) + tuple(float(x) for x in rng.choice(grid, C - 1))
                        fx, fy = (float(x) for x in rng.random(2))
                        f = float(rng.choice([0.6, 0.8, 0.97]))
                        yield (th, C, L, fx, fy, K + f, sh, (a, b, c))


def design2(args):
    r = design(args[:7])
    if r is not None:
        r['tri'] = list(args[7])
    return r


def run2(a):
    from multiprocessing import Pool
    import collections
    sys.path.insert(0, '.')
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'ratgen-rg2', seed=a.seed, args=vars(a), procs=a.procs)
    fr = frontier()
    Ks = range(int(a.K.split('-')[0]), int(a.K.split('-')[1]) + 1)
    T = list(tasks2(Ks, TRIPLES2, a.samples, seed=a.seed))
    print(len(T), 'designs', flush=True)
    kept = collections.defaultdict(dict)               # n -> {round(s, 9): count}
    best = {}
    t0 = time.time()
    with Pool(a.procs) as pool, open(os.path.join(a.out, 'keep.jsonl'), 'w') as fo:
        for k, r in enumerate(pool.imap_unordered(design2, T, chunksize=4)):
            if r is None or r['pen'] > 1e-9:
                continue
            n, s = r['n'], r['s']
            if n not in best or s < best[n]['s'] - 1e-12:
                best[n] = r
            f = fr.get(n)
            if f is not None and s < f + a.gap and abs(s - round(s)) > 1e-9:
                key = round(s, 9)
                if kept[n].get(key, 0) < 3 and (key in kept[n] or len(kept[n]) < 10):
                    kept[n][key] = kept[n].get(key, 0) + 1
                    fo.write(json.dumps(r) + '\n'); fo.flush()
            if k % 20000 == 0:
                print(k, f'{time.time()-t0:.0f}s', flush=True)
    json.dump(best, open(os.path.join(a.out, 'best.json'), 'w'))
    man.finish()


def cmd_import(a):
    from fractions import Fraction as F
    from pk.store import Store
    from pk.packing import Packing
    st = Store(); c = dict(new=0, dup=0, infeasible=0)
    fr = frontier()
    for line in open(os.path.join(a.out, 'keep.jsonl')):
        r = json.loads(line)
        if r['s'] > fr.get(r['n'], 1e9) + a.gap:
            continue
        p = Packing(s=r['s'], sq=np.array(r['sq']))
        tri = r.get('tri')
        tags = ['rational'] + ([f"tri:{tri[0]}-{tri[1]}-{tri[2]}"] if tri else [])
        q = F(r['s']).limit_denominator(100000)
        sx = str(q) if abs(float(q) - r['s']) < 1e-12 else None
        pid, new = st.add(p, source='ours:rational', finder='ours', ref=os.path.basename(a.out.rstrip('/')),
                          kind='rational-design', status='raw', tags=tags, s_exact=sx, commit=False)
        c['new' if new else 'dup'] += 1
    st.db.commit()
    print(c)


def report(a):
    from fractions import Fraction as F
    best = json.load(open(os.path.join(a.out, 'best.json')))
    fr = frontier()
    rows = sorted(((int(n), r) for n, r in best.items()), key=lambda x: x[0])
    for n, r in rows:
        f = fr.get(n)
        if f is None or abs(r['s'] - round(r['s'])) < 1e-9 and r['s'] >= f - 1e-9 and abs(f - round(f)) < 1e-9:
            continue
        d = r['s'] - f if f else float('nan')
        flag = '  <<< BELOW' if d < -1e-9 else ('  = tie' if abs(d) < 1e-9 else ('  near' if d < 1e-3 else ''))
        print(f"{n:4d} s*={r['s']:.10f} ({F(r['s']).limit_denominator(5000)}) front={f:.10f} d={d:+.2e} "
              f"th={r['th']:.2f} C={r['C']} L={r['L']} f=({r['fx']:.2f},{r['fy']:.2f}) s0={r['s0']:.2f}{flag}")


if __name__ == '__main__':
    p = argparse.ArgumentParser(); sub = p.add_subparsers(dest='cmd')
    r = sub.add_parser('run'); r.add_argument('--out', required=True); r.add_argument('--K', default='10-15')
    r.add_argument('--triples', default=''); r.add_argument('--procs', type=int, default=14); r.add_argument('--grid', type=int, default=5)
    r2 = sub.add_parser('run2'); r2.add_argument('--out', required=True); r2.add_argument('--K', default='10-17')
    r2.add_argument('--samples', type=int, default=120); r2.add_argument('--procs', type=int, default=15)
    r2.add_argument('--seed', type=int, default=0); r2.add_argument('--gap', type=float, default=0.03)
    q = sub.add_parser('report'); q.add_argument('out')
    im = sub.add_parser('import'); im.add_argument('out'); im.add_argument('--gap', type=float, default=0.03)
    a = p.parse_args()
    {'run': run, 'run2': run2, 'report': lambda a: report(argparse.Namespace(out=a.out)), 'import': cmd_import}[a.cmd](a)
