#!/usr/bin/env python3
"""Breadth-first basin exploration (10-07, Evan: more explore, less exploit).  Quality-diversity over basins (MAP-Elites
style): an archive of distinct polished minima binned by a content-agnostic descriptor; parents are drawn uniformly over
bins (not by quality), and within a bin favouring rarely expanded basins.  Basins up to S_MAX (above the integer k) are
kept as stepping stones; grid-obstructed states (full lines / staggered axis chains) are discarded.

Quench per proposal: fq screen (ALM + <= 8 SLP iterations, no flips, loosen drawn from 1.0/1.02/1.05).  A screened side
within `same` of a known basin counts as a return (visit, no polish); otherwise full polish (fq --no-alm) and, if new,
an archive entry.

Descriptor (search): (ntilted bin of 2, number of tilt-angle classes, below-k flag).  Evaluation only (n = 110): role
counts (L, B, axis, other) as in census_roles, and novelty against runs/known110.json (48 certified sub-11 minima).

  explore.py --n 110 --starts seeds/rec110.txt ... --minutes 60 --procs 15 --out runs/ex0
  explore.py --report runs/ex0
"""
import argparse, collections, json, math, os, random, subprocess, tempfile, time
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
import mcmin, hop

HERE = os.path.dirname(os.path.abspath(__file__))


class A:                                     # attributes mcmin's moves look up
    alpha = 2.0; reinsert = 'clear'; remove = 'uniform'; tau = 0.0


def mv_bigkick(s, sq, rng, a):
    sig = math.exp(rng.uniform(math.log(0.05), math.log(0.2)))
    return [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), t + rng.gauss(0, 20 * sig)) for x, y, t in sq], f'bigkick {sig:.3f}'


def _axis(t, tol=3.0):
    a = t % 90
    return min(a, 90 - a) < tol


def mv_rowslide(s, sq, rng, a):
    """Slide a contiguous run of axis squares along its row (or column) by 0.2-0.6: a coordinated move kicks never make."""
    ax = [i for i, q in enumerate(sq) if _axis(q[2])]
    if not ax:
        return mv_bigkick(s, sq, rng, a)
    i0 = rng.choice(ax)
    horiz = rng.random() < 0.5
    c = 1 if horiz else 0                    # coordinate shared by the row (y for a horizontal row)
    run = [j for j in ax if abs(sq[j][c] - sq[i0][c]) < 0.35]
    run.sort(key=lambda j: sq[j][1 - c])
    # contiguous component containing i0 (gaps < 0.3 between neighbours along the row)
    k = run.index(i0); lo = hi = k
    while lo > 0 and sq[run[lo]][1 - c] - sq[run[lo - 1]][1 - c] < 1.3: lo -= 1
    while hi < len(run) - 1 and sq[run[hi + 1]][1 - c] - sq[run[hi]][1 - c] < 1.3: hi += 1
    seg = set(run[lo:hi + 1])
    if rng.random() < 0.5:                   # sometimes only part of the run (one end stays)
        seg = set(run[lo:k + 1]) if rng.random() < 0.5 else set(run[k:hi + 1])
    d = rng.choice((-1, 1)) * rng.uniform(0.2, 0.6)
    out = []
    for j, (x, y, t) in enumerate(sq):
        if j in seg:
            x, y = (x + d, y) if horiz else (x, y + d)
        out.append((x, y, t))
    return out, f'rowslide {len(seg)} {d:+.2f} {"h" if horiz else "v"}'


def mv_chainshift(s, sq, rng, a):
    """Vacancy -> hole chain shift: remove a square, move a chain of squares one step each toward the vacancy along the
    line to a clearance hole, re-insert the removed square at the hole."""
    hs = mcmin.holes(s, sq)
    if not hs:
        return mv_bigkick(s, sq, rng, a)
    hx, hy, _ = rng.choices(hs[:10], weights=[h[2] for h in hs[:10]])[0]
    near = sorted(range(len(sq)), key=lambda j: math.hypot(sq[j][0] - hx, sq[j][1] - hy))[:12]
    i = rng.choice(near[2:12]) if len(near) > 3 else near[-1]
    px, py = sq[i][0], sq[i][1]
    L = math.hypot(hx - px, hy - py) or 1.0
    ux, uy = (hx - px) / L, (hy - py) / L
    chain = []
    for j, (x, y, t) in enumerate(sq):       # squares near the segment vacancy -> hole, ordered from the vacancy
        if j == i: continue
        u = (x - px) * ux + (y - py) * uy
        w = abs(-(x - px) * uy + (y - py) * ux)
        if 0 < u < L and w < 0.6:
            chain.append((u, j))
    chain.sort()
    out = list(sq)
    prev = (px, py)
    for _, j in chain:                       # each takes its predecessor's place (keeps its own angle)
        x, y, t = sq[j]
        out[j] = (prev[0], prev[1], t)
        prev = (x, y)
    out[i] = (hx, hy, rng.choice((sq[i][2], 0.0)))
    return out, f'chainshift {len(chain)} L={L:.2f}'


def mv_mirror(s, sq, rng, a):
    """Reflect a local cluster across a horizontal, vertical or diagonal line through a random point (flips junction
    handedness; angles reflect accordingly)."""
    p = (rng.uniform(0, s), rng.uniform(0, s)); r = rng.uniform(1.2, 3.0)
    kind = rng.choice(('h', 'v', 'd'))
    out = []
    n = 0
    for x, y, t in sq:
        if math.hypot(x - p[0], y - p[1]) < r:
            n += 1
            dx, dy = x - p[0], y - p[1]
            if kind == 'h': x, y, t = x, p[1] - dy, -t
            elif kind == 'v': x, y, t = p[0] - dx, y, -t
            else: x, y, t = p[0] + dy, p[1] + dx, 90 - t
        out.append((x, y, t % 90))
    return out, f'mirror {kind} {n}'


MOVES = dict(mcmin.MOVES, bigkick=mv_bigkick, rowslide=mv_rowslide, chainshift=mv_chainshift, mirror=mv_mirror)
WEIGHTS = dict(kick=2, kicksym=1, lkick=2, bigkick=2, crot=1, aswap=2, band=1, reinsert=0.5, rowslide=2, chainshift=2, mirror=2)


def tilt_classes(sq, thr=1.0, gap=3.0):
    a = sorted(t % 90 for _, _, t in sq if min(t % 90, 90 - t % 90) > thr)
    if not a:
        return 0, 0
    cl = 1 + sum(1 for u, v in zip(a, a[1:]) if v - u > gap)
    return len(a), cl


def roles(sq):
    c = collections.Counter()
    for _, _, t in sq:
        a = t % 90; d = min(a, 90 - a)
        c['A' if d < 3 else 'L' if 10 < a < 45 else 'B' if 45 <= a < 80 else 'O'] += 1
    return (c['L'], c['B'], c['A'], c['O'])


def descriptor(s, sq, k):
    nt, ncl = tilt_classes(sq)
    return (nt // 2, ncl, s < k)


def work(job):
    """One proposal: move, screen, (maybe) full polish, grid check.  Runs in a worker process."""
    s, sq, kind, seed, k, known_sides, same, smax = job
    from layout import full_lines
    rng = random.Random(seed)
    prop, desc = MOVES[kind](s, sq, rng, A)
    tmp = tempfile.mkdtemp()
    loosen = rng.choice(('1.0', '1.02', '1.05'))
    t0 = time.time()
    hop.EXTRA[:] = ['--loosen', loosen]
    r = hop.quench(s, prop, tmp, extra=('--pit', '8', '--flip-top', '0'))
    out = dict(kind=kind, desc=desc, loosen=loosen)
    if r is None:
        return dict(out, status='fail', sec=time.time() - t0)
    s1, sq1, _ = r
    if s1 > smax + 1e-3:                     # clearly above the stepping-stone ceiling: no polish
        return dict(out, status='screen-discard', s_screen=s1, sec=time.time() - t0)
    p1 = os.path.join(tmp, 's1.txt'); mcmin.write_deg(p1, s1, sq1)
    if full_lines(p1, k):                    # grid-obstructed already after the screen
        return dict(out, status='screen-grid', s_screen=s1, sec=time.time() - t0)
    near = [x for x in known_sides if abs(x - s1) < same]
    if near:
        return dict(out, status='return', s=min(near, key=lambda x: abs(x - s1)), s_screen=s1, sec=time.time() - t0)
    hop.EXTRA[:] = ['--loosen', '1.0']
    r2 = hop.quench(s1, sq1, tmp, extra=('--no-alm',))
    if r2 is None:
        return dict(out, status='fail', sec=time.time() - t0)
    s2, sq2, _ = r2
    p = os.path.join(tmp, 'q.txt'); mcmin.write_deg(p, s2, sq2)
    lines = bool(full_lines(p, k))
    return dict(out, status='new?', s=s2, sq=sq2, lines=lines, s_screen=s1, sec=time.time() - t0)


class KindBandit:
    """Thompson sampling over move kinds: reward = new below-k archive entry; rate per CPU-second with a Gamma posterior;
    counts decay with half-life `half` seconds (rates drift as the archive fills); each kind keeps a floor probability."""
    def __init__(self, kinds, half=600.0, floor=0.03, a0=1.0, b0=60.0):
        self.kinds = list(kinds); self.half, self.floor, self.a0, self.b0 = half, floor, a0, b0
        self.r = {k: 0.0 for k in kinds}; self.c = {k: 0.0 for k in kinds}; self.last = time.time()

    def decay(self):
        now = time.time(); f = 0.5 ** ((now - self.last) / self.half); self.last = now
        for k in self.kinds:
            self.r[k] *= f; self.c[k] *= f

    def update(self, kind, sec, reward):
        self.decay(); self.c[kind] += sec; self.r[kind] += reward

    def choose(self, rng):
        if rng.random() < self.floor * len(self.kinds):
            return rng.choice(self.kinds)
        draw = {k: rng.gammavariate(self.a0 + self.r[k], 1.0) / (self.b0 + self.c[k]) for k in self.kinds}
        return max(draw, key=draw.get)

    def summary(self):
        return {k: f'{self.r[k]:.1f}/{self.c[k]:.0f}s' for k in self.kinds}


class Archive:
    def __init__(self, out, n, k, smax, known):
        self.out, self.n, self.k, self.smax, self.known = out, n, k, smax, known
        self.E = []                          # entries
        self.log = open(f'{out}/archive.jsonl', 'a')
        self.plog = open(f'{out}/proposals.jsonl', 'a')

    def find(self, s, tol=2e-9):
        for e in self.E:
            if abs(e['s'] - s) < tol:
                return e
        return None

    def add(self, s, sq, parent, kind, t):
        e = self.find(s)
        if e:
            e['visits'] += 1
            return e, False
        i = len(self.E)
        path = f'{self.out}/b{i:05d}.txt'
        mcmin.write_deg(path, s, sq)
        e = dict(i=i, s=s, path=path, desc=descriptor(s, sq, self.k), roles=roles(sq), parent=parent, kind=kind,
                 visits=1, expanded=0, children=0, t=round(t, 1),
                 known=any(abs(s - x) < 2e-9 for x in self.known))
        self.E.append(e)
        self.log.write(json.dumps({k: v for k, v in e.items()}) + '\n'); self.log.flush()
        return e, True

    def dump(self, stats, t):
        tmp = f'{self.out}/state.json.tmp'
        json.dump(dict(t=t, stats=dict(stats), E=self.E), open(tmp, 'w'))
        os.replace(tmp, f'{self.out}/state.json')

    def load(self, path):
        if not os.path.exists(path):         # runs from before state dumps: rebuild from the archive log
            E = [json.loads(l) for l in open(f'{self.out}/archive.jsonl')]
            self.E = E
            return max([e['t'] for e in E] + [0.0]), collections.Counter()
        d = json.load(open(path))
        self.E = d['E']
        return d.get('t', 0.0), collections.Counter(d.get('stats', {}))

    def pick(self, rng, topk=8):
        cells = collections.defaultdict(list)
        for e in self.E:
            if e['s'] < self.smax:
                cells[tuple(e['desc'])].append(e)
        cell = rng.choice(list(cells))
        es = sorted(cells[cell], key=lambda e: e['s'])[:topk]
        return rng.choices(es, weights=[1.0 / (1 + e['expanded']) ** 2 for e in es])[0]


def report(out):
    E = [json.loads(l) for l in open(f'{out}/archive.jsonl')]
    k = math.ceil(math.sqrt(110))
    sub = [e for e in E if e['s'] < k]
    new = [e for e in sub if not e['known']]
    print(f'{out}: {len(E)} basins, {len(sub)} below {k} ({len(new)} not in the known 48); best {min(e["s"] for e in E):.10f}')
    by = collections.Counter(tuple(e['roles']) for e in sub)
    print('  sub-k role groups (L, B, axis, other):', dict(by.most_common()))
    print('  sub-k descriptor cells:', dict(collections.Counter(tuple(e['desc']) for e in sub).most_common()))
    print('  new sub-k basins by kind:', dict(collections.Counter(e['kind'] for e in new)))
    for e in sorted(new, key=lambda e: e['s'])[:12]:
        print(f'    {e["s"]:.10f} roles {tuple(e["roles"])} desc {tuple(e["desc"])} from b{e["parent"]} by {e["kind"]} at {e["t"]}s')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=110); ap.add_argument('--starts', nargs='+')
    ap.add_argument('--minutes', type=float, default=60); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--smax', type=float, default=None, help='stepping-stone ceiling (default k + 0.05)')
    ap.add_argument('--same', type=float, default=1e-6, help='screened side within this of a known basin = return')
    ap.add_argument('--star', action='store_true', help='control: always expand the first start (cen7-style sampling)')
    ap.add_argument('--adapt', action='store_true', help='Thompson sampling over move kinds (reward: new below-k basin)')
    ap.add_argument('--resume', action='store_true', help='continue from <out>/state.json (time offset carried over)')
    ap.add_argument('--out'); ap.add_argument('--report'); ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    if a.report:
        report(a.report); raise SystemExit
    os.makedirs(a.out, exist_ok=True)
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    known = json.load(open(f'{HERE}/runs/known{a.n}.json')) if os.path.exists(f'{HERE}/runs/known{a.n}.json') else []
    ar = Archive(a.out, a.n, k, a.smax or k + 0.05, known)
    rng = random.Random(a.seed)
    tmp = tempfile.mkdtemp()
    t0 = time.time()
    stats = collections.Counter()
    if a.resume:
        toff, stats = ar.load(f'{a.out}/state.json')
        t0 -= toff
        a.starts = []
        print(f'resumed {len(ar.E)} basins at t = {toff:.0f}s', flush=True)
    hop.EXTRA[:] = ['--loosen', '1.0']
    for p in a.starts:                       # starts are polished in place (no ALM, no loosen: keep their basins)
        s, sq = mcmin.load_deg(p)
        r = hop.quench(s, sq, tmp, extra=('--no-alm',))
        if r:
            e, _ = ar.add(r[0], r[1], -1, 'start:' + os.path.basename(p), 0)
            print('start', p, f'{r[0]:.10f}', e['desc'], e['roles'], flush=True)
    kinds, wts = zip(*WEIGHTS.items())
    bandit = KindBandit([kk for kk, w in WEIGHTS.items() if w > 0]) if a.adapt else None
    pend = {}
    t_start = time.time()
    with ProcessPoolExecutor(a.procs) as ex:
        def submit():
            par = ar.E[0] if a.star else ar.pick(rng)
            par['expanded'] += 1
            s, sq = mcmin.load_deg(par['path'])
            kind = bandit.choose(rng) if bandit else rng.choices(kinds, weights=wts)[0]
            f = ex.submit(work, (s, sq, kind, rng.randrange(1 << 30), k, [e['s'] for e in ar.E], a.same, ar.smax))
            pend[f] = par
        for _ in range(a.procs):
            submit()
        last = 0
        while pend:
            done, _ = wait(pend, return_when=FIRST_COMPLETED)
            for f in done:
                par = pend.pop(f)
                try:
                    r = f.result()
                except Exception as exn:
                    stats['error'] += 1; r = None
                if r:
                    stats[r['status']] += 1
                    newi = None
                    if r['status'] == 'return':
                        e = ar.find(r['s'], tol=1e-9)
                        if e: e['visits'] += 1
                    elif r['status'] == 'new?':
                        if r['lines'] or r['s'] >= ar.smax:
                            stats['discard'] += 1
                        else:
                            e, isnew = ar.add(r['s'], r['sq'], par['i'], r['kind'], time.time() - t0)
                            stats['new' if isnew else 'dup'] += 1
                            if isnew:
                                par['children'] += 1
                                newi = e['i']
                    if bandit:
                        bandit.update(r['kind'], r.get('sec', 0.0), 1.0 if (newi is not None and r['s'] < k) else 0.0)
                    ar.plog.write(json.dumps(dict(t=round(time.time() - t0, 1), kind=r['kind'], st=r['status'],
                                                  sec=round(r.get('sec', 0), 2), s=r.get('s'), new=newi, par=par['i'])) + '\n')
                if time.time() - t_start < 60 * a.minutes:
                    submit()
            if time.time() - last > 120:
                last = time.time()
                ar.dump(stats, time.time() - t0)
                sub = [e for e in ar.E if e['s'] < k]
                if bandit: print('  bandit', bandit.summary(), flush=True)
                print(f'{time.time() - t0:7.0f}s basins {len(ar.E)} sub-k {len(sub)} (new vs known {sum(not e["known"] for e in sub)}) '
                      f'best {min(e["s"] for e in ar.E):.10f} cells {len({tuple(e["desc"]) for e in ar.E})} | {dict(stats)}', flush=True)
    ar.dump(stats, time.time() - t0)
    report(a.out)
