#!/usr/bin/env python3
"""Mondrian composer (10-10; after itsnaka's SQUISH "Mondrian composition", jlevy/squares#481).  Log: MONDRIAN.md.

A *piece* is a corner crop of a known packing: the squares lying wholly inside [0, W] x [0, H] at one of the source
box's four corners (W, H tight: the crop's extent).  Pieces that beat the grid (m > floor W * floor H) form the atlas.
Composition puts pieces in corners of a square box and fills the rest with whole rows / columns of upright squares, so
the side is a piece dimension plus an integer (exact).  Then the seams are repaired (fq polish) and the result is
compared with the frontier.

  mondrian.py atlas --out runs/mondrian/atlas_all.npz [--exclude-src pending:#481 ...] [--procs 3]
  mondrian.py t1 --atlas runs/mondrian/atlas_all.npz [--nmin 50 --nmax 400] [--margin 5e-3]
  mondrian.py t3 --atlas ... (two pieces; guillotine cut between them)
  mondrian.py build --atlas ... --layout JSON --out FILE   (coordinates, exact axis fill, fq repair)

Layouts (T1): piece in the lower-left corner, box side s, L of uprights:  variant 'h' = top strip of floor(s-H) full
rows + floor(s-W) columns of height floor(H) right of the piece; variant 'v' = the transpose.
"""
from __future__ import annotations
import argparse, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
KINGBIRD = os.path.join(REPO, 'site', 'data', 'packings.json')
FT = 1e-9          # floor tolerance


def fl(x):
    return np.floor(np.asarray(x) + FT)


# ---------------------------------------------------------------- sources

def kingbird_sources():
    """Ellsworth's catalogue: every packing sighted from source 'kingbird' (all snapshots, revisions included;
    pk.py sync-kingbird).  Falls back to the 09-30 site mirror if the store has no catalogue rows."""
    from pk.store import Store
    st = Store()
    rows = st.db.execute("SELECT p.id, min(g.ref) AS ref FROM packing p JOIN sighting g ON g.packing = p.id "
                         "WHERE g.source = 'kingbird' AND p.status != 'infeasible' GROUP BY p.id").fetchall()
    out = []
    for r in rows:
        p = st.get(r['id'])
        out.append(dict(label='kingbird:' + r['ref'].replace('.svg', ''), n=p.n, s=float(p.s), sq=np.array(p.sq)))
    if out:
        return out
    d = json.load(open(KINGBIRD))
    for name, v in d.items():
        if v['errors'] or len(v['squares']) != int(v['n']):
            continue
        sq = np.array([[float(a) for a in r] for r in v['squares']])
        out.append(dict(label='kingbird:' + name.replace('.svg', ''), n=int(v['n']), s=float(v['s']), sq=sq))
    return out


def store_sources(top_runs=30, exclude=()):
    """Store packings: every non-run sighting (register, pending issues, outside pools), every certified minimum,
    and the best `top_runs` distinct sides per n from our runs.  `exclude`: sighting-source prefixes to drop (a packing
    seen in an excluded source is dropped even if seen elsewhere)."""
    from pk.store import Store
    st = Store()
    db = st.db
    rows = db.execute("SELECT p.id, p.n, p.s, p.status, group_concat(DISTINCT g.source) AS src FROM packing p "
                      "JOIN sighting g ON g.packing = p.id WHERE p.status != 'infeasible' GROUP BY p.id").fetchall()
    keep, per_n = {}, {}
    for r in rows:
        srcs = (r['src'] or '').split(',')
        if any(s.startswith(e) for s in srcs for e in exclude):
            continue
        if srcs and all(s == 'kingbird' for s in srcs):
            continue                                  # catalogue-only rows come in via kingbird_sources()
        ext = any(not s.startswith('run:') and s != 'kingbird' for s in srcs)
        if ext or r['status'] == 'certified':
            keep[r['id']] = ','.join(sorted(set(s for s in srcs if not s.startswith('run:')))[:3]) or r['status']
        else:
            per_n.setdefault(r['n'], []).append((r['s'], r['id']))
    for n, lst in per_n.items():
        lst.sort()
        seen = []
        for s, pid in lst:
            if len(seen) >= top_runs:
                break
            if all(abs(s - t) > 1e-7 for t in seen):
                seen.append(s); keep.setdefault(pid, 'run')
    out = []
    for pid, src in keep.items():
        p = st.get(pid)
        out.append(dict(label=f'store:{pid}:{src}', n=p.n, s=float(p.s), sq=np.array(p.sq)))
    return out


# ---------------------------------------------------------------- pieces

def orient(sq, s, o):
    """o in 0..3: which source corner goes to the origin (bit 0 mirror x, bit 1 mirror y); bit 2 = transpose."""
    x, y, a = sq[:, 0].copy(), sq[:, 1].copy(), sq[:, 2].copy()
    if o & 1: x = s - x; a = -a
    if o & 2: y = s - y; a = -a
    if o & 4: x, y = y, x; a = -a
    return np.c_[x, y, a % 90.0]


def half_ext(a_deg):
    r = np.radians(a_deg)
    return 0.5 * (np.abs(np.cos(r)) + np.abs(np.sin(r)))


def pieces_of(sq, s, min_excess=1):
    """All tight beyond-grid corner crops of one packing: arrays W, H, m, o (W >= H; o has bit 2 when transposed)."""
    out = []
    for o in range(4):
        q = orient(sq, s, o)
        h = half_ext(q[:, 2])
        u, v = q[:, 0] + h, q[:, 1] + h
        us, vs = np.unique(np.round(u, 12)), np.unique(np.round(v, 12))
        iu, iv = np.searchsorted(us, np.round(u, 12)), np.searchsorted(vs, np.round(v, 12))
        C = np.zeros((len(us), len(vs)), np.int32)
        np.add.at(C, (iu, iv), 1)
        M = C.cumsum(0).cumsum(1)
        # tight: growing either threshold to this value adds a square
        Mu = np.vstack([np.zeros((1, len(vs)), np.int32), M[:-1]])
        Mv = np.hstack([np.zeros((len(us), 1), np.int32), M[:, :-1]])
        ex = M - (fl(us)[:, None] * fl(vs)[None, :]).astype(np.int64)
        ok = (M > Mu) & (M > Mv) & (ex >= min_excess)
        I, J = np.nonzero(ok)
        W, H, m = us[I], vs[J], M[I, J]
        oo = np.full(len(I), o)
        sw = W < H
        W2, H2 = np.where(sw, H, W), np.where(sw, W, H)
        oo = np.where(sw, oo | 4, oo)
        out.append((W2, H2, m, oo))
    W, H, m, o = (np.concatenate(z) for z in zip(*out))
    return W, H, m, o


def pareto_prune(W, H, m, *rest):
    """Keep, per (m, floor W, floor H), only (W, H) not dominated by another piece of the same key.  Within one integer
    cell and count a smaller piece is never worse for the composer (its L of uprights is the same or larger)."""
    fW, fH = fl(W).astype(np.int64), fl(H).astype(np.int64)
    order = np.lexsort((H, W, fH, fW, m))
    keep = np.zeros(len(W), bool)
    last_key, minH = None, math.inf
    for i in order:
        key = (m[i], fW[i], fH[i])
        if key != last_key:
            last_key, minH = key, math.inf
        if H[i] < minH - 1e-12:
            keep[i] = True; minH = H[i]
    return [a[keep] for a in (W, H, m) + rest]


def _pieces_job(args):
    k, s, sq = args
    W, H, m, o = pieces_of(sq, s)
    if len(W) == 0:
        return k, None
    W, H, m, o = pareto_prune(W, H, m, o)
    return k, (W, H, m, o)


def cmd_atlas(a):
    t0 = time.time()
    src = []
    if not a.no_kingbird:
        src += kingbird_sources()
    if not a.no_store:
        src += store_sources(a.top_runs, tuple(a.exclude_src))
    print(f'{len(src)} source packings ({time.time() - t0:.0f} s)', flush=True)
    parts = []
    jobs = [(k, p['s'], p['sq']) for k, p in enumerate(src)]
    with ProcessPoolExecutor(a.procs) as ex:
        for k, r in ex.map(_pieces_job, jobs, chunksize=8):
            if r is not None:
                parts.append(r + (np.full(len(r[0]), k),))
    W, H, m, o, sid = (np.concatenate(z) for z in zip(*parts))
    print(f'raw pieces {len(W)}', flush=True)
    W, H, m, o, sid = pareto_prune(W, H, m, o, sid)
    print(f'pruned pieces {len(W)}  ({time.time() - t0:.0f} s)')
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    labels = np.array([p['label'] for p in src])
    np.savez_compressed(a.out, W=W, H=H, m=m, o=o, sid=sid, labels=labels,
                        src_n=np.array([p['n'] for p in src]), src_s=np.array([p['s'] for p in src]),
                        sq=_objarr([p['sq'] for p in src]))


def _objarr(lst):
    out = np.empty(len(lst), object)
    for i, x in enumerate(lst):
        out[i] = x
    return out


def load_atlas(path):
    z = np.load(path, allow_pickle=True)
    return {k: z[k] for k in z.files}


# ---------------------------------------------------------------- T1

def t1_counts(W, H, m, s):
    """Exact squares in the T1 layout (piece in a corner of the s-box, L of uprights), max of the two L variants."""
    fs = fl(s)
    Lh = fl(s - H) * fs + fl(s - W) * fl(H)
    Lv = fl(s - W) * fs + fl(s - H) * fl(W)
    ok = (s >= W - 1e-12)
    return np.where(ok, m + np.maximum(Lh, Lv), -1), np.where(Lh >= Lv, 0, 1)


def frontier_best():
    from pk.store import Store
    st = Store()
    out = {}
    for r in st.db.execute('SELECT f.n, f.s, f.finder, p.s AS ps FROM frontier f LEFT JOIN packing p ON p.id=f.packing'):
        v = min(r['s'], r['ps']) if r['ps'] is not None else r['s']
        if r['n'] not in out or v < out[r['n']][0]:
            out[r['n']] = (v, r['finder'])
    for r in st.db.execute("SELECT n, min(s) AS s FROM packing WHERE status != 'infeasible' AND pen <= 1e-9 GROUP BY n"):
        if r['n'] not in out or r['s'] < out[r['n']][0]:
            out[r['n']] = (r['s'], 'store')
    # beyond the store's frontier: Ellsworth's catalogue (current file names only) and the grid bound
    for name, v in json.load(open(KINGBIRD)).items():
        if v['errors'] or '_' in name or not name.startswith('square-'):
            continue
        n, sv = int(v['n']), float(v['s'])
        if n not in out or sv < out[n][0]:
            out[n] = (sv, 'kingbird-09-30')
    for n in range(1, 4000):
        g = math.ceil(math.sqrt(n) - 1e-12)
        if n not in out or g < out[n][0]:
            out[n] = (float(g), 'grid')
    return out


def cmd_t1(a):
    A = load_atlas(a.atlas)
    W, H, m = A['W'], A['H'], A['m']
    front = frontier_best()
    best = {}                       # n -> (s, piece index, variant, drop)
    for base in ('W', 'H'):
        D = W if base == 'W' else H
        for k in range(a.kmax + 1):
            s = D + k
            cnt, var = t1_counts(W, H, m, s)
            for drop in range(a.drop + 1):
                nn = cnt - drop
                sel = (nn >= a.nmin) & (nn <= a.nmax) & (A['src_n'][A['sid']] != nn)
                for i in np.nonzero(sel)[0]:
                    n = int(nn[i])
                    if n not in best or s[i] < best[n][0] - 1e-12:
                        best[n] = (float(s[i]), int(i), int(var[i]), drop)
    rep = []
    for n in sorted(best):
        s, i, var, drop = best[n]
        f = front.get(n, (math.inf, None))
        rep.append(dict(n=n, s=s, gap=s - f[0], front=f[0], finder=f[1], piece=i, var=var, drop=drop,
                        src=str(A['labels'][A['sid'][i]]), m=int(A['m'][i]), W=float(W[i]), H=float(H[i])))
    out = [r for r in rep if r['gap'] < a.margin and abs(r['gap']) > 1e-9]
    for r in out:
        print(f"{r['n']:4d} s {r['s']:.9f} gap {r['gap']:+.2e} (front {r['front']:.9f} {r['finder']}) "
              f"piece m={r['m']} {r['W']:.4f}x{r['H']:.4f} from {r['src']}")
    print(f'{len(rep)} n composed; {sum(r["gap"] < -1e-9 for r in rep)} below frontier; '
          f'{sum(abs(r["gap"]) <= 1e-9 for r in rep)} ties; {len(out)} within margin')
    if a.json:
        json.dump(rep, open(a.json, 'w'), indent=0)


# ---------------------------------------------------------------- T1 / T2 / T3 by one guillotine cut

class Composer:
    """best(s): most squares in an s x s box made of two columns (widths c, s - c, full height), each column holding at
    most one piece in a corner (bottom of the left column, top of the right one) plus upright rows / columns around it.
    c = s gives T1; a small second piece gives T2 (corner patch); two big ones T3.  best(s) is nondecreasing in s."""

    def __init__(self, A, smax):
        W, H, m = A['W'], A['H'], A['m']
        sel = H <= smax
        idx = np.nonzero(sel)[0]
        # both orientations of every piece: (w, h) = column-width extent, height extent
        self.w = np.concatenate([W[idx], H[idx]])
        self.h = np.concatenate([H[idx], W[idx]])
        self.m = np.concatenate([m[idx], m[idx]]).astype(np.int64)
        self.pid = np.concatenate([idx, idx])
        self.tr = np.concatenate([np.zeros(len(idx), int), np.ones(len(idx), int)])
        sn = A['src_n'][A['sid']]
        self.srcn = np.concatenate([sn[idx], sn[idx]])
        self.skip_n = None                      # pieces cut from a packing of this n are ignored (no self-reproduction)

    def column(self, s):
        """Step function of the best one-piece column of height s: breakpoints c (sorted), value f(c), argmax info."""
        ok = (self.w <= s + 1e-12) & (self.h <= s + 1e-12)
        if self.skip_n is not None:
            ok &= self.srcn != self.skip_n
        w, h, m = self.w[ok], self.h[ok], self.m[ok]
        ids = np.nonzero(ok)[0]
        fs = fl(s)
        cs, vs, js = [], [], []
        kmax = int(fs) + 1
        for i in range(kmax + 1):
            for c in (w + i, np.ceil(w + i - 1e-9)):       # c = w + i, and the next integer at or above it
                good = c <= s + 1e-12
                if not good.any():
                    continue
                cc = c[good]; ww = w[good]; hh = h[good]; mm = m[good]
                v1 = fl(s - hh) * fl(cc) + fl(cc - ww) * fl(hh)
                v2 = fl(cc - ww) * fs + fl(s - hh) * fl(ww)
                cs.append(cc); vs.append(mm + np.maximum(v1, v2).astype(np.int64)); js.append(ids[good])
        ci = np.arange(0, int(fs) + 1, dtype=float)                # piece-free columns
        cs.append(ci); vs.append((fl(ci) * fs).astype(np.int64)); js.append(np.full(len(ci), -1))
        c, v, j = np.concatenate(cs), np.concatenate(vs), np.concatenate(js)
        o = np.lexsort((-v, c))
        c, v, j = c[o], v[o], j[o]
        run = np.maximum.accumulate(v)
        keep = np.r_[True, run[1:] > run[:-1]]                       # first c reaching each new level
        return c[keep], run[keep], j[keep]

    def best(self, s):
        c, v, j = self.column(s)
        d = s - c
        k = np.searchsorted(c, d + 1e-12, side='right') - 1           # best right column of width <= s - c
        tot = v + np.where(k >= 0, v[np.maximum(k, 0)], 0)
        i = int(np.argmax(tot))
        kk = int(k[i])
        return int(tot[i]), dict(c=float(c[i]), left=int(j[i]), right=int(j[kk]) if kk >= 0 else -1,
                                 left_count=int(v[i]), right_count=int(v[kk]) if kk >= 0 else 0)

    def info(self, jj):
        if jj < 0:
            return None
        return dict(piece=int(self.pid[jj]), tr=int(self.tr[jj]), w=float(self.w[jj]), h=float(self.h[jj]),
                    m=int(self.m[jj]))


def _compose_n(args):
    path, n, lo, hi, tol = args
    A = _ATLAS_CACHE.get(path)
    if A is None:
        A = _ATLAS_CACHE[path] = load_atlas(path)
        _ATLAS_CACHE[path + ':C'] = Composer(A, 25.0)
    C = _ATLAS_CACHE[path + ':C']
    C.skip_n = n
    cnt, _ = C.best(hi)
    if cnt < n:
        return n, None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if C.best(mid)[0] >= n:
            hi = mid
        else:
            lo = mid
    cnt, info = C.best(hi)
    info['left'] = C.info(info['left']); info['right'] = C.info(info['right'])
    for side in ('left', 'right'):
        if info[side]:
            info[side]['src'] = str(A['labels'][A['sid'][info[side]['piece']]])
    return n, dict(s=hi, count=cnt, **info)


_ATLAS_CACHE = {}


def cmd_compose(a):
    front = frontier_best()
    ns = [n for n in range(a.nmin, a.nmax + 1) if n in front and (front[n][0] < math.ceil(math.sqrt(n) - 1e-12)
          or not a.subgrid_only)]
    jobs = [(a.atlas, n, front[n][0] - a.below, front[n][0] + a.margin, a.tol) for n in ns]
    rep = []
    t0 = time.time()
    with ProcessPoolExecutor(a.procs) as ex:
        for n, r in ex.map(_compose_n, jobs):
            if r is None:
                continue
            r.update(n=n, front=front[n][0], finder=front[n][1], gap=r['s'] - front[n][0])
            rep.append(r)
            L, R = r['left'], r['right']
            desc = ' + '.join(f"{x['m']}@{x['w']:.4f}x{x['h']:.4f} {x['src']}" for x in (L, R) if x)
            print(f"{n:4d} s {r['s']:.9f} gap {r['gap']:+.2e} count {r['count']} c={r['c']:.4f}  {desc}", flush=True)
    rep.sort(key=lambda r: r['n'])
    print(f'{len(jobs)} n searched, {len(rep)} with a layout within margin, '
          f'{sum(r["gap"] < -1e-9 for r in rep)} below the frontier ({time.time() - t0:.0f} s)')
    if a.json:
        json.dump(rep, open(a.json, 'w'), indent=1)


# ---------------------------------------------------------------- two-piece placement with exact fill (replication)

KB_SVG = os.path.join(REPO, 'site', 'data', 'ellsworth', 'svg')


def kingbird_svg(name):
    """A Kingbird SVG fetched after the 09-30 mirror (svg/new_<name>.svg), parsed by the site's validated parser."""
    sys.path.insert(0, os.path.join(REPO, 'site', 'tools'))
    from parse_svg import parse
    from pk.packing import Packing
    s, sq = parse(os.path.join(KB_SVG, 'new_' + name + '.svg'))
    sq = np.array([[float(a) for a in r] for r in sq]); s = float(s)
    pen = Packing(s, sq).max_pen()
    assert pen < 1e-9, (name, pen)
    return sq, s


def load_src(spec):
    """kingbird:<name> | store:<id> | file:<path>  ->  (sq, s)."""
    kind, _, ref = spec.partition(':')
    if kind == 'kingbird':
        from pk import kingbird as KB
        got = KB.latest(ref)                      # newest catalogue snapshot in the store (pk.py sync-kingbird)
        if got is not None:
            return got
        if os.path.exists(os.path.join(KB_SVG, 'new_' + ref + '.svg')):
            return kingbird_svg(ref)
        v = json.load(open(KINGBIRD))[ref if ref.endswith('.svg') else ref + '.svg']
        return np.array([[float(a) for a in r] for r in v['squares']]), float(v['s'])
    if kind == 'store':
        from pk.store import Store
        p = Store().get(int(ref))
        return np.array(p.sq), float(p.s)
    if kind == 'file':
        from pk.packing import read
        p = read(ref)
        return np.array(p.sq), float(p.s)
    raise ValueError(spec)


def is_tilted(a, eps=1e-6):
    a = np.asarray(a) % 90.0
    return (a > eps) & (a < 90.0 - eps)


def place(sq, s_src, o, corner, i, j, S):
    """Tilted squares of a source packing (orientation o) with its box anchored at a corner of the S-box and pushed
    outward by (i, j) units (i.e. i columns / j rows of the source dropped); only squares wholly inside are kept."""
    q = orient(sq, s_src, o)
    q = q[is_tilted(q[:, 2])].copy()
    d = S - s_src
    dx = {'LL': -i, 'UL': -i, 'UR': d + i, 'LR': d + i}[corner]
    dy = {'LL': -j, 'LR': -j, 'UR': d + j, 'UL': d + j}[corner]
    q[:, 0] += dx; q[:, 1] += dy
    h = half_ext(q[:, 2])
    keep = (q[:, 0] - h >= -1e-9) & (q[:, 0] + h <= S + 1e-9) & (q[:, 1] - h >= -1e-9) & (q[:, 1] + h <= S + 1e-9)
    return q[keep]


def _sat_overlap(P, Q, tol=1e-7):
    """True if some square of P overlaps some square of Q by more than tol (SAT, all four axes)."""
    if len(P) == 0 or len(Q) == 0:
        return False
    from pk.packing import Packing
    near = np.hypot(P[:, None, 0] - Q[None, :, 0], P[:, None, 1] - Q[None, :, 1]) < math.sqrt(2) + 1e-9
    I, J = np.nonzero(near)
    if len(I) == 0:
        return False
    for i, j in zip(I, J):
        pr = Packing(100.0, np.array([P[i], Q[j]]) + [50, 50, 0])
        if pr.max_pen() > tol:
            return True
    return False


def lattice_proxy(T, S):
    """Tilted count + the best single corner lattice of upright cells clear of the tilted squares (a fast lower
    bound on the exact MILP fill, used only for ranking)."""
    k = int(math.floor(S + FT))
    best = 0
    if len(T):
        th = np.radians(T[:, 2]); c, s = np.abs(np.cos(th)), np.abs(np.sin(th)); h = 0.5 * (c + s)
        cs, sn = np.cos(th), np.sin(th)
    for fx in (False, True):
        for fy in (False, True):
            xs = (S - 0.5 - np.arange(k)) if fx else (0.5 + np.arange(k))
            ys = (S - 0.5 - np.arange(k)) if fy else (0.5 + np.arange(k))
            X, Y = np.meshgrid(xs, ys); X = X.ravel(); Y = Y.ravel()
            if len(T) == 0:
                best = max(best, len(X)); continue
            dx = X[:, None] - T[None, :, 0]; dy = Y[:, None] - T[None, :, 1]
            e = 1e-9
            ov = (np.abs(dx) < 0.5 + h - e) & (np.abs(dy) < 0.5 + h - e) \
                & (np.abs(dx * cs + dy * sn) < 0.5 + h - e) & (np.abs(-dx * sn + dy * cs) < 0.5 + h - e)
            best = max(best, int((~ov.any(1)).sum()))
    return len(T) + best


def greedy_proxy(T, S):
    """Tilted count + greedy maximal set of upright cells from all four corner lattices, clear of the tilted squares,
    taken in order of distance from each cell's anchor corner (close to the exact MILP fill; ranking only)."""
    k = int(math.floor(S + FT))
    ar = np.arange(k)
    xs = [(0.5 + ar, 0), (S - 0.5 - ar, 1)]
    cells = []
    for X1, fx in xs:
        for Y1, fy in xs:
            X, Y = np.meshgrid(X1, Y1)
            dist = np.meshgrid(ar, ar)
            cells.append(np.c_[X.ravel(), Y.ravel(), (dist[0] + dist[1]).ravel()])
    C = np.concatenate(cells)
    if len(T):
        th = np.radians(T[:, 2]); cs, sn = np.cos(th), np.sin(th); h = 0.5 * (np.abs(cs) + np.abs(sn))
        dx = C[:, None, 0] - T[None, :, 0]; dy = C[:, None, 1] - T[None, :, 1]; e = 1e-9
        ov = (np.abs(dx) < 0.5 + h - e) & (np.abs(dy) < 0.5 + h - e) \
            & (np.abs(dx * cs + dy * sn) < 0.5 + h - e) & (np.abs(-dx * sn + dy * cs) < 0.5 + h - e)
        C = C[~ov.any(1)]
    # dedupe coincident cells, order by anchor distance
    C = C[np.lexsort((C[:, 1], C[:, 0], C[:, 2]))]
    key = np.round(C[:, :2], 9)
    _, first = np.unique(key, axis=0, return_index=True)
    C = C[np.sort(first)]
    taken = np.zeros(len(C), bool)
    blocked = np.zeros(len(C), bool)
    X, Y = C[:, 0], C[:, 1]
    for i in range(len(C)):
        if blocked[i]:
            continue
        taken[i] = True
        blocked |= (np.abs(X - X[i]) < 1 - 1e-9) & (np.abs(Y - Y[i]) < 1 - 1e-9)
    return len(T) + int(taken.sum())


def exact_fill(T, S, tol=1e-6, time_limit=30):
    import gen
    F = gen.mis_fill2(S, [tuple(r) for r in T], tol=tol, time_limit=time_limit)
    return np.r_[T, np.array(F, float).reshape(-1, 3)]


def two_piece_search(n, S, srcA, srcB, R=2, corners_B=('UR', 'LR'), topk=40, fq_eng=None, verbose=False):
    """All orientations / offsets of piece A at LL and piece B at UR or LR in the S-box; rank by lattice proxy, exact
    fill the top ones; return sorted results [(count, proxy, desc, layout)]."""
    sqA, sA = srcA; sqB, sB = srcB
    PA, PB = {}, {}
    for o in range(8):
        for i in range(R + 1):
            for j in range(R + 1):
                T = place(sqA, sA, o, 'LL', i, j, S)
                PA.setdefault(T.round(9).tobytes(), (T, ('LL', o, i, j)))
                for cb in corners_B:
                    T = place(sqB, sB, o, cb, i, j, S)
                    PB.setdefault(T.round(9).tobytes(), (T, (cb, o, i, j)))
    combos = []
    for TA, da in PA.values():
        for TB, db in PB.values():
            T = np.r_[TA, TB]
            combos.append((lattice_proxy(T, S), da, db, TA, TB))
    combos.sort(key=lambda c: -c[0])
    out = []
    for prox, da, db, TA, TB in combos[:topk * 4]:
        if len(out) >= topk:
            break
        if _sat_overlap(TA, TB):
            continue
        L = exact_fill(np.r_[TA, TB], S)
        out.append((len(L), prox, (da, db), L))
        if verbose:
            print(f'  proxy {prox} exact {len(L)} A{da} B{db}', flush=True)
    out.sort(key=lambda r: -r[0])
    return out


def crop_pieces(src, m, R=3, mtol=0):
    """Rectangular sub-packings of a source: orientation o (8), i columns / j rows dropped at the wall side, tight inner
    cut holding m (+- mtol) squares.  Returns [(T, desc)]: the crop's tilted squares, rectangle corner at the origin."""
    sq, s = src
    out = {}
    for o in range(8):
        q = orient(sq, s, o)
        h = half_ext(q[:, 2])
        for i in range(R + 1):
            for j in range(R + 1):
                x, y = q[:, 0] - i, q[:, 1] - j
                ins = (x - h >= -1e-9) & (y - h >= -1e-9)
                if ins.sum() < m - mtol:
                    continue
                u, v = np.round((x + h)[ins], 9), np.round((y + h)[ins], 9)
                us, vs = np.unique(u), np.unique(v)
                C = np.zeros((len(us), len(vs)), np.int32)
                np.add.at(C, (np.searchsorted(us, u), np.searchsorted(vs, v)), 1)
                M = C.cumsum(0).cumsum(1)
                Mu = np.vstack([np.zeros((1, len(vs)), np.int32), M[:-1]])
                Mv = np.hstack([np.zeros((len(us), 1), np.int32), M[:, :-1]])
                ok = (np.abs(M - m) <= mtol) & (M > Mu) & (M > Mv)
                qq = np.c_[x[ins], y[ins], q[ins, 2]]
                til = is_tilted(qq[:, 2])
                for I, J in zip(*np.nonzero(ok)):
                    sel = (u <= us[I] + 1e-9) & (v <= vs[J] + 1e-9) & til
                    T = qq[sel]
                    out.setdefault(T.round(9).tobytes(), (T, (o, i, j, float(us[I]), float(vs[J]), int(M[I, J]))))
    return list(out.values())


def to_corner(T, corner, S):
    """Move a piece from the lower-left frame to a corner of the S-box (UR: 180-degree turn; LR: mirror in x)."""
    T = T.copy()
    if corner in ('UR', 'LR'):
        T[:, 0] = S - T[:, 0]; T[:, 2] = (-T[:, 2]) % 90.0
    if corner == 'UR':
        T[:, 1] = S - T[:, 1]; T[:, 2] = (-T[:, 2]) % 90.0
    return T


def fits(T, S):
    h = half_ext(T[:, 2])
    return bool(((T[:, 0] - h >= -1e-9) & (T[:, 0] + h <= S + 1e-9) & (T[:, 1] - h >= -1e-9) &
                 (T[:, 1] + h <= S + 1e-9)).all())


def trim(L, n):
    """Keep n squares: every tilted one, then the uprights nearest to a tilted square (surplus dropped far away)."""
    if len(L) <= n:
        return L
    t = is_tilted(L[:, 2])
    T, U = L[t], L[~t]
    if len(T) == 0:
        return L[:n]
    d = np.hypot(U[:, None, 0] - T[None, :, 0], U[:, None, 1] - T[None, :, 1]).min(1)
    return np.r_[T, U[np.argsort(d)[: n - len(T)]]]


def _drop_conflicts(keep, other, tol=1e-7):
    """Tilted squares of `other` that overlap no square of `keep`."""
    if len(keep) == 0 or len(other) == 0:
        return other
    from pk.packing import Packing
    ok = np.ones(len(other), bool)
    near = np.hypot(other[:, None, 0] - keep[None, :, 0], other[:, None, 1] - keep[None, :, 1]) < math.sqrt(2) + 1e-9
    for i, j in zip(*np.nonzero(near)):
        if ok[i] and Packing(100.0, np.array([other[i], keep[j]]) + [50, 50, 0]).max_pen() > tol:
            ok[i] = False
    return other[ok]


def _fill_job(args):
    T, S = args[:2]
    tol = args[2] if len(args) > 2 else 1e-6
    return len(exact_fill(T, S, tol=tol, time_limit=10))


def placements(src, corner, S, R):
    out = {}
    for o in range(8):
        for i in range(R + 1):
            for j in range(R + 1):
                T = place(*src, o, corner, i, j, S)
                if len(T):
                    out.setdefault(T.round(9).tobytes(), (T, (corner, o, i, j)))
    return list(out.values())


def replicate(n, S, piecesA, piecesB, procs=2, top=30, fq_top=6, slack=8, fill_tol=1e-6, insert_max=0, tries=6, log=print):
    """Two pieces (lists from crop_pieces), A at LL and B at UR / LR; conflicts resolved by dropping the other piece's
    tilted squares (both priority orders).  Single-piece MILP excess -> additive estimate -> joint MILP on the best
    pairs -> fq repair of the layouts that reach n.  Returns list of dicts."""
    from pk.engine import FQ, QOpt, POLISH
    from pk.packing import Packing
    grid = int(math.floor(S + FT)) ** 2
    PA = [(to_corner(T, 'LL', S), ('LL',) + d) for T, d in piecesA if fits(T, S)]
    PB = [(to_corner(T, c, S), (c,) + d) for c in ('UR', 'LR') for T, d in piecesB if fits(to_corner(T, c, S), S)]
    if not PA or not PB:
        log(f'  n={n} S={S:.10f}: no fitting pieces'); return []
    # pair screen: rectangle formula (crop totals + guillotine upright count; undercounts by up to ~8)
    WA = np.array([d[4] for _, d in PA]); HA = np.array([d[5] for _, d in PA]); mA = np.array([d[6] for _, d in PA])
    WB = np.array([d[4] for _, d in PB]); HB = np.array([d[5] for _, d in PB]); mB = np.array([d[6] for _, d in PB])

    def col(w, wp, hp, m):          # a full-height column of width w with one piece (wp x hp, m squares) in a corner
        ok = w >= wp - 1e-12
        v = m + np.maximum(fl(S - hp) * fl(w) + fl(w - wp) * fl(hp), fl(w - wp) * fl(S) + fl(S - hp) * fl(wp))
        return np.where(ok, v, -10 ** 6)

    best = np.full((len(PA), len(PB)), -10 ** 6, float)
    for tr in (False, True):        # vertical cut, or horizontal cut (transpose both rectangles)
        a_w, a_h = (HA, WA) if tr else (WA, HA)
        b_w, b_h = (HB, WB) if tr else (WB, HB)
        for k in range(0, 3):
            for left in (True, False):
                c = (a_w[:, None] + k) if left else (S - b_w[None, :] - k)
                cnt = col(c, a_w[:, None], a_h[:, None], mA[:, None]) + col(S - c, b_w[None, :], b_h[None, :], mB[None, :])
                feas = (a_w[:, None] + b_w[None, :] <= S + 1e-12) & (c >= a_w[:, None] - 1e-12) & (S - c >= b_w[None, :] - 1e-12)
                best = np.maximum(best, np.where(feas, cnt, -10 ** 6))
    tA = np.array([len(T) for T, _ in PA]); tB = np.array([len(T) for T, _ in PB])
    I, J = np.nonzero(best >= n - slack)
    order = np.lexsort((-(tA[I] + tB[J]), -best[I, J]))[: 4 * top]
    pairs = [(int(best[I[k], J[k]]), 0, int(I[k]), int(J[k])) for k in order]
    log(f'  n={n} S={S:.10f}: {len(PA)} A x {len(PB)} B placements; formula best {int(best.max())}; '
        f'{len(I)} pairs >= {n} - {slack}')
    jobs, descs = [], []
    for est, _, ia, ib in pairs:
        if len(jobs) >= top:
            break
        TA, da = PA[ia]; TB, db = PB[ib]
        for order in ('A', 'B'):
            if order == 'A':
                T = np.r_[TA, _drop_conflicts(TA, TB)]
            else:
                T = np.r_[_drop_conflicts(TB, TA), TB]
            jobs.append((T, S, fill_tol)); descs.append((est, da, db, order))
    if not jobs:
        log(f'  no pair within slack')
        return []
    t0 = time.time()
    with ProcessPoolExecutor(procs) as ex:
        counts = list(ex.map(_fill_job, jobs, chunksize=2))
    log(f'  {len(jobs)} joint fills in {time.time() - t0:.0f} s')
    res = sorted(zip(counts, range(len(jobs))), key=lambda r: -r[0])
    log(f'  joint fills: best {res[0][0]} (need {n}); {sum(c >= n for c, _ in res)} of {len(jobs)} reach n')
    out = []
    eng = FQ()
    import random, mcmin
    rng = random.Random(n)
    cand = []
    for c, k in res[:fq_top]:
        if c >= n:
            cand.append((c, k, None))
        elif insert_max and c >= n - insert_max:
            for t in range(tries):
                cand.append((c, k, t))
    for c, k, t in cand:
        L = exact_fill(jobs[k][0], S, tol=fill_tol)
        Sx = S
        if t is not None:         # insertion relaxation: dilate a little, add n - c squares at pose holes, let fq squeeze
            f = 1 + rng.uniform(0.002, 0.012) * math.sqrt(n - c)
            Sx = S * f
            sq2 = [(x * f, y * f, d) for x, y, d in L]
            for _ in range(n - c):
                H = mcmin.pose_holes(Sx, sq2)
                if not H:
                    break
                e, x, y, d = H[rng.randrange(min(len(H), 12))]
                sq2.append((x, y, d))
            if len(sq2) != n:
                continue
            L = np.array(sq2)
        p = Packing(Sx, trim(L, n))
        t1 = time.time()
        r = eng.quench(p, QOpt(), timeout=120)
        s1 = r.s if r.ok else math.inf
        r2 = eng.quench(r.p, POLISH, timeout=120) if r.ok else r
        s2 = r2.s if r2.ok else s1
        out.append(dict(count=c, inserted=(n - c if t is not None else 0), est=descs[k][0], A=descs[k][1], B=descs[k][2], order=descs[k][3], pen=p.max_pen(),
                        fq=s1, polish=s2, packing=(r2.p if r2.ok else (r.p if r.ok else p))))
        log(f'    count {c}{"+" + str(n - c) + " ins" if t is not None else ""} A{descs[k][1]} B{descs[k][2]} drop-{descs[k][3]}: composed pen {p.max_pen():.1e} '
            f'-> fq {s1:.10f} -> polish {s2:.10f} ({time.time() - t1:.0f} s)')
    return out


_SQ = os.path.join(HERE, 'runs/ext/itsnaka-squish-certs-d45669b/squish-submission-2026-10-09/sources/squish-s{}-source.cert.txt')
LINEAGE_481 = {   # n: (their S, (source A, m_A), (source B, m_B))  (jlevy/squares#481 lineage table; store ids resolved 10-10)
    131: (11.949659588035860, ('store:79', 26), ('file:' + _SQ.format(155), 33)),
    153: (12.872029849081180, ('kingbird:square-152_r2', 34), ('kingbird:square-152_r2', 34)),
    207: (14.885506308841677, ('kingbird:square-11_r1a', 8), ('store:53', 130)),
    209: (14.946223654487920, ('file:' + _SQ.format(155), 42), ('file:' + _SQ.format(240), 47)),
    232: (15.767428349941842, ('kingbird:square-5', 5), ('store:72', 162)),
    236: (15.863955747159268, ('kingbird:square-152_r2', 34), ('store:33', 88)),
    259: (16.591378145497700, ('kingbird:square-50_r5', 36), ('store:46', 90)),
    263: (16.733166007899882, ('kingbird:square-5', 5), ('kingbird:square-202_r3', 188)),
    269: (16.901513582189132, ('kingbird:square-11_r1a', 8), ('store:338', 180)),
    270: (16.929780126241173, ('store:92', 57), ('store:92', 62)),
    292: (17.591378145497902, ('kingbird:square-50_r5', 43), ('store:46', 90)),
    302: (17.872029849081180, ('kingbird:square-152_r2', 34), ('store:52', 118)),
    303: (17.913065462738533, ('kingbird:square-11_r1a', 8), ('store:91', 208)),
    305: (17.951196144147385, ('store:92', 57), ('file:' + _SQ.format(306), 81)),
    307: (17.980862784054690, ('kingbird:square-241_r3', 64), ('kingbird:square-241_r3', 69)),
}


ANTE_481 = {    # register value before #481 (its "register now" column): the status quo ante each composition had to beat
    131: 11.951150044912, 153: 12.879679373329, 207: 14.887992258303, 209: 14.949617952200, 232: 15.778174593052,
    236: 15.867800839420, 259: 16.602568490493, 263: 16.740419679539, 269: 16.905967058584, 270: 16.937807228446,
    292: 17.597249391156, 302: 17.881306218096, 303: 17.920312372920, 305: 17.952959459016, 307: 17.981030548633,
}


def comp_sides(S, pieces, window=0.02, extra=4):
    """Composition sides: their S (+1e-9), then piece rectangle dims + integer in [S, S + window], nearest first."""
    out = [round(S + 1e-9, 12)]
    cand = set()
    for _, d in pieces:
        for dim in d[3:5]:
            for k in range(0, 20):
                v = dim + k
                if S + 2e-9 <= v <= S + window:
                    cand.add(round(v + 1e-9, 12))
    return out + sorted(cand)[:extra]


def explore_from(n, starts, out, minutes=10, procs=4, seed=1):
    """Hunt-recipe explorer from start files; returns the best side in its archive."""
    import subprocess
    os.makedirs(out, exist_ok=True)
    cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(n), '--starts', *starts, '--procs', str(procs),
           '--minutes', str(minutes), '--seed', str(seed), '--out', os.path.join(out, 'run'), '--adapt',
           '--elite-share', '0.3', '--reward', 'value']
    subprocess.run(cmd, stdout=open(os.path.join(out, 'log'), 'w'), stderr=subprocess.STDOUT)
    best = math.inf
    arch = os.path.join(out, 'run', 'archive.jsonl')
    if os.path.exists(arch):
        for line in open(arch):
            try:
                best = min(best, json.loads(line)['s'])
            except Exception:
                pass
    return best


def cmd_replicate(a):
    ns = a.n or sorted(LINEAGE_481)
    summary = []
    for n in ns:
        S, (sa, ma), (sb, mb) = LINEAGE_481[n]
        A_, B_ = load_src(sa), load_src(sb)
        pa = crop_pieces(A_, ma, a.R, a.mtol); pb = crop_pieces(B_, mb, a.R, a.mtol)
        print(f'n {n}: {len(pa)} pieces of {sa} (m {ma}), {len(pb)} of {sb} (m {mb})', flush=True)
        best, comp, stage = math.inf, [], '-'
        for Sc in comp_sides(S, pa + pb, a.window, a.extra):
            for X, Y in ((pa, pb), (pb, pa)) if (sa, ma) != (sb, mb) else ((pa, pb),):
                res = replicate(n, Sc, X, Y, procs=a.procs, top=a.top, fq_top=a.fq_top)
                for r in res:
                    comp.append((r['polish'], Sc, r['packing']))
                    if r['polish'] < best:
                        best, stage = r['polish'], f'compose+fq (S_c {Sc:.6f})'
                if best <= S + 1e-9:
                    break
            if best <= S + 1e-9:
                break
        if not comp and a.relax:
            print(f'  relaxation fallback: insert up to {a.insert_max} squares into the best composed layouts', flush=True)
            for Sc in comp_sides(S, pa + pb, a.window, a.extra):
                for X, Y in ((pa, pb), (pb, pa)) if (sa, ma) != (sb, mb) else ((pa, pb),):
                    res = replicate(n, Sc, X, Y, procs=a.procs, top=a.top, fq_top=a.fq_top, insert_max=a.insert_max, tries=a.tries)
                    for r in res:
                        comp.append((r['polish'], Sc, r['packing']))
                        if r['polish'] < best:
                            best, stage = r['polish'], f'relax: +{r["inserted"]} inserted at S_c {Sc:.6f} + fq'
        if best > S + 1e-9 and comp and a.explore_min > 0:
            comp.sort(key=lambda c: c[0])
            d = os.path.join(a.out, f'x{n}')
            os.makedirs(d, exist_ok=True)
            starts, seen = [], []
            for sp, Sc, pk in comp:
                if all(abs(sp - t) > 1e-9 for t in seen):
                    seen.append(sp); f = os.path.join(d, f'start{len(starts)}.txt'); pk.write(f); starts.append(f)
                if len(starts) >= 4:
                    break
            print(f'  explorer from {len(starts)} composed starts (best {comp[0][0]:.10f}), '
                  f'{a.explore_procs} procs x {a.explore_min} min', flush=True)
            bx = explore_from(n, starts, d, a.explore_min, a.explore_procs)
            if bx < best:
                best, stage = bx, f'compose+fq+explorer {a.explore_min} min'
        summary.append((n, S, best, stage))
        print(f'== n {n}: theirs {S:.10f}, ours {best:.10f}, diff {best - S:+.2e} [{stage}]', flush=True)
    print('summary:')
    for n, S, b, st in summary:
        ante = ANTE_481[n]
        tag = 'REPLICATED' if b <= S + 1e-9 else ('BEAT-ANTE' if b < ante - 1e-10 else 'no')
        print(f'  {n}: {tag:10s} ours {b:.10f}  theirs {b - S:+.2e}  ante {b - ante:+.2e}  {st}')


def cmd_selftest(a):
    """Positive control: rebuild #481's n = 259 (s(50)_r5 x 36 + ry-xu s(123) x 90) and require fq to land on their side
    within the time budget.  Exercises source loading, crop pieces, the pair screen, the exact fill and fq."""
    t0 = time.time()
    S, (sa, ma), (sb, mb) = LINEAGE_481[259]
    pa = crop_pieces(load_src(sa), ma, 3); pb = crop_pieces(load_src(sb), mb, 3)
    res = replicate(259, round(S + 1e-9, 12), pa, pb, procs=a.procs, top=150, fq_top=2, log=lambda *_: None)
    best = min((r['polish'] for r in res), default=math.inf)
    dt = time.time() - t0
    ok = best <= S + 1e-9 and dt <= a.budget
    print(f'selftest n=259: {len(pa)} x {len(pb)} pieces, {len(res)} layouts, best {best:.10f} (theirs {S:.10f}), '
          f'{dt:.0f} s (budget {a.budget:.0f} s): {"PASS" if ok else "FAIL"}')
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('atlas')
    p.add_argument('--out', required=True)
    p.add_argument('--procs', type=int, default=3)
    p.add_argument('--top-runs', type=int, default=30)
    p.add_argument('--exclude-src', nargs='*', default=[])
    p.add_argument('--no-kingbird', action='store_true')
    p.add_argument('--no-store', action='store_true')
    p = sub.add_parser('t1')
    p.add_argument('--atlas', required=True)
    p.add_argument('--nmin', type=int, default=20)
    p.add_argument('--nmax', type=int, default=400)
    p.add_argument('--kmax', type=int, default=12)
    p.add_argument('--drop', type=int, default=2)
    p.add_argument('--margin', type=float, default=5e-3)
    p.add_argument('--json')
    p = sub.add_parser('compose')
    p.add_argument('--atlas', required=True)
    p.add_argument('--nmin', type=int, default=20)
    p.add_argument('--nmax', type=int, default=324)
    p.add_argument('--subgrid-only', action='store_true', help='skip n whose best is the grid')
    p.add_argument('--margin', type=float, default=3e-3, help='report layouts up to this far above the frontier')
    p.add_argument('--below', type=float, default=0.05)
    p.add_argument('--tol', type=float, default=1e-10)
    p.add_argument('--procs', type=int, default=2)
    p.add_argument('--json')
    p = sub.add_parser('replicate', help='rebuild the 15 #481 compositions from their listed sources')
    p.add_argument('--n', type=int, nargs='*')
    p.add_argument('--R', type=int, default=3)
    p.add_argument('--procs', type=int, default=2)
    p.add_argument('--top', type=int, default=30)
    p.add_argument('--fq-top', type=int, default=6)
    p.add_argument('--window', type=float, default=0.02)
    p.add_argument('--mtol', type=int, default=0)
    p.add_argument('--extra', type=int, default=4, help='composition sides beyond their S to try')
    p.add_argument('--explore-min', type=float, default=0, help='explorer minutes after compose+fq (0 = off)')
    p.add_argument('--explore-procs', type=int, default=4)
    p.add_argument('--insert-max', type=int, default=4)
    p.add_argument('--tries', type=int, default=6, help='insertion tries per layout')
    p.add_argument('--relax', action='store_true', help='if nothing composes within the window, compose at a larger '
                   'side (up to ceil S), trim to n, and let fq + explorer squeeze')
    p.add_argument('--out', default=os.path.join(HERE, 'runs/mondrian/replicate'))
    p = sub.add_parser('selftest', help='positive control: rebuild #481 n = 259 (PASS/FAIL, exit status)')
    p.add_argument('--procs', type=int, default=2)
    p.add_argument('--budget', type=float, default=120)
    a = ap.parse_args()
    {'selftest': cmd_selftest, 'atlas': cmd_atlas, 't1': cmd_t1, 'compose': cmd_compose, 'replicate': cmd_replicate}[a.cmd](a)


if __name__ == '__main__':
    main()
