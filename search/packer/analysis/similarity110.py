#!/usr/bin/env python3
"""Packing similarity / clustering study at n = 110 (10-14).  Analysis only: reads the store READ-ONLY, writes nothing
but runs/similarity110/.  Findings: analysis/similarity110.md.

Question: which cheap measure says "same minimum", and which cluster definition gives "distinct regions" so a search
can sample parents uniformly over clusters?  Ground truth: the store's exact minima at 110 (`minimum`, `copy`).

Measures per packing
  side        stored f64 side
  img         pk.fingerprint image (2 x 24 x 24), distance = min over D4 of L1 / n
  tilt        census3 signature (tilted count, angle-cluster sizes)
  fq          fq polish (--no-alm --loosen 1.0): polished side `sp`, contact fp (FNV of load-bearing (i, j, face) keys:
              depends on square labels), ncontacts
  wl          label-free contact-graph hash: active contacts (SAT gap < 1e-8) on the fq-polished coordinates, rattlers
              stripped (iterated degree < 3 incl. walls), 3 rounds of Weisfeiler-Lehman (node label = tilted?, wall count)
Pair measures
  md          matching distance (rediscover.distance: RMS of centre + angle displacement under the best D4 x assignment)
  mv          fraction of squares that move > 0.1 under that matching (rattler-robust "how many squares differ")
  fq_rl       fq fp of B relabelled + transformed onto A's labels (best-case contact fp: labels made consistent)

  python3 analysis/similarity110.py features   # sample + per-packing features            (~2 min on 8 procs)
  python3 analysis/similarity110.py pairs      # same / different-minimum pairs            (~3 min)
  python3 analysis/similarity110.py reps       # 939 x 939 md / mv / img among representatives (~4 min)
  python3 analysis/similarity110.py edges      # parent -> child transitions from explorer archives
  python3 analysis/similarity110.py analyze    # tables -> runs/similarity110/report.txt + summary.json
  python3 analysis/similarity110.py othern --n 88 237
  python3 analysis/similarity110.py all
"""
from __future__ import annotations
import argparse, glob, hashlib, json, math, os, pickle, random, sqlite3, sys, time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
STORE = os.path.join(PK, 'runs', 'store.sqlite')
OUT = os.path.join(PK, 'runs', 'similarity110')
PROCS = 8
N = 110
COPIES = 6          # sampled copies per minimum (besides the representative)


def db():
    return sqlite3.connect(f'file:{STORE}?mode=ro', uri=True, timeout=60)


def load(con, pid):
    s, c = con.execute('SELECT s, coords FROM packing WHERE id=?', (pid,)).fetchone()
    return float(s), np.frombuffer(c, '<f8').reshape(-1, 3).copy()


def side_key(s):
    return f'{float(s):.9f}'


# ---------------------------------------------------------------- per-packing measures
def image(s, sq):
    from pk import fingerprint as F
    return F.image(s, sq).astype(np.float32)


def tilt_sig(sq):
    angs = sorted(a % 90 for a in sq[:, 2])
    tilt = [a for a in angs if 0.5 < a < 89.5]
    cl = []
    for a in tilt:
        if cl and a - cl[-1][-1] < 2.0:
            cl[-1].append(a)
        else:
            cl.append([a])
    return f'{len(tilt)}:' + ','.join(map(str, sorted((len(c) for c in cl), reverse=True)))


def corners(sq):
    t = np.radians(sq[:, 2]); c, s = np.cos(t), np.sin(t)
    loc = np.array([(0.5, 0.5), (-0.5, 0.5), (-0.5, -0.5), (0.5, -0.5)])
    X = sq[:, 0:1] + c[:, None] * loc[None, :, 0] - s[:, None] * loc[None, :, 1]
    Y = sq[:, 1:2] + s[:, None] * loc[None, :, 0] + c[:, None] * loc[None, :, 1]
    return np.stack([X, Y], -1)          # n x 4 x 2


def contact_graph(s, sq, tol=1e-8):
    """Active contacts: SAT separation gap < tol.  Returns (edges i<j, wall count per square)."""
    n = len(sq)
    C = corners(sq)
    t = np.radians(sq[:, 2])
    ax = np.stack([np.stack([np.cos(t), np.sin(t)], 1), np.stack([-np.sin(t), np.cos(t)], 1)], 1)   # n x 2 x 2
    walls = np.zeros(n, int)
    lo = C.min(1); hi = C.max(1)
    walls += (lo[:, 0] < tol).astype(int) + (lo[:, 1] < tol).astype(int)
    walls += (hi[:, 0] > s - tol).astype(int) + (hi[:, 1] > s - tol).astype(int)
    d = np.hypot(sq[:, None, 0] - sq[None, :, 0], sq[:, None, 1] - sq[None, :, 1])
    I, J = np.nonzero(np.triu(d < math.sqrt(2) + 1e-6, 1))
    E = []
    for i, j in zip(I, J):
        gap = -np.inf
        for a in (ax[i, 0], ax[i, 1], ax[j, 0], ax[j, 1]):
            pi = C[i] @ a; pj = C[j] @ a
            gap = max(gap, pj.min() - pi.max(), pi.min() - pj.max())
        if gap < tol:
            E.append((int(i), int(j)))
    return E, walls


def wl_hash(sq, E, walls, strip=True, rounds=3):
    n = len(sq)
    nb = defaultdict(set)
    for i, j in E:
        nb[i].add(j); nb[j].add(i)
    alive = set(range(n))
    if strip:                        # rattlers: iterated degree (pairs + walls) < 3
        ch = True
        while ch:
            ch = False
            for i in list(alive):
                if len(nb[i] & alive) + walls[i] < 3:
                    alive.discard(i); ch = True
    a = np.mod(sq[:, 2], 90.0); tilted = np.minimum(a, 90 - a) > 1.0
    lab = {i: f'{int(tilted[i])}{walls[i]}' for i in alive}
    for _ in range(rounds):
        lab = {i: hashlib.sha1((lab[i] + '|' + ','.join(sorted(lab[j] for j in nb[i] & alive))).encode()).hexdigest()[:12]
               for i in alive}
    h = hashlib.sha1(','.join(sorted(lab.values())).encode()).hexdigest()[:16]
    return h, len(alive), sum(1 for i, j in E if i in alive and j in alive)


_FQ = None


def polish(s, sq):
    global _FQ
    from pk.engine import FQ, POLISH
    from pk.packing import Packing
    if _FQ is None:
        _FQ = FQ()
    r = _FQ.quench(Packing(s, sq), POLISH, timeout=120)
    if r.p is None:
        return None
    return float(r.p.s), np.asarray(r.p.sq, float), r.fp, r.stats.get('ncontacts')


def features_job(job):
    pid, = job
    con = db()
    s, sq = load(con, pid)
    out = dict(pid=pid, s=s, img=image(s, sq), tilt=tilt_sig(sq))
    r = polish(s, sq)
    if r is None:
        out.update(sp=None)
        return out
    sp, q, fp, nc = r
    E, w = contact_graph(sp, q)
    h, core, ce = wl_hash(q, E, w)
    h0, _, _ = wl_hash(q, E, w, strip=False)
    out.update(sp=sp, fq=fp, nc=nc, wl=h, wl0=h0, core=core, ncon=len(E) + int(w.sum()))
    return out


# ---------------------------------------------------------------- matching
def sym(sq, s, k):
    x, y, a = sq[:, 0].copy(), sq[:, 1].copy(), sq[:, 2].copy()
    if k >= 4:
        x, a = s - x, -a
    for _ in range(k % 4):
        x, y, a = s - y, x, a + 90
    return np.stack([x, y, np.mod(a, 90)], 1)


def match(A, sA, B, sB):
    """rediscover.distance with per-square displacements kept.  -> rms, moved>0.1 fraction, k, perm (B index per A row)."""
    Af = A.copy(); Af[:, :2] *= sB / sA
    best = None
    for k in range(8):
        Bk = sym(B, sB, k)
        dx = Af[:, None, 0] - Bk[None, :, 0]; dy = Af[:, None, 1] - Bk[None, :, 1]
        da = np.mod(Af[:, None, 2] - Bk[None, :, 2] + 45, 90) - 45
        Cm = dx * dx + dy * dy + (0.5 * math.sqrt(2) * np.radians(da)) ** 2
        r, c = linear_sum_assignment(Cm)
        rms = math.sqrt(Cm[r, c].mean())
        if best is None or rms < best[0]:
            best = (rms, np.sqrt(Cm[r, c]), k, c)
    rms, d, k, c = best
    return rms, float((d > 0.1).mean()), k, c


def img_dist(a, b, n=N):
    from pk import fingerprint as F
    return F.dist(a, F.syms(b), n)


def pair_job(job):
    """(i, j, A, sA, B, sB, imA, imB, relabel) -> dict"""
    i, j, sA, A, sB, B, imA, imB, rl = job
    rms, mv, k, c = match(A, sA, B, sB)
    out = dict(i=i, j=j, md=rms, mv=mv, img=img_dist(imA, imB))
    if rl:
        Bk = sym(B, sB, k)[c]                 # B's squares in A's label order, in A's frame
        Bk[:, 2] = A[:, 2] - (np.mod(A[:, 2] - Bk[:, 2] + 45, 90) - 45)       # same angle representative as A
        r = polish(sB, Bk)
        out['fq_rl'] = r[2] if r else None
    return out


_P = None


def _init_reps(P):
    global _P
    _P = P


def rep_row_job(job):
    i, js = job
    P = _P
    sA, A = P[i]
    res = []
    for j in js:
        sB, B = P[j]
        rms, mv, _, _ = match(A, sA, B, sB)
        res.append((j, rms, mv))
    return i, res


# ---------------------------------------------------------------- stages
def minima():
    con = db()
    M = con.execute('SELECT key, CAST(S AS REAL), cls, pid FROM minimum WHERE n=? ORDER BY CAST(S AS REAL)', (N,)).fetchall()
    cp = defaultdict(set)
    for sk, mk in con.execute('SELECT side_key, min_key FROM copy WHERE n=?', (N,)):
        cp[sk].add(mk)
    return M, cp


def stage_features(a):
    os.makedirs(OUT, exist_ok=True)
    rng = random.Random(0)
    M, cp = minima()
    con = db()
    by_min = defaultdict(list)       # min_key -> [(pid, s, status)]
    for pid, s, st in con.execute("SELECT id, s, status FROM packing WHERE n=? AND status IN "
                                  "('polished','certified','weak')", (N,)):
        ks = cp.get(side_key(s))
        if ks and len(ks) == 1:
            by_min[next(iter(ks))].append((pid, s, st))
    collide = {sk: v for sk, v in cp.items() if len(v) > 1}
    sel, rep = [], {}
    for key, S, cls, pid in M:
        L = by_min.get(key, [])
        if pid is None or pid not in {x[0] for x in L}:
            pid = min(L, key=lambda x: abs(x[1] - S))[0] if L else None
        if pid is None:
            r = con.execute("SELECT id FROM packing WHERE n=? AND status!='infeasible' AND s BETWEEN ? AND ? "
                            "ORDER BY abs(s-?) LIMIT 1", (N, S - 1e-8, S + 1e-8, S)).fetchone()
            pid = r[0] if r else None
        rep[key] = pid
        others = [x for x in L if x[0] != pid]
        # prefer unconverged copies (side off by > 1e-9: the hard "same" cases), then random
        far = [x for x in others if abs(x[1] - S) > 1e-9]; near = [x for x in others if abs(x[1] - S) <= 1e-9]
        rng.shuffle(far); rng.shuffle(near)
        pick = (far[:COPIES // 2] + near)[:COPIES]
        if len(pick) < COPIES:
            pick += far[COPIES // 2:COPIES // 2 + COPIES - len(pick)]
        if pid is not None:
            sel.append((pid, key, 'rep'))
        sel += [(x[0], key, 'copy') for x in pick]
    print(f'{len(M)} minima, {sum(len(v) for v in by_min.values())} labelled rows, side-key collisions {len(collide)}; '
          f'sampled {len(sel)} packings', flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(PROCS) as ex:
        F = list(ex.map(features_job, [(p,) for p, _, _ in sel], chunksize=16))
    for f, (p, key, role) in zip(F, sel):
        f['key'] = key; f['role'] = role
    S = {k: S for k, S, _, _ in M}
    cls = {k: c for k, _, c, _ in M}
    pickle.dump(dict(F=F, rep=rep, S=S, cls=cls, collide={k: sorted(v) for k, v in collide.items()},
                     nlabelled={k: len(v) for k, v in by_min.items()}), open(f'{OUT}/features.pkl', 'wb'))
    print(f'features: {len(F)} in {time.time() - t0:.0f} s -> {OUT}/features.pkl')


def coords_of(pids):
    con = db()
    return {p: load(con, p) for p in pids}


def stage_pairs(a):
    D = pickle.load(open(f'{OUT}/features.pkl', 'rb'))
    F = [f for f in D['F'] if f.get('sp') is not None]
    rng = random.Random(1)
    by = defaultdict(list)
    for x, f in enumerate(F):
        by[f['key']].append(x)
    pairs = []
    for key, L in by.items():                     # all within-minimum pairs
        for u in range(len(L)):
            for v in range(u + 1, len(L)):
                pairs.append((L[u], L[v], 'same'))
    order = sorted(range(len(F)), key=lambda x: F[x]['s'])
    for p, x in enumerate(order):                 # hard negatives: the 3 nearest-in-side packings of another minimum
        got = 0
        for q in list(range(p + 1, min(p + 40, len(order)))):
            y = order[q]
            if F[y]['key'] != F[x]['key']:
                pairs.append((x, y, 'near')); got += 1
                if got == 3:
                    break
    reps = [x for x, f in enumerate(F) if f['role'] == 'rep']
    for _ in range(6000):                         # random different-minimum pairs
        x, y = rng.sample(reps, 2)
        pairs.append((x, y, 'rand'))
    C = coords_of({F[x]['pid'] for p in pairs for x in p[:2]})
    jobs = [(x, y, *C[F[x]['pid']], *C[F[y]['pid']], F[x]['img'], F[y]['img'], kind != 'rand') for x, y, kind in pairs]
    t0 = time.time()
    with ProcessPoolExecutor(PROCS) as ex:
        R = list(ex.map(pair_job, jobs, chunksize=32))
    for r, (x, y, kind) in zip(R, pairs):
        r['kind'] = kind
    pickle.dump(R, open(f'{OUT}/pairs.pkl', 'wb'))
    print(f'pairs: {len(R)} ({Counter(k for _, _, k in pairs)}) in {time.time() - t0:.0f} s')


def stage_reps(a):
    D = pickle.load(open(f'{OUT}/features.pkl', 'rb'))
    keys = sorted(D['rep'], key=lambda k: D['S'][k])
    keys = [k for k in keys if D['rep'][k] is not None]
    C = coords_of([D['rep'][k] for k in keys])
    P = [C[D['rep'][k]] for k in keys]
    K = len(keys)
    from pk import fingerprint as Fp
    ims = np.stack([image(*p) for p in P])
    IMG = np.zeros((K, K))
    for i in range(K):
        IMG[i] = Fp.nearest(ims[i], ims, N)
    IMG = np.minimum(IMG, IMG.T)
    MD = np.zeros((K, K)); MV = np.zeros((K, K))
    t0 = time.time()
    chunks = [(i, list(range(i + 1, K))) for i in range(K - 1)]
    with ProcessPoolExecutor(PROCS, initializer=_init_reps, initargs=(P,)) as ex:
        for i, res in ex.map(rep_row_job, chunks, chunksize=4):
            for j, rms, mv in res:
                MD[i, j] = MD[j, i] = rms; MV[i, j] = MV[j, i] = mv
    np.savez_compressed(f'{OUT}/reps.npz', keys=np.array(keys), S=np.array([D['S'][k] for k in keys]),
                        cls=np.array([D['cls'][k] for k in keys]), IMG=IMG, MD=MD, MV=MV)
    print(f'reps: {K} x {K} in {time.time() - t0:.0f} s')


def stage_edges(a):
    """Parent -> child transitions between known minima, from explorer archives (side -> minimum via the copy table)."""
    M, cp = minima()
    Ss = np.array([S for _, S, _, _ in M]); Ks = [k for k, _, _, _ in M]

    def lab(s):
        ks = cp.get(side_key(s))
        if ks and len(ks) == 1:
            return next(iter(ks))
        i = np.searchsorted(Ss, s)
        for j in (i - 1, i):
            if 0 <= j < len(Ss) and abs(Ss[j] - s) < 2e-9:
                return Ks[j]
        return None
    E = Counter(); runs = 0; tot = 0; lab_par = 0
    for d in sorted(glob.glob(f'{PK}/runs/ex*') + glob.glob(f'{PK}/runs/census_v*/*_r*') + glob.glob(f'{PK}/runs/c110_h*')):
        ap = os.path.join(d, 'archive.jsonl')
        if not os.path.exists(ap):
            continue
        A = {}
        for l in open(ap):
            e = json.loads(l)
            A[e['i']] = e
        runs += 1
        for e in A.values():
            p = e.get('parent', -1)
            if p is None or p < 0 or p not in A:
                continue
            tot += 1
            a_, b_ = lab(A[p]['s']), lab(e['s'])
            if a_ and b_:
                lab_par += 1
                if a_ != b_:
                    E[(a_, b_)] += 1
    json.dump(dict(runs=runs, edges=tot, labelled=lab_par, trans=[[a_, b_, c] for (a_, b_), c in E.items()]),
              open(f'{OUT}/edges.json', 'w'))
    print(f'edges: {runs} archives, {tot} parent links, {lab_par} with both ends a known minimum, '
          f'{len(E)} distinct transitions between minima ({sum(E.values())} links)')


# ---------------------------------------------------------------- analysis
def single_linkage(Dm, thr):
    K = len(Dm)
    par = list(range(K))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    I, J = np.nonzero(np.triu(Dm <= thr, 1))
    for i, j in zip(I, J):
        a_, b_ = f(i), f(j)
        if a_ != b_:
            par[a_] = b_
    roots = [f(x) for x in range(K)]
    m = {r: c for c, r in enumerate(dict.fromkeys(roots))}
    return np.array([m[r] for r in roots])


def leader(Dm, order, r):
    """Greedy leader clustering: visit in `order` (deepest first); join the nearest leader within r, else lead."""
    lab = np.full(len(Dm), -1); L = []
    for i in order:
        if L:
            d = Dm[i, L]; j = int(np.argmin(d))
            if d[j] <= r:
                lab[i] = lab[L[j]]; continue
        lab[i] = len(L); L.append(i)
    return lab


def ari(a, b):
    from math import comb
    a = np.asarray(a); b = np.asarray(b)
    ct = Counter(zip(a.tolist(), b.tolist()))
    sij = sum(comb(v, 2) for v in ct.values())
    sa = sum(comb(v, 2) for v in Counter(a.tolist()).values()); sb = sum(comb(v, 2) for v in Counter(b.tolist()).values())
    e = sa * sb / comb(len(a), 2); m = (sa + sb) / 2
    return (sij - e) / (m - e) if m != e else 1.0


def rates(R, F, test, kinds=('same', 'near', 'rand')):
    out = {}
    for k in kinds:
        L = [r for r in R if r['kind'] == k]
        v = [test(r, F[r['i']], F[r['j']]) for r in L]
        v = [x for x in v if x is not None]
        out[k] = (float(np.mean(v)) if v else float('nan'), len(v))
    return out


def stage_analyze(a):
    D = pickle.load(open(f'{OUT}/features.pkl', 'rb'))
    F = [f for f in D['F'] if f.get('sp') is not None]
    R = pickle.load(open(f'{OUT}/pairs.pkl', 'rb'))
    Z = np.load(f'{OUT}/reps.npz')
    lines = []
    P = lambda *x: (lines.append(' '.join(map(str, x))), print(*x))
    J = {}
    # ---- task 2: same-minimum detection, pair rates ("merged" = test says same)
    P('== Task 2: same-minimum detection (pairs of sampled packings; same = same exact minimum) ==')
    P(f"packings {len(F)} ({Counter(f['role'] for f in F)}), pairs {Counter(r['kind'] for r in R)}; "
      f"side-key collisions in copy table: {len(D['collide'])}")
    for f in F:
        f['dS'] = f['s'] - D['S'][f['key']]
    dS = np.array([f['dS'] for f in F if f['role'] == 'copy'])
    P('copy side offset s - S*: ' + ', '.join(f'>{t:g}: {np.mean(dS > t):.3f}' for t in (1e-9, 1e-7, 1e-5, 1e-4)))
    dsp = np.array([f['sp'] - D['S'][f['key']] for f in F])
    P('after fq polish sp - S*: ' + ', '.join(f'>{t:g}: {np.mean(dsp > t):.3f}' for t in (1e-9, 1e-7, 1e-5, 1e-4)))
    tests = {}
    for t in (2e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        tests[f'|ds|<{t:g}'] = lambda r, A, B, t=t: abs(A['s'] - B['s']) < t
    for t in (2e-9, 1e-8, 1e-7, 1e-6):
        tests[f'|dsp|<{t:g} (fq-polished)'] = lambda r, A, B, t=t: abs(A['sp'] - B['sp']) < t
    for t in (0.002, 0.005, 0.01, 0.02):
        tests[f'img<{t:g}'] = lambda r, A, B, t=t: r['img'] < t
    for t in (0.02, 0.05, 0.1):
        tests[f'md<{t:g}'] = lambda r, A, B, t=t: r['md'] < t
    for t in (0.02, 0.05, 0.1):
        tests[f'mv<{t:g}'] = lambda r, A, B, t=t: r['mv'] < t
    tests['tilt equal'] = lambda r, A, B: A['tilt'] == B['tilt']
    tests['fq fp equal (stored labels)'] = lambda r, A, B: A['fq'] == B['fq']
    tests['fq fp equal (B relabelled to A)'] = lambda r, A, B: (A['fq'] == r['fq_rl']) if 'fq_rl' in r and r['fq_rl'] else None
    tests['wl equal (rattlers stripped)'] = lambda r, A, B: A['wl'] == B['wl']
    tests['wl equal (all contacts)'] = lambda r, A, B: A['wl0'] == B['wl0']
    for t in (1e-7, 1e-6, 1e-5):
        for e in (0.01, 0.02):
            tests[f'|ds|<{t:g} & img<{e:g}'] = lambda r, A, B, t=t, e=e: abs(A['s'] - B['s']) < t and r['img'] < e
    for t in (1e-8, 1e-7):
        tests[f'|dsp|<{t:g} & img<0.02'] = lambda r, A, B, t=t: abs(A['sp'] - B['sp']) < t and r['img'] < 0.02
        tests[f'|dsp|<{t:g} & wl equal'] = lambda r, A, B, t=t: abs(A['sp'] - B['sp']) < t and A['wl'] == B['wl']
    tests['|ds|<1e-6 & mv<0.05'] = lambda r, A, B: abs(A['s'] - B['s']) < 1e-6 and r['mv'] < 0.05
    P(f"{'test (says same)':38s} {'false-split (same)':>18s} {'false-merge (near-side diff)':>28s} {'false-merge (random diff)':>26s}")
    J['task2'] = {}
    for name, t in tests.items():
        o = rates(R, F, t)
        J['task2'][name] = o
        P(f"{name:38s} {1 - o['same'][0]:18.4f} {o['near'][0]:28.4f} {o['rand'][0]:26.4f}   (n {o['same'][1]}/{o['near'][1]}/{o['rand'][1]})")
    # side gaps between distinct minima
    S = np.sort(np.array(list(D['S'].values())))
    g = np.diff(S)
    P('gaps between consecutive exact minima S*: ' + ', '.join(f'<{t:g}: {(g < t).sum()}' for t in (1e-12, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)) + f' (of {len(g)})')
    # distance distributions
    for k in ('same', 'near', 'rand'):
        L = [r for r in R if r['kind'] == k]
        q = lambda key: ' '.join(f'{np.percentile([r[key] for r in L], p):.4f}' for p in (50, 90, 99, 100))
        P(f'{k:5s}  img p50/90/99/max {q("img")}   md {q("md")}   mv {q("mv")}')
    # ---- task 3: regions
    P('\n== Task 3: region structure among the 939 representatives ==')
    keys, Sv, cls, IMG, MD, MV = Z['keys'], Z['S'], Z['cls'], Z['IMG'], Z['MD'], Z['MV']
    K = len(keys)
    top = np.argsort(Sv)[:50]
    sub = single_linkage(IMG[np.ix_(top, top)], 0.02)
    rec = top[0]
    A_ = set(top[sub == sub[0]])
    iB = int(np.argmin(np.abs(Sv - 10.996887)))
    B_ = set(top[sub == sub[list(top).index(iB)]]) if iB in top else set()
    P(f'reconstructed A = {len(A_)} minima, B = {len(B_)} (landscape_depth: 11 / 39)')
    dSm = np.abs(Sv[:, None] - Sv[None])
    measures = {'img': IMG, 'md': MD, 'mv': MV}
    for lam in (10, 30, 100):
        measures[f'md+{lam}|dS|'] = MD + lam * dSm
    measures['max(img/0.02, md/0.05)'] = np.maximum(IMG / 0.02, MD / 0.05)
    E = json.load(open(f'{OUT}/edges.json'))
    idx = {k: i for i, k in enumerate(keys)}
    T = [(idx[a_], idx[b_], c) for a_, b_, c in E['trans'] if a_ in idx and b_ in idx]
    tw = np.array([c for _, _, c in T], float)
    P(f"transitions between minima (explorer archives): {len(T)} distinct, {int(tw.sum())} links")
    ti = np.array([t[0] for t in T]); tj = np.array([t[1] for t in T])
    for nm in ('img', 'md', 'mv'):
        d = measures[nm][ti, tj]
        allp = measures[nm][np.triu_indices(K, 1)]
        P(f'  {nm}: transition distance p50/p90/p99 (link-weighted) '
          + ' '.join(f'{np.percentile(np.repeat(d, tw.astype(int)), p):.4f}' for p in (50, 90, 99))
          + '  | all pairs p1/p10/p50 ' + ' '.join(f'{np.percentile(allp, p):.4f}' for p in (1, 10, 50)))
    # transition probability per pair vs distance bin (kernel)
    for nm in ('md', 'img'):
        Mx = measures[nm]; iu = np.triu_indices(K, 1); allp = Mx[iu]
        edges_ = np.quantile(allp, [0, .001, .005, .01, .02, .05, .1, .25, .5, 1])
        cnt = np.zeros((K, K)); cnt[ti, tj] += 1; cnt = np.maximum(cnt, cnt.T)[iu] > 0
        row = []
        for b in range(len(edges_) - 1):
            m = (allp >= edges_[b]) & (allp <= edges_[b + 1])
            row.append(f'{edges_[b + 1]:.3f}:{cnt[m].mean():.3f}')
        P(f'  P(pair ever transitions | {nm} <= edge): ' + ' '.join(row))
    rows = []
    P(f"{'method':8s} {'measure':22s} {'thr':>6s} {'#cl':>5s} {'neff':>6s} {'max':>4s} {'single':>6s} {'A one':>5s} {'A|B':>5s} "
      f"{'ARIimg.02':>9s} {'within':>7s} {'null':>6s} {'lift':>5s} {'Sp(S,cl)':>8s}")
    ref = single_linkage(IMG, 0.02)
    grids = {'img': (0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04), 'md': (0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05),
             'mv': (0.005, 0.01, 0.02, 0.03)}
    order = np.argsort(Sv)
    from scipy.stats import spearmanr
    for meth in ('single', 'leader'):
        for nm, Mx in measures.items():
            base = nm.split('+')[0]
            for t in grids.get(base, (0.5, 0.75, 1, 1.5, 2)):
                lab = single_linkage(Mx, t) if meth == 'single' else leader(Mx, order, t)
                sz = Counter(lab.tolist())
                p_ = np.array(list(sz.values())) / K
                neff = float(np.exp(-(p_ * np.log(p_)).sum()))
                Aone = len({lab[i] for i in A_}) == 1
                split = not ({lab[i] for i in A_} & {lab[i] for i in B_})
                within = tw[lab[ti] == lab[tj]].sum() / tw.sum()
                pa = Counter(); pb = Counter()
                for (i, j, c) in T:
                    pa[lab[i]] += c; pb[lab[j]] += c
                null = sum(pa[c] * pb[c] for c in pa) / tw.sum() ** 2
                cm = np.array([np.mean(Sv[lab == lab[i]]) for i in range(K)])
                sp = spearmanr(Sv, cm)[0] if len(sz) > 1 else float('nan')
                r_ = dict(method=meth, measure=nm, thr=t, ncl=len(sz), neff=neff, max=max(sz.values()),
                          single=sum(v == 1 for v in sz.values()), A_one=Aone, AB_split=split, ari=ari(ref, lab),
                          within=within, null=null, lift=within / max(null, 1e-12), sp=sp)
                rows.append(r_)
                P(f"{meth:8s} {nm:22s} {t:6g} {len(sz):5d} {neff:6.1f} {r_['max']:4d} {r_['single']:6d} {str(Aone)[0]:>5s} "
                  f"{str(split)[0]:>5s} {r_['ari']:9.3f} {within:7.3f} {null:6.3f} {r_['lift']:5.1f} {sp:8.3f}")
    J['task3'] = rows
    json.dump(J, open(f'{OUT}/summary.json', 'w'), indent=1, default=float)
    open(f'{OUT}/report.txt', 'w').write('\n'.join(lines) + '\n')
    print(f'-> {OUT}/report.txt, summary.json')


# ---------------------------------------------------------------- other n
def othern_job(job):
    i, j, sA, A, sB, B = job
    rms, mv, _, _ = match(A, sA, B, sB)
    return i, j, rms, mv


def _md_row(job):
    i, js = job
    sA, A = _P[i]
    return i, [(j, *match(A, sA, _P[j][1], _P[j][0])[:2]) for j in js]


def stage_othern(a):
    """Same protocol without exact minima (side key = minimum proxy): one packing per side key below the grid, a random
    `--sample` of them; md / mv / img matrices, leader + single-linkage cluster counts, and the store's parent -> child
    step size (parent links, both sub-k, different side key)."""
    from pk import fingerprint as Fp
    con = db()
    res = {}
    for n in a.n:
        rng = random.Random(2)
        k = math.ceil(math.sqrt(n) - 1e-12)
        rows = con.execute('SELECT id, s, parent FROM packing WHERE n=? AND s<? ORDER BY s', (n, k)).fetchall()
        seen, U = set(), []
        for pid, s, par in rows:
            if side_key(s) not in seen:
                seen.add(side_key(s)); U.append((pid, s, par))
        nkeys = len(U)
        U = sorted(rng.sample(U, min(a.sample, len(U))), key=lambda u: u[1])
        C = [load(con, p) for p, _, _ in U]
        K = len(U)
        ims = np.stack([image(*c) for c in C]).astype(np.float64)
        IMG = np.zeros((K, K))
        for i in range(K):
            IMG[i] = Fp.nearest(ims[i], ims, n)
        IMG = np.minimum(IMG, IMG.T)
        t0 = time.time()
        MD = np.zeros((K, K)); MV = np.zeros((K, K))
        with ProcessPoolExecutor(PROCS, initializer=_init_reps, initargs=(C,)) as ex:
            for i, R in ex.map(_md_row, [(i, list(range(i + 1, K))) for i in range(K - 1)], chunksize=2):
                for j, rms, mv in R:
                    MD[i, j] = MD[j, i] = rms; MV[i, j] = MV[j, i] = mv
        # parent -> child steps
        side = {pid: s for pid, s, _ in rows}
        pc = [(p, par) for p, s, par in rows if par in side and side_key(side[par]) != side_key(s)]
        rng.shuffle(pc); pc = pc[:300]
        pcm = []
        for p, par in pc:
            sA, A = load(con, par); sB, B = load(con, p)
            pcm.append(match(A, sA, B, sB)[0])
        Dn = MD + np.diag(np.full(K, np.inf)); Dimg = IMG + np.diag(np.full(K, np.inf))
        iu = np.triu_indices(K, 1)
        pct = lambda v: [round(float(np.percentile(v, q)), 4) for q in (10, 50, 90)]
        out = dict(n=n, k=k, side_keys=nkeys, sample=K, md_secs=round(time.time() - t0),
                   parent_child_md=pct(pcm), nn_md=pct(Dn.min(1)), all_md=pct(MD[iu]), nn_img=pct(Dimg.min(1)),
                   all_img=pct(IMG[iu]), nn_mv=pct((MV + np.diag(np.full(K, np.inf))).min(1)))
        order = np.arange(K)
        out['leader_md'] = {r: int(len(set(leader(MD, order, r).tolist()))) for r in (0.02, 0.03, 0.04, 0.05, 0.075, 0.1)}
        out['leader_img'] = {r: int(len(set(leader(IMG, order, r).tolist()))) for r in (0.01, 0.02, 0.03, 0.04)}
        out['single_md'] = {r: max(Counter(single_linkage(MD, r).tolist()).values()) for r in (0.02, 0.03, 0.04, 0.05)}
        out['corr_img_md_spearman'] = float(__import__('scipy.stats').stats.spearmanr(IMG[iu], MD[iu])[0])
        res[n] = out
        print(json.dumps(out), flush=True)
        np.savez_compressed(f'{OUT}/othern_{n}_s{a.sample}.npz', S=np.array([u[1] for u in U]), IMG=IMG, MD=MD, MV=MV)
    json.dump(res, open(f'{OUT}/othern_s{a.sample}_' + '_'.join(map(str, a.n)) + '.json', 'w'), indent=1)


def main():
    global PROCS
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['features', 'pairs', 'reps', 'edges', 'analyze', 'othern', 'all'])
    ap.add_argument('--n', type=int, nargs='+', default=[88, 110, 237])
    ap.add_argument('--procs', type=int, default=PROCS); ap.add_argument('--sample', type=int, default=400)
    a = ap.parse_args()
    PROCS = a.procs
    os.makedirs(OUT, exist_ok=True)
    st = dict(features=stage_features, pairs=stage_pairs, reps=stage_reps, edges=stage_edges, analyze=stage_analyze,
              othern=stage_othern)
    for s in (['features', 'pairs', 'reps', 'edges', 'analyze', 'othern'] if a.stage == 'all' else [a.stage]):
        st[s](a)


if __name__ == '__main__':
    main()
