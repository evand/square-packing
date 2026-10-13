#!/usr/bin/env python3
"""Arm composer (10-13): real tilted arms cut from band packings, varied locally, recombined, exact upright fill.

Why (PACKER 10-13 V-band notes): a record's own tilted squares + mondrian.exact_fill give exactly n at 241, 211, 156, 176
(the fill is expressive for the band family), but idealised arms (vband.py) lose ~20 squares: small boundary deviations
cost whole upright rows.  So keep the arms' real coordinates and change them only where an operator acts.

Library: every store packing in the band family (pk.bands: a tilted band of >= MIN_BAND squares), n in [nmin, nmax],
best few per n; its tilted squares split into arms (vband.split_arms: two components, or two angle clusters, or one
arm); each arm in each of the 8 box symmetries, stored relative to the box corner nearest its centroid (it moves with
that corner when the side changes).

Variation operators (on an arm's real squares, in its own column frame u = (cos, sin), v = (-sin, cos) at its median
angle): extend a column by one square at either end (one unit step along that end square's own v), trim a column end,
add a column beside an edge column (its copy shifted one unit along u), drop an edge column, shift the arm.  Ends are
"promising" (weight W_NEAR) when within NEAR of the other arm or a wall.

Composition for a target n at side S: arm A + arm B anchored at different corners, B squares overlapping A dropped,
squares out of the box dropped, exact fill; S_min by bisection; the best layouts are trimmed to n and polished (fq) and
optionally explored.  Search = random pairs x random ops within a time budget per n.

  python3 arms.py lib --out runs/arms/lib.pkl [--nmin 80 --nmax 324 --per-n 4]
  python3 arms.py run --lib runs/arms/lib.pkl --ns 236-241 --out DIR [--budget 300] [--procs 16] [--exclude-n]
"""
from __future__ import annotations
import argparse, json, math, os, pickle, random, sys, time
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
MIN_BAND = 20
NEAR = 1.6
W_NEAR = 4.0
CORNERS = [(0, 0), (1, 0), (0, 1), (1, 1)]
TOL = 1e-6             # overlap / wall tolerance of the layout (exact_fill's own; touching record arms are ~1e-9 apart)


# ---------------------------------------------------------------- geometry

def sym(sq, S, k):
    """k in 0..7: bit 0 mirror x, bit 1 mirror y, bit 2 transpose (angles transformed so squares stay squares)."""
    X = np.asarray(sq, float).copy()
    if k & 4:
        X[:, [0, 1]] = X[:, [1, 0]]; X[:, 2] = (90.0 - X[:, 2]) % 90.0
    if k & 1:
        X[:, 0] = S - X[:, 0]; X[:, 2] = (-X[:, 2]) % 90.0
    if k & 2:
        X[:, 1] = S - X[:, 1]; X[:, 2] = (-X[:, 2]) % 90.0
    return X


def frame(arm):
    """Median angle, unit vectors u, v, column index per square (rounded u in units), v coordinate."""
    sq = arm['sq']
    a = sq[:, 2] % 90.0
    m = np.median(a)
    a = np.where(a - m > 45, a - 90, np.where(m - a > 45, a + 90, a))
    t = math.radians(float(np.median(a)))
    u, v = np.array([math.cos(t), math.sin(t)]), np.array([-math.sin(t), math.cos(t)])
    pu, pv = sq[:, :2] @ u, sq[:, :2] @ v
    col = np.round(pu - pu.min()).astype(int)
    return t, u, v, col, pv


def place(arm, S):
    """Absolute squares of an arm in a box of side S.  arm['fixed']: keep the source coordinates (the box grows / shrinks
    at the top and right walls only), used when both arms come from one packing so their elbow junction stays intact;
    otherwise the arm moves with its anchor corner."""
    cx, cy = arm['corner']
    X = arm['sq'].copy()
    ref = arm['S'] if arm.get('fixed') else S
    X[:, 0] += cx * ref; X[:, 1] += cy * ref
    return X


# ---------------------------------------------------------------- library

def split_arms(sq):
    """Tilted squares of a band packing -> list of arms (arrays): two components, else two angle clusters (if the two
    cluster medians differ by > 2 deg), else one."""
    from pk.bands import bands
    sq = np.asarray(sq, float)
    _, comps = bands(sq, mn=10)
    if not comps:
        return []
    if len(comps) >= 2:
        return [sq[comps[0]], sq[comps[1]]]
    T = sq[comps[0]]
    a = T[:, 2] % 90.0
    m = np.median(a)
    a = np.where(a - m > 45, a - 90, np.where(m - a > 45, a + 90, a))
    c0, c1 = np.percentile(a, 10), np.percentile(a, 90)
    for _ in range(30):
        lab = np.abs(a - c0) > np.abs(a - c1)
        if lab.all() or (~lab).all():
            break
        c0, c1 = a[~lab].mean(), a[lab].mean()
    if abs(c1 - c0) > 2.0 and 5 <= lab.sum() <= len(T) - 5:
        return [T[~lab], T[lab]]
    return [T]


def build_lib(a):
    from pk.store import Store
    from pk.bands import bands
    st = Store()
    lib = []
    for n in range(a.nmin, a.nmax + 1):
        rows = st.query(n=n, limit=40)
        seen, took = [], 0
        for r in rows:
            if took >= a.per_n:
                break
            if any(abs(r['s'] - s) < 1e-7 for s in seen):
                continue
            p = st.get(r['id'])
            b, _ = bands(p.sq, mn=10)
            if not b or sum(b) < MIN_BAND:
                continue
            seen.append(r['s']); took += 1
            S = float(p.s)
            for k in range(8):
                X = sym(p.sq, S, k)
                for ai, A in enumerate(split_arms(X)):
                    c = A[:, :2].mean(0) / S
                    corner = min(CORNERS, key=lambda q: (c[0] - q[0]) ** 2 + (c[1] - q[1]) ** 2)
                    rel = A.copy(); rel[:, 0] -= corner[0] * S; rel[:, 1] -= corner[1] * S
                    lib.append(dict(src=int(r['id']), n=n, S=S, sym=k, arm=ai, corner=corner, sq=rel, m=len(A)))
        print(f'{n}: {took} packings, lib {len(lib)}', flush=True)
    pickle.dump(lib, open(a.out, 'wb'))
    print(f'{len(lib)} arms from {len({x["src"] for x in lib})} packings -> {a.out}')


# ---------------------------------------------------------------- variation operators

def _ends(arm):
    """[(col, which, index of the end square)] for every column (which: -1 low end, +1 high end)."""
    _, _, _, col, pv = frame(arm)
    out = []
    for c in np.unique(col):
        I = np.nonzero(col == c)[0]
        out.append((int(c), -1, int(I[np.argmin(pv[I])])))
        out.append((int(c), +1, int(I[np.argmax(pv[I])])))
    return out


def mutate(arm, other_abs, S, rng, ops=1):
    """Apply `ops` random operators; promising ends (near the other arm or a wall) weighted W_NEAR."""
    A = dict(arm, sq=arm['sq'].copy())
    log = []
    for _ in range(ops):
        t, u, v, col, pv = frame(A)
        absq = place(A, S)
        ends = _ends(A)
        w = []
        for c, which, i in ends:
            p = absq[i, :2]
            near = (len(other_abs) and np.min(np.hypot(*(other_abs[:, :2] - p).T)) < NEAR + 1) or \
                min(p[0], p[1], S - p[0], S - p[1]) < NEAR
            w.append(W_NEAR if near else 1.0)
        w = np.array(w) / sum(w)
        op = rng.choice(['extend', 'extend', 'trim', 'trim', 'addcol', 'dropcol', 'shift'])
        if op in ('extend', 'trim'):
            c, which, i = ends[rng.choice(len(ends), p=w)]
            if op == 'extend':
                q = A['sq'][i].copy()
                tq = math.radians(q[2])
                vq = np.array([-math.sin(tq), math.cos(tq)])
                # step along the square's own v (or u: pick the axis closest to the arm's v)
                uq = np.array([math.cos(tq), math.sin(tq)])
                step = vq if abs(vq @ v) >= abs(uq @ v) else uq
                step = step * np.sign(step @ v) * which
                q[:2] += step
                A['sq'] = np.r_[A['sq'], q[None]]
            elif (col == c).sum() > 1 or rng.random() < 0.3:
                A['sq'] = np.delete(A['sq'], i, 0)
            log.append(f'{op}:c{c}{"+" if which > 0 else "-"}')
        elif op in ('addcol', 'dropcol'):
            edge = int(rng.choice([col.min(), col.max()]))
            I = np.nonzero(col == edge)[0]
            if op == 'addcol':
                sgn = -1 if edge == col.min() else 1
                N = A['sq'][I].copy(); N[:, :2] += sgn * u
                A['sq'] = np.r_[A['sq'], N]
            elif len(np.unique(col)) > 2:
                A['sq'] = np.delete(A['sq'], I, 0)
            log.append(f'{op}:c{edge}')
        else:
            d = rng.normal(0, 0.15, 2)
            A['sq'] = A['sq'].copy(); A['sq'][:, :2] += d
            log.append(f'shift:{d[0]:+.2f},{d[1]:+.2f}')
    return A, log


# ---------------------------------------------------------------- composition

def layout(A, B, S):
    import vband
    TA = vband._clean(place(A, S), S, tol=TOL)
    if B is None:
        return TA
    TB = vband._clean(place(B, S), S, against=TA, tol=TOL)
    return np.r_[TA, TB]


def smin(A, B, n, lo, hi, steps=9, tl=20):
    import vband
    c, L = vband.fill_count(layout(A, B, hi), hi, tl)
    if c < n:
        return None, c
    best = (hi, L)
    for _ in range(steps):
        mid = (lo + hi) / 2
        c, Lm = vband.fill_count(layout(A, B, mid), mid, tl)
        if c >= n:
            hi, best = mid, (mid, Lm)
        else:
            lo = mid
    return best


LIB = None


def _init(path):
    global LIB
    LIB = pickle.load(open(path, 'rb'))


def job(t):
    """One target n: random arm pairs x random variations within the budget.  Returns the best few layouts."""
    n, f, budget, seed, a_dn, a_ds, exclude_n, max_ops, tl = t
    rng = np.random.default_rng(seed)
    lo, hi = f - 0.03, f + 0.03
    pool = [x for x in LIB if abs(x['n'] - n) <= a_dn and abs(x['S'] - f) <= a_ds and not (exclude_n and x['n'] == n)]
    if not pool:
        return dict(n=n, f=f, tried=0, best=[])
    by_src = {}
    for i, x in enumerate(pool):
        by_src.setdefault((x['src'], x['sym']), []).append(i)
    t0, tried, best, evals = time.time(), 0, [], []
    while time.time() - t0 < budget:
        i = int(rng.integers(len(pool)))
        A = pool[i]
        r = rng.random()
        if r < 0.4:                                    # partner from the same packing + symmetry (its own other arm):
            P = [pool[j] for j in by_src[(A['src'], A['sym'])] if j != i]   # both fixed in the source frame
            B = dict(P[0], fixed=True) if P else None
            A = dict(A, fixed=True)
        elif r < 0.9:                                  # partner from another packing, a different corner
            cand = [x for x in pool if x['corner'] != A['corner']]
            B = cand[int(rng.integers(len(cand)))] if cand else None
        else:
            B = None                                   # one arm alone
        k = int(rng.integers(0, max_ops + 1))
        logs = []
        if k:
            other = place(B, f) if B is not None else np.zeros((0, 3))
            A, la = mutate(A, other, f, rng, ops=int(rng.integers(0, k + 1)))
            logs += ['A.' + s for s in la]
            if B is not None:
                B, lb = mutate(B, place(A, f), f, rng, ops=k - len(la))
                logs += ['B.' + s for s in lb]
        try:
            s, L = smin(A, B, n, lo, hi, tl=tl)
        except Exception:
            continue
        tried += 1
        rec = dict(A=(A['src'], A['n'], A['sym'], A['arm']), B=(B['src'], B['n'], B['sym'], B['arm']) if B is not None else None,
                   ops=logs, S_min=s, count=None if s is not None else L)
        evals.append(rec)
        if s is not None:
            best.append(dict(rec, L=np.asarray(L).tolist()))
            best.sort(key=lambda z: z['S_min'])
            best = best[:6]
    return dict(n=n, f=f, tried=tried, pool=len(pool), best=best, evals=evals)


def run(a):
    from pk import manifest, eventbench as E
    import scratch, vband
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'arms', args=vars(a), procs=a.procs)
    scratch.scope('arms')
    E.set_public(a.target == 'public', a.exclude_src)
    lo, hi = map(int, a.ns.split('-'))
    ns = [n for n in range(lo, hi + 1) if not a.only or n in {int(x) for x in a.only.split(',')}]
    tasks = []
    for n in ns:
        f = E.frontier_now(n)
        if f >= math.ceil(math.sqrt(n) - 1e-12) - 1e-9:
            continue
        for r in range(a.reps):
            tasks.append((n, f - 1e-9, a.budget, a.seed + 1000 * r + n, a.dn, a.ds, a.exclude_n, a.max_ops, a.tl))
    print(f'{len(tasks)} tasks over {len(ns)} n', flush=True)
    res = {}
    with ProcessPoolExecutor(a.procs, initializer=_init, initargs=(a.lib,)) as ex, \
            open(os.path.join(a.out, 'evals.jsonl'), 'w') as fe:
        for r in ex.map(job, tasks):
            for e in r.pop('evals', []):
                fe.write(json.dumps(dict(e, n=r['n'])) + '\n')
            R = res.setdefault(r['n'], dict(n=r['n'], f=r['f'], tried=0, best=[]))
            R['tried'] += r['tried']; R['pool'] = r.get('pool')
            R['best'] = sorted(R['best'] + r['best'], key=lambda z: z['S_min'])[:6]
            b = R['best'][0]['S_min'] if R['best'] else None
            print(f"{r['n']:4d} f {r['f']:.10f} pool {r.get('pool')} tried {r['tried']:5d}  best S_min "
                  f"{'%.10f' % b if b else None}{'  need %+.2e' % (b - r['f']) if b else ''}", flush=True)
    # finish: polish (+ explore) the best layouts per n
    fin = [(n, k, z) for n, R in res.items() for k, z in enumerate(R['best'][:a.finish])]
    out = []
    with ProcessPoolExecutor(a.procs) as ex:
        futs = [ex.submit(vband.finish, n, z['S_min'], z['L'], os.path.join(a.out, f'fin_{n}_{k}'),
                          explore_min=a.explore, seed=k) for n, k, z in fin]
        for (n, k, z), fu in zip(fin, futs):
            side, path = fu.result()
            f = res[n]['f']
            out.append(dict(n=n, k=k, S_min=z['S_min'], side=side, frontier=f, gain=f - side, file=path, A=z['A'],
                            B=z['B'], ops=z['ops']))
            print(f"finish {n} #{k}: S_min {z['S_min']:.8f} -> {side:.10f}  frontier {f:.10f}  gain {f - side:+.2e}"
                  f"{'  BELOW' if f - side > 1e-9 else ''}  {z['ops']}", flush=True)
    json.dump(dict(results=[{k: v for k, v in R.items() if k != 'best'} for R in res.values()], finished=out),
              open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    man.finish(summary=dict(ns=len(res), below=sum(o['gain'] > 1e-9 for o in out)))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    l = sp.add_parser('lib'); l.add_argument('--out', required=True)
    l.add_argument('--nmin', type=int, default=80); l.add_argument('--nmax', type=int, default=324)
    l.add_argument('--per-n', type=int, default=4)
    r = sp.add_parser('run'); r.add_argument('--lib', required=True); r.add_argument('--ns', required=True)
    r.add_argument('--out', required=True); r.add_argument('--only')
    r.add_argument('--budget', type=float, default=300, help='seconds per task'); r.add_argument('--reps', type=int, default=1)
    r.add_argument('--procs', type=int, default=16); r.add_argument('--seed', type=int, default=0)
    r.add_argument('--dn', type=int, default=25, help='arms from packings with |n - n_src| <= dn')
    r.add_argument('--ds', type=float, default=0.6, help='... and |S_src - frontier| <= ds')
    r.add_argument('--exclude-n', action='store_true', help='no arms from packings of the target n itself (control)')
    r.add_argument('--max-ops', type=int, default=3); r.add_argument('--tl', type=float, default=20)
    r.add_argument('--finish', type=int, default=3); r.add_argument('--explore', type=float, default=0.0)
    r.add_argument('--target', default='now'); r.add_argument('--exclude-src', nargs='*', default=[])
    a = ap.parse_args()
    build_lib(a) if a.cmd == 'lib' else run(a)
