#!/usr/bin/env python3
"""V-band constructor (10-13; Ostashev's README step 1, #575): two staircase arms of tilted squares + exact upright fill.

Records near k(k+1) (153-156, 175-182, 236-241, 268-273, 302-308) are aligned corner blocks plus a band of tilted
squares from one wall to another with one bend.  Measured on 241 (pre-#575, store 347): each arm is a staircase block,
C = 5-7 parallel columns of tilted squares along the arm (column i at lattice coordinate u = u0 + i + du_i, squares at
v = v0_i, v0_i + 1, ...; L_i = 7-10 squares; ends tapered by the v0_i shear).  The exact fill (mondrian.exact_fill)
reproduces n from a record's own tilted squares at 241 (both), 211 (both), 156, 176: the fill is expressive enough for
this family, so the band is what the constructor has to get right.

Arm k: angle th_k, columns (du_i, v0_i, L_i) in the arm's own frame u = (cos, sin), v = (-sin, cos), anchored to a box
corner (arm A: top-left, arm B: bottom-right) so its position keeps its meaning while the side is bisected.  Squares of
an arm that overlap an earlier one of the same arm, squares outside the box, and arm-B squares overlapping arm A (the
elbow) are dropped.

Pipeline per parameter set: S_min = smallest side (bisection) at which band + exact fill holds n squares, then trim to n,
fq relax ops (relax.OPS) -> side.  Search = CEM over the continuous parameters (angles, anchors, du_i, v0_i) with the
column lengths fixed per run (`--lengths`: +-1 mutations), started from a fitted packing (`--x0 ID`) or a JSON file.

  python3 vband.py fit ID [--out X.json]            arms read off store packing ID; rebuild + fill check
  python3 vband.py search --n N --x0 ID|X.json --out DIR [--sigma 1] [--pop 32] [--iters 8] [--procs 16]
"""
from __future__ import annotations
import argparse, copy, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ANCHOR = {'A': (0.0, 1.0), 'B': (1.0, 0.0)}        # corner (fraction of S) each arm's coordinates are relative to


def arm_squares(arm, S):
    """Absolute (x, y, th) of an arm: dict(k, th, ax, ay (offset from its anchor corner), cols=[[du, v0, L], ...])."""
    t = math.radians(arm['th'])
    c, s = math.cos(t), math.sin(t)
    cx, cy = ANCHOR[arm['k']][0] * S + arm['ax'], ANCHOR[arm['k']][1] * S + arm['ay']
    out = []
    for i, (du, v0, L) in enumerate(arm['cols']):
        u = i + du
        for j in range(int(L)):
            v = v0 + j
            out.append((cx + u * c - v * s, cy + u * s + v * c, arm['th'] % 90.0))
    return np.array(out, float).reshape(-1, 3)


def _pen(a, b):
    import gen
    return gen.pen(tuple(a), tuple(b))


def _clean(T, S, against=None, tol=1e-9):
    """Keep squares of T inside the box, not overlapping `against`, nor an already kept square of T (greedy, in order)."""
    from gen import wall_viol
    keep = []
    for q in T:
        if wall_viol(tuple(q), S) > tol:
            continue
        near = lambda P: [p for p in P if (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 < 2.0]
        if against is not None and any(_pen(q, p) > tol for p in near(against)):
            continue
        if any(_pen(q, p) > tol for p in near(keep)):
            continue
        keep.append(q)
    return np.array(keep, float).reshape(-1, 3)


def build(V, S):
    """Tilted squares of the band V = dict(A=arm, B=arm) at side S."""
    A = _clean(arm_squares(V['A'], S), S)
    B = _clean(arm_squares(V['B'], S), S, against=A)
    return np.r_[A, B]


def fill_count(T, S, tl=20):
    import mondrian
    L = mondrian.exact_fill(T, S, time_limit=tl)
    return len(L), L


def smin(V, n, lo, hi, steps=9, tl=20):
    """Smallest side in [lo, hi] at which band + fill holds n: (S, L) or (None, best count at hi)."""
    c, L = fill_count(build(V, hi), hi, tl)
    if c < n:
        return None, c
    best = (hi, L)
    for _ in range(steps):
        mid = (lo + hi) / 2
        c, Lm = fill_count(build(V, mid), mid, tl)
        if c >= n:
            hi, best = mid, (mid, Lm)
        else:
            lo = mid
    return best


def write(path, n, S, sq):
    with open(path, 'w') as f:
        f.write(f'{n} {float(S)!r}\n' + ''.join(f'{float(a)!r} {float(b)!r} {float(c)!r}\n' for a, b, c in sq))


def finish(n, S, L, out, ops=('polish', 'soft15c'), seeds=1, explore_min=0.0, seed=0):
    """Trim to n, fq relax ops, then (explore_min > 0) the explorer from the relaxed layouts for that many minutes on
    one core (arrangement-level seeds + explorer recovered Xu's 126 in < 4 min, PACKER 10-09).  (side, file)."""
    import mondrian, relax, subprocess
    X = mondrian.trim(np.asarray(L, float), n)
    inp = out + '.in.txt'
    write(inp, n, S, X)
    best, outs = (math.inf, None), []
    for op in ops:
        for sd in range(1, (1 if op == 'polish' else seeds) + 1):
            r = relax.job((inp, f'{out}.{op}{sd}.txt', op, sd))
            if r.get('ok'):
                outs.append(r['out'])
                if r['s'] < best[0]:
                    best = (r['s'], r['out'])
    if explore_min > 0 and outs:
        d = out + '.ex'
        subprocess.run([sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(n), '--starts', *outs, '--procs', '1',
                        '--minutes', str(explore_min), '--seed', str(seed), '--out', d],
                       cwd=HERE, capture_output=True, text=True)
        try:
            E = json.load(open(os.path.join(d, 'state.json')))['E']
            b = min(E, key=lambda e: e['s'])
            if b['s'] < best[0]:
                p = b['path'] if os.path.isabs(b['path']) else os.path.join(HERE, b['path'])
                best = (b['s'], p)
        except Exception:
            pass
    return best


# ---------------------------------------------------------------- fit (controls, seeds)

def fit(sq, S):
    """Arms read off a packing: the two largest tilted bands (pk.bands, >= 10 squares) if there are two, else the one
    band split by 2-means on the angle; angle outliers (> 3 deg off the arm's median) dropped.  Per arm: columns =
    squares grouped by round(u - u0) in the arm's frame; each column's du (mean fractional u), v0 (lowest v) and L
    (span: gaps inside a column are closed; the relax reopens them).  Arm A = the arm reaching higher."""
    from pk.bands import bands
    sq = np.asarray(sq, float)
    _, comps = bands(sq, mn=10)

    def keep(F):
        a = F[:, 2] % 90.0
        dev = np.abs((a - np.median(a) + 45) % 90 - 45)
        return F[dev <= 3.0]
    if len(comps) >= 2:
        groups = [keep(sq[comps[0]]), keep(sq[comps[1]])]
    else:
        T = sq[comps[0]]
        a = T[:, 2] % 90.0
        m = np.median(a)
        a = np.where(a - m > 45, a - 90, np.where(m - a > 45, a + 90, a))
        c0, c1 = np.percentile(a, 10), np.percentile(a, 90)
        for _ in range(30):
            lab = np.abs(a - c0) > np.abs(a - c1)
            c0, c1 = a[~lab].mean(), a[lab].mean()
        groups = [keep(T[~lab]), keep(T[lab])]
    groups.sort(key=lambda F: F[:, 1].max(), reverse=True)
    V = {}
    for k, F in zip('AB', groups):
        th = float(np.median(F[:, 2] % 90.0))
        t = math.radians(th)
        c, s = math.cos(t), math.sin(t)
        x, y = F[:, 0] - ANCHOR[k][0] * S, F[:, 1] - ANCHOR[k][1] * S      # relative to the anchor corner
        u, v = x * c + y * s, -x * s + y * c
        u0 = float(np.floor(u.min() + 0.5))
        col = np.round(u - u0).astype(int)
        cols = []
        for i in range(col.min(), col.max() + 1):
            m = col == i
            if not m.any():
                cols.append([0.0, 0.0, 0]); continue
            cols.append([float(np.mean(u[m] - u0 - i)), float(v[m].min()), int(round(v[m].max() - v[m].min())) + 1])
        u_start = u0 + col.min()                                 # frame origin at the first column
        V[k] = dict(k=k, th=th, ax=u_start * c, ay=u_start * s, cols=cols)
    return V


# ---------------------------------------------------------------- search

def flat(V):
    """Continuous parameters: per arm th, ax, ay, then du_i, v0_i."""
    out = []
    for k in 'AB':
        a = V[k]
        out += [a['th'], a['ax'], a['ay']] + [c[0] for c in a['cols']] + [c[1] for c in a['cols']]
    return np.array(out, float)


def unflat(V, x):
    W = copy.deepcopy(V)
    i = 0
    for k in 'AB':
        a = W[k]
        C = len(a['cols'])
        a['th'], a['ax'], a['ay'] = (float(z) for z in x[i:i + 3]); i += 3
        for j in range(C):
            a['cols'][j][0] = float(x[i + j]); a['cols'][j][1] = float(x[i + C + j])
        i += 2 * C
    return W


def scales(V, sigma):
    out = []
    for k in 'AB':
        C = len(V[k]['cols'])
        out += [0.5, 0.15, 0.15] + [0.03] * C + [0.15] * C
    return np.array(out) * sigma


def _eval(t):
    V, n, lo, hi, tl = t
    try:
        r = smin(V, n, lo, hi, tl=tl)
    except Exception as e:
        return None, str(e)[:120]
    return r[0], (r[1].tolist() if r[0] is not None else int(r[1]))


def search(a):
    from pk import manifest, eventbench as E
    import scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'vband', args=vars(a), procs=a.procs)
    scratch.scope('vband')
    E.set_public(a.target == 'public', a.exclude_src)
    f = E.frontier_now(a.n)
    lo, hi = f - a.below_span, f + a.above_span
    V0 = load_x0(a.x0)
    if a.grid_cols:                                   # unit column spacing (measured du jitter makes columns overlap)
        for k in 'AB':
            for c in V0[k]['cols']:
                c[0] = 0.0
    if a.lengths:                                     # column lengths +-1 at random (fixed for the run)
        rng0 = np.random.default_rng(a.seed + 7)
        for k in 'AB':
            for c in V0[k]['cols']:
                c[2] = max(0, int(c[2] + rng0.integers(-1, 2)))
    json.dump(V0, open(os.path.join(a.out, 'x0.json'), 'w'), indent=1)
    mu, sd0 = flat(V0), scales(V0, a.sigma)
    sd = sd0.copy()
    rng = np.random.default_rng(a.seed)
    log = open(os.path.join(a.out, 'evals.jsonl'), 'a')
    elite_all = []
    with ProcessPoolExecutor(a.procs) as ex:
        for it in range(a.iters):
            X = ([mu] if it == 0 else []) + [mu + sd * rng.standard_normal(len(mu)) for _ in range(a.pop - (it == 0))]
            Vs = [unflat(V0, x) for x in X]
            res = list(ex.map(_eval, [(V, a.n, lo, hi, a.tl) for V in Vs]))
            scored = []
            for x, V, (s, L) in zip(X, Vs, res):
                log.write(json.dumps(dict(it=it, S=s, count=L if s is None else None, V=V)) + '\n')
                key = s if s is not None else hi + (a.n - (L if isinstance(L, int) else 0)) * 0.01
                scored.append((key, x, V, s, L))
            log.flush()
            scored.sort(key=lambda z: z[0])
            el = scored[:max(2, a.pop // 4)]
            mu = np.mean([z[1] for z in el], 0)
            sd = np.maximum(np.std([z[1] for z in el], 0), sd0 * 0.1)
            elite_all += [z for z in el if z[3] is not None]
            best = scored[0]
            print(f'it {it}: best S_min {best[3]}  (count at hi {best[4] if best[3] is None else n_ok(a.n)})  feasible '
                  f'{sum(z[3] is not None for z in scored)}/{len(scored)}  frontier {f:.10f}', flush=True)
    elite_all.sort(key=lambda z: z[0])
    fin, seen = [], []
    for key, x, V, s, L in elite_all:
        if len(fin) < a.finish and all(abs(s - t) > 1e-6 for t in seen):
            seen.append(s); fin.append((V, s, L))
    out = []
    with ProcessPoolExecutor(a.procs) as ex:
        futs = [ex.submit(finish, a.n, s, L, os.path.join(a.out, f'fin{k}'), explore_min=a.explore, seed=k)
                for k, (V, s, L) in enumerate(fin)]
        for k, (fu, (V, s, L)) in enumerate(zip(futs, fin)):
            side, path = fu.result()
            out.append(dict(k=k, S_min=s, side=side, file=path, gain=f - side, V=V))
            print(f'finish {k}: S_min {s:.6f} -> {side:.10f}  gain {f - side:+.2e}', flush=True)
    json.dump(dict(n=a.n, frontier=f, finished=out), open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    man.finish(summary=dict(n=a.n, frontier=f, best=min((o['side'] for o in out), default=None)))


def n_ok(n):
    return f'>= {n}'


def load_x0(x0):
    if x0.endswith('.json'):
        return json.load(open(x0))
    from pk.store import Store
    p = Store().get(int(x0))
    return fit(p.sq, float(p.s))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    f = sp.add_parser('fit'); f.add_argument('id', type=int); f.add_argument('--out')
    s = sp.add_parser('search'); s.add_argument('--n', type=int, required=True); s.add_argument('--out', required=True)
    s.add_argument('--x0', required=True, help='store id to fit, or a fit JSON')
    s.add_argument('--sigma', type=float, default=1.0); s.add_argument('--pop', type=int, default=32)
    s.add_argument('--iters', type=int, default=8); s.add_argument('--procs', type=int, default=16)
    s.add_argument('--finish', type=int, default=16); s.add_argument('--seed', type=int, default=1)
    s.add_argument('--tl', type=float, default=20); s.add_argument('--lengths', action='store_true',
                                                                   help='mutate column lengths by +-1 (per seed)')
    s.add_argument('--below-span', type=float, default=0.02); s.add_argument('--above-span', type=float, default=0.15)
    s.add_argument('--target', default='now'); s.add_argument('--exclude-src', nargs='*', default=[])
    s.add_argument('--grid-cols', action='store_true', help='du_i = 0 in the start (unit column spacing)')
    s.add_argument('--explore', type=float, default=0.0, help='explorer minutes (1 core) per finished layout')
    a = ap.parse_args()
    if a.cmd == 'fit':
        from pk.store import Store
        from pk.bands import bands
        p = Store().get(a.id)
        V = fit(p.sq, float(p.s))
        for k in 'AB':
            A = V[k]
            print(f"arm {k}: th {A['th']:.3f} anchor offset ({A['ax']:.3f}, {A['ay']:.3f})  columns (L): "
                  + ' '.join(str(c[2]) for c in A['cols']) + f"  = {sum(c[2] for c in A['cols'])}")
        if a.out:
            json.dump(V, open(a.out, 'w'), indent=1)
        for eps in (1e-6, 0.02, 0.05, 0.1):
            T = build(V, p.s + eps)
            c, L = fill_count(T, p.s + eps)
            print(f'  S + {eps:g}: {len(T)} tilted (band of the packing: {sum(bands(p.sq, mn=10)[0])}), tilted + fill = {c}'
                  f' (n = {p.n})')
    else:
        search(a)
