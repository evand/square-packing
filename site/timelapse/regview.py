#!/usr/bin/env python3
"""Display regularisation of every drawing in out/history.json -> out/views/<key>.json (cached).

Per state: load the drawing (site JSON / packing file / store row / grid / derived = base minus squares), orient it,
regularise it, check it.  A *view* only: the side is never changed, and nothing here is a certificate.

  1. angles (degrees, folded to (-45, 45]): angles equal to 1e-6 deg become one exact float (float witnesses carry
     ~1e-12 noise); then near-axis squares (|theta| < 1 deg) are straightened to 0, all at once, if the translation LP
     (search/regularize/reg.py place) stays feasible at the side S + 1e-9 (else left as they are);
  2. orientation: the dihedral image that best matches the next state of the same n (Hungarian assignment on centres
     in units of the side, + an angle-mismatch term), walking back from the last state (which keeps its source
     orientation, as jlevy's atlas does); 'derived' and grid states inherit the rule too;
  3. placement: reg.pretty_place (largest feasible container symmetry, straight tilted chains, tilted clusters closed
     up, axis squares to the walls with gravity to the bottom left), else reg.place with plain gravity, at S + 1e-9;
     accepted if every pair and wall clears -1e-7 (reg.overlap_free), else the angle-cleaned source is drawn;
  4. the last state of an n whose side matches search/regularize/lists/best.json is drawn from that list (the
     certified regularised record, site/www/data/p/best-N.json) instead.

  regview.py [--procs 3] [--ns 1-324] [--force]
"""
import argparse, hashlib, json, math, os, sys, time
from multiprocessing import Pool
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(HERE, 'out')
VIEWS = os.path.join(OUT, 'views')
sys.path.insert(0, os.path.join(ROOT, 'search', 'regularize'))
sys.path.insert(0, os.path.join(ROOT, 'search', 'packer'))
VERSION = 'v3'
PAD = 1e-9              # display side slack for float witnesses


# ------------------------------------------------------------------------------------------------ loading

def fold(t):
    t = t % 90.0
    return t - 90.0 if t > 45.0 else t


def grid(n):
    k = math.ceil(math.sqrt(n) - 1e-12)
    sq = [(c + 0.5, r + 0.5, 0.0) for r in range(k) for c in range(k)][:n]     # rows from the bottom; holes top right
    return float(k), np.array(sq, float)


_DB = None


def store_row(pid):
    global _DB
    import sqlite3
    if _DB is None:
        _DB = sqlite3.connect(f"file:{os.path.join(ROOT, 'search', 'packer', 'runs', 'store.sqlite')}?mode=ro", uri=True)
    s, blob = _DB.execute('SELECT s, coords FROM packing WHERE id = ?', (pid,)).fetchone()
    return float(s), np.frombuffer(blob, '<f8').reshape(-1, 3).copy()


def load(draw, n):
    k = draw['kind']
    if k == 'grid':
        return grid(n)
    if k == 'store':
        return store_row(draw['id'])
    if k == 'file':
        p = os.path.join(ROOT, draw['path'])
        if p.endswith('.json'):
            d = json.load(open(p))
            return float(d['s']), np.array(d['squares'], float)
        from pk.packing import read
        P = read(p)
        return float(P.s), np.array(P.sq, float)
    if k == 'derived':
        S, A = load(draw['base'], draw['n'])
        return S, remove(S, A, draw['remove'])
    raise ValueError(k)


def contacts(S, A, tol=1e-6):
    """rough contact count per square (walls + neighbours within tol), float"""
    import reg
    n = len(A)
    CS = [reg.frame(t) for t in A[:, 2]]
    X, Y = A[:, 0], A[:, 1]
    c = np.zeros(n, int)
    for i in range(n):
        ex, ey = reg.support(*CS[i], 1, 0), reg.support(*CS[i], 0, 1)
        c[i] += (X[i] - ex < tol) + (S - X[i] - ex < tol) + (Y[i] - ey < tol) + (S - Y[i] - ey < tol)
    D2 = (X[:, None] - X[None, :]) ** 2 + (Y[:, None] - Y[None, :]) ** 2
    for i, j in zip(*np.nonzero(np.triu(D2 < 2.05, 1))):
        if reg.best_normal(i, j, X, Y, CS)[3] < tol:
            c[i] += 1; c[j] += 1
    return c


def remove(S, A, m):
    """drop m squares: fewest contacts first, then axis-parallel, then nearest the top-right corner"""
    c = contacts(S, A)
    key = [(c[i], abs(fold(A[i, 2])) > 1e-9, -(A[i, 0] + A[i, 1])) for i in range(len(A))]
    drop = set(sorted(range(len(A)), key=lambda i: key[i])[:m])
    return A[[i for i in range(len(A)) if i not in drop]]


# ------------------------------------------------------------------------------------------------ angles

def merge_noise(A, tol=1e-6):
    T = np.array([fold(t) for t in A[:, 2]])
    order = np.argsort(T)
    groups, cur = [], [order[0]]
    for a, b in zip(order[:-1], order[1:]):
        if T[b] - T[a] < tol:
            cur.append(b)
        else:
            groups.append(cur); cur = [b]
    groups.append(cur)
    # wrap: -45 and 45 are the same orientation
    for g in groups:
        v = 0.0 if any(abs(T[i]) < tol for i in g) else float(np.median(T[g]))
        T[g] = v
    T[np.abs(np.abs(T) - 45.0) < tol] = 45.0
    B = A.copy(); B[:, 2] = T
    return B


def P_of(S, A):
    return dict(n=len(A), S=S + PAD, X=list(A[:, 0]), Y=list(A[:, 1]), T=list(A[:, 2]))


def feasible(S, A):
    import reg
    _, _, st, _ = reg.place(P_of(S, A), [(0, 0)] * len(A), feas_only=True, iters=1)
    return st != 'infeasible'


def straighten(S, A, near=1.0):
    T = A[:, 2]
    idx = np.nonzero((np.abs(T) < near) & (T != 0))[0]
    if not len(idx):
        return A, 0
    B = A.copy(); B[idx, 2] = 0.0
    if feasible(S, B):
        return B, len(idx)
    return A, 0


# ------------------------------------------------------------------------------------------------ orientation

def image(S, A, g):
    X, Y, T = A[:, 0].copy(), A[:, 1].copy(), A[:, 2].copy()
    for _ in range(g & 3):
        X, Y = S - Y, X.copy()
        T = T + 90
    if g & 4:
        X = S - X
        T = -T
    return np.c_[X, Y, [fold(t) for t in T]]


def mismatch(SA, A, SB, B):
    """assignment cost between two drawings (positions in units of the side; angle term per square)"""
    a, b = A[:, :2] / SA, B[:, :2] / SB
    C = ((a[:, None, :] - b[None, :, :]) ** 2).sum(-1) * (max(SA, SB) ** 2)
    dt = np.abs(A[:, 2][:, None] - B[:, 2][None, :]) % 90
    dt = np.minimum(dt, 90 - dt)
    C = C + (dt / 45.0) ** 2 * 0.5
    if len(A) != len(B):          # never happens for one n, but be safe
        return math.inf
    r, c = linear_sum_assignment(C)
    return float(C[r, c].sum())


def orient_like(S, A, ref):
    if ref is None:
        return 0, A
    SR, R = ref
    best = min(range(8), key=lambda g: (round(mismatch(S, image(S, A, g), SR, R), 9), g))
    return best, image(S, A, best)


# ------------------------------------------------------------------------------------------------ placement

def place_gravity(P, rep):
    """reg.py's 'pretty' rules 2-3 without the wall pull and without imposed symmetry: tilted clusters closed up
    (L1 pull to their centroid), straight tilted chains as equalities where feasible (all, else per chain), axis
    squares under gravity to the bottom left.  (The wall pull splits grid blocks in the middle of a wall run and the
    contact shading then flickers between frames; symmetry is left to the source drawings.)"""
    import reg
    weights, absobj = reg.pretty_objective(P, wgrid=0.0)
    align = reg.chain_pairs(P)
    X, Y, st, info = reg.place(P, weights, align=align, absobj=absobj)
    if st == 'infeasible' and align:
        comp = {}
        for i, j, _, _ in align:
            ci, cj = comp.setdefault(i, {i}), comp.setdefault(j, {j})
            if ci is not cj:
                ci |= cj
                for q in cj:
                    comp[q] = ci
        chains = sorted(([a for a in align if a[0] in c] for c in {id(v): v for v in comp.values()}.values()),
                        key=len, reverse=True)
        used = []
        for ch in chains:
            if reg.place(P, weights, align=used + ch, absobj=absobj, iters=1, feas_only=True)[2] != 'infeasible':
                used += ch
        X, Y, st, info = reg.place(P, weights, align=used, absobj=absobj)
        align = used
    rep['align'] = len(align)
    return X, Y, st, info


def regularise(S, A):
    import reg
    rep = {}
    P = P_of(S, A)
    try:
        X, Y, st, info = place_gravity(P, rep)
    except Exception as e:          # noqa: BLE001
        st, rep['error'] = 'error', repr(e)
    if st in ('infeasible', 'error'):
        X, Y, st, info = reg.place(P, [(1e-3, 1.0)] * len(A))
        rep['fallback'] = 'gravity'
    rep['place'] = st
    if st == 'infeasible':
        return A, rep
    B = np.c_[X, Y, A[:, 2]]
    if not reg.overlap_free(dict(n=len(B), S=S + PAD, X=list(X), Y=list(Y), T=list(A[:, 2])), tol=-1e-7):
        rep['place'] = 'overlap: source kept'
        return A, rep
    return B, rep


def view(S, A, ref):
    t0 = time.process_time()
    rep = {}
    if len(A) > 1:
        A = merge_noise(A)
        A, rep['straightened'] = straighten(S, A)
    g, A = orient_like(S, A, ref)
    rep['g'] = g
    if len(A) > 1 and np.any(np.abs(A[:, 2]) > 0) or not np.allclose(A[:, :2] % 1, 0.5):
        A, r = regularise(S, A)
        rep.update(r)
    rep['cpu_s'] = round(time.process_time() - t0, 2)
    return A, rep


def key_of(n, st, ref_key):
    h = hashlib.sha1(json.dumps([VERSION, n, st['draw'], st['s'], ref_key], sort_keys=True).encode()).hexdigest()[:16]
    return f'{n:03d}-{h}'


BEST = os.path.join(ROOT, 'site', 'www', 'data', 'p')


def do_n(args):
    n, states, force = args
    out, ref, ref_key = [], None, None
    log = []
    for i in range(len(states) - 1, -1, -1):        # last state first: it anchors the orientation
        st = states[i]
        if st['draw'] is None:
            out.append(None); continue
        k = key_of(n, st, ref_key)
        path = os.path.join(VIEWS, k + '.json')
        if os.path.exists(path) and not force:
            d = json.load(open(path))
        else:
            S, A = load(st['draw'], n)
            extra = len(A) - n
            if extra > 0:          # catalogue entries drawn as a larger n's packing ("remove a square from s(296)")
                A = remove(S, A, extra)
            elif extra < 0:
                raise ValueError(f'n={n} state {i}: {len(A)} squares')
            bp = os.path.join(BEST, f'best-{n}.json')
            if i == len(states) - 1 and os.path.exists(bp) and abs(float(json.load(open(bp))['s']) - S) < 1e-8:
                b = json.load(open(bp))
                B, rep = np.array(b['squares'], float), dict(src='lists/best (certified regularised record)')
                B[:, 2] = [fold(t) for t in B[:, 2]]
            elif st['draw']['kind'] == 'grid':          # grids stay canonical (holes top right) in every frame
                B, rep = A, dict(src='grid')
            else:
                B, rep = view(S, A, ref)
            if extra > 0:
                rep['removed_from_drawing'] = extra
            d = dict(n=n, s=S, t=st['t'], src=st['src'], squares=B.round(12).tolist(), rep=rep)
            json.dump(d, open(path, 'w'))
        log.append((i, d['rep']))
        out.append(k)
        ref, ref_key = (d['s'], np.array(d['squares'])), k
    return n, out[::-1], log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--procs', type=int, default=3)
    ap.add_argument('--ns', default='1-324')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    lo, hi = (int(x) for x in a.ns.split('-')) if '-' in a.ns else (int(a.ns), int(a.ns))
    H = json.load(open(os.path.join(OUT, 'history.json')))['hist']
    os.makedirs(VIEWS, exist_ok=True)
    idx_path = os.path.join(OUT, 'views.json')
    idx = json.load(open(idx_path)) if os.path.exists(idx_path) else {}
    jobs = [(n, H[str(n)], a.force) for n in range(hi, lo - 1, -1)]
    t0 = time.time()
    with Pool(a.procs) as pool:
        for n, keys, log in pool.imap_unordered(do_n, jobs):
            idx[str(n)] = keys
            cpu = sum(r.get('cpu_s', 0) for _, r in log)
            bad = [(i, r.get('place')) for i, r in log if r.get('place') not in (None, 'ok', 'stalled')]
            print(f'n={n} states={len(keys)} cpu={cpu:.1f}s' + (f' issues={bad}' if bad else ''), flush=True)
            json.dump(idx, open(idx_path, 'w'))
    print(f'done in {time.time() - t0:.0f}s wall')


if __name__ == '__main__':
    main()
