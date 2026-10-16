#!/usr/bin/env python3
"""Composer v3: Mondrian v1 (two corner pieces + exact upright fill) with a per-piece phase table (10-15 night, AD1).

AD1 (PACKER "Composition count structure"): a piece's single-corner excess e(P, S) = fill(S, P at LL) - floor(S)^2 depends
only on frac(S) (300 / 300), and it is nondecreasing in frac(S) (a fill at S fits at any larger S with the same floor;
fillfast's candidate set is not, 5 / 443k table entries dip by 1: the table takes the running max).
The pair sum eA + eB never fell below the joint count on the pairs that reach n (0 false negatives in 2,131; 21 % false
positives), so it is the screen; the joint fill (fillfast, ~60 ms) decides.

Stages per event n:
  1. table: e(P, j / G), j = 0..G (j = G: just below the next integer), computed once per library piece by monotone
     bisection (equal ends -> the interval is constant), cached per library.  Screen bound at S: e at ceil(frac(S) G).
  2. sides: ante + lo .. ante + hi (--sides values), clipped below the next integer k; if the window reaches k, the
     side k - 1e-6 is added.
  3. pairs: A at LL, B at UR (180-degree turn) or LR (mirror), disjoint boxes, i <= j (swap = the whole packing
     rotated / mirrored), bound >= need; ranked by bound then tilted count; joint fill in rank order until --reach-cap
     pairs reach n or --fill-cap fills are spent.
  4. polish: the --polish-k reaching compositions per side with the largest count surplus (then bound) -> trim to n,
     fq quench + polish.  Best --per-n distinct sides -> starts, then the explorer (harden_explore.py) if --minutes > 0.
  Above-integer arm (Evan 10-15: compositions that need more than the integer may still relax below it): when the
  next integer k is within --above-window of ante, the --above-k pairs closest to fitting at the ante (table bound
  there minus need) are composed at k + d for each --above-d (the grid makes the count trivial), trimmed to n
  (uprights nearest the tilted squares kept), quenched + polished; logged with est_ante (calibration of P(< k)), kept
  as starts if below k.  (V1's version, pairs 1-3 short at k - 1e-6, never fired: there every top pair reaches, and
  those polish to k exactly.)

  python3 analysis/compose3.py table --lib runs/mondrian/croplib481t1.npz
  python3 analysis/compose3.py run --set 481 --lib runs/mondrian/croplib481t1.npz --out runs/compose3/v481 [--minutes 45]
  python3 analysis/compose3.py run --live 100-324 --lib LIB --out runs/compose3/live1
"""
import argparse, glob, hashlib, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

FT = 1e-9
G = 128
_P = None          # worker-global piece list [(T, W, H, m)]


def fl(x):
    return int(math.floor(x + FT))


# ---------------------------------------------------------------- pieces

def _store_pieces(source=None, per_cell=3, min_tilted=3):
    """store piece table (rectsearch / rectlib / imported crops): the best per_cell distinct h per (m, d to 1e-4), with
    >= min_tilted squares at a real angle (0.5..89.5 deg); tilted squares only, LL frame."""
    import mondrian as M
    from pk.store import Store
    cnt, seen, out = {}, set(), []
    for p in Store().pieces(mmin=3, source=source):
        if p.meta['source'].endswith(':atlas'):        # atlas crops: croplib covers these (cleaner cuts)
            continue
        m, d, H = p.n, float(p.dx), float(p.s)
        W = H + d
        a = np.mod(p.sq[:, 2], 90.0)
        if int(((a > 0.5) & (a < 89.5)).sum()) < min_tilted:
            continue
        key = (m, round(d, 4))
        if (key, round(H, 7)) in seen or cnt.get(key, 0) >= per_cell:
            continue
        seen.add((key, round(H, 7))); cnt[key] = cnt.get(key, 0) + 1
        T = p.sq[M.is_tilted(p.sq[:, 2])].copy()
        for TT, WW, HH in ((T, W, H), (np.c_[T[:, 1], T[:, 0], (-T[:, 2]) % 90.0], H, W)) if abs(d) > 1e-9 else ((T, W, H),):
            out.append((TT, WW, HH, m, f"piece:{p.meta['source']}#{p.meta['id']}"))
    return out


def load(lib):
    """lib = comma list of croplib .npz files and 'store' / 'store:SOURCE_LIKE' (store piece table).  Deduplicated
    (same box, m and tilted squares to 1e-6) -> list of (T, W, H, m, label)"""
    import l2_screen
    P, seen = [], set()
    for part in lib.split(','):
        if part.startswith('store'):
            src = part.split(':', 1)[1] if ':' in part else None
            items = _store_pieces(src)
        else:
            items = [(T, meta[3], meta[4], meta[5], lab) for T, meta, lab in l2_screen.load_lib(part, 10 ** 6)]
        for T, W, H, m, lab in items:
            T = T[np.lexsort((T[:, 1], T[:, 0]))]
            key = hashlib.sha1(np.round(np.r_[[W, H, m], T.ravel()], 6).tobytes()).hexdigest()
            if key in seen:
                continue
            seen.add(key)
            P.append((T, W, H, m, lab))
    return P


def _init(lib):
    """lib string (worker loads it) or the loaded piece list itself (10-16: the store changes under a running job --
    rectsearch imports at its end -- so workers must use the parent's list, not a fresh load)"""
    global _P
    _P = load(lib) if isinstance(lib, str) else lib


def count(S, T):
    import fillfast
    return len(fillfast.fill(S, T))


# ---------------------------------------------------------------- 1. phase table

def _table_job(k):
    T, W, H, m, _ = _P[k]
    base = math.ceil(max(W, H) - 1e-12)
    S_of = lambda j: base + (j / G if j < G else 1 - 1e-6)
    e = np.full(G + 1, -999, np.int32)
    nf = 0

    def ev(j):
        nonlocal nf
        if e[j] == -999:
            S = S_of(j); e[j] = count(S, T) - fl(S) ** 2; nf += 1
        return e[j]

    stack = [(0, G)]
    while stack:
        lo, hi = stack.pop()
        a, b = ev(lo), ev(hi)
        if a == b:
            e[lo:hi + 1] = a
        elif hi - lo > 1:
            mid = (lo + hi) // 2
            stack += [(lo, mid), (mid, hi)]
    return k, e, nf


def piece_keys(P):
    return [hashlib.sha1(np.round(np.r_[[W, H, m], T.ravel()], 6).tobytes()).hexdigest()[:16] for T, W, H, m, _ in P]


def lib_digest(P):
    h = hashlib.sha1()
    for T, W, H, m, lab in P:
        h.update(np.round(np.r_[[W, H, m], T.ravel()], 6).tobytes())
    return h.hexdigest()[:12]


def table_path(lib, P=None):
    if ',' not in lib and not lib.startswith('store'):
        name = os.path.basename(lib).replace('.npz', f'_G{G}.npz')
    else:          # 10-16: keyed by the pieces too (the store grows; a stale table would misalign rows)
        name = 'lib_' + hashlib.sha1(lib.encode()).hexdigest()[:12] + (f'_{lib_digest(P)}' if P is not None else '') + f'_G{G}.npz'
    return os.path.join(PK, 'runs', 'compose3', 'tables', name)


def table(lib, procs, P=None):
    P = P or load(lib)
    path = table_path(lib, P)
    if os.path.exists(path):
        E = np.load(path)['e']
        assert E.shape[0] == len(P), (path, E.shape, len(P))
        return np.maximum.accumulate(E, axis=1)
    t0 = time.time()
    K = piece_keys(P)
    E = np.zeros((len(P), G + 1), np.int32); nf = 0
    have = {}                                    # 10-16: per-piece rows from earlier tables that saved their keys
    for f in sorted(glob.glob(os.path.join(os.path.dirname(path), '*.npz'))):
        z = np.load(f)
        if 'keys' in z.files:
            for kk, row in zip(z['keys'], z['e']):
                have.setdefault(str(kk), row)
    todo = [i for i in range(len(P)) if K[i] not in have]
    for i in range(len(P)):
        if K[i] in have:
            E[i] = have[K[i]]
    with ProcessPoolExecutor(procs, initializer=_init, initargs=(P,)) as ex:
        for k, e, n_ in ex.map(_table_job, todo, chunksize=8):
            E[k] = e; nf += n_
    bad = (np.diff(E, axis=1) < 0).sum()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, e=E, keys=np.array(K))
    print(f'table: {len(P)} pieces ({len(P) - len(todo)} reused), {nf} fills ({nf / max(1, len(todo)):.1f} per new piece), '
          f'{time.time() - t0:.0f} s; monotonicity violations {bad} (candidate set: wall-anchored lattices shift with S)', flush=True)
    return np.maximum.accumulate(E, axis=1)      # any packing at S' < S (same floor) fits at S


# ---------------------------------------------------------------- 3. pairs

def bound(E, S):
    j = min(G, math.ceil((S - math.floor(S)) * G - 1e-9))
    return E[:, j]


def pairs(n, S, E, W, H, m, t, need, margin=0, cap=200000, score=None, req=None):
    """candidate pairs (score, tilted, i, j, corner), disjoint boxes at S, best first, at most cap.  score = table bound
    at S (default; only pairs with bound >= need - margin) or a given per-piece vector (pair score = sum, no threshold).
    req: bool per piece; if given, at least one of the pair must have it (--new-vs: pairs with a new piece only)."""
    ok = np.nonzero((m <= n - 3) & (W <= S + 1e-12) & (H <= S + 1e-12))[0]
    if len(ok) == 0 or (req is not None and not req[ok].any()):
        return []
    rq = req[ok] if req is not None else None
    u = (bound(E, S) if score is None else np.asarray(score))[ok]; w, h, tt = W[ok], H[ok], t[ok]
    thr = need - margin if score is None else -np.inf
    keys, rows = [], []
    for a in range(0, len(ok), 512):
        sl = slice(a, a + 512)
        est = (u[sl, None] + u[None, :]).astype(float)
        tri = np.arange(a, min(a + 512, len(ok)))[:, None] <= np.arange(len(ok))[None, :]
        good = (est >= thr) & tri
        if rq is not None:
            good &= rq[sl, None] | rq[None, :]
        if not good.any():
            continue
        sw = (w[sl, None] + w[None, :] <= S + 1e-12)
        sh = (h[sl, None] + h[None, :] <= S + 1e-12)
        for ci, dis in ((0, sw | sh), (1, sw)):
            I, J = np.nonzero(good & dis)
            if len(I) == 0:
                continue
            sc, tl = est[I, J], tt[sl][I] + tt[J]
            if len(I) > cap:                                    # keep this chunk's best cap only
                q = np.lexsort((-tl, -sc))[:cap]; I, J, sc, tl = I[q], J[q], sc[q], tl[q]
            rows.append(np.c_[sc, tl, ok[I + a], ok[J], np.full(len(I), ci)])
        if rows and sum(len(r) for r in rows) > 4 * cap:        # compact
            X = np.concatenate(rows); X = X[np.lexsort((-X[:, 1], -X[:, 0]))[:cap]]; rows = [X]
    if not rows:
        return []
    X = np.concatenate(rows); X = X[np.lexsort((-X[:, 1], -X[:, 0]))[:cap]]
    return [(int(r[0]), int(r[1]), int(r[2]), int(r[3]), ('UR', 'LR')[int(r[4])]) for r in X]


def compose(S, a, b, corner):
    import mondrian as M
    A, B = M.to_corner(_P[a][0], 'LL', S), M.to_corner(_P[b][0], corner, S)
    return np.r_[A, M._drop_conflicts(A, B)]


def _joint_job(args):
    S, a, b, corner = args
    return count(S, compose(S, a, b, corner))


def _polish_job(args):
    n, S, a, b, corner = args
    import mondrian as M, fillfast
    from pk.engine import FQ, QOpt, POLISH
    from pk.packing import Packing
    t0 = time.process_time()
    L = fillfast.fill(S, compose(S, a, b, corner))
    if len(L) < n:
        return None
    eng = FQ()
    try:
        r = eng.quench(Packing(S, M.trim(L, n)), QOpt(), timeout=120)
        r2 = eng.quench(r.p, POLISH, timeout=120) if r.ok else r
    finally:
        eng.close()
    if not r2.ok:
        return None
    return dict(s=float(r2.s), sq=r2.p.sq.tolist(), cpu=time.process_time() - t0)


# ---------------------------------------------------------------- driver

def events(a):
    import l2_compose as L2
    if a.live:                                   # non-grid frontier n only (a grid record is not beaten by composing above it)
        return [e for e in L2.live_set(a.live, a.T) if e['ante'] is not None and math.isfinite(e['ante'])
                and abs(e['ante'] - round(e['ante'])) > 1e-9]
    evs = L2.event_set(a.set)
    if a.only:
        ks = {int(x) for x in a.only.split(',')}; evs = [e for e in evs if e['n'] in ks]
    return evs


def run(a):
    from pk import manifest
    from pk.packing import Packing
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'compose3', args=vars(a), procs=a.cores)
    P = load(a.lib)
    E = table(a.lib, a.cores, P)
    W = np.array([p[1] for p in P]); H = np.array([p[2] for p in P]); m = np.array([p[3] for p in P])
    t = np.array([len(p[0]) for p in P])
    REQ = None
    if a.new_vs:                                   # only pairs with a piece absent from that table (new pieces)
        old = set(str(k) for k in np.load(a.new_vs)['keys'])
        REQ = np.array([k not in old for k in piece_keys(P)])
        print(f'--new-vs: {REQ.sum()} new pieces of {len(P)}', flush=True)
    log = open(f'{a.out}/compose.jsonl', 'a')
    ex = ProcessPoolExecutor(a.cores, initializer=_init, initargs=(P,))
    for ev in events(a):
        n, ante = ev['n'], float(ev['ante']); k = fl(ante) + 1
        d = f'{a.out}/starts/n{n}'; os.makedirs(d, exist_ok=True)
        sides = [s for s in np.linspace(ante + a.lo, ante + a.hi, a.sides) if s < k - 1e-7]
        crosses = ante + a.hi >= k - 1e-7
        if crosses:
            sides.append(k - 1e-6)
        comp = []          # (polished side, sq, S, tag)
        polished_n = 0
        for S in sides:
            if polished_n >= a.reach_budget:
                print(f'n={n} S={S:.6f}: skipped (per-n polish budget {a.reach_budget} spent at lower sides)', flush=True)
                continue
            t1 = time.time(); need = n - fl(S) ** 2
            C = pairs(n, S, E, W, H, m, t, need, req=REQ)
            reach, nfill = [], 0
            for b0 in range(0, min(len(C), a.fill_cap), a.cores * 8):
                batch = C[b0:min(b0 + a.cores * 8, a.fill_cap)]
                cnt = list(ex.map(_joint_job, [(S, i, j, c) for _, _, i, j, c in batch], chunksize=2))
                nfill += len(batch)
                reach += [(x - n, c) for x, c in zip(cnt, batch) if x >= n]
                if len(reach) >= a.reach_cap:
                    break
            reach.sort(key=lambda r: (-r[0], -r[1][0], -r[1][1]))
            sel = reach[:min(a.polish_k, a.reach_budget - polished_n)]
            polished_n += len(sel)
            R = list(ex.map(_polish_job, [(n, S, i, j, c) for _, (_, _, i, j, c) in sel], chunksize=1))
            for (sur, (bd, tl, i, j, c)), r in zip(sel, R):
                log.write(json.dumps(dict(n=n, S=S, arm='side', a=P[i][4], b=P[j][4], corner=c, bound=bd, surplus=sur,
                                          polished=r and r['s'], cpu=r and r['cpu'])) + '\n')
                if r:
                    comp.append((r['s'], r['sq'], S, 'side'))
            best = min((r['s'] for r in R if r), default=None)
            print(f'n={n} S={S:.6f} need={need}: {len(C)} pairs pass the bound, {nfill} joint fills, {len(reach)} reach, '
                  f'{len(sel)} polished, best {best} ({time.time() - t1:.0f} s)', flush=True)
        if k - ante <= a.above_window and a.above_k > 0:
            # above-integer arm: compose at k + d (the grid makes the count trivial), polish, see what relaxes below k.
            # Pairs ranked by est_ante = table bound at the ante's phase minus the need there (how close the pair
            # comes to fitting at the target itself); logged for calibration of P(< k | est_ante).
            t1 = time.time(); need_a = n - fl(ante) ** 2
            for dd in a.above_d:
                Sa = k + dd
                C = pairs(n, Sa, E, W, H, m, t, 0, cap=a.above_k, score=bound(E, ante))
                C = [(sc - need_a, tl, i, j, c) for sc, tl, i, j, c in C]
                R = list(ex.map(_polish_job, [(n, Sa, i, j, c) for _, _, i, j, c in C], chunksize=1))
                for (ea, tl, i, j, c), r in zip(C, R):
                    log.write(json.dumps(dict(n=n, S=Sa, arm='above', est_ante=ea, tilted=tl, a=P[i][4], b=P[j][4],
                                              corner=c, polished=r and r['s'], below_k=bool(r and r['s'] < k),
                                              below_ante=bool(r and r['s'] < ante))) + '\n')
                    if r and r['s'] < k:
                        comp.append((r['s'], r['sq'], Sa, 'above'))
                ok = [r['s'] for r in R if r and r['s'] < k]
                print(f'n={n} above k+{dd}: {len(C)} composed + polished, {len(ok)} relax below {k}, best '
                      f'{min(ok) if ok else None} ({time.time() - t1:.0f} s)', flush=True)
        log.flush()
        comp.sort(key=lambda c: c[0])
        seen, kk = [], 0
        for s, sq, S, tag in comp:
            if all(abs(s - x) > 1e-9 for x in seen):
                seen.append(s); Packing(s, np.array(sq)).write(f'{d}/start{kk}.txt'); kk += 1
            if kk >= a.per_n:
                break
        ev['near_best'] = seen[0] if seen else None
        ev['starts'] = seen
        json.dump(ev, open(f'{d}/event.json', 'w'), indent=1, default=str)
        print(f'n={n}: starts {[round(x, 8) for x in seen]}  ante {ante:.10f}', flush=True)
    ex.shutdown()
    if a.minutes > 0:
        import subprocess
        subprocess.run([sys.executable, f'{PK}/harden_explore.py', 'run', '--src', f'{a.out}/starts', '--starts-glob',
                        'start[0-9]*.txt', '--out', f'{a.out}/explore', '--minutes', str(a.minutes), '--cores',
                        str(a.cores)], cwd=PK, check=True)
    man.finish()


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('table'); p.add_argument('--lib', required=True); p.add_argument('--procs', type=int, default=14)
    r = sp.add_parser('run'); r.add_argument('--lib', required=True); r.add_argument('--out', required=True)
    r.add_argument('--set', default='481'); r.add_argument('--only'); r.add_argument('--live'); r.add_argument('--T')
    r.add_argument('--lo', type=float, default=0.04); r.add_argument('--hi', type=float, default=0.16)
    r.add_argument('--sides', type=int, default=4)
    r.add_argument('--new-vs', help='table .npz with keys: only pairs containing a piece not in it')
    r.add_argument('--fill-cap', type=int, default=4000); r.add_argument('--reach-cap', type=int, default=200)
    r.add_argument('--reach-budget', type=int, default=1500, help='per n: skip higher sides once this many polished')
    r.add_argument('--polish-k', type=int, default=8); r.add_argument('--per-n', type=int, default=4)
    r.add_argument('--above-k', type=int, default=0, help='pairs composed per above-integer side')
    r.add_argument('--above-window', type=float, default=0.1, help='above arm when the next integer is within this of ante')
    r.add_argument('--above-d', type=float, nargs='+', default=[0.02, 0.06])
    r.add_argument('--minutes', type=float, default=0); r.add_argument('--cores', type=int, default=14)
    a = ap.parse_args()
    if a.cmd == 'table':
        table(a.lib, a.procs)
    else:
        run(a)


if __name__ == '__main__':
    main()
