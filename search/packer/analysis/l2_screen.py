#!/usr/bin/env python3
"""Blind Mondrian v1 composition (L2; 10-15): the replicate pair screen over every allowed atlas crop instead of the
named pieces.  compose2's bandit is search-bound on 1e8-1e9 arms (L2b); the guillotine count formula is cheap and
vectorised, so it screens every pair, and the exact fill only sees the top.

Allowed sources (time travel): Kingbird 09-30 mirror crops + crops of store packings at the register frontier as of
T_utc (register history) or added to the store before T_local.  Pieces = atlas_v3 crops (corner crops, no dropped
rows / columns: a known gap vs #481, which shifted one piece by a unit row at 259).

Per event n and per side Sc in [S_lo, S_hi] (grid of --sides values from ante to ante + --window): formula best over
all A x B pairs (chunked), the --keep best pairs -> mondrian.replicate on just those pieces (joint exact fill, fq) ->
start files (best distinct composed sides).  Then harden_explore.py (control arm: explorer only) from them.

  python3 analysis/l2_screen.py run --set 481 --out runs/l2/s481 [--window 0.3 --sides 6 --keep 60 --minutes 20]
  python3 analysis/l2_screen.py score --out runs/l2/s481
"""
import argparse, json, math, os, sqlite3, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)
FT = 1e-9
LIB = None                       # run --lib: analysis/croplib.py library (clean crops, dropped rows) instead of the atlas
MIN_EXCESS = None                # run --min-excess: rebuild crops from the allowed sources with this excess floor


def fl(x):
    return np.floor(np.asarray(x) + FT)


def allowed_sids(z, T_utc, T_local):
    from pk import eventbench as EB
    db = sqlite3.connect(f"file:{os.path.join(PK, 'runs', 'store.sqlite')}?mode=ro", uri=True)
    pre = {r[0] for r in db.execute('SELECT id FROM packing WHERE added < ?', (T_local,))}
    ids = {int(str(l).split(':')[1]) for l in z['labels'] if str(l).startswith('store:')}
    fr = {n: EB.frontier_at(n, T_utc) for n in range(2, 401)}
    q = ','.join(map(str, ids))
    for pid, n, s in db.execute(f'SELECT id, n, s FROM packing WHERE id IN ({q})'):   # only the atlas's sources
        f = fr.get(n)
        if f is not None and math.isfinite(f) and abs(s - f) < 2e-9:
            pre.add(pid)
    ok = []
    for k, lab in enumerate(z['labels']):
        lab = str(lab)
        if lab.startswith('kingbird:'):
            ok.append(k)
        elif lab.startswith('store:') and int(lab.split(':')[1]) in pre:
            ok.append(k)
    return set(ok)


def load_lib(lib, nmax):
    """croplib.py library (already restricted to the cutoff's allowed sources): (T, (o, i, j, W, H, m), label)"""
    import mondrian as M
    L = dict(np.load(lib, allow_pickle=True))
    z = dict(np.load(str(L['atlas']), allow_pickle=True))
    P = []
    for W, H, m, o, i, j, sid in zip(L['W'], L['H'], L['m'], L['o'], L['i'], L['j'], L['sid']):
        if m > nmax:
            continue
        q = M.orient(np.asarray(z['sq'][sid], float), float(z['src_s'][sid]), int(o))
        q[:, 0] -= i; q[:, 1] -= j
        h = M.half_ext(q[:, 2])
        inn = (q[:, 0] - h >= -1e-9) & (q[:, 1] - h >= -1e-9) & (q[:, 0] + h <= W + 1e-9) & (q[:, 1] + h <= H + 1e-9)
        T = q[inn & M.is_tilted(q[:, 2])]
        if len(T):
            P.append((T, (int(o), int(i), int(j), float(W), float(H), int(m)), str(z['labels'][sid])))
    return P


def load_pieces(atlas, T_utc, T_local, nmax):
    import mondrian as M
    if LIB:
        return load_lib(LIB, nmax)
    z = dict(np.load(atlas, allow_pickle=True))          # materialise: NpzFile re-reads an array on every access
    A = allowed_sids(z, T_utc, T_local)
    if MIN_EXCESS is not None:                           # rebuild the allowed sources' crops with a relaxed excess
        parts = []                                       # (Mondrian v1 pieces need not beat the grid on their own)
        for k in sorted(A):
            r = M.pieces_of(np.asarray(z['sq'][k], float), float(z['src_s'][k]), min_excess=MIN_EXCESS)
            if len(r[0]):
                r = M.pareto_prune(*r)
                parts.append(tuple(r) + (np.full(len(r[0]), k),))
        W, H, m, o, sid = (np.concatenate(c) for c in zip(*parts))
    else:
        W, H, m, o, sid = z['W'], z['H'], z['m'], z['o'], z['sid']
    keep = np.array([s in A for s in sid]) & (m >= 3) & (m <= nmax)
    P = []
    for k in np.nonzero(keep)[0]:
        q = M.orient(np.asarray(z['sq'][sid[k]], float), float(z['src_s'][sid[k]]), int(o[k]))
        h = M.half_ext(q[:, 2])
        inn = (q[:, 0] + h <= W[k] + 1e-9) & (q[:, 1] + h <= H[k] + 1e-9)
        T = q[inn & M.is_tilted(q[:, 2])]
        if len(T) == 0:
            continue
        P.append((T, (int(o[k]), 0, 0, float(W[k]), float(H[k]), int(m[k])), str(z['labels'][sid[k]])))
    return P


def screen(n, S, W, H, m, keep, chunk=400):
    """formula best count for every (A at LL, B at the opposite corner) pair: the replicate guillotine formula, over
    both cut directions and k = 0..2 offsets; returns the `keep` best (count, ia, ib)."""
    def col(w, wp, hp, mm):
        ok = w >= wp - 1e-12
        v = mm + np.maximum(fl(S - hp) * fl(w) + fl(w - wp) * fl(hp), fl(w - wp) * fl(S) + fl(S - hp) * fl(wp))
        return np.where(ok, v, -10 ** 6)
    fit = (W <= S + 1e-12) & (H <= S + 1e-12)
    idx = np.nonzero(fit)[0]
    Wf, Hf, mf = W[idx], H[idx], m[idx]
    out = []
    for a0 in range(0, len(idx), chunk):
        aw0, ah0, am = Wf[a0:a0 + chunk, None], Hf[a0:a0 + chunk, None], mf[a0:a0 + chunk, None]
        best = np.full((len(aw0), len(idx)), -10 ** 6, float)
        for tr in (False, True):
            a_w, a_h = (ah0, aw0) if tr else (aw0, ah0)
            b_w, b_h = (Hf[None, :], Wf[None, :]) if tr else (Wf[None, :], Hf[None, :])
            for k in range(3):
                for left in (True, False):
                    c = (a_w + k) if left else (S - b_w - k)
                    cnt = col(c, a_w, a_h, am) + col(S - c, b_w, b_h, mf[None, :])
                    feas = (a_w + b_w <= S + 1e-12) & (c >= a_w - 1e-12) & (S - c >= b_w - 1e-12)
                    best = np.maximum(best, np.where(feas, cnt, -10 ** 6))
        I, J = np.nonzero(best >= n - 8)
        if len(I):
            v = best[I, J]
            top = np.argsort(-v)[:keep]
            out += [(float(v[t]), int(idx[a0 + I[t]]), int(idx[J[t]])) for t in top]
        out.sort(key=lambda r: -r[0]); del out[keep:]
    return out


def _ex1(args):
    T, S = args
    import mondrian as M
    if not M.fits(T, S):
        return -99
    return len(M.exact_fill(T, S, time_limit=20)) - int(math.floor(S + FT)) ** 2


def _pair_job(args):
    """joint fill of A at LL + B at `corner` (both priority orders); count >= n -> fq + polish"""
    n, S, TA, TB, corner = args
    import mondrian as M
    from pk.engine import FQ, QOpt, POLISH
    from pk.packing import Packing
    A, B = M.to_corner(TA, 'LL', S), M.to_corner(TB, corner, S)
    out = []
    for T in (np.r_[A, M._drop_conflicts(A, B)], np.r_[M._drop_conflicts(B, A), B]):
        L = M.exact_fill(T, S, time_limit=10)
        if len(L) < n:
            out.append((len(L), None)); continue
        eng = FQ()
        r = eng.quench(Packing(S, M.trim(L, n)), QOpt(), timeout=120)
        r2 = eng.quench(r.p, POLISH, timeout=120) if r.ok else r
        out.append((len(L), (r2.s, r2.p.sq.tolist()) if r2.ok else None))
        eng.close()
    return out


def exact_screen(n, S, P, keep, procs):
    """additive exact excess: e(piece) = single-piece MILP fill at side S minus the grid; pairs with disjoint boxes
    (A at LL; B at UR: wA + wB <= S or hA + hB <= S; B at LR: wA + wB <= S) ranked by eA + eB (then tilted count);
    the top `keep` get the joint fill.  -> list of (side, Packing-like (s, sq))"""
    from concurrent.futures import ProcessPoolExecutor
    g = int(math.floor(S + FT)) ** 2
    idx = [k for k, p in enumerate(P) if p[1][5] <= n - 3 and max(p[1][3], p[1][4]) <= S + 1e-12]
    with ProcessPoolExecutor(procs) as ex:
        e = np.array(list(ex.map(_ex1, [(P[k][0], S) for k in idx], chunksize=8)))
    W = np.array([P[k][1][3] for k in idx]); H = np.array([P[k][1][4] for k in idx]); t = np.array([len(P[k][0]) for k in idx])
    est = (e[:, None] + e[None, :]).astype(float)
    cand = []
    for corner, ok in (('UR', (W[:, None] + W[None, :] <= S + 1e-12) | (H[:, None] + H[None, :] <= S + 1e-12)),
                       ('LR', (W[:, None] + W[None, :] <= S + 1e-12))):
        I, J = np.nonzero(ok & (est >= n - g - 1) & (e[:, None] > -99) & (e[None, :] > -99))
        cand += [(est[i, j], t[i] + t[j], idx[i], idx[j], corner) for i, j in zip(I, J)]
    cand.sort(key=lambda c: (-c[0], -c[1]))
    cand = cand[:keep]
    if not cand:
        return [], 0, max(e.max() * 2, -99) if len(e) else None
    with ProcessPoolExecutor(procs) as ex:
        R = list(ex.map(_pair_job, [(n, S, P[a][0], P[b][0], c) for _, _, a, b, c in cand], chunksize=1))
    comp = [x[1] for r in R for x in r if x[1] is not None]
    reach = sum(1 for r in R for x in r if x[0] >= n)
    return comp, reach, cand[0][0]


def run(a):
    import mondrian as M
    from pk import manifest
    sys.path.insert(0, HERE)
    import l2_compose as L2
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'l2_screen', args=vars(a), procs=a.cores)
    evs = L2.event_set(a.set)
    if a.only:
        ks = {int(x) for x in a.only.split(',')}; evs = [e for e in evs if e['n'] in ks]
    T_utc, T_local = evs[0]['T_utc'], evs[0]['T_local']
    global MIN_EXCESS
    MIN_EXCESS = a.min_excess
    global LIB
    LIB = a.lib
    t0 = time.time()
    P = load_pieces(a.atlas, T_utc, T_local, max(e['n'] for e in evs))
    W = np.array([p[1][3] for p in P]); H = np.array([p[1][4] for p in P]); m = np.array([p[1][5] for p in P])
    print(f'{len(P)} allowed pieces ({time.time() - t0:.0f} s)', flush=True)
    sdir = f'{a.out}/starts'
    for ev in evs:
        n = ev['n']; d = f'{sdir}/n{n}'; os.makedirs(d, exist_ok=True)
        comp = []
        for Sc in np.linspace(ev['ante'] + a.lo, ev['ante'] + a.window, a.sides):
            t1 = time.time()
            if a.rank == 'exact':
                from pk.packing import Packing
                cs, reach, top = exact_screen(n, float(Sc), P, a.keep, a.cores)
                comp += [(s_, Packing(s_, np.array(sq)), float(Sc)) for s_, sq in cs]
                print(f"n={n} Sc={Sc:.6f}: best estimate {top}, {reach} joint fills reach n -> {len(cs)} composed "
                      f"({time.time() - t1:.0f} s)", flush=True)
                continue
            pairs = screen(n, float(Sc), W, H, m, a.keep)
            if not pairs:
                continue
            ia = sorted({p[1] for p in pairs}); ib = sorted({p[2] for p in pairs})
            PA = [P[i][:2] for i in ia]; PB = [P[i][:2] for i in ib]
            res = M.replicate(n, float(Sc), PA, PB, procs=a.cores, top=a.keep, fq_top=4, log=lambda *x: None)
            for r in res:
                comp.append((r['polish'], r['packing'], float(Sc)))
            print(f"n={n} Sc={Sc:.6f}: formula best {pairs[0][0]:.0f}, {len(ia)} x {len(ib)} pieces -> "
                  f"{len(res)} composed ({time.time() - t1:.0f} s)", flush=True)
        comp.sort(key=lambda c: c[0])
        seen, k = [], 0
        for s, pk, Sc in comp:
            if all(abs(s - t) > 1e-9 for t in seen):
                seen.append(s); pk.write(f'{d}/start{k}.txt'); k += 1
            if k >= a.per_n:
                break
        ev['near_best'] = seen[0] if seen else None
        json.dump(ev, open(f'{d}/event.json', 'w'), indent=1)
        print(f"n={n}: {len(comp)} compositions, best {ev['near_best']}  ante {ev['ante']:.10f} posted {ev['S']:.10f}", flush=True)
    if a.minutes > 0:
        subprocess.run([sys.executable, f'{PK}/harden_explore.py', 'run', '--src', sdir, '--starts-glob', 'start[0-9]*.txt',
                        '--out', f'{a.out}/explore', '--minutes', str(a.minutes), '--cores', str(a.cores)], cwd=PK, check=True)
    man.finish()


def score(a):
    subprocess.run([sys.executable, f'{PK}/harden_explore.py', 'score', '--out', f'{a.out}/explore'], cwd=PK)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--set', default='481'); r.add_argument('--out', required=True)
    r.add_argument('--only'); r.add_argument('--atlas', default=os.path.join(PK, 'runs', 'mondrian', 'atlas_v3.npz'))
    r.add_argument('--window', type=float, default=0.3); r.add_argument('--sides', type=int, default=8)
    r.add_argument('--keep', type=int, default=100); r.add_argument('--per-n', type=int, default=4)
    r.add_argument('--minutes', type=float, default=20); r.add_argument('--cores', type=int, default=16)
    r.add_argument('--min-excess', type=int, help='rebuild crops with excess >= this (e.g. -3) instead of the atlas (>= 1)')
    r.add_argument('--lo', type=float, default=1e-9, help='first side = ante + lo')
    r.add_argument('--rank', choices=['formula', 'exact'], default='formula',
                   help='exact: additive single-piece MILP excess (10-15: the formula ranks #481 pairs thousands deep)')
    r.add_argument('--lib', help='analysis/croplib.py library npz (built for this event set\'s cutoff)')
    s = sp.add_parser('score'); s.add_argument('--out', required=True)
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
