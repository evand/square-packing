#!/usr/bin/env python3
"""Composition count structure (10-15 night, Evan: "isn't it just filling a rectangle?").

Two checks on the croplib481t1 library (the L2e library), counts by fillfast (exact; = gen.mis_fill2):
  period:  single-piece excess e(P, S) = fill(S, piece at LL) - floor(S)^2 at S = S0 + k, k = 0..3, fixed frac(S).
           Is e a function of frac(S) only (given the piece fits)?
  pairs:   at the 15 #481 events, S = ante + DS: e for every allowed piece, then joint fill of disjoint (UR / LR) pairs
           (top by eA + eB as in L2e, plus random disjoint pairs).  delta = (joint - floor(S)^2) - (eA + eB), against the
           box gaps gx = S - wA - wB, gy = S - hA - hB.

  python3 analysis/additivity.py period --n 300 --out runs/additivity/period.jsonl
  python3 analysis/additivity.py pairs --ds 0.08 --top 300 --rand 300 --out runs/additivity/pairs.jsonl
  python3 analysis/additivity.py report runs/additivity
"""
import argparse, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

LIB = os.path.join(PK, 'runs', 'mondrian', 'croplib481t1.npz')
FT = 1e-9
_P = None


def pieces():
    global _P
    if _P is None:
        import l2_screen
        _P = l2_screen.load_lib(LIB, 400)
    return _P


def g(S):
    return int(math.floor(S + FT)) ** 2


def count(S, T):
    import fillfast
    return len(fillfast.fill(S, T))


def _period(args):
    k, phi = args
    import mondrian as M
    T, meta, lab = pieces()[k]
    W, H = meta[3], meta[4]
    big = max(W, H)
    S0 = math.floor(big) + phi
    if S0 < big - 1e-12:
        S0 += 1
    out = []
    for j in range(4):
        S = S0 + j
        if not M.fits(T, S):
            out.append(None); continue
        out.append(count(S, T) - g(S))
    return dict(k=k, label=lab, W=W, H=H, m=meta[5], tilted=len(T), phi=phi, S0=S0, slack0=S0 - big, e=out)


def cmd_period(a):
    rng = np.random.default_rng(a.seed)
    P = pieces()
    ks = rng.choice(len(P), size=min(a.n, len(P)), replace=False)
    jobs = [(int(k), float(rng.uniform(0.02, 0.98))) for k in ks]
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.procs) as ex, open(a.out, 'w') as f:
        for r in ex.map(_period, jobs, chunksize=4):
            f.write(json.dumps(r) + '\n')
    print(f'period: {len(jobs)} pieces, {time.time() - t0:.0f} s', flush=True)


def _single(args):
    k, S = args
    import mondrian as M
    T = pieces()[k][0]
    if not M.fits(T, S):
        return None
    return count(S, T) - g(S)


def _joint(args):
    n, S, a, b, corner = args
    import mondrian as M
    P = pieces()
    A, B = M.to_corner(P[a][0], 'LL', S), M.to_corner(P[b][0], corner, S)
    c = []
    for T in (np.r_[A, M._drop_conflicts(A, B)], np.r_[M._drop_conflicts(B, A), B]):
        c.append(count(S, T) - g(S))
    return c


def cmd_pairs(a):
    import l2_compose as L2
    P = pieces()
    W = np.array([p[1][3] for p in P]); H = np.array([p[1][4] for p in P]); m = np.array([p[1][5] for p in P])
    rng = np.random.default_rng(a.seed)
    evs = L2.event_set('481')
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    f = open(a.out, 'w')
    with ProcessPoolExecutor(a.procs) as ex:
        for ev in evs:
            n = ev['n']; S = float(ev['ante']) + a.ds; t0 = time.time()
            idx = [k for k in range(len(P)) if m[k] <= n - 3 and max(W[k], H[k]) <= S + 1e-12]
            e = np.array([-99 if x is None else x for x in ex.map(_single, [(k, S) for k in idx], chunksize=16)])
            ok = e > -99
            idx = np.array(idx)[ok]; e = e[ok]
            w, h = W[idx], H[idx]
            cand = []
            for corner, dis in (('UR', (w[:, None] + w[None, :] <= S + 1e-12) | (h[:, None] + h[None, :] <= S + 1e-12)),
                                ('LR', (w[:, None] + w[None, :] <= S + 1e-12))):
                I, J = np.nonzero(dis)
                cand.append((corner, I, J))
            allp = [(c, i, j) for c, I, J in cand for i, j in zip(I, J)]
            est = np.array([e[i] + e[j] for _, i, j in allp])
            need = n - g(S)
            order = np.argsort(-est, kind='stable')
            top = [allp[q] for q in order[:a.top]]
            rnd = [allp[q] for q in rng.choice(len(allp), size=min(a.rand, len(allp)), replace=False)]
            sel = [(x, 'top') for x in top] + [(x, 'rand') for x in rnd]
            jobs = [(n, S, int(idx[i]), int(idx[j]), c) for (c, i, j), _ in sel]
            R = list(ex.map(_joint, jobs, chunksize=4))
            for ((c, i, j), kind), r in zip(sel, R):
                ka, kb = int(idx[i]), int(idx[j])
                f.write(json.dumps(dict(n=n, S=S, need=need, kind=kind, corner=c, a=ka, b=kb,
                                        la=P[ka][2], lb=P[kb][2], eA=int(e[i]), eB=int(e[j]),
                                        wA=W[ka], hA=H[ka], wB=W[kb], hB=H[kb],
                                        gx=S - W[ka] - W[kb], gy=S - H[ka] - H[kb], joint=r)) + '\n')
            f.flush()
            print(f'n={n} S={S:.5f} need={need}: {len(idx)} pieces, max e {e.max()}, {len(allp)} disjoint pairs, '
                  f'top est {est[order[0]] if len(est) else None}, {len(jobs)} joint ({time.time() - t0:.0f} s)', flush=True)


def cmd_report(a):
    import collections
    pf = os.path.join(a.dir, 'period.jsonl')
    if os.path.exists(pf):
        R = [json.loads(l) for l in open(pf)]
        const = from1 = n_ok = 0
        diffs = collections.Counter()
        for r in R:
            e = [x for x in r['e'] if x is not None]
            if len(e) < 2:
                continue
            n_ok += 1
            const += len(set(e)) == 1
            from1 += len(set(e[1:])) <= 1
            diffs[tuple(x - e[0] for x in e)] += 1
        print(f'period: {n_ok} pieces; e constant over S0..S0+3: {const}; constant from S0+1 on: {from1}')
        for d, c in diffs.most_common(8):
            print(f'   e(S0+k) - e(S0) = {d}: {c}')
        bad = [r for r in R if len(set(x for x in r['e'] if x is not None)) > 1]
        sl = np.array([r['slack0'] for r in bad]); print(f'   non-constant: slack0 (S0 - max(W, H)) quartiles {np.percentile(sl, [0, 25, 50, 75, 100]).round(3) if len(sl) else "-"}')
    qf = os.path.join(a.dir, 'pairs.jsonl')
    if os.path.exists(qf):
        R = [json.loads(l) for l in open(qf)]
        for kind in ('top', 'rand'):
            X = [r for r in R if r['kind'] == kind]
            if not X:
                continue
            d = np.array([max(r['joint']) - r['eA'] - r['eB'] for r in X])
            sep = np.array([max(r['gx'], r['gy']) for r in X])
            print(f'\npairs [{kind}]: {len(X)}; delta = joint excess (better order) - (eA + eB)')
            print('   delta histogram:', dict(sorted(collections.Counter(d.tolist()).items())))
            for lo, hi in ((-1e9, 0.5), (0.5, 1), (1, 2), (2, 3), (3, 1e9)):
                s = (sep >= lo) & (sep < hi)
                if s.sum():
                    print(f'   sep in [{lo:g}, {hi:g}): {s.sum():5d} pairs, exact additive {np.mean(d[s] == 0):.1%}, '
                          f'delta>0 {np.mean(d[s] > 0):.1%}, delta<0 {np.mean(d[s] < 0):.1%}, range [{d[s].min()}, {d[s].max()}]')
            reach = [r for r in X if max(r['joint']) >= r['need']]
            fp = [r for r in X if r['eA'] + r['eB'] >= r['need'] and max(r['joint']) < r['need']]
            fn = [r for r in X if r['eA'] + r['eB'] < r['need'] and max(r['joint']) >= r['need']]
            print(f'   reach n: {len(reach)}; additive says reach but joint does not: {len(fp)}; joint reaches, additive '
                  f'says no: {len(fn)}')


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('period'); p.add_argument('--n', type=int, default=300); p.add_argument('--seed', type=int, default=1)
    p.add_argument('--procs', type=int, default=14); p.add_argument('--out', required=True)
    p = sp.add_parser('pairs'); p.add_argument('--ds', type=float, default=0.08); p.add_argument('--top', type=int, default=300)
    p.add_argument('--rand', type=int, default=300); p.add_argument('--seed', type=int, default=1)
    p.add_argument('--procs', type=int, default=14); p.add_argument('--out', required=True)
    p = sp.add_parser('report'); p.add_argument('dir')
    a = ap.parse_args()
    dict(period=cmd_period, pairs=cmd_pairs, report=cmd_report)[a.cmd](a)


if __name__ == '__main__':
    main()
