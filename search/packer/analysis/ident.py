"""How many distinct minima are in a pool, and can a cheap test tell copies apart from new minima?  (10-20; Evan:
"only keeping the novel stuff ... dedupe can be a hard problem at this scale when a proper polish is compute intensive")

  python3 analysis/ident.py cluster --n 274 --smax 17.01              # cheap: side window + fingerprint clusters
  python3 analysis/ident.py calibrate --n 274 --smax 17.01 --pairs 40 # expensive check on a sample: long polish both

Cheap test: a row is a copy of an earlier (lower) row if its side is within W above it and the fingerprint distance
(pk/fingerprint: smoothed centre images, min over 8 symmetries, L1 / n) is <= tau.  Unconverged copies of one minimum sit
up to ~1e-3 above it (pk/identity.py), hence the one-sided window.
Expensive check: nearmiss's long fq polish (loosen 1.0, no ALM, flip-top 64, pit 600) on both rows; same minimum if the
polished sides agree to DISTINCT (1e-7, nearmiss.py's rule).  Not proof either way (a polish can stall), but it is the
store's current working definition.
"""
from __future__ import annotations
import argparse, json, os, random, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from pk import fingerprint as FP

W = 1e-3
TAUS = (0.002, 0.005, 0.01, 0.02, 0.04)


def load(n, smax, status=('polished', 'certified', 'weak', 'not-min', 'unresolved')):
    import sqlite3
    db = sqlite3.connect('file:runs/store.sqlite?mode=ro', uri=True)
    q = (f"SELECT id, s, coords FROM packing WHERE n=? AND s<? AND status IN ({','.join('?' * len(status))}) "
         "AND (grid IS NULL OR grid = 0) ORDER BY s")
    R = db.execute(q, (n, smax, *status)).fetchall()
    ids = np.array([r[0] for r in R]); S = np.array([r[1] for r in R])
    SQ = [np.frombuffer(r[2], '<f8').reshape(-1, 3) for r in R]
    return ids, S, SQ


def images(S, SQ):
    return np.stack([FP.image(s, sq) for s, sq in zip(S, SQ)])


def nearest_lower(S, IM, n, w=W):
    """for each row: (fingerprint distance, index) of the nearest row with side in [s - w, s) (min over symmetries)"""
    out = np.full(len(S), np.inf); arg = np.full(len(S), -1)
    j0 = 0
    for i in range(len(S)):
        while S[j0] < S[i] - w:
            j0 += 1
        if j0 < i:
            d = FP.nearest(IM[i], IM[j0:i], n)
            k = int(np.argmin(d)); out[i] = d[k]; arg[i] = j0 + k
    return out, arg


def cluster(a):
    ids, S, SQ = load(a.n, a.smax)
    n = SQ[0].shape[0]
    IM = images(S, SQ)
    d, arg = nearest_lower(S, IM, n)
    keys = len(set(np.round(S / 2e-9).astype(np.int64).tolist()))
    print(f'n {a.n}: {len(S):,} polished non-grid rows below {a.smax} (best {S[0]:.9f}); distinct side keys (2e-9) {keys:,}')
    print('nearest lower neighbour (side window 1e-3) fingerprint distance quantiles 10/25/50/75/90: '
          + ' / '.join(f'{v:.4f}' for v in np.quantile(d[np.isfinite(d)], [.1, .25, .5, .75, .9])))
    for t in TAUS:
        print(f'  tau {t:<6g}: {np.sum(d > t):>7,} rows are not a copy of a lower row  ({np.mean(d > t) * 100:4.1f} %)')
    np.savez(f'runs/ident_{a.n}.npz', ids=ids, S=S, d=d, arg=arg)


def _polish(args):
    s, sq = args
    from pk.engine import QOpt, fq_quench
    from pk.nearmiss import RECIPES
    with tempfile.TemporaryDirectory() as t:
        q = fq_quench((s, sq), QOpt.from_args(RECIPES['polish'] + ['--seed', '1']), tmp=t, timeout=7200)
        try:
            d = json.loads(q.stdout.strip().splitlines()[-1])
        except Exception:
            return None
        ok = d.get('min_gap', 0) > -1e-9 and d.get('min_wall', 0) > -1e-9
        return float(d['s']) if ok else None


def calibrate(a):
    from concurrent.futures import ProcessPoolExecutor
    ids, S, SQ = load(a.n, a.smax)
    z = np.load(f'runs/ident_{a.n}.npz'); d, arg = z['d'], z['arg']
    rng = random.Random(a.seed)
    # stratified by fingerprint distance: pairs (row, its nearest lower row)
    edges = (0, 0.002, 0.005, 0.01, 0.02, 0.04, np.inf)
    P = []
    for lo, hi in zip(edges, edges[1:]):
        idx = [i for i in range(len(S)) if lo <= d[i] < hi and arg[i] >= 0]
        P += [(i, int(arg[i]), lo) for i in rng.sample(idx, min(a.pairs, len(idx)))]
    rows = sorted({k for i, j, _ in P for k in (i, j)})
    with ProcessPoolExecutor(a.procs) as ex:
        pol = dict(zip(rows, ex.map(_polish, [(float(S[k]), SQ[k]) for k in rows])))
    out = []
    print(f'{len(P)} pairs, {len(rows)} polishes')
    for lo, hi in zip(edges, edges[1:]):
        B = [(i, j) for i, j, l in P if l == lo]
        ok = [(i, j) for i, j in B if pol[i] is not None and pol[j] is not None]
        same = sum(abs(pol[i] - pol[j]) < 1e-7 for i, j in ok)
        drop = [S[i] - pol[i] for i, j in ok]
        print(f'  fp dist [{lo:g}, {hi:g}): {len(B)} pairs, polished both {len(ok)}, same minimum {same} '
              f'({same / max(1, len(ok)) * 100:.0f} %); polish moved the upper row by median {np.median(drop) if drop else 0:.2e}')
        out += [dict(i=int(ids[i]), j=int(ids[j]), d=float(d[i]), si=float(S[i]), sj=float(S[j]), pi=pol[i], pj=pol[j])
                for i, j in B]
    json.dump(out, open(f'runs/ident_{a.n}_cal.json', 'w'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    for c in ('cluster', 'calibrate'):
        p = sp.add_parser(c); p.add_argument('--n', type=int, required=True); p.add_argument('--smax', type=float, required=True)
        p.add_argument('--pairs', type=int, default=20); p.add_argument('--procs', type=int, default=15)
        p.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    dict(cluster=cluster, calibrate=calibrate)[a.cmd](a)
