#!/usr/bin/env python3
"""s(110) landscape explainer data (10-08).

  build.py reps    --minima M.json --dist D.npy --k 24 --out runs/landscape110
  build.py barrier --out runs/landscape110 --procs 15        # all pairs among the reps
  build.py cycle   --out runs/landscape110                   # tour + stored leg paths

reps: k-medoids on the corner-movement matrix (certified minima only), both records forced in.
barrier(A, B): smallest side s (bisection, res 2e-4) at which morph.walk_carrot connects A to B (spread A to s, follow
the straight-line morph with feasible LP steps, shrink to B): an upper bound on the minimal expansion for A -> B.
"""
import argparse, json, os, sys, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import morph
from jobpool import run_jobs

FORCE = ('10.99678339663', '10.99679327395')          # Couzo's record, our earlier record


def kmedoids(D, k, force, iters=100, seed=0):
    force = list(dict.fromkeys(force))
    rng = np.random.default_rng(seed)
    n = len(D)
    med = list(force) + [int(i) for i in rng.choice([i for i in range(n) if i not in force], k - len(force), replace=False)]
    for _ in range(iters):
        lab = np.argmin(D[:, med], axis=1)
        new = list(force)
        for c in range(len(force), k):
            mem = np.where(lab == c)[0]
            if len(mem) == 0:
                new.append(med[c]); continue
            new.append(int(mem[np.argmin(D[np.ix_(mem, mem)].sum(1))]))
        if new == med:
            break
        med = new
    return med


def cmd_reps(a):
    M = json.load(open(a.minima)); D = np.load(a.dist)
    ok = json.load(open(a.ok)) if a.ok else None
    idx = [i for i in range(len(M)) if ok is None or ok[i]]
    Dk = D[np.ix_(idx, idx)]
    force = [idx.index(next(i for i in idx if M[i]['S'].startswith(f))) for f in FORCE]
    if a.keep:                                       # grow an existing set: its members stay medoids
        old = [r['census_index'] for r in json.load(open(a.keep))]
        force = [idx.index(i) for i in old]
    med = kmedoids(Dk, a.k, force)
    reps = sorted((idx[m] for m in med), key=lambda i: float(M[i]['S']))
    os.makedirs(a.out, exist_ok=True)
    json.dump([dict(M[i], census_index=i) for i in reps], open(f'{a.out}/reps.json', 'w'), indent=0)
    print(len(reps), 'reps:', [M[i]['S'][:12] for i in reps])


def _barrier(job):
    i, j, pa, pb = job
    A, B = morph.load(pa), morph.load(pb)
    perm = morph.match2(A, B)
    base = max(A[0], B[0])
    lo, hi, h = base, None, 0.01
    while h <= 0.32:
        if morph.walk_carrot(base + h, A, B, perm)[0]:
            hi = base + h; break
        lo = base + h; h *= 2
    if hi is None:
        return dict(i=i, j=j, s=None)
    while hi - lo > 2e-4:
        mid = (lo + hi) / 2
        if morph.walk_carrot(mid, A, B, perm)[0]:
            hi = mid
        else:
            lo = mid
    return dict(i=i, j=j, s=hi, lo=lo, base=base)


def _barrier_rev(job):
    i, j, pa, pb = job                                 # walk j -> i; reported under (i, j) with rev = True
    r = _barrier((j, i, pb, pa))
    return dict(r, i=i, j=j, rev=True)


def cmd_barrier_rev(a):
    """retry the failed pairs walking the other way (the walk is not symmetric); merged into barriers.json."""
    R = json.load(open(f'{a.out}/reps.json'))
    B = json.load(open(f'{a.out}/barriers.json'))
    fail = [b for b in B if b['s'] is None and not b.get('rev_tried')]
    res = run_jobs(_barrier_rev, [(b['i'], b['j'], R[b['i']]['path'], R[b['j']]['path']) for b in fail], procs=a.procs,
                   timeout=3600, on_timeout=lambda j: None, on_error=lambda j, x: None)
    got = {(r['i'], r['j']): dict(r, S_i=R[r['i']]['S'], S_j=R[r['j']]['S']) for r in res if r and r['s'] is not None}
    B = [got.get((b['i'], b['j']), dict(b, rev_tried=True) if b['s'] is None else b) for b in B]
    json.dump(B, open(f'{a.out}/barriers.json', 'w'), indent=0)
    print(f'{len(got)} of {len(fail)} failed pairs connect the other way')


def cmd_barrier(a):
    R = json.load(open(f'{a.out}/reps.json'))
    pairs = [(i, j) for i in range(len(R)) for j in range(i + 1, len(R))]
    if a.near:                                       # only each rep's `near` nearest reps by corner movement
        D = np.load(a.dist)
        ci = [r['census_index'] for r in R]
        nb = {i: set(np.argsort([D[ci[i], ci[j]] if j != i else 1e9 for j in range(len(R))])[:a.near].tolist()) for i in range(len(R))}
        pairs = [(i, j) for i, j in pairs if j in nb[i] or i in nb[j]]
    old = {}
    if os.path.exists(f'{a.out}/barriers.json'):
        for b in json.load(open(f'{a.out}/barriers.json')):
            old[(b['S_i'], b['S_j'])] = b
    key = lambda i, j: (R[i]['S'], R[j]['S'])
    jobs = [(i, j, R[i]['path'], R[j]['path']) for i, j in pairs if key(i, j) not in old]
    print(f'{len(pairs)} pairs, {len(pairs) - len(jobs)} already computed, {len(jobs)} to run', flush=True)
    t0 = time.time()
    res = run_jobs(_barrier, jobs, procs=a.procs, timeout=3600, on_timeout=lambda j: None, on_error=lambda j, x: None)
    res = [dict(r, S_i=R[r['i']]['S'], S_j=R[r['j']]['S']) for r in res if r]
    keep = [dict(b, i=[k for k in range(len(R)) if R[k]['S'] == b['S_i']][0], j=[k for k in range(len(R)) if R[k]['S'] == b['S_j']][0])
            for b in old.values() if any(R[k]['S'] == b['S_i'] for k in range(len(R))) and any(R[k]['S'] == b['S_j'] for k in range(len(R)))]
    res = keep + res
    json.dump(res, open(f'{a.out}/barriers.json', 'w'), indent=0)
    print(f'{len(res)}/{len(jobs)} pairs in {time.time() - t0:.0f} s; failed {sum(r["s"] is None for r in res)}')


def tour(W):
    """closed tour over all nodes: nearest neighbour from every start + 2-opt; cost = sum W + 3 * max W (keeps the
    worst leg low as well as the total)."""
    n = len(W)
    cost = lambda t: sum(W[t[k], t[(k + 1) % n]] for k in range(n)) + 3 * max(W[t[k], t[(k + 1) % n]] for k in range(n))
    best = None
    for st in range(n):
        t, left = [st], set(range(n)) - {st}
        while left:
            j = min(left, key=lambda j: W[t[-1], j]); t.append(j); left.remove(j)
        imp = True
        while imp:
            imp = False
            for a in range(1, n - 1):
                for b in range(a + 1, n):
                    t2 = t[:a] + t[a:b + 1][::-1] + t[b + 1:]
                    if cost(t2) < cost(t) - 1e-12:
                        t, imp = t2, True
        if best is None or cost(t) < cost(best):
            best = t
    return best


def thin(path, k):
    """<= k frames, evenly spaced in cumulative movement (first and last kept)."""
    if len(path) <= k:
        return path
    d = [0.0]
    for p, q in zip(path, path[1:]):
        d.append(d[-1] + morph.dist(p, q))
    tg = np.linspace(0, d[-1], k)
    idx = sorted(set(int(np.searchsorted(d, x)) for x in tg) | {0, len(path) - 1})
    return [path[min(i, len(path) - 1)] for i in idx]


def _leg(job):
    pa, pb, s, flip = job
    if flip:
        r = _leg((pb, pa, s, False))
        return r and dict(r, frames=r['frames'][::-1], reversed=True)
    A, B = morph.load(pa), morph.load(pb)
    perm = morph.match2(A, B)
    ok, d, path = morph.walk_carrot(s, A, B, perm)
    if not ok:                                         # bisection's last success was at this s; retry slightly above
        s = s + 2e-4
        ok, d, path = morph.walk_carrot(s, A, B, perm)
    fa, fb = s / A[0], s / B[0]
    frames = []
    for f in np.linspace(1, fa, 6):                    # spread A to fill the box
        frames.append((A[0] * f, [(x * f, y * f, t) for x, y, t in A[1]]))
    frames += [(s, q) for q in thin(path, 40)[1:]]
    tg = morph.target(A, B, perm)
    for f in np.linspace(fb, 1, 6)[1:]:                # shrink to B (angles continuous with the walk)
        frames.append((B[0] * f, [(x * f, y * f, t) for x, y, t in tg]))
    return dict(ok=bool(ok), s=float(s), perm=[int(p) for p in perm], frames=frames)


def cmd_cycle(a):
    R = json.load(open(f'{a.out}/reps.json'))
    n = len(R)
    S = [float(r['S']) for r in R]
    W = np.full((n, n), 10.0)
    Sb = np.full((n, n), np.nan)
    fwd = {}                                           # (a, b) -> True if the walk that worked goes a -> b
    for r in json.load(open(f'{a.out}/barriers.json')):
        if r['s'] is not None:
            i, j = r['i'], r['j']
            Sb[i, j] = Sb[j, i] = r['s']
            W[i, j] = W[j, i] = r['s'] - max(S[i], S[j])
            a_, b_ = (j, i) if r.get('rev') else (i, j)
            fwd[(a_, b_)] = True
    np.fill_diagonal(W, 0)
    deg = (W < 10).sum(1) - 1
    keep = [i for i in range(n) if deg[i] >= 2]
    if len(keep) < n:
        print('left out of the tour (fewer than 2 connections):', [i for i in range(n) if i not in keep])
    t = [keep[k] for k in tour(W[np.ix_(keep, keep)])]
    n = len(t)
    legs = [(t[k], t[(k + 1) % n]) for k in range(n)]
    print('tour', t, 'excess per leg', [round(float(W[i, j]), 4) for i, j in legs])
    res = run_jobs(_leg, [(R[i]['path'], R[j]['path'], float(Sb[i, j]), (i, j) not in fwd) for i, j in legs], procs=a.procs, timeout=1800,
                   on_timeout=lambda j: None, on_error=lambda j, x: None)
    out = []
    for (i, j), r in zip(legs, res):
        out.append(dict(i=i, j=j, s_barrier=float(Sb[i, j]), **(r or dict(ok=False))))
        print(i, j, 'ok' if r and r['ok'] else 'FAILED', r and len(r['frames']))
    json.dump(dict(tour=t, W=W.tolist(), Sb=np.nan_to_num(Sb).tolist(), legs=out), open(f'{a.out}/cycle.json', 'w'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'runs/landscape110'))
    ap.add_argument('--minima'); ap.add_argument('--dist'); ap.add_argument('--ok'); ap.add_argument('--k', type=int, default=24)
    ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--keep'); ap.add_argument('--near', type=int, default=0)
    a = ap.parse_args()
    {'reps': cmd_reps, 'barrier': cmd_barrier, 'barrier-rev': cmd_barrier_rev, 'cycle': cmd_cycle}[a.cmd](a)
