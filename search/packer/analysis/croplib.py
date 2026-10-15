#!/usr/bin/env python3
"""Crop library for blind Mondrian v1 composition (L2; 10-15).  The atlas keeps only corner crops that beat the grid on
their own (excess >= 1), and #481's pieces mostly don't (s(11)_r1a's 8-square crop, s(50)_r5's 36-square crop ...):
their value is in how the tilted squares meet the upright fill.  Relaxing the excess floor alone explodes (1.07M crops
over the 1,553 sources allowed at the #481 cutoff).  This library instead keeps crops that look like pieces:

  * every orientation (8) and up to R columns / rows dropped at the wall side (259 / 292 use one);
  * tight corner cut (growing either threshold adds a square), m in [3, nmax];
  * no tilted square cut by the crop's boundary (its bounding box either inside the box or clear of it): a piece is a
    whole sub-structure, not a slice through one;
  * >= --min-tilted tilted squares inside; excess m - floor(W) floor(H) >= --min-excess;
  * Pareto prune per (m, floor W, floor H) globally across sources (a smaller box with the same count and floors is never
    worse for the composer's count formula).

  python3 analysis/croplib.py build --set 481 --out runs/mondrian/croplib481.npz [--R 2 --min-tilted 3 --min-excess -3]
  python3 analysis/croplib.py cover --lib runs/mondrian/croplib481.npz        # are #481's named pieces in it?
"""
import argparse, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
EPS = 1e-9


def crops(sq, s, R=2, min_tilted=3, min_excess=-3, nmax=400):
    """-> W, H, m, o, i, j arrays (W >= H not enforced: orientation o in 0..7 carries the transpose)."""
    import mondrian as M
    out = []
    for o in range(8):
        q = M.orient(sq, s, o)
        h = M.half_ext(q[:, 2])
        til = M.is_tilted(q[:, 2])
        for i in range(R + 1):
            for j in range(R + 1):
                x, y = q[:, 0] - i, q[:, 1] - j
                lo_x, lo_y, hi_x, hi_y = x - h, y - h, x + h, y + h
                # a tilted square cut by the dropped lines at the wall: then no crop from this (i, j) is clean near it;
                # handled per candidate below through the generic partial-overlap test
                ins = (lo_x >= -EPS) & (lo_y >= -EPS)
                if ins.sum() < 3:
                    continue
                u, v = np.round(hi_x[ins], 9), np.round(hi_y[ins], 9)
                us, vs = np.unique(u), np.unique(v)
                C = np.zeros((len(us), len(vs)), np.int32)
                np.add.at(C, (np.searchsorted(us, u), np.searchsorted(vs, v)), 1)
                Mc = C.cumsum(0).cumsum(1)
                Ct = np.zeros_like(C)
                np.add.at(Ct, (np.searchsorted(us, u[til[ins]]), np.searchsorted(vs, v[til[ins]])), 1)
                Mt = Ct.cumsum(0).cumsum(1)
                Mu = np.vstack([np.zeros((1, len(vs)), np.int32), Mc[:-1]])
                Mv = np.hstack([np.zeros((len(us), 1), np.int32), Mc[:, :-1]])
                ex = Mc - (np.floor(us + EPS)[:, None] * np.floor(vs + EPS)[None, :]).astype(np.int64)
                ok = (Mc > Mu) & (Mc > Mv) & (Mc >= 3) & (Mc <= nmax) & (Mt >= min_tilted) & (ex >= min_excess)
                I, J = np.nonzero(ok)
                if len(I) == 0:
                    continue
                W, H = us[I], vs[J]
                # partial overlap of any tilted square's bbox with [0, W] x [0, H]
                tx0, tx1, ty0, ty1 = lo_x[til], hi_x[til], lo_y[til], hi_y[til]
                clean = np.ones(len(I), bool)
                for a in range(0, len(I), 4096):
                    w, hh = W[a:a + 4096, None], H[a:a + 4096, None]
                    meet = (tx0 < w - EPS) & (tx1 > EPS) & (ty0 < hh - EPS) & (ty1 > EPS)
                    inside = (tx0 >= -EPS) & (tx1 <= w + EPS) & (ty0 >= -EPS) & (ty1 <= hh + EPS)
                    clean[a:a + 4096] = ~(meet & ~inside).any(1)
                k = np.nonzero(clean)[0]
                n_ = len(k)
                out.append((W[k], H[k], Mc[I[k], J[k]], np.full(n_, o), np.full(n_, i), np.full(n_, j)))
    if not out:
        return [np.zeros(0)] * 6
    return [np.concatenate(c) for c in zip(*out)]


def prune(W, H, m, *rest):
    """global Pareto per (m, floor W, floor H) after canonicalising W >= H (the composer tries both transposes)."""
    import mondrian as M
    sw = W < H
    W2, H2 = np.where(sw, H, W), np.where(sw, W, H)
    W2, H2, m, W, H, *rest = M.pareto_prune(W2, H2, m, W, H, *rest)
    return [W, H, m, *rest]


def _job(args):
    k, sq, s, kw = args
    r = crops(np.asarray(sq, float), float(s), **kw)
    if len(r[0]) == 0:
        return k, None
    r = prune(*r)          # per source first (cheap), globally after
    return k, r


def build(a):
    import l2_screen as LS
    import l2_compose as L2
    from multiprocessing import Pool
    z = dict(np.load(a.atlas, allow_pickle=True))
    if a.set == 'all':
        A = range(len(z['labels']))
    else:
        ev = L2.event_set(a.set)[0]
        A = sorted(LS.allowed_sids(z, ev['T_utc'], ev['T_local']))
    kw = dict(R=a.R, min_tilted=a.min_tilted, min_excess=a.min_excess, nmax=a.nmax)
    t0 = time.time()
    parts = []
    with Pool(a.procs) as P:
        for k, r in P.imap_unordered(_job, [(k, z['sq'][k], z['src_s'][k], kw) for k in A], chunksize=4):
            if r is not None:
                parts.append(tuple(r) + (np.full(len(r[0]), k),))
    Wo, Ho, m, o, i, j, sid = (np.concatenate(c) for c in zip(*parts))
    n0 = len(Wo)
    Wo, Ho, m, o, i, j, sid = prune(Wo, Ho, m, o, i, j, sid)
    W = Wo
    print(f'{len(A)} sources -> {n0} crops (per-source pruned) -> {len(W)} (global Pareto)  {time.time() - t0:.0f} s', flush=True)
    np.savez(a.out, W=Wo, H=Ho, m=m, o=o, i=i, j=j, sid=sid, atlas=a.atlas)
    import collections
    c = collections.Counter(np.minimum(m // 25 * 25, 300))
    print('by m (bin 25):', dict(sorted(c.items())))


def cover(a):
    """for each #481 lineage piece (named source, m): the named crops (mondrian.crop_pieces, R) vs the library: is
    there a library piece with the same m, the same floors and a box no larger (in either transpose)?"""
    import mondrian as M
    L = dict(np.load(a.lib, allow_pickle=True))
    z = dict(np.load(str(L['atlas']), allow_pickle=True))
    lab = [str(x).split('@')[0] for x in z['labels']]
    W = np.maximum(L['W'], L['H']); H = np.minimum(L['W'], L['H']); m = L['m']
    for n, (S, pa, pb) in sorted(M.LINEAGE_481.items()):
        row = []
        for src, mm in (pa, pb):
            P = M.crop_pieces(M.load_src(src), mm, a.R, 0)
            best = None
            for T, d in P:
                w, hh = max(d[3], d[4]), min(d[3], d[4])
                hit = (m == mm) & (W <= w + 1e-9) & (H <= hh + 1e-9) & (np.floor(W + EPS) == np.floor(w + EPS)) & \
                      (np.floor(H + EPS) == np.floor(hh + EPS))
                if hit.any():
                    k = np.nonzero(hit)[0][0]
                    best = f"yes ({lab[L['sid'][k]]})"; break
            name = src.split('/')[-1].replace('-source.cert.txt', '')
            row.append(f"{name}:{mm} {best or 'NO'} [{len(P)} named crops]")
        print(n, ' | '.join(row))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    b = sp.add_parser('build'); b.add_argument('--set', default='481'); b.add_argument('--out', required=True)
    b.add_argument('--atlas', default=os.path.join(PK, 'runs', 'mondrian', 'atlas_v3.npz'))
    b.add_argument('--R', type=int, default=2); b.add_argument('--min-tilted', type=int, default=3)
    b.add_argument('--min-excess', type=int, default=-3); b.add_argument('--nmax', type=int, default=330)
    b.add_argument('--procs', type=int, default=16)
    c = sp.add_parser('cover'); c.add_argument('--lib', required=True); c.add_argument('--R', type=int, default=2)
    a = ap.parse_args()
    {'build': build, 'cover': cover}[a.cmd](a)
