#!/usr/bin/env python3
"""t3_chain: the H-lemma (one wall-to-wall chain of T + the transverse chains that bound its
rise) as a RESTRICTED FARKAS DUAL, tested over the whole angle space at T = 3, n = 6.

Task tasks/t3-chain/README.md.  Nothing in search/ is modified; s6skel, s6local and bandcut_k
are imported.

Semantics: the wall-carrying LP of search/BANDCUT_K.md sec 1.1 = search/S6_SKELETON.md sec 3.1,

    variables  (delta, c_1..c_n) ,   every row of the form   a . c  -  delta  >=  b ,
      walls   lo-x(i):  +x_i >= P_i + delta            hi-x(i):  -x_i >= P_i - T + delta
              lo-y(i):  +y_i >= P_i + delta            hi-y(i):  -y_i >= P_i - T + delta
      pairs   (i,j,o,kind,s):  s * n_{o,kind} . (c_j - c_i)  >=  m_ij + delta
              n_{o,0} = (cos th_o, sin th_o),  n_{o,1} = (-sin th_o, cos th_o),
              m_ij = 1/2 + (|cos D| + |sin D|)/2,  D = th_j - th_i,   P_i = (|cos|+|sin|)/2 .

A Farkas certificate is w >= 0 with  sum_r w_r a_r = 0  and  sum_r w_r = 1; it proves
delta <= -sum_r w_r b_r for every configuration satisfying the rows in its support.  Restricting
the support to an "H" -- one wall-to-wall chain of T plus transverse rows -- is the H-lemma.
"""
import argparse
import itertools
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                     # noqa: E402
import s6local                                                    # noqa: E402
import bandcut_k as BK                                            # noqa: E402

TOL = 1e-9


# ===================================================================== the row system
def geom(z, n):
    X, Y, TH = np.array(z[1::3]), np.array(z[2::3]), np.array(z[3::3])
    C, S = np.cos(TH), np.sin(TH)
    P = 0.5 * (np.abs(C) + np.abs(S))
    return X, Y, TH, C, S, P


def all_rows(z, n, T, lvl, tol=1e-7):
    """every row of the LP that the configuration z SATISFIES at level lvl (= its own margin).

    -> (A, b, tags) with A the (nrows, 2n) coefficient matrix on the centres, b the rhs.
    A dual supported on any subset of these rows is valid for z.  Using the full set reproduces
    delta* exactly when z is the global optimum of its angle vector (adding rows that z satisfies
    cannot cut z off, and the leaf assignment's rows are all present)."""
    X, Y, TH, C, S, P = geom(z, n)
    A, b, tags = [], [], []
    for i in range(n):
        for (ax, v) in (('x', X[i]), ('y', Y[i])):
            off = 0 if ax == 'x' else n
            r = np.zeros(2 * n); r[off + i] = 1.0
            A.append(r); b.append(P[i]); tags.append(('lo', ax, i))
            r = np.zeros(2 * n); r[off + i] = -1.0
            A.append(r); b.append(P[i] - T); tags.append(('hi', ax, i))
    for i in range(n):
        for j in range(i + 1, n):
            D = TH[j] - TH[i]
            m = 0.5 + 0.5 * (abs(math.cos(D)) + abs(math.sin(D)))
            dd = np.array([X[j] - X[i], Y[j] - Y[i]])
            for o in (i, j):
                for kind, nrm in enumerate((np.array([C[o], S[o]]),
                                            np.array([-S[o], C[o]]))):
                    pr = float(nrm @ dd)
                    for sg in (1, -1):
                        if sg * pr - m >= lvl - tol:
                            r = np.zeros(2 * n)
                            r[j] += sg * nrm[0]; r[i] -= sg * nrm[0]
                            r[n + j] += sg * nrm[1]; r[n + i] -= sg * nrm[1]
                            A.append(r); b.append(m)
                            tags.append(('pair', i, j, o, kind, sg))
    return np.array(A), np.array(b), tags


def row_axis(tag, TH):
    """'x' or 'y': which global axis this row's normal is nearer.  Walls are their own axis."""
    if tag[0] in ('lo', 'hi'):
        return tag[1]
    _, i, j, o, kind, sg = tag
    c, s = math.cos(TH[o]), math.sin(TH[o])
    d0, d1 = ((c, s), (-s, c))[kind]
    return 'x' if abs(d0) > abs(d1) else 'y'


# ===================================================================== restricted dual LP
def dual_bound(A, b, keep, want_w=False):
    """min -sum w_r b_r over w >= 0 supported on `keep`, sum w_r a_r = 0, sum w_r = 1.
    -> (bound, w) with bound = +inf when no cancelling combination exists on that support."""
    from scipy.optimize import linprog
    idx = np.asarray(sorted(keep), dtype=int)
    if len(idx) == 0:
        return math.inf, None
    Ai, bi = A[idx], b[idx]
    Aeq = np.vstack([Ai.T, np.ones((1, len(idx)))])
    beq = np.zeros(Aeq.shape[0]); beq[-1] = 1.0
    res = linprog(-bi, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * len(idx),
                  method='highs',
                  options={'primal_feasibility_tolerance': 1e-11,
                           'dual_feasibility_tolerance': 1e-11})
    if not res.success:
        return math.inf, None
    w = np.zeros(len(b)); w[idx] = res.x
    return float(res.fun), (w if want_w else None)


# ===================================================================== chains
def dag(z, n, lvl, tol=1e-7):
    """-> (gx, dirx, gy, diry) of bandcut_k.normal_gaps, over ALL n squares (tilted or not)."""
    return BK.normal_gaps(np.asarray(z, dtype=float), n)


def chains_of(z, n, lvl, L, tol=1e-7):
    """every path on exactly L vertices of the x-DAG and of the y-DAG at level lvl.
    -> list of (axis, path)."""
    X, Y = np.array(z[1::3]), np.array(z[2::3])
    gx, dxr, gy, dyr = BK.normal_gaps(np.asarray(z, dtype=float), n)
    out = []
    for (g, dr, co, ax) in ((gx, dxr, X, 'x'), (gy, dyr, Y, 'y')):
        for ch in BK.all_chains(list(range(n)), g, dr, co, lvl - tol, L):
            out.append((ax, ch))
    return out


def link_rows(tags, TH, ax, i, j, coord):
    """indices of pair rows for the pair {i,j} whose normal is `ax`-type and which separate with
    the higher-`coord` square on the + side (i.e. rows usable as a link of an ax-chain i -> j)."""
    out = []
    for q, tg in enumerate(tags):
        if tg[0] != 'pair':
            continue
        _, a, bq, o, kind, sg = tg
        if {a, bq} != {i, j}:
            continue
        if row_axis(tg, TH) != ax:
            continue
        # the row reads  sg * n.(c_b - c_a) >= m + delta ; it is the link a->b when sg=+1, b->a
        # when sg=-1.  We want the link  i -> j.
        head = bq if sg > 0 else a
        tail = a if sg > 0 else bq
        if (tail, head) == (i, j):
            out.append(q)
    return out


MODES = ('C', 'CW', 'H', 'HP', 'F')


def h_support(tags, TH, ax, path, mode):
    """the H-restricted support for a main chain `path` along axis `ax`.  `tr` is the transverse
    axis.  Walls are container rows, not chain rows, so every mode but 'C' keeps all 4n of them.

    'C'   chain only: the chain's ax-links and the two ax-wall rows at its ends.
    'CW'  chain + every wall row (no pair row off the chain at all).
    'H'   THE H-LEMMA: chain + every wall row + every TRANSVERSE-type pair row.  Excluded: every
          ax-type (main-direction) pair row that is not a link of the chain.
    'HP'  H plus the ax-type pair rows between squares OF THE CHAIN (staircase closure).
    'F'   everything (reproduces delta* when the configuration is the global optimum)."""
    keep = set()
    for q, tg in enumerate(tags):
        if tg[0] in ('lo', 'hi'):
            if mode == 'C':
                if tg[1] == ax and ((tg[0] == 'lo' and tg[2] == path[0])
                                    or (tg[0] == 'hi' and tg[2] == path[-1])):
                    keep.add(q)
            else:
                keep.add(q)
    for s in range(len(path) - 1):
        keep |= set(link_rows(tags, TH, ax, path[s], path[s + 1], None))
    if mode in ('C', 'CW'):
        return keep
    if mode == 'F':
        return set(range(len(tags)))
    ps = set(path)
    tr = 'y' if ax == 'x' else 'x'
    for q, tg in enumerate(tags):
        if tg[0] != 'pair':
            continue
        a = row_axis(tg, TH)
        if a == tr:
            keep.add(q)
        elif mode == 'HP' and {tg[1], tg[2]} <= ps:
            keep.add(q)
    return keep


def best_h(z, n, T, lvl, modes=('H',), maxchains=400, rows=None):
    """-> dict mode -> (bound, axis, path): the best restricted dual bound over all wall-to-wall
    chains of T at level lvl."""
    A, b, tags = rows if rows is not None else all_rows(z, n, T, lvl)
    TH = np.array(z[3::3])
    chs = chains_of(z, n, lvl, int(round(T)))
    out = {'nchains': len(chs)}
    for mode in modes:
        best = (math.inf, None, None)
        for (ax, path) in chs[:maxchains]:
            keep = h_support(tags, TH, ax, path, mode)
            v, _ = dual_bound(A, b, keep)
            if v < best[0]:
                best = (v, ax, list(path))
        out[mode] = dict(bound=best[0], axis=best[1], path=best[2])
    return out


def full_dual(z, n, T, lvl, wtol=1e-6):
    """the unrestricted dual at z: value (= delta* when z is the global optimum) and support."""
    A, b, tags = all_rows(z, n, T, lvl)
    v, w = dual_bound(A, b, range(len(b)), want_w=True)
    TH = np.array(z[3::3])
    supp = [(float(w[q]), tags[q], row_axis(tags[q], TH)) for q in range(len(b)) if w[q] > wtol]
    supp.sort(key=lambda t: -t[0])
    return v, supp


# ===================================================================== delta* at one theta
def delta_star(theta, T, rng, limit=400, jitters=3, rho=0.08, extra=400, seeds=(), keep=1):
    """max over centres of the wall-carrying margin, by LP ascent from tiling starts with random
    label dealings, `extra` uniform random starts, and any `seeds` (X, Y) supplied.
    A LOWER bound on delta* (feasible points only).  -> (value, z, [top-`keep` z's])."""
    n = len(theta)
    F = BK.FixedBand(list(theta), T, n, None, 1e-6)
    pool = list(s6local.sampled_starts(n, T, rng, limit, jitters, rho))
    for _ in range(extra):
        pool.append((rng.uniform(0.4, T - 0.4, n), rng.uniform(0.4, T - 0.4, n)))
    pool += list(seeds)
    found = []
    for (X, Y) in pool:
        v, bx, by = F.ascent(np.array(X, dtype=float), np.array(Y, dtype=float))
        if bx is not None:
            found.append((v, bx.copy(), by.copy()))
    found.sort(key=lambda t: -t[0])
    out = []
    for (v, bx, by) in found[:max(keep, 1)]:
        z = np.zeros(1 + 3 * n)
        z[1::3], z[2::3], z[3::3] = bx, by, theta
        z[0] = v
        out.append(z)
    return found[0][0], out[0], out


def classify_support(supp, n, TH):
    """name the shape of a dual support: how many heavy rows, on which axis, chain or not."""
    heavy = [s for s in supp if s[0] >= 0.05]
    ax = {}
    for (w, tg, a) in supp:
        ax[a] = ax.get(a, 0.0) + w
    walls = sum(w for (w, tg, a) in supp if tg[0] in ('lo', 'hi'))
    pairs = sum(w for (w, tg, a) in supp if tg[0] == 'pair')
    sq = set()
    for (w, tg, a) in heavy:
        sq |= ({tg[2]} if tg[0] in ('lo', 'hi') else {tg[1], tg[2]})
    return dict(nsupp=len(supp), nheavy=len(heavy), wall=walls, pair=pairs,
                wx=ax.get('x', 0.0), wy=ax.get('y', 0.0), heavy_squares=sorted(sq),
                heavy_tags=[str(s[1]) for s in heavy])


# ===================================================================== angle-space families
def sample_theta(fam, rng, n=6):
    """-> (theta, label).  Families covering the whole of [0,90 deg)^6 at T = 3.
    Tilt = distance to 0 mod 90 deg; the sign of a tilt is the sign of norm_tilt."""
    R, d = rng, math.radians

    def band(lo, hi, k, sign=None):
        """k tilts drawn from [lo, hi] degrees with a random (or fixed) sign."""
        out = []
        for _ in range(k):
            t = d(R.uniform(lo, hi))
            s = sign if sign is not None else (1 if R.random() < 0.5 else -1)
            out.append(s * t)
        return out

    if fam.startswith('Z'):                     # j angles exactly 0, rest uniform generic
        j = int(fam[1:])
        th = [0.0] * j + list(R.uniform(0, math.pi / 2, n - j))
        return th, f'{fam} j={j}'
    if fam.startswith('near'):                  # all six inside a band of eps degrees
        eps = float(fam[4:])
        th = band(0, eps, n)
        return th, f'near eps={eps}'
    if fam.startswith('coh'):                   # all six small, ONE sign (the coherent cone)
        eps = float(fam[3:])
        th = band(0, eps, n, sign=1)
        return th, f'coherent eps={eps}'
    if fam.startswith('hole'):                  # k near-axis (<= eps), 6-k far (>= eps)
        k, eps = fam[4:].split('_')
        k, eps = int(k), float(eps)
        th = band(0, eps, k) + band(eps, 45, n - k)
        return th, f'hole k={k} eps={eps}'
    if fam.startswith('unif'):                  # all six at one common tilt
        t = float(fam[4:])
        return [d(t)] * n, f'uniform {t} deg'
    if fam.startswith('split'):                 # a at +t, 6-a at -t
        a, t = fam[5:].split('_')
        a, t = int(a), float(t)
        return [d(t)] * a + [-d(t)] * (n - a), f'split {a}/+{t} {6-a}/-{t}'
    if fam == 'generic':
        return list(R.uniform(0, math.pi / 2, n)), 'generic'
    if fam.startswith('bigrand'):               # all tilts >= eps, random signs
        eps = float(fam[7:])
        return band(eps, 45, n), f'bigrand eps={eps}'
    raise SystemExit('unknown family ' + fam)


DEFAULT_FAMS = (['Z3', 'Z4', 'Z5', 'Z6']
                + ['near0.5', 'near1', 'near2', 'near5', 'near10', 'near20']
                + ['coh0.5', 'coh1', 'coh2', 'coh5', 'coh10', 'coh20']
                + [f'hole{k}_{e}' for k in (3, 4) for e in (1, 2, 5, 10, 20)]
                + [f'hole{k}_{e}' for k in (0, 1, 2) for e in (1, 5, 20)]
                + ['unif1', 'unif5', 'unif10', 'unif20', 'unif30', 'unif40', 'unif44', 'unif45']
                + [f'split{a}_{t}' for a in (1, 2, 3) for t in (2, 5, 20, 45)]
                + ['generic', 'bigrand5', 'bigrand10', 'bigrand20', 'bigrand30'])

_CFG = None


def _init(cfg):
    global _CFG
    _CFG = cfg


def _job(arg):
    fam, seed = arg
    cfg = _CFG
    n, T = cfg['n'], cfg['T']
    rng = np.random.default_rng(seed)
    th, label = sample_theta(fam, rng, n)
    v, z, _ = delta_star(th, T, rng, limit=cfg['limit'], extra=cfg['extra'])
    rows = all_rows(z, n, T, v)
    r = best_h(z, n, T, v, modes=cfg['modes'], rows=rows)
    A, b, tags = rows
    fv, w = dual_bound(A, b, range(len(b)), want_w=True)
    TH = np.array(z[3::3])
    supp = [(float(w[q]), tags[q], row_axis(tags[q], TH)) for q in range(len(b))
            if w[q] > 1e-6]
    supp.sort(key=lambda t: -t[0])
    out = dict(fam=fam, label=label, seed=seed, theta=[float(t) for t in th],
               delta=float(v), fulldual=float(fv), nchains=r['nchains'],
               z=[float(q) for q in z],
               supp=[[s[0], list(map(_js, s[1])), s[2]] for s in supp],
               shape=classify_support(supp, n, TH))
    for m in cfg['modes']:
        out[m] = float(r[m]['bound'])
        out[m + '_path'] = r[m]['path']
        out[m + '_axis'] = r[m]['axis']
    return out


def _js(x):
    return int(x) if isinstance(x, (int, np.integer)) else x


def cmd_scan(a):
    import multiprocessing as mp
    fams = a.fams.split(',') if a.fams else DEFAULT_FAMS
    cfg = dict(n=a.n, T=a.T, limit=a.limit, extra=a.extra, modes=list(a.modes.split(',')))
    jobs = [(f, a.seed0 + 1009 * i + j) for i, f in enumerate(fams) for j in range(a.reps)]
    t0 = __import__('time').time()
    with open(a.out, 'w') as fh:
        if a.nproc > 1:
            with mp.Pool(a.nproc, initializer=_init, initargs=(cfg,)) as pool:
                for k, rec in enumerate(pool.imap_unordered(_job, jobs, chunksize=1)):
                    fh.write(json.dumps(rec) + '\n'); fh.flush()
                    if k % 20 == 0:
                        print(f'{k}/{len(jobs)}  {__import__("time").time()-t0:.0f}s',
                              flush=True)
        else:
            _init(cfg)
            for k, j in enumerate(jobs):
                rec = _job(j)
                fh.write(json.dumps(rec) + '\n'); fh.flush()
                print(k, rec['fam'], rec['delta'], {m: rec[m] for m in cfg['modes']}, flush=True)
    print('done', len(jobs), f'{__import__("time").time()-t0:.0f}s')


def cmd_report(a):
    recs = []
    for f in a.file:
        recs += [json.loads(l) for l in open(f)]
    modes = [m for m in MODES if m in recs[0]]
    fams = {}
    for r in recs:
        fams.setdefault(r['fam'], []).append(r)
    print(f'# {len(recs)} samples, {len(fams)} families.  "H<=0" = the H-restricted dual bound is')
    print('# <= 1e-9, i.e. the H-lemma certifies delta <= 0 at that sampled optimum.')
    print(f'# {"family":14s} {"N":>3s} {"max delta*":>12s} {"median H":>12s} {"max H":>12s} '
          f'{"H<=0":>7s} {"HP<=0":>7s} {"CW<=0":>7s} {"nochain":>7s}')
    for fam in sorted(fams, key=lambda f: (recs and 0)):
        rs = fams[fam]
        H = sorted(r['H'] for r in rs)
        nz = sum(1 for r in rs if r['nchains'] == 0)
        ok = sum(1 for r in rs if r['H'] <= 1e-9)
        okp = sum(1 for r in rs if r.get('HP', 9) <= 1e-9)
        okw = sum(1 for r in rs if r.get('CW', 9) <= 1e-9)
        print(f'  {fam:14s} {len(rs):3d} {max(r["delta"] for r in rs):+12.4e} '
              f'{H[len(H)//2]:+12.4e} {max(H):+12.4e} {ok:3d}/{len(rs):<3d} '
              f'{okp:3d}/{len(rs):<3d} {okw:3d}/{len(rs):<3d} {nz:3d}')
    bad = [r for r in recs if r['H'] > 1e-9]
    print(f'\n# H fails at {len(bad)} / {len(recs)} sampled optima.')
    if bad:
        print('# worst 15 by H - delta*:')
        bad.sort(key=lambda r: -(r['H'] - r['delta']))
        for r in bad[:15]:
            tl = ' '.join(f'{math.degrees(s6skel.norm_tilt(t)):+.1f}' for t in r['theta'])
            print(f'  {r["fam"]:12s} d*={r["delta"]:+.4e} H={r["H"]:+.4e} '
                  f'nch={r["nchains"]:2d} tilts[{tl}]')
            print('      dual: ' + ' '.join(f'{w:.3f}{tg}' for w, tg, ax in r['supp'][:6]))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('scan')
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--T', type=float, default=3.0)
    p.add_argument('--fams', default='')
    p.add_argument('--reps', type=int, default=8)
    p.add_argument('--limit', type=int, default=300)
    p.add_argument('--extra', type=int, default=400)
    p.add_argument('--modes', default='CW,H,HP')
    p.add_argument('--seed0', type=int, default=12345)
    p.add_argument('--nproc', type=int, default=8)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_scan)
    p = sub.add_parser('report')
    p.add_argument('file', nargs='+')
    p.set_defaults(fn=cmd_report)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
