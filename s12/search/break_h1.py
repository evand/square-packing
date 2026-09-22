#!/usr/bin/env python3
"""break_h1: adversarial hunt for a configuration where the H1 restricted dual FAILS, i.e.

    H1(z) > 0   while   delta*(z) < 0 ,

at T = 4, n = 12, plus the ladder-depth / margin-exponent measurement at T = 3, 4, 5.

Task tasks/break-h1/README.md.  NOTHING in search/ is modified: s6skel, s6local, bandcut_k,
t3_chain, t3_chain_shape and t4_cycles are imported.

Semantics and row system exactly as search/t3_chain.py (= BANDCUT_K.md sec 1.1): the
wall-carrying LP in the 2n centre coordinates at fixed angles, every row  a . c - delta >= b.
A restricted dual (Farkas) bound on a support S is

    min  -sum_r w_r b_r   over  w >= 0 supported on S, sum_r w_r a_r = 0, sum_r w_r = 1 ,

= +inf when no cancelling combination exists on S.  Every such bound is >= delta*, so it can
never be negative at a margin-0 point (the filter).

The ladder of supports, each one main-direction row further out than the last:

    CW  chain of T + every wall row                          (empirically O(eps))
    H   + every transverse pair row      (the H-lemma)       (empirically O(eps^2))
    H1  + ONE further main-direction pair row                (empirically O(eps^3))
    H2  + TWO further main-direction pair rows               (O(eps^4)?)

"certificate depth" d(z) = the smallest level of that ladder whose bound is <= 0.
"""
import argparse
import itertools
import json
import math
import os
import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                     # noqa: E402
import s6local                                                    # noqa: E402
import t3_chain as TC                                             # noqa: E402
import t3_chain_shape as TS                                       # noqa: E402
import t4_cycles as TCY                                           # noqa: E402
import bandcut_k as BK                                            # noqa: E402


# ===================================================================== the ladder of bounds
def h1_h2(A, b, tags, TH, chains, maxchains=400, nch2=2, cap2=4000):
    """H1 (best over chains x one extra main-direction row) and H2 (best over pairs of extra
    rows, on the `nch2` chains that did best at H1 / H -- a RESTRICTED search, so the H2 value
    reported is an upper bound on the true best-pair value; that is the conservative direction
    for a `H2 fails` claim and the optimistic one for `H2 holds`).

    -> dict(H, H1, H2, H1_extra, H2_extra, H_path, H1_path, H2_path, nlp)
    """
    nlp = 0
    perch = []                       # (Hval, H1val, ax, path, base, extras, bestq)
    for (ax, path) in chains[:maxchains]:
        base = TC.h_support(tags, TH, ax, path, 'H')
        hv, _ = TC.dual_bound(A, b, base); nlp += 1
        extras = [q for q, tg in enumerate(tags)
                  if tg[0] == 'pair' and q not in base and TC.row_axis(tg, TH) == ax]
        best1, bq = math.inf, None
        for q in extras:
            v, _ = TC.dual_bound(A, b, base | {q}); nlp += 1
            if v < best1:
                best1, bq = v, q
        perch.append((hv, best1, ax, list(path), base, extras, bq))
    if not perch:
        return dict(H=math.inf, H1=math.inf, H2=math.inf, H_path=None, H1_path=None,
                    H2_path=None, H1_extra=None, H2_extra=None, nlp=nlp)
    bh = min(perch, key=lambda t: t[0])
    b1 = min(perch, key=lambda t: t[1])
    # H2: all pairs of extra rows, on a couple of promising chains
    if nch2 <= 0:
        return dict(H=bh[0], H1=b1[1], H2=math.inf, nlp=nlp,
                    H_path=(bh[2], bh[3]), H1_path=(b1[2], b1[3]),
                    H1_extra=None if b1[6] is None else list(map(TC._js, tags[b1[6]])),
                    H2_path=None, H2_extra=None)
    cands = []
    for t in sorted(perch, key=lambda t: t[1])[:nch2]:
        if t not in cands:
            cands.append(t)
    if bh not in cands:
        cands.append(bh)
    best2, b2 = math.inf, None
    for (hv, h1v, ax, path, base, extras, bq) in cands:
        npairs = len(extras) * (len(extras) - 1) // 2
        ex = extras
        if npairs > cap2:            # keep the best `m` singles, m(m-1)/2 <= cap2
            sc = []
            for q in extras:
                v, _ = TC.dual_bound(A, b, base | {q}); nlp += 1
                sc.append((v, q))
            sc.sort()
            m = int((1 + math.sqrt(1 + 8 * cap2)) / 2)
            ex = [q for _, q in sc[:m]]
        for q1, q2 in itertools.combinations(ex, 2):
            v, _ = TC.dual_bound(A, b, base | {q1, q2}); nlp += 1
            if v < best2:
                best2, b2 = v, (ax, list(path), q1, q2)
    return dict(H=bh[0], H1=b1[1], H2=best2, nlp=nlp,
                H_path=(bh[2], bh[3]), H1_path=(b1[2], b1[3]),
                H1_extra=None if b1[6] is None else list(map(TC._js, tags[b1[6]])),
                H2_path=None if b2 is None else (b2[0], b2[1]),
                H2_extra=None if b2 is None else [list(map(TC._js, tags[b2[2]])),
                                                 list(map(TC._js, tags[b2[3]]))])


def weight_levels(w, ratio, wtol=1e-13, band=3.0):
    """cluster the dual weights into geometric levels of ratio `ratio` (= tan eps).

    -> (depth, levels) with levels = [(level index, count, mean weight)] and depth the number of
    CONSECUTIVE levels 0, 1, ... that are populated."""
    ws = np.asarray([q for q in w if q > wtol], dtype=float)
    if len(ws) == 0 or ratio <= 0 or ratio >= 1:
        return 0, []
    wmax = ws.max()
    lg = np.log(ws / wmax) / np.log(ratio)
    lev = np.round(lg).astype(int)
    ok = np.abs(lg - lev) < 0.5 * band       # band = 3 -> accept anything within ratio^1.5
    out = defaultdict(list)
    for q, l, o in zip(ws, lev, ok):
        out[int(l) if o else -999].append(float(q))
    levels = sorted((l, len(v), float(np.mean(v))) for l, v in out.items() if l != -999)
    depth = 0
    have = {l for l, _, _ in levels}
    while depth in have:
        depth += 1
    return depth, levels


def bounds(z, n, T, lvl, tol=1e-7, maxchains=400, kmax=9, do_cyc=True, eps_rad=None,
           nch2=2, cap2=4000, do_h2=True):
    """every bound of the ladder at the configuration z, whose own margin is `lvl`."""
    t0 = time.time()
    A, b, tags = TC.all_rows(z, n, T, lvl, tol=tol)
    TH = np.array(z[3::3])
    chains = TC.chains_of(z, n, lvl, int(round(T)), tol=tol)
    out = {'nchains': len(chains), 'nrows': len(b)}
    for mode in ('CW', 'HP'):
        best = (math.inf, None)
        for (ax, path) in chains[:maxchains]:
            v, _ = TC.dual_bound(A, b, TC.h_support(tags, TH, ax, path, mode))
            if v < best[0]:
                best = (v, (ax, list(path)))
        out[mode] = best[0]
        out[mode + '_path'] = best[1]
    r = h1_h2(A, b, tags, TH, chains, maxchains=maxchains,
              nch2=(nch2 if do_h2 else 0), cap2=cap2)
    out.update({k: r[k] for k in ('H', 'H1', 'H2', 'H_path', 'H1_path', 'H2_path',
                                  'H1_extra', 'H2_extra')})
    out['nlp'] = r['nlp']
    if do_cyc:
        tight, _ = TCY.tight_indices(A, b, z, lvl, tol=tol)
        pg = TCY.tight_pair_graph(tags, tight)
        cv, cyc, cw, ncyc = TCY.best_cycle_bound(A, b, tags, pg, kmax=kmax)
        out['CYC'] = cv
        out['ncycles'] = ncyc
    fv, w = TC.dual_bound(A, b, range(len(b)), want_w=True)
    out['F'] = fv
    supp = [(float(w[q]), tags[q], TC.row_axis(tags[q], TH)) for q in range(len(b))
            if w[q] > 1e-13]
    supp.sort(key=lambda t: -t[0])
    out['supp'] = [[s[0], list(map(TC._js, s[1])), s[2]] for s in supp]
    g = TS.support_graph(supp, wtol=1e-9)
    out['mu'] = g['mu']; out['E'] = g['nedge']; out['V'] = g['nvert']
    out['cls'] = TCY.classify(dict(g, walls={f'{k[0]}-{k[1]}': float(v)
                                             for k, v in g['walls'].items()}))
    if eps_rad:
        d, levels = weight_levels([s[0] for s in supp], math.tan(eps_rad))
        out['wdepth'] = d
        out['wlevels'] = levels
    # certificate depth: the first level of the ladder that reaches <= 0
    dep = None
    for d, key in ((1, 'CW'), (2, 'H'), (3, 'H1'), (4, 'H2')):
        if out.get(key, math.inf) <= 1e-12:
            dep = d
            break
    out['cdepth'] = dep
    out['secs'] = time.time() - t0
    return out


# ===================================================================== structured angle sets
def perm_hole_cells(T, perm):
    """the T^2 tiling cell centres minus the permutation hole set {(i, perm[i])}."""
    Ti = int(round(T))
    holes = {(i, perm[i]) for i in range(Ti)}
    return [(a + 0.5, b + 0.5) for a in range(Ti) for b in range(Ti) if (a, b) not in holes]


def tilt_subsets(T, cells, kind):
    """-> list of indices into `cells` to tilt, for the named structured choice."""
    Ti = int(round(T))
    lo, hi = (Ti - 2) / 2.0, Ti / 2.0                # the central 2x2 block of cell centres
    lo3, hi3 = (Ti - 3) / 2.0, (Ti + 1) / 2.0
    idx = {}
    idx['central2'] = [q for q, (x, y) in enumerate(cells) if lo < x < hi + 1 and lo < y < hi + 1]
    idx['central3'] = [q for q, (x, y) in enumerate(cells)
                       if lo3 < x < hi3 + 1 and lo3 < y < hi3 + 1]
    idx['diag'] = [q for q, (x, y) in enumerate(cells) if abs(x - y) < 1e-9]
    idx['antidiag'] = [q for q, (x, y) in enumerate(cells) if abs(x + y - T) < 1e-9]
    idx['border'] = [q for q, (x, y) in enumerate(cells)
                     if x < 1 or y < 1 or x > T - 1 or y > T - 1]
    idx['all'] = list(range(len(cells)))
    return idx.get(kind, [])


def sign_pattern(cells, sub, kind):
    """+-1 per tilted square.  'coh' = all +, 'alt' = checkerboard, 'pin' = pinwheel swirl
    (sign by which quadrant of the block the square is in, rotationally)."""
    out = []
    for q in sub:
        x, y = cells[q]
        if kind == 'coh':
            out.append(1.0)
        elif kind == 'neg':
            out.append(-1.0)
        elif kind == 'alt':
            out.append(1.0 if (int(x) + int(y)) % 2 == 0 else -1.0)
        elif kind == 'pin':
            out.append(1.0 if (int(x) - int(y)) % 2 == 0 else -1.0)
        else:
            out.append(1.0)
    return out


def structured_T4(eps, n=12, T=4.0):
    """every structured (theta, seed centres, label) at tilt eps, T = 4, n = 12.

    1. the 4x4 tiling minus each permutation hole set, with a structured subset of the 12
       remaining squares tilted by +-eps (coherent / alternating / pinwheel signs);
    2. central 2x2 and 3x3 pinwheels at tilt eps on the tiling-minus-hole-set;
    3. the (4,1) band stack perturbed (added by the caller, it needs a file).
    """
    Ti = int(round(T))
    out = []
    for perm in itertools.permutations(range(Ti)):
        cells = perm_hole_cells(T, perm)
        if len(cells) != n:
            continue
        for tk in ('central2', 'central3', 'diag', 'antidiag', 'border', 'all'):
            sub = tilt_subsets(T, cells, tk)
            if not sub:
                continue
            for sk in ('coh', 'alt', 'pin'):
                sg = sign_pattern(cells, sub, sk)
                if sk != 'coh' and len(set(sg)) == 1:
                    continue
                th = [0.0] * n
                for q, s in zip(sub, sg):
                    th[q] = s * eps
                X = np.array([c[0] for c in cells])
                Y = np.array([c[1] for c in cells])
                out.append((th, (X, Y), f'perm{"".join(map(str, perm))}_{tk}_{sk}'))
    return out


def _d4_canon(cellpos_tilt, Ti):
    """canonical form of a (cells, tilt) structure under the dihedral group of the TxT grid.
    `cellpos_tilt` = frozenset of ((col, row), signed tilt in units of eps)."""
    best = None
    for (fx, fy, tr) in itertools.product((0, 1), (0, 1), (0, 1)):
        for sgn in (1, -1):
            t = []
            for ((c, r), s) in cellpos_tilt:
                cc, rr = (Ti - 1 - c if fx else c), (Ti - 1 - r if fy else r)
                # a reflection reverses the sign of a tilt, a transposition does too
                sg = s * sgn * ((-1) ** (fx + fy + tr))
                t.append(((rr, cc) if tr else (cc, rr), sg))
            key = tuple(sorted(t))
            if best is None or key < best:
                best = key
    return best


def dedup_theta(items):
    """drop repeats of the same angle MULTISET-with-position-independent structure: two starts
    whose sorted tilt vector agrees and whose seed centres agree up to the dihedral group."""
    seen, out = set(), []
    for (th, seed, lab) in items:
        X, Y = seed
        Ti = int(round(max(max(X), max(Y)) + 0.5))
        mx = max(abs(q) for q in th) or 1.0
        st = frozenset(((int(x), int(y)), int(round(t / mx)))
                       for x, y, t in zip(X, Y, th))
        key = _d4_canon(st, Ti)
        if key in seen:
            continue
        seen.add(key)
        out.append((th, seed, lab))
    return out


# ============================================ the top-of-hole chain-free family at any T
def merged_cells(T, a, bq):
    """the map from the T x T tiling cells to the (T-1)^2 LEVEL grid obtained by merging the
    adjacent tiling columns a, a+1 and the adjacent tiling rows bq, bq+1.

    -> dict (level cell) -> [tiling cells mapping to it]."""
    Ti = int(round(T))
    out = defaultdict(list)
    for c in range(Ti):
        for r in range(Ti):
            lc = c if c <= a else (a if c == a + 1 else c - 1)
            mr = r if r <= bq else (bq if r == bq + 1 else r - 1)
            out[(lc, mr)].append((c, r))
    return out


def tophole_family(T, eps, rng, nsamp, drop=False, sgn_mix=False):
    """the `BANDCUT_K.md` sec 3(b) construction, written out for general `T`.

    A chain-free configuration of `n = T(T-1)` squares needs its `k = (T-1)^2` chain-cut squares
    to sit ONE PER CELL of a (T-1) x (T-1) level grid.  Take that grid to be the tiling's, with
    one adjacent column pair and one adjacent row pair merged (verified to be the level map of
    the recorded `T = 4`, `k = 8, 9` optimum).  Then

      * the merged-column x merged-row level cell covers FOUR tiling cells: put a pinwheel of
        four squares at tilt eps on them, one of which is the cell's cut-set representative;
      * each of the other `2(T-2)` cross level cells covers TWO tiling cells: one carries an
        axis-parallel cut-set square, the other is a HOLE;
      * the `(T-2)^2` remaining level cells cover one tiling cell each: an axis-parallel
        cut-set square.

    Counting: cut-set squares `(T-1)^2`, far squares `T - 1` (three of the pinwheel plus
    `T - 4` more placed on removed cells), holes `T`, total `T(T-1) = n`.  At `T = 4` this is
    exactly the recorded optimum (holes a permutation set, central 2x2 pinwheel); at `T = 5` the
    holes are NO LONGER a permutation set -- see the note.

    `drop=True` drops one cut-set square to a far square, giving `k = (T-1)^2 - 1`.
    -> list of (theta, (X, Y), pattern, perm, label)."""
    Ti = int(round(T))
    n = Ti * (Ti - 1)
    out = []
    for _ in range(nsamp):
        a, bq = int(rng.integers(Ti - 1)), int(rng.integers(Ti - 1))
        mc = merged_cells(T, a, bq)
        centre = (a, bq)
        pin = mc[centre]                                   # the four central tiling cells
        surv = {}
        for key, cs in mc.items():
            if key == centre:
                surv[key] = pin[int(rng.integers(4))]
            else:
                surv[key] = cs[int(rng.integers(len(cs)))]
        removed = [c for key, cs in mc.items() for c in cs if c != surv[key]]
        cen = [c for c in pin if c != surv[centre]]
        others = [c for c in removed if c not in pin]
        nfar = Ti - 1                                  # = n - k far squares
        if nfar <= len(cen):
            # T = 3: only `nfar` of the three central removed cells carry a square; the rest
            # are holes, so the central block is a pinwheel of `nfar + 1` and not of four.
            far = [cen[q] for q in rng.permutation(len(cen))[:nfar]]
        else:
            nextra = nfar - len(cen)
            if nextra > len(others):
                continue
            far = cen + [others[q] for q in rng.permutation(len(others))[:nextra]]
        holes = [c for c in removed if c not in far]
        if len(holes) != Ti:
            continue
        cut = [(key, surv[key]) for key in sorted(mc)]
        dropped = None
        if drop:
            cand = [q for q in range(len(cut)) if cut[q][1] not in pin]
            dropped = cut[cand[int(rng.integers(len(cand)))]]
            cut = [c for c in cut if c != dropped]
            far = far + [dropped[1]]
        pattern = [list(key) for (key, _c) in cut]
        cells = [c for (_k, c) in cut] + list(far)
        if len(cells) != n or len(set(cells)) != n:
            continue
        theta, X, Y = [], [], []
        gs = 1.0 if (not sgn_mix or rng.random() < 0.5) else -1.0
        for (c, r) in cells:
            X.append(c + 0.5); Y.append(r + 0.5)
            if (c, r) in set(far) | ({surv[centre]} if surv[centre] in pin else set()):
                sg = gs if not sgn_mix else (1.0 if rng.random() < 0.5 else -1.0)
                theta.append(sg * eps)
            else:
                theta.append(0.0)
        lab = f'tophole T{Ti} m{a}{bq} k{len(cut)}'
        out.append((theta, (np.array(X), np.array(Y)), pattern, list(range(n)), lab))
    return out


def _th_init(cfg):
    _HG['cfg'] = cfg


def _th_job(arg):
    (theta, seedXY, pattern, perm, lab, seed) = arg
    cfg = _HG['cfg']
    T, n = cfg['T'], len(theta)
    k = len(pattern)
    rng = np.random.default_rng(seed)
    X, Y = seedXY
    seeds = [(X, Y)]
    for _ in range(cfg['njit']):
        seeds.append((X + rng.normal(0, cfg['rho'], n), Y + rng.normal(0, cfg['rho'], n)))
    v, z, hit = delta_star_at(theta, T, rng, seeds=seeds, limit=cfg['limit'], jitters=2,
                              extra=cfg['extra'], pattern=[tuple(c) for c in pattern], k=k)
    if z is None:
        return dict(lab=lab, err='infeasible')
    return dict(lab=lab, k=k, eps=cfg['epsdeg'], delta=float(v), hit=int(hit),
                cfree=bool(chainfree_ok(z, n, T, v, range(k))),
                pattern=pattern, z=[float(q) for q in z])


def band_solve(z0, n, T, k, eps, pattern, eta=1e-6, cycles=3, rounds=8, maxiter=250):
    """the full bandcut_k pipeline at one structured start: LP in the centres (angles fixed),
    then SLSQP over centres AND angles with the separating-axis assignment fixed, alternately.

    The angles matter: at the recorded `T = 4` optimum the eight "axis-parallel" squares carry a
    tilt of `1e-4 deg`, and it is exactly that tilt that stops the diagonal pair of the middle
    x-level from separating.  A fixed-angle search cannot find the cell."""
    pr = BK.BandProblem(n, T, k, eps, BK.start_signs(z0, n, k), pattern=pattern, eta=eta)
    best = (-1e18, None)
    z, v = BK.lp_polish(z0, n, T, k, pattern, eta)
    if v > best[0]:
        best = (v, z.copy())
    for _c in range(cycles):
        v, z = BK.solve_band(pr, z, rounds=rounds, maxiter=maxiter)
        if z is None:
            break
        if v > best[0]:
            best = (v, z.copy())
        z, v = BK.lp_polish(z, n, T, k, pattern, eta)
        if v > best[0]:
            best = (v, z.copy())
    v, z = best
    if z is None:
        return None, None, None
    z = np.array(z); z[0] = v
    ok, info = BK.verify(z, n, T, k, eps, v, need_chainfree=True)
    return (v if ok else None), z, info


def _th2_job(arg):
    (theta, seedXY, pattern, lab, seed, epsdeg, carry) = arg
    cfg = _HG['cfg']
    T = cfg['T']
    n = len(theta)
    k = len(pattern)
    eps = math.radians(epsdeg)
    rng = np.random.default_rng(seed)
    pat = [tuple(c) for c in pattern]
    starts = []
    if carry is not None:
        z0 = BK.clamp_band(np.array(carry, dtype=float), n, k, eps)
        starts.append(z0)
        for (rho, al) in ((0.0, 0.0), (0.03, 0.3 * eps), (0.1, 0.6 * eps)):
            for _ in range(2):
                starts.append(BK.nudge(z0, rng, n, T, k, eps, rho, al))
    else:
        X, Y = seedXY
        for q in range(cfg['njit'] + 1):
            z = np.zeros(1 + 3 * n)
            z[1::3] = X + (0 if q == 0 else rng.normal(0, cfg['rho'], n))
            z[2::3] = Y + (0 if q == 0 else rng.normal(0, cfg['rho'], n))
            z[3::3] = [t + (0 if q == 0 else rng.normal(0, 0.15 * eps)) for t in theta]
            for i in range(n):
                z[3 + 3 * i] = (min(max(z[3 + 3 * i], -eps), eps) if i < k else
                                (1 if z[3 + 3 * i] >= 0 else -1) * min(
                                    max(abs(z[3 + 3 * i]), eps), math.pi / 4))
            starts.append(z)
    best = (-1e18, None, None)
    for z0 in starts:
        v, z, info = band_solve(z0, n, T, k, eps, pat, eta=cfg['eta'], cycles=cfg['cycles'])
        if v is not None and v > best[0]:
            best = (v, z, info)
    if best[1] is None:
        return dict(lab=lab, eps=epsdeg, k=k, err='no chain-free optimum')
    v, z, info = best
    return dict(lab=lab, eps=epsdeg, k=k, delta=float(v), pattern=pattern,
                cfree=True, info={q: (list(w) if isinstance(w, list) else w)
                                  for q, w in info.items()},
                z=[float(q) for q in z])


def _th2_star(pk):
    return _th2_job(pk)


def cmd_tophole2(a):
    """the top-of-hole chain-free cell, structured starts + the FULL bandcut_k pipeline (angles
    free inside the band), with continuation in eps from the largest down."""
    import multiprocessing as mp
    T = a.T
    rng = np.random.default_rng(a.seed)
    epss = sorted([float(s) for s in a.eps.split(',')], reverse=True)
    fam = []
    for drop in ((False, True) if a.drop else (False,)):
        for mix in ((False, True) if a.mix else (False,)):
            fam += tophole_family(T, math.radians(epss[0]), rng, a.nsamp,
                                  drop=drop, sgn_mix=mix)
    seen, ded = set(), []
    for (t, sd, pt, pm, lb) in fam:
        key = (tuple(round(q, 12) for q in t),
               tuple(sorted((round(x, 6), round(y, 6)) for x, y in zip(sd[0], sd[1]))))
        if key in seen:
            continue
        seen.add(key); ded.append((t, sd, pt, pm, lb))
    print(f'# {len(fam)} samples, {len(ded)} distinct structures', flush=True)
    cfg = dict(T=T, njit=a.njit, rho=a.rho, eta=a.eta, cycles=a.cycles)
    out = open(a.out, 'w')
    carry = {}
    best = {}
    with mp.Pool(a.nproc, initializer=_th_init, initargs=(cfg,)) as pool:
        for epsdeg in epss:
            eps = math.radians(epsdeg)
            jobs = []
            for i, (t, sd, pt, pm, lb) in enumerate(ded):
                th = [(eps if q > 0 else (-eps if q < 0 else 0.0)) for q in t]
                lab = lb.split(' e')[0]
                jobs.append((th, sd, pt, lab, a.seed + 7919 * i, epsdeg,
                             carry.get((len(pt), i))))
            print(f'#  eps={epsdeg}: {len(jobs)} jobs', flush=True)
            res = []
            for i, rec in enumerate(pool.imap(_th2_star, jobs)):
                out.write(json.dumps(rec) + '\n'); out.flush()
                res.append((i, rec))
                if 'delta' in rec:
                    key = (rec['k'], epsdeg)
                    if key not in best or rec['delta'] > best[key]['delta']:
                        best[key] = rec
            for (i, rec) in res:
                if 'delta' in rec:
                    carry[(rec['k'], i)] = rec['z']
            for key in sorted(k for k in best if k[1] == epsdeg):
                r = best[key]
                e = math.radians(epsdeg)
                nh = sum(1 for (_i, q) in res
                         if 'delta' in q and q['k'] == key[0] and q['delta'] >= r['delta'] - 1e-9)
                print(f"   k={key[0]} eps={epsdeg:g}: delta* = {r['delta']:+.8e}  "
                      f"/eps^2 {-r['delta']/e**2:.5f}  /eps^3 {-r['delta']/e**3:.4f}  "
                      f"/eps^4 {-r['delta']/e**4:.3f}   hits {nh}/{len(res)}   {r['lab']}",
                      flush=True)
    out.close()


def cmd_tophole(a):
    """the top-of-hole chain-free cell at T = 4 or 5, built from the structured family above."""
    import multiprocessing as mp
    T = a.T
    Ti = int(round(T))
    jobs = []
    rng = np.random.default_rng(a.seed)
    for epsdeg in [float(s) for s in a.eps.split(',')]:
        eps = math.radians(epsdeg)
        fam = []
        for drop in ((False, True) if a.drop else (False,)):
            for mix in ((False, True) if a.mix else (False,)):
                fam += [(t, s, p, pm, lb + (' mix' if mix else '') + f' e{epsdeg:g}')
                        for (t, s, p, pm, lb) in
                        tophole_family(T, eps, rng, a.nsamp, drop=drop, sgn_mix=mix)]
        seen, ded = set(), []
        for f in fam:
            key = (f[4].split(' e')[0], tuple(round(q, 12) for q in f[0]),
                   tuple(sorted((round(x, 6), round(y, 6)) for x, y in zip(f[1][0], f[1][1]))))
            if key in seen:
                continue
            seen.add(key); ded.append(f)
        print(f'#  eps={epsdeg}: {len(fam)} samples, {len(ded)} distinct', flush=True)
        for i, (t, s, p, pm, lb) in enumerate(ded):
            jobs.append((t, s, p, pm, lb, a.seed + 7919 * i + int(epsdeg * 101)))
    cfg = dict(T=T, limit=a.limit, extra=a.extra, njit=a.njit, rho=a.rho, epsdeg=0.0)
    print(f'# {len(jobs)} jobs', flush=True)
    best = {}
    with mp.Pool(a.nproc, initializer=_th_init, initargs=(cfg,)) as pool, open(a.out, 'w') as f:
        for i, rec in enumerate(pool.imap_unordered(_th_star, [(j, float(j[4].split(' e')[-1]))
                                                               for j in jobs])):
            f.write(json.dumps(rec) + '\n'); f.flush()
            if 'delta' not in rec or not rec['cfree']:
                continue
            key = (rec['k'], rec['eps'])
            if key not in best or rec['delta'] > best[key]['delta']:
                best[key] = rec
                print(f"  [{i+1}/{len(jobs)}] NEW BEST {rec['lab']:30s} k={rec['k']} "
                      f"delta*={rec['delta']:+.6e}", flush=True)
    print('\n# best chain-free configuration per (k, eps):', flush=True)
    for key in sorted(best):
        r = best[key]
        print(f"  k={key[0]} eps={key[1]:g}  delta* = {r['delta']:+.8e}  "
              f"/eps^2 {r['delta']/math.radians(key[1])**2:+.5f}  "
              f"/eps^3 {r['delta']/math.radians(key[1])**3:+.5f}  "
              f"/eps^4 {r['delta']/math.radians(key[1])**4:+.5f}   {r['lab']}", flush=True)


def _th_star(pk):
    j, e = pk
    _HG['cfg'] = dict(_HG['cfg'], epsdeg=e)
    return _th_job(j)


# ===================================================================== delta* at fixed angles
def delta_star_at(theta, T, rng, seeds=(), limit=200, jitters=2, rho=0.08, extra=200,
                  pattern=None, k=0, eta=1e-6, keep=1):
    """max margin over centres at the fixed angle vector, by LP ascent from tiling starts, the
    supplied structured seeds, and uniform random starts.  A LOWER bound on delta*.

    With `pattern` given, the chain rows of bandcut_k.FixedBand forbid a chain of T among the
    first `k` squares (relabel the angle vector to choose WHICH squares those are)."""
    n = len(theta)
    F = BK.FixedBand(list(theta), T, k, pattern, eta)
    pool = list(seeds)
    pool += list(s6local.sampled_starts(n, T, rng, limit, jitters, rho))
    for _ in range(extra):
        pool.append((rng.uniform(0.4, T - 0.4, n), rng.uniform(0.4, T - 0.4, n)))
    found = []
    for (X, Y) in pool:
        v, bx, by = F.ascent(np.array(X, dtype=float), np.array(Y, dtype=float))
        if bx is not None:
            found.append((v, bx.copy(), by.copy()))
    if not found:
        return None, None, 0
    found.sort(key=lambda t: -t[0])
    hit = sum(1 for f in found if f[0] > found[0][0] - 1e-9)
    z = np.zeros(1 + 3 * n)
    z[1::3], z[2::3], z[3::3] = found[0][1], found[0][2], theta
    z[0] = found[0][0]
    return found[0][0], z, hit


def chainfree_ok(z, n, T, lvl, near, tol=1e-7):
    """re-check chain-freeness from the geometry alone: no path on T vertices of the x- or
    y-DAG among the squares `near`, at the configuration's own margin."""
    X, Y = np.array(z[1::3]), np.array(z[2::3])
    gx, dxr, gy, dyr = BK.normal_gaps(np.asarray(z, dtype=float), n)
    for (g, dr, co) in ((gx, dxr, X), (gy, dyr, Y)):
        h, _p = BK.longest_path(list(near), g, dr, co, lvl - tol)
        if h >= int(round(T)):
            return False
    return True


# ===================================================================== structured hunt (T = 4)
def merge_level_map(cellpos, S, T):
    """every chain-free level pattern of the cut set S obtained from the tiling by MERGING one
    adjacent pair of tiling columns and one adjacent pair of tiling rows.

    That is the combinatorial form of the T = 4 top-of-hole chain-free optimum of
    `BANDCUT_K.md` sec 3(b): its near-axis squares sit on tiling cells whose columns 1, 2 and
    whose rows 1, 2 are NOT x- (resp. y-) separated, so the x- and y-heights drop from T to
    T - 1.  Verified on that configuration (the map below is injective there).

    -> list of (pattern, (a, b)) with pattern[q] the level of S[q]."""
    Ti = int(round(T))
    out = []
    for a in range(Ti - 1):
        for bq in range(Ti - 1):
            pat = []
            for i in S:
                c, r = cellpos[i]
                lc = c if c <= a else (a if c == a + 1 else c - 1)
                mr = r if r <= bq else (bq if r == bq + 1 else r - 1)
                pat.append((lc, mr))
            if len(set(pat)) != len(pat):
                continue
            if max(q[0] for q in pat) > Ti - 2 or max(q[1] for q in pat) > Ti - 2:
                continue
            out.append(([list(q) for q in pat], (a, bq)))
    return out


def cutset_jobs(rng, theta, cellpos, T, n, ks, npat, patterns_by_k, nrand=1):
    """(k, pattern, perm, ntilted, label) jobs.  `perm` relabels the squares so that the
    chain-cut set is perm[:k] -- which may contain TILTED squares, unlike bandcut_k's scan,
    where the cut set is by construction the near-axis band.

    Cut sets, in order of how close they are to `BANDCUT_K.md` sec 3(b)'s optimum:
      the whole near-axis set;  near-axis + one tilted square (its `k = 9` labelling);
      near-axis + two tilted;  near-axis minus one square (the k below the top);
      and `nrand` random ones per k as a control.
    Patterns come from merge_level_map (one adjacent column pair and one adjacent row pair of
    the tiling merged) and, for the random cut sets, from bandcut_k.level_patterns."""
    L2 = (int(round(T)) - 1) ** 2
    out = []
    tilted = [i for i in range(n) if abs(s6skel.norm_tilt(theta[i])) > 1e-9]
    flat = [i for i in range(n) if i not in tilted]
    Ss = []
    if len(flat) <= L2:
        Ss.append((list(flat), 0, 'near'))
        for t in tilted:
            if len(flat) + 1 <= L2:
                Ss.append(([t] + list(flat), 1, 'near+1til'))
        for t1, t2 in itertools.combinations(tilted, 2):
            if len(flat) + 2 <= L2:
                Ss.append(([t1, t2] + list(flat), 2, 'near+2til'))
        for d in flat:
            Ss.append(([q for q in flat if q != d], 0, 'near-1'))
    for _ in range(npat):
        k = int(rng.integers(min(ks), max(ks) + 1))
        ntil = int(rng.integers(0, min(len(tilted), k) + 1))
        if k - ntil > len(flat):
            continue
        S = [int(q) for q in rng.permutation(tilted)[:ntil]] + \
            [int(q) for q in rng.permutation(flat)[:k - ntil]]
        Ss.append((S, ntil, 'mix'))
    seen = set()
    for (S, ntil, nm) in Ss:
        if len(S) not in ks or tuple(sorted(S)) in seen:
            continue
        seen.add(tuple(sorted(S)))
        rest = [i for i in range(n) if i not in S]
        for (pat, ab) in merge_level_map(cellpos, S, T):
            out.append((len(S), pat, list(S) + rest, ntil, f'{nm}m{ab[0]}{ab[1]}'))
    for _ in range(nrand):
        k = int(rng.integers(min(ks), max(ks) + 1))
        pats = patterns_by_k[k]
        pat = pats[int(rng.integers(len(pats)))]
        ntil = int(rng.integers(0, min(len(tilted), k) + 1))
        if k - ntil > len(flat):
            continue
        S = [int(q) for q in rng.permutation(tilted)[:ntil]] + \
            [int(q) for q in rng.permutation(flat)[:k - ntil]]
        rest = [i for i in range(n) if i not in S]
        out.append((k, [list(c) for c in pat], S + rest, ntil, 'rnd'))
    return out


def _hunt_init(cfg):
    _HG['cfg'] = cfg


_HG = {}


def _hunt_job(arg):
    (idx, theta, seed, lab, seedXY, k, pat, perm, ntil) = arg
    cfg = _HG['cfg']
    T, n = cfg['T'], len(theta)
    rng = np.random.default_rng(seed)
    th = [theta[p] for p in perm] if perm is not None else list(theta)
    seeds = []
    if seedXY is not None:
        X, Y = seedXY
        if perm is not None:
            X, Y = np.array([X[p] for p in perm]), np.array([Y[p] for p in perm])
        seeds = [(X, Y)]
        for _ in range(cfg['njit']):
            seeds.append((X + rng.normal(0, 0.07, n), Y + rng.normal(0, 0.07, n)))
    try:
        v, z, hit = delta_star_at(th, T, rng, seeds=seeds, limit=cfg['limit'],
                                  jitters=2, extra=cfg['extra'],
                                  pattern=pat, k=(k or 0))
    except Exception as e:                                        # noqa: BLE001
        return dict(lab=lab, err=str(e))
    if z is None:
        return dict(lab=lab, err='no feasible start')
    rec = dict(lab=lab, k=k, ntil=ntil, delta=float(v), hit=int(hit), idx=idx,
               eps=cfg['epsdeg'], pat=None if pat is None else [list(c) for c in pat],
               perm=perm, z=[float(q) for q in z])
    if pat is not None:
        rec['cfree'] = bool(chainfree_ok(z, n, T, v, range(k)))
    if v > cfg['dmax']:
        rec['skip'] = 'delta >= 0'
        return rec
    o = bounds(z, n, T, v, tol=cfg['tol'], maxchains=cfg['maxchains'],
               eps_rad=math.radians(cfg['epsdeg']), kmax=cfg['kmax'],
               do_cyc=cfg['docyc'], nch2=cfg['nch2'], cap2=cfg['cap2'])
    for key in ('CW', 'H', 'HP', 'H1', 'H2', 'CYC', 'F', 'cdepth', 'wdepth', 'wlevels',
                'mu', 'E', 'V', 'cls', 'nchains', 'ncycles', 'H1_extra', 'secs'):
        if key in o:
            rec[key] = o[key]
    rec['supp'] = o['supp'][:14]
    return rec


def cmd_hunt(a):
    import multiprocessing as mp
    T, n = 4.0, 12
    patterns_by_k = {k: BK.level_patterns(T, k) for k in range(1, 10)}
    jobs = []
    idx = 0
    for epsdeg in [float(s) for s in a.eps.split(',')]:
        eps = math.radians(epsdeg)
        fams = dedup_theta(structured_T4(eps, n=n, T=T))
        if a.maxfam and len(fams) > a.maxfam:
            rr = np.random.default_rng(a.seed)
            fams = [fams[q] for q in rr.permutation(len(fams))[:a.maxfam]]
        for (th, seedXY, lab) in fams:
            rng = np.random.default_rng(a.seed + 7919 * idx)
            if a.free:
                jobs.append((idx, th, a.seed + idx, f'{lab}|e{epsdeg:g}|free', seedXY,
                             None, None, None, 0))
                idx += 1
            cellpos = {q: (int(seedXY[0][q]), int(seedXY[1][q])) for q in range(n)}
            cj = cutset_jobs(
                rng, th, cellpos, T, n, [int(s) for s in a.ks.split(',')],
                a.npat, patterns_by_k, nrand=a.nrand)
            if a.jpf and len(cj) > a.jpf:
                cj = [cj[q] for q in rng.permutation(len(cj))[:a.jpf]]
            for (k, pat, perm, ntil, pl) in cj:
                jobs.append((idx, th, a.seed + idx, f'{lab}|e{epsdeg:g}|k{k}t{ntil}{pl}',
                             seedXY, k, pat, perm, ntil))
                idx += 1
    cfg = dict(T=T, limit=a.limit, extra=a.extra, njit=a.njit, tol=a.tol,
               maxchains=a.maxchains, kmax=a.kmax, docyc=not a.nocyc, nch2=a.nch2,
               cap2=a.cap2, dmax=a.dmax, epsdeg=0.0)
    print(f'# {len(jobs)} jobs', flush=True)
    rng = np.random.default_rng(a.seed)
    if a.limitjobs and len(jobs) > a.limitjobs:
        jobs = [jobs[q] for q in rng.permutation(len(jobs))[:a.limitjobs]]
        print(f'# sampled down to {len(jobs)}', flush=True)
    # eps is per-job (encoded in the label); carry it by rebuilding cfg per job is overkill --
    # the label holds it and bounds() only uses it for the weight-level ratio, so pass it through
    out = open(a.out, 'w')
    nfail = nneg = 0
    packed = []
    for j in jobs:
        e = float(j[3].split('|e')[1].split('|')[0])
        packed.append((j, e))
    with mp.Pool(a.nproc, initializer=_hunt_init, initargs=(cfg,)) as pool:
        for i, rec in enumerate(pool.imap_unordered(_hunt_star, packed)):
            out.write(json.dumps(rec) + '\n'); out.flush()
            if 'H1' in rec:
                nneg += 1
                if rec['H1'] > 1e-9:
                    nfail += 1
                    print(f"  ** H1 FAILS {rec['lab']} delta*={rec['delta']:+.4e} "
                          f"H1={rec['H1']:+.4e} H2={rec['H2']:+.4e}", flush=True)
            if (i + 1) % 25 == 0:
                print(f'  [{i+1}/{len(jobs)}] delta<0: {nneg}, H1 failures: {nfail}', flush=True)
    print(f'# done: {nneg} optima with delta* < 0, {nfail} H1 failures', flush=True)
    out.close()


def _hunt_star(pk):
    j, e = pk
    _HG['cfg'] = dict(_HG['cfg'], epsdeg=e)
    return _hunt_job(j)


# ===================================================================== evaluate recorded optima
def _eval_init(cfg):
    _HG['cfg'] = cfg


def _eval_job(rec):
    cfg = _HG['cfg']
    z = np.array(rec['z'], dtype=float)
    n = (len(z) - 1) // 3
    T = cfg['T']
    lvl = float(rec['delta'])
    out = dict({k: rec[k] for k in rec if k != 'z'}, n=n)
    if lvl > cfg['dmax']:
        out['skip'] = 'delta >= 0'
    o = bounds(z, n, T, lvl, tol=cfg['tol'], maxchains=cfg['maxchains'],
               eps_rad=(math.radians(rec['eps']) if rec.get('eps') else None),
               kmax=cfg['kmax'], do_cyc=cfg['docyc'], nch2=cfg['nch2'], cap2=cfg['cap2'])
    for key in ('CW', 'H', 'HP', 'H1', 'H2', 'CYC', 'F', 'cdepth', 'wdepth', 'wlevels',
                'mu', 'E', 'V', 'cls', 'nchains', 'ncycles', 'H1_extra', 'H2_extra', 'secs'):
        if key in o:
            out[key] = o[key]
    out['supp'] = o['supp'][:16]
    out['z'] = [float(q) for q in z]
    return out


def cmd_eval(a):
    """bounds at every optimum recorded in a bandcut_k-style jsonl (or the T = 4 census)."""
    import multiprocessing as mp
    recs = []
    if a.census:
        for r in TCY.load_recorded(a.root):
            recs.append(dict(lab=r['cell'], src=r['src'], k=r['k'], eps=r['eps'],
                             delta=r['delta'], z=[float(q) for q in r['z']]))
    for fn in (a.file or []):
        for line in open(fn):
            r = json.loads(line)
            if r.get('value') is None or r.get('z') is None:
                continue
            recs.append(dict(lab=f"k={r['k']} eps={r['eps']:g} "
                             f"{'chainfree' if r['pidx'] >= 0 else 'free'} pat{r['pidx']}",
                             src=os.path.basename(fn), k=r['k'], eps=r['eps'], pidx=r['pidx'],
                             delta=float(r['value']), z=r['z'],
                             nhit=r.get('nhit'), nstart=r.get('nstart')))
    if a.only:
        recs = [r for r in recs if any(s in r['lab'] for s in a.only)]
    if a.maxn and len(recs) > a.maxn:
        rr = np.random.default_rng(a.seed)
        recs = [recs[q] for q in rr.permutation(len(recs))[:a.maxn]]
    cfg = dict(T=a.T, tol=a.tol, maxchains=a.maxchains, kmax=a.kmax, docyc=not a.nocyc,
               nch2=a.nch2, cap2=a.cap2, dmax=a.dmax)
    print(f'# {len(recs)} recorded optima', flush=True)
    nfail = nneg = nzero = nbug = 0
    with mp.Pool(a.nproc, initializer=_eval_init, initargs=(cfg,)) as pool, \
            open(a.out, 'w') as f:
        for i, rec in enumerate(pool.imap_unordered(_eval_job, recs)):
            f.write(json.dumps(rec) + '\n'); f.flush()
            if rec['delta'] > a.dmax:
                nzero += 1
                worst = min(rec.get(k, math.inf) for k in ('CW', 'H', 'HP', 'H1', 'H2', 'CYC'))
                if worst < -1e-8:
                    nbug += 1
                    print(f"  !! FILTER BREACH {rec['lab']} delta={rec['delta']:+.3e} "
                          f"worst bound {worst:+.4e}", flush=True)
            else:
                nneg += 1
                if rec.get('H1', -1) > 1e-9:
                    nfail += 1
                    print(f"  ** H1 FAILS {rec['lab']} delta*={rec['delta']:+.4e} "
                          f"H1={rec['H1']:+.4e} H2={rec.get('H2'):+.4e}", flush=True)
            if (i + 1) % 25 == 0:
                print(f'  [{i+1}/{len(recs)}] delta<0: {nneg} (H1 fails {nfail}), '
                      f'margin-0: {nzero} (filter breaches {nbug})', flush=True)
    print(f'# done: {nneg} with delta*<0, {nfail} H1 failures; {nzero} margin-0 points, '
          f'{nbug} filter breaches', flush=True)


# ===================================================================== reporting
def fit_exponent(eps_deg, vals):
    """least squares p, gamma in  value = -gamma eps^p  (eps in radians)."""
    xs = np.log([math.radians(e) for e in eps_deg])
    ys = np.log([-v for v in vals])
    if len(xs) < 2:
        return None, None
    p, c = np.polyfit(xs, ys, 1)
    return float(p), float(math.exp(c))


def cmd_cells(a):
    """exponent fit and ladder depth per cell, from an `eval` jsonl."""
    recs = [json.loads(l) for l in open(a.file)]
    cells = defaultdict(dict)
    for r in recs:
        if r.get('delta') is None or r.get('eps') is None:
            continue
        key = a.key(r) if callable(a.key) else (r.get('k'), 'chainfree' if
                                                r.get('pidx', -1) >= 0 else 'free')
        e = r['eps']
        if e not in cells[key] or r['delta'] > cells[key][e]['delta']:
            cells[key][e] = r
    print(f"# {'cell':22s} " + ' '.join(f"{'e=%g' % e:>13s}" for e in a.epss)
          + f"  {'p':>6s} {'gamma':>9s}  depth(by eps)")
    for key in sorted(cells, key=lambda t: (str(t))):
        row = cells[key]
        es = [e for e in a.epss if e in row and row[e]['delta'] < -1e-14]
        if len(es) < 2:
            continue
        p, g = fit_exponent(es, [row[e]['delta'] for e in es])
        dd = ' '.join(f"{e:g}:{row[e].get('cdepth')}/{row[e].get('wdepth')}" for e in es)
        print(f"  {str(key):22s} "
              + ' '.join(f"{row[e]['delta']:+13.5e}" if e in row else ' ' * 13
                         for e in a.epss)
              + f"  {p:6.3f} {g:9.4f}  {dd}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    def common(p):
        p.add_argument('--tol', type=float, default=1e-7)
        p.add_argument('--maxchains', type=int, default=200)
        p.add_argument('--kmax', type=int, default=9)
        p.add_argument('--nocyc', action='store_true')
        p.add_argument('--nch2', type=int, default=2)
        p.add_argument('--cap2', type=int, default=3000)
        p.add_argument('--dmax', type=float, default=-1e-12)
        p.add_argument('--nproc', type=int, default=2)
        p.add_argument('--seed', type=int, default=20260921)
        p.add_argument('--out', required=True)

    p = sub.add_parser('hunt')
    p.add_argument('--eps', default='1,5')
    p.add_argument('--ks', default='8,9')
    p.add_argument('--npat', type=int, default=2)
    p.add_argument('--nrand', type=int, default=1)
    p.add_argument('--jpf', type=int, default=0)
    p.add_argument('--maxfam', type=int, default=0)
    p.add_argument('--limitjobs', type=int, default=0)
    p.add_argument('--limit', type=int, default=120)
    p.add_argument('--extra', type=int, default=120)
    p.add_argument('--njit', type=int, default=12)
    p.add_argument('--free', action='store_true')
    common(p)
    p.set_defaults(fn=cmd_hunt)

    p = sub.add_parser('eval')
    p.add_argument('file', nargs='*')
    p.add_argument('--census', action='store_true')
    p.add_argument('--root', default='.')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--only', nargs='*', default=None)
    p.add_argument('--maxn', type=int, default=0)
    common(p)
    p.set_defaults(fn=cmd_eval)

    p = sub.add_parser('tophole')
    p.add_argument('--T', type=float, default=5.0)
    p.add_argument('--eps', default='1,2,5')
    p.add_argument('--nsamp', type=int, default=400)
    p.add_argument('--drop', action='store_true')
    p.add_argument('--mix', action='store_true')
    p.add_argument('--limit', type=int, default=0)
    p.add_argument('--extra', type=int, default=0)
    p.add_argument('--njit', type=int, default=8)
    p.add_argument('--rho', type=float, default=0.05)
    p.add_argument('--nproc', type=int, default=2)
    p.add_argument('--seed', type=int, default=4242)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_tophole)

    p = sub.add_parser('tophole2')
    p.add_argument('--T', type=float, default=5.0)
    p.add_argument('--eps', default='1,2,5')
    p.add_argument('--nsamp', type=int, default=200)
    p.add_argument('--drop', action='store_true')
    p.add_argument('--mix', action='store_true')
    p.add_argument('--njit', type=int, default=4)
    p.add_argument('--rho', type=float, default=0.05)
    p.add_argument('--eta', type=float, default=1e-6)
    p.add_argument('--cycles', type=int, default=2)
    p.add_argument('--nproc', type=int, default=2)
    p.add_argument('--seed', type=int, default=515)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_tophole2)

    p = sub.add_parser('cells')
    p.add_argument('file')
    p.add_argument('--epss', type=float, nargs='*', default=[0.5, 1.0, 2.0, 5.0, 10.0])
    p.set_defaults(fn=cmd_cells, key=None)

    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
