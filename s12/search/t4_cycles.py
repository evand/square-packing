#!/usr/bin/env python3
"""t4_cycles: which CYCLE certificates exist in the LP duals at T = 4, n = 12.

Task tasks/t4-cycles/README.md.  Nothing in search/ is modified; t3_chain, t3_chain_shape,
bandcut_k, s6skel, s6local are imported.

Semantics and row system: search/t3_chain.py (= BANDCUT_K.md sec 1.1 = S6_SKELETON.md sec 3.1),
the wall-carrying LP, every row  a . c - delta >= b.  A Farkas certificate is w >= 0 with
sum w_r a_r = 0 and sum w_r = 1, proving delta <= -sum w_r b_r.

Objects measured here
  * the FULL dual at a recorded/fresh optimum and its support graph (cyclomatic number mu,
    wall rows per wall, turns) -- t3_chain_shape.support_graph;
  * the best H-restricted bound -- t3_chain.best_h (mode 'H');
  * the best PURE-CYCLE bound: a restricted dual supported on the pair rows of one simple cycle
    of the TIGHT separation graph plus the wall rows, nothing else;
  * the closed form of notes/t3-chain.md sec 5.3 evaluated on an explicit cycle, together with
    the Farkas residual that says whether its hypotheses actually hold.

A restricted-dual LP bound is a Farkas bound and is therefore >= delta* by construction: it can
never be < 0 at a margin-0 point.  The closed FORM can be, if its hypotheses are not checked;
that is what the filter (cmd_filter) tests.
"""
import argparse
import itertools
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s6skel                                                     # noqa: E402
import t3_chain as TC                                             # noqa: E402
import t3_chain_shape as TS                                       # noqa: E402
import bandcut_k as BK                                            # noqa: E402

TOL = 1e-9


# ===================================================================== loading recorded optima
def _z_from_sq(sq):
    n = len(sq)
    z = np.zeros(1 + 3 * n)
    for i, (x, y, th) in enumerate(sq):
        z[1 + 3 * i], z[2 + 3 * i], z[3 + 3 * i] = x, y, th
    return z


def load_recorded(root='.'):
    """-> list of dicts(src, cell, z, delta).  Every recorded T = 4, n = 12 optimum."""
    out = []
    for fn, tag in (('runs/bandcut_k_T4a.jsonl', 'T4a'),
                    ('runs/bandcut_k_T4free.jsonl', 'T4free'),
                    ('runs/bandcut_k_T4free_pilot.jsonl', 'T4pilot'),
                    ('runs/bandcut_k_T4_hunt.jsonl', 'T4hunt')):
        p = os.path.join(root, fn)
        if not os.path.exists(p):
            continue
        for line in open(p):
            r = json.loads(line)
            if r.get('value') is None or r.get('z') is None:
                continue
            z = np.array(r['z'], dtype=float)
            if (len(z) - 1) // 3 != 12:
                continue
            free = 'chainfree' if r['pidx'] >= 0 else 'free'
            out.append(dict(src=tag, cell=f"k={r['k']} eps={r['eps']:g} {free}",
                            k=r['k'], eps=r['eps'], free=free,
                            z=z, delta=float(r['value'])))
    p = os.path.join(root, 'runs/bandcut_scan_pT_T4.json')
    if os.path.exists(p):
        d = json.load(open(p))
        out.append(dict(src='bandstack', cell='(4,1) band stack', k=8, eps=None, free='free',
                        z=_z_from_sq(d['sq']), delta=float(d['delta'])))
    for p3 in ('runs/bandcut_scan_bestbar_T4_p3.json', 'runs/bandcut_scan_bestbar_T4_p4.json'):
        p = os.path.join(root, p3)
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        sq = d.get('sq') or d.get('squares')
        if sq is None:
            continue
        out.append(dict(src='bestbar', cell=os.path.basename(p3)[:-5].split('_')[-1] + ' bar',
                        k=None, eps=None, free='free',
                        z=_z_from_sq(sq), delta=float(d.get('delta', 0.0))))
    return out


# ===================================================================== the tight separation graph
def tight_indices(A, b, z, lvl, tol=1e-7):
    """indices of rows whose slack at z is <= tol.

    Complementary slackness: the optimal dual of the LP at z is supported on tight rows only, so
    cycles that can carry a certificate of value delta* live in the graph of tight PAIR rows."""
    c = np.concatenate([z[1::3], z[2::3]])
    slack = A @ c - lvl - b
    return [q for q in range(len(b)) if abs(slack[q]) <= tol], slack


def tight_rows(z, n, T, lvl, tol=1e-7):
    A, b, tags = TC.all_rows(z, n, T, lvl, tol=tol)
    tight, slack = tight_indices(A, b, z, lvl, tol=tol)
    return A, b, tags, tight, slack


def tight_pair_graph(tags, tight):
    """-> dict (i,j) -> [row indices], over tight PAIR rows only."""
    g = defaultdict(list)
    for q in tight:
        tg = tags[q]
        if tg[0] == 'pair':
            g[(min(tg[1], tg[2]), max(tg[1], tg[2]))].append(q)
    return dict(g)


def simple_cycles(edges, kmax=12, cap=200000):
    """every simple cycle (as a vertex tuple, canonical rotation/reflection) of the undirected
    graph given by `edges` (iterable of (i,j)), of length 3..kmax.  Stops at `cap` cycles."""
    adj = defaultdict(set)
    for (i, j) in edges:
        adj[i].add(j); adj[j].add(i)
    seen, out = set(), []

    def canon(cyc):
        k = len(cyc)
        best = None
        for r in range(k):
            for d in (1, -1):
                t = tuple(cyc[(r + d * s) % k] for s in range(k))
                if best is None or t < best:
                    best = t
        return best

    def walk(start, path, pset):
        last = path[-1]
        for nb in adj[last]:
            if len(out) >= cap:
                return
            if nb == start and len(path) >= 3:
                c = canon(path)
                if c not in seen:
                    seen.add(c); out.append(c)
            elif nb not in pset and nb > start and len(path) < kmax:
                walk(start, path + [nb], pset | {nb})
    for s in sorted(adj):
        if len(out) >= cap:
            break
        walk(s, [s], {s})
    return out


# ===================================================================== restricted dual bounds
def wall_rows(tags):
    return [q for q, tg in enumerate(tags) if tg[0] in ('lo', 'hi')]


def cycle_support(tags, cyc, pairrows, walls=True, tight_only=None):
    """support = every pair row on a cycle edge (+ all wall rows).  `pairrows` maps (i,j)->rows."""
    keep = set()
    k = len(cyc)
    for s in range(k):
        i, j = cyc[s], cyc[(s + 1) % k]
        keep |= set(pairrows.get((min(i, j), max(i, j)), []))
    if walls:
        keep |= set(wall_rows(tags))
    return keep


def best_cycle_bound(A, b, tags, pairrows, kmax=9, maxcyc=4000, want=False, seed=7):
    """min over simple cycles of the tight pair graph of the restricted dual on
    (that cycle's pair rows) + (all wall rows).  -> (bound, cycle, w, ncycles).

    When there are more than `maxcyc` cycles a uniform random sample of that many is taken
    (fixed seed), so the number reported is a LOWER bound on what a cycle can do -- the honest
    direction for a certificate search, which can only fail to find one."""
    cycs = simple_cycles(list(pairrows), kmax=kmax)
    use = cycs
    if len(cycs) > maxcyc:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(cycs), size=maxcyc, replace=False)
        use = [cycs[i] for i in idx]
    best = (math.inf, None, None)
    for cyc in use:
        keep = cycle_support(tags, cyc, pairrows)
        v, w = TC.dual_bound(A, b, keep, want_w=True)
        if v < best[0]:
            best = (v, cyc, w)
    return best[0], best[1], best[2], len(cycs)


# ===================================================================== the sec 5.3 closed form
def _cycle_choices(z, cyc, rows_by_edge, tags, b, cap=200000):
    """-> list per cycle edge of (unit normal as the link cyc[s] -> cyc[s+1], m, row index)."""
    TH = z[3::3]
    k = len(cyc)
    cand = []
    for s in range(k):
        i, j = cyc[s], cyc[(s + 1) % k]
        lst = []
        for q in rows_by_edge.get((min(i, j), max(i, j)), []):
            _, aa, bb, o, kind, sg = tags[q]
            nrm = np.array([math.cos(TH[o]), math.sin(TH[o])]) if kind == 0 else \
                np.array([-math.sin(TH[o]), math.cos(TH[o])])
            # the row reads  sg * nrm . (c_bb - c_aa) >= m + delta
            e = sg * nrm if (aa, bb) == (i, j) else -sg * nrm
            lst.append((e, float(b[q]), q))
        if not lst:
            return None
        cand.append(lst)
    if float(np.prod([len(c) for c in cand])) > cap:
        return None
    return cand


def cycle_closed_form(z, n, T, cyc, rows_by_edge, tags, b, A=None, axtol=1e-7):
    """evaluate notes/t3-chain.md sec 5.3 on every choice of one row per cycle edge.

    For a cyclic sequence of links with normals n_0..n_{k-1} the Farkas residual at vertex
    cyc[s] is  n_{s-1} - n_s , so ONE wall row of weight lambda_s = |n_s - n_{s-1}| cancels it
    exactly when that difference is axis-aligned (the lemma's hypothesis).  The value returned is
    the exact  -sum w b / sum w  of that weight vector; it is a genuine Farkas bound iff
    `axis_ok`.  When no choice conforms, the number the naive formula WOULD give is still
    reported (`naive`), because that is what the filter has to catch.

    -> dict(best=<best conforming, or None>, naive=<smallest value over all choices>, ...)"""
    P = 0.5 * (np.abs(np.cos(z[3::3])) + np.abs(np.sin(z[3::3])))
    k = len(cyc)
    cand = _cycle_choices(z, cyc, rows_by_edge, tags, b)
    if cand is None:
        return None
    wallrow = {}
    for q, tg in enumerate(tags):
        if tg[0] in ('lo', 'hi'):
            wallrow[(tg[0], tg[1], tg[2])] = q
    best, naive = None, None
    for choice in itertools.product(*cand):
        ns = [c[0] for c in choice]
        ms = [c[1] for c in choice]
        lam, ok, whichwall = [], True, []
        for s in range(k):
            d = ns[s] - ns[s - 1]
            nd = float(np.hypot(*d))
            if nd < 1e-9:
                continue                    # straight vertex: no wall row needed
            ax = 0 if abs(d[0]) >= abs(d[1]) else 1
            if abs(d[1 - ax]) > axtol:
                ok = False                  # not axis-aligned: one wall row cannot cancel it
            lam.append(nd)
            whichwall.append(('lo' if d[ax] > 0 else 'hi', 'xy'[ax], int(cyc[s])))
        W = k + sum(lam)
        wb = sum(ms) + sum(l * (P[v] if lh == 'lo' else P[v] - T)
                           for l, (lh, _, v) in zip(lam, whichwall))
        val = float(-wb / W)
        lo = sum(l for l, (lh, _, _) in zip(lam, whichwall) if lh == 'lo')
        hi = sum(l for l, (lh, _, _) in zip(lam, whichwall) if lh == 'hi')
        rec = dict(value=val, q=len(lam), k=k, lam=[float(x) for x in lam],
                   axis_ok=bool(ok), lo=float(lo), hi=float(hi),
                   walls=whichwall, rows=[int(c[2]) for c in choice])
        if A is not None:
            w = np.zeros(len(b))
            for c in choice:
                w[c[2]] += 1.0
            miss = False
            for l, key in zip(lam, whichwall):
                q = wallrow.get(key)
                if q is None:
                    miss = True
                    break
                w[q] += l
            if not miss:
                rec['resid'], _, _ = farkas_residual(A, b, w)
            else:
                rec['resid'] = None
                rec['axis_ok'] = False      # the needed wall row is not even in the system
        if naive is None or val < naive['value']:
            naive = rec
        if rec['axis_ok'] and (best is None or val < best['value']):
            best = rec
    return dict(best=best, naive=naive)


def cycle_closed_form2(z, n, T, cyc, rows_by_edge, tags, b, A=None):
    """the sec 5.3 cycle bound GENERALISED: a turn whose normal jump is not axis-aligned is
    cancelled by TWO wall rows (one per axis) at the same square instead of one.

    The turn weight is then the L1 norm |d_x| + |d_y| of the jump, not its Euclidean norm, and
    the hypothesis 'every turn is axis-aligned' disappears -- so this version applies to every
    cycle, and it is still an exact Farkas certificate (`resid` is printed).

        delta  <=  [ sum_j ( lam_jx (T - u_j)/2 + lam_jy (T - u_j)/2 ) - sum_i m_i ]
                   / ( k + sum_j (lam_jx + lam_jy) )        when the lo/hi split is even,
    which it always is: summing the link coefficients round a closed cycle telescopes to 0, so
    the wall coefficients sum to 0 axis by axis.  [proved]"""
    P = 0.5 * (np.abs(np.cos(z[3::3])) + np.abs(np.sin(z[3::3])))
    k = len(cyc)
    cand = _cycle_choices(z, cyc, rows_by_edge, tags, b)
    if cand is None:
        return None
    wallrow = {}
    for q, tg in enumerate(tags):
        if tg[0] in ('lo', 'hi'):
            wallrow[(tg[0], tg[1], tg[2])] = q
    best = None
    for choice in itertools.product(*cand):
        ns = [c[0] for c in choice]
        ms = [c[1] for c in choice]
        walls, lam, q = [], [], 0
        for s in range(k):
            d = ns[s] - ns[s - 1]
            if float(np.hypot(*d)) < 1e-9:
                continue
            q += 1
            for ax in (0, 1):
                if abs(d[ax]) < 1e-12:
                    continue
                walls.append(('lo' if d[ax] > 0 else 'hi', 'xy'[ax], int(cyc[s])))
                lam.append(abs(float(d[ax])))
        W = k + sum(lam)
        wb = sum(ms) + sum(l * (P[v] if lh == 'lo' else P[v] - T)
                           for l, (lh, _, v) in zip(lam, walls))
        val = float(-wb / W)
        rec = dict(value=val, q=q, nwall=len(lam), k=k, lam=[float(x) for x in lam],
                   walls=walls, rows=[int(c[2]) for c in choice])
        if A is not None:
            w = np.zeros(len(b))
            for c in choice:
                w[c[2]] += 1.0
            okw = True
            for l, key in zip(lam, walls):
                qq = wallrow.get(key)
                if qq is None:
                    okw = False
                    break
                w[qq] += l
            rec['resid'] = farkas_residual(A, b, w)[0] if okw else None
        if best is None or val < best['value']:
            best = rec
    return best


def farkas_residual(A, b, w):
    """(max |sum w a|, sum w, bound) -- the honest check that w is a Farkas certificate."""
    r = A.T @ w
    return float(np.max(np.abs(r))), float(np.sum(w)), float(-(b @ w))


# ===================================================================== classification
def classify(g):
    """the t3-chain.md sec 3.1 classes.  `g` = t3_chain_shape.support_graph output.

    chain      mu = 0, the pair rows form one path, wall rows on ONE axis only (the bare
               wall-to-wall chain of T: T-1 links + 2 walls)
    tree       mu = 0 otherwise (an H, or a chain with transverse legs)
    chain+rung mu = 1 but fewer than one wall row per wall
    cycle      mu = 1 and all four walls carry a row (the sec 3.3 object)
    other      mu >= 2
    """
    mu, nw = g['mu'], len(g['walls'])
    axes = {(k[1] if isinstance(k, tuple) else k.split('-')[1]) for k in g['walls']}
    if mu == 0:
        return 'chain' if (nw <= 2 and len(axes) <= 1
                           and g['nedge'] == g['nvert'] - 1) else 'tree'
    if mu == 1:
        return 'cycle' if nw >= 4 else 'chain+rung'
    return 'other'


# ===================================================================== commands
def analyse(rec, T=4.0, kmax=9, do_h=True, do_cycle=True, maxchains=120):
    z, n = rec['z'], (len(rec['z']) - 1) // 3
    lvl = rec['delta']
    A, b, tags, tight, slack = tight_rows(z, n, T, lvl)
    fv, w = TC.dual_bound(A, b, range(len(b)), want_w=True)
    TH = np.array(z[3::3])
    supp = [(float(w[q]), tags[q], TC.row_axis(tags[q], TH)) for q in range(len(b))
            if w[q] > 1e-6]
    supp.sort(key=lambda t: -t[0])
    g = TS.support_graph(supp)
    cls = classify(g)
    g['walls'] = {f'{k[0]}-{k[1]}': float(v) for k, v in g['walls'].items()}
    out = dict(src=rec['src'], cell=rec['cell'], k=rec.get('k'), eps=rec.get('eps'),
               delta=float(lvl), fulldual=float(fv), shape=g, cls=cls,
               tilts=[float(math.degrees(s6skel.norm_tilt(t))) for t in TH],
               z=[float(q) for q in z],
               maxtilt=float(max(abs(math.degrees(s6skel.norm_tilt(t))) for t in TH)),
               supp=[[s[0], list(map(TC._js, s[1])), s[2]] for s in supp],
               ntight=len(tight))
    pg = tight_pair_graph(tags, tight)
    out['ntight_edges'] = len(pg)
    if do_h:
        r = TC.best_h(z, n, T, lvl, modes=('H',), maxchains=maxchains, rows=(A, b, tags))
        out['H'] = float(r['H']['bound']); out['H_path'] = r['H']['path']
        out['H_axis'] = r['H']['axis']; out['nchains'] = r['nchains']
    if do_cycle:
        cv, cyc, cw, ncyc = best_cycle_bound(A, b, tags, pg, kmax=kmax)
        out['cycle'] = float(cv); out['cycle_path'] = list(cyc) if cyc else None
        out['ncycles'] = ncyc
        if cyc is not None:
            out['cf'] = cycle_closed_form(z, n, T, cyc, pg, tags, b, A=A)
    return out


# ================================================== the eps^3 cell: which object carries it?
def ladder(z, n, T, lvl, tol=1e-7, maxchains=400, kmax=9):
    """the bound as the support is widened, one geometric fact at a time.

    C     the wall-to-wall chain of T alone (links + its two end walls)
    CW    that chain + every wall row
    H     the H-lemma support: chain + every wall + every TRANSVERSE-type pair row
    HP    H + the main-direction rows BETWEEN chain squares (staircase closure)
    H1    H + ONE further main-direction pair row (the t3-chain.md sec 2.5 repair), best choice
    CYC   the best pure-cycle support (one simple cycle's pair rows + every wall row)
    F     everything the configuration satisfies  (= delta* at a global optimum)
    """
    A, b, tags = TC.all_rows(z, n, T, lvl, tol=tol)
    TH = np.array(z[3::3])
    chs = TC.chains_of(z, n, lvl, int(round(T)), tol=tol)
    out = {'nchains': len(chs)}
    for mode in ('C', 'CW', 'H', 'HP'):
        best = (math.inf, None)
        for (ax, path) in chs[:maxchains]:
            keep = TC.h_support(tags, TH, ax, path, mode)
            v, _ = TC.dual_bound(A, b, keep)
            if v < best[0]:
                best = (v, (ax, list(path)))
        out[mode] = best[0]
        out[mode + '_path'] = best[1]
    # H + one extra main-direction pair row
    best = (math.inf, None, None)
    for (ax, path) in chs[:maxchains]:
        base = TC.h_support(tags, TH, ax, path, 'H')
        extra = [q for q, tg in enumerate(tags)
                 if tg[0] == 'pair' and q not in base and TC.row_axis(tg, TH) == ax]
        for q in extra:
            v, _ = TC.dual_bound(A, b, base | {q})
            if v < best[0]:
                best = (v, (ax, list(path)), tags[q])
    out['H1'] = best[0]
    out['H1_path'] = best[1]
    out['H1_extra'] = None if best[2] is None else list(map(TC._js, best[2]))
    tight, _ = tight_indices(A, b, z, lvl, tol=tol)
    pg = tight_pair_graph(tags, tight)
    cv, cyc, cw, ncyc = best_cycle_bound(A, b, tags, pg, kmax=kmax)
    out['CYC'] = cv
    out['CYC_path'] = list(cyc) if cyc else None
    out['ncycles'] = ncyc
    fv, w = TC.dual_bound(A, b, range(len(b)), want_w=True)
    out['F'] = fv
    supp = [(float(w[q]), tags[q], TC.row_axis(tags[q], TH)) for q in range(len(b))
            if w[q] > 1e-9]
    supp.sort(key=lambda t: -t[0])
    out['supp'] = [[s[0], list(map(TC._js, s[1])), s[2]] for s in supp]
    out['shape'] = {}
    for wt in (1e-6, 1e-9):
        g = TS.support_graph(supp, wtol=wt)
        gg = dict(g, walls={f'{k[0]}-{k[1]}': float(v) for k, v in g['walls'].items()})
        out['shape'][f'{wt:g}'] = dict(mu=g['mu'], E=g['nedge'], V=g['nvert'],
                                       walls=gg['walls'], cls=classify(gg))
    # does the HEAVY part of the dual alone still certify?  (re-solve on w > 1e-6 rows only)
    heavy = [q for q in range(len(b)) if w[q] > 1e-6]
    out['heavy'] = TC.dual_bound(A, b, heavy)[0]
    out['nheavy'] = len(heavy)
    hs = [(float(w[q]), tags[q], TC.row_axis(tags[q], TH)) for q in heavy]
    gh = TS.support_graph(hs, wtol=0.0)
    out['heavy_mu'] = gh['mu']
    return out


def cmd_cell89(a):
    recs = load_recorded(a.root)
    sel = [r for r in recs if r.get('k') in (a.kk or (8, 9)) and r['free'] == 'chainfree']
    cells = {}
    for r in sel:
        key = (r['k'], r['eps'])
        if key not in cells or r['delta'] > cells[key]['delta']:
            cells[key] = r
    for key in sorted(cells):
        r = cells[key]
        z, n = r['z'], 12
        print(f"\n### k={key[0]} eps={key[1]:g} chain-free   delta* = {r['delta']:+.6e}")
        tl = [math.degrees(s6skel.norm_tilt(t)) for t in z[3::3]]
        print('    tilts (deg): ' + ' '.join(f'{q}:{tl[q]:+.3f}' for q in range(n)))
        for tol in a.tols:
            L = ladder(z, n, a.T, r['delta'], tol=tol, kmax=a.kmax)
            s6, s9 = L['shape']['1e-06'], L['shape']['1e-09']
            print(f"  tol={tol:g}  nchains={L['nchains']} ncycles={L['ncycles']}  "
                  f"dual support @w>1e-6: mu={s6['mu']} E={s6['E']} V={s6['V']} "
                  f"cls={s6['cls']} | @w>1e-9: mu={s9['mu']} E={s9['E']} V={s9['V']} "
                  f"cls={s9['cls']}")
            print(f"     C  {L['C']:+.6e}   CW {L['CW']:+.6e}   H  {L['H']:+.6e}   "
                  f"HP {L['HP']:+.6e}")
            print(f"     H1 {L['H1']:+.6e} (extra row {L['H1_extra']})   "
                  f"CYC {L['CYC']:+.6e} {L['CYC_path']}   F {L['F']:+.6e}")
            print(f"     heavy-only (w>1e-6, {L['nheavy']} rows, mu={L['heavy_mu']}): "
                  f"{L['heavy']:+.6e}")
            if tol == a.tols[0]:
                for w, tg, ax in L['supp'][:18]:
                    sq = tg[2:] if tg[0] in ('lo', 'hi') else tg[1:3]
                    print(f"       {w:11.8f} {str(tg):<32s} {ax}  tilts "
                          + ' '.join(f'{tl[q]:+.3f}' for q in sq))


# ============================================ the dichotomy failures and the ladder that fixes
def cmd_dich(a):
    """every sampled optimum at which neither H <= 0 nor cycle <= 0 nor delta* >= 0, re-tested
    with the LADDER: H + one extra main-direction pair row (t3-chain.md sec 2.5's repair)."""
    recs = []
    for f in a.file:
        recs += [json.loads(l) for l in open(f)]
    recs = [r for r in recs if 'err' not in r]
    byz = {}
    for r in load_recorded(a.root):
        byz[(r['src'], r['cell'], round(r['delta'], 15))] = r['z']
    bad = [r for r in recs if r.get('H', 9) > 1e-9 and r.get('cycle', 9) > 1e-9
           and r['delta'] < -1e-12]
    print(f'# {len(bad)} / {len(recs)} sampled optima with H > 0 AND cycle > 0 AND delta* < 0.')
    print(f'# {"cell":24s} {"delta*":>12s} {"H":>12s} {"CYC":>12s} {"H1":>12s} '
          f'{"F":>12s} {"mu":>3s} {"H1==F?":>7s}')
    nfix = 0
    for r in bad:
        z = r.get('z')
        if z is None:
            z = byz.get((r['src'], r['cell'], round(r['delta'], 15)))
        if z is None:
            print(f"  {r['cell']:24s} {r['delta']:+12.4e}  (z not recoverable)")
            continue
        z = np.array(z, dtype=float)
        L = ladder(z, 12, a.T, r['delta'], kmax=a.kmax)
        ok = 'yes' if L['H1'] <= 1e-9 else 'NO'
        nfix += L['H1'] <= 1e-9
        print(f"  {r['cell']:24s} {r['delta']:+12.4e} {r['H']:+12.4e} {r['cycle']:+12.4e} "
              f"{L['H1']:+12.4e} {L['F']:+12.4e} {r['shape']['mu']:3d} {ok:>7s}")
    print(f'\n# H + one extra main-direction link certifies delta <= 0 at {nfix} / {len(bad)}.')


# ===================================================================== the filter
def filter_point(z, n, T, lvl, kmax=9, maxcyc=4000, maxchains=400, tol=1e-7, seed=7):
    """every chain and cycle bound available at one configuration.

    CW/H/HP   restricted duals on a wall-to-wall chain of T  (Farkas, so always >= delta*)
    CYC       restricted dual on one simple cycle + walls    (Farkas, so always >= delta*)
    CF        the sec 5.3 closed form, minimised over cycles and row choices, restricted to
              choices whose hypotheses HOLD (axis-aligned turns, exact cancellation)
    CFnaive   the same formula with the hypotheses NOT checked -- not a bound, and the number
              the filter has to catch."""
    A, b, tags = TC.all_rows(z, n, T, lvl, tol=tol)
    TH = np.array(z[3::3])
    out = {}
    chs = TC.chains_of(z, n, lvl, int(round(T)), tol=tol)
    out['nchains'] = len(chs)
    for mode in ('CW', 'H', 'HP'):
        best = math.inf
        for (ax, path) in chs[:maxchains]:
            v, _ = TC.dual_bound(A, b, TC.h_support(tags, TH, ax, path, mode))
            best = min(best, v)
        out[mode] = best
    tight, _ = tight_indices(A, b, z, lvl, tol=tol)
    pg = tight_pair_graph(tags, tight)
    cycs = simple_cycles(list(pg), kmax=kmax)
    out['ncycles'] = len(cycs)
    use = cycs
    if len(cycs) > maxcyc:
        rng = np.random.default_rng(seed)
        use = [cycs[i] for i in rng.choice(len(cycs), maxcyc, replace=False)]
    out['ncyc_eval'] = len(use)
    bl, bc, bn = math.inf, math.inf, math.inf
    arg = {}
    for c in use:
        v, _ = TC.dual_bound(A, b, cycle_support(tags, c, pg))
        if v < bl:
            bl, arg['CYC'] = v, list(c)
        cf = cycle_closed_form(z, n, T, c, pg, tags, b, A=A)
        if cf is None:
            continue
        if cf['best'] is not None and cf['best']['value'] < bc:
            bc, arg['CF'] = cf['best']['value'], (list(c), cf['best']['q'],
                                                  cf['best'].get('resid'))
        if cf['naive'] is not None and cf['naive']['value'] < bn:
            bn, arg['CFnaive'] = cf['naive']['value'], (list(c), cf['naive']['q'],
                                                        cf['naive']['axis_ok'])
    out['CYC'], out['CF'], out['CFnaive'] = bl, bc, bn
    out['arg'] = arg
    out['F'] = TC.dual_bound(A, b, range(len(b)))[0]
    return out


def cmd_filter(a):
    recs = load_recorded(a.root)
    pts = [r for r in recs if r['delta'] >= -a.zero]
    if a.only:
        pts = [r for r in pts if any(s in r['cell'] for s in a.only)]
    seen, uniq = set(), []
    for r in pts:
        key = tuple(np.round(r['z'], 9))
        if key in seen:
            continue
        seen.add(key); uniq.append(r)
    print(f'# {len(uniq)} distinct T = 4 optima with delta* >= -{a.zero:g}.')
    print('# every bound below must be >= 0.  CW/H/HP/CYC are Farkas bounds (restricted duals),')
    print('# so they cannot be < 0; CF is the sec 5.3 closed form WITH its hypotheses checked;')
    print('# CFnaive is the same formula with the hypotheses NOT checked.')
    print(f'# {"point":26s} {"delta*":>11s} {"F":>11s} {"CW":>11s} {"H":>11s} {"HP":>11s} '
          f'{"CYC":>11s} {"CF":>11s} {"CFnaive":>11s} {"ncyc":>7s}')
    worst = {}
    for r in uniq:
        f = filter_point(r['z'], 12, a.T, r['delta'], kmax=a.kmax, maxcyc=a.maxcyc)
        print(f"  {r['cell']:26s} {r['delta']:+11.3e} {f['F']:+11.3e} {f['CW']:+11.3e} "
              f"{f['H']:+11.3e} {f['HP']:+11.3e} {f['CYC']:+11.3e} {f['CF']:+11.3e} "
              f"{f['CFnaive']:+11.3e} {f['ncycles']:7d}")
        for kk in ('CW', 'H', 'HP', 'CYC', 'CF', 'CFnaive'):
            if kk not in worst or f[kk] < worst[kk][0]:
                worst[kk] = (f[kk], r['cell'], f['arg'].get(kk))
    print('\n# tightest (most negative) value of each bound over all margin-0 points:')
    for kk in ('CW', 'H', 'HP', 'CYC', 'CF', 'CFnaive'):
        v, cell, ag = worst[kk]
        flag = '   <-- BELOW ZERO' if v < -1e-8 else ''
        print(f'  {kk:8s} {v:+.6e}  at {cell}  {ag}{flag}')


# ===================================================================== fresh samples
T4_FAMS = (['unif1', 'unif5', 'unif10', 'unif20', 'unif28.0725', 'unif30', 'unif45']
           + [f'hole{k}_{e}' for k in (0, 1, 2, 3) for e in (1, 5, 20)]      # the far field
           + [f'hole{k}_{e}' for k in (4, 6, 8, 9) for e in (1, 5, 10)]      # the hole
           + ['near1', 'near5', 'near20', 'coh1', 'coh5', 'coh20']
           + ['Z4', 'Z8', 'Z9', 'Z12']
           + ['generic', 'bigrand5', 'bigrand10', 'bigrand20', 'bigrand30', 'bigrand40'])

_FCFG = None


def _finit(cfg):
    global _FCFG
    _FCFG = cfg


def _fjob(arg):
    fam, seed = arg
    c = _FCFG
    rng = np.random.default_rng(seed)
    th, label = TC.sample_theta(fam, rng, c['n'])
    v, z, _ = TC.delta_star(th, c['T'], rng, limit=c['limit'], extra=c['extra'])
    rec = dict(src='fresh', cell=fam, k=None, eps=None, z=z, delta=float(v))
    o = analyse(rec, T=c['T'], kmax=c['kmax'])
    o['fam'] = fam
    o['seed'] = seed
    o['label'] = label
    o['z'] = [float(q) for q in z]
    # how many squares are near-axis at 1 deg / 5 deg (the k of BANDCUT_K)
    tl = o['tilts']
    o['k1'] = sum(1 for t in tl if abs(t) <= 1.0 + 1e-9)
    o['k5'] = sum(1 for t in tl if abs(t) <= 5.0 + 1e-9)
    o['maxtilt'] = max(abs(t) for t in tl)
    return o


def cmd_fresh(a):
    import multiprocessing as mp
    import time
    fams = a.fams.split(',') if a.fams else T4_FAMS
    jobs = [(f, a.seed0 + 1009 * i + j) for i, f in enumerate(fams) for j in range(a.reps)]
    cfg = dict(T=a.T, n=a.n, limit=a.limit, extra=a.extra, kmax=a.kmax)
    t0 = time.time()
    print(f'# {len(jobs)} fresh samples, n={a.n} T={a.T}', file=sys.stderr, flush=True)
    with open(a.out, 'w') as fh:
        with mp.Pool(a.nproc, initializer=_finit, initargs=(cfg,)) as pool:
            for i, o in enumerate(pool.imap_unordered(_fjob, jobs, chunksize=1)):
                fh.write(json.dumps(o) + '\n'); fh.flush()
                if i % 10 == 0:
                    print(f'  {i}/{len(jobs)}  {time.time()-t0:.0f}s',
                          file=sys.stderr, flush=True)
    print(f'done {time.time()-t0:.0f}s', file=sys.stderr)


def cmd_report(a):
    recs = [json.loads(l) for l in open(f) for f in ()] if False else []
    for f in a.file:
        recs += [json.loads(l) for l in open(f)]
    recs = [r for r in recs if 'err' not in r]
    cells = defaultdict(list)
    for r in recs:
        cells[r['cell']].append(r)
    print(f'# {len(recs)} T=4 optima.  mu = cyclomatic number of the dual support graph.')
    print(f'# {"cell":26s} {"N":>3s} {"max delta*":>12s} {"mu=0":>5s} {"mu=1":>5s} {"mu>=2":>5s}'
          f'  {"chain":>5s} {"tree":>5s} {"c+rung":>6s} {"cycle":>5s} {"other":>5s}'
          f'  {"H<=0":>7s} {"cyc<=0":>7s}')
    tot = defaultdict(int)
    for cell in sorted(cells):
        rs = cells[cell]
        mu = [r['shape']['mu'] for r in rs]
        cl = [r['cls'] for r in rs]
        hok = sum(1 for r in rs if r.get('H', 9) <= 1e-9)
        cok = sum(1 for r in rs if r.get('cycle', 9) <= 1e-9)
        for r in rs:
            tot['n'] += 1
            tot['mu%d' % min(r['shape']['mu'], 2)] += 1
            tot[r['cls']] += 1
            if r.get('H', 9) <= 1e-9:
                tot['H'] += 1
            if r.get('cycle', 9) <= 1e-9:
                tot['C'] += 1
            if r.get('H', 9) <= 1e-9 or r.get('cycle', 9) <= 1e-9 or r['delta'] >= -1e-12:
                tot['dich'] += 1
        print(f'  {cell:26s} {len(rs):3d} {max(r["delta"] for r in rs):+12.4e} '
              f'{mu.count(0):5d} {mu.count(1):5d} {sum(1 for x in mu if x>=2):5d}  '
              f'{cl.count("chain"):5d} {cl.count("tree"):5d} {cl.count("chain+rung"):6d} '
              f'{cl.count("cycle"):5d} {cl.count("other"):5d}  '
              f'{hok:3d}/{len(rs):<3d} {cok:3d}/{len(rs):<3d}')
    print(f"\n# totals {tot['n']}: mu=0 {tot['mu0']}, mu=1 {tot['mu1']}, mu>=2 {tot['mu2']}")
    print(f"# classes: chain {tot['chain']}, tree {tot['tree']}, chain+rung {tot['chain+rung']},"
          f" cycle {tot['cycle']}, other {tot['other']}")
    print(f"# H<=0 {tot['H']}/{tot['n']}, cycle<=0 {tot['C']}/{tot['n']}, "
          f"dichotomy (H<=0 or cycle<=0 or delta*>=0) {tot['dich']}/{tot['n']}")

    print('\n# by max tilt of the configuration')
    bins = defaultdict(lambda: [0, 0, 0])
    for r in recs:
        b = int(max(abs(t) for t in r['tilts']) // 5) * 5
        bins[b][min(r['shape']['mu'], 2)] += 1
    print(f'# {"maxtilt":>9s} {"mu=0":>6s} {"mu=1":>6s} {"mu>=2":>6s}')
    for b in sorted(bins):
        print(f'  {b:3d}-{b+5:<5d}' + ''.join(f'{x:6d}' for x in bins[b]))

    print('\n# the cyclic duals: up to 15 examples')
    ex = [r for r in recs if r['shape']['mu'] >= 1]
    ex.sort(key=lambda r: -r['delta'])
    for r in ex[:15]:
        print(f"  {r['cell']:24s} d*={r['delta']:+.4e} mu={r['shape']['mu']} "
              f"E={r['shape']['nedge']} V={r['shape']['nvert']} walls={len(r['shape']['walls'])}"
              f" cls={r['cls']} H={r.get('H',float('nan')):+.4e} "
              f"cyc={r.get('cycle',float('nan')):+.4e}")
        print('      ' + '  '.join(f"{w:.3f}{tg[0]}{tuple(tg[1:3])}" for w, tg, ax in r['supp'][:9]))

    print('\n# filter: any bound < 0 at a delta* >= -1e-12 point?')
    bad = [r for r in recs if r['delta'] >= -1e-12
           and (r.get('H', 0) < -1e-9 or r.get('cycle', 0) < -1e-9)]
    print(f'#   {len(bad)} violations among '
          f'{sum(1 for r in recs if r["delta"] >= -1e-12)} margin-0 optima')
    for r in bad[:20]:
        print(f"    {r['cell']:24s} d*={r['delta']:+.3e} H={r.get('H'):+.4e} "
              f"cyc={r.get('cycle'):+.4e}")


_CFG = None


def _init(cfg):
    global _CFG
    _CFG = cfg


def _job(r):
    try:
        return analyse(r, T=_CFG['T'], kmax=_CFG['kmax'])
    except Exception as e:                                           # noqa: BLE001
        return dict(src=r['src'], cell=r['cell'], delta=r['delta'], err=repr(e))


def cmd_census(a):
    import multiprocessing as mp
    recs = load_recorded(a.root)
    if a.limit:
        recs = recs[:a.limit]
    print(f'# {len(recs)} recorded T=4 optima', file=sys.stderr, flush=True)
    cfg = dict(T=a.T, kmax=a.kmax)
    with open(a.out, 'w') as fh:
        if a.nproc > 1:
            with mp.Pool(a.nproc, initializer=_init, initargs=(cfg,)) as pool:
                for i, o in enumerate(pool.imap_unordered(_job, recs, chunksize=1)):
                    fh.write(json.dumps(o) + '\n'); fh.flush()
                    if i % 25 == 0:
                        print(f'  {i}/{len(recs)}', file=sys.stderr, flush=True)
        else:
            _init(cfg)
            for i, r in enumerate(recs):
                fh.write(json.dumps(_job(r)) + '\n'); fh.flush()
                if i % 25 == 0:
                    print(f'  {i}/{len(recs)}', file=sys.stderr, flush=True)
    print('done', file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('census')
    p.add_argument('--root', default='.')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--kmax', type=int, default=9)
    p.add_argument('--limit', type=int, default=0)
    p.add_argument('--nproc', type=int, default=12)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_census)
    p = sub.add_parser('report')
    p.add_argument('file', nargs='+')
    p.set_defaults(fn=cmd_report)
    p = sub.add_parser('dich')
    p.add_argument('file', nargs='+')
    p.add_argument('--root', default='.')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--kmax', type=int, default=9)
    p.set_defaults(fn=cmd_dich)
    p = sub.add_parser('filter')
    p.add_argument('--root', default='.')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--kmax', type=int, default=9)
    p.add_argument('--maxcyc', type=int, default=3000)
    p.add_argument('--zero', type=float, default=1e-12)
    p.add_argument('--only', nargs='*', default=None)
    p.set_defaults(fn=cmd_filter)
    p = sub.add_parser('cell89')
    p.add_argument('--root', default='.')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--kk', type=int, nargs='*', default=None)
    p.add_argument('--kmax', type=int, default=9)
    p.add_argument('--tols', type=float, nargs='*', default=[1e-7, 1e-10])
    p.set_defaults(fn=cmd_cell89)
    p = sub.add_parser('fresh')
    p.add_argument('--T', type=float, default=4.0)
    p.add_argument('--n', type=int, default=12)
    p.add_argument('--fams', default='')
    p.add_argument('--reps', type=int, default=6)
    p.add_argument('--limit', type=int, default=200)
    p.add_argument('--extra', type=int, default=200)
    p.add_argument('--kmax', type=int, default=9)
    p.add_argument('--seed0', type=int, default=20260921)
    p.add_argument('--nproc', type=int, default=12)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_fresh)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
