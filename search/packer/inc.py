#!/usr/bin/env python3
"""Incidence contact model for the SLP squeeze (f64, numpy-vectorised); the f64 twin of exact/find_contacts.

Constraint g >= 0 per feature incidence (same functions as exact/geom.py):
  wall:  corner a of square i vs wall L/R/B/T          g = Px, S - Px, Py, S - Py
  pair:  corner a of square t on the line of side k of square o      g = n_ok . (P_ta - c_o) - H
A pair contributes incidences only on side lines that separate it (every corner of t outside, up to `sep`), and only for
corners that project onto the side (|tau| <= H + reach).  So a side-side contact is its two segment-end incidences (one
corner of each square in the offset case), not "both corners of t outside o's line": the old rigid.pair_rows model,
which forbids feasible relative rotations and jams falsely.
Pure corner-corner touches (every incidence at a side end, no side line holding two corners) are disjunctive: the pair is
disjoint iff one of the candidate side lines separates.  They are returned separately, one alternative per line.

All gaps/gradients are first order; the LP is  min dS  s.t.  g + grad . z >= 0,  |z| <= R.
"""
import math
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix, vstack as svstack

H = 0.5
CU = np.array([(1, 1), (-1, 1), (-1, -1), (1, -1)], float)          # corner a = c + H(su u + sv v)
EK = np.array([(1, 0), (0, 1), (-1, 0), (0, -1)], float)            # side k normal in the local frame


def frames(sq):
    q = np.asarray(sq, float)
    X, Y, T = q[:, 0], q[:, 1], q[:, 2]
    C, S = np.cos(T), np.sin(T)
    # offsets of corners: (n, 4, 2);  normals: (n, 4, 2)
    ox = H * (C[:, None] * CU[None, :, 0] - S[:, None] * CU[None, :, 1])
    oy = H * (S[:, None] * CU[None, :, 0] + C[:, None] * CU[None, :, 1])
    nx = C[:, None] * EK[None, :, 0] - S[:, None] * EK[None, :, 1]
    ny = S[:, None] * EK[None, :, 0] + C[:, None] * EK[None, :, 1]
    return X, Y, ox, oy, nx, ny


def model(s, sq, delta, reach=None, sep=None, etol=None):
    """Returns (rows, disj):
      rows: (G, I, J, V) -> gaps G (m,), sparse grad entries (row I, col J, value V); col 3n = side S
      disj: list of (pair, [alt]), alt = (G, I, J, V) with rows numbered from 0 within the alternative.
    delta: incidence range (g < delta).  reach: how far past a side end a corner still counts (default delta).
    sep: a line separates if every corner of t has g > -sep (default delta).  etol: 'at the side end' (default delta)."""
    reach = delta if reach is None else reach
    sep = delta if sep is None else sep
    etol = delta if etol is None else etol
    n = len(sq)
    X, Y, ox, oy, nx, ny = frames(sq)
    PX, PY = X[:, None] + ox, Y[:, None] + oy
    G, I, J, V = [], [], [], []
    m = 0
    DS = 3 * n
    # walls (vectorised): g for each corner and wall
    for w, (g, gx, gy, gs, dth) in {
            'L': (PX, 1.0, 0.0, 0.0, -oy), 'R': (s - PX, -1.0, 0.0, 1.0, oy),
            'B': (PY, 0.0, 1.0, 0.0, ox), 'T': (s - PY, 0.0, -1.0, 1.0, -ox)}.items():
        ii, aa = np.nonzero(g < delta)
        for i, a in zip(ii.tolist(), aa.tolist()):
            G.append(g[i, a])
            ent = [(3 * i + 2, dth[i, a])]
            if gx: ent.append((3 * i, gx))
            if gy: ent.append((3 * i + 1, gy))
            if gs: ent.append((DS, gs))
            for c, v in ent:
                I.append(m); J.append(c); V.append(v)
            m += 1
    # near pairs
    d2 = (X[:, None] - X[None, :]) ** 2 + (Y[:, None] - Y[None, :]) ** 2
    pi, pj = np.nonzero(np.triu(d2 < (math.sqrt(2) + delta) ** 2, 1))
    disj = []
    for i, j in zip(pi.tolist(), pj.tolist()):
        inc = []                                   # (o, t, k, a, g, tau)
        for o, t in ((i, j), (j, i)):
            dx = PX[t][None, :] - X[o]             # (1, 4) corners of t relative to centre of o
            dy = PY[t][None, :] - Y[o]
            g = nx[o][:, None] * dx + ny[o][:, None] * dy - H          # (4 sides, 4 corners)
            tau = -ny[o][:, None] * dx + nx[o][:, None] * dy
            ok_line = g.min(axis=1) > -sep
            for k in np.nonzero(ok_line)[0].tolist():
                for a in range(4):
                    if g[k, a] < delta and abs(tau[k, a]) <= H + reach:
                        inc.append((o, t, k, a, g[k, a], tau[k, a]))
        if not inc:
            continue
        lines = {}
        for e in inc:
            lines.setdefault((e[0], e[2]), []).append(e)
        at_end = all(abs(e[5]) > H - etol for e in inc)
        if at_end and max(len(v) for v in lines.values()) < 2 and len(lines) >= 2:
            # alternative = "this side line separates": every corner of the other square outside it, including corners
            # past the end of the side (a nearly parallel neighbour's far corner can dip below the line)
            alts = []
            for (o, k) in lines:
                t = j if o == i else i
                dx, dy = PX[t] - X[o], PY[t] - Y[o]
                g = nx[o, k] * dx + ny[o, k] * dy - H
                tau = -ny[o, k] * dx + nx[o, k] * dy
                aG, aI, aJ, aV = [], [], [], []
                r = 0
                for a in range(4):
                    if g[a] < delta:
                        _row((o, t, k, a, g[a], tau[a]), X, Y, ox, oy, nx, ny, r, aG, aI, aJ, aV)
                        r += 1
                alts.append((aG, aI, aJ, aV))
            disj.append(((i, j), alts))
        else:
            for e in inc:
                if abs(e[5]) <= H + 1e-9:          # corners past a side end only matter for corner-corner pairs
                    _row(e, X, Y, ox, oy, nx, ny, m, G, I, J, V)
                    m += 1
    return (np.array(G), np.array(I, int), np.array(J, int), np.array(V)), disj



def model_fast(s, sq, delta, reach=None, sep=None, etol=None):
    """Vectorised model(): identical rows (same order, same floats) and disjunctions.  Checked against model() by
    test_inc_fast.py; model() stays as the reference."""
    reach = delta if reach is None else reach
    sep = delta if sep is None else sep
    etol = delta if etol is None else etol
    n = len(sq)
    X, Y, ox, oy, nx, ny = frames(sq)
    PX, PY = X[:, None] + ox, Y[:, None] + oy
    G, I, J, V = [], [], [], []
    m = 0
    DS = 3 * n
    for w, (g, gx, gy, gs, dth) in {
            'L': (PX, 1.0, 0.0, 0.0, -oy), 'R': (s - PX, -1.0, 0.0, 1.0, oy),
            'B': (PY, 0.0, 1.0, 0.0, ox), 'T': (s - PY, 0.0, -1.0, 1.0, -ox)}.items():
        ii, aa = np.nonzero(g < delta)
        for i, a in zip(ii.tolist(), aa.tolist()):
            G.append(g[i, a])
            ent = [(3 * i + 2, dth[i, a])]
            if gx: ent.append((3 * i, gx))
            if gy: ent.append((3 * i + 1, gy))
            if gs: ent.append((DS, gs))
            for c, v in ent:
                I.append(m); J.append(c); V.append(v)
            m += 1
    G = np.array(G, float); I = np.array(I, int); J = np.array(J, int); V = np.array(V, float)
    d2 = (X[:, None] - X[None, :]) ** 2 + (Y[:, None] - Y[None, :]) ** 2
    pi, pj = np.nonzero(np.triu(d2 < (math.sqrt(2) + delta) ** 2, 1))
    disj = []
    if len(pi) == 0:
        return (G, I, J, V), disj
    O = np.stack([pi, pj], 1)                       # (M, 2): o for dir 0 / 1
    T = np.stack([pj, pi], 1)
    dx = PX[T][:, :, None, :] - X[O][:, :, None, None]           # (M, 2, 1, 4)
    dy = PY[T][:, :, None, :] - Y[O][:, :, None, None]
    nxo, nyo = nx[O][:, :, :, None], ny[O][:, :, :, None]        # (M, 2, 4, 1)
    g = nxo * dx + nyo * dy - H                                  # (M, 2, 4 sides, 4 corners)
    tau = -nyo * dx + nxo * dy
    ok_line = g.min(axis=3) > -sep
    inc = ok_line[..., None] & (g < delta) & (np.abs(tau) <= H + reach)
    has = inc.any(axis=(1, 2, 3))
    cnt = inc.sum(axis=3)                                        # incidences per (dir, side line)
    at_end = ~(inc & (np.abs(tau) <= H - etol)).any(axis=(1, 2, 3))
    isdj = has & at_end & (cnt.max(axis=(1, 2)) < 2) & ((cnt > 0).sum(axis=(1, 2)) >= 2)
    # regular incidences: non-disjunctive pairs, corners on the side proper
    keep = inc & ~isdj[:, None, None, None] & (np.abs(tau) <= H + 1e-9)
    q, d, k, a = np.nonzero(keep)                                # C order = old loop order (pair, dir, side, corner)
    if len(q):
        o, t = O[q, d], T[q, d]
        n_x, n_y = nx[o, k], ny[o, k]
        ddx, ddy = PX[t, a] - X[o], PY[t, a] - Y[o]
        rows = m + np.arange(len(q))
        Jn = np.stack([3 * o, 3 * o + 1, 3 * t, 3 * t + 1, 3 * o + 2, 3 * t + 2], 1)
        Vn = np.stack([-n_x, -n_y, n_x, n_y, -n_y * ddx + n_x * ddy, -n_x * oy[t, a] + n_y * ox[t, a]], 1)
        G = np.concatenate([G, g[q, d, k, a]]); I = np.concatenate([I, np.repeat(rows, 6)])
        J = np.concatenate([J, Jn.ravel()]); V = np.concatenate([V, Vn.ravel()])
    for qq in np.nonzero(isdj)[0].tolist():
        i, j = int(pi[qq]), int(pj[qq])
        lines = []
        for dd in (0, 1):
            for kk in range(4):
                if cnt[qq, dd, kk] > 0:
                    lines.append((int(O[qq, dd]), kk))
        # old code: dict insertion order = order of first incidence = (dir, side) order
        alts = []
        for (o, kk) in lines:
            t = j if o == i else i
            ddx, ddy = PX[t] - X[o], PY[t] - Y[o]
            gg = nx[o, kk] * ddx + ny[o, kk] * ddy - H
            tt = -ny[o, kk] * ddx + nx[o, kk] * ddy
            aG, aI, aJ, aV = [], [], [], []
            r = 0
            for a_ in range(4):
                if gg[a_] < delta:
                    _row((o, t, kk, a_, gg[a_], tt[a_]), X, Y, ox, oy, nx, ny, r, aG, aI, aJ, aV)
                    r += 1
            alts.append((aG, aI, aJ, aV))
        disj.append(((i, j), alts))
    return (G, I, J, V), disj

def _row(e, X, Y, ox, oy, nx, ny, r, G, I, J, V):
    o, t, k, a, g, tau = e
    n_x, n_y = nx[o, k], ny[o, k]
    dx, dy = X[t] + ox[t, a] - X[o], Y[t] + oy[t, a] - Y[o]
    G.append(g)
    for c, v in ((3 * o, -n_x), (3 * o + 1, -n_y), (3 * t, n_x), (3 * t + 1, n_y),
                 (3 * o + 2, -n_y * dx + n_x * dy), (3 * t + 2, -n_x * oy[t, a] + n_y * ox[t, a])):
        I.append(r); J.append(c); V.append(v)


def _mat(rows, nz):
    G, I, J, V = rows
    m = len(G)
    return csr_matrix((V, (I, J)), shape=(m, nz)), np.asarray(G, float)


def _lp(A, g, nz, R, exact, fixed, time_limit=60.0, clip=True):
    """min dS s.t. g + A z >= 0, |z| <= R  (in units of R).  clip=False keeps negative gaps (overlaps must be opened)."""
    b = (np.zeros_like(g) if exact else (np.maximum(g, 0.0) if clip else g)) / R
    c = np.zeros(nz); c[-1] = 1.0
    bnd = np.array([(-1.0, 1.0)] * nz)
    for v in fixed:
        bnd[v] = (0.0, 0.0)
    res = linprog(c, A_ub=-A, b_ub=b, bounds=bnd, method='highs',
                  options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10,
                           'time_limit': time_limit})
    if res.status != 0:
        return None, None
    return res.fun * R, res.x * R


def _alt_slack(alt, z):
    G, I, J, V = alt
    if not G:
        return 0.0
    v = np.asarray(G, float).copy()
    np.add.at(v, np.asarray(I), np.asarray(V) * z[np.asarray(J)])
    return v.min()


def shrink(s, sq, delta, R, fixed=(), exact=False, switch=True, rounds=3):
    """One SLP subproblem.  Disjunctive pairs: relax (drop them) -> z; round each to the alternative z violates least;
    restricted LP; repeat.  Returns (ds, z, n_rows, n_disj)."""
    n = len(sq)
    nz = 3 * n + 1
    rows, disj = model_fast(s, sq, delta)
    A, g = _mat(rows, nz)
    if not disj:
        ds, z = _lp(A, g, nz, R, exact, fixed)
        return ds, z, len(g), 0

    def restricted(choice):
        # one COO assembly (was one csr per alternative, then vstack: ~1/4 of a descent)
        Gs, Is, Js, Vs = [np.asarray(rows[0], float)], [np.asarray(rows[1], int)], [np.asarray(rows[2], int)], [np.asarray(rows[3], float)]
        off = len(rows[0])
        for (p, alts), q in zip(disj, choice):
            aG, aI, aJ, aV = alts[q]
            if not aG:
                continue
            Gs.append(np.asarray(aG, float)); Is.append(np.asarray(aI, int) + off); Js.append(np.asarray(aJ, int)); Vs.append(np.asarray(aV, float))
            off += len(aG)
        M = csr_matrix((np.concatenate(Vs), (np.concatenate(Is), np.concatenate(Js))), shape=(off, nz))
        return _lp(M, np.concatenate(Gs), nz, R, exact, fixed)

    # default branch: the line with the largest current minimum gap (as the old fixed-axis choice)
    choice = [max(range(len(alts)), key=lambda q: min(alts[q][0])) for _, alts in disj]
    best = restricted(choice)
    if switch:
        _, z = _lp(A, g, nz, R, exact, fixed)             # relaxation: corner-corner pairs dropped
        for _ in range(rounds):
            if z is None:
                break
            ch = [max(range(len(alts)), key=lambda q: _alt_slack(alts[q], z)) for _, alts in disj]
            if ch == choice:
                break
            choice = ch
            r = restricted(choice)
            if r[0] is None:
                break
            if best[0] is None or r[0] < best[0] - 1e-18:
                best = r
            z = r[1]
    return best[0], best[1], len(g), len(disj)


def shrink_milp(s, sq, delta=1e-7, R=1e-3, fixed=(), exact=True, time_limit=60.0):
    """Exhaustive version: one binary per candidate line of each corner-corner pair (big-M).  For the final jam test."""
    n = len(sq)
    nz = 3 * n + 1
    rows, disj = model_fast(s, sq, delta)
    A, g = _mat(rows, nz)
    nb = sum(len(a) for _, a in disj)
    b0 = (np.zeros_like(g) if exact else np.maximum(g, 0.0)) / R
    Afull = [np.hstack([A.toarray(), np.zeros((A.shape[0], nb))])]
    lo = [-b0]
    ones = []
    col = nz
    for _, alts in disj:
        first = col
        for alt in alts:
            Aq, gq = _mat(alt, nz)
            gq = (np.zeros_like(gq) if exact else gq) / R
            M = np.abs(gq) + np.abs(Aq.toarray()).sum(axis=1) + 1.0
            blk = np.hstack([Aq.toarray(), np.zeros((Aq.shape[0], nb))])
            blk[:, col] = -M                       # gq + Aq z' + M(1 - y) >= 0
            Afull.append(blk); lo.append(-gq - M)
            col += 1
        r = np.zeros(nz + nb); r[first:col] = 1.0
        ones.append(r)
    Amat = np.vstack(Afull + ([np.array(ones)] if ones else []))
    L = np.concatenate(lo + ([np.ones(len(ones))] if ones else []))
    lb = np.r_[-np.ones(nz), np.zeros(nb)]
    ub = np.r_[np.ones(nz), np.ones(nb)]
    for v in fixed:
        lb[v] = ub[v] = 0.0
    c = np.zeros(nz + nb); c[nz - 1] = 1.0
    res = milp(c, constraints=LinearConstraint(Amat, L, np.inf), integrality=np.r_[np.zeros(nz), np.ones(nb)],
               bounds=Bounds(lb, ub), options={'time_limit': time_limit, 'mip_rel_gap': 1e-4, 'mip_abs_gap': 1e-12,
                                                          'mip_feasibility_tolerance': 1e-9})
    if res.x is None:
        return None, None, len(disj)
    return res.fun * R, res.x[:nz] * R, len(disj)


def correct(s, sq, R, fixed=()):
    """Second-order correction: after a step, restore the linearised contacts at the new point with the side fixed
    (overlaps from the step's curvature are O(R^2); this removes them to O(R^4) instead of paying for them by scaling the
    whole packing).  Corner-corner pairs take their currently best line.  Returns the corrected configuration or None."""
    n = len(sq)
    nz = 3 * n + 1
    rows, disj = model_fast(s, sq, 3 * R)
    A, g = _mat(rows, nz)
    blocks, gs = [A], [g]
    for _, alts in disj:
        q = max(range(len(alts)), key=lambda q: min(alts[q][0]))
        Aq, gq = _mat(alts[q], nz)
        blocks.append(Aq); gs.append(gq)
    g = np.concatenate(gs)
    if g.min() >= 0:
        return sq
    # the side is fixed, so the LP has no objective: box the correction at ~10x the overlap (an arbitrary vertex of a
    # box of size R would move everything by R and create new O(R^2) overlaps)
    Af = svstack(blocks).tocsr()
    r = 10.0 * -g.min()
    z = None
    while z is None and r <= 10 * R:
        _, z = _lp(Af, g, nz, r, False, tuple(fixed) + (nz - 1,), clip=False)
        r *= 10.0
    if z is None:
        return None
    return [(float(x + z[3 * i]), float(y + z[3 * i + 1]), float(t + z[3 * i + 2])) for i, (x, y, t) in enumerate(sq)]
