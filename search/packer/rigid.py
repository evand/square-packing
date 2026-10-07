#!/usr/bin/env python3
"""First-order (LP) shrink / jam / rigidity analysis of a packing (Donev-Torquato-Stillinger-Connelly style).

Contacts: pairs whose separating gap < delta, walls with clearance < delta.  Each is linearized in
z = (dx_i, dy_i, dtheta_i)_i and ds (side change).  Pair gap uses the currently best separating axis (a
restriction of the true disjunctive constraint, so "jammed" means jammed for that axis choice).  Kinks of
|cos|, |sin| at 0 (axis-aligned squares) are split into one constraint per sign, which is exact for them.

  shrink LP:   min ds  s.t. linearized contacts, |z| <= R      -> ds* < 0: shrink direction; ds* = 0: jammed
  rigidity:    with ds = 0, max/min each coordinate of each square -> squares with zero range are rigid

Usage: rigid.py file.txt [--delta 1e-7] [--rigidity]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')   # BLAS threads spin (4 s CPU per import, 32 threads); parallelism is across jobs
import math, sys, time
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

H = 0.5


def load(path):
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1])
    return s, [(float(a), float(b), math.radians(float(c))) for a, b, c, *_ in L[1:n + 1]]


def _sgns(x, tiny=1e-9):
    return (-1.0, 1.0) if abs(x) < tiny else ((1.0,) if x > 0 else (-1.0,))


def _axis_rows(dx, dy, owner, cp, sp, proj, ca, sa, ktol):
    sgp = 1.0 if proj >= 0 else -1.0
    rows = []
    for s1 in _sgns(ca, ktol):
        for s2 in _sgns(sa, ktol):
            dr = H * (-s1 * sa + s2 * ca)            # d r_other / d alpha
            dproj_dphi = sgp * (-dx * sp + dy * cp)    # d |proj| / d phi
            g_owner = dproj_dphi - dr                  # alpha = phi - theta_other
            g_other = dr
            gdx, gdy = sgp * cp, sgp * sp              # w.r.t. dx = xj - xi
            gi_t, gj_t = (g_owner, g_other) if owner == 0 else (g_other, g_owner)
            rows.append((-gdx, -gdy, gi_t, gdx, gdy, gj_t))
    return rows


def pair_axes(qi, qj):
    """The 4 candidate separating axes (side normals of either square): list of (gap, owner, cp, sp, proj, ca, sa).
    The pair is disjoint iff some gap > 0 (SAT); the true gap is the max."""
    dx, dy = qj[0] - qi[0], qj[1] - qi[1]
    out = []
    for owner in (0, 1):
        to, tot = (qi[2], qj[2]) if owner == 0 else (qj[2], qi[2])
        for k in (0, 1):
            phi = to + k * math.pi / 2
            cp, sp = math.cos(phi), math.sin(phi)
            proj = dx * cp + dy * sp
            al = phi - tot
            ca, sa = math.cos(al), math.sin(al)
            out.append((abs(proj) - H - H * (abs(ca) + abs(sa)), owner, cp, sp, proj, ca, sa))
    return out


def pair_rows(qi, qj, ktol=1e-9):
    """Best separating axis gap and list of gradient variants w.r.t. (dxi,dyi,dti,dxj,dyj,dtj)."""
    ax = pair_axes(qi, qj)
    gap, owner, cp, sp, proj, ca, sa = max(ax, key=lambda a: a[0])
    return gap, _axis_rows(qj[0] - qi[0], qj[1] - qi[1], owner, cp, sp, proj, ca, sa, ktol)


def pair_alternatives(qi, qj, reach, ktol=1e-9):
    """Axes whose gap is within `reach` of 0 (could separate after a step of size ~reach): [(gap, rows), ...],
    best axis first.  More than one entry = a disjunctive (corner-corner type) contact."""
    dx, dy = qj[0] - qi[0], qj[1] - qi[1]
    ax = sorted(pair_axes(qi, qj), key=lambda a: -a[0])
    return [(a[0], _axis_rows(dx, dy, *a[1:], ktol)) for a in ax if a[0] > -reach]


def contacts(s, sq, delta, ktol=1e-9):
    """List of (gap, {var: coef}) linear constraints gap + coef.z >= 0.  var 3i.. for squares, 'ds' = 3n."""
    n = len(sq)
    DS = 3 * n
    cons = []
    for i, (x, y, t) in enumerate(sq):
        ct, st = math.cos(t), math.sin(t)
        w = H * (abs(ct) + abs(st))
        for s1 in _sgns(ct, ktol):
            for s2 in _sgns(st, ktol):
                dw = H * (-s1 * st + s2 * ct)
                for gap, cx, cy, cs in ((x - w, 1, 0, 0), (s - x - w, -1, 0, 1), (y - w, 0, 1, 0), (s - y - w, 0, -1, 1)):
                    if gap < delta:
                        c = {3 * i + 2: -dw}
                        if cx: c[3 * i] = cx
                        if cy: c[3 * i + 1] = cy
                        if cs: c[DS] = 1.0
                        cons.append((gap, c))
    for i in range(n):
        for j in range(i + 1, n):
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 > (math.sqrt(2) + delta) ** 2:
                continue
            gap, rows = pair_rows(sq[i], sq[j], ktol)
            if gap < delta:
                for r in rows:
                    c = {3 * i: r[0], 3 * i + 1: r[1], 3 * i + 2: r[2], 3 * j: r[3], 3 * j + 1: r[4], 3 * j + 2: r[5]}
                    cons.append((gap, c))
    return cons


def disjunctive(s, sq, delta, reach, ktol=1e-9, parallel=True):
    """Pairs in contact (best gap < delta) with >= 2 distinct axes within `reach` of separating.
    Returns {(i, j): [(gap, [coef dicts]), ...]} (best axis first)."""
    out = {}
    n = len(sq)
    for i in range(n):
        for j in range(i + 1, n):
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 > (math.sqrt(2) + delta) ** 2:
                continue
            alts = pair_alternatives(sq[i], sq[j], reach, ktol)
            if len(alts) < 2 or alts[0][0] >= delta:
                continue
            # parallel=False: only pairs with candidate axes in different directions (corner-corner type, as in
            # exact/find_contacts).  Parallel sides (one line, seen from either square) are disjunctive too, for rotations.
            dirs = [(a[2], a[3]) for a in pair_axes(sq[i], sq[j]) if a[0] > -reach]
            if not parallel and all(abs(d[0] * dirs[0][1] - d[1] * dirs[0][0]) < 1e-6 for d in dirs):
                continue
            seen, uniq = set(), []
            for gap, rows in alts:
                key = (round(gap, 12), tuple(sorted(tuple(round(v, 12) for v in r) for r in rows)))
                if key in seen:
                    continue
                seen.add(key)
                uniq.append((gap, [{3 * i: r[0], 3 * i + 1: r[1], 3 * i + 2: r[2], 3 * j: r[3], 3 * j + 1: r[4],
                                    3 * j + 2: r[5]} for r in rows]))
            if len(uniq) >= 2:
                out[(i, j)] = uniq
    return out


def shrink_milp(s, sq, delta=1e-7, R=1e-3, ktol=None, fixed=(), exact=False, time_limit=60.0, tie=1e-6, only=None, parallel=True):
    """shrink_lp with every disjunctive pair (corner-corner type) as an OR over its candidate separating axes
    (big-M, one binary per axis).  Returns (ds, z, n_disj); n_disj = 0 means no disjunctions (caller can skip)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = len(sq)
    kt = 3 * R if ktol is None else ktol
    disj = disjunctive(s, sq, min(delta, tie), reach=tie, ktol=kt, parallel=parallel)
    if only is not None:
        disj = {p: a for p, a in disj.items() if p in only}
    if not disj:
        return None, None, 0
    cons = [c for c in contacts(s, sq, delta, kt)]
    # drop the fixed-axis rows of disjunctive pairs (contacts() put the best axis in); keep everything else
    cons = [(g, c) for g, c in cons if pair_of(c) not in disj]
    nz = 3 * n + 1
    nb = sum(len(a) for a in disj.values())
    rows, lo, hi = [], [], []
    for gap, c in cons:                                   # gap/R + c.z' >= 0
        r = np.zeros(nz + nb)
        for k, v in c.items():
            r[k] = v
        rows.append(r); lo.append(-(0.0 if exact else max(gap, 0.0)) / R); hi.append(np.inf)
    b = nz
    for (i, j), alts in disj.items():
        first = b
        for a, (gap, cl) in enumerate(alts):
            g = (0.0 if exact else max(gap, 0.0)) if a == 0 else gap
            for c in cl:                                  # g/R + c.z' + M(1 - y) >= 0
                M = abs(g) / R + sum(abs(v) for v in c.values()) + 1.0
                r = np.zeros(nz + nb)
                for k, v in c.items():
                    r[k] = v
                r[b] = -M
                rows.append(r); lo.append(-g / R - M); hi.append(np.inf)
            b += 1
        r = np.zeros(nz + nb); r[first:b] = 1.0          # at least one axis separates
        rows.append(r); lo.append(1.0); hi.append(np.inf)
    lb = np.r_[-np.ones(nz), np.zeros(nb)]
    ub = np.r_[np.ones(nz), np.ones(nb)]
    for v in fixed:
        lb[v] = ub[v] = 0.0
    cvec = np.zeros(nz + nb); cvec[nz - 1] = 1.0
    integ = np.r_[np.zeros(nz), np.ones(nb)]
    res = milp(cvec, constraints=LinearConstraint(np.array(rows), lo, hi), integrality=integ, bounds=Bounds(lb, ub),
               options={'time_limit': time_limit, 'mip_rel_gap': 1e-4, 'mip_abs_gap': 1e-12})
    if res.x is None:
        return None, None, len(disj)
    return res.fun * R, res.x[:nz] * R, len(disj)


def build(n, cons, exact=False):
    A = lil_matrix((len(cons), 3 * n + 1))
    b = np.zeros(len(cons))
    for r, (gap, c) in enumerate(cons):
        for k, v in c.items():
            A[r, k] = -v          # -coef.z <= gap
        b[r] = 0.0 if exact else max(gap, 0.0)
    return A.tocsr(), b


def pair_of(c):
    """(i, j) for a pair constraint dict (6 keys), None for a wall constraint."""
    if len(c) != 6:
        return None
    k = sorted(c)
    return k[0] // 3, k[3] // 3


def shrink_lp(s, sq, delta=1e-7, R=1e-3, exact=False, ktol=None, fixed=(), duals=None):
    """duals: if a dict is passed, it is filled with {(i, j): total |multiplier|} of the pair rows."""
    n = len(sq)
    cons = contacts(s, sq, delta, 3 * R if ktol is None else ktol)
    A, b = build(n, cons, exact)
    c = np.zeros(3 * n + 1); c[-1] = 1.0
    # solve in units of R (z = R z'), so the solver's absolute feasibility tolerance scales with the trust radius
    bnd = [(-1.0, 1.0)] * (3 * n + 1)
    for v in fixed:
        bnd[v] = (0.0, 0.0)
    res = linprog(c, A_ub=A, b_ub=b / R, bounds=bnd, method='highs',
                  options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10,
                           'time_limit': 60.0})
    if res.status != 0:
        return None, None, len(cons)
    if duals is not None:
        for (g, c), m in zip(cons, res.ineqlin.marginals):
            p = pair_of(c)
            if p is not None and abs(m) > 0:
                duals[p] = duals.get(p, 0.0) + abs(m)
    return res.fun * R, res.x * R, len(cons)


def _lp(n, cons, R, exact, fixed, time_limit=60.0):
    A, b = build(n, cons, exact)
    c = np.zeros(3 * n + 1); c[-1] = 1.0
    bnd = [(-1.0, 1.0)] * (3 * n + 1)
    for v in fixed:
        bnd[v] = (0.0, 0.0)
    res = linprog(c, A_ub=A, b_ub=b / R, bounds=bnd, method='highs',
                  options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10,
                           'time_limit': time_limit})
    if res.status != 0:
        return None, None
    return res.fun * R, res.x * R


def shrink_switch(s, sq, delta=1e-7, R=1e-3, ktol=None, fixed=(), exact=False, tie=1e-6, rounds=3):
    """LP-only branch search at corner-corner touches (cheap stand-in for shrink_milp).
    Relax: drop the rows of every disjunctive pair (any axis allowed) -> z.  Round: give each such pair the candidate axis
    that z violates least, solve the restricted LP (a valid restriction: its steps are feasible to first order).
    Repeat from the new z.  Returns (ds, z, n_disj) of the best restricted LP, or (None, None, n) if none solved."""
    n = len(sq)
    kt = 3 * R if ktol is None else ktol
    disj = disjunctive(s, sq, min(delta, tie), reach=tie, ktol=kt)
    if not disj:
        return None, None, 0
    base = [(g, c) for g, c in contacts(s, sq, delta, kt) if pair_of(c) not in disj]
    ds, z = _lp(n, base, R, exact, fixed)
    best = (None, None)
    for _ in range(rounds):
        if z is None:
            break
        chosen = []
        for p, alts in disj.items():
            def slack(a):
                gap, cl = alts[a]
                g = (0.0 if exact else max(gap, 0.0)) if a == 0 else gap
                return min(g + sum(v * z[k] for k, v in c.items()) for c in cl)
            a = max(range(len(alts)), key=slack)
            gap, cl = alts[a]
            g = (0.0 if exact else max(gap, 0.0)) if a == 0 else gap
            chosen += [(g, c) for c in cl]
        ds2, z2 = _lp(n, base + chosen, R, exact, fixed)
        if ds2 is None:
            break
        if best[0] is None or ds2 < best[0] - 1e-18:
            best = (ds2, z2)
            z = z2
        else:
            break
    return best[0], best[1], len(disj)


def rigidity(s, sq, delta=1e-7, R=1e-3, tol=1e-7):
    """Per-coordinate motion ranges with ds = 0.  Returns array (n, 3): max |motion| of x, y, theta per square
    (0 = that coordinate is pinned to first order)."""
    n = len(sq)
    cons = contacts(s, sq, delta)
    A, b = build(n, cons, exact=True)
    bounds = [(-R, R)] * (3 * n) + [(0, 0)]
    rng = np.zeros((n, 3))
    for i in range(n):
        for k in range(3):
            for sign in (1.0, -1.0):
                c = np.zeros(3 * n + 1); c[3 * i + k] = -sign
                res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
                if res.status == 0:
                    rng[i, k] = max(rng[i, k], -res.fun)
                if rng[i, k] > tol:
                    break
    return rng, len(cons)


def classify(rng, tol=1e-7):
    pinned = (rng <= tol).sum(axis=1)
    return int((pinned == 3).sum()), int(((pinned > 0) & (pinned < 3)).sum()), int((pinned == 0).sum())


if __name__ == '__main__':
    path = sys.argv[1]
    delta = float(sys.argv[sys.argv.index('--delta') + 1]) if '--delta' in sys.argv else 1e-7
    s, sq = load(path)
    t0 = time.time()
    ds, z, nc = shrink_lp(s, sq, delta)
    ds0, _, _ = shrink_lp(s, sq, delta, exact=True)
    print(f'{path}: n={len(sq)} s={s:.12f} contacts={nc} shrink ds*={ds:.3e} (gaps as slack), {ds0:.3e} (exact contacts: <0 = not first-order jammed) ({time.time() - t0:.2f}s)')
    if '--rigidity' in sys.argv:
        t0 = time.time()
        rng, _ = rigidity(s, sq, delta)
        r, p, f = classify(rng)
        print(f'  rigid {r}  partly pinned {p}  free {f}  (of {len(sq)}; {time.time() - t0:.1f}s)')
