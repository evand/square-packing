#!/usr/bin/env python3
"""chains.py -- what holds a zero-margin configuration shut: first-order rigidity cores, by type.

Why (2026-09-19).  `search/RANK8.md` §3 linearised leaf A about the 4x4 tiling and read the
infeasibility certificate as wall-to-wall chains; §1 showed the zero set is a plateau reaching 32
degrees of tilt, so the statement a proof needs is about EVERY interior-disjoint configuration, not
about the tiling.  This file does the first-order analysis at an arbitrary configuration, with NO
pattern constraints (patterns are bookkeeping; the geometric fact should not need them):

  unknowns   v = (dx_k, dy_k, dphi_k) per square, and delta
  walls      every active wall/vertex row of square k:                d(slack) . v >= delta
  pairs      every touching pair (gap < tol):   d(gap; v) >= delta, where
             gap = max over the 4 edge normals of  min over the sign branches of |cos|+|sin|,
             so d(gap; v) = max over ACTIVE normals of min over ACTIVE branches: a disjunction
             (binary per active normal) of conjunctions.
  maximise   delta, |v| <= 1.

  delta* > 0  => z + eps v is a genuine closed packing for small eps: a COUNTEREXAMPLE to
                 s(12) = 4 (every active constraint strictly improves, the rest have room).
                 Checked by an explicit line search on the exact gap; never expected.
  delta* = 0  => first-order rigid.  The CORE is an irreducible set of contacts (pair contacts
                 and wall contacts as units) that is still rigid: drop any one and a direction
                 opens.  Found by a deletion filter that tries to delete contacts of tilted
                 squares first, so an axis-parallel core is found whenever one exists.

The taxonomy is over cores: which squares, how tilted, and -- when the core is axis-parallel --
its decomposition into x-chains and y-chains and whether each is anchored on both walls.  A core
that is one straight wall-to-wall row of four axis-parallel squares is killed by the chord lemma
(`notes/chord-lemma.md`); anything else is what a uniform lemma has to be about.

Floating point throughout; a measurement.

  analyse FILE...   FILE is a rank8.py configuration file or a census JSONL
"""
import argparse
import json
import math
import os
import sys
from collections import Counter

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rank8                                                         # noqa: E402

T = 4.0
SGN = ((1, 1), (1, -1), (-1, 1), (-1, -1))
WALLN = ('L', 'R', 'B', 'T')


def norm_tilt(th):
    return (th + math.pi / 4) % (math.pi / 2) - math.pi / 4


# ===================================================================== active structure
def pair_structure(z, i, j, tol):
    """-> (gap, [ (axis_vec, [row gradient over (xi, yi, thi, xj, yj, thj)] ) for active axes ])"""
    X, Y, TH = z[1::3], z[2::3], z[3::3]
    dcx, dcy = X[j] - X[i], Y[j] - Y[i]
    D = TH[j] - TH[i]
    cD, sD = math.cos(D), math.sin(D)
    axes = []
    for o in (0, 1):
        th = TH[i] if o == 0 else TH[j]
        c, s = math.cos(th), math.sin(th)
        for (d0, d1, p0, p1) in ((c, s, -s, c), (-s, c, -c, -s)):
            pr = d0 * dcx + d1 * dcy
            sg = 1.0 if pr >= 0 else -1.0
            combos = []
            for (s1, s2) in SGN:
                val = sg * pr - 0.5 - 0.5 * (s1 * cD + s2 * sD)
                g = np.zeros(6)
                g[0], g[1] = -sg * d0, -sg * d1
                g[3], g[4] = sg * d0, sg * d1
                g[2 if o == 0 else 5] += sg * (p0 * dcx + p1 * dcy)
                w = -0.5 * (-s1 * sD + s2 * cD)
                g[5] += w
                g[2] -= w
                combos.append((val, g))
            aval = min(v for v, _g in combos)
            axes.append((aval, (sg * d0, sg * d1), [g for v, g in combos if v <= aval + tol]))
    gap = max(a[0] for a in axes)
    act = [(d, rows) for (v, d, rows) in axes if v >= gap - tol]
    return gap, act


def wall_structure(z, k, tol):
    """-> [(side, [gradients over (xk, yk, thk)])] for the active walls of square k"""
    x, y, th = z[1 + 3 * k], z[2 + 3 * k], z[3 + 3 * k]
    c, s = math.cos(th), math.sin(th)
    out = []
    for side, (base, gx, gy) in enumerate(((x, 1, 0), (T - x, -1, 0), (y, 0, 1), (T - y, 0, -1))):
        rows = []
        for (s1, s2) in SGN:
            val = base - 0.5 * (s1 * c + s2 * s)
            if val <= tol:
                rows.append(np.array([gx, gy, -0.5 * (-s1 * s + s2 * c)], dtype=float))
        if rows:
            out.append((side, rows))
    return out


def structure(z, n, tol):
    pairs = {}
    for i in range(n):
        for j in range(i + 1, n):
            gap, act = pair_structure(z, i, j, tol)
            if gap <= tol:
                pairs[(i, j)] = act
    walls = {}
    for k in range(n):
        for (side, rows) in wall_structure(z, k, tol):
            walls[(k, side)] = rows
    return pairs, walls


# ===================================================================== the first-order MILP
def first_order(n, pairs, walls, use=None):
    """max delta over the contacts in `use` (default: all).  -> (delta, v)"""
    from scipy.optimize import milp, LinearConstraint, Bounds
    nv = 1 + 3 * n
    BIG = 20.0
    rowsA, lo = [], []
    nb = 0
    binrows = []
    for key, rows in walls.items():
        if use is not None and ('w',) + key not in use:
            continue
        k = key[0]
        for g in rows:
            r = np.zeros(nv)
            r[0] = -1.0
            r[1 + 3 * k: 4 + 3 * k] = g
            rowsA.append((r, None))
            lo.append(0.0)
    for key, act in pairs.items():
        if use is not None and ('p',) + key not in use:
            continue
        i, j = key
        single = len(act) == 1
        bs = []
        for (_d, rows) in act:
            b = None if single else nb
            if not single:
                bs.append(nb)
                nb += 1
            for g in rows:
                r = np.zeros(nv)
                r[0] = -1.0
                r[1 + 3 * i: 4 + 3 * i] += g[:3]
                r[1 + 3 * j: 4 + 3 * j] += g[3:]
                rowsA.append((r, b))
                lo.append(0.0 if single else -BIG)
        if bs:
            binrows.append(bs)
    if not rowsA:
        return 5.0, None
    A = np.zeros((len(rowsA) + len(binrows), nv + nb))
    lob = np.zeros(len(rowsA) + len(binrows))
    for r, (row, b) in enumerate(rowsA):
        A[r, :nv] = row
        if b is not None:
            A[r, nv + b] = -BIG
        lob[r] = lo[r]
    for q, bs in enumerate(binrows):
        A[len(rowsA) + q, [nv + b for b in bs]] = 1.0
        lob[len(rowsA) + q] = 1.0
    cost = np.zeros(nv + nb)
    cost[0] = -1.0
    lb = np.concatenate([[-1.0], -np.ones(nv - 1), np.zeros(nb)])
    ub = np.concatenate([[5.0], np.ones(nv - 1), np.ones(nb)])
    integ = np.concatenate([np.zeros(nv), np.ones(nb)])
    res = milp(c=cost, constraints=LinearConstraint(A, lob, np.inf), integrality=integ,
               bounds=Bounds(lb, ub))
    if not res.success:
        return float('nan'), None
    return -res.fun, res.x[1:nv]


def true_value(z, n):
    return min(rank8.true_objective(z, n)[0], rank8.admissibility_slack(z, n))


def line_search(z, v, n):
    best, arg = true_value(z, n), 0.0
    for e in np.geomspace(1e-9, 0.3, 60):
        w = z.copy()
        w[1:] += e * v
        val = true_value(w, n)
        if val > best:
            best, arg = val, e
    return best, arg


# ===================================================================== the core
def core(z, n, pairs, walls, rng, ntry, zero):
    tilt = np.abs(norm_tilt(z[3::3]))
    allc = [('p',) + k for k in pairs] + [('w',) + k for k in walls]

    def weight(c):
        return max(tilt[c[1]], tilt[c[2]]) if c[0] == 'p' else tilt[c[1]]

    best = None
    for t in range(ntry):
        order = list(allc)
        rng.shuffle(order)
        # tilted contacts are tried for deletion first; among equals, random
        order.sort(key=lambda c: -round(float(weight(c)), 6))
        keep = set(allc)
        for c in order:
            trial = keep - {c}
            d, _v = first_order(n, pairs, walls, trial)
            if d <= zero:
                keep = trial
        score = (round(float(max([weight(c) for c in keep], default=0.0)), 6), len(keep))
        if best is None or score < best[0]:
            best = (score, keep)
    return best[1]


def signature(z, n, keep, pairs, names):
    tilt = np.abs(norm_tilt(z[3::3]))
    sq = sorted({c[1] for c in keep} | {c[2] for c in keep if c[0] == 'p'})
    mt = max(tilt[k] for k in sq) if sq else 0.0
    npair = sum(1 for c in keep if c[0] == 'p')
    nwall = sum(1 for c in keep if c[0] == 'w')
    if mt > 1e-6:
        return (f'TILTED squares={len(sq)} pairs={npair} walls={nwall}',
                dict(squares=[names[k] for k in sq], max_tilt_deg=math.degrees(mt)))
    # axis-parallel core: split the pair contacts by separating direction
    edges = {'x': [], 'y': [], 'xy': []}
    for c in keep:
        if c[0] != 'p':
            continue
        dirs = set()
        for (d, _rows) in pairs[(c[1], c[2])]:
            dirs.add('x' if abs(d[0]) > abs(d[1]) else 'y')
        edges[''.join(sorted(dirs))].append((c[1], c[2]))
    anchors = {'x': {}, 'y': {}}
    for c in keep:
        if c[0] == 'w':
            ax = 'x' if c[2] < 2 else 'y'
            anchors[ax].setdefault(c[1], set()).add(WALLN[c[2]])
    parts, detail = [], []
    for ax in ('x', 'y'):
        adj = {}
        for (i, j) in edges[ax]:
            adj.setdefault(i, set()).add(j)
            adj.setdefault(j, set()).add(i)
        for k in anchors[ax]:
            adj.setdefault(k, set())
        seen = set()
        for v0 in sorted(adj):
            if v0 in seen:
                continue
            comp, st = set(), [v0]
            while st:
                v = st.pop()
                if v not in comp:
                    comp.add(v)
                    st += list(adj[v])
            seen |= comp
            an = ''.join(sorted(set().union(*[anchors[ax].get(v, set()) for v in comp])))
            parts.append(f'{ax}{len(comp)}{an}')
            detail.append(f'{ax}:{"-".join(names[v] for v in sorted(comp, key=lambda q: z[(1 if ax == 'x' else 2) + 3 * q]))}[{an}]')
    if edges['xy']:
        parts.append(f'corner{len(edges["xy"])}')
        detail.append('corner:' + ','.join(f'{names[i]}/{names[j]}' for i, j in edges['xy']))
    return ' + '.join(sorted(parts)), dict(chains=detail)


# ===================================================================== io
def load_any(path, zero_only, near):
    out = []
    if path.endswith('.jsonl'):
        for line in open(path):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if zero_only and abs(r['best']) > near:
                continue
            z = np.array(r['z'])
            out.append((f'class{r["id"]}:{r["tag"]}', z, [f'S{k}' for k in range(12)]))
    else:
        for (tag, val, rows) in rank8.load(path):
            if zero_only and abs(val) > near:
                continue
            z = rank8.z_from([(x, y, th) for (_n, x, y, th) in rows])
            out.append((tag, z, [r[0] for r in rows]))
    return out


_A = {}


def _init(a):
    _A['a'] = a


def _one(item):
    a = _A['a']
    tag, z, names = item
    n = (len(z) - 1) // 3
    rng = np.random.default_rng(abs(hash(tag)) % (2 ** 32))
    val = true_value(z, n)
    pairs, walls = structure(z, n, a.tol)
    d, v = first_order(n, pairs, walls)
    rec = dict(tag=tag, value=val, npairs=len(pairs), nwalls=len(walls), delta=d,
               max_tilt_deg=math.degrees(float(np.max(np.abs(norm_tilt(z[3::3]))))))
    if not (d <= a.zero):
        ls, eps = line_search(z, v, n) if v is not None else (val, 0.0)
        rec.update(sig='FIRST-ORDER FEASIBLE', line_search=ls, eps=eps)
        return rec
    keep = core(z, n, pairs, walls, rng, a.ntry, a.zero)
    sig, det = signature(z, n, keep, pairs, names)
    rec.update(sig=sig, core_pairs=sum(1 for c in keep if c[0] == 'p'),
               core_walls=sum(1 for c in keep if c[0] == 'w'), **det)
    return rec


def cmd_analyse(a):
    import multiprocessing as mp
    items = []
    for p in a.FILES:
        items += load_any(p, not a.all, a.near)
    if a.limit:
        items = items[:a.limit]
    print(f'# {len(items)} configurations (|value| <= {a.near}) from {len(a.FILES)} file(s); '
          f'active tol {a.tol}', flush=True)
    with mp.Pool(a.nproc, initializer=_init, initargs=(a,)) as pool:
        recs = pool.map(_one, items, chunksize=4)
    if a.out:
        with open(a.out, 'w') as f:
            for r in recs:
                f.write(json.dumps(r) + '\n')
    feas = [r for r in recs if r['sig'] == 'FIRST-ORDER FEASIBLE']
    print(f'\nfirst-order FEASIBLE at {len(feas)} of {len(recs)} configurations')
    for r in sorted(feas, key=lambda r: -r['line_search'])[:15]:
        print(f'  {r["tag"]:<22} delta {r["delta"]:.3e}  value {r["value"]:+.3e} -> line search '
              f'{r["line_search"]:+.3e} at eps {r["eps"]:.1e}')
    cnt = Counter(r['sig'] for r in recs)
    print('\ncore signatures (xN = N squares chained along x; L/R/B/T = wall anchors):')
    for sig, k in cnt.most_common():
        ex = next(r for r in recs if r['sig'] == sig)
        tl = [r['max_tilt_deg'] for r in recs if r['sig'] == sig]
        print(f'  {k:6d}  {sig:<44} config tilt up to {max(tl):6.2f} deg   e.g. {ex["tag"]}: '
              + '; '.join(ex.get('chains', ex.get('squares', [])) or []))
    rigid = [r for r in recs if r['sig'] != 'FIRST-ORDER FEASIBLE']
    simple = [r for r in rigid if r['sig'] in ('x4LR', 'y4BT')]
    tilted = [r for r in rigid if r['sig'].startswith('TILTED')]
    print(f'\nof {len(rigid)} rigid configurations: {len(simple)} have a single straight '
          f'wall-to-wall row of four as core,\n  {len(tilted)} have no axis-parallel core at '
          f'all, {len(rigid) - len(simple) - len(tilted)} have a compound axis-parallel core.')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('analyse')
    c.add_argument('FILES', nargs='+')
    c.add_argument('--tol', type=float, default=1e-7, help='a constraint is active below this')
    c.add_argument('--zero', type=float, default=1e-6)
    c.add_argument('--near', type=float, default=1e-7)
    c.add_argument('--all', action='store_true')
    c.add_argument('--ntry', type=int, default=3)
    c.add_argument('--limit', type=int, default=0)
    c.add_argument('--nproc', type=int, default=4)
    c.add_argument('--out', default=None)
    a = ap.parse_args()
    cmd_analyse(a)


if __name__ == '__main__':
    main()
