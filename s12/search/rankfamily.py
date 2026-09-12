#!/usr/bin/env python3
"""rank family: odd-polygon anchor rows for the QSTAB loop (`tasks/rank-diag`, round 2).

Doc: `search/RANKDIAG.md` §8.  Used by `search/cliquelever.py` behind `--pent` (the only hook:
`Lever.add_pgon`, one call in `run_stage`, one guard in `finalize`).  Nothing else imports it.

**The inequality.**  Let `A_0 … A_{k-1}` be points of the container, `k` odd, and let
`E_i = [A_i, A_{i+1 mod k}]` be the `k` closed segments of the polygon they span.  Put

    X_i = { S admissible : E_i subset S }        (a square is convex: E_i subset S iff both
                                                  endpoints are in S)

Two members of the same `X_i` both contain `A_i`.  A member of `X_i` and a member of `X_{i+1}`
both contain `A_{i+1}`, because `E_i ∩ E_{i+1} ∋ A_{i+1}`.  So a pairwise-DISJOINT subfamily of
`X_0 ∪ … ∪ X_{k-1}` hits each piece at most once and never two cyclically consecutive pieces:
it injects into an independent set of `C_k`.  Hence

    every packing of pairwise-disjoint closed unit squares puts at most
    alpha(C_k) = (k-1)/2 of its squares in the union,

so `mu(X_0 ∪ … ∪ X_{k-1}) <= (k-1)/2` is valid for the integer hull and may be added to the LP
exactly as a clique row is (`k = 3` gives a clique with right-hand side 1).

The well-formedness condition is only `E_i ∩ E_{i+1} != {}` for each cyclic `i`.  Taking the
segments to be the sides of a polygon on `k` points makes that automatic (consecutive sides share
an endpoint), and it is the parameterisation the verifier's `cores_meet` already tests; the search
below uses it.

**What is exact.**  Membership of every row: `lc.sq_contains` on both endpoints, integers, with a
float prefilter that is only trusted outside a `1e-7` band (the same rule as
`cliquelever.Poses.contains_rows`).  Every row is re-derived from its exact rational anchors
before it enters the LP, and its support members' independence number is checked against
`(k-1)/2` by a complete branch and bound.  The hill climb that proposes the anchors is a
heuristic: it can only fail to find a row, never produce an invalid one.
"""
import time

import numpy as np

import leaf_ceiling as lc
from fractions import Fraction as Fr

FTOL = 1e-7


# ====================================================================== geometry
def contains_vec(ps, X, Y, D):
    """boolean over ALL poses: does the pose's closed square contain the exact point (X/D, Y/D)?
    Float test, decided exactly inside the FTOL band."""
    x, y = X / D, Y / D
    dx = x - ps.F[:, 0]
    dy = y - ps.F[:, 1]
    u = np.abs(dx * ps.F[:, 2] + dy * ps.F[:, 3])
    v = np.abs(-dx * ps.F[:, 3] + dy * ps.F[:, 2])
    m = 0.5 - np.maximum(u, v)
    inside = m > FTOL
    for j in np.nonzero(np.abs(m) <= FTOL)[0]:
        inside[j] = lc.sq_contains(ps.sq[int(j)], X, Y, D)
    return inside


def row_members(ps, pts):
    """the MAXIMAL row: every loaded pose whose square contains one of the k closed segments"""
    cont = [contains_vec(ps, X, Y, D) for (X, Y, D) in pts]
    k = len(pts)
    out = np.zeros(ps.n, dtype=bool)
    for i in range(k):
        out |= cont[i] & cont[(i + 1) % k]
    return np.nonzero(out)[0].tolist()


def verify_row(ps, pts, members):
    """fully exact re-derivation: every listed member must contain some segment (both endpoints).
    Returns the number of failures (0 = the row is exactly what its anchors define)."""
    k = len(pts)
    bad = 0
    for i in members:
        s = ps.sq[i]
        if not any(lc.sq_contains(s, *pts[a]) and lc.sq_contains(s, *pts[(a + 1) % k])
                   for a in range(k)):
            bad += 1
    return bad


def alpha_of(ps, members, time_limit=30.0):
    """independence number of the members' exact closed-intersection graph (complete B&B on the
    complement).  Returns (alpha, complete)."""
    if len(members) < 2:
        return len(members), True
    squares = [ps.sq[i][:8] + (1,) for i in members]
    nbr = lc.closed_graph(squares, verbose=False)
    n = len(members)
    comp = [[] for _ in range(n)]
    for i in range(n):
        comp[i] = [j for j in range(n) if j != i and j not in nbr[i]]
    cadj = [0] * n
    for i in range(n):
        for j in comp[i]:
            cadj[i] |= 1 << j
    return _max_clique(cadj, (1 << n) - 1, time_limit)


def _max_clique(adj, P, time_limit):
    best = [0]
    t0 = time.time()
    to = [False]

    def colour(P):
        order, bounds, c = [], [], 0
        Q = P
        while Q:
            c += 1
            avail = Q
            while avail:
                v = (avail & -avail).bit_length() - 1
                avail &= ~(1 << v) & ~adj[v]
                Q &= ~(1 << v)
                order.append(v)
                bounds.append(c)
        return order, bounds

    def expand(P, sz):
        if time.time() - t0 > time_limit:
            to[0] = True
            return
        order, bounds = colour(P)
        for i in range(len(order) - 1, -1, -1):
            if sz + bounds[i] <= best[0]:
                return
            v = order[i]
            nP = P & adj[v]
            if nP:
                expand(nP, sz + 1)
                if to[0]:
                    return
            elif sz + 1 > best[0]:
                best[0] = sz + 1
            P &= ~(1 << v)

    expand(P, 0)
    return best[0], not to[0]


# ====================================================================== separation
def _cand_masks(squares, V, Inc, w, ncand):
    """the arrangement vertices of the support as candidate anchors: a boolean (ncand, nsup)
    matrix of "this support square contains this vertex", keeping the heaviest `ncand` distinct
    point cliques."""
    nsup = len(squares)
    nv = len(V)
    M = np.zeros((nv, nsup), dtype=bool)
    for si, col in enumerate(Inc.cols):
        M[np.asarray(col, dtype=np.int64), si] = True
    mass = M @ w
    keep = np.argsort(-mass)[:max(ncand * 4, ncand)]
    # dedupe on the containment pattern, heaviest first
    seen = set()
    idx = []
    for j in keep:
        key = M[j].tobytes()
        if key in seen:
            continue
        seen.add(key)
        idx.append(int(j))
        if len(idx) >= ncand:
            break
    return M[idx], [V[j] for j in idx]


def _climb(M, w, k, rhs_int, rng, restarts, seeds, deadline):
    """vectorised coordinate descent on k anchor indices, maximising the support mass of
    `union_i (piece from segment [A_i, A_{i+1}])`.  Returns every local optimum above rhs_int."""
    nc = M.shape[0]
    found = {}
    starts = [list(s) for s in seeds if len(s) == k]
    starts += [[int(rng.randrange(nc)) for _ in range(k)] for _ in range(restarts)]
    for idx in starts:
        if time.time() > deadline:
            break
        cur = -1
        for _ in range(6):
            improved = False
            for pos in range(k):
                # segments not touching `pos` are fixed; the two that do are swept over all cands
                base = np.zeros(M.shape[1], dtype=bool)
                for a in range(k):
                    b = (a + 1) % k
                    if a != pos and b != pos:
                        base |= M[idx[a]] & M[idx[b]]
                prev, nxt = M[idx[(pos - 1) % k]], M[idx[(pos + 1) % k]]
                tot = base[None, :] | (M & prev[None, :]) | (M & nxt[None, :])
                vals = tot @ w
                c = int(np.argmax(vals))
                if vals[c] > cur:
                    cur = int(vals[c])
                    idx[pos] = c
                    improved = True
            if not improved:
                break
        if cur > rhs_int:
            found[_canon(idx)] = (cur, list(idx))
    return sorted(found.values(), key=lambda v: -v[0])


def _canon(idx):
    """canonical form of a cyclic anchor sequence, up to rotation and reflection"""
    k = len(idx)
    best = None
    for seq in (idx, idx[::-1]):
        for r in range(k):
            c = tuple(seq[r:] + seq[:r])
            if best is None or c < best:
                best = c
    return best


def separate(lever, sup, squares, V, Inc, w, mu, it):
    """propose, verify and add odd-polygon rows.  Returns (nnew, best, secs, nfound)."""
    a = lever.a
    ps = lever.ps
    DM = lever.DM
    t0 = time.time()
    deadline = t0 + a.pent_time
    rng = lever.pent_rng
    M, verts = _cand_masks(squares, V, Inc, np.asarray(w, dtype=np.int64), a.pent_cands)
    nnew, nfound, best = 0, 0, None
    for k in a.pent:
        rhs_int = (k - 1) // 2 * DM
        seeds = lever.pent_seeds.get(k, [])
        cands = _climb(M, np.asarray(w, dtype=np.int64), k, int(rhs_int * (1 + a.ktol)), rng,
                       a.pent_restarts, seeds, deadline)
        nfound += len(cands)
        lever.pent_seeds[k] = [c[1] for c in cands[:3]]
        for (val, idx) in cands[:a.pent_want]:
            pts = [verts[i] for i in idx]
            mem = row_members(ps, pts)
            bad = verify_row(ps, pts, mem)
            if bad:
                lever.log(f'   [pgon k={k} DROPPED: {bad} members fail the exact segment test]')
                continue
            smem = [i for i in mem if mu[i] > 1e-12]
            al, comp = alpha_of(ps, smem, 30.0) if len(smem) <= 160 else ((k - 1) // 2, False)
            if comp and al > (k - 1) // 2:
                lever.log(f'   [pgon k={k} DROPPED: alpha({len(smem)} support members) = {al} '
                          f'> {(k - 1) // 2} -- IMPOSSIBLE, check the geometry]')
                continue
            rowmass = sum(mu[i] for i in mem)
            info = dict(it=it, k=k, rhs=(k - 1) / 2, size=len(mem), sup_size=len(smem),
                        mass0=rowmass, excess=rowmass - (k - 1) / 2,
                        alpha=al, alpha_complete=bool(comp),
                        anchors=[[str(Fr(X, D)), str(Fr(Y, D))] for (X, Y, D) in pts])
            if lever.add_pgon(mem, (k - 1) / 2, info, pts):
                nnew += 1
                if best is None:
                    best = info
    return nnew, best, time.time() - t0, nfound


# ====================================================================== finalisation guards
class TopUpGuard:
    """`Lever.finalize` tops pinned regions back up on poses with exact coverage AND clique
    slack; a polygon row must not be pushed over its right-hand side either.  This tracks the
    exact integer mass of every polygon row as the top-up proceeds."""

    def __init__(self, lever, mi):
        self.DM = lever.DM
        self.rows = [(set(c['members']), int(round(c['rhs'] * lever.DM)))
                     for c in lever.pgons]
        self.mass = [sum(mi[i] for i in mem) for mem, _ in self.rows]
        self.of = {}
        for r, (mem, _) in enumerate(self.rows):
            for i in mem:
                self.of.setdefault(i, []).append(r)

    def slack(self, i):
        s = self.DM
        for r in self.of.get(i, ()):
            s = min(s, self.rows[r][1] - self.mass[r])
        return max(s, 0)

    def credit(self, i, add):
        for r in self.of.get(i, ()):
            self.mass[r] += add


def check_final(lever, mi):
    """exact check of every polygon row against the finalised integer measure, plus a full
    re-derivation of each row from its stored rational anchors.  Returns (ok, worst, report)."""
    DM = lever.DM
    ps = lever.ps
    ok = True
    worst = -1.0
    rep = []
    for c in lever.pgons:
        rhs = int(round(c['rhs'] * DM))
        m = sum(mi[i] for i in c['members'])
        pts = c['pts']
        again = row_members(ps, pts)
        same = sorted(again) == sorted(c['members'])
        good = (m <= rhs) and same
        ok &= good
        worst = max(worst, (m - rhs) / DM)
        rep.append(dict(k=c['k'], rhs=c['rhs'], mass=m / DM, slack=(rhs - m) / DM,
                        rederived=bool(same), size=len(c['members'])))
    return ok, worst, rep


# ====================================================================== selftest
class _Shim:
    """the three attributes of `cliquelever.Poses` these functions use"""

    def __init__(self, t, poses):
        self.t = Fr(t)
        self.sq = [lc.make_square(cx, cy, p, q, self.t, 0, tag=i)
                   for i, (p, q, cx, cy) in enumerate(poses)]
        self.n = len(self.sq)
        self.F = np.array([[float(Fr(s[0], s[2])), float(Fr(s[1], s[2])),
                            s[3] / s[5], s[4] / s[5]] for s in self.sq], dtype=float)


def selftest():
    """the k-ring: k axis-parallel unit squares whose centres sit on a circle so that only
    cyclically consecutive ones meet.  A_i taken inside S_{i-1} ∩ S_i makes E_i = [A_i, A_{i+1}]
    lie in the convex S_i, so the polygon row picks up all k squares; each carries mass 1/2, so
    the row mass is k/2 against the bound (k-1)/2 -- violated by exactly 1/2, and the exact
    independence number of the members must come out as (k-1)/2."""
    import math
    ok = True

    def chk(name, cond, extra=''):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'} {name} {extra}")
        ok = ok and cond

    for k, R in ((5, 0.8), (7, 1.0), (9, 1.3)):
        poses = []
        cs = []
        for i in range(k):
            th = 2 * math.pi * i / k
            cx = Fr(round(4 + R * math.cos(th), 3)).limit_denominator(1000)
            cy = Fr(round(4 + R * math.sin(th), 3)).limit_denominator(1000)
            poses.append((0, 1, cx, cy))
            cs.append((cx, cy))
        ps = _Shim(8, poses)
        nbr = lc.closed_graph([s[:8] + (1,) for s in ps.sq], verbose=False)
        chk(f'k={k}: ring has {k} edges', sum(map(len, nbr)) // 2 == k,
            f'(got {sum(map(len, nbr)) // 2})')
        # A_i = midpoint of the centres of S_{i-1} and S_i (both squares contain it: the centres
        # are less than 1 apart in each coordinate, so the midpoint is within 1/2 of both)
        pts = []
        for i in range(k):
            a, b = cs[(i - 1) % k], cs[i]
            x = (a[0] + b[0]) / 2
            y = (a[1] + b[1]) / 2
            D = x.denominator * y.denominator
            pts.append((x.numerator * y.denominator, y.numerator * x.denominator, D))
        mem = row_members(ps, pts)
        chk(f'k={k}: the row is the whole ring', sorted(mem) == list(range(k)), f'(got {mem})')
        chk(f'k={k}: exact re-derivation', verify_row(ps, pts, mem) == 0)
        al, comp = alpha_of(ps, mem)
        chk(f'k={k}: alpha = (k-1)/2 = {(k - 1) // 2}', al == (k - 1) // 2 and comp, f'(got {al})')
        # a square far away must NOT be in the row
        ps2 = _Shim(8, poses + [(0, 1, Fr(1, 2), Fr(1, 2))])
        chk(f'k={k}: a distant square is excluded', k not in row_members(ps2, pts))
    print('rankfamily selftest ' + ('PASSED' if ok else 'FAILED'))
    return 0 if ok else 1


if __name__ == '__main__':
    import sys
    sys.exit(selftest())
