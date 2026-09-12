#!/usr/bin/env python3
"""rank-diag: what non-clique structure carries the residual excess of the QSTAB measures?

Task: `tasks/rank-diag/README.md`.  Doc: `search/RANKDIAG.md`.

On a converged QSTAB measure (a `leaf_ceiling.py` measure file) this computes, with the exact
geometry of `search/leaf_ceiling.py` (`read_measure`, `build_squares`, `closed_graph`: integer
separating-axis SAT with `<=`, so touching counts):

  alpha     the independence number of the support's closed-intersection graph `G`, exactly,
            and the gap `mu(V) - alpha(G)` that no clique inequality can explain;
  cycles    every chordless odd cycle of length 5 and 7 with mass `> (k-1)/2`, every odd
            antihole on 7 vertices with mass `> 2`, every 5-wheel (hub adjacent to a chordless
            C5) with hub mass + rim mass `> 2`, with anatomy;
  search    a heuristic hunt for subsets `X` with `mu(X) - alpha(G[X])` large, and for
            *minimal* violated rank subsets (no single member removable) -- the question being
            whether a small certifiable witness exists at all; everything found goes onto a
            size-vs-excess Pareto front;
  anchors   step 4: the separation problem for the certifiable rank-2 family of five anchor
            points `A_0..A_4` and the five segments `E_i = [A_i, A_{i+1}]` of the pentagon they
            span -- the pieces `X_i = {S : E_i in S}` have guaranteed-meet graph `C5`, so
            `mu(union X_i) <= alpha(C5) = 2` by Lemma 0 of `notes/clique-family.md`.

All mass arithmetic is on the integer masses over the measure's common denominator `DM`, so
every comparison (`mu(X) > alpha`, `mu(C) > (k-1)/2`, `mu(pentagon) > 2`) is exact.  The only
floats are the reported decimals, the heuristics' choice of which subsets to try (every set
they return is verified exactly afterwards), the centre-distance prefilter inside
`closed_graph` (which only decides which pairs get the exact test), and the printed degrees.

Reproduce: see `search/RANKDIAG.md` 6.
"""
import argparse
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import leaf_ceiling as lc                                             # noqa: E402

GRID = [(Fr(i), Fr(j)) for i in (1, 2, 3) for j in (1, 2, 3)]


def log(msg=''):
    print(msg, flush=True)


# ============================================================================ loading
class Support:
    """the support of a measure: exact squares, integer masses, closed-intersection graph"""

    def __init__(self, path, r='1', t_override=None):
        self.path = path
        self.name = os.path.basename(path).replace('_exact.txt', '')
        labels = []
        t, sym, poses = lc.read_measure(path, t_override=t_override, labels=labels)
        self.t, self.sym = t, sym
        squares, DM, pose_of = lc.build_squares(t, sym, poses)
        self.squares, self.DM, self.pose_of = squares, DM, pose_of
        self.n = n = len(squares)
        self.w = [s[8] for s in squares]                       # integer mass over DM
        self.total = sum(self.w)
        self.cx = [Fr(s[0], s[2]) for s in squares]
        self.cy = [Fr(s[1], s[2]) for s in squares]
        self.deg = [math.degrees(math.atan2(s[4], s[3])) % 90.0 for s in squares]
        boxes = lc.region_boxes(t, Fr(r))
        self.region = [','.join(lc.regions_of(self.cx[i], self.cy[i], boxes)) for i in range(n)]
        nbr = lc.closed_graph(squares, verbose=False)
        self.nbr = nbr
        self.adj = [0] * n
        for i in range(n):
            m = 0
            for j in nbr[i]:
                m |= 1 << j
            self.adj[i] = m
        self.full = (1 << n) - 1
        self.cadj = [(~self.adj[i]) & self.full & ~(1 << i) for i in range(n)]   # complement
        self.nedge = sum(len(s) for s in nbr) // 2
        # which grid vertices each square contains (exact)
        self.grid = []
        for i, s in enumerate(squares):
            g = 0
            for k, (gx, gy) in enumerate(GRID):
                if lc.sq_contains(s, gx.numerator * gy.denominator, gy.numerator * gx.denominator,
                                  gx.denominator * gy.denominator):
                    g |= 1 << k
            self.grid.append(g)

    def mass(self, bits):
        return sum(self.w[v] for v in bits_of(bits))

    def f(self, m_int):
        return m_int / self.DM


def bits_of(x):
    while x:
        b = x & -x
        yield b.bit_length() - 1
        x ^= b


def popcount(x):
    return bin(x).count('1')


# ============================================================================ max clique / alpha
def colour_sort(adj, P):
    """greedy colouring of P; returns (order, bounds) with bounds non-decreasing (Tomita)"""
    order, bounds = [], []
    colour = 0
    Q = P
    while Q:
        colour += 1
        avail = Q
        while avail:
            v = (avail & -avail).bit_length() - 1
            avail &= ~(1 << v)
            avail &= ~adj[v]
            Q &= ~(1 << v)
            order.append(v)
            bounds.append(colour)
    return order, bounds


def max_clique_bits(adj, P, time_limit=300.0, lb=0):
    """max cardinality clique of the graph `adj` restricted to the vertex mask `P`.
    Returns (size, members_mask, complete)."""
    best = [lb, 0]
    t0 = time.time()
    state = {'timeout': False}

    def expand(P, Rsz, R):
        if time.time() - t0 > time_limit:
            state['timeout'] = True
            return
        order, bounds = colour_sort(adj, P)
        for i in range(len(order) - 1, -1, -1):
            if Rsz + bounds[i] <= best[0]:
                return
            v = order[i]
            newP = P & adj[v]
            if newP:
                expand(newP, Rsz + 1, R | (1 << v))
                if state['timeout']:
                    return
            elif Rsz + 1 > best[0]:
                best[0] = Rsz + 1
                best[1] = R | (1 << v)
            P &= ~(1 << v)

    expand(P, 0, 0)
    return best[0], best[1], not state['timeout']


def alpha(sup, P=None, time_limit=300.0):
    """independence number of G[P] (exact): max clique in the complement"""
    if P is None:
        P = sup.full
    return max_clique_bits(sup.cadj, P, time_limit=time_limit)


# ============================================================================ odd structures
def chordless_odd_cycles(adj, n, k, w, thresh_int, order=None):
    """every chordless cycle of length k (odd) whose integer mass exceeds `thresh_int`.

    Canonical form: `v0` is the cycle member of smallest rank in `order` (order is by mass
    descending, so v0 is the heaviest member), and rank(v1) < rank(v_{k-1}) picks one of the
    two traversal directions.  Prune: the partial path can gain at most (k - len) * w[v0].
    Returns a list of vertex tuples in cycle order."""
    if order is None:
        order = sorted(range(n), key=lambda v: (-w[v], v))
    rank = [0] * n
    for i, v in enumerate(order):
        rank[v] = i
    out = []

    for v0 in order:
        if k * w[v0] <= thresh_int:
            break                                     # order is mass-descending: none left
        allowed = 0
        for v in order[rank[v0] + 1:]:
            allowed |= 1 << v
        n0 = adj[v0] & allowed

        def rec(path, cur, blocked):
            i = len(path) - 1                          # index of the last placed vertex
            last = path[-1]
            if i == k - 2:                             # place the closing vertex
                cands = adj[last] & allowed & ~blocked & n0
                lo = rank[path[1]]
                for u in bits_of(cands):
                    if rank[u] <= lo:
                        continue
                    if cur + w[u] > thresh_int:
                        out.append(tuple(path) + (u,))
                return
            if cur + (k - 1 - i) * w[v0] <= thresh_int:
                return
            cands = adj[last] & allowed & ~blocked & ~n0
            for u in bits_of(cands):
                # the next vertex must miss N(v_1..v_i): add the neighbourhood of the OLD last
                rec(path + [u], cur + w[u], blocked | adj[last] | (1 << u))

        for v1 in bits_of(n0):
            rec([v0, v1], w[v0] + w[v1], (1 << v0) | (1 << v1))
    return out


def wheels(sup, c5s, thresh_int):
    """5-wheels: hub h adjacent to every vertex of a chordless C5, hub + rim mass > thresh"""
    out = []
    for cyc in c5s:
        common = sup.full
        rim = 0
        for v in cyc:
            common &= sup.adj[v]
            rim |= 1 << v
        rm = sum(sup.w[v] for v in cyc)
        for h in bits_of(common & ~rim):
            if rm + sup.w[h] > thresh_int:
                out.append((h, cyc))
    return out


# ============================================================================ anatomy
def anatomy(sup, members, kind, bound_int, cyc=None, hub=None):
    """members: iterable of vertex indices"""
    ms = list(members)
    m_int = sum(sup.w[v] for v in ms)
    gx = sup.full and 0
    common_grid = None
    for v in ms:
        common_grid = sup.grid[v] if common_grid is None else (common_grid & sup.grid[v])
    grid_all = [f'({int(GRID[k][0])},{int(GRID[k][1])})' for k in bits_of(common_grid or 0)]
    # grid vertices covered by a majority
    cnt = {}
    for v in ms:
        for k in bits_of(sup.grid[v]):
            cnt[k] = cnt.get(k, 0) + 1
    best_grid = sorted(cnt.items(), key=lambda kv: -kv[1])[:2]
    regs = {}
    for v in ms:
        regs[sup.region[v]] = regs.get(sup.region[v], 0) + sup.w[v]
    angles = sorted({round(sup.deg[v], 2) for v in ms})
    return {
        'kind': kind,
        'members': ms,
        'size': len(ms),
        'mass': sup.f(m_int),
        'mass_int': m_int,
        'bound': sup.f(bound_int),
        'excess': sup.f(m_int - bound_int),
        'cycle': list(cyc) if cyc else None,
        'hub': hub,
        'centres': [(float(sup.cx[v]), float(sup.cy[v])) for v in ms],
        'angles_deg': angles,
        'masses': [sup.f(sup.w[v]) for v in ms],
        'regions': {k: sup.f(v) for k, v in sorted(regs.items(), key=lambda kv: -kv[1])},
        'common_grid_vertex': grid_all,
        'grid_vertex_counts': [(f'({int(GRID[k][0])},{int(GRID[k][1])})', c) for k, c in best_grid],
    }


def fmt_struct(sup, a):
    c = ' '.join(f'({x:.3f},{y:.3f})' for x, y in a['centres'])
    return (f"    {a['kind']} |X|={a['size']} mass={a['mass']:.6f} bound={a['bound']:.1f} "
            f"excess={a['excess']:+.6f}\n"
            f"      members {a['members']}  masses "
            f"{['%.4f' % m for m in a['masses']]}\n"
            f"      centres {c}\n"
            f"      angles  {a['angles_deg']}\n"
            f"      regions {{{', '.join('%s %.3f' % (k, v) for k, v in a['regions'].items())}}}"
            f"  common grid vertex: {a['common_grid_vertex'] or 'none'}"
            f"  top grid: {a['grid_vertex_counts']}")


# ============================================================================ commands
def cmd_alpha(sup, args, res):
    t0 = time.time()
    a, mask, complete = alpha(sup, time_limit=args.time)
    res['n'] = sup.n
    res['edges'] = sup.nedge
    res['density'] = 2 * sup.nedge / (sup.n * (sup.n - 1))
    res['mass'] = sup.f(sup.total)
    res['alpha'] = a
    res['alpha_complete'] = complete
    res['alpha_witness'] = list(bits_of(mask))
    res['gap'] = sup.f(sup.total - a * sup.DM)
    res['alpha_secs'] = time.time() - t0
    log(f"  n = {sup.n} poses, {sup.nedge} edges (density {res['density']:.3f}), "
        f"mass = {res['mass']:.6f}")
    log(f"  alpha(G) = {a}  (complete B&B: {complete}, {res['alpha_secs']:.1f} s)   "
        f"gap = mass - alpha = {res['gap']:.6f}")
    wit = res['alpha_witness']
    log(f"    a maximum independent set: {wit}")
    log(f"      centres " + ' '.join(f'({float(sup.cx[v]):.3f},{float(sup.cy[v]):.3f})'
                                     for v in wit))
    log(f"      masses  " + ' '.join(f'{sup.f(sup.w[v]):.3f}' for v in wit))
    return a


def descend(gen, bound_int, cap=200):
    """run `gen(threshold)` at the rank bound, then at 0.9, 0.8, ... of it, stopping at the
    first threshold that yields something (so the near misses are visible when nothing is
    violated).  Returns (threshold_used, results)."""
    fracs = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
    for f in fracs:
        thr = int(bound_int * f)
        r = gen(thr)
        if r:
            return thr, (r if len(r) <= cap else r)
    return 0, []


def cmd_cycles(sup, args, res):
    DM = sup.DM
    out = {}
    near = {}
    # ---- chordless odd holes C5 (bound 2) and C7 (bound 3)
    for k in (5, 7):
        t0 = time.time()
        bnd = (k - 1) // 2 * DM
        cyc = chordless_odd_cycles(sup.adj, sup.n, k, sup.w, bnd)
        an = [anatomy(sup, c, f'C{k}', bnd, cyc=c) for c in cyc]
        an.sort(key=lambda a: -a['excess'])
        out[f'C{k}'] = an
        log(f"  chordless C{k} with mass > {(k - 1) // 2}: {len(an)} ({time.time() - t0:.1f} s)")
        for a in an[:10]:
            log(fmt_struct(sup, a))
        if not an:                                   # how close does the family get?
            thr, r = descend(lambda T: chordless_odd_cycles(sup.adj, sup.n, k, sup.w, T), bnd)
            if r:
                best = max(sum(sup.w[v] for v in c) for c in r)
                near[f'C{k}'] = sup.f(best)
                log(f"    heaviest chordless C{k} anywhere: mass {sup.f(best):.6f} "
                    f"(bound {(k - 1) // 2}; {len(r)} above {sup.f(thr):.3f})")
    # ---- odd antiholes on 7 vertices (alpha = 2)
    t0 = time.time()
    bnd = 2 * DM
    ah = chordless_odd_cycles(sup.cadj, sup.n, 7, sup.w, bnd)
    an = [anatomy(sup, c, 'antihole7', bnd, cyc=c) for c in ah]
    an.sort(key=lambda a: -a['excess'])
    out['antihole7'] = an
    log(f"  odd antiholes on 7 vertices with mass > 2: {len(an)} ({time.time() - t0:.1f} s)")
    for a in an[:10]:
        log(fmt_struct(sup, a))
    if not an:
        thr, r = descend(lambda T: chordless_odd_cycles(sup.cadj, sup.n, 7, sup.w, T), bnd)
        if r:
            best = max(sum(sup.w[v] for v in c) for c in r)
            near['antihole7'] = sup.f(best)
            log(f"    heaviest 7-antihole anywhere: mass {sup.f(best):.6f} (bound 2; "
                f"{len(r)} above {sup.f(thr):.3f})")
    # ---- 5-wheels (alpha = 2): need ALL chordless C5 whose rim can still reach 2 with a hub
    t0 = time.time()
    wmax = max(sup.w)
    rim_thr = max(0, 2 * DM - wmax)
    all5 = chordless_odd_cycles(sup.adj, sup.n, 5, sup.w, rim_thr)
    wh = wheels(sup, all5, 2 * DM)
    an = [anatomy(sup, list(c) + [h], 'wheel5', 2 * DM, cyc=c, hub=h) for (h, c) in wh]
    an.sort(key=lambda a: -a['excess'])
    out['wheel5'] = an
    log(f"  chordless C5 with rim mass > {rim_thr / DM:.4f} (wheel candidates): {len(all5)}; "
        f"5-wheels with hub+rim > 2: {len(an)} ({time.time() - t0:.1f} s)")
    for a in an[:10]:
        log(fmt_struct(sup, a))
    if not an:
        best = 0
        for c in all5:
            common = sup.full
            for v in c:
                common &= sup.adj[v]
            common &= ~sum(1 << v for v in c)
            if common:
                best = max(best, sum(sup.w[v] for v in c) + max(sup.w[h]
                                                                for h in bits_of(common)))
        if best:
            near['wheel5'] = sup.f(best)
            log(f"    heaviest 5-wheel anywhere: hub+rim {sup.f(best):.6f} (bound 2)")
    res['cycles'] = {k: [strip(a) for a in v[:10]] for k, v in out.items()}
    res['cycle_counts'] = {k: len(v) for k, v in out.items()}
    res['cycle_best_excess'] = {k: (v[0]['excess'] if v else 0.0) for k, v in out.items()}
    res['cycle_total_excess'] = {k: sum(a['excess'] for a in v) for k, v in out.items()}
    res['cycle_near_miss'] = near
    return out


def strip(a):
    b = dict(a)
    b.pop('centres', None)
    return {**b, 'centres': [[round(x, 6), round(y, 6)] for x, y in a['centres']]}


class Front:
    """the size-vs-excess Pareto front: for every |X| seen, the largest `mu(X) - alpha(G[X])`.
    Every violated set the search touches is offered to it, so the front answers the question
    the brief really asks -- how small can a witness be and still carry real excess."""

    def __init__(self, sup):
        self.sup = sup
        self.best = {}                       # size -> (excess_int, mass_int, alpha, mask)

    def offer(self, mask, m_int, a):
        k = popcount(mask)
        ex = m_int - a * self.sup.DM
        if ex <= 0:
            return
        cur = self.best.get(k)
        if cur is None or ex > cur[0]:
            self.best[k] = (ex, m_int, a, mask)

    def rows(self):
        """the monotone front: the sizes at which the achievable excess strictly improves"""
        out, seen = [], -1
        for k in sorted(self.best):
            ex, m, a, mask = self.best[k]
            if ex > seen:
                seen = ex
                out.append({'size': k, 'excess': self.sup.f(ex), 'mass': self.sup.f(m),
                            'alpha': a, 'members': list(bits_of(mask))})
        return out

    def smallest_with(self, ex_int):
        for k in sorted(self.best):
            if self.best[k][0] >= ex_int:
                return k, self.best[k]
        return None, None


def cmd_search(sup, args, res):
    """heuristic hunt for X maximising mu(X) - alpha(G[X]), and for MINIMAL violated X.

    The key economy: if `W` is a maximum independent set of `G[X]` and `v` not in `W`, then
    `W ⊆ X\\v` so `alpha(X\\v) = alpha(X)`.  Only the `alpha(X)` members of `W` can make alpha
    fall, so each round needs one alpha call for the witness plus `|W| <= 11` tests -- not
    `|X|` of them."""
    DM = sup.DM
    tl = args.time
    front = Front(sup)

    def al(mask, lb=0):
        return max_clique_bits(sup.cadj, mask, time_limit=tl, lb=lb)

    full = sup.full
    a0, w0, _ = al(full)
    front.offer(full, sup.total, a0)
    log(f"  start: |X| = {sup.n}, mu = {sup.f(sup.total):.6f}, alpha = {a0}, "
        f"excess = {sup.f(sup.total - a0 * DM):+.6f}")

    # (a) ascent: remove the lightest vertex whose deletion drops alpha (excess gains 1 - w[v])
    cur, curA, curM, wit = full, a0, sup.total, w0
    steps, trace = 0, []
    while True:
        cand = [(sup.w[v], v) for v in bits_of(wit) if sup.w[v] < DM and drops(sup, cur, v, curA, tl)]
        if not cand:
            break
        wv, v = min(cand)
        cur &= ~(1 << v)
        curM -= wv
        curA, wit, _ = al(cur)
        front.offer(cur, curM, curA)
        steps += 1
        trace.append((v, sup.f(wv), curA, sup.f(curM - curA * DM)))
    asc = {'size': popcount(cur), 'mass': sup.f(curM), 'alpha': curA,
           'excess': sup.f(curM - curA * DM), 'members': list(bits_of(cur)), 'steps': steps,
           'trace': trace}
    log(f"  ascent (remove the lightest vertex whose deletion drops alpha): "
        f"|X| = {asc['size']}, mu = {asc['mass']:.6f}, alpha = {asc['alpha']}, "
        f"excess = {asc['excess']:+.6f} after {steps} removals")
    for (v, wv, a, ex) in trace:
        log(f"      -{v} (mass {wv:.4f}) -> alpha {a}, excess {ex:+.6f}")

    # (b) minimalisation: keep deleting while the excess stays > 0, to a minimal violated set
    rng = random.Random(args.seed)
    mins, seen = [], set()
    starts = [(full, True)] + [(full, False)] * args.restarts
    for st, greedy in starts:
        X, M, A = minimalise(sup, st, rng, greedy, tl, front)
        if X in seen:
            continue
        seen.add(X)
        mins.append({'size': popcount(X), 'mass': sup.f(M), 'alpha': A,
                     'excess': sup.f(M - A * DM), 'members': list(bits_of(X))})
    mins.sort(key=lambda d: (d['size'], -d['excess']))
    log(f"  minimal violated rank subsets (no single vertex removable without killing the "
        f"violation), {len(mins)} distinct from {len(starts)} runs:")
    for d in mins[:6]:
        log(f"    |X| = {d['size']:3d}  mu = {d['mass']:.6f}  alpha = {d['alpha']}  "
            f"excess = {d['excess']:+.6f}")
    smallest = min(d['size'] for d in mins)
    sm = min(mins, key=lambda d: (d['size'], -d['excess']))
    log(f"  SMALLEST violated rank subset found: |X| = {smallest}, alpha = {sm['alpha']}, "
        f"mu = {sm['mass']:.6f}, excess = {sm['excess']:+.6f} "
        f"(the excess floor is one mass quantum)")

    # (c) the alpha-level sweep: heaviest X with alpha(G[X]) <= k, for every k.  A violated
    # rank inequality of rank k exists iff this exceeds k.
    sweep = alpha_sweep(sup, args, tl, rng, front)

    # (d) exhaustive check of the small end
    t0 = time.time()
    small = exhaustive_small(sup, args.small, tl)
    log(f"  exhaustive: X with |X| <= {args.small}, alpha(G[X]) <= 2 and mu > 2: "
        f"{len(small)} ({time.time() - t0:.1f} s)")
    for a in small[:5]:
        log(fmt_struct(sup, a))

    # (e) the answer: how small can a witness be and still carry excess?
    rows = front.rows()
    gap = sup.total - a0 * DM
    log(f"  size-vs-excess front (largest mu(X) - alpha(G[X]) found at each |X|):")
    for r in rows:
        log(f"    |X| = {r['size']:3d}  alpha = {r['alpha']:2d}  mu = {r['mass']:9.6f}  "
            f"excess = {r['excess']:+.6f}  ({100 * r['excess'] / sup.f(gap):.0f} % of the gap)")
    marks = {}
    for label, want in (('0.1', int(0.1 * DM)), ('half', gap // 2), ('all', gap)):
        k, _ = front.smallest_with(want)
        marks[label] = k
        log(f"    smallest |X| with excess >= {sup.f(want):.6f}: {k}")
    # independent re-verification of every reported set: mass resummed, alpha recomputed by a
    # fresh unseeded complete B&B
    bad = 0
    for r in rows:
        mask = sum(1 << v for v in r['members'])
        a, _, comp = al(mask)
        if not comp or a != r['alpha'] or abs(sup.mass(mask) / DM - r['mass']) > 1e-12:
            bad += 1
    log(f"    re-verified {len(rows)} front sets (mass resummed, alpha by a fresh complete "
        f"B&B): {len(rows) - bad} ok, {bad} MISMATCH")
    res['front_verified'] = (bad == 0)
    # the heaviest alpha = 2 set this search found (a LOWER bound on what a rank-2 anchor
    # family could capture; step 4's pentagon search usually beats it)
    r2 = max((v for v in front.best.values() if v[2] == 2), default=None,
             key=lambda v: v[1])
    if r2:
        log(f"    heaviest alpha = 2 subset found by this search: mu = {sup.f(r2[1]):.6f} on "
            f"|X| = {popcount(r2[3])} (excess {sup.f(r2[0]):+.6f}) -- a heuristic lower bound "
            f"on what a rank-2 family can capture, not a ceiling")
        res['rank2_best'] = {'mass': sup.f(r2[1]), 'excess': sup.f(r2[0]),
                                'size': popcount(r2[3]), 'members': list(bits_of(r2[3]))}
    best = max(front.best.items(), key=lambda kv: kv[1][0])
    log(f"  BEST excess anywhere: {sup.f(best[1][0]):+.6f} on |X| = {best[0]} "
        f"(alpha = {best[1][2]}, mu = {sup.f(best[1][1]):.6f})")
    ba = anatomy(sup, list(bits_of(best[1][3])), f"best rank-{best[1][2]}", best[1][2] * DM)
    log(fmt_struct(sup, ba))

    res['ascent'] = asc
    res['minimal'] = mins[:10]
    res['minimal_smallest'] = smallest
    res['alpha_sweep'] = sweep
    res['small_exhaustive'] = [strip(a) for a in small[:10]]
    res['front'] = rows
    res['front_marks'] = marks
    res['best_excess'] = {'excess': sup.f(best[1][0]), 'size': best[0], 'alpha': best[1][2],
                          'mass': sup.f(best[1][1]), 'anatomy': strip(ba)}
    return mins


def drops(sup, mask, v, curA, tl):
    """does removing v lower alpha(G[mask])?  (seeded B&B: only looks for alpha-sized sets)"""
    a, _, _ = max_clique_bits(sup.cadj, mask & ~(1 << v), time_limit=tl, lb=curA - 1)
    return a < curA


def minimalise(sup, start, rng, greedy, tl, front=None):
    """delete vertices from `start` while `mu(X) > alpha(G[X])` survives, until no single
    vertex can be removed: a MINIMAL violated rank subset.  Only the members of a maximum
    independent set can make alpha fall, so one alpha call per round plus `alpha` tests."""
    DM = sup.DM
    cur = start
    curA, wit, _ = max_clique_bits(sup.cadj, cur, time_limit=tl)
    curM = sup.mass(cur)
    while True:
        if front is not None:
            front.offer(cur, curM, curA)
        ex = curM - curA * DM
        cand = []
        for v in bits_of(cur):
            if sup.w[v] < ex:                       # alpha cannot rise; excess stays > 0
                cand.append((v, curM - sup.w[v], curA))
            elif (1 << v) & wit and drops(sup, cur, v, curA, tl):
                if curM - sup.w[v] - (curA - 1) * DM > 0:
                    cand.append((v, curM - sup.w[v], curA - 1))
        if not cand:
            break
        if greedy:                                  # leave the largest excess
            v, nm, na = max(cand, key=lambda c: c[1] - c[2] * DM)
        else:
            v, nm, na = cand[rng.randrange(len(cand))]
        cur &= ~(1 << v)
        curM = nm
        curA, wit, _ = max_clique_bits(sup.cadj, cur, time_limit=tl)
    return cur, curM, curA


def alpha_sweep(sup, args, tl, rng, front=None):
    """for k = 1..alpha(G): a heuristic maximum of mu(X) over X with alpha(G[X]) <= k.
    Greedy growth (heaviest first) with randomised restarts and a re-add pass; each level is
    also seeded with the best sets of the levels below (which satisfy alpha <= k too), making
    the sweep monotone in k and much less noisy.

    `mu > k` means a violated rank inequality of rank k; `mu <= k` for every k below alpha(G)
    would mean the whole gap sits in the one full-support inequality.  Each best set has its
    alpha recomputed exactly by an unseeded complete B&B and is then minimalised, so the
    reported inequality is the one actually certified."""
    DM = sup.DM
    n = sup.n
    order = sorted(range(n), key=lambda v: (-sup.w[v], v))
    out = {}
    a0, _, _ = max_clique_bits(sup.cadj, sup.full, time_limit=tl)
    log(f"  alpha-level sweep (heuristic maxima of mu(X) subject to alpha(G[X]) <= k; "
        f"{args.sweep + 1} restarts):")

    def ok_add(X, v, k):
        """alpha(X + v) <= k, given alpha(X) <= k: any bigger independent set must use v"""
        if k == 0:
            return False
        sub = X & sup.cadj[v]
        a, _, _ = max_clique_bits(sup.cadj, sub, time_limit=tl, lb=k - 1)
        return a <= k - 1

    def grow(seed, seq, k):
        X, M = seed, sup.mass(seed)
        for v in seq:
            if not (X >> v) & 1 and ok_add(X, v, k):
                X |= 1 << v
                M += sup.w[v]
        improved = True
        while improved:                                        # add anything still addable
            improved = False
            for v in order:
                if not (X >> v) & 1 and ok_add(X, v, k):
                    X |= 1 << v
                    M += sup.w[v]
                    improved = True
        return M, X

    seeds = []
    for k in range(1, a0 + 1):
        best = (0, 0)
        seqs = [order]
        for _ in range(args.sweep):
            seq = order[:]
            for i in range(len(seq) - 1, 0, -1):               # mass-biased shuffle
                j = rng.randrange(i + 1) if rng.random() < 0.35 else i
                seq[i], seq[j] = seq[j], seq[i]
            seqs.append(seq)
        for sd in [0] + seeds:                                 # monotone in k
            for seq in seqs:
                M, X = grow(sd, seq, k)
                if M > best[0]:
                    best = (M, X)
        M, X = best
        seeds = [X] + seeds[:2]
        ak, _, comp = max_clique_bits(sup.cadj, X, time_limit=tl)      # exact, unseeded
        if front is not None:
            front.offer(X, M, ak)
        rec = {'mass': sup.f(M), 'size': popcount(X), 'alpha_exact': ak, 'complete': comp,
               'excess': sup.f(M - ak * DM), 'members': list(bits_of(X))}
        flag = 'VIOLATED' if M > ak * DM else ''
        log(f"    k = {k:2d}: mu <= {sup.f(M):9.6f} on |X| = {popcount(X):3d}, "
            f"alpha = {ak}  (excess {sup.f(M - ak * DM):+.6f}) {flag}")
        if M > ak * DM:                                        # shrink to a minimal witness
            Y, MY, AY = minimalise(sup, X, rng, True, tl, front)
            rec['minimal'] = {'size': popcount(Y), 'mass': sup.f(MY), 'alpha': AY,
                              'excess': sup.f(MY - AY * DM), 'members': list(bits_of(Y))}
            a = anatomy(sup, list(bits_of(X)), f'rank-{ak} |X|={popcount(X)}', ak * DM)
            rec['anatomy'] = strip(a)
            log(f"         minimal witness: |X| = {popcount(Y)}, mu = {sup.f(MY):.6f}, "
                f"alpha = {AY}, excess = {sup.f(MY - AY * DM):+.6f}")
            if popcount(X) <= args.show:
                log(fmt_struct(sup, a))
        out[k] = rec
    return out


def exhaustive_small(sup, kmax, tl):
    """every X with alpha(G[X]) = 2 and mu(X) > 2, of size <= kmax, found exhaustively over
    the heavy vertices: a set with alpha <= 2 has a complement graph with no triangle, so we
    grow sets keeping the complement triangle-free.  (alpha(G[X]) = 1 is a clique -- already
    handled by the QSTAB constraints, and every measure here has max clique = 1.)"""
    DM = sup.DM
    n = sup.n
    out = []
    order = sorted(range(n), key=lambda v: (-sup.w[v], v))
    rank = [0] * n
    for i, v in enumerate(order):
        rank[v] = i
    seen = set()

    def rec(cur, curmask, curw, start):
        if len(cur) >= 3 and curw > 2 * DM:
            key = curmask
            if key not in seen:
                seen.add(key)
                a, _, _ = max_clique_bits(sup.cadj, curmask, time_limit=tl)
                if a <= 2:
                    out.append(anatomy(sup, sorted(cur), f'rank|X|={len(cur)}', 2 * DM))
        if len(cur) >= kmax:
            return
        rem = kmax - len(cur)
        for i in range(start, n):
            v = order[i]
            if curw + rem * sup.w[v] <= 2 * DM:
                break                                   # mass-descending: nothing left
            # keep alpha <= 2: no independent triple.  v must not form an independent triple
            # with two current members.
            nonadj = sup.cadj[v] & curmask
            bad = False
            for u in bits_of(nonadj):
                if sup.cadj[u] & nonadj:
                    bad = True
                    break
            if bad:
                continue
            rec(cur + [v], curmask | (1 << v), curw + sup.w[v], i + 1)

    rec([], 0, 0, 0)
    out.sort(key=lambda a: -a['excess'])
    return out


def cmd_anchors(sup, args, res):
    """step 4: the certifiable anchor-graph version of a rank-2 inequality.

    Take five points `A_0..A_4` and the five closed segments `E_i = [A_i, A_{i+1}]` of the
    pentagon they span.  Let `X_i = {S : E_i in S}`.  For `i ~ i+1` any `S in X_i` and
    `S' in X_{i+1}` both contain `A_{i+1}`, so they meet -- the guaranteed-meet graph of the
    five pieces contains `C5`, whose independence number is 2, and Lemma 0 of
    `notes/clique-family.md` gives

        mu(X_0 u ... u X_4) <= alpha(C5) = 2

    for EVERY fractional packing measure, with no geometry beyond "a square is convex, so it
    contains a segment iff it contains both endpoints".  This is the cheapest certifiable rank
    inequality that is not a clique inequality.  The question the brief asks is whether any
    such pentagon carries mass > 2 on these measures: that is the separation problem for the
    family, solved here by hill-climbing on the five anchors over the exact arrangement
    vertices of the support (plus the square centres and the grid points `{1,2,3}^2`)."""
    t0 = time.time()
    V = lc.enumerate_vertices(sup.squares, sup.t, verbose=False)
    Inc = lc.incidences(sup.squares, V, verbose=False)
    cm = [0] * len(V)                                  # vertex -> mask of squares containing it
    for si, col in enumerate(Inc.cols):
        b = 1 << si
        for vi in col:
            cm[int(vi)] |= b
    uniq = {}
    for vi, m in enumerate(cm):
        if m and (m not in uniq or sup.mass(m) > sup.mass(uniq[m])):
            uniq[m] = vi
    cands = sorted(uniq.items(), key=lambda kv: -sup.mass(kv[0]))[:args.cands]
    masks = [m for m, _ in cands]
    verts = [V[vi] for _, vi in cands]
    log(f"  step 4: {len(V)} arrangement vertices -> {len(uniq)} distinct point cliques, "
        f"top {len(masks)} kept ({time.time() - t0:.1f} s)")

    memo = {}

    def value(idx):
        cov = (masks[idx[0]] & masks[idx[1]]) | (masks[idx[1]] & masks[idx[2]]) | \
              (masks[idx[2]] & masks[idx[3]]) | (masks[idx[3]] & masks[idx[4]]) | \
              (masks[idx[4]] & masks[idx[0]])
        m = memo.get(cov)
        if m is None:
            m = memo[cov] = sup.mass(cov)
        return m, cov

    rng = random.Random(args.seed)
    best = (0, None, 0)
    seeds = [[rng.randrange(len(masks)) for _ in range(5)] for _ in range(args.pent)]
    heavy = list(range(min(8, len(masks))))            # also start from the heaviest points
    for a in heavy:
        for b in heavy:
            seeds.append([a, b, a, b, a])
    for idx in seeds:
        cur, cov = value(idx)
        improved = True
        while improved:
            improved = False
            for pos in range(5):
                old = idx[pos]
                for c in range(len(masks)):
                    idx[pos] = c
                    m, cv = value(idx)
                    if m > cur:
                        cur, cov, old, improved = m, cv, c, True
                idx[pos] = old
        if cur > best[0]:
            best = (cur, list(idx), cov)
    m, idx, cov = best
    A = [verts[i] for i in idx] if idx else []
    pts = [[float(Fr(a[0], a[2])), float(Fr(a[1], a[2]))] for a in A]
    gap = res.get('gap', 0) or 1e-18
    log(f"  best pentagon anchor family: mu = {sup.f(m):.6f} against the rank bound 2 "
        f"(excess {sup.f(m - 2 * sup.DM):+.6f} = {100 * sup.f(m - 2 * sup.DM) / gap:.0f} % of "
        f"the gap mass - alpha) on {popcount(cov)} poses")
    log(f"    anchors A_0..A_4 = {pts}")
    # ---- independent re-derivation, from the exact rational anchors, sharing nothing with the
    # Incidence machinery above: a convex square contains the segment iff it contains both ends
    cov2 = 0
    for si, s in enumerate(sup.squares):
        for i in range(5):
            P, Q = A[i], A[(i + 1) % 5]
            if lc.sq_contains(s, P[0], P[1], P[2]) and lc.sq_contains(s, Q[0], Q[1], Q[2]):
                cov2 |= 1 << si
                break
    ap, _, comp = max_clique_bits(sup.cadj, cov, time_limit=args.time)
    log(f"    re-derived membership from the exact anchors: "
        f"{'MATCHES' if cov2 == cov else 'MISMATCH'}; "
        f"alpha(G[X]) = {ap} (complete {comp}; Lemma 0 predicts <= 2), "
        f"mu = {sup.f(sup.mass(cov)):.6f}")
    if cov:
        log(fmt_struct(sup, anatomy(sup, list(bits_of(cov)), 'pentagon rank-2', 2 * sup.DM)))
    res['pentagon'] = {'mass': sup.f(m), 'excess': sup.f(m - 2 * sup.DM),
                       'share_of_gap': sup.f(m - 2 * sup.DM) / gap,
                       'size': popcount(cov), 'anchors': pts, 'alpha': ap,
                       'anchors_exact': [[str(Fr(a[0], a[2])), str(Fr(a[1], a[2]))] for a in A],
                       'rederived_ok': cov2 == cov, 'members': list(bits_of(cov))}
    return res['pentagon']


# ============================================================================ driver
def run_one(path, args):
    log(f"\n=== {os.path.basename(path)}")
    t0 = time.time()
    sup = Support(path, r=args.r)
    res = {'file': path, 'name': sup.name}
    log(f"  loaded in {time.time() - t0:.1f} s; DM = {sup.DM}")
    cmd_alpha(sup, args, res)
    if args.step >= 2:
        out = cmd_cycles(sup, args, res)
    else:
        out = {}
    if args.step >= 3:
        cmd_search(sup, args, res)
    if args.step >= 4:
        cmd_anchors(sup, args, res)
    return res


def selftest():
    """the 5-cycle of unit squares that pairwise-meet only consecutively: mass 5 * 0.5 = 2.5
    against the C5 rank bound 2 -- checks the hole finder, alpha and the exhaustive small
    search on a synthetic measure."""
    import tempfile
    # five unit squares arranged so that S_i meets S_{i+1} only (a ring), t = 8
    t = Fr(8)
    R = 0.8       # ring radius: consecutive centres 0.94 apart (meet), the others 1.52 (do not)
    lines = [f'# t = {t}', 'sym 1']
    cs = []
    for i in range(5):
        th = 2 * math.pi * i / 5
        cx = Fr(round(4 + R * math.cos(th), 3)).limit_denominator(1000)
        cy = Fr(round(4 + R * math.sin(th), 3)).limit_denominator(1000)
        cs.append((cx, cy))
        lines.append(f'pose 0 1 {cx} {cy} 1/2')
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as fh:
        fh.write('\n'.join(lines) + '\n')
        p = fh.name
    sup = Support(p)
    ok = True

    def chk(name, cond, extra=''):
        nonlocal ok
        log(f"  {'ok  ' if cond else 'FAIL'} {name} {extra}")
        ok = ok and cond

    chk('5 squares', sup.n == 5)
    chk('ring graph: 5 edges', sup.nedge == 5, f'(got {sup.nedge})')
    a, mask, comp = alpha(sup)
    chk('alpha(C5) = 2', a == 2 and comp)
    chk('mass 2.5', abs(sup.f(sup.total) - 2.5) < 1e-12)
    cyc = chordless_odd_cycles(sup.adj, sup.n, 5, sup.w, 2 * sup.DM)
    chk('one violated chordless C5', len(cyc) == 1, f'(got {len(cyc)})')
    if cyc:
        chk('  excess 0.5', abs(sup.f(sum(sup.w[v] for v in cyc[0])) - 2.5) < 1e-12)
    c7 = chordless_odd_cycles(sup.adj, sup.n, 7, sup.w, 3 * sup.DM)
    chk('no C7', len(c7) == 0)
    ah = chordless_odd_cycles(sup.cadj, sup.n, 7, sup.w, 2 * sup.DM)
    chk('no 7-antihole', len(ah) == 0)
    sm = exhaustive_small(sup, 5, 60.0)
    chk('exhaustive small finds the C5', any(len(x['members']) == 5 for x in sm),
        f'(got {[x["size"] for x in sm]})')
    os.unlink(p)

    # a 7-ring: consecutive centres 0.87 apart (meet), the next 1.56 and 1.95 (do not)
    lines = [f'# t = {t}', 'sym 1']
    for i in range(7):
        th = 2 * math.pi * i / 7
        cx = Fr(round(4 + math.cos(th), 3)).limit_denominator(1000)
        cy = Fr(round(4 + math.sin(th), 3)).limit_denominator(1000)
        lines.append(f'pose 0 1 {cx} {cy} 1/2')
    with tempfile.NamedTemporaryFile('w', suffix='.txt', delete=False) as fh:
        fh.write('\n'.join(lines) + '\n')
        p7 = fh.name
    s7 = Support(p7)
    chk('7-ring: 7 edges', s7.nedge == 7, f'(got {s7.nedge})')
    a7, _, _ = alpha(s7)
    chk('alpha(C7) = 3', a7 == 3, f'(got {a7})')
    c7 = chordless_odd_cycles(s7.adj, s7.n, 7, s7.w, 3 * s7.DM)
    chk('one violated chordless C7 (mass 3.5)', len(c7) == 1, f'(got {len(c7)})')
    c5 = chordless_odd_cycles(s7.adj, s7.n, 5, s7.w, 0)
    chk('no chordless C5 in the 7-ring', len(c5) == 0, f'(got {len(c5)})')
    os.unlink(p7)
    log('  selftest ' + ('PASSED' if ok else 'FAILED'))
    return 0 if ok else 1


def recheck(path):
    """re-verify the pentagons of a results JSON from scratch: rebuild each measure from its
    file, take the five anchors as exact rationals, recompute the membership (a convex square
    contains the segment iff it contains both endpoints), resum the mass, and check Lemma 0 by
    brute force -- no triple of the union is pairwise disjoint.  Shares no state with the run
    that produced the JSON."""
    from itertools import combinations
    out = json.load(open(path))
    allok = True
    for r in out:
        pent = r.get('pentagon')
        if not pent:
            continue
        t, sym, poses = lc.read_measure(r['file'])
        squares, DM, _ = lc.build_squares(t, sym, poses)
        pts = []
        for (sx, sy) in pent['anchors_exact']:
            x, y = Fr(sx), Fr(sy)
            D = x.denominator * y.denominator
            pts.append((x.numerator * y.denominator, y.numerator * x.denominator, D))
        members = []
        for si, s in enumerate(squares):
            for i in range(5):
                P, Q = pts[i], pts[(i + 1) % 5]
                if lc.sq_contains(s, *P) and lc.sq_contains(s, *Q):
                    members.append(si)
                    break
        mass = sum(squares[i][8] for i in members)
        meets = {(a, b): lc.sq_meets_sq(squares[a], squares[b])
                 for a, b in combinations(members, 2)}
        bad = sum(1 for a, b, c in combinations(members, 3)
                  if not meets[(a, b)] and not meets[(a, c)] and not meets[(b, c)])
        same = sorted(members) == sorted(pent['members'])
        ok = (same and abs(mass / DM - pent['mass']) < 1e-12 and bad == 0 and mass > 2 * DM)
        allok = allok and ok
        log(f"{r['name']:<12} |X| = {len(members):3d}  mu = {mass / DM:.6f}  bound 2  "
            f"excess {mass / DM - 2:+.6f}  members match: {same}  "
            f"pairwise-disjoint triples: {bad}  -> {'OK' if ok else 'FAIL'}")
    log('ALL PENTAGONS VERIFIED' if allok else 'FAILURES')
    return 0 if allok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='*', help='measure files (leaf_ceiling format)')
    ap.add_argument('--step', type=int, default=3, help='1 alpha, 2 +cycles, 3 +search, 4 +anchors')
    ap.add_argument('--time', type=float, default=300.0, help='B&B time limit per alpha call')
    ap.add_argument('--r', default='1', help='region inset (region_boxes)')
    ap.add_argument('--restarts', type=int, default=12, help='randomised minimalisation restarts')
    ap.add_argument('--small', type=int, default=8, help='exhaustive small-X size cap')
    ap.add_argument('--cands', type=int, default=300, help='step 4: candidate anchor points')
    ap.add_argument('--pent', type=int, default=40, help='step 4: pentagon hill-climb restarts')
    ap.add_argument('--show', type=int, default=40, help='print full anatomy up to this size')
    ap.add_argument('--sweep', type=int, default=6, help='randomised restarts per alpha level')
    ap.add_argument('--seed', type=int, default=20260911)
    ap.add_argument('--json', default=None)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--recheck', default=None,
                    help='re-verify the pentagons of a results JSON from scratch and exit')
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if args.recheck:
        return recheck(args.recheck)
    if not args.files:
        ap.error('give at least one measure file (or --selftest)')
    out = []
    for p in args.files:
        out.append(run_one(p, args))
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(out, fh, indent=1)
        log(f"\nwrote {args.json}")
    log('\n=== summary')
    log(f"{'measure':<12} {'n':>4} {'mass':>10} {'a':>3} {'gap':>9} "
        f"{'C5':>3} {'C7':>3} {'AH7':>4} {'W5':>3} {'|X|>=.1':>8} {'best ex':>9} "
        f"{'pentagon':>9}")
    for r in out:
        cc = r.get('cycle_counts', {})
        fm = r.get('front_marks', {})
        be = r.get('best_excess', {})
        pg = r.get('pentagon', {})
        log(f"{r['name']:<12} {r['n']:>4} {r['mass']:>10.6f} {r['alpha']:>3} {r['gap']:>9.6f} "
            f"{cc.get('C5', '-'):>3} {cc.get('C7', '-'):>3} {cc.get('antihole7', '-'):>4} "
            f"{cc.get('wheel5', '-'):>3} {str(fm.get('0.1', '-')):>8} "
            f"{be.get('excess', float('nan')):>9.6f} "
            f"{pg.get('mass', float('nan')):>9.6f}")
    return 0


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    sys.exit(main())
