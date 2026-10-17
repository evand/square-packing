#!/usr/bin/env python3
"""Basin radius of a certified minimum: when is a float posting provably "the same local minimum"?  (10-17 prototype)

Evan: "a provably correct way to say two packings are the same when they're numerically solved … thresholds derived
from numerical analysis, rather than magic constants".  For a certified minimum M (an exactsolve output directory) and
postings P:

  distance  d(P, M) on M's rigid squares (load-bearing, minus squares that move in M's exact flat modes): squares are
            interchangeable (assignment on centres), angles mod 90, the box's 8 symmetries; per-square metric
            max(|dx|, |dy|, |dθ|/√2), so every corner moves by at most 2d.
  margins   where M's contact structure could change: the smallest gap of non-touching pairs / walls (exactsolve's
            margin_separated, min_wall_other), of corner-on-side incidences of touching pairs that aren't touching, and
            the distance of each touching corner from its side's ends (corner-corner touches are excepted: they are
            disjunctions the certifier branches on).
  growth    φ(t) = min over rigid coordinates c and signs of the side slp2 reaches (trust radius t, from M) with c held
            at displacement t.  Every point on the radius-t sphere (in the metric above) has such a coordinate, so φ(t)
            is the empirical lower bound of S - S* on that sphere.
  test      P is in M's basin if d(P, M) < t and S(P) - S* < φ(t) for some tested t, with φ > 0 at every smaller tested
            radius: a continuous descent from P cannot cross the sphere, so it ends inside the ball.

Not yet a proof: φ comes from a local solver (an upper estimate of each face's minimum), and the test doesn't exclude a
second local minimum inside the ball (that needs a uniqueness radius from Lipschitz bounds on the constraint
derivatives).  The prototype measures whether the scales make such a proof worth building.

  basin.py growth MINDIR [--t 1e-6,1e-5,1e-4,1e-3,1e-2] [--procs 15] [--sample K]   -> MINDIR/basin_growth.json
  basin.py report CASES.yaml                                                        distances, margins, growth, test
"""
import argparse, itertools, json, math, os, random, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
sys.path.insert(0, os.path.join(os.path.dirname(PK), 'exact'))
from fractions import Fraction
from multiprocessing import Pool
import numpy as np


def read_pack(path):
    """(S, [(x, y, theta_rad)]) from an 'ours' file; accepts fractions and long decimals"""
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n, S = int(L[0][0]), float(Fraction(L[0][1]))
    sq = []
    for r in L[1:n + 1]:
        d = float(Fraction(r[2]))
        d = (d + 45.0) % 90.0 - 45.0
        sq.append((float(Fraction(r[0])), float(Fraction(r[1])), math.radians(d)))
    return S, sq


def load_min(d):
    """a certified minimum from an exactsolve output dir: S*, configuration, rigid / flat / free squares, report"""
    stem = 'desc' if os.path.exists(os.path.join(d, 'desc.exact.txt')) else 'in'
    rep = json.load(open(os.path.join(d, stem + '.json')))
    S, sq = read_pack(os.path.join(d, stem + '.exact.txt'))
    so = rep.get('second_order') or {}
    free = set(json.loads(rep['free_squares']) if isinstance(rep.get('free_squares'), str) else rep.get('free_squares') or [])
    flat = set(so.get('flat_squares') or [])
    rigid = [i for i in range(len(sq)) if i not in free and i not in flat]
    cts = json.load(open(os.path.join(d, stem + '.contacts.json')))
    return dict(dir=d, S=S, S_exact=rep['S_exact'], sq=sq, rigid=rigid, free=sorted(free), flat=sorted(flat), rep=rep, cts=cts)


# ------------------------------------------------------------------------------------------------ distance
def _d4(S, sq, k):
    """the k-th of the box's 8 symmetries applied to a configuration"""
    out = []
    for x, y, t in sq:
        if k & 4:
            x, y, t = y, x, -t                          # reflect in the diagonal
        if k & 1:
            x, t = S - x, -t
        if k & 2:
            y, t = S - y, -t
        out.append((x, y, (t + math.pi / 4) % (math.pi / 2) - math.pi / 4))
    return out


def distance(P, M):
    """d(P, M) = max over M's rigid squares of max(|dx|, |dy|, |dθ|/√2) after the best box symmetry and assignment"""
    from scipy.optimize import linear_sum_assignment
    SP, sqP = P
    best = None
    B = np.array([(x, y) for x, y, _ in M['sq']])
    for k in range(8):
        A = _d4(SP, sqP, k)
        Ac = np.array([(x, y) for x, y, _ in A])
        C = np.abs(Ac[:, None, :] - B[None, :, :]).max(axis=2)
        r, c = linear_sum_assignment(C)
        to = dict(zip(c, r))                             # M's square -> P's square
        ds = []
        for i in M['rigid']:
            x, y, t = A[to[i]]
            X, Y, T = M['sq'][i]
            dt = (t - T + math.pi / 4) % (math.pi / 2) - math.pi / 4
            ds.append(max(abs(x - X), abs(y - Y), abs(dt) / math.sqrt(2)))
        dall = float(C[r, c].max())
        cand = (max(ds), k, dall)
        if best is None or cand[0] < best[0]:
            best = cand
    return dict(d=best[0], sym=best[1], d_centres_all=best[2])


# ------------------------------------------------------------------------------------------------ margins
def _local(p, sq):
    x, y, t = sq
    c, s = math.cos(t), math.sin(t)
    dx, dy = p[0] - x, p[1] - y
    return c * dx + s * dy, -s * dx + c * dy


def _corners(sq):
    x, y, t = sq
    c, s = math.cos(t), math.sin(t)
    return [(x + 0.5 * (a * c - b * s), y + 0.5 * (a * s + b * c)) for a, b in ((1, 1), (-1, 1), (-1, -1), (1, -1))]


def margins(M):
    """contact-structure margins of M (see the docstring)"""
    sq = M['sq']
    vv = {tuple(sorted(p)) for p in M['cts'].get('vv') or []}
    touching = set()
    for e in M['cts']['equations']:
        if e[0] == 'C':
            touching.add(tuple(sorted((e[1], e[3]))))
    inact, ends = [], []
    for i, j in touching:
        for a, b in ((i, j), (j, i)):
            for p in _corners(sq[a]):
                u, v = _local(p, sq[b])
                gap = max(abs(u), abs(v)) - 0.5                  # >0: corner outside b
                along = 0.5 - min(abs(u), abs(v))                # distance from the nearest end of that side
                if gap < 1e-12:                                  # a touching corner
                    if (i, j) not in vv:
                        ends.append(along)
                elif along > -gap:                               # a corner facing a side of b (not around its end)
                    inact.append(gap)
    rep = M['rep']
    f = lambda k: float(rep[k]) if rep.get(k) not in (None, 'None') else None
    return dict(gap_nontouching=f('margin_separated'), gap_wall=f('min_wall_other'),
                gap_inactive_incidence=min(inact) if inact else None, end_margin=min(ends) if ends else None,
                n_vv=len(vv), n_rigid=len(M['rigid']), n_flat=len(M['flat']), n_free=len(M['free']),
                eig_min_positive=(rep.get('second_order') or {}).get('eig_min_positive'))


# ------------------------------------------------------------------------------------------------ growth
def _step_lp(s, sq, R, c, delta, rounds=3):
    """one SLP step of  min dS  s.t. linearised contacts (inc's model), |z| <= R, z_c = delta, dS free.  Corner-corner
    pairs are disjunctions: as inc.shrink, solve with them dropped, put each on the line that solution violates least,
    re-solve, repeat; keep the best.  Returns (ds, z) or (None, None)."""
    import inc
    from scipy.optimize import linprog
    from scipy.sparse import vstack as svstack
    nz = 3 * len(sq) + 1
    rows, disj = inc.model_fast(s, sq, 3 * R)
    A0, g0 = inc._mat(rows, nz)
    obj = np.zeros(nz); obj[-1] = 1.0
    bnd = [(-1.0, 1.0)] * nz
    bnd[c] = (delta / R, delta / R)
    bnd[-1] = (None, None)                           # the side responds freely (a push may need dS > R)

    def lp(choice):
        blocks, gs = [A0], [g0]
        for (_, alts), q in zip(disj, choice or ()):
            Aq, gq = inc._mat(alts[q], nz)
            blocks.append(Aq); gs.append(gq)
        A = svstack(blocks).tocsr(); g = np.concatenate(gs)
        res = linprog(obj, A_ub=-A, b_ub=np.maximum(g, 0.0) / R, bounds=bnd, method='highs',
                      options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10, 'time_limit': 60})
        return (res.fun * R, res.x * R) if res.status == 0 else (None, None)

    if not disj:
        return lp(None)
    choice = [max(range(len(alts)), key=lambda q: min(alts[q][0]) if alts[q][0] else 1e9) for _, alts in disj]
    best = lp(choice)
    _, z = lp([])                                    # relaxation: corner-corner pairs dropped
    for _ in range(rounds):
        if z is None:
            break
        ch = [max(range(len(alts)), key=lambda q: inc._alt_slack(alts[q], z)) for _, alts in disj]
        if ch == choice:
            break
        choice = ch
        r = lp(choice)
        if r[0] is None:
            break
        if best[0] is None or r[0] < best[0] - 1e-18:
            best = r
        z = r[1]
    return best


def push_min(s, sq, c, target, R0, iters=600, budget=300):
    """Push coordinate c continuously by `target` (trust-region steps; every other coordinate and the side respond),
    then minimise the side with c held.  Overlaps from each step's curvature are removed by inc.correct at fixed side
    with c held (never by rescaling, which would move c).  Returns (s, sq, max overlap left, steps)."""
    import inc, slp2
    t0 = time.time()
    R, moved, it = R0, 0.0, 0
    sq = list(sq)
    while it < iters and time.time() - t0 < budget and R > 1e-13:
        it += 1
        rem = target - moved
        delta = math.copysign(min(abs(rem), R), rem) if abs(rem) > 1e-15 else 0.0
        ds, z = _step_lp(s, sq, R, c, delta)
        if z is None:
            R /= 3; continue
        new = [(x + z[3 * i], y + z[3 * i + 1], t + z[3 * i + 2]) for i, (x, y, t) in enumerate(sq)]
        s2 = s + ds
        if slp2.worst(s2, new) <= 1e-13:
            fix = new                                # already feasible (correct's LP is then degenerate)
        else:
            fix = inc.correct(s2, new, R, fixed=(c,))
            if fix is None or slp2.worst(s2, fix) > 1e-13:
                R /= 3; continue                     # e.g. the step ran into a contact the model didn't include
        if delta == 0.0 and ds > -1e-15:
            R /= 3                                   # held phase: no progress at this radius
            if ds > -1e-15 and R < 1e-11:
                break
            continue
        s, sq, moved = s2, fix, moved + delta
        if delta == 0.0:
            R = min(R * 2, R0)
    return s, sq, slp2.worst(s, sq), it, moved


def _face(args):
    d, c, sign, t = args
    M = load_min(d)
    i, k = divmod(c, 3)
    target = sign * (t * math.sqrt(2) if k == 2 else t)
    t0 = time.time()
    try:
        s2, sq2, w, it, moved = push_min(M['S'], M['sq'], c, target, R0=t / 4)
    except Exception as e:
        return dict(c=c, sign=sign, t=t, err=repr(e)[:80])
    held = abs(sq2[i][k] - M['sq'][i][k]) / (math.sqrt(2) if k == 2 else 1.0)
    return dict(c=c, sign=sign, t=t, dS=s2 - M['S'], held=held, overlap=w, steps=it, reached=abs(moved - target) < 1e-15,
                sec=time.time() - t0)


def growth(a):
    M = load_min(a.mindir)
    T = [float(x) for x in a.t.split(',')]
    coords = [3 * i + k for i in M['rigid'] for k in range(3)]
    if a.sample and a.sample < len(coords):
        coords = sorted(random.Random(0).sample(coords, a.sample))
    jobs = [(a.mindir, c, sg, t) for t in T for c in coords for sg in (1, -1)]
    print(f'{a.mindir}: n={len(M["sq"])}, {len(M["rigid"])} rigid squares, {len(coords)} coordinates x 2 signs x {len(T)} radii '
          f'= {len(jobs)} slp2 runs', flush=True)
    with Pool(a.procs) as pool:
        R = list(pool.imap_unordered(_face, jobs, chunksize=1))
    out = dict(S_exact=M['S_exact'], radii=T, runs=R, phi={})
    for t in T:
        v = [r for r in R if r['t'] == t and 'dS' in r]
        if v:
            w = min(v, key=lambda r: r['dS'])
            out['phi'][repr(t)] = dict(phi=w['dS'], at=w['c'], sign=w['sign'], held=w['held'], n=len(v),
                                       median=float(np.median([r['dS'] for r in v])))
            print(f'  t={t:.0e}: phi={w["dS"]:.3e} (coord {w["c"]}, held {w["held"]:.2e}), median {out["phi"][repr(t)]["median"]:.3e}, {len(v)} runs')
    json.dump(out, open(os.path.join(a.mindir, 'basin_growth.json'), 'w'), indent=1)


def report(a):
    import yaml
    C = yaml.safe_load(open(a.cases))
    for case in C['cases']:
        print(f"== {case['name']}")
        mins = {k: load_min(v) for k, v in case['minima'].items()}
        for k, M in mins.items():
            g = os.path.join(M['dir'], 'basin_growth.json')
            phi = json.load(open(g))['phi'] if os.path.exists(g) else {}
            print(f"  minimum {k}: S* = {M['S_exact'][:22]}  margins {json.dumps(margins(M))}")
            if phi:
                print('    phi(t): ' + ', '.join(f"{float(t):.0e}: {v['phi']:.2e}" for t, v in phi.items()))
            M['phi'] = phi
        for pk_, path in case['postings'].items():
            P = read_pack(path)
            for k, M in mins.items():
                D = distance(P, M)
                dS = P[0] - M['S']
                ok = None
                ts = sorted((float(t), v['phi']) for t, v in M['phi'].items())
                for idx, (t, ph) in enumerate(ts):
                    if any(p <= 0 for _, p in ts[:idx + 1]):
                        break
                    if D['d'] < t and dS < ph:
                        ok = t
                        break
                verdict = (f'IN basin (t = {ok:.0e})' if ok else ('not shown' if ts else 'no growth data'))
                print(f"  posting {pk_:<14} vs {k:<10} d = {D['d']:.3e} (sym {D['sym']}, all centres {D['d_centres_all']:.2e}), "
                      f"S(P) - S* = {dS:+.3e}  -> {verdict}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    g = sp.add_parser('growth'); g.add_argument('mindir'); g.add_argument('--t', default='1e-6,1e-5,1e-4,1e-3,1e-2')
    g.add_argument('--procs', type=int, default=15); g.add_argument('--sample', type=int)
    r = sp.add_parser('report'); r.add_argument('cases')
    a = ap.parse_args()
    {'growth': growth, 'report': report}[a.cmd](a)
