#!/usr/bin/env python3
"""Composer v2: two pieces in opposite corners + the exact MILP upright fill, at a target side (10-12 night).

Why: the formula composer (mondrian.Composer, guillotine columns) undercounts; #528's s(175) = two copies of a
38-square (h + 2) x h piece + 99 uprights is not guillotine-separable (170 by formula, 175 by exact fill).

For each target n: S = target side (default frontier_now(n) - 1e-9, i.e. a strict improvement).  Pieces come from the
store (`python3 -m pk.pieces import ...`).  A (m, d) cell is in when one of its pieces has a side D with S - D within
TIGHT above an integer (the uprights then fill whole rows / columns); all its pieces that fit are used.  Arms =
(pair, variant): A at LL in 4 variants (transpose x mirror-x), B at UR in 4; area screen mA + mB + floor(S^2 - WA HA -
WB HB) >= n; tilted squares of A and B must not overlap; then exact_fill re-fills the uprights.  Arms are pulled within
a time budget per n: pk.bandit.CellBandit over (cell A, cell B) groups (reward 2^-(n - count)), prior weight within a
group.  Every fill is logged to fills.jsonl (to calibrate the prior).  count >= n = success; the layout is fq-polished.

  python3 compose2.py run --ns 100-324 --out runs/compose2/c1 [--target now|cutoff:ISO] [--budget 120] [--procs 14]
"""
from __future__ import annotations
import argparse, json, math, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FQ = os.path.join(HERE, 'target', 'release', 'fq')
PIECES = None
PRIOR = (1.0, 0.01)          # weight = (1 + area slack)^PRIOR[0] * exp(-(gap_A + gap_B) / PRIOR[1])


def load_pieces(mmin, per_cell=10 ** 9, source=None, cf_eps=0.0):
    """Pieces from the store (pk.pieces): those that beat the grid in their own box, one per (cell, h to 1e-7) (copies
    and rattler moves of one minimum), up to per_cell per (m, d), best h first; gap = h / best h in the cell - 1."""
    from pk.store import Store
    keep, seen, cnt = [], set(), {}
    for p in Store().pieces(mmin=mmin, source=source):
        m, d, H = p.n, p.dx, p.s
        W = H + d
        if m <= math.floor(W + 1e-9) * math.floor(H + 1e-9):
            continue                                    # does not beat the grid in its own box
        k = (m, d, round(H, 7))
        cnt_k = cnt.get((m, d), 0)
        if k in seen or cnt_k >= per_cell:
            continue
        seen.add(k); cnt[(m, d)] = cnt_k + 1
        sq = p.sq.copy()
        if cf_eps and d > 0:                            # counterfactual: the same piece with h smaller by eps h
            #                                             (rectangle pieces only: a better square piece = a better s(m))
            sq[:, :2] *= 1 - cf_eps; W *= 1 - cf_eps; H *= 1 - cf_eps
        keep.append(dict(m=m, W=W, H=H, sq=sq, cell=(m, d), label=f"{p.meta['source']}:{m}:{d:g}#{p.meta['id']}"))
    hb = {}
    for q in keep:
        hb[q['cell']] = min(hb.get(q['cell'], math.inf), q['H'])
    for q in keep:
        q['gap'] = q['H'] / hb[q['cell']] - 1
    return keep


def variants(q):
    """(squares, w, h, tag) for transpose x mirror-x (mirror-y is covered by the 180-degree turn of the UR piece)."""
    out = []
    for tr in (0, 1):
        P, w, h = q['sq'], q['W'], q['H']
        if tr:
            P = np.c_[P[:, 1], P[:, 0], (-P[:, 2]) % 90.0]; w, h = h, w
        for mx in (0, 1):
            R = P.copy()
            if mx:
                R[:, 0] = w - R[:, 0]; R[:, 2] = (-R[:, 2]) % 90.0
            out.append((R, w, h, (tr, mx)))
    return out


NONE = dict(m=0, W=0.0, H=0.0, sq=np.zeros((0, 3)), cell=(0, 0), label='-', gap=0.0)


def place(A, B, v, S):
    """Tilted squares of A (variant v // 4) at LL and B (variant v % 4, turned 180 degrees) at UR of the S-box; None if
    they overlap."""
    import mondrian
    if max(A['W'], A['H'], B['W'], B['H']) > S + 1e-12:   # a box sticking out of the S-box (fill does not check walls)
        return None
    RA = variants(A)[v // 4][0]
    RB = variants(B)[v % 4][0]
    TA = RA[mondrian.is_tilted(RA[:, 2])]
    TB = RB.copy(); TB[:, 0] = S - TB[:, 0]; TB[:, 1] = S - TB[:, 1]
    TB = TB[mondrian.is_tilted(TB[:, 2])]
    return None if len(TB) and mondrian._sat_overlap(TA, TB) else np.r_[TA, TB]


def tight(D, S, tol):
    r = S - D
    return r >= -1e-9 and (r - math.floor(r + 1e-9)) < tol


def job(args):
    n, S, tol, tl, budget, max_fills, seed, probe, probe_steps = args
    S = S + probe                                    # probe: fill at the relaxed side, bisect each fitting arm
    near = []
    import mondrian
    P = PIECES
    cell = {}
    for i, q in enumerate(P):
        cell.setdefault(q['cell'], []).append(i)
    # A cell is in if any of its pieces is tight; then every piece of it that fits is tried (a better piece leaves
    # more slack than a tight worse one, which only makes the screen non-monotone, not the piece useless).
    cand = [i for c, I in cell.items() if any(tight(P[i]['W'], S, tol + probe) or tight(P[i]['H'], S, tol + probe) for i in I)
            for i in I if max(P[i]['W'], P[i]['H']) <= S + 1e-9 and P[i]['m'] < n]
    P = P + [NONE]                                      # B = NONE: one piece + uprights (e.g. s(327) = s(123) + 7)
    pairs, w = [], []
    for ia in cand:
        for ib in cand + [len(P) - 1]:
            A, B = P[ia], P[ib]
            sl = A['m'] + B['m'] + math.floor(S * S - A['W'] * A['H'] - B['W'] * B['H'] + 1e-9) - n
            if sl >= 0:
                pairs.append((ia, ib, sl))
                # prior weight (uncalibrated; fit it from fills.jsonl): more spare area, closer to the cell's best h
                w.append((1 + sl) ** PRIOR[0] * math.exp(-(A['gap'] + B['gap']) / PRIOR[1]))
    # Arms: groups = (cell of A, cell of B), chosen by pk.bandit.CellBandit (Thompson, reward per fill-second; the
    # reward is graded, 2^-(n - count), since hits are too rare to learn from); within a group the next untried
    # (pair, variant) is drawn by the prior weight.  Budget: seconds per n.
    import random
    from pk.bandit import CellBandit
    rng, nrng = random.Random(seed), np.random.default_rng(seed)
    w = np.array(w, float)
    groups = {}
    for k, (ia, ib, sl) in enumerate(pairs):
        groups.setdefault((P[ia]['cell'], P[ib]['cell']), []).append(k)
    left = {k: list(range(16)) if P[pairs[k][1]]['m'] else [0, 4, 8, 12] for k in range(len(pairs))}  # untried variants
    bandit = CellBandit(b0=5.0)                                   # fills take ~0.1-10 s, not explorer minutes
    best, tried, fills, t0 = None, 0, [], time.time()
    while groups and time.time() - t0 < budget and tried < max_fills:
        g = bandit.choose(sorted(groups), rng)
        ks = np.array(groups[g])
        k = int(nrng.choice(ks, p=w[ks] / w[ks].sum()))
        v = left[k].pop(int(nrng.integers(len(left[k]))))
        if not left[k]:
            del left[k]; groups[g].remove(k)
            if not groups[g]:
                del groups[g]
        ia, ib, sl = pairs[k]
        A, B = P[ia], P[ib]
        T = place(A, B, v, S)
        if T is None:
            continue
        tried += 1
        t1 = time.time()
        L = mondrian.exact_fill(T, S, time_limit=tl)
        c = len(L)
        bandit.update((A['cell'], B['cell']), time.time() - t1, 2.0 ** -(n - c) if c < n else 1.0)
        rec = dict(n=n, A=A['label'], B=B['label'], cA=A['cell'], cB=B['cell'], hA=A['H'], hB=B['H'], v=v, sl=sl, gA=A['gap'], gB=B['gap'], c=c,
                   sec=round(time.time() - t1, 2))
        if probe and c >= n:                    # near-miss probe: smallest side at which this arm still fits n
            lo, hi, Lh = S - 2 * probe, S, L        # bisect over [frontier - probe, frontier + probe]
            for _ in range(probe_steps):
                mid = (lo + hi) / 2
                Tm = place(A, B, v, mid)
                Lm = mondrian.exact_fill(Tm, mid, time_limit=tl) if Tm is not None else []
                if len(Lm) >= n:
                    hi, Lh = mid, Lm
                else:
                    lo = mid
            rec.update(S_min=hi, need=hi - (S - probe))     # need = S_min - frontier target (> 0: miss by that much)
            near.append(dict(rec, sq=np.asarray(Lh).tolist()))
        fills.append(rec)
        if best is None or c > best[0]:
            best = (c, A['label'], (v // 4), B['label'], (v % 4), L)
        if c >= n and not probe:
            break
    near.sort(key=lambda r: r['S_min'])
    base = dict(n=n, S=S - probe, near=near[:20], cand=len(cand), pairs=len(pairs), arms=16 * len(pairs), tried=tried,
                sec=round(time.time() - t0, 1), fills=fills)
    if best is None:
        return base
    c, la, ta, lb, tb, L = best
    if probe and near and near[0]['need'] < -1e-9:          # an arm fits below the frontier: a hit; polish its layout
        x = near[0]
        return dict(base, count=n, A=x['A'], B=x['B'], v=x['v'], S=x['S_min'], sq=x['sq'])
    return dict(base, count=c, A=la, ta=ta, B=lb, tb=tb, sq=L.tolist() if c >= n and not probe else None)


def polish(n, S, sq, out, tag=''):
    """Composition layout -> trimmed to n at side S -> fq polish (all squares free: the seam relaxes).  {s_pol | err}."""
    import mondrian
    X = mondrian.trim(np.array(sq, float), n)
    raw = os.path.join(out, f'raw-{n}{tag}.txt'); outp = os.path.join(out, f'n-{n}{tag}.txt')
    with open(raw, 'w') as g:
        g.write(f'{n} {float(S)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in X))
    q = subprocess.run([FQ, 'quench', '--in', raw, '--out', outp, '--loosen', '1.0', '--no-alm',
                        '--flip-top', '64', '--pit', '600'], capture_output=True, text=True)
    try:
        return dict(s_pol=json.loads(q.stdout.strip().splitlines()[-1])['s'], out=outp)
    except Exception:
        return dict(err=q.stderr[-200:])


def _polish_job(t):
    return t[0], t[1], polish(*t[2:])


def polish_near(a):
    """Seam polish of the probe's tightest arms (need <= MAX_NEED: ties within the fill tolerance, or near misses):
    the same side as the frontier from different pieces is a different packing, and freeing every square may go lower.
    Success = s_pol < frontier_now - 1e-9, then candidate.py."""
    from pk import eventbench as E
    R = [json.loads(l) for d in a.dirs for l in open(os.path.join(d, 'near.jsonl'))]
    by = {}
    for r in sorted(R, key=lambda r: r['need']):
        if r['need'] <= a.max_need:
            L = by.setdefault(r['n'], [])
            if len(L) < a.per_n and all((r['A'], r['B']) != (x['A'], x['B']) for x in L):
                L.append(r)
    os.makedirs(a.out, exist_ok=True)
    tasks = [(n, r, n, r['S_min'], r['sq'], a.out, f'_{i}') for n, L in sorted(by.items()) for i, r in enumerate(L)]
    print(f'{len(tasks)} polishes over {len(by)} n', flush=True)
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'polish.jsonl'), 'w') as f:
        for n, r, res in ex.map(_polish_job, tasks):
            fr = E.frontier_now(n)
            gain = fr - res.get('s_pol', math.inf)
            f.write(json.dumps(dict(n=n, A=r['A'], B=r['B'], need=r['need'], frontier=fr, gain=gain, **res)) + '\n')
            print(f"{n:4d} {r['A'][:28]:28s} {r['B'][:28]:28s} need {r['need']:.1e} s_pol {res.get('s_pol', math.nan):.10f}"
                  f" gain {gain:.2e}{'  BELOW' if gain > 1e-9 else ''}", flush=True)


def place3(A, B, C, va, vb, vc, S):
    """Three corners: A at LL (variant va), B at UR (vb, turned 180 degrees), C at UL (vc, mirrored in y).  Tilted squares
    only; None if any box sticks out or two pieces' tilted squares overlap."""
    import mondrian
    if max(q[k] for q in (A, B, C) for k in ('W', 'H')) > S + 1e-12:
        return None
    RA, RB, RC = variants(A)[va][0], variants(B)[vb][0], variants(C)[vc][0]
    TA = RA[mondrian.is_tilted(RA[:, 2])]
    TB = RB.copy(); TB[:, 0] = S - TB[:, 0]; TB[:, 1] = S - TB[:, 1]; TB = TB[mondrian.is_tilted(TB[:, 2])]
    TC = RC.copy(); TC[:, 1] = S - TC[:, 1]; TC[:, 2] = (-TC[:, 2]) % 90.0; TC = TC[mondrian.is_tilted(TC[:, 2])]
    for X, Y in ((TA, TB), (TA, TC), (TB, TC)):
        if len(X) and len(Y) and mondrian._sat_overlap(X, Y):
            return None
    return np.r_[TA, TB, TC]


def job3(args):
    """Three-piece arms, weighted random sampling within the budget (cand^3 is too many to enumerate): A, B, C drawn by
    prior weight among tight pieces (C may also be none), area screen, overlap check, exact fill.  Probe as in job."""
    n, S, tol, tl, budget, max_fills, seed, probe, probe_steps = args
    import mondrian
    S = S + probe
    P = PIECES
    cell = {}
    for i, q in enumerate(P):
        cell.setdefault(q['cell'], []).append(i)
    cand = [i for c, I in cell.items() if any(tight(P[i]['W'], S, tol + probe) or tight(P[i]['H'], S, tol + probe) for i in I)
            for i in I if max(P[i]['W'], P[i]['H']) <= S + 1e-9 and P[i]['m'] < n]
    rng = np.random.default_rng(seed)
    if not cand:
        return dict(n=n, S=S - probe, cand=0, pairs=0, arms=0, tried=0, sec=0.0, fills=[], near=[])
    wc = np.array([math.exp(-P[i]['gap'] / PRIOR[1]) for i in cand]); wc /= wc.sum()
    best, tried, fills, near, seen, t0, draws = None, 0, [], [], set(), time.time(), 0
    while time.time() - t0 < budget and tried < max_fills and draws < 200000:
        draws += 1
        ia, ib, ic = (int(x) for x in rng.choice(cand, 3, p=wc))
        A, B, C = P[ia], P[ib], P[ic]
        if A['m'] + B['m'] + C['m'] + math.floor(S * S - A['W'] * A['H'] - B['W'] * B['H'] - C['W'] * C['H'] + 1e-9) < n:
            continue
        va, vb, vc = (int(x) for x in rng.integers(4, size=3))
        key = (ia, ib, ic, va, vb, vc)
        if key in seen:
            continue
        seen.add(key)
        T = place3(A, B, C, va, vb, vc, S)
        if T is None:
            continue
        tried += 1
        t1 = time.time()
        L = mondrian.exact_fill(T, S, time_limit=tl)
        c = len(L)
        rec = dict(n=n, A=A['label'], B=B['label'], C=C['label'], v=[va, vb, vc], c=c, sec=round(time.time() - t1, 2))
        if probe and c >= n:
            lo, hi, Lh = S - 2 * probe, S, L
            for _ in range(probe_steps):
                mid = (lo + hi) / 2
                Tm = place3(A, B, C, va, vb, vc, mid)
                Lm = mondrian.exact_fill(Tm, mid, time_limit=tl) if Tm is not None else []
                if len(Lm) >= n:
                    hi, Lh = mid, Lm
                else:
                    lo = mid
            rec.update(S_min=hi, need=hi - (S - probe))
            near.append(dict(rec, sq=np.asarray(Lh).tolist()))
        fills.append(rec)
        if best is None or c > best[0]:
            best = (c, L)
    near.sort(key=lambda r: r['S_min'])
    base = dict(n=n, S=S - probe, near=near[:20], cand=len(cand), pairs=draws, arms=len(seen), tried=tried,
                sec=round(time.time() - t0, 1), fills=fills)
    if near and near[0]['need'] < -1e-9:
        x = near[0]
        return dict(base, count=n, A=x['A'], B=x['B'], S=x['S_min'], sq=x['sq'])
    return dict(base, count=best[0] if best else None)


def init(mmin, per_cell, source, cf_eps=0.0):
    global PIECES
    PIECES = load_pieces(mmin, per_cell, source, cf_eps)


def run(a):
    from pk import manifest, eventbench as E
    import scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'compose2', args=vars(a), procs=a.procs)
    scratch.scope('compose2')
    lo, hi = map(int, a.ns.split('-'))
    tasks = []
    only = {int(x) for x in a.only.split(',')} if a.only else None
    for n in range(lo, hi + 1):
        if only and n not in only:
            continue
        if a.target == 'now':
            f = E.frontier_now(n)
        else:
            f = E.frontier_at(n, a.target.split(':', 1)[1])
        if f >= math.ceil(math.sqrt(n) - 1e-12) - 1e-9 and not a.grid:
            continue
        tasks.append((n, f - a.below, a.tight, a.time_limit, a.budget, a.max_fills, a.seed + n, a.probe, a.probe_steps))
    print(f'{len(tasks)} targets, {len(load_pieces(a.mmin, a.per_cell, a.source))} pieces', flush=True)
    res = []
    with ProcessPoolExecutor(a.procs, initializer=init, initargs=(a.mmin, a.per_cell, a.source, a.cf_eps)) as ex, \
            open(os.path.join(a.out, 'results.jsonl'), 'w') as f, open(os.path.join(a.out, 'fills.jsonl'), 'w') as fl:
        for r in ex.map(job3 if a.three else job, tasks):
            sq = r.pop('sq', None)
            fl.write(''.join(json.dumps(x) + '\n' for x in r.pop('fills', []))); fl.flush()
            near = r.pop('near', [])
            if near:
                with open(os.path.join(a.out, 'near.jsonl'), 'a') as g:
                    g.write(''.join(json.dumps(x) + '\n' for x in near))
                r['need'] = near[0]['need']
            if sq is not None and not a.cf_eps:          # counterfactual layouts are not packings
                r.update(polish(r['n'], r['S'], sq, a.out))
            f.write(json.dumps(r) + '\n'); f.flush()
            print(f"{r['n']:4d} S {r['S']:.10f} cand {r['cand']:4d} arms {r['arms']:7d} tried {r['tried']:5d} {r['sec']:6.0f}s best count {r.get('count')}"
                  f"{'  need %.2e' % r['need'] if 'need' in r else ''}"
                  f"{'  HIT ' + str(r.get('s_pol')) if 's_pol' in r else ''}", flush=True)
    man.finish(summary=dict(targets=len(tasks)))


def rank(a):
    """Probe output -> cells worth searching.  For each n keep its best arms (need = smallest side above the frontier
    target at which the arm fits n); a cell's score = sum over n of exp(-need / tau), tau = the h improvement a ~1 CPU-h
    search plausibly buys (TAU x h; the (38, 2) control went 5.7996 -> 5.7642, 0.6 %, from 3 min to 1 CPU-h).  Both
    cells of an arm are credited (improving either piece closes the gap only partly: a coarse proxy)."""
    R = [json.loads(l) for d in a.dirs for l in open(os.path.join(d, 'near.jsonl'))]
    best = {}
    for r in R:                                         # per (n, cell pair): smallest need
        k = (r['n'], tuple(r['cA']), tuple(r['cB']))
        if k not in best or r['need'] < best[k]['need']:
            best[k] = r
    score, ns, tgt = {}, {}, {}
    for (n, cA, cB), r in best.items():
        for c, h in ((cA, r['hA']), (cB, r['hB'])):
            if c[0] == 0:                               # B = none (one-piece arm)
                continue
            if r['need'] > 0:                           # this cell's piece shrunk by need closes n (proxy): target h
                tgt[c] = max(tgt.get(c, 0.0), h - r['need'])
            w = math.exp(-max(r['need'], 0) / (a.tau * h))
            score[c] = score.get(c, 0) + w
            ns.setdefault(c, {})
            ns[c][n] = min(ns[c].get(n, 9), r['need'])
    hits = sorted({n for (n, *_), r in best.items() if r['need'] < -1e-9})
    print(f'{len(R)} near arms, {len(best)} (n, cell pair); n with need < 0 (hits): {hits}')
    out = []
    for c, sc in sorted(score.items(), key=lambda x: -x[1])[:a.top]:
        nn = sorted(ns[c].items(), key=lambda x: x[1])[:6]
        print(f'  cell m={c[0]:3d} d={c[1]:g}  score {sc:6.2f}  target h {tgt.get(c, 0):.6f}  n:need '
              + ' '.join(f'{n}:{v:.1e}' for n, v in nn))
        out.append(dict(m=c[0], d=c[1], score=sc, ns=nn, target=tgt.get(c, 0.0)))
    if a.out:
        json.dump(out, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--source', help='store piece source filter (SQL LIKE), default all'); r.add_argument('--ns', default='100-324')
    r.add_argument('--out', required=True); r.add_argument('--target', default='now')
    r.add_argument('--below', type=float, default=1e-9); r.add_argument('--tight', type=float, default=0.03)
    r.add_argument('--mmin', type=int, default=10); r.add_argument('--time-limit', type=float, default=10)
    r.add_argument('--procs', type=int, default=14); r.add_argument('--grid', action='store_true')
    r.add_argument('--budget', type=float, default=120, help='seconds of fills per n (weighted random arms)')
    r.add_argument('--max-fills', type=int, default=10 ** 9); r.add_argument('--seed', type=int, default=0)
    r.add_argument('--probe', type=float, default=0.0,
                   help='near-miss probe: fill at frontier + PROBE, bisect fitting arms to their smallest side (near.jsonl)')
    r.add_argument('--probe-steps', type=int, default=10)
    r.add_argument('--three', action='store_true', help='three-piece arms (LL, UR, UL corners)')
    r.add_argument('--only', help='comma list of n (within --ns)')
    r.add_argument('--cf-eps', type=float, default=0.0,
                   help='counterfactual: every piece scaled to h (1 - eps); hits say which cells are worth searching')
    r.add_argument('--per-cell', type=int, default=10 ** 9, help='pieces per (m, d) cell, best h first (seam diversity)')
    k = sp.add_parser('rank'); k.add_argument('dirs', nargs='+'); k.add_argument('--tau', type=float, default=0.003)
    k.add_argument('--top', type=int, default=40); k.add_argument('--out')
    q = sp.add_parser('polish'); q.add_argument('dirs', nargs='+'); q.add_argument('--out', required=True)
    q.add_argument('--max-need', type=float, default=2e-3); q.add_argument('--per-n', type=int, default=3)
    q.add_argument('--procs', type=int, default=14)
    a = ap.parse_args()
    {'run': run, 'rank': rank, 'polish': polish_near}[a.cmd](a)
