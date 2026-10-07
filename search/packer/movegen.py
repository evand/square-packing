#!/usr/bin/env python3
"""Structural moves on the tilted bands of a two-channel packing (count-pair table, step 2 of README "Next session").

A start packing (e.g. the s(110) record) has two tilted bands, L (10-45 deg mod 90) and B (45-80), near-axis junction squares
J (0.05-3 deg), and axis squares.  A proposal edits only the tilted set T = L + B + J, on each band's own lattice
(basis u = (cos a, sin a), v = (-sin a, cos a) at the band's median angle):
  add k / del k      grow / shrink a band by a connected cluster of k lattice sites (at the junction, the far end, or anywhere)
  transfer k         del k from one band at the junction, add k to the other at the junction
  row_add / row_del  a whole lattice line along one side of a band (thicker / thinner, longer / shorter)
  shift / rotate     the band rigidly (along u or v; by a few degrees about its centroid)
  drop_J             remove the junction squares (the axis refill decides what goes there)
Evaluation (--eval): 'fill' = ftmc.evaluate with reasons (relax T alone, exact axis fill at the smallest side holding n, penalty
squeeze, reject full lines); 'keep' (default) = keep the start's axis squares, fix the count locally, penalty squeeze from the
overlapping state.  Then slp2 finish -> side, jam, full lines, role counts (L, B, axis, other).

  movegen.py start.txt [start2.txt ...] --n 110 --k 11 --trials 200 --out runs/mv1/m [--procs 15] [--seed 1]   (--trials per start)
"""
import argparse, collections, json, math, os, random, tempfile
from ftmc import load, run_packer
from layout import full_lines
from gen import pen, wall_viol, write, mis_fill2
from jobpool import run_jobs

JOB_BUDGET = 300


FTOL = 0.03


def evaluate(T, n, k, s_cur, tmp, mu0, smax=0.6):
    """ftmc.evaluate with reasons: relax T alone, exact axis fill at the smallest scanned side holding n, penalty squeeze."""
    a, b, c = (os.path.join(tmp, x) for x in ('t.txt', 'tr.txt', 'f.txt'))
    sb = s_cur + 0.3
    write(a, sb, T)
    E, _ = run_packer(['relax', '--in', a, '--s', repr(sb), '--out', b, '--maxit', '5000'])
    if E > 1e-20:
        return dict(status='tilted overlap'), None
    _, T = load(b)
    s = s_cur - 0.03
    best = 0
    while s < s_cur + smax:
        if all(wall_viol(q, s) <= 1e-9 for q in T):
            ax = mis_fill2(s, T, FTOL)
            best = max(best, len(ax))
            if len(T) + len(ax) >= n:
                write(c, s, (T + ax)[:n])
                _, sq = run_packer(['relax', '--in', c, '--squeeze-pen', '--mu0', repr(mu0), '--out', b])
                if sq == float('inf'):
                    return dict(status='squeeze failed', s_fill=s), None
                if full_lines(b, k) > 0:
                    return dict(status='full line after squeeze', s_fill=s, soft=sq), None
                return dict(status='ok', s_fill=s, soft=sq, extra_axis=len(T) + len(ax) - n), b
        s += 0.04
    return dict(status='fill short', short=n - len(T) - best), None


def evaluate_keep(T, A, n, k, s0, tmp, mu0, loosen=1.01):
    """Keep the start's axis squares A: drop the ones most overlapped by T (too many squares), or insert axis squares at the
    least-overlap grid positions (too few); then penalty squeeze from the overlapping state at side s0 * loosen."""
    A = list(A)
    over = lambda q, S: sum(max(0.0, pen(q, o)) for o in S)
    while len(T) + len(A) > n:
        A.remove(max(A, key=lambda q: over(q, T)))
    added = 0
    while len(T) + len(A) < n:
        S = T + A
        grid = [(0.5 + i * (s0 - 1) / 60, 0.5 + j * (s0 - 1) / 60, 0.0) for i in range(61) for j in range(61)]
        A.append(min(grid, key=lambda q: over(q, S)))
        added += 1
    a, b = (os.path.join(tmp, x) for x in ('k.txt', 'kr.txt'))
    s = s0 * loosen
    sq = [(min(max(x, 0.5), s - 0.5), min(max(y, 0.5), s - 0.5), t) for x, y, t in T + A]
    write(a, s, sq)
    _, ss = run_packer(['relax', '--in', a, '--squeeze-pen', '--mu0', repr(mu0), '--out', b])
    if ss == float('inf'):
        return dict(status='squeeze failed', added=added), None
    if full_lines(b, k) > 0:
        return dict(status='full line after squeeze', soft=ss, added=added), None
    return dict(status='ok', soft=ss, added=added), b


def role(a):
    a %= 90
    d = min(a, 90 - a)
    if d < 3:
        return 'J' if d > 0.05 else 'A'
    # bands at ~25-32 and ~62-66 deg (s(110) record); 40-50 deg junction squares are 'O', not band members
    return 'L' if 12 < a < 40 else ('B' if 50 < a < 78 else 'O')


def split(C):
    out = collections.defaultdict(list)
    for q in C:
        out[role(q[2])].append(q)
    return out


def basis(band):
    a = sorted(q[2] for q in band)[len(band) // 2]
    t = math.radians(a)
    return a, (math.cos(t), math.sin(t)), (-math.sin(t), math.cos(t))


def occupied(p, band, r=0.6):
    return any(math.hypot(p[0] - q[0], p[1] - q[1]) < r for q in band)


def frontier(band, others, s, a=None):
    """Empty lattice-neighbour sites of the band, inside the box, not overlapping other tilted squares much."""
    a0, u, v = basis(band)
    a = a0 if a is None else a
    out = []
    for q in band:
        for d in (u, v, (-u[0], -u[1]), (-v[0], -v[1])):
            p = (q[0] + d[0], q[1] + d[1], a)
            if wall_viol(p, s) > 0.05 or occupied(p, band) or occupied(p, out, 0.3):
                continue
            if any(pen(p, o) > 0.15 for o in others):
                continue
            out.append(p)
    return out


def cluster(cands, k, anchor):
    """k candidates nearest to anchor (a point), greedy-connected."""
    if not cands:
        return []
    pick = [min(cands, key=lambda p: math.hypot(p[0] - anchor[0], p[1] - anchor[1]))]
    rest = [p for p in cands if p is not pick[0]]
    while len(pick) < k and rest:
        nxt = min(rest, key=lambda p: min(math.hypot(p[0] - q[0], p[1] - q[1]) for q in pick))
        pick.append(nxt); rest.remove(nxt)
    return pick


def anchor_point(band, where, junction, rng):
    if where == 'junction':
        return junction
    if where == 'far':
        q = max(band, key=lambda q: math.hypot(q[0] - junction[0], q[1] - junction[1]))
        return q[:2]
    q = band[rng.randrange(len(band))]
    return q[:2]


def snap_band(band, a, anchor, s, offset=(0.0, 0.0)):
    """Re-place a band on one lattice at angle a (deg) anchored at `anchor` (+ offset); merge collisions."""
    t = math.radians(a)
    u, v = (math.cos(t), math.sin(t)), (-math.sin(t), math.cos(t))
    ax, ay = anchor[0] + offset[0], anchor[1] + offset[1]
    out = []
    for x, y, _ in band:
        dx, dy = x - ax, y - ay
        i, j = round(dx * u[0] + dy * u[1]), round(dx * v[0] + dy * v[1])
        p = (ax + i * u[0] + j * v[0], ay + i * u[1] + j * v[1], a)
        if wall_viol(p, s) <= 0.05 and not occupied(p, out, 0.5):
            out.append(p)
    return out


def rebuild_band(band, a, s, offset, others, dil=0.3):
    """Fresh lattice patch at angle a over the band's region (centres within dil of the old band's squares)."""
    t = math.radians(a)
    u, v = (math.cos(t), math.sin(t)), (-math.sin(t), math.cos(t))
    ax, ay = band[0][0] + offset[0], band[0][1] + offset[1]
    out = []
    R = int(2 * s)
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            p = (ax + i * u[0] + j * v[0], ay + i * u[1] + j * v[1], a)
            if not (0 < p[0] < s and 0 < p[1] < s) or wall_viol(p, s) > 0.05:
                continue
            if min(math.hypot(p[0] - q[0], p[1] - q[1]) for q in band) > 0.5 + dil:
                continue
            if any(pen(p, o) > 0.15 for o in others):
                continue
            out.append(p)
    return out


def move_focus(G, s, rng, junction, b):
    """Moves aimed at one band b: the generic ones on b, plus 'uniform' (one angle, one lattice) and 'rebuild'."""
    kind = rng.choice(['uniform', 'uniform', 'rebuild', 'rebuild', 'generic'])
    others = [q for n2 in 'LBJ' if n2 != b for q in G[n2]]
    if kind == 'generic':
        return move(G, s, rng, junction, band=b)
    lo, hi = (20, 34) if b == 'L' else (56, 70)
    a = rng.uniform(lo, hi)
    if kind == 'uniform':
        anchor = min(G[b], key=lambda q: q[0] + q[1]) if b == 'L' else min(G[b], key=lambda q: q[1] - q[0])
        G[b] = snap_band(G[b], a, anchor, s, (rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)))
        return f'uniform {b} {a:.1f} -> {len(G[b])}'
    G[b] = rebuild_band(G[b], a, s, (rng.uniform(0, 1), rng.uniform(0, 1)), others, dil=rng.choice([0.0, 0.3, 0.6]))
    return f'rebuild {b} {a:.1f} -> {len(G[b])}'


def move(G, s, rng, junction, band=None):
    """Apply one random move to G = {'L': [...], 'B': [...], 'J': [...]}; returns a description."""
    kind = rng.choice(['add', 'del', 'transfer', 'transfer', 'row_add', 'row_del', 'shift', 'rotate', 'drop_J'])
    b = band or rng.choice('LB')
    o = 'B' if b == 'L' else 'L'
    others = lambda name: [q for n2 in 'LBJ' if n2 != name for q in G[n2]]
    if kind in ('add', 'del'):
        k = rng.choice([1, 1, 2, 3])
        where = rng.choice(['junction', 'far', 'any'])
        anc = anchor_point(G[b], where, junction, rng)
        if kind == 'add':
            new = cluster(frontier(G[b], others(b), s), k, anc)
            G[b] += new
            return f'add {b} {len(new)} {where}'
        rm = cluster(list(G[b]), k, anc)
        G[b] = [q for q in G[b] if q not in rm]
        return f'del {b} {len(rm)} {where}'
    if kind == 'transfer':
        k = rng.choice([1, 1, 2, 3])
        rm = cluster(list(G[b]), k, junction)
        G[b] = [q for q in G[b] if q not in rm]
        new = cluster(frontier(G[o], others(o), s), len(rm), junction)
        G[o] += new
        return f'transfer {b}>{o} {len(rm)}/{len(new)}'
    if kind in ('row_add', 'row_del'):
        a, u, v = basis(G[b])
        d = rng.choice([u, v, (-u[0], -u[1]), (-v[0], -v[1])])
        proj = lambda q: q[0] * d[0] + q[1] * d[1]
        if kind == 'row_add':
            m = max(proj(q) for q in G[b])
            edge = [q for q in G[b] if proj(q) > m - 0.5]
            new = [(q[0] + d[0], q[1] + d[1], a) for q in edge]
            new = [p for p in new if wall_viol(p, s) <= 0.05 and not any(pen(p, x) > 0.15 for x in others(b))]
            G[b] += new
            return f'row_add {b} {len(new)}'
        m = max(proj(q) for q in G[b])
        rm = [q for q in G[b] if proj(q) > m - 0.5]
        if len(rm) < len(G[b]) - 3:
            G[b] = [q for q in G[b] if q not in rm]
        return f'row_del {b} {len(rm)}'
    if kind == 'shift':
        a, u, v = basis(G[b])
        d = rng.choice([u, v]); t = rng.choice([-1, 1]) * rng.uniform(0.15, 1.0)
        G[b] = [(x + t * d[0], y + t * d[1], an) for x, y, an in G[b]]
        return f'shift {b} {t:+.2f}'
    if kind == 'rotate':
        ph = rng.choice([-1, 1]) * rng.uniform(1, 5)
        cx = sum(q[0] for q in G[b]) / len(G[b]); cy = sum(q[1] for q in G[b]) / len(G[b])
        c, sn = math.cos(math.radians(ph)), math.sin(math.radians(ph))
        G[b] = [(cx + c * (x - cx) - sn * (y - cy), cy + sn * (x - cx) + c * (y - cy), an + ph) for x, y, an in G[b]]
        return f'rotate {b} {ph:+.1f}'
    n = len(G['J'])
    G['J'] = []
    return f'drop_J {n}'


def job(args):
    start, n, k, seed, out, mode, focus = args
    from slp2 import slp2
    from slp import save
    rng = random.Random(seed)
    s0, C = load(start)
    G = split(C)
    L, B = G['L'], G['B']
    # junction: midpoint of the closest L-B pair
    p, q = min(((p, q) for p in L for q in B), key=lambda t: math.hypot(t[0][0] - t[1][0], t[0][1] - t[1][1]))
    junction = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
    A0 = list(G['A'])
    G = {'L': list(L), 'B': list(B), 'J': list(G['J']) + list(G['O'])}      # 'O' squares stay in the tilted set
    if focus:
        desc = [move_focus(G, s0, rng, junction, focus) for _ in range(rng.choice([1, 1, 2]))]
    else:
        desc = [move(G, s0, rng, junction) for _ in range(rng.choice([1, 1, 2]))]
    T = G['L'] + G['B'] + G['J']
    mu0 = rng.choice([1e2, 1e3])
    tmp = tempfile.mkdtemp()
    if mode == 'fill':
        ev, path = evaluate(T, n, k, s0, tmp, mu0)
    else:
        ev, path = evaluate_keep(T, A0, n, k, s0, tmp, mu0)
    rec = dict(start=start, seed=seed, moves=desc, mu0=mu0, nT=(len(G['L']), len(G['B']), len(G['J'])), **ev)
    if path is None:
        return rec
    from rigid import load as load_rad                 # slp2 works in radians (ftmc.load returns degrees)
    s2, sq = load_rad(path)
    info = {}
    s, sq2, _ = slp2(s2, sq, R=1e-3, rmin=1e-8, budget=JOB_BUDGET, info=info)
    fn = f'{out}_{seed}.txt'
    save(fn, s, sq2)
    c = collections.Counter(role(math.degrees(q[2])) for q in sq2)
    rec.update(status='ok', s=s, jam=info['jammed'], lines=full_lines(fn, k), path=fn,
               counts=dict(L=c['L'], B=c['B'], A=c['A'] + c['J'], O=c['O']))
    return rec


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('start', nargs='+'); ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--trials', type=int, default=100); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--eval', choices=['fill', 'keep'], default='keep'); ap.add_argument('--focus', choices=['L', 'B'], default=None); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or '.', exist_ok=True)
    # trials per start; seeds unique across starts
    jobs = [(st, a.n, a.k, a.seed * 100000 + i * a.trials + t, a.out, a.eval, a.focus)
            for i, st in enumerate(a.start) for t in range(a.trials)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=3 * JOB_BUDGET,
                   on_timeout=lambda j: dict(seed=j[3], status='killed'),
                   on_error=lambda j, e: dict(seed=j[3], status=f'failed: {e}'))
    res = [r if r is not None else dict(seed=j[3], status='failed') for r, j in zip(res, jobs)]
    json.dump(dict(start=a.start, n=a.n, k=a.k, eval=a.eval, trials=res), open(a.out + '.json', 'w'), indent=0)
    ok = [r for r in res if r.get('status') == 'ok']
    print(f'{len(ok)}/{len(res)} evaluated; statuses {dict(collections.Counter(r["status"][:30] for r in res))}')
    tab = collections.defaultdict(list)
    for r in ok:
        if r['lines'] == 0:
            tab[(r['counts']['L'], r['counts']['B'], r['counts']['O'])].append(r)
    print('(L, B, other): trials, best side, best moves  [no full lines]')
    for key in sorted(tab, key=lambda kk: min(r['s'] for r in tab[kk])):
        R = tab[key]; b = min(R, key=lambda r: r['s'])
        print(f'  {key}: {len(R):3d}  {b["s"]:.7f}  {b["moves"]}')
