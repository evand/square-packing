#!/usr/bin/env python3
"""Rotated lattice block + exact axis fill (10-09; the 126-trio generator from TODO): an a x b block of unit squares on a
lattice at angle theta, placed in the box of side s, rest filled with axis squares by gen.mis_fill2 (exact MILP).
Random search over (theta, a, b, position) for fills reaching n; each hit is fq-quenched (screening, f64).

  blockgen.py --n N --s S [--theta-set 45,30,...] [--evals K] [--procs P] [--out DIR] [--seed k]

Positive control: n = 126, s = 11.7427 must find Ryan Xu's 45-degree 6 x 6 diamond (s = 15/2 + 3*sqrt 2).
"""
import argparse, math, os, random, tempfile
from jobpool import run_jobs
import mcmin, hop
from gen import mis_fill2
DEFICIT, LOOSEN = 0, '1.0'


def block(s, th, a, b, ox, oy, rows=None):
    """Squares of an a x b lattice block at angle th (deg); (ox, oy) = centre of the block as a fraction of the free range.
    rows = [(length, shift)] per row (staircase block: row j has its own length and offset along the row direction)."""
    t = math.radians(th); c, si = math.cos(t), math.sin(t)
    rows = rows or [(a, 0.0)] * b
    b = len(rows)
    P = [((i + d - (L - 1) / 2) * c - (j - (b - 1) / 2) * si, (i + d - (L - 1) / 2) * si + (j - (b - 1) / 2) * c)
         for j, (L, d) in enumerate(rows) for i in range(L)]
    w = 0.5 * (abs(c) + abs(si))                            # half extent of one rotated square
    xs = [p[0] for p in P]; ys = [p[1] for p in P]
    lo_x, hi_x = w - min(xs), s - w - max(xs)
    lo_y, hi_y = w - min(ys), s - w - max(ys)
    if lo_x > hi_x or lo_y > hi_y:
        return None
    if isinstance(ox, str):                                  # 'c+dx,dy': block centre at box centre + (dx, dy) (absolute)
        dx, dy = map(float, ox[2:].split(','))
        cx, cy = s / 2 + dx, s / 2 + dy
        if not (lo_x - 1e-9 <= cx <= hi_x + 1e-9 and lo_y - 1e-9 <= cy <= hi_y + 1e-9):
            return None
    else:
        cx, cy = lo_x + ox * (hi_x - lo_x), lo_y + oy * (hi_y - lo_y)
    return [(cx + x, cy + y, th % 90) for x, y in P]


def evaluate(j):
    n, s, th, a, b, ox, oy, out = j[:8]
    rows = j[8] if len(j) > 8 else None
    T = block(s, th, a, b, ox, oy, rows)
    if T is None:
        return None
    try:
        ax = mis_fill2(s, T, tol=1e-6, time_limit=20)
    except Exception:
        return None
    tot = len(T) + len(ax)
    r = None
    if tot >= n - DEFICIT:
        allsq = T + [(q[0], q[1], 0.0) for q in ax][:n - len(T)]
        if len(allsq) < n:                                   # deficit: missing squares at the largest clearance holes
            hs = mcmin.holes(s, allsq)
            while len(allsq) < n and hs:
                x, y, _ = hs.pop(0); allsq.append((x, y, 0.0))
            if len(allsq) < n:
                return (tot, None, th, a, b, ox, oy, rows)
        tmp = tempfile.mkdtemp()
        hop.EXTRA[:] = ['--loosen', LOOSEN]
        q = hop.quench(s, allsq, tmp)
        if q:
            r = q[0]
            mcmin.write_deg(f'{out}/bg_n{n}_th{th:.3f}_{a}x{b}{"s" if rows else ""}_{r:.9f}.txt', q[0], q[1])
    return (tot, r, th, a, b, ox, oy, rows)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--s', type=float, required=True)
    ap.add_argument('--theta-set', default=''); ap.add_argument('--evals', type=int, default=400)
    ap.add_argument('--procs', type=int, default=2); ap.add_argument('--out', default='runs/bg'); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--amax', type=int, default=0, help='max block edge (default floor(s) - 1)')
    ap.add_argument('--amin', type=int, default=2); ap.add_argument('--theta-range', default='')
    ap.add_argument('--stair', type=float, default=0.0, help='share of staircase blocks (per-row length a - {0,1,2}, shift)')
    ap.add_argument('--deficit', type=int, default=0); ap.add_argument('--loosen', default='1.0')
    ap.add_argument('--fixed-ab', default='', help='a,b: fixed block size')
    a = ap.parse_args()
    DEFICIT, LOOSEN = a.deficit, a.loosen
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    ths = [float(x) for x in a.theta_set.split(',')] if a.theta_set else None
    amax = a.amax or int(a.s) - 1
    jobs = []
    for _ in range(a.evals):
        if a.theta_range:
            lo, hi = map(float, a.theta_range.split(',')); th = rng.uniform(lo, hi)
        else:
            th = rng.choice(ths) if ths else rng.uniform(5, 45)
        A, B = (map(int, a.fixed_ab.split(',')) if a.fixed_ab else (rng.randint(a.amin, amax), rng.randint(a.amin, amax)))
        if rng.random() < 0.5:                                # centre offsets on the quarter lattice (45-degree records)
            ox, oy = f'c+{rng.randint(-8, 8) / 4},{rng.randint(-8, 8) / 4}', None
        else:
            ox, oy = rng.choice([0.0, 0.5, 1.0, rng.random()]), rng.choice([0.0, 0.5, 1.0, rng.random()])
        rows = None
        if rng.random() < a.stair:
            rows = [(A - rng.choice([0, 0, 1, 2]), rng.choice([0.0, rng.uniform(-1, 1)])) for _ in range(B)]
        jobs.append((a.n, a.s, th, A, B, ox, oy, a.out, rows))
    res = [r for r in run_jobs(evaluate, jobs, procs=a.procs, timeout=300) if isinstance(r, tuple)]
    res.sort(key=lambda r: (-r[0], r[1] if r[1] is not None else 99))
    print(f'{len(res)} fills; best totals: ' + ', '.join(str(r[0]) for r in res[:10]))
    q = sorted(r[1] for r in res if r[1] is not None)
    print(f'quenched {len(q)}; best sides {q[:8]}')
    for r in res[:15]:
        print(f'  total {r[0]:4d} quench {r[1]}  th {r[2]:.2f} block {r[3]}x{r[4]} pos {r[5]} {r[6]} rows {r[7]}')
