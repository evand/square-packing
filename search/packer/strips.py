#!/usr/bin/env python3
"""Diagonal strips / staircase band + exact axis fill (10-09; Evan: test on s(110), the landscape we know).

Band = K parallel strips of face-touching unit squares at angle theta, strips 1 apart (faces touching) along the normal,
the band's centre line at signed offset t from the box centre, strip k shifted along its own direction by d_k (a
staircase when the shifts differ).  Each strip keeps only squares fully inside the box (optionally minus `trim` squares
at either end).  The rest is filled with axis squares (gen.mis_fill2, exact MILP).  Fills reaching n are fq-quenched.

  strips.py --n N --s S [--evals K] [--theta-range 20,45] [--kmax 4] [--procs P] [--out DIR] [--seed k]
"""
import argparse, math, os, random, tempfile
from jobpool import run_jobs
import mcmin, hop
from gen import mis_fill2
DEFICIT, LOOSEN, CUT = 0, '1.0', False


def band(s, th, t, shifts, trims):
    r = math.radians(th); u = (math.cos(r), math.sin(r)); nv = (-math.sin(r), math.cos(r))
    w = 0.5 * (abs(u[0]) + abs(u[1]))
    K = len(shifts); out = []
    for k in range(K):
        off = t + (k - (K - 1) / 2)
        cx, cy = s / 2 + off * nv[0], s / 2 + off * nv[1]
        row = []
        for j in range(-2 * int(s) - 2, 2 * int(s) + 3):
            x, y = cx + (j + shifts[k]) * u[0], cy + (j + shifts[k]) * u[1]
            if w <= x <= s - w and w <= y <= s - w:
                row.append((x, y, th % 90))
        a, b = trims[k]
        row = row[a:len(row) - b] if len(row) > a + b else []
        out += row
    return out


def cuts_all(s, T):
    """no unit-wide axis band (row or column) is free of tilted squares, so no full line of axis squares can form"""
    from gen import corners
    for ax_ in (0, 1):
        iv = sorted((min(p[ax_] for p in corners(*q)), max(p[ax_] for p in corners(*q))) for q in T)
        reach = 0.0
        for lo, hi in iv:
            if lo - reach >= 1.0:
                return False
            reach = max(reach, hi)
        if s - reach >= 1.0:
            return False
    return True


def evaluate(j):
    n, s, th, t, shifts, trims, out, quench = j
    T = band(s, th, t, shifts, trims)
    try:
        ax = mis_fill2(s, T, tol=1e-6, time_limit=20)
    except Exception:
        return None
    tot = len(T) + len(ax)
    r = None
    if CUT and not cuts_all(s, T):
        return (tot, None, th, t, shifts, trims, len(T))
    if tot >= n - DEFICIT and quench:
        allsq = T + [(q[0], q[1], 0.0) for q in ax][:n - len(T)]
        rr = random.Random(hash((th, t)) & 0xffff)
        hs = mcmin.holes(s, allsq) if len(allsq) < n else []
        while len(allsq) < n:                    # missing squares: largest clearance holes first, then random
            if hs:
                x, y, _ = hs.pop(0); allsq.append((x, y, rr.uniform(0, 90)))
            else:
                allsq.append((rr.uniform(0.5, s - 0.5), rr.uniform(0.5, s - 0.5), rr.uniform(0, 90)))
        tmp = tempfile.mkdtemp()
        hop.EXTRA[:] = ['--loosen', LOOSEN]
        q = hop.quench(s, allsq, tmp)
        if q:
            r = q[0]
            mcmin.write_deg(f'{out}/st_n{n}_th{th:.2f}_K{len(shifts)}_{r:.9f}.txt', q[0], q[1])
    return (tot, r, th, t, shifts, trims, len(T))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--s', type=float, required=True)
    ap.add_argument('--evals', type=int, default=400); ap.add_argument('--theta-range', default='20,45')
    ap.add_argument('--kmax', type=int, default=4); ap.add_argument('--procs', type=int, default=2)
    ap.add_argument('--out', default='runs/st'); ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--no-quench', action='store_true')
    ap.add_argument('--deficit', type=int, default=0, help='also quench fills short by up to this many (holes, then random)')
    ap.add_argument('--loosen', default='1.0'); ap.add_argument('--cut', action='store_true', help='require the band to cut every row and column')
    ap.add_argument('--tmax', type=float, default=None, help='max |band offset| (default s/2 - 1.5)')
    a = ap.parse_args()
    DEFICIT, LOOSEN, CUT = a.deficit, a.loosen, a.cut
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    lo, hi = map(float, a.theta_range.split(','))
    jobs = []
    for _ in range(a.evals):
        th = rng.choice([lo, hi, rng.uniform(lo, hi)])
        K = rng.randint(1, a.kmax)
        stair = rng.random() < 0.5
        d = rng.uniform(-0.5, 0.5)
        shifts = [(d + (k * rng.uniform(-0.5, 0.5) if stair else 0.0)) for k in range(K)]
        trims = [(rng.choice([0, 0, 1, 2]), rng.choice([0, 0, 1, 2])) for _ in range(K)]
        s_half = a.s / 2 - 1.5 if a.tmax is None else a.tmax
        t = rng.uniform(-s_half, s_half) if s_half > 0 and rng.random() < 0.5 else 0.0
        jobs.append((a.n, a.s, th, t, shifts, trims, a.out, not a.no_quench))
    res = [r for r in run_jobs(evaluate, jobs, procs=a.procs, timeout=300) if isinstance(r, tuple)]
    res.sort(key=lambda r: (-r[0], r[1] if r[1] is not None else 99))
    print(f'{len(res)} fills; totals: ' + ', '.join(str(r[0]) for r in res[:12]))
    for r in res[:12]:
        print(f'  total {r[0]:4d} (tilted {r[6]}) quench {r[1]}  th {r[2]:.2f} t {r[3]:+.2f} K {len(r[4])} shifts '
              + ' '.join(f'{x:+.2f}' for x in r[4]) + f' trims {r[5]}')
    hits = sorted(r[1] for r in res if r[1] is not None)
    print(f'quenched: {len(hits)}; best {hits[:8]}')
