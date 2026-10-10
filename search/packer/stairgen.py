#!/usr/bin/env python3
"""Staircase blocks with a continuous optimiser (10-09).  Tilted part = C columns of unit squares at angle theta:
column i (i = 0..C-1) runs along v = (-sin, cos) at u-position i, holds L_i squares starting at offset v0_i:
    centre(i, j) = origin + i * u + (v0_i + j) * v,     u = (cos theta, sin theta).
Plain rectangles (all v0_i equal), strips (C = 1 or L_i = 1) and Xu-type diamonds (126: 6 columns of 6, offsets
0 / 0 / .77 / .68 / .71 / .81) are special cases.  The axis part is the exact fill (gen.mis_fill2).

Search: cross-entropy over the continuous parameters (placement fractions, v0_i) for fixed (theta, C, L).  Objective:
the structure in a box of side s + exact axis fill, completed to n squares at clearance holes, fq-quenched; score =
quenched side (results jammed at an integer side are rejected).  Always n squares, scored by side (Evan 10-09).
The exact fill alone cannot express the new-style records (ry-xu 126: axis squares in pockets against the block and
against near-axis rotated squares), so the quench does the last step; seeds then go to the explorer.

  stairgen.py --n N --s S --theta T --cols L1,L2,... [--pop 40 --iters 15 --restarts R --procs P]
"""
import argparse, math, os, random, tempfile
import numpy as np
from jobpool import run_jobs
import mcmin, hop
from gen import mis_fill2, corners


def stair(theta, origin, v0, L):
    t = math.radians(theta); u = (math.cos(t), math.sin(t)); v = (-math.sin(t), math.cos(t))
    o = (float(origin[0]), float(origin[1])); v0 = [float(q) for q in v0]       # plain floats: repr() goes into files
    return [(o[0] + i * u[0] + (v0[i] + j) * v[0], o[1] + i * u[1] + (v0[i] + j) * v[1], float(theta % 90))
            for i in range(len(L)) for j in range(L[i])]


def place(s, theta, f, v0, L):
    """origin from fractions f = (fx, fy) of the feasible range (block inside the box); None if the block is too big"""
    T0 = stair(theta, (0.0, 0.0), v0, L)
    P = [p for q in T0 for p in corners(*q)]
    x0, x1 = min(p[0] for p in P), max(p[0] for p in P); y0, y1 = min(p[1] for p in P), max(p[1] for p in P)
    if x1 - x0 > s or y1 - y0 > s:
        return None
    fx, fy = (min(max(v, 0.0), 1.0) for v in f)
    return stair(theta, (-x0 + fx * (s - (x1 - x0)), -y0 + fy * (s - (y1 - y0))), v0, L)


def inside(s, T):
    m = min(min(min(p[0], p[1], s - p[0], s - p[1]) for p in corners(*q)) for q in T)
    return m


def score(s, T, reach=4):
    if not T or inside(s, T) < -1e-9:
        return -1, None
    ax = mis_fill2(s, T, tol=1e-6, reach=reach, time_limit=20)
    return len(T) + len(ax), ax


def finish(n, s, T, ax, loosen='1.0', tries=8, seed=0):
    """complete to n (missing squares at random choices among the 10 largest holes, else random) and fq-quench; several
    variants, best non-grid result (a side within 1e-6 of an integer is a jammed axis line: useless as a seed)"""
    base = T + [(q[0], q[1], 0.0) for q in ax][:n - len(T)]
    hs0 = mcmin.holes(s, base) if len(base) < n else []
    rr = random.Random(seed)
    best = None
    for k in range(tries if len(base) < n else 1):
        allsq = list(base); hs = list(hs0[:10]); rr.shuffle(hs) if k else None
        while len(allsq) < n and hs:
            x, y, _ = hs.pop(0); allsq.append((x, y, rr.choice([0.0, rr.uniform(0, 90)]) if k else 0.0))
        while len(allsq) < n:
            allsq.append((rr.uniform(0.5, s - 0.5), rr.uniform(0.5, s - 0.5), rr.uniform(0, 90)))
        hop.EXTRA[:] = ['--loosen', loosen if k % 2 == 0 else '1.02']
        try:
            r = hop.quench(s, allsq, tempfile.mkdtemp())
        except Exception:                               # fq crash / empty output
            r = None
        if r and abs(r[0] - round(r[0])) > 1e-6 and (best is None or r[0] < best[0]):
            best = r
    return best


def smin(n, theta, f, v0, L, lo, hi, tol=2e-4):
    """smallest side (to tol) in [lo, hi] at which the structure (placed by fractions f) + exact axis fill holds n
    squares; returns (s, packing) or (None, None).  count(s) is not monotone in general (placement scales with s);
    bisection finds a crossing, good enough for a seed."""
    def full(s):
        T = place(s, theta, f, v0, L)
        if T is None:
            return None
        k, ax = score(s, T)
        return (T + [(q[0], q[1], 0.0) for q in ax])[:n] if k >= n else None
    P = full(hi)
    if P is None:
        return None, None
    while hi - lo > tol:
        m = (lo + hi) / 2
        Q = full(m)
        if Q is None:
            lo = m
        else:
            hi, P = m, Q
    return hi, P


def qside(n, s, theta, x, L, tries=3, seed=0):
    """objective: the structure placed in a box of side s, exact axis fill, completed to n squares (holes, then random;
    a few variants), fq-quenched: best quenched side, ignoring results jammed at an integer side (a full axis line).
    Always an n-square packing scored by its side (Evan 10-09)."""
    T = place(s, theta, x[:2], x[2:], L)
    if T is None:
        return np.inf, None
    k, ax = score(s, T)
    if k < 0:
        return np.inf, None
    r = finish(n, s, T, ax, tries=tries, seed=seed)
    return (r[0], r) if r else (np.inf, None)


def cem(j):
    """one cross-entropy restart over (placement fractions, offsets) minimising qside; returns (side, params, path)"""
    n, s, theta, L, pop, iters, seed, out = j
    rng = np.random.default_rng(seed)
    C = len(L)
    mu = np.concatenate([rng.uniform(0, 1, 2), rng.uniform(-0.5, 0.5, C)])
    sd = np.concatenate([[0.4, 0.4], np.full(C, 0.5)])
    best = (np.inf, None, None)
    for it in range(iters):
        X = mu + sd * rng.standard_normal((pop, 2 + C))
        sc = []
        for q, x in enumerate(X):
            v, r = qside(n, s, theta, x, L, seed=seed * 7919 + it * 101 + q)
            sc.append(v)
            if v < best[0]:
                best = (v, x.copy(), r)
        el = X[np.argsort(np.array(sc))[:max(4, pop // 6)]]
        mu, sd = el.mean(0), np.maximum(el.std(0), 0.01 * (0.7 ** it))
    v, x, r = best
    if r is None:
        return (None, None, None)
    path = f'{out}/sg_n{n}_th{theta:g}_C{C}_{v:.9f}_{seed}.txt'
    mcmin.write_deg(path, r[0], r[1])
    return (float(v), [round(float(q), 4) for q in x], path)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--s', type=float, required=True, help='box side for the construction (around the target)')
    ap.add_argument('--theta', type=float, default=45.0); ap.add_argument('--cols', required=True)
    ap.add_argument('--pop', type=int, default=40); ap.add_argument('--iters', type=int, default=15)
    ap.add_argument('--restarts', type=int, default=4); ap.add_argument('--deficit', type=int, default=1)
    ap.add_argument('--procs', type=int, default=2); ap.add_argument('--out', default='runs/sg'); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    L = [int(v) for v in a.cols.split(',')]
    jobs = [(a.n, a.s, a.theta, L, a.pop, a.iters, a.seed * 1000 + r, a.out) for r in range(a.restarts)]
    res = [r for r in run_jobs(cem, jobs, procs=a.procs, timeout=7200) if isinstance(r, tuple)]
    for r in sorted(res, key=lambda r: r[0] if r[0] is not None else 1e9):
        print(f'side {r[0]} params {r[1]} {r[2] or ""}')
