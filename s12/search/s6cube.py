#!/usr/bin/env python3
"""s6cube.py -- delta*(eps*v) over many directions v of the small-angle cube (2026-09-20).

E1: for each direction v (signed tilts, max-norm 1) measure delta*(eps*v) at several eps, giving the
    local order and the coefficient.  Lower bounds on delta* (s6local ascent; structured starts at
    the middle eps, continuation to the others).
E2: at the optimum, the dual of the fixed-assignment LP (s6exact.rows_numeric): support and weights.

    setsid nohup python3 -u search/s6cube.py run --n 6 --T 3 --out runs/s6cube_n6.jsonl > runs/s6cube_n6.log 2>&1 &
    python3 search/s6cube.py report runs/s6cube_n6.jsonl
"""
import argparse
import itertools
import json
import math
import os
import sys

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6local                                                       # noqa: E402
import s6exact                                                       # noqa: E402


def directions(n, rng, ndense, nsparse):
    out = []
    for k in range(n // 2 + 1):                       # sign patterns, k minus signs
        out.append(('sign%d' % k, [-1.0] * k + [1.0] * (n - k)))
    for j in range(1, n - 1):                         # j zeros, the rest +1
        out.append(('zeros%d' % j, [0.0] * j + [1.0] * (n - j)))
    for _ in range(ndense):
        v = rng.normal(size=n); v /= np.max(np.abs(v))
        out.append(('dense', list(v)))
    for j in (1, 2, 3, 4):
        for _ in range(nsparse):
            v = rng.normal(size=n); v[:j] = 0.0; v /= np.max(np.abs(v))
            out.append(('sparse%d' % j, list(v)))
    return out


def _job(arg):
    kind, v, n, T, epsdeg, limit, jitters, seed = arg
    rng = np.random.default_rng(seed)
    res = {}
    mid = epsdeg[len(epsdeg) // 2]
    order = [mid] + [e for e in epsdeg if e > mid] + [e for e in reversed(epsdeg) if e < mid]
    carried, best_cfg = [], {}
    for e in order:
        theta = [math.radians(e) * x for x in v]
        F = s6local.Fixed(theta, T)
        starts = s6local.sampled_starts(n, T, rng, limit, jitters, 0.08) if e == mid else []
        starts += carried if e != mid else []
        if e < mid and mid in best_cfg:
            starts += best_cfg[mid]
        runs = []
        for (X, Y) in starts:
            val, bX, bY, asg = F.ascent(X, Y)
            if bX is not None:
                runs.append((val, bX, bY, asg))
        runs.sort(key=lambda r: -r[0])
        top = runs[:20]
        best_cfg[e] = [(r[1], r[2]) for r in top]
        carried = best_cfg[e]
        val, bX, bY, asg = top[0]
        nhit = sum(1 for r in runs if r[0] >= val - 1e-9)
        rows = s6exact.rows_numeric(F, asg)
        w, dlt = s6exact.solve_dual_numeric(rows, n)
        supp = [(rows[k][3], float(w[k])) for k in range(len(rows)) if w[k] > 1e-7]
        res[str(e)] = {'value': float(val), 'lp': float(dlt), 'nhit': nhit, 'nruns': len(runs),
                       'X': [float(x) for x in bX], 'Y': [float(y) for y in bY],
                       'supp': [[list(map(lambda q: q if isinstance(q, str) else int(q), tag)), ww]
                                for (tag, ww) in supp]}
    return {'kind': kind, 'v': [float(x) for x in v], 'eps': epsdeg, 'res': res}


def cmd_run(a):
    from multiprocessing import Pool
    rng = np.random.default_rng(a.seed)
    dirs = directions(a.n, rng, a.ndense, a.nsparse)
    epsdeg = [float(s) for s in a.eps.split(',')]
    jobs = [(k, v, a.n, a.T, epsdeg, a.limit, a.jitters, a.seed * 100003 + i)
            for i, (k, v) in enumerate(dirs)]
    print(f"# {len(jobs)} directions, eps={epsdeg}, limit={a.limit} x jitters={a.jitters}", flush=True)
    with Pool(a.nproc) as pool, open(a.out, 'w') as f:
        for i, r in enumerate(pool.imap_unordered(_job, jobs)):
            f.write(json.dumps(r) + '\n'); f.flush()
            print(f"  done {i + 1}/{len(jobs)}  {r['kind']}", flush=True)


def chain_shape(supp, n):
    """main chain = rows with weight >= 0.2 of the total pair weight (which is 1)."""
    big = [(tag, w) for (tag, w) in supp if w >= 0.2]
    walls = [t for (t, w) in big if t[0] in ('lo', 'hi')]
    pairs = [t for (t, w) in big if t[0] == 'pair']
    axes = set([t[1] for t in walls])
    return len(walls), len(pairs), ''.join(sorted(axes))


def cmd_report(a):
    recs = [json.loads(l) for l in open(a.file)]
    print(f"# {len(recs)} directions")
    print("# kind      s3     s4    |  delta*/eps^2 at each eps  |  local order (last two eps)  | main chain (walls,pairs,axis) | support")
    rows = []
    for r in recs:
        v = np.array(r['v']); s = np.sort(np.abs(v))
        eps = r['eps']
        vals = [r['res'][str(e)]['value'] for e in eps]
        q = [vals[k] / math.radians(eps[k]) ** 2 for k in range(len(eps))]
        if vals[0] < -1e-12 and vals[1] < -1e-12:
            order = math.log(vals[1] / vals[0]) / math.log(eps[1] / eps[0])
        else:
            order = float('nan')
        sup = r['res'][str(eps[0])]['supp']
        sup = [((t[0],) + tuple(t[1:]), w) for (t, w) in sup]
        shape = chain_shape(sup, len(v))
        rows.append((r['kind'], s, q, order, shape, len(sup), r))
    rows.sort(key=lambda x: (x[0], -x[2][0]))
    for (kind, s, q, order, shape, ns, _r) in rows:
        print(f"  {kind:8s} {s[2]:.3f} {s[3]:.3f} | " + ' '.join(f"{x:+9.4f}" for x in q)
              + f" | {order:5.2f} | {shape} | {ns}")
    # candidate laws
    print("\n# candidate: delta*/eps^2 (smallest eps) against simple invariants")
    Q = np.array([x[2][0] for x in rows])
    for name, fn in (('s3^2', lambda s, v: s[2] ** 2), ('s4^2', lambda s, v: s[3] ** 2),
                     ('mean v^2', lambda s, v: float(np.mean(v ** 2))),
                     ('s1^2+..', lambda s, v: float(np.sum(s[:3] ** 2)))):
        I = np.array([fn(x[1], np.array(x[6]['v'])) for x in rows])
        ok = I > 1e-9
        ratio = Q[ok] / I[ok]
        print(f"    -Q / {name:9s}: min {(-ratio).min():8.4f}   median {np.median(-ratio):8.4f}   "
              f"max {(-ratio).max():8.4f}   (n={ok.sum()})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('run')
    p.add_argument('--n', type=int, default=6)
    p.add_argument('--T', type=float, default=3.0)
    p.add_argument('--eps', default='0.5,1,2,4')
    p.add_argument('--ndense', type=int, default=50)
    p.add_argument('--nsparse', type=int, default=10)
    p.add_argument('--limit', type=int, default=3000)
    p.add_argument('--jitters', type=int, default=2)
    p.add_argument('--nproc', type=int, default=14)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser('report')
    p.add_argument('file')
    p.set_defaults(fn=cmd_report)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
