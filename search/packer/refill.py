#!/usr/bin/env python3
"""Refill test (10-13): keep a packing's tilted squares, re-fill the uprights exactly (mondrian.exact_fill, MILP) at the
same side.  count > n means s(n + k) <= s(n) for k = count - n: an insertion record wherever the frontier at n + k is
higher.  Inputs: the store's best packing per n.

  python3 refill.py --ns 20-340 --out runs/refill/r1 [--procs 14]
"""
from __future__ import annotations
import argparse, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def job(n):
    import mondrian
    from pk.store import Store
    from pk import eventbench as E
    p = Store().best(n)
    if p is None:
        return dict(n=n)
    T = p.sq[mondrian.is_tilted(p.sq[:, 2])]
    S = float(p.s)
    L = mondrian.exact_fill(T, S, time_limit=60)
    c = len(L)
    out = dict(n=n, s=S, tilted=len(T), count=c)
    if c > n:
        gains = {k: E.frontier_now(n + k) - S for k in range(1, c - n + 1)}
        out['gains'] = {k: g for k, g in gains.items() if g > 1e-9}
        if out['gains']:
            out['sq'] = L.tolist()
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', default='20-340'); ap.add_argument('--out', required=True)
    ap.add_argument('--procs', type=int, default=14)
    a = ap.parse_args()
    from pk import manifest
    import scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'refill', args=vars(a), procs=a.procs)
    scratch.scope('refill')
    lo, hi = map(int, a.ns.split('-'))
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'results.jsonl'), 'w') as f:
        for r in ex.map(job, range(lo, hi + 1)):
            f.write(json.dumps(r) + '\n'); f.flush()
            if r.get('count', 0) > r['n']:
                print(f"{r['n']} count {r['count']} at {r['s']:.10f}  gains {r.get('gains')}", flush=True)
    man.finish()
