#!/usr/bin/env python3
"""Phase 1 of the Final5 generator (2026-10-21): collect all polynomial leaves needed by the case theorems
(gen5.emit_case), dedupe, certify (cert5.certify_split) in parallel.  Cache: runs/rectfront/r5cert/plans.pkl."""
import itertools, os, pickle, sys, time
from multiprocessing import Pool
import sympy as sp
from gen5 import emit_case
from cert5 import certify_split

CACHE = '/home/evand/math/square-packing/public/runs/rectfront/r5cert/plans.pkl'
CASES = [(e1, e2, b1, b2, c) for e1, e2, b1, b2 in itertools.product((False, True), (False, True), 'UP', 'UP')
         for c in ('up0', 'up90', 'a1', 'a2')]

def collect():
    need = {}
    for case in CASES:
        certs = {}
        emit_case(case, certs, 'f5p')
        for lname, P, Gh, c, sub, Hint in certs['need']:
            key = (sp.srepr(P), tuple((g[0], sp.srepr(g[1])) for g in Gh))
            need.setdefault(key, (P, Gh, []))[2].append(lname)
    return need

def work(item):
    key, (P, Gh, names) = item
    t = time.time()
    plan = certify_split(P, Gh)
    return key, plan, names, time.time() - t

if __name__ == '__main__':
    need = collect()
    print('distinct leaves:', len(need), ' total uses:', sum(len(v[2]) for v in need.values()), flush=True)
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    todo = [(k, v) for k, v in need.items() if k not in cache]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 10) as pool:
        for i, (key, plan, names, dt) in enumerate(pool.imap_unordered(work, todo)):
            cache[key] = plan
            pickle.dump(cache, open(CACHE, 'wb'))
            print(f'{i + 1}/{len(todo)} {names[0]} {"OK " + str(len(plan)) + " box(es)" if plan else "FAIL"} {dt:.0f}s', flush=True)
