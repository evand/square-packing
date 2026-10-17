#!/usr/bin/env python3
"""MX1 (10-16): is the pair interface where composition loses squares, and does keeping both pieces whole
("meshing": overlaps left for the quench to resolve) beat compose3's drop-conflicts?  Pre-registered in PACKER.md.

  python3 analysis/mx1.py --pairs ONESHORT.json --out runs/mx1 [--n-max 40 --per-n 8 --cores 14]"""
import argparse, json, math, os, sys, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path[:0] = [PK, HERE]
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import compose3 as C

LIB = 'runs/mondrian/croplib_live1016.npz,store'


def top_up(L, n, S, seed):
    import mondrian as M
    from rectsearch import holes
    rng = np.random.default_rng(seed)
    L = np.array(L, float)
    while len(L) < n:
        x, y = holes(L, S, S, rng)
        L = np.r_[L, [[x, y, 0.0]]]
    return M.trim(L, n)


def quench(n, S, L):
    from pk.engine import FQ, QOpt, POLISH
    from pk.packing import Packing
    eng = FQ()
    try:
        r = eng.quench(Packing(S, L), QOpt(), timeout=180)
        r2 = eng.quench(r.p, POLISH, timeout=180) if r.ok else r
    finally:
        eng.close()
    return float(r2.s) if r2.ok else None


def job(args):
    n, S, need, a, b, corner = args
    import mondrian as M, fillfast
    t0 = time.process_time()
    P = C._P
    A, B = M.to_corner(P[a][0], 'LL', S), M.to_corner(P[b][0], corner, S)
    Bk = M._drop_conflicts(A, B)
    fl2 = C.fl(S) ** 2
    eA = C.count(S, P[a][0]) - fl2; eB = C.count(S, P[b][0]) - fl2
    LR = fillfast.fill(S, np.r_[A, Bk]); kR = len(LR)
    LM = fillfast.fill(S, np.r_[A, B]); kM = len(LM)
    sR = quench(n, S, top_up(LR, n, S, 1))
    sM = quench(n, S, top_up(LM, n, S, 1))
    return dict(n=n, S=S, need=need, a=P[a][4], b=P[b][4], corner=corner, eA=int(eA), eB=int(eB), jointR=int(kR - fl2),
                jointM=int(kM - fl2), dropped=int(len(B) - len(Bk)), tiltA=len(A), tiltB=len(B), sR=sR, sM=sM,
                cpu=round(time.process_time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pairs', required=True); ap.add_argument('--out', required=True)
    ap.add_argument('--n-max', type=int, default=40); ap.add_argument('--per-n', type=int, default=8)
    ap.add_argument('--cores', type=int, default=14)
    a = ap.parse_args()
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'mx1', args=vars(a), procs=a.cores)
    rows = json.load(open(a.pairs))['rows']
    by = collections.OrderedDict()
    for n, S, need, i, j, c in rows:
        L = by.setdefault(int(n), [])
        if len(L) < a.per_n and (int(i), int(j), c) not in [(x[3], x[4], x[5]) for x in L]:
            L.append((int(n), float(S), int(need), int(i), int(j), c))
    ns = sorted(by)
    pick = [ns[round(k * (len(ns) - 1) / (a.n_max - 1))] for k in range(a.n_max)] if len(ns) > a.n_max else ns
    jobs = [x for n in sorted(set(pick)) for x in by[n]]
    print(f'{len(jobs)} pairs over {len(set(pick))} n', flush=True)
    out = open(f'{a.out}/mx1.jsonl', 'w')
    with ProcessPoolExecutor(a.cores, initializer=C._init, initargs=(LIB,)) as ex:
        for r in ex.map(job, jobs):
            out.write(json.dumps(r) + '\n'); out.flush()
            print(f"n={r['n']} need {r['need']} eA+eB {r['eA'] + r['eB']} jointR {r['jointR']} jointM {r['jointM']} "
                  f"dropped {r['dropped']}  sR {r['sR']}  sM {r['sM']}", flush=True)
    man.finish()


if __name__ == '__main__':
    main()
