#!/usr/bin/env python3
"""GS1 (10-16): slide a piece along its H(W) curve (warm-started fixed-d polishes) and re-fill the pair; do the
MX1 one-short pairs (inter-box gaps off-integer) recover?  Pre-registered in PACKER.md.

  python3 analysis/gs1.py --out runs/gs1 [--span 0.4 --step 0.02 --cores 14]"""
import argparse, json, os, sys, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path[:0] = [PK, HERE]
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import compose3 as C
LIB = 'runs/mondrian/croplib_live1016.npz,store'
_ST = None


def store_piece(lab):
    global _ST
    if _ST is None:
        from pk.store import Store
        _ST = {f"piece:{p.meta['source']}#{p.meta['id']}": p for p in Store().pieces(mmin=3)}
    return _ST.get(lab)


def fq_quench(fq, h, sq, d):
    from pk.engine import QOpt
    from pk.packing import Packing
    r = fq.quench(Packing(float(h), sq), QOpt(loosen=1.0, rect_dx=float(d), flip_top=0), timeout=300)
    return (r.s, np.array(r.p.sq, float)) if r.ok else (None, None)


def trace(pid, span, step):
    """traced shapes of table row pid: [(d, h, T in the row's LL frame)]"""
    import mondrian as M
    from pk.engine import FQ
    sp = store_piece(C._P[pid][4])
    h0, d0, sq0 = float(sp.s), float(sp.dx), np.array(sp.sq, float)
    fq = FQ(); pts = [(d0, h0, sq0)]
    for sgn in (-1, 1):
        h, sq, d = h0, sq0, d0
        for _ in range(int(round(span / step))):
            dn = d + sgn * step
            if h + dn <= 1: break
            s2 = sq.copy(); s2[:, 0] *= (h + dn) / (h + d)
            hn, sn = fq_quench(fq, h * 1.002, s2, dn)
            if hn is None: break
            pts.append((dn, hn, sn)); h, sq, d = hn, sn, dn
    flip = abs(C._P[pid][1] - (h0 + d0)) > 1e-6
    out = []
    for d, h, sq in pts:
        T = sq[M.is_tilted(sq[:, 2])].copy()
        if flip:
            T = np.c_[T[:, 1], T[:, 0], (-T[:, 2]) % 90.0]
        out.append((d, h, T))
    return out


def polish(n, S, A, B):
    import mondrian as M, fillfast
    from pk.engine import FQ, QOpt, POLISH
    from pk.packing import Packing
    L = fillfast.fill(S, np.r_[A, M._drop_conflicts(A, B)])
    eng = FQ()
    try:
        r = eng.quench(Packing(S, M.trim(L, n)), QOpt(), timeout=180)
        r2 = eng.quench(r.p, POLISH, timeout=180) if r.ok else r
    finally:
        eng.close()
    return float(r2.s) if r2.ok else None


def job(args):
    pid, cases, span, step = args          # cases: [(n, S, need, i, j, corner)] all containing pid
    import mondrian as M
    t0 = time.process_time()
    pts = trace(pid, span, step)
    res = []
    for n, S, need, i, j, c in cases:
        fl2 = C.fl(S) ** 2
        base = C.count(S, C.compose(S, i, j, c)) - fl2
        best, bd, hit = -999, None, None
        for d, h, T in pts:
            A = M.to_corner(T if i == pid else C._P[i][0], 'LL', S)
            B = M.to_corner(T if j == pid else C._P[j][0], c, S)
            if not (M.fits(A, S) and M.fits(B, S)):
                continue
            k = C.count(S, np.r_[A, M._drop_conflicts(A, B)]) - fl2
            if k > best: best, bd = k, d
            if k >= need and hit is None:
                hit = dict(d=d, s=polish(n, S, A, B))
        res.append(dict(n=n, S=S, need=need, a=C._P[i][4], b=C._P[j][4], corner=c, traced=C._P[pid][4],
                        d0=pts[0][0], base=int(base), best=int(best), best_d=bd, hit=hit, npts=len(pts)))
    return res, time.process_time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--pairs', default='runs/mx1/oneshort.json')
    ap.add_argument('--n-max', type=int, default=40); ap.add_argument('--per-n', type=int, default=8)
    ap.add_argument('--span', type=float, default=0.4); ap.add_argument('--step', type=float, default=0.02)
    ap.add_argument('--cores', type=int, default=14)
    a = ap.parse_args()
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'gs1', args=vars(a), procs=a.cores)
    C._init(LIB)
    rows = json.load(open(a.pairs))['rows']
    by = collections.OrderedDict()                       # the same selection as mx1.py
    for n, S, need, i, j, c in rows:
        L = by.setdefault(int(n), [])
        if len(L) < a.per_n and (int(i), int(j), c) not in [(x[3], x[4], x[5]) for x in L]:
            L.append((int(n), float(S), int(need), int(i), int(j), c))
    ns = sorted(by)
    pick = [ns[round(k * (len(ns) - 1) / (a.n_max - 1))] for k in range(a.n_max)] if len(ns) > a.n_max else ns
    pairs = [x for n in sorted(set(pick)) for x in by[n]]
    per = collections.defaultdict(list)
    for x in pairs:
        for p in {x[3], x[4]}:
            if C._P[p][4].startswith('piece:') and store_piece(C._P[p][4]) is not None:
                per[p].append(x)
    jobs = [(p, cs, a.span, a.step) for p, cs in per.items()]
    print(f'{len(pairs)} pairs, {len(jobs)} traceable pieces, {sum(len(j[1]) for j in jobs)} (piece, pair) cases', flush=True)
    out = open(f'{a.out}/gs1.jsonl', 'w'); cpu = 0
    with ProcessPoolExecutor(a.cores, initializer=C._init, initargs=(LIB,)) as ex:
        for res, c in ex.map(job, jobs):
            cpu += c
            for r in res:
                out.write(json.dumps(r) + '\n'); out.flush()
                print(f"n={r['n']} need {r['need']} base {r['base']} traced best {r['best']} (d {r['best_d']} vs {r['d0']})"
                      f"{'  REACH ' + json.dumps(r['hit']) if r['hit'] else ''}", flush=True)
    man.finish(summary=dict(cpu_h=cpu / 3600))


if __name__ == '__main__':
    main()
