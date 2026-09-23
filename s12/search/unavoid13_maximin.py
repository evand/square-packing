#!/usr/bin/env python3
"""unavoid13: continuous maximin continuation of a k-point set (upper-bound search, [heuristic]).

    python3 search/unavoid13_maximin.py P.npy [P2.npy ...] --m 4 [--iters 30] [--out runs/unavoid13_maximin]

For each start set P: repeat { W = worst poses of P (grid + Nelder-Mead local minima, all with
max_p depth < +0.03, accumulated over iterations); assign each pose of W to its deepest point;
LP: maximise t s.t. depth(p_a(Q), Q) >= t for all Q in W (trust region 0.05); } and report the
margin  min over all admissible poses of max_p depth  after each iteration (re-measured by a fresh
search).  A margin >= 0 (up to 1e-9) means the set is unavoidable up to float precision and is
written as a certificate for unavoid13_check.py.  This is the brief's "multistart on
max_Q min_p dist(p,Q)" with the IP solutions as the starts.
"""
import sys, os, math, time, argparse, json
from fractions import Fraction as F
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unavoid13_lib as L


def continue_set(P, m, fam, iters=30, want=60, trust=0.05, verbose=True):
    P = np.asarray(P, float).copy()
    sqF = L.Squares(fam) if fam else None
    W = []          # accumulated poses (cx, cy, th)
    hist = []
    best = (-9, None)
    for it in range(iters):
        viol, gmin = L.find_violations(P, m, dth_deg=0.5, pitch=0.01, top=600, want=want, cutoff=0.03, verbose=False)
        margin = viol[0][0] if viol else gmin
        hist.append(float(margin))
        if margin > best[0]: best = (float(margin), P.copy())
        if verbose:
            print(f"    iter {it}: margin {margin:+.5f} (grid min {gmin:+.5f}), |W| = {len(W)}", flush=True)
        if margin >= -1e-9:
            break
        for (f, cx, cy, th) in viol:
            W.append((cx, cy, th))
        # squares object from W (float poses; polish only needs cx, cy, cos, sin)
        sq = L.Squares([(F(math.tan(t / 2)).limit_denominator(10 ** 6), F(cx).limit_denominator(10 ** 6),
                         F(cy).limit_denominator(10 ** 6)) for (cx, cy, t) in W])
        Pn, tstar = L.polish_positions(P, m, ([sqF] if sqF else []) + [sq], rounds=4, delta=0.03, trust=trust, move_pen=0.0, verbose=False)
        moved = np.abs(Pn - P).max()
        P = Pn
        if moved < 1e-10 and it > 2:
            break
    return best[0], best[1], hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('starts', nargs='+')
    ap.add_argument('--m', type=str, default='4')
    ap.add_argument('--iters', type=int, default=30)
    ap.add_argument('--out', default='runs/unavoid13_maximin')
    ap.add_argument('--family', default=None, help='family file whose squares are always enforced in the polish')
    a = ap.parse_args()
    m = F(a.m); os.makedirs(a.out, exist_ok=True)
    fam = L.read_family(a.family)[1] if a.family else []
    results = []
    for path in a.starts:
        P0 = np.load(path)
        print(f"start {path}: {len(P0)} points", flush=True)
        t0 = time.time()
        best, Pb, hist = continue_set(P0, m, fam, iters=a.iters)
        tag = os.path.basename(path).replace('.npy', '')
        np.save(f"{a.out}/{tag}_best.npy", Pb)
        print(f"  -> best margin {best:+.5f} after {len(hist)} iters, {time.time()-t0:.0f}s", flush=True)
        results.append(dict(start=path, best=best, hist=hist))
        if best >= -1e-9:
            Pex = L.snap_points(Pb, 10 ** 4)
            cert = f"{a.out}/{tag}_cert.txt"; L.write_cert(cert, m, Pex)
            print(f"  margin >= 0: certificate candidate written to {cert}", flush=True)
    json.dump(results, open(f"{a.out}/results.json", 'w'), indent=1)


if __name__ == '__main__':
    main()
