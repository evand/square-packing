#!/usr/bin/env python3
"""s21_cover_eval.py -- HEURISTIC evaluation of a certificate-format cover of [0,s]^2 (search/S21_COVER.md).

Three float checks, none of them a proof (the exact word is search/zeromargin.py):

  hardscan CERT [--nproc 7] [--pitch 0.004]
      the S21_OVERHEAD.md sec 2 protocol: lattice scan (closed4.scan_angle, pitch 0.004 + wall bands) at
      theta = 0, 0.02..3 deg step 0.02 (near-tile), 3..45 deg step 0.25, then closed4.polish from the
      1400 worst poses.  Stricter than `closed4.py stress` (which samples only integer degrees and a few
      tiny angles, and misses the 0.1-0.5 deg and non-integer tilted dips).  Prints total / min = cost.
      Calibration: certificates/rung2/s13_closed_cover_4.txt -> min 1.0200245 = 1.05 x 0.97145,
      reproducing S21_OVERHEAD.md's 0.9715 for the r5 LP weights.
  family CERT [--pitch 0.0005]
      the ZEROMARGIN.md sec 4 item-3 family (family_rows.item3: c_x or c_y on a half-integer line, the
      other coordinate on a fine pitch, 11 small angles 0 .. 10^-1.5 rad), which the lattice steps over.
      Calibration: runs/closed4_best.txt -> 0.9420217, the known exact violation (RUNG2.md sec 4.3).
  regions CERT ...
      where the weight sits: corner cells / wall bands / interior, per unit cell (a point on an
      interior grid line split equally between the cells it bounds, as s21_overhead.py regions),
      per distance to the wall, on / near / off the integer grid lines.

The covers are D4-symmetric (closed4.Model orbits), so angles in [0, 45] suffice for hardscan.
"""
import sys, os, math, time, argparse
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import closed4 as C
import family_rows as F


def _scan(args):
    P, w, s, th, pitch = args
    v, p = C.scan_angle(P, w, s, th, pitch)
    o = np.argsort(v)[:30]
    return v[o], p[o]


def _pol(args):
    P, w, s, seeds, seed = args
    fp, fv, cur, curv = C.polish(P, w, s, seeds, rounds=9, nper=64, rng=np.random.default_rng(seed))
    return cur, curv


def _cap(args):
    P, w, Q = args
    return C.captured(P, w, Q)


def hardscan(path, nproc=7, pitch=0.004):
    s, P, w = C.load_points(path); t0 = time.time()
    degs = [0.0] + list(np.arange(0.02, 3.0, 0.02)) + list(np.arange(3.0, 45.0 + 1e-9, 0.25))
    with Pool(nproc) as pool:
        outs = pool.map(_scan, [(P, w, s, math.radians(d), pitch) for d in degs])
        V = np.concatenate([o[0] for o in outs]); Q = np.concatenate([o[1] for o in outs])
        lat = float(V.min()); seeds = Q[np.argsort(V)[:1400]]
        res = pool.map(_pol, [(P, w, s, ch, i) for i, ch in enumerate(np.array_split(seeds, 4 * nproc))])
    cur = np.concatenate([r[0] for r in res]); curv = np.concatenate([r[1] for r in res])
    j = int(curv.argmin()); m = min(lat, float(curv[j]))
    print(f"{path}: total {w.sum():.6f}; hardscan {len(degs)} angles pitch {pitch}: lattice min {lat:.7f}; "
          f"polished min {curv[j]:.7f} at ({cur[j,0]:.5f}, {cur[j,1]:.5f}, {math.degrees(cur[j,2]):.4f} deg); "
          f"cost total/min = {w.sum()/m:.6f}; {time.time()-t0:.0f}s")
    for k in np.argsort(curv)[:6]:
        print(f"   {curv[k]:.7f} at ({cur[k,0]:.5f}, {cur[k,1]:.5f}, {math.degrees(cur[k,2]):.4f} deg)")
    return m


def family(path, pitch=0.0005, nproc=7):
    s, P, w = C.load_points(path)
    angles = [0.0] + [10 ** e for e in (-7, -6, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5)]
    Q = np.array(F.item3(s, pitch, angles))
    with Pool(nproc) as pool:
        v = np.concatenate(pool.map(_cap, [(P, w, ch) for ch in np.array_split(Q, 8 * nproc)]))
    j = int(v.argmin())
    print(f"{path}: total {w.sum():.6f}; item-3 family ({len(Q)} poses, pitch {pitch}): min {v[j]:.7f} at "
          f"({Q[j,0]:.5f}, {Q[j,1]:.5f}, {math.degrees(Q[j,2]):.6f} deg); cost {w.sum()/v[j]:.6f}; "
          f"poses < 1: {(v < 1 - 1e-9).sum()}")
    return float(v[j])


def regions(path):
    sfr, D, WD, rows = C.read_cert(path)
    K = sfr.numerator * D // sfr.denominator; n = int(sfr)
    X, Y, W = rows[:, 0], rows[:, 1], rows[:, 2].astype(float) / WD
    tot = W.sum()
    dx = np.minimum(X, K - X); dy = np.minimum(Y, K - Y); d = np.minimum(dx, dy); e = np.maximum(dx, dy)
    print(f"== {path}: s={float(sfr)} points={len(W)} total={tot:.6f}")
    for name, msk in (('closed corner squares  max(dx,dy)<=1', e <= D), ('wall bands  min(dx,dy)<=1<max', (d <= D) & (e > D)),
                      ('interior  min(dx,dy)>1', d > D)):
        print(f"   {name:38s} {W[msk].sum():9.4f}  {100*W[msk].sum()/tot:5.1f}%")
    print("   by distance d to the nearest wall:")
    for name, msk in (('0<d<1', (d > 0) & (d < D)), ('d=1', d == D), ('1<d<2', (d > D) & (d < 2 * D)), ('d=2', d == 2 * D),
                      ('d>2', d > 2 * D), ('d=0 (on a wall)', d == 0)):
        print(f"      {name:18s} {W[msk].sum():9.4f}  {100*W[msk].sum()/tot:5.1f}%")
    online = ((X % D == 0) & (X > 0) & (X < K)) | ((Y % D == 0) & (Y > 0) & (Y < K))
    near = ~online & ((np.abs(((X + D // 2) % D) - D // 2) <= D // 20) | (np.abs(((Y + D // 2) % D) - D // 2) <= D // 20))
    print(f"   on integer grid lines {W[online].sum():.4f} ({100*W[online].sum()/tot:.1f}%), within 0.05 of one "
          f"{W[near].sum():.4f}, elsewhere {W[~online & ~near].sum():.4f}")
    print("   per line: " + "  ".join(f"x={k}: {W[(X == k * D) & (Y > 0) & (Y < K)].sum():.4f}" for k in range(1, n)))
    T = np.zeros((n, n))
    def split(Z):
        return [[(q - 1, .5), (q, .5)] if (r == 0 and 0 < q < n) else [(min(q, n - 1), 1.0)] for q, r in zip(Z // D, Z % D)]
    for sx, sy, wt in zip(split(X), split(Y), W):
        for i, a in sx:
            for j, b in sy: T[i, j] += wt * a * b
    print("   unit-cell weights (grid-line points split equally between the cells they bound):")
    for j in range(n - 1, -1, -1): print("      " + " ".join(f"{T[i, j]:6.3f}" for i in range(n)))
    ne = np.array([[(min(i, n - 1 - i) == 0) + (min(j, n - 1 - j) == 0) for j in range(n)] for i in range(n)])
    cm, em, im = T[ne == 2].mean(), T[ne == 1].mean(), T[ne == 0].mean()
    print(f"   class means: corner {cm:.4f}  edge {em:.4f}  interior {im:.4f}   (4a + 4(m-2)b + (m-2)^2 c = "
          f"{4*cm + 4*(n-2)*em + (n-2)**2*im:.4f})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['hardscan', 'family', 'regions'])
    ap.add_argument('certs', nargs='+')
    ap.add_argument('--nproc', type=int, default=7)
    ap.add_argument('--pitch', type=float, default=None)
    a = ap.parse_args()
    for c in a.certs:
        if a.mode == 'hardscan': hardscan(c, a.nproc, a.pitch or 0.004)
        elif a.mode == 'family': family(c, a.pitch or 0.0005, a.nproc)
        else: regions(c)


if __name__ == '__main__':
    main()
