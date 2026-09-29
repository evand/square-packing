#!/usr/bin/env python3
"""qx2_family_check.py -- independent exact check of a qx2 family file against its box cover (task
quadrant-exact-w2).  Reads ONLY the human-readable family file (qx2_exact.py family) and a FORMAT.md v1 box cover;
shares no code with qx2_lp / qx2_exact / quadrant_lp.

Checks (all Fraction arithmetic):
  * sigma = 0: one period of the profile carries mass w (atoms + the Lebesgue layer [a, w]);
  * mirror symmetry of the profile (phase p -> 1 - p) and diagonal symmetry of the corner module, as multisets;
  * D = R^2 - nu(C) - m_v, with nu(C) including the Lebesgue piece [a, R]^2 and m_v = mass of the phase-0 vertical
    pieces; 4 D > 3;
  * the box measure mu_k built here from the family (corner module at the 4 corners, profile on [R, k-R] along the
    4 walls, Lebesgue on [a, k-a]^2) equals the box cover file as a multiset of (segment, mass) and polygon; total
    = k^2 - 4 D.
usage: python3 qx2_family_check.py FAMILY.txt BOX.txt [--k 7 --R 2 --w 2]
"""
import sys, argparse
from fractions import Fraction as F
from collections import Counter


def read_family(path):
    prof, corner = [], []
    a = None
    for ln in open(path):
        if ln.startswith('#'):
            if 'Lebesgue on U = [a, inf)^2 with a =' in ln: a = F(ln.split('a =')[1].split(';')[0].strip())
            continue
        t = ln.split()
        if not t: continue
        if t[0] in ('h', 'v'): prof.append((t[0],) + tuple(F(v) for v in t[1:]))
        elif t[0] in ('H', 'V'): corner.append((t[0],) + tuple(F(v) for v in t[1:]))
        else: raise ValueError(ln)
    return a, prof, corner


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('family'); ap.add_argument('box')
    ap.add_argument('--k', type=int, default=7); ap.add_argument('--R', type=int, default=2); ap.add_argument('--w', default='2')
    A = ap.parse_args()
    k, R, w = A.k, A.R, F(A.w)
    a, prof, corner = read_family(A.family)
    ok = True
    # ---- sigma
    per = sum(e[-1] for e in prof) + (w - a)
    print(f"a = {a}; one period of the profile: atoms {sum(e[-1] for e in prof)} + layer {w - a} = {per}  "
          f"(sigma = {w - per})")
    ok &= (per == w)
    # ---- mirror symmetry of the profile
    def mir(e):
        if e[0] == 'h': return ('h', e[1], 1 - e[3], 1 - e[2], e[4])
        return ('v', (1 - e[1]) % 1, e[2], e[3], e[4])
    sym_p = Counter(prof) == Counter(mir(e) for e in prof)
    # ---- diagonal symmetry of the corner module
    def diag(e): return ('V' if e[0] == 'H' else 'H',) + e[1:]
    sym_c = Counter(corner) == Counter(diag(e) for e in corner)
    print(f"profile mirror-symmetric: {sym_p}; corner module diagonal-symmetric: {sym_c}")
    ok &= sym_p and sym_c
    # ---- supports inside the right regions; nothing in the open layer / open U
    for e in prof:
        if e[0] == 'h': assert 0 < e[1] <= a and 0 <= e[2] < e[3] <= 1, e
        else: assert 0 <= e[1] < 1 and 0 <= e[2] < e[3] <= a, e
    for e in corner:
        p, t0, t1 = e[1], e[2], e[3]
        assert 0 <= p <= R and 0 <= t0 < t1 <= R, e
        assert not (p > a and t0 >= a), ("atom inside U", e)
    # ---- D
    nuC = sum(e[-1] for e in corner) + (R - a) ** 2
    mv = sum(e[-1] for e in prof if e[0] == 'v' and e[1] == 0)
    D = R * R - nuC - mv
    print(f"nu(C) = {nuC} = {float(nuC):.12f}; m_v = {mv} = {float(mv):.12f}; D = {D} = {float(D):.12f}; 4D > 3: {4 * D > 3}")
    ok &= 4 * D > 3
    # ---- the box measure from the family
    segs = Counter()
    K = F(k)
    maps = [lambda x, y: (x, y), lambda x, y: (K - x, y), lambda x, y: (x, K - y), lambda x, y: (K - x, K - y)]
    for e in corner:
        o, p, t0, t1, ms = e
        for f in maps:
            P, Q = (f(t0, p), f(t1, p)) if o == 'H' else (f(p, t0), f(p, t1))
            segs[(tuple(sorted([P, Q])), ms)] += 1
    rots = [lambda x, y: (x, y), lambda x, y: (K - y, x), lambda x, y: (K - x, K - y), lambda x, y: (y, K - x)]
    for e in prof:
        for j in range(R, k - R + 1):
            if e[0] == 'h':
                _, y, p0, p1, ms = e
                if j + p1 > k - R: continue
                for f in rots: segs[(tuple(sorted([f(j + p0, y), f(j + p1, y)])), ms)] += 1
            else:
                _, p, y0, y1, ms = e
                if j + p > k - R: continue
                for f in rots: segs[(tuple(sorted([f(j + p, y0), f(j + p, y1)])), ms)] += 1
    # ---- the box cover file
    tok = []
    for ln in open(A.box):
        ln = ln.split('#', 1)[0]; tok += ln.split()
    it = iter(tok)
    assert next(it) == 'mixed' and next(it) == '1'
    sn, sd = int(next(it)), int(next(it)); Dd = int(next(it)); Wd = int(next(it))
    assert F(sn, sd) == K
    npnt = int(next(it)); assert npnt == 0, "points present"
    ns = int(next(it))
    fsegs = Counter()
    for _ in range(ns):
        X0, Y0, X1, Y1, wm = (int(next(it)) for _ in range(5))
        fsegs[(tuple(sorted([(F(X0, Dd), F(Y0, Dd)), (F(X1, Dd), F(Y1, Dd))])), F(wm, Wd))] += 1
    npg = int(next(it)); assert npg == 1
    kk = int(next(it)); wm = F(int(next(it)), Wd); V = [(F(int(next(it)), Dd), F(int(next(it)), Dd)) for _ in range(kk)]
    same = (segs == fsegs)
    lebok = (sorted(V) == sorted([(a, a), (K - a, a), (K - a, K - a), (a, K - a)]) and wm == (K - 2 * a) ** 2)
    total = sum(ms * c for (_, ms), c in fsegs.items()) + wm
    print(f"box k = {k}: {sum(segs.values())} segments built from the family; file has {sum(fsegs.values())}; "
          f"identical multisets: {same}; Lebesgue square [{a},{K-a}]^2 with density 1: {lebok}")
    print(f"total = {total} = {float(total):.12f}; k^2 - 4D = {K*K - 4*D}; equal: {total == K*K - 4*D}; < k^2 - 3: {total < K*K - 3}")
    ok &= same and lebok and total == K * K - 4 * D and total < K * K - 3
    print("FAMILY CHECK", "OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
