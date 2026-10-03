#!/usr/bin/env python3
"""Re-run QXChecker (the certificate's checker, same args as run V2) on a few root boxes near U's boundary with
dump=True, and stress every PIECE leaf whose bound uses the polygon: exact independent mass at sampled admissible
rational poses (must be >= 1 and >= the leaf's bound), polygon part vs area(Q∩U), K's vertices in Q.
Usage: taskset -c 6 python3 test_leaves.py cx0 cy0 iu  (root = [cx0, cx0+1/10] x [cy0, cy0+1/10] x u-bin iu of 8)
"""
import sys, os, time
from fractions import Fraction as F
import test_poly as TP
from test_poly import S, ZM, zm, MC
sys.path.insert(0, S)
import qx2_zm as QX


def main():
    cx0, cy0, iu = F(sys.argv[1]), F(sys.argv[2]), int(sys.argv[3])
    root = (cx0, cx0 + F(1, 10), cy0, cy0 + F(1, 10), F(iu, 16), F(iu + 1, 16))
    path = os.path.join(S, 'qx2_data/L4_k02_box7.txt')
    cv = MC.load(path); cov = ZM.Cover(cv)
    a, b = QX.u_square(cov); m = cov.m
    chk = QX.QXChecker(cov, exact_umax=F(1, 2), exact_from=3, use_exact=True, max_depth=18, use_chain=False,
                       cert_mode=False, dump=True)
    chk0 = ZM.MixedChecker(ZM.Cover(dict(cv, polygons=[])), max_depth=18, use_chain=False, cert_mode=True)
    t0 = time.time()
    st, unc, leaves = chk.run_box(root)
    print('root', [str(v) for v in root], {k: v for k, v in st.items() if v}, 'unc', len(unc),
          '%.0fs' % (time.time() - t0), flush=True)
    npiece = nposes = nbad = 0; minmass = None; minpolygap = None
    for box, kind, wit in leaves:
        if kind != 'PIECE': continue
        B = zm.bin_data(box[4], box[5])
        L, _ = chk.piece_bound(box, B)
        L0, _ = chk0.piece_bound(box, B)
        polyL = L - L0
        if polyL == 0: continue
        npiece += 1
        K, _ = TP.region_K(chk, box, B)
        poses = [p for p in TP.sample_poses(box, m, nu=6, nc=3) if TP.admissible(m, *p)]
        for cx, cy, u in poses:
            nposes += 1
            A = TP.area_QU(cx, cy, u, a, b)
            M = TP.exact_mass_indep(cv, cx, cy, u, a, b)
            if minmass is None or M < minmass[0]: minmass = (M, box, (cx, cy, u))
            g = A - polyL
            if minpolygap is None or g < minpolygap: minpolygap = g
            bad = M < 1 or polyL > A or M < L or any(not TP.in_square(v, cx, cy, u) for v in K)
            if bad:
                nbad += 1
                print('BAD', [str(v) for v in box], (str(cx), str(cy), str(u)), float(M), float(L), float(polyL),
                      float(A), flush=True)
    print(f'PIECE leaves with polygon part > 0: {npiece}; poses {nposes}; bad {nbad}; '
          f'min mass {float(minmass[0]) if minmass else None}; min area-polybound gap '
          f'{float(minpolygap) if minpolygap is not None else None}', flush=True)


if __name__ == '__main__':
    main()
