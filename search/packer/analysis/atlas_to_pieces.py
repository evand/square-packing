#!/usr/bin/env python3
"""Import Mondrian atlas crops (corner crops of known packings, any aspect W x H) into the store's piece table as
pieces with d = W - H (non-integer in general), so compose2 can use off-diagonal pieces (10-15: L2-481 found the piece
table held only (h + d) x h pieces with integer d: 0 off-diagonal of 45k).  Source label: '<atlas source label>:atlas'
(store:<id>:... keeps piece_parent working for --blind).

  python3 analysis/atlas_to_pieces.py runs/mondrian/atlas_v3.npz [--mmin 10 --mmax 250]
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('atlas'); ap.add_argument('--mmin', type=int, default=10)
    ap.add_argument('--mmax', type=int, default=250); a = ap.parse_args()
    import mondrian as M
    from pk.packing import Packing
    from pk.store import Store
    z = np.load(a.atlas, allow_pickle=True)
    W, H, m, o, sid, labels, src_s, SQ = z['W'], z['H'], z['m'], z['o'], z['sid'], z['labels'], z['src_s'], z['sq']
    S = Store(); new = seen = bad = 0
    for k in np.nonzero((m >= a.mmin) & (m <= a.mmax))[0]:
        q = M.orient(np.asarray(SQ[sid[k]], float), float(src_s[sid[k]]), int(o[k]))
        h = M.half_ext(q[:, 2])
        ok = (q[:, 0] + h <= W[k] + 1e-9) & (q[:, 1] + h <= H[k] + 1e-9)
        if ok.sum() != m[k]:
            bad += 1; continue
        lab = str(labels[sid[k]]); parts = lab.split(':')
        src = (':'.join(parts[:2]) if parts[0] == 'store' else lab.split('@')[0]) + ':atlas'
        _, c = S.add_piece(Packing(float(H[k]), q[ok], dx=float(W[k] - H[k])), src, ref=os.path.basename(a.atlas), commit=False)
        new += c; seen += 1
    S.db.commit()
    print(f'seen {seen}, new {new}, count mismatch {bad}')


if __name__ == '__main__':
    main()
