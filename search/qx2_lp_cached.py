#!/usr/bin/env python3
"""qx2_lp_cached.py -- run qx2_lp.py with the exact germ-limit rows cached on disk (2026-10-02).

GermRows depends only on the element model (R, w, delta, hc, hp, points), not on kappa or the rows, and costs 30-60 min
at R >= 4 on 2 cores; every warm restart after a HiGHS failure recomputed it.  This wrapper replaces qx2_lp.GermRows by
a subclass that loads/saves (A, rhs) under runs/qx2_germcache/<key>.npz.  Same command line as qx2_lp.py.
"""
import os, sys, hashlib, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
import qx2_lp as X

_Orig = X.GermRows
CACHE = os.path.join(X.Q.RUNS, 'qx2_germcache')


class CachedGermRows(_Orig):
    def __init__(self, m, nproc):
        key = json.dumps([m.R, m.w, m.a, m.nvar, sorted((k, len(v)) for k, v in m.E.items()),
                          [list(map(float, l[1:])) if isinstance(l, tuple) else l for l in m.label]], default=str)
        h = hashlib.sha256(key.encode()).hexdigest()[:16]
        os.makedirs(CACHE, exist_ok=True); f = os.path.join(CACHE, h + '.npz')
        if os.path.exists(f):
            z = np.load(f)
            self.A = sp.csr_matrix((z['data'], z['indices'], z['indptr']), shape=tuple(z['shape']))
            self.rhs = z['rhs']; print(f"germ rows: loaded cache {f}", flush=True)
        else:
            super().__init__(m, nproc)
            A = self.A.tocsr()
            np.savez(f, data=A.data, indices=A.indices, indptr=A.indptr, shape=np.array(A.shape), rhs=self.rhs)
            print(f"germ rows: saved cache {f}", flush=True)
        self.added = np.zeros(self.A.shape[0], bool)


X.GermRows = CachedGermRows

if __name__ == '__main__':
    X.run(X.parser().parse_args())
