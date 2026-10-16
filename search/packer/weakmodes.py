#!/usr/bin/env python3
"""Run exactsolve on a weak minimum and capture its second-order data (P exact point, A active set, lambda, V, weak
contacts, the zero modes of the reduced Hessian) for the higher-order analysis (weakslice.py).  Also describes the modes:
which squares move, rotation share, upright / tilted.

  python3 weakmodes.py DIR        (DIR/in.txt; writes DIR/modes.pkl, prints a one-line JSON summary)
"""
import contextlib, io, json, math, os, pickle, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'exact')); sys.path.insert(0, HERE)
import numpy as np
import exactsolve as E

CAP = {}
_so = E.second_order


def _cap(P, A, lamA, V, squares, weak=()):
    r = _so(P, A, lamA, V, squares, weak)
    CAP.clear(); CAP.update(P=P.copy(), A=list(A), lam=[float(x) for x in lamA], V=list(V), squares=list(squares),
                            weak=list(weak), modes=r.get('_zero_modes'), status=r['status'], eig=r.get('eig'))
    return r


E.second_order = _cap


def main(d):
    with contextlib.redirect_stdout(io.StringIO()):
        rep = E.run(f'{d}/in.txt', outdir=f'{d}/modes.ex', quiet=True)[0]
    P = CAP['P']; V = CAP['V']; M = CAP['modes']
    with open(f'{d}/modes.pkl', 'wb') as f:
        pickle.dump(dict(CAP, S=str(P.S), n=P.n, X=[str(x) for x in P.X], Y=[str(y) for y in P.Y],
                         T=[str(t) for t in P.T]), f)
    idx = {v: k for k, v in enumerate(V)}
    out = dict(dir=os.path.basename(d), status=CAP['status'][:40], nA=len(CAP['A']), nweak=len(CAP['weak']),
               lam_min=min(CAP['lam']) if CAP['lam'] else None)
    if M is not None and M.shape[1]:
        rot = np.array([idx.get(3 * i + 2) for i in range(P.n)], dtype=object)
        movers, rshare = set(), []
        for k in range(M.shape[1]):
            m = M[:, k]
            th = np.array([m[idx[3 * i + 2]] if (3 * i + 2) in idx else 0.0 for i in range(P.n)])
            rshare.append(float(np.linalg.norm(th) / np.linalg.norm(m)))
            for i in range(P.n):
                if any((3 * i + c) in idx and abs(m[idx[3 * i + c]]) > 1e-6 for c in range(3)):
                    movers.add(i)
        tilt = lambda i: abs(math.remainder(math.degrees(float(P.T[i])), 90.0)) > 1e-6
        out.update(n_modes=M.shape[1], movers=len(movers), movers_tilted=sum(tilt(i) for i in movers),
                   rot_share_max=max(rshare), rot_share_med=float(np.median(rshare)))
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
