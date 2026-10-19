#!/usr/bin/env python3
"""Funnel width, free vs Md-symmetric soft-to-rigid (10-18).  Start = the record, loosened x 1.05 about the centre,
then a Gaussian kick sigma on x, y (20 sigma degrees on angles); the sym arm's kick is symmetrised (projected onto the
Md subspace: pairs averaged with their images, diagonal squares on the diagonal at their own angle).  Both: harden
reharden (r0 0.1, box free while soft; sym arm + --sym-md 1), fq polish.  Return = polished side within 1e-7 of the
record.
  symfunnel.py --ns 27,52,65,82,101,125 --sig 0.05,0.1,0.2,0.4 --k 12 --out runs/symsub/funnel1
"""
import sys, os, json, argparse, subprocess, tempfile, time, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symscan import invariant
from symharden import write_cfg, HARDEN
from scipy.spatial import cKDTree


def partners(sq):
    tr = cKDTree(sq[:, :2]); _, k = tr.query(sq[:, [1, 0]]); return k


def sym_kick(sq, sig, rng):
    k = partners(sq); n = len(sq)
    d = np.c_[rng.normal(0, sig, (n, 2)), rng.normal(0, 20 * sig, n)]
    out = sq + d
    for i in range(n):
        j = k[i]
        if j == i:
            m = 0.5 * (out[i, 0] + out[i, 1]); out[i, :2] = m; out[i, 2] = sq[i, 2]
        elif i < j:
            a = 0.5 * (out[i, 0] + out[j, 1]); b = 0.5 * (out[i, 1] + out[j, 0])
            th = out[i, 2]
            out[i, :2] = a, b; out[j, :2] = b, a; out[j, 2] = (-th) % 90
    return out


def job(t):
    from pk.engine import fq_quench, POLISH
    from pk.packing import read
    arm, n, s, sq, sig, seed = t
    sq = np.array(sq); rng = np.random.default_rng(seed)
    L = 1.05; c = s / 2
    sq = sq.copy(); sq[:, :2] = c * L + L * (sq[:, :2] - c)
    sq = sym_kick(sq, sig, rng) if arm == 'sym' else sq + np.c_[rng.normal(0, sig, (len(sq), 2)), rng.normal(0, 20 * sig, len(sq))]
    with tempfile.TemporaryDirectory() as tmp:
        a, b = os.path.join(tmp, 'a'), os.path.join(tmp, 'b')
        write_cfg(a, s * L, sq)
        args = [HARDEN, 'reharden', '--in', a, '--out', b, '--seed', str(seed), '--hold-move', '1', '--r0', '0.1']
        if arm == 'sym':
            args += ['--sym-md', '1']
        subprocess.run(args, capture_output=True, text=True, timeout=900)
        run = fq_quench(b, POLISH, tmp=tmp, timeout=900)
        p = read(run.out) if run.out else None
    if p is None:
        return dict(arm=arm, n=n, sig=sig, seed=seed, s=None)
    return dict(arm=arm, n=n, sig=sig, seed=seed, s=float(p.s), gap=float(p.s - s), pen=float(p.max_pen()))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--ns', required=True); ap.add_argument('--sig', default='0.05,0.1,0.2,0.4')
    ap.add_argument('--k', type=int, default=12); ap.add_argument('--out', required=True); ap.add_argument('--procs', type=int, default=15)
    a = ap.parse_args()
    from multiprocessing import Pool
    from pk.store import Store
    os.makedirs(a.out, exist_ok=True)
    sc = json.load(open('runs/symsub/symscan.json')); st = Store()
    T = []
    for n in map(int, a.ns.split(',')):
        p = st.get(sc[str(n)]['pid'])
        for sig in map(float, a.sig.split(',')):
            for i in range(a.k):
                for arm in ('free', 'sym'):
                    T.append((arm, n, float(p.s), p.sq.tolist(), sig, 100 * n + i))
    with Pool(a.procs) as pool:
        R = pool.map(job, T)
    json.dump(R, open(os.path.join(a.out, 'runs.json'), 'w'))
    g = collections.defaultdict(list)
    for r in R:
        if r['s'] is not None and r['pen'] <= 1e-9:
            g[(r['n'], r['sig'], r['arm'])].append(r['gap'])
    print('n sig: free return / median gap | sym return / median gap')
    for n in map(int, a.ns.split(',')):
        for sig in map(float, a.sig.split(',')):
            f, s_ = g[(n, sig, 'free')], g[(n, sig, 'sym')]
            fr = sum(abs(x) < 1e-7 for x in f); sr = sum(abs(x) < 1e-7 for x in s_)
            print(n, sig, f'{fr}/{len(f)} {np.median(f):+.3f} | {sr}/{len(s_)} {np.median(s_):+.3f}')
    for arm in ('free', 'sym'):
        v = [x for k, vv in g.items() if k[2] == arm for x in vv]
        print('TOTAL', arm, sum(abs(x) < 1e-7 for x in v), '/', len(v))
