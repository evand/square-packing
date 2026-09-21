#!/usr/bin/env python3
"""bandcut_khunt.py -- stage 2 of `tasks/bandcut-k`: hunt for a CHAIN-FREE configuration of
margin 0 (2026-09-20).

The stage-1 scan (`search/bandcut_k.py scan`) searches each (k, eps, pattern) cell from its own
starts, so the chain-free runs never see the zero-margin configurations that the unrestricted runs
find.  This stage closes that loop, which is the sharpest test available for a counterexample:

  * take every configuration stage 1 stored (in particular the delta* = 0 ones, which DO contain a
    wall-to-wall chain of T among their near-axis squares), relabel it so that the k least-tilted
    squares are the near-axis ones, clamp it into the band, and re-solve it under EVERY chain-free
    level pattern.  If a zero-margin configuration can be deformed into a chain-free one, this is
    what finds it;
  * plus synthetic Cleemann-type starts: a bent band of squares at 2 arctan(1/T) (the T = 4 window
    angle, `search/ARCH_TLEDGER.md` sec 5.5) cutting the container, with axis-parallel blocks on
    either side -- the structure that at T = 17 has T <= k <= (T-1)^2 and no chain.

    python3 search/bandcut_khunt.py --n 12 --T 4 --eps 1,5,10 --seeds runs/bandcut_k_T4.jsonl \\
        --out runs/bandcut_k_T4_hunt.jsonl --nproc 9
"""
import argparse
import json
import math
import os
import sys

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bandcut_k as B                                            # noqa: E402
import s6skel                                                    # noqa: E402


def relabel(z, n, k):
    """reorder the squares so that the k least tilted come first (they become the near-axis set)."""
    order = sorted(range(n), key=lambda q: abs(s6skel.norm_tilt(z[3 + 3 * q])))
    w = z.copy()
    for new, old in enumerate(order):
        w[1 + 3 * new:4 + 3 * new] = z[1 + 3 * old:4 + 3 * old]
    return w


def band_start(rng, n, T, k, eps, tiltdeg):
    """a Cleemann-type start: a diagonal band of n-k squares at +-tiltdeg cutting the container,
    the k near-axis squares on tiling cells on either side of it."""
    Ti = int(round(T))
    t = math.radians(tiltdeg) * (1.0 if rng.random() < 0.5 else -1.0)
    z = np.zeros(1 + 3 * n)
    c = rng.uniform(-0.6, 0.6)
    # band: centres along the anti-diagonal x + y = T + c, spaced by ~1 along it
    nb = n - k
    d = np.array([math.cos(t + math.pi * 0.75), math.sin(t + math.pi * 0.75)])
    mid = np.array([T / 2 + c / 2, T / 2 + c / 2])
    for q in range(nb):
        p = mid + (q - (nb - 1) / 2.0) * 1.05 * d
        z[1 + 3 * (k + q)] = p[0] + rng.normal(0, 0.05)
        z[2 + 3 * (k + q)] = p[1] + rng.normal(0, 0.05)
        z[3 + 3 * (k + q)] = t
    cells = [(a + 0.5, b + 0.5) for a in range(Ti) for b in range(Ti)]
    far = sorted(cells, key=lambda cc: -abs(cc[0] + cc[1] - (T + c)))
    for q in range(k):
        cc = far[q % len(far)]
        z[1 + 3 * q] = cc[0] + rng.normal(0, 0.05)
        z[2 + 3 * q] = cc[1] + rng.normal(0, 0.05)
        z[3 + 3 * q] = rng.uniform(-1, 1) * eps * 0.2
    return z


_G = {}


def _init(a, seeds):
    _G['a'], _G['seeds'] = a, seeds


def _job(arg):
    (k, epsdeg, pidx, pattern) = arg
    a = _G['a']
    n, T = a.n, a.T
    eps = math.radians(epsdeg)
    rng = np.random.default_rng(hash((k, pidx, int(epsdeg * 10))) % (2 ** 31))
    tw = 2 * math.degrees(math.atan(1.0 / T))            # the straight-chain window angle

    def run(z0):
        signs = B.start_signs(z0, n, k)
        pr = B.BandProblem(n, T, k, eps, signs, pattern=pattern, eta=a.eta)
        best = (-1e18, None)
        z, v = B.lp_polish(z0, n, T, k, pattern, a.eta)
        if v > best[0]:
            best = (v, z.copy())
        for _cyc in range(2):
            v, z = B.solve_band(pr, z, rounds=8, maxiter=200)
            if z is None:
                break
            if v > best[0]:
                best = (v, z.copy())
            z, v = B.lp_polish(z, n, T, k, pattern, a.eta)
            if v > best[0]:
                best = (v, z.copy())
        v, z = best
        if z is None:
            return None
        z[0] = v
        ok, info = B.verify(z, n, T, k, eps, v, tol_edge=a.tol_edge, need_chainfree=True)
        return (v, z, info) if ok else None

    starts = [B.make_start(rng, n, T, k, eps, pattern,
                           ['spread', 'compact', 'jitgrid', 'tilingall', 'random'][s % 5])
              for s in range(a.nstart)]
    for z in _G['seeds']:
        w = B.clamp_band(relabel(np.array(z), n, k), n, k, eps)
        starts.append(w)
        for (rho, al) in ((0.03, 0.3 * eps), (0.12, 0.7 * eps)):
            starts.append(B.nudge(w, rng, n, T, k, eps, rho, al))
    for _ in range(a.bands):
        for td in (tw, 45.0, max(tw, math.degrees(eps))):
            if math.radians(td) >= eps - 1e-12:
                starts.append(band_start(rng, n, T, k, eps, td))
    keep, vals = [], []
    for z0 in starts:
        r = run(z0)
        if r is None:
            continue
        vals.append(r[0])
        keep.append(r)
    keep.sort(key=lambda r: -r[0])
    keep = keep[:4]
    for _rd in range(a.polish):
        newk = list(keep)
        for (v, z, _i) in keep:
            for (rho, al) in ((0.02, 0.2 * eps), (0.08, 0.5 * eps)):
                for _ in range(3):
                    r = run(B.nudge(z, rng, n, T, k, eps, rho, al))
                    if r is not None:
                        vals.append(r[0])
                        newk.append(r)
        newk.sort(key=lambda r: -r[0])
        keep = newk[:4]
    pat = [list(c) for c in pattern]
    if not keep:
        return dict(k=k, eps=epsdeg, pidx=pidx, value=None, nfeas=0, nstart=len(starts), nhit=0,
                    pattern=pat, z=None, info=None)
    v, z, info = keep[0]
    return dict(k=k, eps=epsdeg, pattern=pat, pidx=pidx, value=float(v), nfeas=len(vals),
                nstart=len(vals), nhit=sum(1 for u in vals if u >= v - 1e-9),
                z=[float(q) for q in z],
                info={kk: (vv if not isinstance(vv, list) else list(vv))
                      for kk, vv in info.items()})


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--T', type=float, required=True)
    ap.add_argument('--eps', default='1,5,10')
    ap.add_argument('--ks', default=None)
    ap.add_argument('--seeds', nargs='*', default=[])
    ap.add_argument('--maxseeds', type=int, default=60)
    ap.add_argument('--bands', type=int, default=12)
    ap.add_argument('--nstart', type=int, default=40)
    ap.add_argument('--polish', type=int, default=2)
    ap.add_argument('--eta', type=float, default=1e-6)
    ap.add_argument('--tol-edge', dest='tol_edge', type=float, default=1e-7)
    ap.add_argument('--nproc', type=int, default=9)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    seeds = []
    for f in a.seeds:
        for line in open(f):
            r = json.loads(line)
            if r.get('z'):
                seeds.append((r['value'], r['z']))
    seeds.sort(key=lambda s: -s[0])
    seeds = [z for (_v, z) in seeds[:a.maxseeds]]
    print(f"# {len(seeds)} seed configurations", flush=True)

    ks = ([int(s) for s in a.ks.split(',')] if a.ks
          else list(range(1, (int(round(a.T)) - 1) ** 2 + 1)))
    jobs = [(k, float(e), pidx, pat) for e in a.eps.split(',') for k in ks
            for pidx, pat in enumerate(B.level_patterns(a.T, k))]
    print(f"# {len(jobs)} jobs", flush=True)
    with mp.Pool(a.nproc, initializer=_init, initargs=(a, seeds)) as pool, open(a.out, 'w') as f:
        for i, r in enumerate(pool.imap_unordered(_job, jobs)):
            f.write(json.dumps(r) + '\n')
            f.flush()
            print(f"  [{i + 1}/{len(jobs)}] k={r['k']:<2d} eps={r['eps']:<4g} pat{r['pidx']:<3d} "
                  f"value={r['value'] if r['value'] is None else round(r['value'], 9)} "
                  f"hit {r['nhit']}/{r['nstart']}", flush=True)


if __name__ == '__main__':
    main()
