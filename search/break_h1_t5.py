#!/usr/bin/env python3
"""break_h1_t5: the T = 5, n = 20 top-of-hole chain-free cell (k = 15, 16) -- continuation of
the structured winner of `break_h1.py tophole` in eps, the k = 15 -> 16 extension, and the
CW / H / H1 / H2 ladder + weight-level depth at each eps.

Imports only; nothing in search/ is modified."""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import break_h1 as B                                             # noqa: E402
import s6skel                                                    # noqa: E402

T, n = 5.0, 20


def extend_to_16(rec):
    """the k = 15 winner has one level cell of the 4 x 4 grid empty and a TILTED square sitting
    where it would be (exactly the T = 4 `k = 8 -> k = 9` move: the extra cut-set square is a far
    square at tilt exactly eps).  -> (theta, pattern, perm) for k = 16, or None."""
    z = np.array(rec['z'])
    pat = [tuple(c) for c in rec['pattern']]
    k = len(pat)
    missing = [(a, b) for a in range(4) for b in range(4) if (a, b) not in set(pat)]
    if len(missing) != 1:
        return None
    X, Y = z[1::3], z[2::3]
    # the free square nearest the centroid of the missing level cell's neighbours
    xs = {}
    ys = {}
    for q, (a, b) in enumerate(pat):
        xs.setdefault(a, []).append(X[q])
        ys.setdefault(b, []).append(Y[q])
    ma, mb = missing[0]
    tx = np.mean(xs[ma]) if ma in xs else 2.5
    ty = np.mean(ys[mb]) if mb in ys else 2.5
    free = list(range(k, n))
    j = min(free, key=lambda q: (X[q] - tx) ** 2 + (Y[q] - ty) ** 2)
    perm = list(range(k)) + [j] + [q for q in free if q != j]
    return perm, [list(c) for c in pat] + [list(missing[0])]


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else 'runs/break_h1_t5_cont.jsonl'
    epss = [0.125, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
    rs = [json.loads(line) for line in open(src)]
    best = {}
    for r in rs:
        if not r.get('cfree'):
            continue
        key = (r['k'], r['eps'])
        if key not in best or r['delta'] > best[key]['delta']:
            best[key] = r
    base = min((r for r in best.values() if r['k'] == 15),
               key=lambda r: r['delta'] / math.radians(r['eps']) ** 2)
    print(f"# seed: {base['lab']}  k={base['k']} eps={base['eps']}g "
          f"delta*={base['delta']:+.6e}", flush=True)
    ext = extend_to_16(base)
    rng = np.random.default_rng(11)
    out = open('runs/break_h1_t5_cells.jsonl', 'w')
    for (k, pat, perm, tag) in ([(15, [tuple(c) for c in base['pattern']],
                                 list(range(n)), 'k15')]
                                + ([(16, [tuple(c) for c in ext[1]], ext[0], 'k16')]
                                   if ext else [])):
        z0 = np.array(base['z'])
        X0, Y0 = z0[1::3][perm], z0[2::3][perm]
        TH0 = z0[3::3][perm]
        print(f'\n## {tag}  (cut set of {k}, pattern {pat})', flush=True)
        prev = None
        for epsdeg in sorted(epss, reverse=True):
            eps = math.radians(epsdeg)
            th = [(0.0 if abs(s6skel.norm_tilt(t)) < 1e-9 else (eps if t > 0 else -eps))
                  for t in TH0]
            seeds = [(X0.copy(), Y0.copy())]
            if prev is not None:
                seeds.append((prev[0].copy(), prev[1].copy()))
            for _ in range(40):
                bx, by = (prev if prev is not None else (X0, Y0))
                seeds.append((bx + rng.normal(0, 0.03, n), by + rng.normal(0, 0.03, n)))
            v, z, hit = B.delta_star_at(th, T, rng, seeds=seeds, limit=0, jitters=0, extra=0,
                                        pattern=pat, k=k)
            if z is None:
                print(f'  eps={epsdeg:g}: infeasible', flush=True)
                continue
            cf = B.chainfree_ok(z, n, T, v, range(k))
            prev = (z[1::3].copy(), z[2::3].copy())
            o = B.bounds(z, n, T, v, tol=1e-10, maxchains=int(os.environ.get('MAXCH', 120)),
                         eps_rad=eps, kmax=7, do_cyc=False, nch2=1, cap2=1500)
            print(f"  eps={epsdeg:6g}  delta*={v:+.8e}  /e^2 {-v/eps**2:.5f} "
                  f"/e^3 {-v/eps**3:.4f} /e^4 {-v/eps**4:.3f}  cfree={cf} hit={hit}", flush=True)
            print(f"      CW {o['CW']:+.4e}  H {o['H']:+.4e}  HP {o['HP']:+.4e}  "
                  f"H1 {o['H1']:+.4e}  H2 {o['H2']:+.4e}  F {o['F']:+.6e}  "
                  f"cdepth={o['cdepth']} wdepth={o['wdepth']} mu={o['mu']} "
                  f"nchains={o['nchains']} ({o['secs']:.0f}s)", flush=True)
            print(f"      wlevels {[(l, c, float('%.3e' % w)) for l, c, w in o['wlevels']]}",
                  flush=True)
            rec = dict(lab=f'T5 tophole {tag}', k=k, eps=epsdeg, delta=float(v),
                       cfree=bool(cf), hit=int(hit), pattern=[list(c) for c in pat],
                       z=[float(q) for q in z])
            for key in ('CW', 'H', 'HP', 'H1', 'H2', 'F', 'cdepth', 'wdepth', 'wlevels',
                        'mu', 'E', 'V', 'cls', 'nchains', 'supp'):
                rec[key] = o[key]
            out.write(json.dumps(rec) + '\n'); out.flush()
    out.close()


if __name__ == '__main__':
    main()
