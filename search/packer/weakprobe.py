#!/usr/bin/env python3
"""Is a weak (degenerate) exact minimum a local minimum?  Radius probe (10-15).

m(t) = S* - min { S : packing, |v - v*|_inf <= t }   (v = all x, y, theta; S free), estimated by IPOPT (nlpdesc's
fixed-line NLP, lines picked at each start) from K random starts v* + t xi, xi uniform in the box; the result is repaired
to an exactly-checked f64 packing (scaling about the origin) before S is read.  A local minimum gives m(t) <= noise
(~1e-13) for all small t; a non-minimum gives m(t) ~ c t^p (p = 3, 4 along a degenerate zero mode).  Fixed lines: a
corner-corner pair only gets the branch its start point prefers, so random starts are also the branch sampler.

  ./ipenv.sh python3 weakprobe.py EXACT.txt OUT.json [--ts 1e-2,3e-3,1e-3,3e-4] [--K 12] [--seed 0]
"""
import argparse, json, math, sys, time
import numpy as np
import cyipopt
import nlpdesc as N


def load_exact(path):
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n = int(L[0][0]); S = L[0][1]
    sq = np.array([(float(a), float(b), math.radians(float(c))) for a, b, c, *_ in L[1:n + 1]])
    return S, sq


def box_solve(center_s, center_sq, start_s, start_sq, t, tol=1e-13):
    n = len(center_sq)
    cons, near = N.build(start_s, start_sq, 0.0)
    P = N.Prob(n, cons)
    vc = np.concatenate([center_sq.ravel(), [center_s]])
    v0 = np.concatenate([start_sq.ravel(), [start_s]])
    lb, ub = vc - t, vc + t
    lb[-1], ub[-1] = center_s - 1.0, center_s + 1.0
    nlp = cyipopt.Problem(n=P.nv, m=len(cons), problem_obj=P, lb=lb, ub=ub, cl=np.zeros(len(cons)), cu=np.full(len(cons), 1e20))
    for k, val in (('tol', tol), ('constr_viol_tol', 1e-14), ('acceptable_tol', 1e-11), ('max_iter', 500),
                   ('print_level', 0), ('mu_strategy', 'adaptive'), ('sb', 'yes'), ('nlp_scaling_method', 'none'),
                   ('honor_original_bounds', 'yes')):
        nlp.add_option(k, val)
    v, info = nlp.solve(v0)
    s1, sq1 = v[-1], v[:-1].reshape(n, 3)
    s2, sq2, p = N.repair(s1, sq1, N.near_pairs(sq1))
    return float(s2), sq2, info['status']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('exact'); ap.add_argument('out')
    ap.add_argument('--ts', default='1e-2,3e-3,1e-3,3e-4'); ap.add_argument('--K', type=int, default=12)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--cont', type=int, default=10)
    a = ap.parse_args()
    Sstr, sq = load_exact(a.exact); S = float(Sstr)
    rng = np.random.default_rng(a.seed)
    res = dict(S=Sstr, runs=[])
    best = (0.0, None, None)
    for t in [float(x) for x in a.ts.split(',')]:
        ms = []
        for k in range(a.K):
            t0 = time.time()
            st_sq = sq + t * rng.uniform(-1, 1, sq.shape)
            s2, sq2, status = box_solve(S, sq, S + t, st_sq, t)
            m = S - s2
            ms.append(m)
            if m > best[0]:
                best = (m, s2, sq2)
            res['runs'].append(dict(t=t, k=k, m=m, status=int(status), secs=round(time.time() - t0, 1)))
        print(f't={t:.0e}: m max {max(ms):.3e}  median {np.median(ms):.3e}  (#m>1e-12: {sum(x > 1e-12 for x in ms)}/{a.K})',
              flush=True)
    res['m_by_t'] = {str(t): max(r['m'] for r in res['runs'] if r['t'] == t) for t in sorted({r['t'] for r in res['runs']})}
    # continuation: shrink the deepest descent found toward v* (start v* + (t / t0)(v_best - v*)), t0 / 3, t0 / 9, ...
    res['cont'] = []
    if best[2] is not None and best[0] > 1e-9:
        vb = np.concatenate([best[2].ravel(), [best[1]]]); vs = np.concatenate([sq.ravel(), [S]])
        t0 = float(np.abs(vb[:-1] - vs[:-1]).max())
        prev = vb
        for j in range(1, a.cont + 1):
            t = t0 / 3 ** j
            w = vs + (t / t0) * (vb - vs) if j == 1 else vs + (prev - vs) / 3
            s2, sq2, status = box_solve(S, sq, float(w[-1]), w[:-1].reshape(-1, 3), t)
            m = S - s2
            res['cont'].append(dict(t=t, m=m, status=int(status)))
            print(f'  cont t={t:.2e}: m {m:.3e}', flush=True)
            if m <= 1e-10:
                break
            prev = np.concatenate([sq2.ravel(), [s2]])
        res['t0'] = t0
    json.dump(res, open(a.out, 'w'))
    if best[2] is not None:
        N.save(a.out.replace('.json', '.best.txt'), best[1], best[2])


if __name__ == '__main__':
    main()
