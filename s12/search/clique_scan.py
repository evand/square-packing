#!/usr/bin/env python3
"""Best segment-anchored anchor clique on a measure, with the anchor OFFSET CAP LIFTED.

`clique_family.anchor_local` caps the offset at

    clear = min(p_x, p_y, t-p_x, t-p_y);  base = max(1e-4, min(0.6, 1-clear if clear < 1 else 0.05))

so at an INTERIOR anchor point (`clear >= 1`) it only ever tries `eps <= 0.05`.  The cliques these
measures actually violate hardest are interior with `eps` of 0.2-0.45, so the cap is what hides
them.  This scans `eps` up to `--epsmax` at every anchor point.
"""
import sys, os, math, time, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clique_family as CF

path = sys.argv[1]
t = float(sys.argv[2]) if len(sys.argv) > 2 else 3.99
maxpts = int(sys.argv[3]) if len(sys.argv) > 3 else 200
epsmax = float(sys.argv[4]) if len(sys.argv) > 4 else 0.9
NDIR, NEPS, NRHO = 24, 12, 10

imgs, poses = CF.load_measure(path, t)
IM = np.array([(a, b, c, d) for (a, b, c, d) in imgs if d > 1e-12])
print(f"# {path}: {len(IM)} images, mass {IM[:,3].sum():.9f}, epsmax {epsmax}")
cand = CF.tight_points(imgs, t, 0.01, 0.95)
cand.sort(key=lambda z: -z[1])
step = max(1, len(cand) // maxpts)
cand = cand[::step][:maxpts]
dirs = [(math.cos(a), math.sin(a)) for a in np.linspace(0, 2 * math.pi, NDIR, endpoint=False)]
best = []
t0 = time.time()
for k, (p, cv) in enumerate(cand):
    thru = CF.contains_np(IM, [p])
    for d in dirs:
        dp = (-d[1], d[0])
        for ie in range(1, NEPS + 1):
            eps = epsmax * ie / NEPS
            mpt = (p[0] + eps * d[0], p[1] + eps * d[1])
            for ir in range(1, NRHO + 1):
                rho = 0.49 * ir / NRHO
                A = [(mpt[0] - rho * dp[0], mpt[1] - rho * dp[1]),
                     (mpt[0] + rho * dp[0], mpt[1] + rho * dp[1])]
                mm, n1, n2 = CF.kmass_np(IM, thru, A)
                if mm > 1.0 + 1e-9:
                    best.append(dict(mass=mm, p=p, A=A, eps=eps, rho=rho, cov=cv,
                                     n_point=n1, n_anchor=n2))
    if k % 25 == 0:
        print(f"  {k}/{len(cand)} best {max((z['mass'] for z in best), default=0):.6f} "
              f"({time.time()-t0:.0f}s)", flush=True)
best.sort(key=lambda z: -z['mass'])
print(f"{'K mass':>10} {'cov(p)':>9} {'walldist':>9} {'eps':>8} {'|A|':>8} {'#thru':>6} {'#anch':>6}  p")
seen = set()
shown = 0
for r in best:
    key = (round(r['p'][0], 5), round(r['p'][1], 5))
    if key in seen:
        continue
    seen.add(key)
    wd = min(r['p'][0], r['p'][1], t - r['p'][0], t - r['p'][1])
    print(f"{r['mass']:10.6f} {r['cov']:9.6f} {wd:9.4f} {r['eps']:8.5f} {2*r['rho']:8.5f} "
          f"{r['n_point']:6d} {r['n_anchor']:6d}  ({r['p'][0]:8.4f},{r['p'][1]:8.4f})")
    shown += 1
    if shown >= 15:
        break
if best:
    r = best[0]
    wd = min(r['p'][0], r['p'][1], t - r['p'][0], t - r['p'][1])
    print(f"\n# BEST SEGMENT anchor-clique mass: {r['mass']:.9f} (violation {r['mass']-1:+.6f}) "
          f"at wall distance {wd:.4f}, eps {r['eps']:.4f}, |A| {2*r['rho']:.4f}")
    nint = sum(1 for b in best
               if min(b['p'][0], b['p'][1], t - b['p'][0], t - b['p'][1]) >= 1.0)
    ncap = sum(1 for b in best if b['eps'] > 0.05
               and min(b['p'][0], b['p'][1], t - b['p'][0], t - b['p'][1]) >= 1.0)
    print(f"# {len(best)} violated cliques; {nint} at wall distance >= 1; "
          f"{ncap} of those with eps > 0.05 (invisible to clique_family.anchor_local)")
    json.dump([dict(mass=b['mass'], p=list(b['p']), A=[list(v) for v in b['A']],
                    eps=b['eps'], rho=b['rho']) for b in best[:100]],
              open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'scan_') +
                   os.path.basename(path).replace('.txt', '') + '.json', 'w'), indent=1)
