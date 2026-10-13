#!/usr/bin/env python3
"""Is depth regional? (10-13, Evan: "do regions have good draws, or are draws random and we just need diversity?")

Population: the certified exact minima at n (store `minimum` table), each represented by its nearest store packing.
Distance: pk.fingerprint (smoothed centre images, min over 8 symmetries; calibrated against matching distance).

Part 1, the landscape (no search data):
  variogram      gamma(d) = mean (s_i - s_j)^2 / 2 over pairs in distance bins, vs a side-shuffled null: rising gamma
                 = nearby minima have similar depth (regions); flat = depth is spatially random.
  kNN            Spearman(side, mean side of the k nearest minima), vs shuffled.
  concentration  top-K deepest: mean pairwise distance and the share whose nearest other minimum is also top-K, vs
                 random K-subsets.
Part 2, the chain (explorer proposals with parent links: census v2-v4 archives + proposals.jsonl):
  heritability   child side vs parent side over all proposals with a sub-k child; P(child in the deep tail | parent
                 depth decile); proposal yields per decile = the weights a parent-selection rule should use.

  python3 landscape_depth.py [--n 110] [--runs 'runs/census_v*/*_r*'] [--out runs/landscape_depth.json]
"""
from __future__ import annotations
import argparse, glob, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def population(n):
    from pk.store import Store
    st = Store()
    M = st.db.execute("SELECT CAST(S AS REAL) AS S FROM minimum WHERE n=? AND cls='certified' ORDER BY S", (n,)).fetchall()
    out = []
    for m in M:
        r = st.db.execute("SELECT id FROM packing WHERE n=? AND status != 'infeasible' AND s BETWEEN ? AND ? "
                          "ORDER BY abs(s - ?) LIMIT 1", (n, m['S'] - 1e-8, m['S'] + 1e-8, m['S'])).fetchone()
        if r:
            p = st.get(r['id'])
            out.append((float(m['S']), p))
    return out


def distances(pop, n):
    from pk import fingerprint as F
    ims = np.stack([F.image(float(p.s), p.sq) for _, p in pop])
    N = len(ims)
    D = np.zeros((N, N))
    for i in range(N):
        D[i] = F.nearest(ims[i], ims, n)
    return np.minimum(D, D.T)


def variogram(s, D, bins, rng, reps=200):
    iu = np.triu_indices(len(s), 1)
    d, ds2 = D[iu], (s[:, None] - s[None])[iu] ** 2 / 2
    q = np.quantile(d, np.linspace(0, 1, bins + 1))
    idx = np.clip(np.searchsorted(q, d, side='right') - 1, 0, bins - 1)
    g = np.array([ds2[idx == b].mean() for b in range(bins)])
    null = []
    for _ in range(reps):
        t = rng.permutation(s)
        ds2n = (t[:, None] - t[None])[iu] ** 2 / 2
        null.append([ds2n[idx == b].mean() for b in range(bins)])
    null = np.array(null)
    return q, g, null.mean(0), np.quantile(null, 0.025, 0), np.quantile(null, 0.975, 0)


def spearman(a, b):
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def knn(s, D, k, rng, reps=200):
    Dn = D + np.diag(np.full(len(s), np.inf))
    nb = np.argsort(Dn, 1)[:, :k]
    obs = spearman(s, s[nb].mean(1))
    null = [spearman(t, t[nb].mean(1)) for t in (rng.permutation(s) for _ in range(reps))]
    return obs, float(np.mean(null)), float(np.quantile(null, 0.975))


def concentration(s, D, K, rng, reps=500):
    order = np.argsort(s)
    top = order[:K]
    def stats(I):
        sub = D[np.ix_(I, I)]
        mpd = sub[np.triu_indices(K, 1)].mean()
        Dn = D[I] + 0.0
        Dn[np.arange(K), I] = np.inf
        nn = np.argmin(Dn, 1)                               # nearest other minimum (any)
        return mpd, float(np.isin(nn, I).mean())
    o = stats(top)
    null = np.array([stats(rng.choice(len(s), K, replace=False)) for _ in range(reps)])
    return dict(K=K, mean_pair_dist=o[0], null_mpd=float(null[:, 0].mean()), p_mpd=float((null[:, 0] <= o[0]).mean()),
                nn_in_top=o[1], null_nn=float(null[:, 1].mean()), p_nn=float((null[:, 1] >= o[1]).mean()))


def chain(run_glob, k, s_tail):
    """(parent side, child side or None (grid / discard)) for every logged proposal."""
    pairs = []
    for d in sorted(glob.glob(run_glob)):
        ap, pp = os.path.join(d, 'archive.jsonl'), os.path.join(d, 'proposals.jsonl')
        if not (os.path.exists(ap) and os.path.exists(pp)):
            continue
        A = {}
        for l in open(ap):
            e = json.loads(l); A[e['i']] = e['s']
        for l in open(pp):
            p = json.loads(l)
            par = p.get('par')
            if par is None or par not in A:
                continue
            s = p.get('s')
            ok = s is not None and s < k - 1e-7 and p.get('st') not in ('screen-grid', 'discard', 'screen-discard')
            pairs.append((A[par], s if ok else None, p.get('cpu') or p.get('sec') or 0.0))
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=110); ap.add_argument('--runs', default='runs/census_v*/*_r*')
    ap.add_argument('--out', default='runs/landscape_depth.json'); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    pop = population(a.n)
    s = np.array([x[0] for x in pop])
    print(f'{len(s)} certified minima, sides {s.min():.10f} .. {s.max():.10f}', flush=True)
    D = distances(pop, a.n)
    out = dict(n=a.n, N=len(s))
    q, g, gn, lo, hi = variogram(s, D, 10, rng)
    print('variogram (distance-decile upper edge: gamma vs shuffled 95% band), gamma in 1e-8 units:')
    for b in range(10):
        flag = ' BELOW null' if g[b] < lo[b] else (' above' if g[b] > hi[b] else '')
        print(f'  d <= {q[b + 1]:.3f}: {g[b] * 1e8:8.3f}  null {gn[b] * 1e8:8.3f} [{lo[b] * 1e8:.3f}, {hi[b] * 1e8:.3f}]{flag}')
    out['variogram'] = dict(edges=q.tolist(), gamma=g.tolist(), null=gn.tolist(), lo=lo.tolist(), hi=hi.tolist())
    out['knn'] = {}
    for kk in (1, 5, 20):
        o, m, h = knn(s, D, kk, rng)
        out['knn'][kk] = (o, m, h)
        print(f'kNN k={kk:2d}: Spearman(side, neighbours\' mean side) = {o:.3f}  (shuffled mean {m:.3f}, 97.5 % {h:.3f})')
    out['concentration'] = []
    for K in (5, 10, 25, 50, 100):
        c = concentration(s, D, K, rng)
        out['concentration'].append(c)
        print(f"top {K:3d}: mean pair distance {c['mean_pair_dist']:.3f} vs random {c['null_mpd']:.3f} (p {c['p_mpd']:.3f});"
              f" nearest neighbour also top-{K}: {c['nn_in_top']:.2f} vs random {c['null_nn']:.2f} (p {c['p_nn']:.3f})")
    # Part 2
    P = chain(a.runs, k, None)
    ps = np.array([p[0] for p in P]); cs = np.array([np.nan if p[1] is None else p[1] for p in P])
    cpu = np.array([p[2] for p in P])
    print(f'\nchain: {len(P)} proposals with a known parent; {np.isfinite(cs).sum()} sub-{k} children')
    ok = np.isfinite(cs) & (ps < k)
    print(f'Spearman(parent side, child side) over sub-{k} children: {spearman(ps[ok], cs[ok]):.3f}')
    tails = {f'top{K}': np.sort(s)[K - 1] for K in (1, 10, 50)}
    dec = np.quantile(ps[ps < k], np.linspace(0, 1, 11))
    pd = np.clip(np.searchsorted(dec, ps, side='right') - 1, 0, 9)
    out['chain'] = dict(n=len(P), deciles=dec.tolist(), rows=[])
    print('parent depth decile (side range): proposals, sub-k rate, child in top-1 / top-10 / top-50 per 1000 proposals,'
          ' per CPU-h (top-50)')
    for b in range(10):
        m = (pd == b) & (ps < k)
        if not m.any():
            continue
        c = cs[m]
        row = dict(decile=b, lo=float(dec[b]), hi=float(dec[b + 1]), props=int(m.sum()), subk=float(np.isfinite(c).mean()),
                   cpu_h=float(cpu[m].sum() / 3600))
        for name, t in tails.items():
            row[name] = float(np.nansum(c <= t + 1e-9) / m.sum() * 1000)
        row['top50_per_cpuh'] = float(np.nansum(c <= tails['top50'] + 1e-9) / max(row['cpu_h'], 1e-9))
        out['chain']['rows'].append(row)
        print(f"  {b} ({dec[b]:.7f}-{dec[b + 1]:.7f}): {row['props']:6d}  {row['subk']:.2f}  "
              f"{row['top1']:6.2f} {row['top10']:6.2f} {row['top50']:6.2f}   {row['top50_per_cpuh']:7.2f}")
    json.dump(out, open(a.out, 'w'), indent=1)
    print(a.out)


if __name__ == '__main__':
    main()
