#!/usr/bin/env python3
"""Do role counts separate the n = 110 sub-11 minima into distance clusters? (10-11; Evan: region-count holdouts.)

Sample K minima of a refset (representatives), junction-aware role counts (explore.roles_j: L, B, axis, J ~45 deg, other),
all pairwise matching distances (pk.crossing.match_distance, ~0.2 s each) on --procs.  Report: within-group vs
between-group distance distributions (group = full count tuple, and = axis count only), nearest-neighbour group agreement,
and the group table.  Output json: pairs + roles, reusable for distance-rule calibration.

  count_sep.py --ref census110_1011 --k 150 --procs 8 --out runs/count_sep/k150.json
"""
import argparse, collections, itertools, json, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from jobpool import run_jobs


def job(arg):
    i, j, a, b = arg
    from pk.store import Store
    from pk.crossing import match_distance
    st = Store()
    return i, j, match_distance(st.get(a), st.get(b))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', default='census110_1011'); ap.add_argument('--k', type=int, default=150)
    ap.add_argument('--procs', type=int, default=8); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default=f'{HERE}/runs/count_sep/k150.json')
    a = ap.parse_args()
    import scratch; scratch.scope('count_sep')
    from pk import refset, manifest
    from pk.store import Store
    import explore
    st = Store(); R = refset.load(a.ref)
    M = [m for m in R['minima'] if m['rep'] and m['cls'] == 'certified']
    M = random.Random(a.seed).sample(M, min(a.k, len(M)))
    roles = [explore.roles_j([tuple(map(float, r)) for r in st.get(m['rep']).sq]) for m in M]
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    man = manifest.start(os.path.dirname(a.out), 'count_sep', seed=a.seed, args=vars(a), procs=a.procs)
    jobs = [(i, j, M[i]['rep'], M[j]['rep']) for i, j in itertools.combinations(range(len(M)), 2)]
    res = [r for r in run_jobs(job, jobs, procs=a.procs, timeout=60) if r]
    man.finish()
    json.dump(dict(ref=a.ref, ref_sha256=R['sha256'], keys=[m['key'] for m in M], S=[m['S'][:20] for m in M], roles=roles,
                   pairs=res), open(a.out, 'w'))
    report(a.out)


def report(path):
    import numpy as np
    d = json.load(open(path))
    K = len(d['keys']); D = np.full((K, K), np.nan)
    for i, j, x in d['pairs']: D[i, j] = D[j, i] = x
    ro = [tuple(r) for r in d['roles']]; S = np.array([float(s) for s in d['S']])
    print(f'{K} minima, {len(d["pairs"])} pairs; groups (L, B, axis, J, other): '
          f'{collections.Counter(ro).most_common(12)}')
    for name, g in (('full tuple', ro), ('axis count', [r[2] for r in ro]), ('J count', [r[3] for r in ro])):
        w = [D[i, j] for i, j in itertools.combinations(range(K), 2) if g[i] == g[j]]
        b = [D[i, j] for i, j in itertools.combinations(range(K), 2) if g[i] != g[j]]
        nn = [g[int(np.nanargmin(np.where(np.arange(K) == i, np.nan, D[i])))] == g[i] for i in range(K)]
        base = sum(v * (v - 1) for v in collections.Counter(g).values()) / (K * (K - 1))
        q = lambda v: ' '.join(f'{x:.3f}' for x in np.percentile(v, [10, 50, 90])) if v else '-'
        print(f'  {name:10}: within p10/50/90 {q(w)} ({len(w)}) | between {q(b)} ({len(b)}) | '
              f'nearest neighbour same group {np.mean(nn):.0%} (chance {base:.0%})')
    rec = int(np.argmin(S))
    print(f'  distance to the best sampled ({S[rec]:.10f}): p10/50/90 {np.percentile(np.delete(D[rec], rec), [10, 50, 90]).round(3)}')


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--report':
        report(sys.argv[2])
    else:
        main()
