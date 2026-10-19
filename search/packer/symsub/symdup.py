#!/usr/bin/env python3
"""Blend vs duplicate-half symmetrisation (10-18, Evan: "blend them, or duplicate one half and see if that relaxes?").
Per TT-S1 start: blend (TT-S1's average), dup-lo, dup-hi; each relaxed once by symmetric soft-to-rigid (r0 0.05 / 0.1
/ 0.2, 2 seeds) + free polish.  Report the released side vs the start and today's record.
  symdup.py --out runs/symsub/dup1"""
import sys, os, json, argparse, tempfile, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import symtt


def job(t):
    n, rec, pid, s, sq, how, el, r0, seed = t
    sq = np.array(sq); rng = np.random.default_rng(seed)
    x = symtt.to_md(sq, s, el)
    if how == 'blend':
        y, ok = symtt.symmetrize(x, s)
    else:
        y = symtt.duplicate_half(x, s, n, how[4:], rng)
    with tempfile.TemporaryDirectory() as tmp:
        r = symtt.relax_once('sym', s, y, seed, tmp, r0)
    if r is None:
        return dict(n=n, pid=pid, how=how, r0=r0, seed=seed, s=None)
    hs, ps, psq = r
    return dict(n=n, pid=pid, how=how, el=el, r0=r0, seed=seed, rec=rec, s_start=s, hs=hs,
                s=ps if psq is not None else None, sq=None if psq is None else np.asarray(psq).tolist())


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--procs', type=int, default=15)
    a = ap.parse_args()
    from multiprocessing import Pool
    os.makedirs(a.out, exist_ok=True)
    T = []
    for n, rec, pid, s, sq in symtt.starts():
        el, _ = symtt.best_diag_mirror(sq, s)
        for how in ('blend', 'dup-lo', 'dup-hi'):
            for r0 in (0.05, 0.1, 0.2):
                for seed in (1, 2):
                    T.append((n, rec, pid, s, sq.tolist(), how, el, r0, seed))
    with Pool(a.procs) as pool:
        R = pool.map(job, T)
    json.dump(R, open(os.path.join(a.out, 'runs.json'), 'w'))
    g = collections.defaultdict(list)
    for r in R:
        if r['s'] is not None:
            g[(r['n'], r['pid'], r['how'])].append(r['s'] - r['rec'])
    print('n pid start-rec | blend best | dup-lo best | dup-hi best   (released side - record)')
    for (n, pid) in sorted({(k[0], k[1]) for k in g}):
        st = next(r for r in R if r['n'] == n and r['pid'] == pid)
        cells = [min(g[(n, pid, h)], default=float('nan')) for h in ('blend', 'dup-lo', 'dup-hi')]
        print(n, pid, f"{st['s_start'] - st['rec']:+.4f} |", ' | '.join(f'{c:+.4f}' for c in cells),
              '<<< RECORD' if min(cells) < 1e-7 else '')
