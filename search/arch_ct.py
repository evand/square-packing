#!/usr/bin/env python3
"""arch_ct.py -- the second-order constant c_T of the uniform tilt at the T x T tiling (2026-09-20).

`notes/proof-architecture.md` §3.2 records  delta*(t,...,t) = -c_T t^2 + O(t^3)  at
n = T^2 - T squares in [0,T]^2, with c_2 = 1, c_3 = 3/4, c_4 = 1/3 exact on their leaves and
c_5 <= 0.3745 measured (non-monotone => suspected under-optimisation).  This script attacks the
T = 5, 6 cases with *hole-set* starts instead of `s6local.py`'s random n-subsets:

    n = T^2 - T squares on the T x T tiling  <=>  T of the T^2 cells are EMPTY.

The hole set is what the configuration is; the brief's guess is that the optimum has one hole per
row and one per column (a permutation matrix), of which there are only T! -- 120 at T = 5, 720 at
T = 6 -- so the search space collapses from C(25,20) = 53,130 / C(36,30) = 1,947,792 to something
exhaustive, and every hole set gets many more jitters than a blanket multistart can afford.

Everything printed is a LOWER bound on delta* (a feasible point of the disjunctive LP of
`search/S6_SKELETON.md` §3.1, evaluated geometrically by `s6skel.value`), hence an UPPER bound on
c_T.  `s6exact.py` then turns the winning leaf into a closed form.

    python3 search/arch_ct.py holes --T 4 --family all  --tdeg 0.8                 # 1,820 hole sets
    python3 search/arch_ct.py holes --T 5 --family perm --tdeg 0.4,0.8,1.6 -j 24
    python3 search/arch_ct.py holes --T 5 --family all  --tdeg 0.8 -j 2
    python3 search/arch_ct.py window                                               # §3.2 straight-chain window
"""
import argparse
import itertools
import math
import os
import sys
import time

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np                                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6local                                                        # noqa: E402


# --------------------------------------------------------------------------- hole sets
def dihedral(holes, T):
    """the 8 images of a hole set under the symmetry group of the T x T grid."""
    out = []
    for tr in (False, True):
        hs = [(c, r) for (r, c) in holes] if tr else list(holes)
        for fx in (False, True):
            for fy in (False, True):
                out.append(frozenset(((T - 1 - r) if fy else r, (T - 1 - c) if fx else c)
                                     for (r, c) in hs))
    return out


def canon(holes, T):
    return min(tuple(sorted(h)) for h in dihedral(holes, T))


def perm_type(holes, T):
    """if the hole set is a permutation matrix, its cycle type (a partition of T); else None."""
    if len({r for (r, _c) in holes}) != T or len({c for (_r, c) in holes}) != T:
        return None
    sigma = {r: c for (r, c) in holes}
    seen, cyc = set(), []
    for r in range(T):
        if r in seen:
            continue
        ln, x = 0, r
        while x not in seen:
            seen.add(x); x = sigma[x]; ln += 1
        cyc.append(ln)
    return tuple(sorted(cyc, reverse=True))


def shape(holes, T):
    """a coarse invariant: sorted row counts and column counts of the holes."""
    rc = tuple(sorted(sum(1 for (r, _c) in holes if r == i) for i in range(T)))
    cc = tuple(sorted(sum(1 for (_r, c) in holes if c == j) for j in range(T)))
    return (rc, cc)


def hole_family(T, family, rng, limit):
    cells = [(r, c) for r in range(T) for c in range(T)]
    if family == 'perm':
        hs = [frozenset((r, p[r]) for r in range(T)) for p in itertools.permutations(range(T))]
    elif family == 'all':
        hs = [frozenset(s) for s in itertools.combinations(cells, T)]
    elif family == 'sample':
        hs = []
        seen = set()
        while len(hs) < limit:
            s = frozenset(cells[k] for k in rng.choice(len(cells), size=T, replace=False))
            if s not in seen:
                seen.add(s); hs.append(s)
    else:
        raise SystemExit(f'unknown family {family}')
    # keep one representative per dihedral class
    keep, seen = [], set()
    for s in hs:
        k = canon(s, T)
        if k not in seen:
            seen.add(k); keep.append(s)
    return keep


def picture(holes, T):
    return '\n'.join('      ' + ''.join('.' if (r, c) in holes else '#' for c in range(T))
                     for r in range(T))


def final_holes(X, Y, T):
    """the hole set of a RELAXED configuration: round each centre to its tiling cell.
    Returns None if the rounding is not injective (the configuration left the tiling)."""
    occ = set()
    for x, y in zip(X, Y):
        r = int(round(y - 0.5)); c = int(round(x - 0.5))
        if not (0 <= r < T and 0 <= c < T) or (r, c) in occ:
            return None
        occ.add((r, c))
    return frozenset((r, c) for r in range(T) for c in range(T) if (r, c) not in occ)


def describe(h, T):
    if h is None:
        return 'off-tiling'
    return f"perm={perm_type(h, T)} rows={shape(h, T)[0]} cols={shape(h, T)[1]} holes={tuple(sorted(h))}"


# --------------------------------------------------------------------------- workers
_F = None


def _init(theta, T):
    global _F
    _F = s6local.Fixed(theta, T)


def _job(arg):
    k, X, Y = arg
    v, bX, bY, _asg = _F.ascent(X, Y)
    return k, v, bX, bY


def starts_for(holes, T, n, rng, jitters, rho):
    cells = [(c + 0.5, r + 0.5) for r in range(T) for c in range(T) if (r, c) not in holes]
    assert len(cells) == n
    out = []
    for _ in range(jitters):
        X = np.array([c[0] for c in cells]) + rng.uniform(-rho, rho, n)
        Y = np.array([c[1] for c in cells]) + rng.uniform(-rho, rho, n)
        out.append((X, Y))
    return out


def cmd_holes(a):
    T = a.T
    n = T * T - T
    rng = np.random.default_rng(a.seed)
    fam = hole_family(T, a.family, rng, a.limit)
    tdegs = [float(s) for s in a.tdeg.split(',')]
    print(f"# arch_ct holes: T={T} n={n} family={a.family} -> {len(fam)} hole sets "
          f"(dihedral classes), jitters={a.jitters} rho={a.rho} nproc={a.nproc} seed={a.seed}")
    sys.stdout.flush()

    carry = {}                                   # hole index -> list of (X, Y) to re-use
    for td in sorted(tdegs, reverse=True):
        t = math.radians(td)
        theta = [t] * n
        jobs = []
        for k, h in enumerate(fam):
            for (X, Y) in starts_for(h, T, n, rng, a.jitters, a.rho):
                jobs.append((k, X, Y))
            for (X, Y) in carry.get(k, []):
                jobs.append((k, X, Y))
        t0 = time.time()
        from multiprocessing import Pool
        best = {}                                # start hole index -> (v, X, Y)
        byfinal = {}                             # canonical FINAL hole set -> (v, X, Y, count)
        with Pool(a.nproc, initializer=_init, initargs=(theta, T)) as pool:
            for k, v, X, Y in pool.imap_unordered(_job, jobs,
                                                  chunksize=max(1, len(jobs) // (a.nproc * 16))):
                if X is None:
                    continue
                if k not in best or v > best[k][0]:
                    best[k] = (v, X, Y)
                fh = final_holes(X, Y, T)
                key = canon(fh, T) if fh is not None else None
                old = byfinal.get(key)
                if old is None or v > old[0]:
                    byfinal[key] = (v, X, Y, (0 if old is None else old[3]) + 1)
                else:
                    byfinal[key] = (old[0], old[1], old[2], old[3] + 1)
        dt = time.time() - t0
        carry = {k: [(X, Y)] for k, (v, X, Y) in best.items()}
        rank = sorted(best.items(), key=lambda kv: -kv[1][0])
        print(f"\n## t = {td} deg   ({len(jobs)} ascents, {dt:.1f} s)")
        vbest = rank[0][1][0]
        print(f"  best delta* >= {vbest:+.12f}   /t^2 = {vbest / t**2:+.6f}   "
              f"(=> c_{T} <= {-vbest / t**2:.6f})")
        print(f"  --- best by FINAL hole set ({len(byfinal)} distinct, dihedral classes) ---")
        for key, (v, X, Y, cnt) in sorted(byfinal.items(), key=lambda kv: -kv[1][0])[:a.topk]:
            h = frozenset(key) if key is not None else None
            print(f"  {v:+.12f}  /t^2={v / t**2:+.6f}  hits={cnt}  {describe(h, T)}")
        print(f"  --- best by STARTING hole set ---")
        bytype = {}
        for k, (v, _X, _Y) in best.items():
            key = (perm_type(fam[k], T), shape(fam[k], T)[0])
            if key not in bytype or v > bytype[key][0]:
                bytype[key] = (v, k)
        for key, (v, k) in sorted(bytype.items(), key=lambda kv: -kv[1][0])[:a.topk]:
            print(f"    {v:+.12f}  /t^2={v / t**2:+.6f}  cycle={key[0]}  rows={key[1]}  "
                  f"n_class={sum(1 for kk in best if (perm_type(fam[kk], T), shape(fam[kk], T)[0]) == key)}")
        if a.show:
            _k, (v, X, Y) = rank[0]
            fh = final_holes(X, Y, T)
            print(f"  winning configuration at t={td} (value {v:+.6e}); FINAL hole set:")
            if fh is not None:
                print(picture(fh, T))
                print(f"      {describe(fh, T)}")
            else:
                print("      (off-tiling)")
            s6local.grid_print(X, Y, theta, T)
        sys.stdout.flush()


# --------------------------------------------------------------------------- §3.2 window
def cmd_window(a):
    """A straight row of T unit squares at a common tilt t has bounding box
       (T cos t + sin t) x (cos t + T sin t).  It fits in [0,T]^2 iff both are <= T."""
    print("# straight-chain window:  T cos t + sin t <= T   <=>   tan(t/2) >= 1/T")
    print("#   T cos t + sin t <= T  <=>  sin t <= T (1 - cos t)")
    print("#                        <=>  2 sin(t/2) cos(t/2) <= 2 T sin^2(t/2)")
    print("#                        <=>  cos(t/2) <= T sin(t/2)   (sin(t/2) > 0 for t in (0, 180))")
    print("#                        <=>  tan(t/2) >= 1/T .                          QED")
    print("# the transverse row  cos t + T sin t <= T  is the same with t -> 90 - t,")
    print("# so the window is  t in [t_T, 90 - t_T],  t_T = 2 arctan(1/T);  empty iff t_T > 45,")
    print("# i.e. iff 1/T > tan 22.5 = sqrt2 - 1, i.e. iff T < 1 + sqrt2 = 2.4142...  (only T = 2).")
    print()
    print("  T     t_T = 2 atan(1/T)   window [t_T, 90-t_T]   T cos t_T + sin t_T   2/T (rad, deg)")
    for T in (2, 3, 4, 5, 6, 7, 10, 17, 20, 100):
        tT = 2 * math.atan(1.0 / T)
        d = math.degrees(tT)
        w = T * math.cos(tT) + math.sin(tT)
        ok = 'empty' if d > 45 else f"[{d:.4f}, {90 - d:.4f}]"
        print(f" {T:3d}      {d:9.4f} deg     {ok:>22s}   {w:.12f}   "
              f"{2.0 / T:.6f} = {math.degrees(2.0 / T):.4f} deg")
    print()
    print("# numeric check of the equivalence on a grid (t in (0, 90) deg):")
    bad = 0
    for T in (2, 3, 4, 5, 6, 17):
        for i in range(1, 9000):
            t = math.radians(i / 100.0)
            lhs = (T * math.cos(t) + math.sin(t) <= T + 1e-12)
            rhs = (math.tan(t / 2) >= 1.0 / T - 1e-12)
            if lhs != rhs and abs(math.tan(t / 2) - 1.0 / T) > 1e-9:
                bad += 1
    print(f"  mismatches over 6 x 8999 samples: {bad}")
    print()
    print("# sanity: at t just above t_T the straight row of T fits, so a *row* is not what")
    print("# forbids a tilt; only stacking T-1 such rows does.  Slack of the transverse extent")
    print("# at t = t_T:  T - (cos t_T + T sin t_T)")
    for T in (3, 4, 5, 17):
        tT = 2 * math.atan(1.0 / T)
        print(f"  T={T:3d}:  transverse extent {math.cos(tT) + T * math.sin(tT):.9f} "
              f"(slack {T - (math.cos(tT) + T * math.sin(tT)):+.9f})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('holes')
    p.add_argument('--T', type=int, required=True)
    p.add_argument('--family', default='perm', choices=('perm', 'all', 'sample'))
    p.add_argument('--tdeg', default='0.8')
    p.add_argument('-j', '--jitters', type=int, default=8)
    p.add_argument('--rho', type=float, default=0.08)
    p.add_argument('--limit', type=int, default=4000)
    p.add_argument('--topk', type=int, default=10)
    p.add_argument('--nproc', type=int, default=8)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--show', action='store_true')
    p.set_defaults(fn=cmd_holes)
    p = sub.add_parser('window')
    p.set_defaults(fn=cmd_window)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
