#!/usr/bin/env python3
"""Diverse near-record donors at any n (10-19; Evan: "spend more effort generating additional starting packings, even if
they're not new records themselves").  d146.py generalized: starts by fingerprint farthest-point sampling, not tilt band
(308's pool is all 100-119 tilted).

  bank:   fingerprints (pk/fingerprint.py) of every distinct store packing at n below --smax (status != infeasible)
          -> runs/dgen/n<n>/bank.npz.
  run:    --runs explorer runs x 1 proc x --minutes; starts = greedy farthest-point sample (fingerprint distance, seeded
          with the record) over the best --pool distinct, --starts-per-run each, round-robin; --adapt --adapt-cells
          --sparse-share 0.5 --fp-bank, --smax.
  report: archive packings below --smax: distinct sides, novel (nnd to the bank > bank median nnd) -> donors.json
          (path, s, tilt, nnd; best side first), consumed by hborder.py --donor-file / lcut.py --donor-files.

  dgen.py bank --n 308 --smax 18; dgen.py run --n 308 --smax 18 [--runs 8 --minutes 40]; dgen.py report --n 308 --smax 18
"""
import argparse, collections, glob, json, os, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def root(a):
    return f'{HERE}/runs/dgen/n{a.n}{getattr(a, "tag", "")}'


def ntilt(sq):
    t = np.mod(np.asarray([x[2] for x in sq], float), 90.0)
    return int((np.minimum(t, 90 - t) > 3.0).sum())


def distinct(st, n, smax):
    rows = st.db.execute("SELECT id, s, ntilted FROM packing WHERE n=? AND s < ? AND status != 'infeasible' ORDER BY s",
                         (n, smax)).fetchall()
    seen, D = set(), []
    for pid, s, t in rows:
        if round(s, 7) not in seen:
            seen.add(round(s, 7)); D.append((pid, s, t or 0))
    return D


def images(st, D, tmp):
    from pk import fingerprint as FP
    import mcmin
    ims = []
    for pid, s, t in D:
        st.get(pid).write(tmp)
        s2, sq = mcmin.load_deg(tmp)
        ims.append(FP.image(s2, sq))
    return np.stack(ims).astype(np.float32)


def bank(a):
    from pk.store import Store
    R0 = root(a); os.makedirs(R0, exist_ok=True)
    st = Store()
    D = distinct(st, a.n, a.smax)
    ims = images(st, D, f'{R0}/tmp.txt')
    np.savez_compressed(f'{R0}/bank.npz', images=ims, ids=np.array([d[0] for d in D]), s=np.array([d[1] for d in D]))
    c = collections.Counter(t // 5 * 5 for _, _, t in D)
    print(f'bank n {a.n}: {len(D)} distinct below {a.smax}; by tilt {sorted(c.items())}')


def fps(ims, n, m):
    """greedy farthest-point sample of m indices (fingerprint distance), seeded with index 0 (the best side)"""
    from pk import fingerprint as FP
    pick = [0]
    d = FP.nearest(ims[0], ims, n).astype(float)
    while len(pick) < min(m, len(ims)):
        j = int(np.argmax(d)); pick.append(j)
        d = np.minimum(d, FP.nearest(ims[j], ims, n))
    return pick


def run(a):
    from pk.store import Store
    from pk import manifest
    import scratch
    R0 = root(a)
    scratch.scope(f'dgen{a.n}')
    st = Store()
    bk = np.load(f'{R0}/bank.npz')
    ids, sides = list(bk['ids']), list(bk['s'])
    P0 = min(a.pool, len(ids))
    if a.band:                      # SW (10-19): starts from [F + lo, F + hi), record excluded, farthest from it first
        from pk import fingerprint as FP
        F = frontier(st, a.n)
        lo, hi = map(float, a.band.split(','))
        C = [j for j in range(len(ids)) if F + lo <= sides[j] < F + hi][:a.pool]
        ims = bk['images']
        d = FP.nearest(ims[0], ims[C], a.n).astype(float)    # bank index 0 = the best stored packing (the record)
        sel = []
        while len(sel) < min(a.runs * a.starts_per_run, len(C)):
            j = int(np.argmax(d)); sel.append(C[j]); d[j] = -1
            d = np.minimum(d, FP.nearest(ims[C[j]], ims[C], a.n)); d[[C.index(x) for x in sel]] = -1
        print(f'n {a.n}: frontier {F:.9f}; band {lo:g}..{hi:g}: {len(C)} distinct; {len(sel)} starts', flush=True)
    else:
        sel = fps(bk['images'][:P0], a.n, a.runs * a.starts_per_run)
    man = manifest.start(R0, f'dgen-run-{a.n}', seed=a.seed, args=vars(a), procs=a.runs, fname='manifest_run.json')
    os.makedirs(f'{R0}/starts', exist_ok=True)
    print(f'n {a.n}: {len(sel)} farthest-point starts from the best {P0} (sides {min(sides[j] for j in sel):.6f}..'
          f'{max(sides[j] for j in sel):.6f})', flush=True)
    P = []
    for r in range(a.runs):
        S = sel[r::a.runs]
        paths = []
        for j in S:
            p = f'{R0}/starts/{ids[j]}.txt'; st.get(int(ids[j])).write(p); paths.append(p)
        out = f'{R0}/x{r}'
        os.makedirs(out, exist_ok=True)
        json.dump([dict(id=int(ids[j]), s=float(sides[j])) for j in S], open(f'{out}/starts.json', 'w'))
        cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(a.n), '--starts', *paths,
               *(['--cpu-h', str(a.cpu_h / a.runs)] if a.cpu_h else ['--minutes', str(a.minutes)]),
               '--procs', '1', '--seed', str(a.seed + r), '--adapt', '--adapt-cells', '--sparse-share', '0.5',
               '--fp-bank', f'{R0}/bank.npz', '--smax', str(a.smax), '--frontier-s', str(sides[0]),
               '--final-polish', '8', '--out', out]
        P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log.txt', 'w'), stderr=subprocess.STDOUT))
    for p in P:
        p.wait()
    man.finish()
    report(a)


def frontier(st, n):
    return st.db.execute('SELECT min(s) FROM frontier WHERE n=?', (n,)).fetchone()[0]


def sweep(a):
    """SW: bank + band run + record check for one n (smax = frontier + band hi + 0.002)"""
    from pk.store import Store
    F = frontier(Store(), a.n)
    a.smax = F + float(a.band.split(',')[1]) + 0.002
    bank(a); run(a)
    best = []
    for r in range(a.runs):
        for f in glob.glob(f'{root(a)}/x{r}/*.txt'):
            h = open(f).readline().split()
            if len(h) == 2 and h[0] == str(a.n):
                best.append((float(h[1]), f))
    best.sort()
    print(f'SW n {a.n}: frontier {F:.9f}; best {best[0][0]:.9f} ({best[0][0] - F:+.3e}) {best[0][1]}' if best else 'none')
    print('SW-RECORD-CANDIDATE' if best and best[0][0] < F - 1e-9 else 'SW-no-record', flush=True)


def report(a):
    from pk import fingerprint as FP
    import mcmin
    R0 = root(a)
    bk = np.load(f'{R0}/bank.npz')['images']
    rng = np.random.default_rng(0)
    smp = rng.choice(len(bk), size=min(200, len(bk)), replace=False)
    med = float(np.median([np.partition(FP.nearest(bk[j], bk, a.n), 1)[1] for j in smp]))
    seen, rows = set(), []
    for d in sorted(glob.glob(f'{R0}/x*')):
        for f in glob.glob(f'{d}/*.txt'):
            try:
                h = open(f).readline().split()
                if len(h) != 2 or h[0] != str(a.n):
                    continue
                s = float(h[1])
            except Exception:
                continue
            if s >= a.smax or round(s, 7) in seen:
                continue
            seen.add(round(s, 7))
            s2, sq = mcmin.load_deg(f)
            rows.append(dict(path=f, s=s, tilt=ntilt(sq), nnd=float(FP.nearest(FP.image(s2, sq), bk, a.n).min()),
                             run=os.path.basename(d)))
    nov = sorted((r for r in rows if r['nnd'] > med), key=lambda r: r['s'])
    print(f'n {a.n}: bank median nnd {med:.4f}; archive distinct below {a.smax}: {len(rows)}, novel {len(nov)}; '
          f'novel best {[round(r["s"], 6) for r in nov[:5]]}')
    q = np.quantile([r['s'] for r in rows], [0.1, 0.5]) if rows else []
    print(f'  archive side 10 / 50 % quantiles {np.round(q, 6).tolist()}; novel by run '
          f'{sorted(collections.Counter(r["run"] for r in nov).items())}')
    json.dump(dict(n=a.n, median_nnd=med, donors=nov), open(f'{R0}/donors.json', 'w'), indent=0)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['bank', 'run', 'report', 'sweep'])
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--smax', type=float, default=None)
    ap.add_argument('--band', help='SW: start band above the store frontier, "lo,hi" (record excluded)')
    ap.add_argument('--tag', default='', help='output dir suffix: runs/dgen/n<n><tag>')
    ap.add_argument('--runs', type=int, default=8); ap.add_argument('--minutes', type=float, default=40)
    ap.add_argument('--cpu-h', type=float, default=None, help='total CPU budget, split evenly over --runs (replaces --minutes)')
    ap.add_argument('--pool', type=int, default=400, help='farthest-point sample over the best POOL distinct')
    ap.add_argument('--starts-per-run', type=int, default=6); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    a.seed = a.seed or a.n
    {'bank': bank, 'run': run, 'report': report, 'sweep': sweep}[a.cmd](a)
