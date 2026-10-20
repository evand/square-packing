#!/usr/bin/env python3
"""D146 (10-19; Evan: "spend some effort on more diverse near-record 146 packings").  HB's 147 seeds: best seed per
donor tracks structure more than donor side, and the store's sub-12.66 146s are 9,186 / 10,482 the Xu 35-39-tilt family.

  bank:   fingerprint images (pk/fingerprint.py) of every distinct store 146 below --smax (status != infeasible) ->
          runs/d146/bank.npz; known configurations for explore.py --sparse-share.
  run:    explorer runs at n = 146 from the minority tilt bands (runs per band below), 1 proc each, --adapt --adapt-cells
          --sparse-share 0.5 --fp-bank bank.npz, --smax 12.652 (HB's d1 donor ceiling T(147) - 0.005).
  report: archive 146s below --smax per run: distinct sides, tilt band, fingerprint distance to the bank (nnd); writes
          runs/d146/donors.json (the novel ones, for hborder.py seeds --donor-file).

  d146.py bank; d146.py run [--minutes 40]; d146.py report
"""
import argparse, collections, glob, json, os, random, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = f'{HERE}/runs/d146'
REC = 12.583782190152
BANDS = {'t40': ((40, 45), 6), 't30': ((30, 35), 6), 'tedge': ((20, 30), 2), 'thigh': ((45, 55), 2)}   # band: (tilt range, runs)


def ntilt(sq):
    a = np.mod(np.asarray([t for _, _, t in sq], float), 90.0)
    return int((np.minimum(a, 90 - a) > 3.0).sum())


def distinct(st, smax):
    rows = st.db.execute("SELECT id, s, ntilted FROM packing WHERE n=146 AND s < ? AND status != 'infeasible' ORDER BY s",
                         (smax,)).fetchall()
    seen, D = set(), []
    for pid, s, t in rows:
        if round(s, 7) not in seen:
            seen.add(round(s, 7)); D.append((pid, s, t or 0))
    return D


def bank(a):
    from pk.store import Store
    from pk import fingerprint as FP
    import mcmin
    os.makedirs(f'{ROOT}/tmp', exist_ok=True)
    st = Store()
    D = distinct(st, a.smax)
    ims, ids = [], []
    for pid, s, t in D:
        p = f'{ROOT}/tmp/b.txt'; st.get(pid).write(p)
        s2, sq = mcmin.load_deg(p)
        ims.append(FP.image(s2, sq)); ids.append(pid)
    np.savez_compressed(f'{ROOT}/bank.npz', images=np.stack(ims).astype(np.float32), ids=np.array(ids))
    c = collections.Counter(t // 5 * 5 for _, _, t in D)
    print(f'bank: {len(ims)} distinct 146s below {a.smax}; by tilt {sorted(c.items())}')


def run(a):
    from pk.store import Store
    from pk import manifest
    import scratch
    scratch.scope('d146')
    st = Store(); rng = random.Random(a.seed)
    D = distinct(st, a.smax)
    nruns = sum(r for _, r in BANDS.values())
    man = manifest.start(ROOT, 'd146-run', seed=a.seed, args=vars(a), procs=nruns, fname='manifest_run.json')
    os.makedirs(f'{ROOT}/starts', exist_ok=True)
    P, i = [], 0
    for band, ((lo, hi), runs) in BANDS.items():
        B = [d for d in D if lo <= d[2] < hi][:runs * a.starts_per_run]
        print(f'{band}: {len(B)} starts ({B[0][1]:.6f}..{B[-1][1]:.6f})' if B else f'{band}: none', flush=True)
        for r in range(runs):
            R = B[r::runs]
            if not R:
                continue
            paths = []
            for pid, s, t in R:
                p = f'{ROOT}/starts/{pid}.txt'; st.get(pid).write(p); paths.append(p)
            out = f'{ROOT}/x{i}_{band}'; i += 1
            os.makedirs(out, exist_ok=True)
            json.dump([dict(id=pid, s=s, tilt=t) for pid, s, t in R], open(f'{out}/starts.json', 'w'))
            cmd = [sys.executable, f'{HERE}/explore.py', '--n', '146', '--starts', *paths, '--minutes', str(a.minutes),
                   '--procs', '1', '--seed', str(a.seed + i), '--adapt', '--adapt-cells', '--sparse-share', '0.5',
                   '--fp-bank', f'{ROOT}/bank.npz', '--smax', str(a.smax), '--frontier-s', str(REC),
                   '--final-polish', '8', '--out', out]
            P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log.txt', 'w'), stderr=subprocess.STDOUT))
    for p in P:
        p.wait()
    man.finish()
    report(a)


def report(a):
    from pk import fingerprint as FP
    import mcmin
    bk = np.load(f'{ROOT}/bank.npz')['images']
    # nnd scale: bank-internal nearest distances of a random sample
    rng = np.random.default_rng(0)
    smp = rng.choice(len(bk), size=min(200, len(bk)), replace=False)
    ref = np.array([np.partition(FP.nearest(bk[j], bk, 146), 1)[1] for j in smp])
    med = float(np.median(ref))
    seen, rows = set(), []
    for d in sorted(glob.glob(f'{ROOT}/x*')):
        for f in glob.glob(f'{d}/*.txt'):
            try:
                h = open(f).readline().split()
                if len(h) != 2 or h[0] != '146':
                    continue
                s = float(h[1])
            except Exception:
                continue
            if s >= a.smax or round(s, 7) in seen:
                continue
            seen.add(round(s, 7))
            s2, sq = mcmin.load_deg(f)
            nnd = float(FP.nearest(FP.image(s2, sq), bk, 146).min())
            rows.append(dict(path=f, s=s, tilt=ntilt(sq), nnd=nnd, run=os.path.basename(d)))
    print(f'bank internal nnd median {med:.4f}; archive 146s below {a.smax}, distinct sides: {len(rows)}')
    by = collections.defaultdict(list)
    for r in rows:
        by[r['tilt'] // 5 * 5].append(r)
    for t in sorted(by):
        R = by[t]
        print(f'  tilt {t:3d}-{t + 4}: {len(R):5d}, < 12.62: {sum(r["s"] < 12.62 for r in R):4d}, < 12.60: '
              f'{sum(r["s"] < 12.60 for r in R):4d}, novel (nnd > median) {sum(r["nnd"] > med for r in R):4d}; best '
              f'{min(r["s"] for r in R):.6f}')
    nov = sorted((r for r in rows if r['nnd'] > med), key=lambda r: r['s'])
    json.dump(dict(median_nnd=med, donors=nov), open(f'{ROOT}/donors.json', 'w'), indent=0)
    print(f'novel donors (nnd > {med:.4f}): {len(nov)}; best {[round(r["s"], 6) for r in nov[:5]]} -> {ROOT}/donors.json')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['bank', 'run', 'report'])
    ap.add_argument('--smax', type=float, default=12.652); ap.add_argument('--minutes', type=float, default=40)
    ap.add_argument('--starts-per-run', type=int, default=6); ap.add_argument('--seed', type=int, default=146)
    a = ap.parse_args()
    {'bank': bank, 'run': run, 'report': report}[a.cmd](a)
