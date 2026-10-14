#!/usr/bin/env python3
"""TT2 census analysis (pre-registered in PACKER.md "TT2 pre-registration", census bullet, 10-14 ~12:30, before scoring).

Reference per n: store rows added before T (as-of, local time), sub-k (s < ceil(sqrt n)), grid-free, feasible.
New (per arm) = a polished, grid-free, sub-k row sighted in that arm's rounds (`run:runs/tt2/<arm>/...`) whose side is
> 2e-9 (pk.cluster.SAME) from every reference side; distinct new minima = distinct side keys.
  (1) distinct new sub-k minima
  (2) new regions: new minima that start their own leader ball (pk.cluster.Leaders, r = 0.42 / sqrt n) against the
      reference leaders (built here on the FULL as-of reference, deepest first, uncapped; then the arm's new minima
      deepest first)
  (3) near-frontier new minima: within 1e-3 of frontier_at(n, T)
  (4) S / C overlap: on (1) side keys within 2e-9; on (2) S's new leaders that fall in a ball of a C new leader (and
      the reverse)
  per measured CPU-h (run manifests' cpu_s).
Certification sample: 40 random new minima per arm, sub-k and within 1e-2 of frontier_at(n, T), through exactsolve
(pk/identity verdicts) for the false-new rate.

  python3 analysis/tt2_census.py census [--procs 15]      -> runs/tt2/census.json
  python3 analysis/tt2_census.py certify [--procs 14]     -> runs/tt2/census_cert/ + census_cert.json
  python3 analysis/tt2_census.py report
"""
from __future__ import annotations
import argparse, bisect, copy, glob, json, math, os, random, sqlite3, sys
from multiprocessing import Pool
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):    # one thread per worker (process-parallel)
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
from pk import cluster as CL, eventbench as EB, parents as PA      # noqa: E402

T = '2026-10-10T06:00'
ARMS = ('S', 'C')
OUT = os.path.join(PK, 'runs', 'tt2')
DB = os.path.join(PK, 'runs', 'store.sqlite')
NEAR, SAMPLE_GAP, SAMPLE_N = 1e-3, 1e-2, 40


def lookup(L, s, sq):
    """index of the leader whose ball holds (s, sq) (as Leaders.assign), else -1; L is not changed"""
    if not L.ids:
        return -1
    im = CL.FP.image(s, sq)
    d = CL.FP.nearest(im, L.bank, L.n)
    cand = np.nonzero(d < 2 * L.r)[0]
    if not len(cand):
        return -1
    cand = cand[np.argsort(d[cand])][:8]
    dm, j = min((CL.md(L.sq[j], L.s[j], sq, s), j) for j in cand)
    return int(j) if dm <= L.r else -1


def rows(db, q, a):
    return [(i, s, np.frombuffer(c, '<f8').reshape(-1, 3)) for i, s, c in db.execute(q, a)]


def dedupe(rs):
    """one row per side key (2e-9), deepest first"""
    out, keys = [], []
    for r in sorted(rs, key=lambda r: (r[1], r[0])):
        i = bisect.bisect_left(keys, r[1] - CL.SAME)
        if i < len(keys) and abs(keys[i] - r[1]) < CL.SAME:
            continue
        keys.insert(i, r[1]); out.append(r)
    return out


def near_key(sorted_sides, s):
    i = bisect.bisect_left(sorted_sides, s - CL.SAME)
    return i < len(sorted_sides) and abs(sorted_sides[i] - s) < CL.SAME


def one_n(n):
    db = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    k = math.ceil(math.sqrt(n) - 1e-12)
    asof = PA.asof_local(T)
    fT = EB.frontier_at(n, T)
    ref = dedupe(rows(db, "SELECT id, s, coords FROM packing WHERE n=? AND s<? AND status!='infeasible' AND added<?"
                          " AND (grid IS NULL OR grid=0)", (n, k, asof)))
    rsides = [r[1] for r in ref]
    L = CL.Leaders(n)
    for pid, s, sq in ref:
        L.assign(pid, s, sq)
    nref = len(L.ids)
    res = dict(n=n, k=k, frontier_T=fT, ref_keys=len(ref), ref_leaders=nref)
    new, Ls = {}, {}
    for a in ARMS:
        got = rows(db, "SELECT DISTINCT p.id, p.s, p.coords FROM packing p JOIN sighting g ON g.packing=p.id WHERE p.n=?"
                       " AND p.s<? AND p.status IN ('polished','certified') AND (p.grid IS NULL OR p.grid=0)"
                       f" AND g.source LIKE 'run:runs/tt2/{a}/%'", (n, k))
        nw = [r for r in dedupe(got) if not near_key(rsides, r[1])]
        La = copy.deepcopy(L)
        reg = []
        for pid, s, sq in nw:
            j = La.assign(pid, s, sq)
            if j >= nref and La.ids[j] == pid:
                reg.append(j)
        new[a], Ls[a] = nw, (La, reg)
        res[a] = dict(arm_keys=len(dedupe(got)), new=len(nw), regions=len(reg),
                      near=sum(1 for r in nw if r[1] - fT < NEAR),
                      sample_pool=[(r[0], r[1], r[1] - fT) for r in nw if r[1] - fT < SAMPLE_GAP],
                      best_new_gap=min((r[1] - fT for r in nw), default=None))
    sS = sorted(r[1] for r in new['S'])
    res['overlap_new'] = sum(1 for r in new['C'] if near_key(sS, r[1]))
    for a, b in (('S', 'C'), ('C', 'S')):         # a's new regions that fall in a ball of one of b's new leaders
        Lb, regb = Ls[b]
        La, rega = Ls[a]
        hit = 0
        regb = set(regb)
        for j in rega:
            hit += lookup(Lb, La.s[j], La.sq[j]) in regb
        res[f'overlap_regions_{a}in{b}'] = hit
    db.close()
    return res


def cpu_h(a):
    c = 0.0
    for m in glob.glob(os.path.join(OUT, a, '**', 'manifest.json'), recursive=True):
        c += json.load(open(m)).get('cpu_s') or 0
    return c / 3600


def census(procs):
    ns = json.load(open(os.path.join(PK, 'runs', 'night1014_ns.json')))
    ns = ns if isinstance(ns, list) else ns['ns']
    ns = sorted(ns, key=lambda n: -n)                 # big n first (slowest)
    R = []
    with Pool(procs) as P:
        for r in P.imap_unordered(one_n, ns):
            R.append(r)
            print(f"n={r['n']:3d} ref {r['ref_keys']:5d}/{r['ref_leaders']:4d}  " + '  '.join(
                f"{a}: new {r[a]['new']:4d} reg {r[a]['regions']:3d} near {r[a]['near']:3d}" for a in ARMS), flush=True)
    R.sort(key=lambda r: r['n'])
    json.dump(dict(T=T, cpu_h={a: cpu_h(a) for a in ARMS}, rows=R), open(os.path.join(OUT, 'census.json'), 'w'))


def _solve(args):
    pid, path, outdir = args
    import contextlib, io
    sys.path.insert(0, os.path.join(PK, '..', 'exact'))
    import census_exact as CE
    import exactsolve as E
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rep = E.run(path, outdir=outdir, quiet=True)[0]
        return pid, CE.classify(rep)
    except Exception as e:                            # noqa: BLE001
        return pid, dict(cls='unresolved', status=f'error: {e}'[:120])


def certify(procs, seed=0):
    from pk.store import Store
    D = json.load(open(os.path.join(OUT, 'census.json')))
    wd = os.path.join(OUT, 'census_cert'); os.makedirs(wd, exist_ok=True)
    st = Store(); jobs, meta = [], {}
    for a in ARMS:
        pool = [(r['n'], *x) for r in D['rows'] for x in r[a]['sample_pool']]
        pick = random.Random(seed).sample(pool, min(SAMPLE_N, len(pool)))
        for n, pid, s, gap in pick:
            path = os.path.join(wd, f'{a}_n{n}_{pid}.txt')
            if not os.path.exists(path):
                st.get(pid).write(path)
            jobs.append((pid, path, os.path.join(wd, 'exact')))
            meta.setdefault(pid, []).append(dict(arm=a, n=n, s=s, gap=gap))
    st.close()
    os.makedirs(os.path.join(wd, 'exact'), exist_ok=True)
    out = {}
    with Pool(procs) as P:
        for pid, v in P.imap_unordered(_solve, jobs):
            out[pid] = v
            print(pid, meta[pid][0]['arm'], meta[pid][0]['n'], v['cls'], v.get('status', '')[:60], flush=True)
    json.dump(dict(meta={str(k): v for k, v in meta.items()}, verdicts={str(k): v for k, v in out.items()}),
              open(os.path.join(OUT, 'census_cert.json'), 'w'))


def report():
    D = json.load(open(os.path.join(OUT, 'census.json')))
    R, H = D['rows'], D['cpu_h']
    print(f"cutoff {D['T']}; measured CPU-h S {H['S']:.1f}, C {H['C']:.1f}; {len(R)} n")
    tot = {a: {m: sum(r[a][m] for r in R) for m in ('new', 'regions', 'near')} for a in ARMS}
    for m, lab in (('new', '(1) distinct new sub-k minima'), ('regions', '(2) new regions'), ('near', '(3) new within 1e-3')):
        s, c = tot['S'][m], tot['C'][m]
        print(f"{lab:32s} S {s:6d} ({s / H['S']:7.2f}/CPU-h)  C {c:6d} ({c / H['C']:7.2f}/CPU-h)  ratio S/C per CPU-h "
              f"{(s / H['S']) / (c / H['C']) if c else float('nan'):.2f}")
    print(f"(4) overlap: new side keys in both {sum(r['overlap_new'] for r in R)}; S regions inside a C new region "
          f"{sum(r['overlap_regions_SinC'] for r in R)}, C inside S {sum(r['overlap_regions_CinS'] for r in R)}")
    for a in ARMS:
        nn = sum(1 for r in R if r[a]['regions'])
        print(f"  {a}: n with >= 1 new region {nn}; n with >= 1 near-frontier new minimum {sum(1 for r in R if r[a]['near'])}")
    p = os.path.join(OUT, 'census_cert.json')
    if os.path.exists(p):
        C = json.load(open(p))
        from collections import Counter
        for a in ARMS:
            v = Counter(C['verdicts'][pid]['cls'] for pid, ms in C['meta'].items() for m in ms if m['arm'] == a
                        and pid in C['verdicts'])
            print(f"  certification sample {a}: {dict(v)}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['census', 'certify', 'report'])
    ap.add_argument('--procs', type=int, default=15)
    a = ap.parse_args()
    {'census': lambda: census(a.procs), 'certify': lambda: certify(a.procs), 'report': report}[a.cmd]()
