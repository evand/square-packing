#!/usr/bin/env python3
"""Insertion test at 110 (10-10 evening; Evan): do 110 packings built by inserting a square into non-regular 109 packings
add a pile of new sub-11 minima (the "every new technique adds a pile of sub-11 packings" theory), and do any lie outside
the record's family (census: all 754 certified minima within matching distance 0.104 of the record)?

  seeds:   control = 20 certified 110 census minima -> one square out (chain.seed_job k = 1, must shrink) -> one square
           back in (k = -1: dilate 0.2-1.2 %, one of the 12 best pose holes, fq quench); donors = distinct non-regular
           sub-11 109 packings from the store, stratified by side and tilted count, 3 insertion seeds each.
  explore: seeds that quench below --smax-seed go to explorer runs (novelty reward against every known sub-11 110 side,
           no elite share).
  report:  quench outcomes per arm; explorer sub-11 entries -> census_certify / descend_census -> new distinct minima and
           their matching distance from the record (pk.crossing.match_distance vs store id 40600).

  insert110.py seeds   [--donors 40 --per 3 --control 20 --tries 4 --procs 10]
  insert110.py explore [--runs 8 --procs-per 2 --minutes 20]
  insert110.py report
"""
import argparse, collections, glob, json, math, os, random, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = f'{HERE}/runs/insert110'
REC = 40600
REG109 = 6 + 7 / math.sqrt(2)                     # best(109): the regular 45-degree construction


def pick_control(st, m, rng):
    rows = st.db.execute("SELECT key, S FROM basin WHERE n=110 AND cls='certified' AND S IS NOT NULL").fetchall()
    rows = sorted(rows, key=lambda r: float(r[0]))
    step = max(1, len(rows) // m)
    out = []
    for key, S in rows[::step][:m]:
        s = float(key)
        r = st.db.execute('SELECT id FROM packing WHERE n=110 AND s BETWEEN ? AND ? ORDER BY pen LIMIT 1', (s - 5e-10, s + 5e-10)).fetchone()
        if r:
            out.append(r[0])
    return out


def pick_donors(st, m, rng):
    rows = st.db.execute('SELECT id, s, ntilted FROM packing WHERE n=109 AND s < 11 ORDER BY s').fetchall()
    seen, D = set(), []
    for pid, s, t in rows:
        k = round(s, 7)
        if k in seen or abs(s - REG109) < 1e-6:
            continue
        seen.add(k); D.append((pid, s, t))
    strata = collections.defaultdict(list)
    lo, hi = D[0][1], 11.0
    for d in D:
        strata[(min(4, int(5 * (d[1] - lo) / (hi - lo))), d[2] // 2)].append(d)
    keys = sorted(strata)
    out = []
    while len(out) < m and any(strata.values()):          # round-robin over strata
        for k in keys:
            if strata[k] and len(out) < m:
                out.append(strata[k].pop(rng.randrange(len(strata[k]))))
    return out


def _seed(j):
    import chain
    chain.TRIES = j[-1]
    return chain.seed_job(j[:-1])


def seeds(a):
    from pk.store import Store
    from jobpool import run_jobs
    st = Store(); rng = random.Random(a.seed)
    os.makedirs(f'{ROOT}/control', exist_ok=True); os.makedirs(f'{ROOT}/donors', exist_ok=True); os.makedirs(f'{ROOT}/seeds', exist_ok=True)
    ctrl = pick_control(st, a.control, rng)
    donors = pick_donors(st, a.donors, rng)
    for pid in ctrl:
        st.get(pid).write(f'{ROOT}/control/c{pid}.txt')
    for pid, s, t in donors:
        st.get(pid).write(f'{ROOT}/donors/d{pid}.txt')
    t0 = time.time()
    # control: 110 -> 109 (removal)
    J = [(f'{ROOT}/control/c{pid}.txt', 1, 109, f'{ROOT}/control/r{pid}.txt', rng.randrange(10 ** 9), a.tries) for pid in ctrl]
    R = run_jobs(_seed, J, procs=a.procs, timeout=1200)
    rem = [(pid, r) for pid, r in zip(ctrl, R) if isinstance(r, float)]
    print(f'control removals: {len(rem)}/{len(ctrl)} shrank (109 sides {min(r for _, r in rem):.6f} .. {max(r for _, r in rem):.6f}) '
          f'[{time.time() - t0:.0f}s]', flush=True)
    # insertions: control (one seed each, as many as donors' per) and donors
    J, meta = [], []
    for pid, _ in rem:
        for q in range(a.per):
            J.append((f'{ROOT}/control/r{pid}.txt', -1, 110, f'{ROOT}/seeds/c{pid}_{q}.txt', rng.randrange(10 ** 9), a.tries))
            meta.append(dict(arm='control', src=pid, q=q))
    for pid, s, t in donors:
        for q in range(a.per):
            J.append((f'{ROOT}/donors/d{pid}.txt', -1, 110, f'{ROOT}/seeds/d{pid}_{q}.txt', rng.randrange(10 ** 9), a.tries))
            meta.append(dict(arm='donor', src=pid, src_s=s, src_tilt=t, q=q))
    R = run_jobs(_seed, J, procs=a.procs, timeout=1200)
    with open(f'{ROOT}/seeds.jsonl', 'w') as fo:
        for j, m, r in zip(J, meta, R):
            m.update(path=j[3], s=r if isinstance(r, float) else None)
            fo.write(json.dumps(m) + '\n')
    for arm in ('control', 'donor'):
        S = [m['s'] for m in meta if m['arm'] == arm]
        ok = [s for s in S if s is not None]
        print(f'{arm}: {len(S)} seeds, {len(ok)} feasible; below 11: {sum(s < 11 - 1e-9 for s in ok)}, '
              f'exactly 11 (grid): {sum(abs(s - 11) < 1e-9 for s in ok)}, 11-11.02: {sum(11 + 1e-9 < s < 11.02 for s in ok)}, '
              f'above: {sum(s >= 11.02 for s in ok)}; min {min(ok):.7f}' if ok else f'{arm}: none feasible', flush=True)
    print(f'seeds done [{time.time() - t0:.0f}s]')


def explore(a):
    from pk.store import Store
    st = Store()
    M = [json.loads(l) for l in open(f'{ROOT}/seeds.jsonl')]
    ref = sorted(set(float(k) for (k,) in st.db.execute('SELECT key FROM basin WHERE n=110')) |
                 set(s for (s,) in st.db.execute('SELECT s FROM packing WHERE n=110 AND s < 11')))
    json.dump(ref, open(f'{ROOT}/novel_ref.json', 'w'))
    use = [m for m in M if m['s'] is not None and m['s'] < a.smax_seed and abs(m['s'] - 11) > 1e-9]
    by = collections.defaultdict(list)                    # group by source so a run's starts share no lineage
    for m in use:
        by[(m['arm'], m['src'])].append(m)
    groups = sorted(by.values(), key=lambda g: (g[0]['arm'], g[0]['src']))
    ctrl = [g for g in groups if g[0]['arm'] == 'control']; don = [g for g in groups if g[0]['arm'] == 'donor']
    # baseline (old technique): explorer from the control's source census minima, same flags and minutes
    base = sorted(glob.glob(f'{ROOT}/control/c*.txt'))
    runs = [[dict(arm='baseline', src=os.path.basename(p), path=p, s=float(open(p).readline().split()[1]))
             for p in base[i::a.baseline]] for i in range(a.baseline)]
    runs += [[] for _ in range(a.control_runs)]
    for i, g in enumerate(ctrl):
        runs[a.baseline + i % a.control_runs].extend(g)
    nd = a.runs - a.baseline - a.control_runs
    runs += [[] for _ in range(nd)]
    for i, g in enumerate(don):
        runs[a.baseline + a.control_runs + i % nd].extend(g)
    P = []
    for i, R in enumerate(runs):
        if not R:
            continue
        out = f'{ROOT}/x{i}_{R[0]["arm"]}'
        os.makedirs(out, exist_ok=True)
        json.dump(R, open(f'{out}/starts.json', 'w'))
        cmd = [sys.executable, f'{HERE}/explore.py', '--n', '110', '--starts', *[m['path'] for m in R], '--minutes', str(a.minutes),
               '--procs', str(a.procs_per), '--seed', str(500 + i), '--adapt', '--reward', 'legacy', '--novel-ref',
               f'{ROOT}/novel_ref.json', '--out', out]
        P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
        print(f'x{i} ({R[0]["arm"]}): {len(R)} starts, sides {min(m["s"] for m in R):.5f} .. {max(m["s"] for m in R):.5f}', flush=True)
    for p in P:
        p.wait()
    print('explore done')


def report(a):
    """Per arm: sub-11 sides (excluding starts) -> exactsolve class; equal random sample of not-min per arm descended;
    new = exact S not in known_before.json (census frozen before the explorers); d = matching distance from the record."""
    from pk.store import Store
    from pk.packing import read
    from pk.crossing import match_distance
    from jobpool import run_jobs
    import census_certify as CC, descend_census as DC
    st = Store(); rec = st.get(REC)
    known = set(json.load(open(f'{ROOT}/known_before.json')))
    if not a.no_solve:
        subprocess.run([sys.executable, f'{HERE}/pk.py', 'import-run', *sorted(glob.glob(f'{ROOT}/x*'))], cwd=HERE)
        CC.certify(argparse.Namespace(n=110, k=11, procs=a.procs, timeout=900, limit=None))
    B = {k: (c, S) for k, c, S in st.db.execute('SELECT key, cls, S FROM basin WHERE n=110')}
    rng = random.Random(1)
    arms = collections.defaultdict(lambda: dict(runs=0, ph=0.0, keys={}))
    for d in sorted(glob.glob(f'{ROOT}/x*')):
        arm = d.rsplit('_', 1)[1]
        E = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
        A = arms[arm]; A['runs'] += 1; A['ph'] += a.minutes / 60 * a.procs_per
        for e in E:
            if e['s'] < 11 and e['parent'] >= 0:
                A['keys'].setdefault(f"{e['s']:.9f}", e['path'])
    root = f'{ROOT}/descend'; os.makedirs(root, exist_ok=True)
    jobs, owner = [], []
    for arm, A in sorted(arms.items()):
        nm = sorted(k for k in A['keys'] if B.get(k, ('?',))[0] == 'not-min')
        A['nm'] = len(nm)
        A['sample'] = rng.sample(nm, min(a.sample, len(nm)))
        for k in A['sample']:
            r = st.db.execute('SELECT id FROM packing WHERE n=110 AND s BETWEEN ? AND ? ORDER BY pen LIMIT 1',
                              (float(k) - 5e-10, float(k) + 5e-10)).fetchone()
            if r:
                jobs.append((k, r[0], 110, root, 1e-5, 20, 1e-9)); owner.append(arm)
    cache = f'{ROOT}/descend.json'
    D = json.load(open(cache)) if os.path.exists(cache) else {}
    todo = [(j, o) for j, o in zip(jobs, owner) if j[0] not in D]
    if todo and not a.no_solve:
        R = run_jobs(DC.job, [j for j, _ in todo], procs=a.procs, timeout=900)
        for (j, _), r in zip(todo, R):
            D[j[0]] = r or dict(cls='unresolved', status='timeout')
        json.dump(D, open(cache, 'w'))
    def dist(path):
        return match_distance(read(path), rec)
    print(f"{'arm':9s} {'proc-h':>6s} {'sub11':>5s} {'cert':>4s} {'new':>4s} {'new d>0.104':>11s} | {'not-min':>7s} {'sampled':>7s} "
          f"{'->new':>5s} {'->known':>7s} {'still':>5s} | {'est new/proc-h':>14s}")
    for arm, A in sorted(arms.items()):
        cert = [k for k in A['keys'] if B.get(k, ('?',))[0] == 'certified']
        newS = {}
        for k in cert:
            S = (B[k][1] or '')[:23]
            if S and S not in known:
                newS.setdefault(S, A['keys'][k])
        far = sum(dist(p) > 0.104 for p in newS.values())
        smp = [D[k] for k in A['sample'] if k in D]
        to_new = {r['S'][:23] for r in smp if r.get('cls') == 'certified' and r['S'][:23] not in known and r['S'][:23] not in newS}
        to_known = sum(1 for r in smp if r.get('cls') == 'certified' and r['S'][:23] in known)
        still = sum(1 for r in smp if r.get('cls') != 'certified')
        frac = len(to_new) / len(smp) if smp else 0.0
        est = len(newS) + frac * A['nm']
        A.update(cert=len(cert), new=len(newS), far=far)
        print(f"{arm:9s} {A['ph']:6.2f} {len(A['keys']):5d} {len(cert):4d} {len(newS):4d} {far:11d} | {A['nm']:7d} {len(smp):7d} "
              f"{len(to_new):5d} {to_known:7d} {still:5d} | {len(newS) / A['ph']:5.1f} .. {est / A['ph']:5.1f}")
        if newS:
            best = min(newS, key=lambda S: float(S[:20]))
            print(f"          best new {best[:16]}; new minima d from record: "
                  f"{sorted(round(dist(p), 3) for p in newS.values())[:12]}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('seeds'); s.add_argument('--donors', type=int, default=40); s.add_argument('--per', type=int, default=3)
    s.add_argument('--control', type=int, default=20); s.add_argument('--tries', type=int, default=4)
    s.add_argument('--procs', type=int, default=10); s.add_argument('--seed', type=int, default=0)
    e = sp.add_parser('explore'); e.add_argument('--runs', type=int, default=8); e.add_argument('--procs-per', type=int, default=2)
    e.add_argument('--minutes', type=float, default=20); e.add_argument('--smax-seed', type=float, default=11.02)
    e.add_argument('--baseline', type=int, default=2); e.add_argument('--control-runs', type=int, default=1)
    r = sp.add_parser('report'); r.add_argument('--procs', type=int, default=12); r.add_argument('--sample', type=int, default=15)
    r.add_argument('--minutes', type=float, default=30); r.add_argument('--procs-per', type=int, default=2); r.add_argument('--no-solve', action='store_true')
    a = ap.parse_args()
    {'seeds': seeds, 'explore': explore, 'report': report}[a.cmd](a)
