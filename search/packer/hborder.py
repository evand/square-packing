#!/usr/bin/env python3
"""H183's recipe at other targets (10-19; Evan: "let's try it" on 147, 157, 274).  h183.py, parameterized; h183.py is
kept as run for reproducibility.

  target n with frontier side T (the grid k at 157 / 274; the 45-degree tie 7 + 4 sqrt 2 at 147).
  seeds:   arm 'record' = one square inserted into best(n - 1), many seeds; 'd1' = one square into distinct
           (n - 1)-packings (best half + stratified by side band x tilted count, s < T - gap1); 'd2' = two squares into
           distinct (n - 2)-packings (s < T - gap2).  chain.seed_job (dilate 0.2-1.2 %, pose hole, fq quench).
           Donors with store status 'infeasible' are skipped (h183.py did not filter; e.g. 308 id 212401 kingbird).
  explore: explorer runs per arm from the best runs x starts-per-run distinct line-free (layout.full_lines at k = ceil
           sqrt n), non-tie seeds (cap T + --seed-gap), round-robin; equal minutes per run.  Explorer --smax = max(T + 0.06,
           worst start + 0.01) (h183: fixed cap 14.04 / smax 14.06; 147's seeds quench ~0.12 above its tie).
           Arm 'lcut' joins when runs/hborder/n<n>/lcut/seeds.jsonl exists (lcut.py's line-free 'lcut' seeds).
  report:  seed side distribution per arm; explorer best per run; every run best below T.

  hborder.py seeds   --n 157 --target 13 [--donors 120 --per 4 --tries 4 --procs 15]
  hborder.py explore --n 157 --target 13 [--runs-per-arm 4 --procs-per 1 --minutes 40]
  hborder.py report  --n 157 --target 13
"""
import argparse, collections, glob, json, math, os, random, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ARMS = ('record', 'd1', 'd2')


def root(a):
    return f'{HERE}/runs/hborder/n{a.n}{a.tag}'


def _seed(j):
    import chain
    chain.TRIES = j[-1]
    return chain.seed_job(j[:-1])


def pick(st, n, m, smax, rng):
    """the best m // 2 distinct sides + m // 2 stratified by (side band, tilted count) over the rest below smax"""
    rows = st.db.execute("SELECT id, s, ntilted FROM packing WHERE n=? AND s < ? AND status != 'infeasible' ORDER BY s",
                         (n, smax)).fetchall()
    seen, D = set(), []
    for pid, s, t in rows:
        k = round(s, 7)
        if k in seen:
            continue
        seen.add(k); D.append((pid, s, t or 0))
    out = D[:m // 2]
    rest = D[m // 2:]
    strata = collections.defaultdict(list)
    if rest:
        lo, hi = rest[0][1], smax
        for d in rest:
            strata[(min(4, int(5 * (d[1] - lo) / (hi - lo + 1e-12))), d[2] // 6)].append(d)
    keys = sorted(strata)
    while len(out) < m and any(strata.values()):
        for k in keys:
            if strata[k] and len(out) < m:
                out.append(strata[k].pop(rng.randrange(len(strata[k]))))
    return out


def seeds(a):
    from pk.store import Store
    from pk import manifest
    from jobpool import run_jobs
    import scratch
    R0 = root(a)
    scratch.scope(f'hborder{a.n}')
    os.makedirs(R0, exist_ok=True)
    man = manifest.start(R0, f'hborder-seeds-{a.n}', seed=a.seed, args=vars(a), procs=a.procs, fname='manifest_seeds.json')
    st = Store(); rng = random.Random(a.seed)
    for d in ('donors', 'seeds'):
        os.makedirs(f'{R0}/{d}', exist_ok=True)
    if a.donor_file:
        return seeds_ab(a, st, rng, R0, man)
    rec = st.db.execute("SELECT id FROM packing WHERE n=? AND status != 'infeasible' ORDER BY s LIMIT 1",
                        (a.n - 1,)).fetchone()[0]
    d1 = [x for x in pick(st, a.n - 1, a.donors + 1, a.target - a.gap1, rng) if x[0] != rec][:a.donors]
    d2 = pick(st, a.n - 2, a.donors2, a.target - a.gap2, rng)
    print(f'n {a.n} target {a.target}: record {rec}; d1 {len(d1)} donors ({d1[0][1]:.6f}..{d1[-1][1]:.6f}); '
          f'd2 {len(d2)} donors', flush=True)
    J, meta = [], []
    st.get(rec).write(f'{R0}/donors/d{rec}.txt')
    for q in range(a.per_record):
        J.append((f'{R0}/donors/d{rec}.txt', -1, a.n, f'{R0}/seeds/r{rec}_{q}.txt', rng.randrange(10 ** 9), a.tries))
        meta.append(dict(arm='record', src=rec, q=q))
    for arm, D, k in (('d1', d1, -1), ('d2', d2, -2)):
        for pid, s, t in D:
            st.get(pid).write(f'{R0}/donors/d{pid}.txt')
            for q in range(a.per):
                J.append((f'{R0}/donors/d{pid}.txt', k, a.n, f'{R0}/seeds/{arm}_{pid}_{q}.txt', rng.randrange(10 ** 9), a.tries))
                meta.append(dict(arm=arm, src=pid, src_s=s, src_tilt=t, q=q))
    t0 = time.time()
    R = run_jobs(_seed, J, procs=a.procs, timeout=1800)
    from layout import full_lines
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    with open(f'{R0}/seeds.jsonl', 'w') as fo:
        for j, m, r in zip(J, meta, R):
            s = r if isinstance(r, float) else None
            m.update(path=j[3], s=s, lines=(full_lines(j[3], k) if s is not None else None))
            fo.write(json.dumps(m) + '\n')
    summary(a, meta)
    print(f'seeds done [{time.time() - t0:.0f}s]')
    man.finish()


def seeds_ab(a, st, rng, R0, man):
    """HB2 (10-19): arm 'new' = one-square insertions into --donor-file's packings (d146.py report), arm 'ctl' = as many
    fresh d1-rule donors from the store, excluding the donors of --exclude (an earlier seeds.jsonl); same per / tries."""
    from jobpool import run_jobs
    from layout import full_lines
    new = json.load(open(a.donor_file))['donors'][:a.donors]
    used = {m['src'] for m in map(json.loads, open(a.exclude))} if a.exclude else set()
    ctl = [x for x in pick(st, a.n - 1, 4 * a.donors, a.target - a.gap1, rng) if x[0] not in used]
    ctl = [] if a.no_ctl else rng.sample(ctl, min(len(new), len(ctl)))
    print(f'n {a.n}: new {len(new)} donors, ctl {len(ctl)} donors', flush=True)
    J, meta = [], []
    for i, d in enumerate(new):
        for q in range(a.per):
            J.append((d['path'], -1, a.n, f'{R0}/seeds/new_{i}_{q}.txt', rng.randrange(10 ** 9), a.tries))
            meta.append(dict(arm='new', src=d['path'], src_s=d['s'], src_tilt=d['tilt'], src_nnd=d['nnd'], q=q))
    for pid, s, t in ctl:
        st.get(pid).write(f'{R0}/donors/d{pid}.txt')
        for q in range(a.per):
            J.append((f'{R0}/donors/d{pid}.txt', -1, a.n, f'{R0}/seeds/ctl_{pid}_{q}.txt', rng.randrange(10 ** 9), a.tries))
            meta.append(dict(arm='ctl', src=pid, src_s=s, src_tilt=t, q=q))
    R = run_jobs(_seed, J, procs=a.procs, timeout=1800)
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    with open(f'{R0}/seeds.jsonl', 'w') as fo:
        for j, m, r in zip(J, meta, R):
            sv = r if isinstance(r, float) else None
            m.update(path=j[3], s=sv, lines=(full_lines(j[3], k) if sv is not None else None))
            fo.write(json.dumps(m) + '\n')
    summary(a, meta)
    man.finish()


def arms_of(meta):
    """arms in order of first appearance (record, d1, d2 / new, ctl / any composed seeds.jsonl, e.g. LC274b)"""
    return list(dict.fromkeys(m['arm'] for m in meta))


def summary(a, meta):
    T = a.target
    for arm in arms_of(meta):
        S = [m for m in meta if m['arm'] == arm]
        ok = [m for m in S if m['s'] is not None]
        if not ok:
            print(f'{arm}: none feasible'); continue
        free = sorted(m['s'] for m in ok if not m.get('lines'))
        print(f'{arm}: {len(S)} seeds, {len(ok)} feasible, line-free {len(free)}; line-free sides: < T {sum(s < T - 1e-9 for s in free)}, '
              f'tie {sum(abs(s - T) < 1e-9 for s in free)}, T..+.01 {sum(T + 1e-9 < s < T + .01 for s in free)}, '
              f'+.01..+.02 {sum(T + .01 <= s < T + .02 for s in free)}, +.02..+.05 {sum(T + .02 <= s < T + .05 for s in free)}, '
              f'above {sum(s >= T + .05 for s in free)}; best {[round(s, 6) for s in free[:3]]}', flush=True)


def explore(a):
    from pk import manifest
    R0 = root(a)
    M = [json.loads(l) for l in open(f'{R0}/seeds.jsonl')]
    M = [m for m in M if m['s'] is not None and not m['lines']]      # line-blocked seeds can't go below the grid
    arms = arms_of(M)
    lc = f'{R0}/lcut/seeds.jsonl'                 # lcut.py --out runs/hborder/n<n>/lcut: arm 'lcut' (its 'random' arm unused)
    if os.path.exists(lc):
        arms.append('lcut')
        M += [dict(arm='lcut', s=r['s'], path=r['path'], src=r['src'], cut=r['cut']) for r in map(json.loads, open(lc))
              if r.get('arm') == 'lcut' and r.get('s') and not r['grid']]
    man = manifest.start(R0, f'hborder-explore-{a.n}', args=vars(a), procs=a.runs_per_arm * len(arms) * a.procs_per,
                         fname='manifest_explore.json')
    P, i = [], 0
    for arm in arms:
        use = sorted((m for m in M if m['arm'] == arm and m['s'] < a.target + a.seed_gap and abs(m['s'] - a.target) > 1e-9),
                     key=lambda m: m['s'])
        seen, U = set(), []
        for m in use:
            if round(m['s'], 8) not in seen:
                seen.add(round(m['s'], 8)); U.append(m)
        U = U[:a.runs_per_arm * a.starts_per_run]
        if not U:
            print(f'{arm}: no seeds below {a.target + a.seed_gap}'); continue
        for r in range(a.runs_per_arm):
            R = U[r::a.runs_per_arm]
            if not R:
                continue
            out = f'{R0}/x{i}_{arm}'; i += 1
            os.makedirs(out, exist_ok=True)
            json.dump(R, open(f'{out}/starts.json', 'w'))
            smax = max(a.target + 0.06, max(m['s'] for m in R) + 0.01)     # 147's seeds sit ~0.12 above the tie
            cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(a.n), '--starts', *[m['path'] for m in R],
                   *(['--cpu-h', str(a.cpu_h / max(1, len(arms) * a.runs_per_arm))] if a.cpu_h else ['--minutes', str(a.minutes)]), '--procs', str(a.procs_per), '--seed', str(a.seed + 700 + i), '--adapt', '--smax',
                   str(smax), '--frontier-s', str(a.target), '--final-polish', '8', '--out', out]
            P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log.txt', 'w'), stderr=subprocess.STDOUT))
    for p in P:
        p.wait()
    man.finish()
    report(a)


def report(a):
    R0 = root(a)
    M = [json.loads(l) for l in open(f'{R0}/seeds.jsonl')]
    summary(a, M)
    allbest = []
    for d in sorted(glob.glob(f'{R0}/x*')):
        best = None
        for f in glob.glob(f'{d}/**/*.txt', recursive=True):
            try:
                h = open(f).readline().split()
                if len(h) == 2 and h[0] == str(a.n):
                    s = float(h[1])
                    if best is None or s < best[0]:
                        best = (s, f)
            except Exception:
                pass
        print(f'{os.path.basename(d)}: best {best}')
        if best:
            allbest.append(best)
    sub = sorted(b for b in allbest if b[0] < a.target - 1e-9)
    print(f'run bests below {a.target}: {len(sub)}', sub[:5])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    for c in ('seeds', 'explore', 'report'):
        p = sp.add_parser(c)
        p.add_argument('--n', type=int, required=True); p.add_argument('--target', type=float, required=True)
        p.add_argument('--seed', type=int, default=0); p.add_argument('--tag', default='')
        if c == 'seeds':
            p.add_argument('--donors', type=int, default=120); p.add_argument('--donors2', type=int, default=30)
            p.add_argument('--per', type=int, default=4); p.add_argument('--per-record', type=int, default=120)
            p.add_argument('--tries', type=int, default=4); p.add_argument('--procs', type=int, default=15)
            p.add_argument('--gap1', type=float, default=0.005); p.add_argument('--gap2', type=float, default=0.03)
            p.add_argument('--donor-file', help="HB2: arms 'new' (this json's donors) vs 'ctl' (fresh store donors)")
            p.add_argument('--exclude', help='seeds.jsonl whose donors the ctl arm skips')
            p.add_argument('--no-ctl', action='store_true', help="HB2c: arm 'new' only (baseline = an earlier ctl)")
        if c == 'explore':
            p.add_argument('--runs-per-arm', type=int, default=4); p.add_argument('--procs-per', type=int, default=1)
            p.add_argument('--starts-per-run', type=int, default=6); p.add_argument('--minutes', type=float, default=40)
            p.add_argument('--cpu-h', type=float, default=None, help='total CPU budget over all explorer runs (replaces --minutes)')
            p.add_argument('--seed-gap', type=float, default=0.3, help='sanity cap; starts = the best distinct seeds')
    a = ap.parse_args()
    if not a.seed:
        a.seed = a.n
    {'seeds': seeds, 'explore': explore, 'report': report}[a.cmd](a)
