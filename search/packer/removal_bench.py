#!/usr/bin/env python3
"""Removal baseline, level 1 (method fidelity; Evan 10-15: get to baseline on the major tools in order).

For every family-F record event at n <= 324 with a named parent "X s(m) - k" (bench/events.json): take the named parent
packing (store row with the parent side; the parent is the frontier of m at T or a packing posted in the same update),
do the removals (k = 1: every square; k = 2: --pairs random pairs), quench each with fq (chain's seed quench), keep the
--keep best distinct, then run the explorer from them (current config: gate 2e-3, value reward, elite 0.3, frontier_s =
frontier_at(n, T)) for --minutes on --procs.  Success tiers: REPLICATED (<= the posted side + 1e-9), BEAT-ANTE
(< frontier_at(n, T) - 1e-9), no.

  python3 removal_bench.py list
  python3 removal_bench.py run [--only 179,208] [--minutes 20 --procs 2 --par 8] --out runs/removal_bench/b1
  python3 removal_bench.py score runs/removal_bench/b1
"""
import argparse, glob, json, math, os, random, re, sqlite3, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import mcmin                                                           # noqa: E402
from pk import eventbench as EB                                        # noqa: E402

DB = os.path.join(HERE, 'runs', 'store.sqlite')


def events():
    E = json.load(open(os.path.join(HERE, 'bench', 'events.json')))
    E = E['events'] if isinstance(E, dict) else E
    out = []
    for e in E:
        if e['n'] > 324 or e['family'] != 'F' or e['gain'] < 1e-8:
            continue
        m = re.search(r's\((\d+)\)\s*-\s*(\d+)', str(e['parents']))
        if not m:
            continue
        pm, k = int(m[1]), int(m[2])
        T = e['date']
        same = [x['new'] for x in E if x['n'] == pm and x['date'] == T]        # parent posted in the same update
        ps = min([EB.frontier_at(pm, T)] + same)
        out.append(dict(n=e['n'], T=T, author=e['author'], issue=e['issue'], new=e['new'], ante=EB.frontier_at(e['n'], T),
                        pm=pm, k=k, parent_side=ps, parents=e['parents']))
    return out


def parent_row(pm, ps):
    db = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    r = db.execute("SELECT id, s FROM packing WHERE n=? AND abs(s-?)<2e-9 AND status!='infeasible' ORDER BY abs(s-?) LIMIT 1",
                   (pm, ps, ps)).fetchone()
    return r


def quench_one(j):
    s, sq, idx, tmpd = j
    import hop
    from chain import SEED_Q
    sq2 = [q for i, q in enumerate(sq) if i not in idx]
    with tempfile.TemporaryDirectory(dir=tmpd) as tmp:
        r = hop.quench(s, sq2, tmp, opt=SEED_Q)
    if not r:
        return idx, None, None
    return idx, r[0], r[1]


def edir(ev, out):
    return f'{out}/n{ev["n"]}_{ev["T"][:10]}_{ev["pm"]}m{ev["k"]}'


def prepare(ev, a, out):
    from pk.store import Store
    d = edir(ev, out); os.makedirs(d, exist_ok=True)
    pr = parent_row(ev['pm'], ev['parent_side'])
    if not pr:
        return []
    st = Store(); p = st.get(pr[0]); st.close()
    s, sq = float(p.s), [tuple(map(float, q)) for q in p.sq]
    ev['parent_id'] = pr[0]
    rng = random.Random(ev['n'])
    if ev['k'] == 1:
        idxs = [(i,) for i in range(len(sq))]
    else:
        idxs = list({tuple(sorted(rng.sample(range(len(sq)), ev['k']))) for _ in range(a.pairs)})
    ev['removals'] = len(idxs); ev['parent_s'] = s
    return [(s, sq, set(ix), d) for ix in idxs]


def pick_starts(ev, res, a, out):
    d = edir(ev, out)
    s0 = res[0][0] if res else None
    ok = sorted((r for r in res if r[1] is not None), key=lambda r: r[1])
    keep, seen = [], []
    for idx, s2, sq2 in ok:
        if any(abs(s2 - x) < 2e-9 for x in seen):
            continue
        seen.append(s2); f = f'{d}/start_{len(keep):02d}.txt'; mcmin.write_deg(f, s2, sq2); keep.append(f)
        if len(keep) >= a.keep:
            break
    ev['removal_best'] = ok[0][1] if ok else None
    return keep


def finish(ev, out):
    d = edir(ev, out); best = ev.get('removal_best')
    for ad in glob.glob(f'{d}/explore*/archive.jsonl'):
        for l in open(ad):
            e = json.loads(l)
            if not e.get('unpolished') and (best is None or e['s'] < best):
                best = e['s']
    tier = 'REPLICATED' if best is not None and best <= ev['new'] + 1e-9 else \
        'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
    json.dump(dict(ev, best=best, tier=tier, need=(best - ev['ante']) if best is not None else None),
              open(f'{d}/result.json', 'w'), indent=1)


def run_pooled(evs, a):
    """phase 1: every removal quench of every event in one --cores queue; phase 2: every event's explorer at once
    (1 proc each), spare cores to extra replicates (Evan 10-15: no waves, no per-event pools)."""
    from jobpool import run_jobs
    jobs, owner = [], []
    for i, ev in enumerate(evs):
        J = prepare(ev, a, a.out); jobs += J; owner += [i] * len(J)
    print(f'{len(evs)} events, {len(jobs)} removal quenches', flush=True)
    res = run_jobs(quench_one, jobs, procs=a.cores, timeout=900)
    per = {i: [] for i in range(len(evs))}
    for i, r in zip(owner, res):
        if isinstance(r, tuple) and r[1] is not None and r[1] < evs[i]['parent_s'] - 1e-7:
            per[i].append(r)
    runs = []
    for i, ev in enumerate(evs):
        keep = pick_starts(ev, per[i], a, a.out)
        if keep and a.minutes > 0:
            runs.append((ev, keep))
    slots = [(ev, keep, 0) for ev, keep in runs]
    r = 1
    while len(slots) < a.cores and runs:
        slots += [(ev, keep, r) for ev, keep in runs][:a.cores - len(slots)]; r += 1
    P = []
    for ev, keep, rep in slots:
        d = edir(ev, a.out)
        cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(ev['n']), '--adapt', '--reward', 'value', '--elite-share', '0.3',
               '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(ev['ante']), '--starts', *keep,
               '--minutes', str(a.minutes), '--procs', '1', '--seed', str(ev['n'] * 10 + rep), '--out', f'{d}/explore{rep}']
        P.append(subprocess.Popen(cmd, stdout=open(f'{d}/explore{rep}.log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
    for p in P:
        p.wait()
    for ev in evs:
        finish(ev, a.out)


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    sp.add_parser('list')
    r = sp.add_parser('run'); r.add_argument('--only'); r.add_argument('--minutes', type=float, default=20)
    r.add_argument('--cores', type=int, default=16)
    r.add_argument('--keep', type=int, default=8); r.add_argument('--pairs', type=int, default=400); r.add_argument('--out', required=True)
    s = sp.add_parser('score'); s.add_argument('out')
    a = ap.parse_args()
    if a.cmd == 'list':
        for ev in events():
            pr = parent_row(ev['pm'], ev['parent_side'])
            print(f"{ev['n']:4d} {ev['T'][:16]} {ev['author'][:12]:12s} s({ev['pm']})-{ev['k']} parent {ev['parent_side']:.10f} "
                  f"store {pr[0] if pr else 'MISSING'}  ante {ev['ante']:.10f} posted {ev['new']:.10f}")
    elif a.cmd == 'run':
        from pk import manifest
        os.makedirs(a.out, exist_ok=True)
        man = manifest.start(a.out, 'removal_bench', args=vars(a), procs=a.cores)
        evs = events()
        if a.only:
            keep = {int(x) for x in a.only.split(',')}; evs = [e for e in evs if e['n'] in keep]
        run_pooled(evs, a)
        man.finish()
    else:
        R = [json.load(open(f)) for f in sorted(glob.glob(f'{a.out}/*/result.json'))]
        from collections import Counter
        for r in R:
            print(f"{r['n']:4d} s({r['pm']})-{r['k']}  removal {r.get('removal_best') or float('nan'):.10f}  final "
                  f"{r.get('best') or float('nan'):.10f}  need {r.get('need') or float('nan'):+.2e}  {r.get('tier', r.get('status'))}")
        print(Counter(r.get('tier', r.get('status')) for r in R))


if __name__ == '__main__':
    main()
