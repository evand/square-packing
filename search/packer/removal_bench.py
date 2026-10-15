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


def run_event(ev, a, out):
    from pk.store import Store
    d = f'{out}/n{ev["n"]}_{ev["T"][:10]}_{ev["pm"]}m{ev["k"]}'
    os.makedirs(d, exist_ok=True)
    pr = parent_row(ev['pm'], ev['parent_side'])
    if not pr:
        json.dump(dict(ev, status='no-parent'), open(f'{d}/result.json', 'w')); return
    st = Store(); p = st.get(pr[0]); st.close()
    s, sq = float(p.s), [tuple(map(float, q)) for q in p.sq]
    rng = random.Random(ev['n'])
    if ev['k'] == 1:
        idxs = [(i,) for i in range(len(sq))]
    else:
        idxs = list({tuple(sorted(rng.sample(range(len(sq)), ev['k']))) for _ in range(a.pairs)})
    jobs = [(s, sq, set(ix), d) for ix in idxs]
    res = []
    with ProcessPoolExecutor(a.procs) as ex:
        for idx, s2, sq2 in ex.map(quench_one, jobs, chunksize=4):
            if s2 is not None and s2 < s - 1e-7:
                res.append((s2, idx, sq2))
    res.sort(key=lambda r: r[0])
    keep, seen = [], []
    for s2, idx, sq2 in res:
        if any(abs(s2 - x) < 2e-9 for x in seen):
            continue
        seen.append(s2); f = f'{d}/start_{len(keep):02d}.txt'; mcmin.write_deg(f, s2, sq2); keep.append(f)
        if len(keep) >= a.keep:
            break
    rem_best = res[0][0] if res else None
    best = rem_best
    if keep and a.minutes > 0:
        cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(ev['n']), '--adapt', '--reward', 'value', '--elite-share', '0.3',
               '--polish-gate', '2e-3', '--frontier-s', repr(ev['ante']), '--starts', *keep, '--minutes', str(a.minutes),
               '--procs', str(a.procs), '--seed', str(ev['n']), '--out', f'{d}/explore']
        subprocess.run(cmd, stdout=open(f'{d}/explore.log', 'w'), stderr=subprocess.STDOUT, cwd=HERE)
        try:
            E = [json.loads(l) for l in open(f'{d}/explore/archive.jsonl')]
            b = min((e['s'] for e in E if not e.get('unpolished')), default=None)
            best = min(x for x in (best, b) if x is not None)
        except FileNotFoundError:
            pass
    tier = 'REPLICATED' if best is not None and best <= ev['new'] + 1e-9 else \
        'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
    json.dump(dict(ev, parent_id=pr[0], removals=len(idxs), removal_best=rem_best, best=best, tier=tier,
                   need=(best - ev['ante']) if best is not None else None), open(f'{d}/result.json', 'w'), indent=1)
    print(f"n={ev['n']} {ev['T'][:16]} parent {ev['pm']}-{ev['k']}: removal best {rem_best}, final {best}, ante {ev['ante']:.10f}, "
          f"posted {ev['new']:.10f} -> {tier}", flush=True)


def main():
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    sp.add_parser('list')
    r = sp.add_parser('run'); r.add_argument('--only'); r.add_argument('--minutes', type=float, default=20)
    r.add_argument('--procs', type=int, default=2); r.add_argument('--par', type=int, default=8)
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
        man = manifest.start(a.out, 'removal_bench', args=vars(a), procs=a.par * a.procs)
        evs = events()
        if a.only:
            keep = {int(x) for x in a.only.split(',')}; evs = [e for e in evs if e['n'] in keep]
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(a.par) as tx:
            list(tx.map(lambda ev: run_event(ev, a, a.out), evs))
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
