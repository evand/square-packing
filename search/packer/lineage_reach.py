#!/usr/bin/env python3
"""Lineage reach test (10-10): does chaining (explorer lineages) reach basins that single moves cannot?

Single-step reach is ~0.05 matching distance for every move (c110r, battery 10-09: 0 withheld hits from parents at
0.05-0.1 in 1,144 proposals).  Here explorers start only from certified s(110) basins at matching distance [dlo, dhi)
from the record; everything closer is unknown to the run (novelty reference = known basins >= dlo from the record, so
closer basins count as novel, as they would in a real search).  Measured per run: time to the first sub-11 archive
entry within 0.05 / 0.03 / 0.02 of the record, the record itself, and the distance frontier over time.  Variants run
concurrently with the same starts and paired seeds.

  lineage_reach.py freeze  [--dlo 0.05 --dhi 0.10 --per-band 6]     -> runs/lineage/starts/, novel_ref.json
  lineage_reach.py run     --name L1 --minutes 60 --procs 4 --reps 2   (variants: new = cross2 + ashallow, old = without)
  lineage_reach.py report  --name L1 [--certify]
"""
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = f'{HERE}/runs/lineage'
REC = 40600                                           # store id of the s(110) record (battery suite_110 'record')
VARIANTS = {'new': [], 'old': ['--cross2', '0', '--ashallow', '0']}
RECIPE = ['--adapt', '--elite-share', '0.3', '--reward', 'value']    # hunt recipe (chain.py) + value reward
BANDS = (0.05, 0.03, 0.02)


def freeze(a):
    from pk.store import Store
    S = json.load(open(f'{HERE}/runs/battery/suite_110.json'))
    st = Store()
    D = sorted((d, int(i)) for i, d in S['reach_d'].items())
    edges = [a.dlo, (a.dlo + a.dhi) / 2, a.dhi]
    os.makedirs(f'{ROOT}/starts', exist_ok=True)
    picked = []
    for lo, hi in zip(edges, edges[1:]):
        band = [(d, i) for d, i in D if lo <= d < hi]
        step = max(1, len(band) // a.per_band)
        picked += band[::step][:a.per_band]
    for d, i in picked:
        st.get(i).write(f'{ROOT}/starts/b{i}.txt')
    ref = sorted(st.get(i).s for d, i in D if d >= a.dlo)
    json.dump(ref, open(f'{ROOT}/novel_ref.json', 'w'))
    json.dump(dict(dlo=a.dlo, dhi=a.dhi, starts=[dict(id=i, d=d) for d, i in picked], nref=len(ref),
                   frozen=time.strftime('%F %T')), open(f'{ROOT}/meta.json', 'w'), indent=1)
    print(f'{len(picked)} starts, d = {picked[0][0]:.3f} .. {picked[-1][0]:.3f}; novelty ref {len(ref)} sides')


GROUPS = {                                    # (gap = s - s_record range, d range); 10-10 after L1 (Evan: near / far in side)
    'near': ((0, 5e-4), (0.05, 0.07)),       # the d ~ 0.062 family, almost as good as the record
    'far1': ((1.5e-3, 1), (0.05, 0.07)),
    'far2': ((1.5e-3, 1), (0.07, 0.10)),
}


def freeze_groups(a):
    from pk.store import Store
    S = json.load(open(f'{HERE}/runs/battery/suite_110.json'))
    st = Store(); rec_s = st.get(REC).s
    D = sorted((d, int(i), st.get(int(i)).s - rec_s) for i, d in S['reach_d'].items())
    meta = {}
    for g, ((glo, ghi), (dlo, dhi)) in GROUPS.items():
        pool = [(d, i, gap) for d, i, gap in D if glo <= gap < ghi and dlo <= d < dhi]
        step = max(1, len(pool) // a.per_group)
        pick = pool[::step][:a.per_group]
        os.makedirs(f'{ROOT}/groups/{g}', exist_ok=True)
        for d, i, gap in pick:
            st.get(i).write(f'{ROOT}/groups/{g}/b{i}.txt')
        meta[g] = [dict(id=i, d=round(d, 4), gap=gap) for d, i, gap in pick]
        print(f'{g}: pool {len(pool)}, picked {len(pick)}: d {pick[0][0]:.3f}..{pick[-1][0]:.3f}, '
              f'gap {min(x[2] for x in pick):.1e}..{max(x[2] for x in pick):.1e}')
    json.dump(dict(groups=meta, frozen=time.strftime('%F %T')), open(f'{ROOT}/groups/meta.json', 'w'), indent=1)


def run(a):
    sets = {g: sorted(f'{ROOT}/groups/{g}/{f}' for f in os.listdir(f'{ROOT}/groups/{g}') if f.endswith('.txt'))
            for g in a.groups} if a.groups else {'': sorted(f'{ROOT}/starts/{f}' for f in os.listdir(f'{ROOT}/starts'))}
    P = []
    for rep in range(a.reps):
      for g, starts in sets.items():
        for v, extra in VARIANTS.items():
            out = f'{ROOT}/{a.name}/{g + "_" if g else ""}{v}{a.seed_base + rep}'
            os.makedirs(out, exist_ok=True)
            cmd = [sys.executable, f'{HERE}/explore.py', '--n', '110', '--starts', *starts, '--minutes', str(a.minutes),
                   '--procs', str(a.procs), '--seed', str(100 + a.seed_base + rep), '--novel-ref', f'{ROOT}/novel_ref.json',
                   '--out', out] + RECIPE + extra
            P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
    json.dump(dict(minutes=a.minutes, procs=a.procs, reps=a.reps, seed_base=a.seed_base, groups=a.groups, variants=VARIANTS, recipe=RECIPE,
                   started=time.strftime('%F %T')), open(f'{ROOT}/{a.name}/meta{a.seed_base}.json', 'w'), indent=1)
    for p in P:
        p.wait()


def _dist(path):
    from pk.packing import read
    from pk.store import Store
    from pk.crossing import match_distance
    global _REC
    if '_REC' not in globals():
        _REC = Store().get(REC)
    return path, match_distance(read(path), _REC)


def report(a):
    base = f'{ROOT}/{a.name}'
    cache_p = f'{base}/dist_cache.json'
    cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}
    runs = sorted(d for d in os.listdir(base) if os.path.isdir(f'{base}/{d}'))
    E = {}
    for r in runs:
        p = f'{base}/{r}/archive.jsonl'
        E[r] = [json.loads(l) for l in open(p)] if os.path.exists(p) else []
    todo = [e['path'] for r in runs for e in E[r] if e['s'] < 11 and e['path'] not in cache]
    if todo:
        with ProcessPoolExecutor(a.procs) as ex:
            for pth, d in ex.map(_dist, todo, chunksize=8):
                cache[pth] = d
        json.dump(cache, open(cache_p, 'w'))
    rec_s = 10.996783396634157
    print(f'{"run":6s} {"sub11":>6s} {"min d":>6s}  ' + '  '.join(f't(<{b})' for b in BANDS) + '  t(rec)   n<0.05 n<0.02')
    hits = []
    for r in runs:
        sub = [e for e in E[r] if e['s'] < 11 and e['parent'] >= 0]
        ds = [(e['t'], cache[e['path']], e) for e in sub if e['path'] in cache]
        if not ds:
            print(f'{r:6s} {len(sub):6d}   (no data)'); continue
        first = []
        for b in BANDS:
            ts = [t for t, d, _ in ds if d < b]
            first.append(f'{min(ts) / 60:6.1f}m' if ts else '     - ')
        tr = [t for t, d, e in ds if abs(e['s'] - rec_s) < 5e-9]
        print(f'{r:6s} {len(sub):6d} {min(d for _, d, _ in ds):6.3f}  ' + '  '.join(first) +
              f'  {(min(tr) / 60 if tr else float("nan")):6.1f}m  {sum(d < 0.05 for _, d, _ in ds):6d} {sum(d < 0.02 for _, d, _ in ds):6d}')
        hits += [(r, e, d) for _, d, e in ds if d < 0.05]
        if a.trace:                                    # distance frontier over time (10-min bins)
            tmax = max(t for t, _, _ in ds)
            row = []
            for m in range(10, int(tmax / 60) + 11, 10):
                dd = [d for t, d, _ in ds if t < 60 * m]
                row.append(f'{m}m:{min(dd):.3f}' if dd else f'{m}m:-')
            print('       ' + ' '.join(row))
    # lineage depth of the hits: number of explorer steps from a start
    for r, e, d in sorted(hits, key=lambda h: h[2])[:a.show]:
        idx = {x['i']: x for x in E[r]}
        depth, x, kinds = 0, e, []
        bypath = {y['path']: y for y in E[r]}
        root = e.get('root')
        while x['parent'] >= 0:
            kd = x['kind']
            if x.get('mate') in bypath:            # crossover mate: same root (own lineage) or another start's lineage
                m = bypath[x['mate']]
                kd += f'[{"own" if m.get("root") == root else "borrow:" + idx[m["root"]]["kind"].split(":")[-1]}' \
                      f'{"," + format(cache[m["path"]], ".3f") if m["path"] in cache else ""}]'
            kinds.append(kd); x = idx[x['parent']]; depth += 1
        print(f'  {r} s={e["s"]:.10f} d={d:.4f} t={e["t"] / 60:.1f}m depth {depth} from {x["kind"]}: {" <- ".join(kinds[:a.chain])}')
    if a.certify and hits:
        import census_exact, tempfile
        from pk.packing import read
        for r, e, d in sorted(hits, key=lambda h: h[2])[:a.show]:
            with tempfile.TemporaryDirectory() as tmp:
                read(e['path']).write(f'{tmp}/in.txt')
                c = census_exact.solve((f'{tmp}/in.txt', tmp))
            print(f'  cert {r} s={e["s"]:.10f} d={d:.4f}: {c.get("cls")} {c.get("status", "")}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    f = sp.add_parser('freeze'); f.add_argument('--dlo', type=float, default=0.05); f.add_argument('--dhi', type=float, default=0.10)
    f.add_argument('--per-band', type=int, default=6)
    fg = sp.add_parser('freeze-groups'); fg.add_argument('--per-group', type=int, default=4)
    r = sp.add_parser('run'); r.add_argument('--name', required=True); r.add_argument('--minutes', type=float, default=60)
    r.add_argument('--procs', type=int, default=4); r.add_argument('--reps', type=int, default=2)
    r.add_argument('--groups', nargs='*', help='start groups (runs/lineage/groups/<g>), one explorer per group x variant')
    r.add_argument('--seed-base', type=int, default=0, help='replicate offset (a later batch into the same --name)')
    p = sp.add_parser('report'); p.add_argument('--name', required=True); p.add_argument('--procs', type=int, default=8)
    p.add_argument('--certify', action='store_true'); p.add_argument('--trace', action='store_true')
    p.add_argument('--show', type=int, default=10); p.add_argument('--chain', type=int, default=16)
    a = ap.parse_args()
    {'freeze': freeze, 'freeze-groups': freeze_groups, 'run': run, 'report': report}[a.cmd](a)
