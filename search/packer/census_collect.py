#!/usr/bin/env python3
"""Census of certified local minima at one n, from any set of explore.py runs (explorer, chain rounds, benchmark runs).
Every archive entry with side < --below (default: k = ceil(sqrt n)) goes through exactsolve (census_exact.solve), cached
per run in <run>/exact.json (same cache as explore_exact.py).  Distinct certified minima are keyed by the exact side
(census_exact key, 1e-20); output runs/census<n>.json: one row per minimum (S, roles, occurrences, first run / entry / t),
plus the union with an optional list of already-known certified sides (--known, 2e-9).

  census_collect.py --n 110 --runs 'runs/ex*' 'runs/ch*/r*_n110' --known runs/known110_all.json --procs 15
  census_collect.py --n 129 --runs 'runs/ch1/r*_n129' --dry       # count what would be solved
"""
import argparse, collections, glob, json, math, os
from jobpool import run_jobs
import census_exact

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--below', type=float); ap.add_argument('--known'); ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--out'); ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    below = a.below if a.below is not None else k
    runs = sorted({d for g in a.runs for d in glob.glob(g) if os.path.exists(f'{d}/archive.jsonl')})
    # only runs at this n: the archive's start entries are packings of n squares
    todo, per = [], {}
    for d in runs:
        E = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
        if not E:
            continue
        with open(E[0]['path']) as fh:
            if int(fh.readline().split()[0]) != a.n:
                continue
        C = json.load(open(f'{d}/exact.json')) if os.path.exists(f'{d}/exact.json') else {}
        sub = [e for e in E if e['s'] < below]
        per[d] = (E, sub, C)
        todo += [(d, e) for e in sub if str(e['i']) not in C]
    print(f'n = {a.n}: {len(per)} runs, {sum(len(v[1]) for v in per.values())} archive entries below {below:g}, '
          f'{len(todo)} to solve', flush=True)
    if a.dry:
        return
    for d in {d for d, _ in todo}:
        os.makedirs(f'{d}/exact', exist_ok=True)
    res = run_jobs(census_exact.solve, [(e['path'], f'{d}/exact') for d, e in todo], procs=a.procs, timeout=900,
                   on_timeout=lambda j: dict(cls='unresolved', status='killed'),
                   on_error=lambda j, x: dict(cls='unresolved', status=str(x)[:60]))
    for (d, e), r in zip(todo, res):
        per[d][2][str(e['i'])] = r or dict(cls='unresolved', status='none')
    for d, (E, sub, C) in per.items():
        json.dump(C, open(f'{d}/exact.json', 'w'), indent=0)
    rows = {}
    cls = collections.Counter()
    for d, (E, sub, C) in per.items():
        for e in sub:
            r = C[str(e['i'])]
            cls[r['cls']] += 1
            if r['cls'] != 'certified':
                continue
            row = rows.setdefault(r['key'], dict(S=r['key'], roles=e.get('roles'), count=0, first=dict(run=d, i=e['i'], t=e.get('t'), path=e['path'])))
            row['count'] += 1
    known = json.load(open(a.known)) if a.known else []
    new = [r for r in rows.values() if not any(abs(float(r['S']) - x) < 2e-9 for x in known)]
    union = sorted(set(round(float(r['S']), 9) for r in rows.values()) | set(round(x, 9) for x in known))
    out = a.out or f'{HERE}/runs/census{a.n}.json'
    json.dump(dict(n=a.n, below=below, runs=list(per), classes=dict(cls), minima=sorted(rows.values(), key=lambda r: float(r['S'])),
                   known=a.known, union_sides=union), open(out, 'w'), indent=0)
    print(f'  archive entries by class: {dict(cls)}; distinct certified minima {len(rows)}; not in --known: {len(new)}; '
          f'union {len(union)}  -> {out}')
    best = sorted(rows.values(), key=lambda r: float(r['S']))[:5]
    for r in best:
        print(f'    {float(r["S"]):.12f} roles {r["roles"]} x{r["count"]} first {r["first"]["run"]}#{r["first"]["i"]}')


if __name__ == '__main__':
    main()
