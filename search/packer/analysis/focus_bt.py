#!/usr/bin/env python3
"""FB1 (10-16): do fresh run bests deserve focused expansion?  Offline backtest over explorer proposal logs.

Each proposal (proposals.jsonl: completion time t, sec, cpu, par, new, s, unp) is classed by its parent at submission
time (t - sec): rank among the run's polished archive entries (1 = run best), start or descendant, and the parent's
expansions so far.  Outcome: the proposal created a polished archive entry below the run best at completion - 1e-9.

  python3 analysis/focus_bt.py [--glob 'runs/**/proposals.jsonl'] [--json OUT]
"""
import argparse, collections, glob, json, math, os, sys
import numpy as np


def rank_bin(r):
    return '1' if r == 1 else '2-5' if r <= 5 else '6-20' if r <= 20 else '>20'


def exp_bin(x):
    return '0-4' if x < 5 else '5-19' if x < 20 else '20+'


def one_run(d):
    try:
        A = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
        P = [json.loads(l) for l in open(f'{d}/proposals.jsonl')]
        man = json.load(open(f'{d}/manifest.json'))
    except Exception:
        return None
    if not A or not P:
        return None
    fs = (man.get('args') or {}).get('frontier_s')
    E = {e['i']: e for e in A}
    pol = [e for e in A if not e.get('unpolished')]
    if not pol:
        return None
    T = np.array([e.get('t', 0.0) for e in pol]); S = np.array([e['s'] for e in pol])
    start_s = min((e['s'] for e in A if e.get('parent', -1) < 0), default=None)
    if start_s is None:
        return None
    need = None if fs is None else start_s - fs
    # submission order for expansion counts
    Q = sorted(range(len(P)), key=lambda q: P[q]['t'] - (P[q].get('sec') or 0))
    nexp = collections.Counter(); rows = []
    for q in Q:
        p = P[q]; par = E.get(p.get('par'))
        if par is None:
            continue
        ts = p['t'] - (p.get('sec') or 0)
        if par.get('unpolished'):
            rk = None
        else:
            m = T <= ts + 1e-9
            rk = int((S[m] < par['s'] - 1e-12).sum()) + 1
        ex = nexp[par['i']]; nexp[par['i']] += 1
        # run best at completion, excluding the proposal's own entry
        m2 = T <= p['t'] + 1e-9
        if p.get('new') is not None:
            m2 &= np.array([e['i'] != p['new'] for e in pol])
        best = S[m2].min() if m2.any() else start_s
        hit = (p.get('new') is not None and not p.get('unp') and p.get('s') is not None and p['s'] < best - 1e-9)
        rows.append(dict(rank=rk, desc=par.get('parent', -1) >= 0, exp=ex, cpu=(p.get('cpu') or p.get('sec') or 0.0),
                         hit=bool(hit), gain=(best - p['s']) if hit else 0.0))
    return dict(dir=d, n=(man.get('args') or {}).get('n'), need=need, rows=rows)


def table(rows, key):
    C = collections.defaultdict(lambda: [0, 0.0, 0])
    for r in rows:
        k = key(r)
        if k is None:
            continue
        c = C[k]; c[0] += 1; c[1] += r['cpu']; c[2] += r['hit']
    return {k: dict(props=v[0], cpu_h=round(v[1] / 3600, 2), hits=v[2], per_cpu_h=round(v[2] / (v[1] / 3600), 3) if v[1] else None)
            for k, v in sorted(C.items(), key=lambda kv: str(kv[0]))}


def cls(r):
    if r['rank'] is None:
        return 'unpolished'
    return f"rank {rank_bin(r['rank'])} {'desc' if r['desc'] else 'start'}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--glob', default='runs/**/proposals.jsonl'); ap.add_argument('--json')
    a = ap.parse_args()
    dirs = sorted({os.path.dirname(f) for f in glob.glob(a.glob, recursive=True)})
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(12) as ex:
        R = [r for r in ex.map(one_run, dirs, chunksize=4) if r]
    out = dict(runs=len(R), dirs=len(dirs))
    print(f'{len(R)} runs with logs ({len(dirs)} dirs)')
    for split, f in (('all', lambda r: True), ('at record', lambda r: r['need'] is not None and r['need'] <= 1e-9),
                     ('weak start', lambda r: r['need'] is not None and r['need'] > 1e-9),
                     ('no frontier_s', lambda r: r['need'] is None)):
        rs = [r for r in R if f(r)]
        rows = [x for r in rs for x in r['rows']]
        print(f'\n== {split}: {len(rs)} runs, {len(rows)} proposals, {sum(x["hit"] for x in rows)} new-best hits')
        t1 = table(rows, cls)
        for k, v in t1.items():
            print(f'  {k:22s} {v}')
        # within runs that had at least one descendant run best (the case the heuristic acts on)
        rs2 = [r for r in rs if any(x['hit'] for x in r['rows'])]
        rows2 = [x for r in rs2 for x in r['rows']]
        t2 = table(rows2, cls)
        print(f'  -- within the {len(rs2)} runs with >= 1 new-best hit:')
        for k, v in t2.items():
            print(f'  {k:22s} {v}')
        t3 = table([x for x in rows2 if x['rank'] == 1], lambda x: f"rank 1 {'desc' if x['desc'] else 'start'} exp {exp_bin(x['exp'])}")
        print('  -- rank 1 parents by expansions so far (same runs):')
        for k, v in t3.items():
            print(f'  {k:22s} {v}')
        r1 = t2.get('rank 1 desc', {}); mids = [t2.get(f'rank {b} {s}') for b in ('2-5', '6-20') for s in ('desc', 'start')]
        h = sum(m['hits'] for m in mids if m); c = sum(m['cpu_h'] for m in mids if m)
        ratio = (r1.get('per_cpu_h') or 0) / (h / c) if c and h else None
        print(f'  ** rank-1 descendant per CPU-h {r1.get("per_cpu_h")} vs ranks 2-20 {round(h / c, 3) if c else None}: ratio {ratio}')
        out[split] = dict(runs=len(rs), table=t1, within=t2, rank1_exp=t3, ratio=ratio)
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
