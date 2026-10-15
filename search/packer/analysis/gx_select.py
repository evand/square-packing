#!/usr/bin/env python3
"""Grid-crossing seeds at one n (10-15, Evan: "grid-crossing descents are of interest in general"): corner signatures
(analysis/cornersig.py) of every distinct polished archive minimum in the given runs, split below / above the grid;
above-grid corner cells that never occur below the grid ("foreign" cells), best member per cell -> seeds.

  python3 analysis/gx_select.py --n 110 --runs 'runs/c110_h*' 'runs/h110_pm' 'runs/bench110' --out runs/gx110 [--within 2e-3]
"""
import argparse, glob, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def _sig(args):
    s, path = args
    import cornersig as CS
    try:
        with open(path) as fh:
            n, s2 = fh.readline().split()[:2]
            sq = np.array([list(map(float, x.split()[:3])) for x in fh if x.strip()])
        g = CS.signature(float(s2), sq)
    except Exception:
        return None
    so = g['sorted']
    return dict(s=s, path=path, cell=[so[0], so[1], so[2] + so[3]], sorted=so, T=g['T'], U=g['U'],
                slid=g['slid_rows'] + g['slid_cols'], cutoff=g['cutoff'], free=g['free'], L=g['L'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--within', type=float, default=2e-3)
    ap.add_argument('--procs', type=int, default=8)
    a = ap.parse_args()
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    E, seen = [], set()
    for g in a.runs:
        for f in glob.glob(f'{g}/**/archive.jsonl', recursive=True):
            for l in open(f):
                e = json.loads(l)
                if e.get('unpolished') or e['s'] >= k + 0.05:
                    continue
                key = round(e['s'], 9)
                if key in seen or abs(e['s'] - k) < 1e-9:
                    continue
                seen.add(key); E.append((e['s'], e['path']))
    with ProcessPoolExecutor(a.procs) as ex:
        R = [r for r in ex.map(_sig, E, chunksize=64) if r]
    below = [r for r in R if r['s'] < k]; above = [r for r in R if r['s'] > k]
    bcells = {tuple(r['cell']) for r in below}
    best = {}
    for r in above:
        c = tuple(r['cell']) + (r['cutoff'] >= 3,)
        if tuple(r['cell']) in bcells and r['cutoff'] < 3:
            continue                                       # a cell the sub-grid family already occupies
        if c not in best or r['s'] < best[c]['s']:
            best[c] = r
    seeds = sorted([r for r in best.values() if r['s'] - k < a.within], key=lambda r: r['s'])
    os.makedirs(a.out, exist_ok=True)
    with open(f'{a.out}/above.jsonl', 'w') as fa:
        for r in above:
            fa.write(json.dumps(r) + '\n')
    json.dump(dict(n=a.n, below_cells=sorted(map(list, bcells)), n_below=len(below), n_above=len(above),
                   foreign_cells=len(best), seeds=seeds), open(f'{a.out}/seeds.json', 'w'), indent=0)
    print(f'n={a.n}: {len(below)} distinct below-grid minima ({len(bcells)} corner cells), {len(above)} above; '
          f'{len(best)} foreign above-grid cells; {len(seeds)} with best member within {a.within:g} of the grid')


if __name__ == '__main__':
    main()
