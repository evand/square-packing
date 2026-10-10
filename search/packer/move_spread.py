#!/usr/bin/env python3
"""Per-move landing spread (10-10): matching distance child -> parent for polished sub-k non-grid children in the idea
battery's trial files (bench + reach suites, n = 110), with the share of proposals lost (screen-grid / discard / grid)
and of children that are the parent's own basin (d < 0.005) or a different certified basin.

  move_spread.py [--per 250] [--procs 4]   -> runs/battery/move_spread.json + table
"""
import argparse, json, os, random, sys, collections
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = f'{HERE}/runs/battery'


def _d(job):
    pid, s, fin = job
    from pk.store import Store; from pk.packing import Packing; from pk.crossing import match_distance
    import numpy as np
    global _ST
    if '_ST' not in globals(): _ST = Store()
    return match_distance(Packing(s, np.array(fin)), _ST.get(pid))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--per', type=int, default=250); ap.add_argument('--procs', type=int, default=4)
    a = ap.parse_args()
    rng = random.Random(1)
    G = collections.defaultdict(list); N = collections.Counter(); L = collections.Counter()
    for run in ('baseline', 'ashallow', 'cross2', 'crossP5near', 'crossP5far'):
        for l in open(f'{ROOT}/{run}/trials.jsonl'):
            r = json.loads(l)
            if r['suite'] not in ('bench', 'reach'): continue
            arm = r['arm'] if run == 'baseline' else run
            N[arm] += 1
            if 'fin' in r: G[arm].append((r['pid'], r['s'], r['fin']))
            else: L[arm] += 1
    out = {}
    with ProcessPoolExecutor(a.procs) as ex:
        for arm, rows in sorted(G.items()):
            samp = rng.sample(rows, min(a.per, len(rows)))
            ds = sorted(ex.map(_d, samp, chunksize=4))
            q = lambda f: ds[min(len(ds) - 1, int(f * len(ds)))]
            out[arm] = dict(props=N[arm], lost=L[arm] / N[arm], kept=len(rows) / N[arm], n=len(ds), same=sum(d < 0.005 for d in ds) / len(ds),
                            q25=q(.25), q50=q(.5), q75=q(.75), q90=q(.9), far=sum(d > 0.05 for d in ds) / len(ds))
            o = out[arm]
            print(f"{arm:12s} props {o['props']:5d} lost {o['lost']:.2f} | kept children: same-basin {o['same']:.2f}  d q25/50/75/90 "
                  f"{o['q25']:.3f} {o['q50']:.3f} {o['q75']:.3f} {o['q90']:.3f}  >0.05 {o['far']:.2f}", flush=True)
    json.dump(out, open(f'{ROOT}/move_spread.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
