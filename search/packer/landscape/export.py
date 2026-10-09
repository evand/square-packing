#!/usr/bin/env python3
"""Pack the explainer data into one JSON (10-08).

  export.py --census-rerun RERUN.json --old OLD.json --out runs/landscape110/landscape110.json

Coordinates are integers: x, y in 1e-4 units, theta in 1e-2 degrees.  Sides as strings (12 decimals)."""
import argparse, glob, json, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
import morph, forces

RUN = os.path.join(os.path.dirname(HERE), 'runs/landscape110')


def enc(sq):
    out = []
    for x, y, t in sq:
        out += [round(x * 1e4), round(y * 1e4), round(math.degrees(t) * 100)]
    return out


def nflat(st):
    m = re.search(r'modulo (\d+)', st or '') or re.search(r'PSD with (\d+) zero', st or '')
    return int(m.group(1)) if m else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--census-rerun', required=True); ap.add_argument('--old', required=True)
    ap.add_argument('--out', default=f'{RUN}/landscape110.json')
    a = ap.parse_args()
    cen = []
    for o in json.load(open(a.census_rerun)):
        if o['new']['cls'] != 'certified':
            continue
        s, sq = morph.load(o['old']['path'])
        cen.append(dict(S=o['new']['S'][:14], roles=o['old'].get('roles'), flat=nflat(o['new']['status']), sq=enc(sq)))
    for o in json.load(open(a.old)):
        s, sq = morph.load(o['path'])
        cen.append(dict(S=o['S'][:14], roles=None, flat=nflat(o['status']), sq=enc(sq)))
    cen.sort(key=lambda r: r['S'])
    R = json.load(open(f'{RUN}/reps.json'))
    reps = []
    for r in R:
        f = forces.forces(r['path'])
        reps.append(dict(S=r['S'][:14], roles=r.get('roles'), sq=enc(f['sq']), free=f['free'], flat=f['flat'],
                         contacts=[[c[0], c[1], c[2], round(c[3] * 1e4), round(c[4] * 1e4), round(c[5], 5)] for c in f['contacts']]))
    cyc = json.load(open(f'{RUN}/cycle.json'))
    legs = [dict(i=l['i'], j=l['j'], ok=l.get('ok'), s=l.get('s'),
                 frames=[[round(s_, 6), enc(q)] for s_, q in l.get('frames', [])]) for l in cyc['legs']]
    settle = []
    for f in sorted(glob.glob(f'{RUN}/settle/*.json')):
        d = json.load(open(f))
        tr = d['traj']
        # thin to <= 160 frames, evenly in log(side excess over the final side)
        k = 160
        if len(tr) > k:
            idx = sorted(set(np.linspace(0, len(tr) - 1, k).round().astype(int).tolist()))
            tr = [tr[i] for i in idx]
        settle.append(dict(kind=d['kind'], seed=d['seed'], sigma=d.get('sigma'), s_final=d['s_final'], jammed=d['jammed'],
                           frames=[[s_, enc(q)] for s_, q in tr]))
    out = dict(n=110, census=cen, reps=reps, tour=cyc['tour'], W=cyc['W'], Sb=cyc['Sb'], legs=legs, settle=settle)
    json.dump(out, open(a.out, 'w'), separators=(',', ':'))
    print(a.out, os.path.getsize(a.out) // 1024, 'KB;', len(cen), 'census,', len(reps), 'reps,', len(legs), 'legs,', len(settle), 'settle')


if __name__ == '__main__':
    main()
