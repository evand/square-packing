#!/usr/bin/env python3
"""Where did archived melt moves melt?  (10-21; Evan: "Of the successful melt moves we've seen ... is it the junction?")

Per archived melt entry (kind 'melt', explorer archives): parent and child packings; displacement per child square =
distance to the nearest parent square of similar angle (mod 90, 3 deg); melted region = the 12 most displaced squares,
centre = their centroid.  Junction point (two-band packings only: >= 10 squares tilted > +5 deg and >= 10 < -5 deg) =
midpoint of the closest (+, -) tilted pair in the parent.  Statistic: rank of the melt centre's distance to the junction
among all parent squares' distances (0 = nearest square to the junction, 1 = farthest); a uniform random centre gives
U(0, 1).  Also the shim cluster (0.5-5 deg) distance.

  python3 analysis/melt_junction.py RUNDIR... [--out runs/melt_junction.jsonl]
"""
import argparse, glob, json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from pk.packing import read
from pk import archive


def tilt(sq):
    a = sq[:, 2] % 90.0
    return np.where(a < 45, a, a - 90)


def junction(sq):
    t = tilt(sq)
    P, N = sq[t > 5, :2], sq[t < -5, :2]
    if len(P) < 10 or len(N) < 10:
        return None
    d = np.linalg.norm(P[:, None] - N[None], axis=-1)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return (P[i] + N[j]) / 2


def melt_centre(par, ch, k=12):
    tp, tc = par[:, 2] % 90.0, ch[:, 2] % 90.0
    d = np.linalg.norm(ch[:, None, :2] - par[None, :, :2], axis=-1)
    da = np.abs(((tc[:, None] - tp[None]) + 45) % 90 - 45)
    d = np.where(da < 3, d, d + 5)
    disp = d.min(1)
    top = np.argsort(-disp)[:k]
    return ch[top, :2].mean(0), float(np.median(disp[top])), float(np.median(disp))


def load(path):
    try:
        archive.ensure(path)
    except Exception:
        pass
    p = read(path)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+'); ap.add_argument('--out', default=None); ap.add_argument('--kind', default='melt')
    a = ap.parse_args()
    rows = []
    for f in a.runs:
        for arch in glob.glob(os.path.join(f, '**', 'archive.jsonl'), recursive=True):
            E = {}
            for l in open(arch):
                try:
                    e = json.loads(l); E[e['i']] = e
                except Exception:
                    pass
            for e in E.values():
                if e.get('kind') != a.kind or e.get('parent', -1) not in E:
                    continue
                pe = E[e['parent']]
                try:
                    ch, par = load(e['path']), load(pe['path'])
                except Exception:
                    continue
                if ch is None or par is None or ch.n != par.n:
                    continue
                J = junction(par.sq)
                if J is None:
                    continue
                c, dtop, dmed = melt_centre(par.sq, ch.sq)
                dist = np.linalg.norm(par.sq[:, :2] - J, axis=1)
                dc = float(np.linalg.norm(c - J))
                rank = float((dist < dc).mean())
                rows.append(dict(run=os.path.dirname(arch), n=par.n, i=e['i'], ds=float(e['s'] - pe['s']),
                                 s=float(e['s']), dc=dc, rank=rank, dtop=dtop, side=float(par.s)))
    if a.out:
        with open(a.out, 'w') as fo:
            for r in rows:
                fo.write(json.dumps(r) + '\n')
    R = np.array([r['rank'] for r in rows]); D = np.array([r['ds'] for r in rows])
    print(f'{len(rows)} {a.kind} entries in two-band packings, n: {sorted(set(r["n"] for r in rows))[:40]}')
    for name, m in (('all archived', np.ones(len(rows), bool)), ('improved (ds < -1e-9)', D < -1e-9),
                    ('improved > 1e-4', D < -1e-4), ('worse / equal', D >= -1e-9)):
        if m.sum():
            q = np.quantile(R[m], [.25, .5, .75])
            print(f'{name:24s} {m.sum():5d}  rank quartiles {q.round(2)}  share in nearest 10 % {np.mean(R[m] < .1):.2f}  '
                  f'nearest 25 % {np.mean(R[m] < .25):.2f}  (uniform: 0.10 / 0.25)')


if __name__ == '__main__':
    main()
