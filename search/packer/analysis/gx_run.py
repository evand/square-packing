#!/usr/bin/env python3
"""GX: grid-crossing descents from structurally foreign above-grid minima (10-15; seeds from analysis/gx_select.py).

run:   one explorer per seed (explore.py, 1 proc, --minutes; frontier = the grid side k, so value reward = distance below
       the grid), --cores at a time.
score: per seed: best polished side; CROSSED if < k; for every sub-k archive entry its (rattler-snapped) corner cell:
       KNOWN (a cell the below-grid census family occupies) or NEW (a sub-grid cell never seen below the grid before);
       plus whether the 40-50 degree "junction" count (explore.roles_j J) of the seed survives in the best descendant.

  python3 analysis/gx_run.py run --seeds runs/gx110/seeds.json --extra runs/c110_h2/b02882.txt --out runs/gx110/x1 [--minutes 15]
  python3 analysis/gx_run.py score --out runs/gx110/x1
"""
import argparse, glob, json, math, os, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)


def run(a):
    from pk import manifest
    d = json.load(open(a.seeds)); n = d['n']; k = math.ceil(math.sqrt(n) - 1e-12)
    seeds = [r['path'] for r in d['seeds']] + list(a.extra or [])
    os.makedirs(a.out, exist_ok=True)
    json.dump(dict(n=n, k=k, seeds=seeds, below_cells=d['below_cells']), open(f'{a.out}/plan.json', 'w'), indent=0)
    man = manifest.start(a.out, 'gx_run', args=vars(a), procs=a.cores)
    Q, P = list(enumerate(seeds)), []
    while Q or P:
        P = [p for p in P if p.poll() is None]
        while Q and len(P) < a.cores:
            i, s = Q.pop(0)
            cmd = [sys.executable, f'{PK}/explore.py', '--n', str(n), '--adapt', '--reward', 'value', '--elite-share', '0.3',
                   '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(float(k)), '--starts', s,
                   '--minutes', str(a.minutes), '--procs', '1', '--seed', str(1000 + i), '--out', f'{a.out}/s{i:03d}']
            P.append(subprocess.Popen(cmd, stdout=open(f'{a.out}/s{i:03d}.log', 'w'), stderr=subprocess.STDOUT, cwd=PK))
        time.sleep(5)
    man.finish()
    score(a)


def _load(p):
    with open(p) as fh:
        n, s = fh.readline().split()[:2]
        return float(s), np.array([list(map(float, x.split()[:3])) for x in fh if x.strip()])


def score(a):
    import cornersig as CS
    import explore as X
    plan = json.load(open(f'{a.out}/plan.json')); k = plan['k']; known = {tuple(c) for c in plan['below_cells']}
    tally = dict(seeds=0, crossed=0, new_cell=0, junction_kept=0)
    for i, seed in enumerate(plan['seeds']):
        ad = f'{a.out}/s{i:03d}/archive.jsonl'
        if not os.path.exists(ad):
            continue
        tally['seeds'] += 1
        E = [json.loads(l) for l in open(ad)]
        E = [e for e in E if not e.get('unpolished')]
        s0, q0 = _load(seed); g0 = CS.signature(s0, q0); j0 = X.roles_j(q0)[3]
        best = min(E, key=lambda e: e['s']) if E else None
        sub = sorted([e for e in E if e['s'] < k - 1e-9], key=lambda e: e['s'])
        cells, newc = {}, []
        for e in sub[:40]:
            s1, q1 = _load(e['path']); g = CS.signature(s1, q1)
            c = (g['sorted'][0], g['sorted'][1], g['sorted'][2] + g['sorted'][3])
            cells[c] = min(cells.get(c, 9e9), e['s'])
            if c not in known:
                newc.append((e['s'], c, e['path']))
        jb = X.roles_j(_load(best['path'])[1])[3] if best else None
        crossed = bool(sub)
        tally['crossed'] += crossed; tally['new_cell'] += bool(newc); tally['junction_kept'] += bool(crossed and jb is not None and abs(jb - j0) <= 1 and j0 >= 3)
        print(f"seed {i:2d} s0 {s0:.6f} cell {g0['sorted']} J{j0} -> best {best['s'] if best else float('nan'):.6f} J{jb} "
              f"{'CROSSED' if crossed else 'above'}  sub-grid cells {len(cells)} new {len(newc)}"
              + (f"  NEW: {newc[0][1]} at {newc[0][0]:.6f} {newc[0][2]}" if newc else ''))
    print(tally)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--seeds', required=True); r.add_argument('--extra', nargs='*')
    r.add_argument('--out', required=True); r.add_argument('--minutes', type=float, default=15); r.add_argument('--cores', type=int, default=16)
    s = sp.add_parser('score'); s.add_argument('--out', required=True)
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
