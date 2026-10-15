#!/usr/bin/env python3
"""Reach test (10-15): soft-to-rigid re-hardening rounds (harden.py) from the removal-bench starts, per F event, target =
frontier_at(n, T) (the ante).  All events concurrently, --procs each.  Score with `removal`-style tiers.

  python3 reach_harden.py run --bench runs/removal_bench/b1 --out runs/reach/h1 [--procs 1]
  python3 reach_harden.py score --out runs/reach/h1
"""
import argparse, glob, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))

def run(a):
    P = []
    for rj in sorted(glob.glob(f'{a.bench}/*/result.json') + glob.glob(f'{a.bench}/*/event.json')):
        ev = json.load(open(rj)); d = os.path.dirname(rj)
        if 'new' not in ev and 'S' in ev:                  # comp_l1 event.json
            ev['new'] = ev['S']; ev.setdefault('T', ev.get('T', ''))
        starts = sorted(glob.glob(f'{d}/start_*.txt') + glob.glob(f'{d}/start[0-9]*.txt'))[:a.starts]
        if not starts: continue
        out = f'{a.out}/{os.path.basename(d)}'; os.makedirs(out, exist_ok=True)
        cmd = [sys.executable, f'{HERE}/harden.py', 'run', '--files', *starts, '--out', out, '--target', repr(ev['ante']),
               '--grid', a.grid, '--seeds', str(a.seeds), '--rounds', str(a.rounds), '--keep', '4', '--procs', str(a.procs)]
        json.dump(ev, open(f'{out}/event.json', 'w'))
        P.append(subprocess.Popen(cmd, stdout=open(f'{out}/log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
    for p in P: p.wait()

def score(a):
    for ej in sorted(glob.glob(f'{a.out}/*/event.json')):
        ev = json.load(open(ej)); d = os.path.dirname(ej); best = None
        for f in glob.glob(f'{d}/**/*.jsonl', recursive=True) + glob.glob(f'{d}/**/results*.json', recursive=True):
            for l in open(f):
                try: r = json.loads(l)
                except Exception: continue
                for rr in (r if isinstance(r, list) else [r]):
                    s = rr.get('s_pol') or rr.get('s') if isinstance(rr, dict) else None
                    if isinstance(s, (int, float)) and (best is None or s < best): best = s
        tier = 'REPLICATED' if best is not None and best <= ev['new'] + 1e-9 else 'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
        print(f"n={ev['n']:3d} {ev['T'][:10]} harden best {best}  explorer best {ev.get('best')}  ante {ev['ante']:.8f} posted {ev['new']:.8f} {tier}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--bench', required=True); r.add_argument('--out', required=True)
    r.add_argument('--procs', type=int, default=1); r.add_argument('--starts', type=int, default=4)
    r.add_argument('--grid', default='r0=0.1,0.2,0.3;noise=0.004,0.008'); r.add_argument('--seeds', type=int, default=3)
    r.add_argument('--rounds', type=int, default=3)
    s = sp.add_parser('score'); s.add_argument('--out', required=True)
    a = ap.parse_args(); {'run': run, 'score': score}[a.cmd](a)
