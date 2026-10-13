#!/usr/bin/env python3
"""Removal positive control (PLAN-2026-10-12 item 4; 10-12).

Event: s(179) from Couzo's s(180) (10-05 register, 13.9278881405) minus one square, then "shrink-hopping" (SQUISH
#401: 13.8913126; frontier at the cutoff 13.8953410700).  Our time-travel chain gave 179 one 4-min round; its removal
seed quenched to 13.9224.  Here: every single-square removal x completion modes, all fq:
  plain   quench --loosen 1.0
  shrink  quench --shrink (fixed-side overlap relax + shrink hops, then polish)
  soft    soft-to-rigid: --round-r0 0.15 --round-compress
  kick    quench --kick 0.02 --loosen 1.01
Success: any output below the cutoff frontier by > 1e-8 (secondary: below SQUISH's 13.8913126).

  python3 removal_ctl.py run --src ../exact/batch/inputs/n-180.txt --target 13.89534106997649 --out runs/removal_ctl/r179
  python3 removal_ctl.py report runs/removal_ctl/r179 --target 13.89534106997649
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
FQ = os.path.join(HERE, 'target', 'release', 'fq')
from relax import OPS                    # one operator table (10-13)
MODES = dict(OPS, plain=['--loosen', '1.0'], soft=['--loosen', '1.0', '--round-r0', '0.15', '--round-compress'],
             kick=['--loosen', '1.01', '--kick', '0.02'])


def read(path):
    L = [l for l in open(path).read().split('\n') if l.strip()]
    from fractions import Fraction
    n, s = L[0].split()[:2]
    sq = [[float(t) for t in l.split()[:3]] for l in L[1:1 + int(n)]]
    return float(Fraction(s)), sq


def job(t):
    inp, outp, mode, k = t
    r = subprocess.run([FQ, 'quench', '--in', inp, '--out', outp, '--seed', str(k + 1)] + MODES[mode],
                       capture_output=True, text=True, timeout=1200)
    try:
        d = json.loads(r.stdout.strip().splitlines()[-1])
        return dict(k=k, mode=mode, s=d['s'], ok=d['min_gap'] > -1e-9 and d['min_wall'] > -1e-9, out=outp, sec=d['sec'])
    except Exception:
        return dict(k=k, mode=mode, err=r.stderr[-300:])


def run(a):
    from pk import manifest
    import scratch
    os.makedirs(a.out, exist_ok=True)
    m = manifest.start(a.out, 'removal_ctl', args=vars(a), procs=a.procs)
    scratch.scope('removal_ctl')
    s, sq = read(a.src)
    tasks = []
    for k in range(len(sq)):
        inp = os.path.join(a.out, f'rm{k:03d}.txt')
        rest = [q for i, q in enumerate(sq) if i != k]
        with open(inp, 'w') as f:
            f.write(f'{len(rest)} {s!r}\n' + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in rest))
        for mode in a.modes.split(','):
            tasks.append((inp, os.path.join(a.out, f'out_{mode}_{k:03d}.txt'), mode, k))
    with ProcessPoolExecutor(a.procs) as ex, open(os.path.join(a.out, 'results.jsonl'), 'w') as f:
        for r in ex.map(job, tasks, chunksize=1):
            f.write(json.dumps(r) + '\n'); f.flush()
    m.finish(summary=dict(tasks=len(tasks)))


def report(a):
    R = [json.loads(l) for l in open(os.path.join(a.dir, 'results.jsonl'))]
    R = [r for r in R if r.get('ok')]
    for mode in sorted({r['mode'] for r in R}):
        rs = sorted((r for r in R if r['mode'] == mode), key=lambda r: r['s'])
        print(f'{mode:7s} {len(rs)} runs; best {rs[0]["s"]:.10f} (k={rs[0]["k"]}); below target: '
              f'{sum(r["s"] < a.target - 1e-8 for r in rs)}; best 5: ' + ' '.join(f'{r["s"]:.6f}' for r in rs[:5]))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--src', required=True); r.add_argument('--out', required=True)
    r.add_argument('--target', type=float, required=True); r.add_argument('--procs', type=int, default=8)
    r.add_argument('--modes', default='plain,shrink,soft,kick')
    q = sp.add_parser('report'); q.add_argument('dir'); q.add_argument('--target', type=float, required=True)
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else report(a)
