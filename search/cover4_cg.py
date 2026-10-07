#!/usr/bin/env python3
"""Exact restricted-master column generation for the packing lower bound at a container side t.

`search/dual_exact.py build` snaps a float pose list to rationals, solves ONE LP over *all*
arrangement vertices of that pose set, and certifies the result exactly.  Its value is only as
good as the columns it is given, and the arrangement (hence the LP) grows quadratically in the
number of poses, so a large float pose pool cannot be fed to it in one piece.

This script runs the obvious restricted master: keep a current support, add a batch of candidate
poses from the pool, call `dual_exact.py build` on the union, keep the poses that come out with
positive mass, repeat.  Every step's `L` is an exact Fraction produced by `dual_exact.py`; adding
columns can *lower* it (the new squares' arrangement vertices are extra, genuine `cov <= 1`
constraints that the previous LP had not seen), so the best step is kept and each certified value
stands on its own.  Nothing here is load-bearing for the proof: it only decides which poses the
exact certification is run on.

Usage
    python3 search/cover4_cg.py pool OUT [--Q 10000000] [--Dc 100000000] [--scale LAM] FILE...
        union of the pose lists of several support files (float `pose cx cy theta_deg mu` or exact
        `pose p q cx cy mass`), deduplicated on the snapped pose, centres optionally scaled by LAM
        (to reuse a support found at another container side)
    python3 search/cover4_cg.py cg POOL START TAG [--t 4] [--batch 130] [--rounds 10] [--Q ...]
        the restricted master; writes runs/dual_exact_<TAG>_best_support.txt (the best step) and
        prints the exact L of every step
"""
import sys, os, math, argparse, subprocess, shutil
from math import gcd
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RUNS = os.path.join(REPO, 'runs')
PY = sys.executable


def read_poses(path):
    """(cx, cy, theta_deg, mass) from a float or an exact support file"""
    out = []
    for line in open(path):
        q = line.split()
        if not q or q[0] != 'pose':
            continue
        if len(q) == 5:                                        # float: pose cx cy theta_deg mu
            out.append((float(q[1]), float(q[2]), float(q[3]), float(q[4])))
        else:                                                  # exact: pose p q cx cy mass
            p, qq, cx, cy, m = int(q[1]), int(q[2]), Fr(q[3]), Fr(q[4]), Fr(q[5])
            out.append((float(cx), float(cy), math.degrees(2 * math.atan2(p, qq)), float(m)))
    return out


def write_poses(path, rows, note=''):
    with open(path, 'w') as f:
        f.write(f"# {len(rows)} poses {note}\n")
        for cx, cy, th, m in rows:
            f.write(f"pose {cx:.13f} {cy:.13f} {th:.11f} {m:.12f}\n")


def key(r, Q, Dc):
    """the pose as dual_exact.py will snap it, so the pool holds exactly the distinct columns"""
    p = round(math.tan(math.radians(r[2]) / 2) * Q)
    g = gcd(abs(p), Q)
    return (p // g, Q // g, round(r[0] * Dc), round(r[1] * Dc))


def cmd_pool(a):
    seen = {}
    for path in a.files:
        for cx, cy, th, m in read_poses(path):
            cx *= a.scale; cy *= a.scale
            k = key((cx, cy, th, m), a.Q, a.Dc)
            if k not in seen or seen[k][3] < m:
                seen[k] = (cx, cy, th, m)
    rows = sorted(seen.values(), key=lambda r: -r[3])
    write_poses(a.out, rows, f"(Q={a.Q}, Dc={a.Dc}, scale={a.scale}) from: {' '.join(a.files)}")
    print(f"{len(rows)} distinct poses -> {a.out}")


def cmd_cg(a):
    cur = read_poses(a.start)
    pool = read_poses(a.pool)
    have = {key(r, a.Q, a.Dc) for r in cur}
    rest = [r for r in pool if key(r, a.Q, a.Dc) not in have]
    print(f"start {len(cur)} poses, pool {len(pool)}, {len(rest)} new candidates", flush=True)
    best = None
    for rd in range(a.rounds):
        if not rest:
            break
        take, rest = rest[:a.batch], rest[a.batch:]
        src = os.path.join(RUNS, f'cg_{a.tag}_{rd}_src.txt')
        write_poses(src, cur + take)
        r = subprocess.run([PY, os.path.join(HERE, 'dual_exact.py'), 'build', '--t', a.t,
                            '--tag', f'cg{a.tag}{rd}', '--src', src, '--Q', str(a.Q),
                            '--Dc', str(a.Dc), '--procs', str(a.procs)],
                           capture_output=True, text=True, cwd=REPO)
        line = [l for l in (r.stdout + r.stderr).splitlines() if l.strip().startswith('EXACT:')]
        if not line:
            print(f"round {rd}: build failed\n{(r.stdout + r.stderr)[-2000:]}", flush=True)
            break
        L = Fr(line[0].split('L = mass/M = ')[1].split(' = ')[0])
        sup = os.path.join(RUNS, f'dual_exact_cg{a.tag}{rd}_support.txt')
        cur = read_poses(sup)
        print(f"round {rd}: +{len(take)} candidates -> L = {L} = {float(L):.9f}, support {len(cur)}", flush=True)
        if best is None or L > best:
            best = L
            shutil.copy(sup, os.path.join(RUNS, f'dual_exact_{a.tag}_best_support.txt'))
    print(f"BEST L = {best} = {float(best):.9f}" if best else "no result", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('pool')
    p.add_argument('out'); p.add_argument('files', nargs='+')
    p.add_argument('--Q', type=int, default=10 ** 7); p.add_argument('--Dc', type=int, default=10 ** 8)
    p.add_argument('--scale', type=float, default=1.0, help='multiply the centres by this (reuse another side)')
    c = sub.add_parser('cg')
    c.add_argument('pool'); c.add_argument('start'); c.add_argument('tag')
    c.add_argument('--t', default='4'); c.add_argument('--batch', type=int, default=130)
    c.add_argument('--rounds', type=int, default=10); c.add_argument('--procs', type=int, default=4)
    c.add_argument('--Q', type=int, default=10 ** 7); c.add_argument('--Dc', type=int, default=10 ** 8)
    a = ap.parse_args()
    (cmd_pool if a.cmd == 'pool' else cmd_cg)(a)


if __name__ == '__main__':
    main()
