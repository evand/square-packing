#!/usr/bin/env python3
"""Is the polish worth it per hop step? (10-07, Evan's question.)  Hop proposals (hop.propose) from the register packing;
per quench record the polished side, the ALM-only side (raw, and repaired to feasible), and the time split.

  polish_value.py --n 130 180 237 263 --props 40 --procs 15 --out runs/pv0
"""
import argparse, json, os, random, subprocess, tempfile
from jobpool import run_jobs
import mcmin, hop


def job(j):
    n, k, loosen = j
    rng = random.Random(n * 7919 + k)
    s0, sq0 = mcmin.load_deg(f'{hop.BATCH}/n-{n}.txt')
    prop, kind = hop.propose(sq0, rng)
    tmp = tempfile.mkdtemp()
    a, b = os.path.join(tmp, 'a.txt'), os.path.join(tmp, 'b.txt')
    mcmin.write_deg(a, s0, prop)
    r = subprocess.run([hop.FQ, 'quench', '--in', a, '--out', b, '--loosen', str(loosen)], capture_output=True, text=True, timeout=900)
    d = json.loads(r.stdout)
    d.update(n=n, k=k, kind=kind, loosen=loosen, start=s0)
    return d


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+'); ap.add_argument('--props', type=int, default=40)
    ap.add_argument('--loosen', type=float, nargs='+', default=[1.02])
    ap.add_argument('--procs', type=int, default=15); ap.add_argument('--out', default='runs/pv0')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    res = [r for r in run_jobs(job, [(n, k, l) for n in a.n for l in a.loosen for k in range(a.props)], procs=a.procs, timeout=1200) if r]
    json.dump(res, open(f'{a.out}/pv.json', 'w'), indent=0)
