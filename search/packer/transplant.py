#!/usr/bin/env python3
"""Gadget transplant (10-09): tilted squares of a source packing (any n) into a box of side s for target n, axis part
filled exactly (gen.mis_fill2).  Scans s and the 8 D4 placements (gadget anchored at the lower-left corner of the source
box, after the symmetry), reports the axis + tilted count per s, and quenches every fill reaching n with fq.

  transplant.py --src FILE --n N --smin A --smax B [--steps K] [--out DIR] [--procs P]

Motivation: the 2026 records at 88 (SQUISH) and 146 (ry-xu) use ~37 near-45-degree squares where the old 45-degree
constructions at 89 / 109 / 147 use 49-64; their gadgets at n + 1 test whether the hint carries over.
"""
import argparse, os, tempfile
import numpy as np
from jobpool import run_jobs
import mcmin, hop
from gen import is_axis, mis_fill2
from mkseed import D4, D4a


def job(j):
    src, n, s, sym, out = j
    ss, sq = mcmin.load_deg(src)
    t = []
    for x, y, a in sq:
        if is_axis((x, y, a), 1e-3):
            continue
        x, y = D4[sym](x, y, ss)
        t.append((x, y, (a * D4a[sym]) % 90))
    t = [q for q in t if 0.5 <= q[0] <= s - 0.5 and 0.5 <= q[1] <= s - 0.5]
    try:
        ax = mis_fill2(s, t, tol=1e-6)
    except Exception as e:
        return (s, sym, len(t), -1, None, repr(e)[:80])
    tot = len(t) + len(ax)
    if tot < n:
        return (s, sym, len(t), tot, None, '')
    allsq = t + [tuple(q[:3]) if len(q) >= 3 else (q[0], q[1], 0.0) for q in ax][:n - len(t)]
    tmp = tempfile.mkdtemp()
    hop.EXTRA[:] = ['--loosen', '1.0']
    r = hop.quench(s, allsq, tmp)
    if r:
        p = f'{out}/tp_n{n}_sym{sym}_s{s:.5f}_{r[0]:.9f}.txt'
        mcmin.write_deg(p, r[0], r[1])
        return (s, sym, len(t), tot, r[0], p)
    return (s, sym, len(t), tot, None, 'quench failed')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True); ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--smin', type=float, required=True); ap.add_argument('--smax', type=float, required=True)
    ap.add_argument('--steps', type=int, default=8); ap.add_argument('--out', default='runs/tp'); ap.add_argument('--procs', type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = [(a.src, a.n, float(s), sym, a.out) for s in np.linspace(a.smin, a.smax, a.steps) for sym in range(8)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=600)
    for r in sorted((r for r in res if isinstance(r, tuple)), key=lambda r: (r[4] is None, r[4] or 0, -r[3])):
        print(f's={r[0]:.5f} sym={r[1]} tilted={r[2]} total={r[3]} quench={r[4]} {r[5]}')
