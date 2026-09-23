#!/usr/bin/env python3
"""break_h1_filter: the filter of `tasks/break-h1/README.md` item 4.

Every restricted-dual bound (CW, H, HP, H1, H2) is a Farkas bound, so it is >= delta* for EVERY
configuration satisfying its rows; at a margin-0 configuration it must therefore be >= 0.
Anything below -1e-8 is a bug in the rows, not a discovery.

Points tested: the (4,1) band stack (T = 4, delta = 0 exactly), the T = 4 tiling-minus-
permutation-hole-set zeros with a tilted central block, and the T = 5 tiling-minus-hole-set
zeros produced by the unrestricted k = 15, 16 runs of `bandcut_k.py scan --T 5`.

Imports only; nothing in search/ is modified."""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import break_h1 as B                                             # noqa: E402
import t4_cycles as TCY                                          # noqa: E402

TOLBUG = -1e-8


def report(lab, z, n, T, lvl, eps=None, maxchains=120):
    o = B.bounds(np.asarray(z, dtype=float), n, T, lvl, tol=1e-9, maxchains=maxchains,
                 eps_rad=eps, kmax=7, do_cyc=(n <= 12), nch2=1, cap2=1200)
    vals = {k: o[k] for k in ('CW', 'H', 'HP', 'H1', 'H2', 'F') if k in o}
    if 'CYC' in o:
        vals['CYC'] = o['CYC']
    worst = min(v for v in vals.values() if v == v and v != math.inf)
    bug = worst < TOLBUG
    print(f"  {lab:46s} delta={lvl:+.3e}  " +
          '  '.join(f'{k} {v:+.3e}' for k, v in vals.items()) +
          f"   worst {worst:+.3e}  {'*** FILTER BREACH ***' if bug else 'ok'}", flush=True)
    return bug, worst, vals


def main():
    nbug = 0
    print('# T = 4 margin-0 points', flush=True)
    for r in TCY.load_recorded('.'):
        if r['delta'] < -1e-12:
            continue
        if r['src'] not in ('bandstack', 'bestbar') and r.get('eps') not in (1.0,):
            continue
        # the `p3 bar` is recorded with delta = -1.2e-02: it is NOT a margin-0 point and the
        # filter does not apply to it (its bounds equal delta*, which is correct and negative).
        b, _w, _v = report(f"{r['src']}: {r['cell']}", r['z'], 12, 4.0, r['delta'],
                           eps=(math.radians(r['eps']) if r.get('eps') else None))
        nbug += b
    print('\n# T = 4 tiling minus a permutation hole set, central block tilted (structured)',
          flush=True)
    rng = np.random.default_rng(5)
    for epsdeg in (1.0, 5.0):
        eps = math.radians(epsdeg)
        fams = B.dedup_theta(B.structured_T4(eps))
        done = 0
        for (th, seed, lab) in fams:
            if 'central2' not in lab:
                continue
            v, z, _hit = B.delta_star_at(th, 4.0, rng, seeds=[seed], limit=40, jitters=2,
                                         extra=40)
            if z is None or v < -1e-12:
                continue
            b, _w, _v = report(f'{lab} e{epsdeg:g}', z, 12, 4.0, v, eps=eps)
            nbug += b
            done += 1
            if done >= 4:
                break
    print('\n# T = 5 tiling minus a hole set, unrestricted (delta* = 0)', flush=True)
    p = 'runs/break_h1_t5_scan.jsonl'
    if os.path.exists(p):
        seen = 0
        for line in open(p):
            r = json.loads(line)
            if r.get('value') is None or r['pidx'] >= 0 or r['value'] < -1e-12:
                continue
            b, _w, _v = report(f"T5 free k={r['k']} eps={r['eps']:g}", r['z'], 20, 5.0,
                               r['value'], eps=math.radians(r['eps']), maxchains=80)
            nbug += b
            seen += 1
            if seen >= 4:
                break
    print(f'\n# filter breaches (any bound < {TOLBUG:g} at a margin-0 point): {nbug}', flush=True)


if __name__ == '__main__':
    main()
