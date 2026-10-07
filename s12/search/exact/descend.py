#!/usr/bin/env python3
"""Escape false jams of the f64 squeeze at corner-corner touches, then hand over to exactsolve.

slp2 linearises each pair on one fixed separating axis.  At a corner-corner touch the true constraint is a disjunction
(either side line may separate), so slp2 can stall at a point that is not a local minimum.  Loop:
  MILP (exactsolve.vv_branch_test): min first-order dS with corner-corner touches as disjunctions;
  if dS* < -thr: kick along the MILP direction (step t, then slp2.repair), then re-run slp2 from there;
  until dS* >= -thr or no progress.

  descend.py in.txt out.txt [--kick 1e-5] [--rounds 20]
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, math, argparse, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'packer'))   # slp2 (local optimizer)
import numpy as np
import scipy.optimize as so
from mpmath import mp
import exactsolve as E
import slp2
from slp import save

_cap = {}
_orig_milp = so.milp


def _milp(*a, **k):
    r = _orig_milp(*a, **k)
    _cap['x'] = r.x
    return r


so.milp = _milp


def milp_dir(s, sq, tol=1e-6):
    P = E.Pack(len(sq), mp.mpf(s), [mp.mpf(x) for x, _, _ in sq], [mp.mpf(y) for _, y, _ in sq],
               [mp.mpf(t) for _, _, t in sq])
    c, vv = E.find_contacts(P, tol)
    dS, nvv, _ = E.vv_branch_test(P, c, vv, tol)
    return dS, (_cap['x'][:3 * len(sq) + 1] if _cap.get('x') is not None else None)


def descend(s, sq, kick=1e-5, rounds=20, thr=1e-6, log=print, tol=1e-6):
    mp.dps = 30
    for rd in range(rounds):
        dS, z = milp_dir(s, sq, tol if rd == 0 else 1e-6)
        log(f'round {rd}: S = {s:.13f}, MILP first-order dS* = {dS:.3e}')
        if dS is None or dS > -thr:
            break
        best = None
        for t in (kick, kick / 10, kick / 100):
            sq2 = [(x + t * z[3 * i], y + t * z[3 * i + 1], th + t * z[3 * i + 2]) for i, (x, y, th) in enumerate(sq)]
            s2, sq2 = slp2.repair(s + t * z[-1], sq2)
            if s2 < s - 1e-15:
                best = (s2, sq2)
                break
        if best is None:
            log('  kick does not decrease S after repair; stop')
            break
        s1, sq1 = best
        t0 = time.time()
        s, sq, ds0 = slp2.slp2(s1, sq1)
        log(f'  kick -> {s1:.13f}, slp2 -> {s:.13f} ({time.time() - t0:.0f}s)')
    return s, sq


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('input'); ap.add_argument('output')
    ap.add_argument('--kick', type=float, default=1e-5)
    ap.add_argument('--rounds', type=int, default=20)
    ap.add_argument('--tol', type=float, default=1e-6, help='contact tolerance of the first MILP (1e-10 for an exact input)')
    a = ap.parse_args()
    s, sq = slp2.load(a.input)
    s2, sq2 = descend(s, sq, a.kick, a.rounds, tol=a.tol)
    save(a.output, s2, sq2)
    print(f'{a.input}: {s:.13f} -> {s2:.13f}')
