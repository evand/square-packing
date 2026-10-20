"""Mechanism A: noisy greedy best fit (Evan's sketch, 10-16).

Each placement: rasterise the covered region (squares + everything outside the box) at step h; for every angle in a
small palette (0 plus the dominant tilts of the squares near the uncovered region) compute by FFT convolution
    overlap(x, y)  = covered pixels inside the unit square at (x, y, angle)
    contact(x, y)  = covered pixels in a thin ring (width RING) around it
and the cheap score  ov - LAM * contact  (low = tight fit: no overlap, hugging walls / squares).
Then draw K random poses (uniform over uncovered pixels as centres, uniform angle from the palette) and keep the best
by the score ("best of K random options": K = 1 is pure random, large K is near-greedy).  With SNAP on, the chosen
pose slides to the local score minimum within SNAP_R at its angle (pushes it against its neighbours).  Paint it,
repeat until m squares.  Small overlaps are allowed (the harness's fq quench relaxes them); when nothing fits the
least-bad pose is taken.  Returns C candidates (different random draws) per call.

Knobs (env vars, read at call time): SHOVE_K (comma list, cycled over the candidates; default '64,16,16'), SHOVE_C (3), SHOVE_LAM (0.05), SHOVE_SNAP (1),
SHOVE_SNAP_R (0.3), SHOVE_H (0.04), SHOVE_TILT (1: use neighbourhood tilts; 0: uprights only),
SHOVE_PUSH (1: continuous SAT push-out of each placed pose), SHOVE_RANDANG (0: extra random angles per candidate).
"""
from __future__ import annotations
import math, os
import numpy as np
import harness as H  # noqa: F401  (kept: callers may use fill_greedy.H)
from pk.fill import (BESTOFK_DEFAULTS, Grid, palette as _palette, bestofk_candidate, _push, _descend,  # noqa: F401
                     _kernel, _corners, _sep, RING, PAD)

# 10-20: the engine moved to pk/fill.py (placer 'bestofK'); this file is the shove-bench adapter (env knobs, C
# candidates per call), same draws and output as the 10-16 version.


def _knobs():
    e = os.environ.get
    return dict(K=[int(v) for v in e('SHOVE_K', '64,16,16').split(',')], C=int(e('SHOVE_C', 3)), LAM=float(e('SHOVE_LAM', 0.05)),
                SNAP=int(e('SHOVE_SNAP', 1)), SNAP_R=float(e('SHOVE_SNAP_R', 0.3)), h=float(e('SHOVE_H', 0.04)),
                TILT=int(e('SHOVE_TILT', 1)), PUSH=int(e('SHOVE_PUSH', 1)),
                RANDANG=int(e('SHOVE_RANDANG', 0)))


def palette(S, sq, meta, tilt):
    return _palette(S, sq, tilt)


def one_candidate(S, sq, meta, m, rng, kn, angs):
    return bestofk_candidate(S, sq, m, rng, kn, angs)[0]


def fill(S, sq, meta, m, rng):
    kn = _knobs()
    angs = palette(S, sq, meta, kn['TILT'])
    out = []
    for c in range(kn['C']):
        kc = dict(kn, K=kn['K'][c % len(kn['K'])])
        A = angs + [rng.uniform(-45, 45) for _ in range(kn['RANDANG'])]
        out.append(one_candidate(S, sq, meta, m, rng, kc, A))
    return out
