#!/usr/bin/env python3
"""model_fast == model (rows bitwise, order, disjunctions) on kicked / overlapping states from saved minima."""
import glob, math, random, sys
import numpy as np
from rigid import load
import inc

files = sorted(glob.glob('data/s110_minima/*.txt'))[:8] + ['seeds/rec132.txt', 'data/s90_best_10.0095668.txt', 'seeds/rec71.txt']
rng = random.Random(1)
bad = 0; tot = 0
for f in files:
    s, sq = load(f)
    for sig in (0.0, 1e-4, 1e-3, 1e-2, 5e-2):
        q = [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), t + rng.gauss(0, sig)) for x, y, t in sq]
        for delta in (1e-7, 3e-3, 3e-2):
            r0, d0 = inc.model(s, q, delta)
            r1, d1 = inc.model_fast(s, q, delta)
            ok = all(np.array_equal(np.asarray(a), np.asarray(b)) for a, b in zip(r0, r1))
            ok = ok and [p for p, _ in d0] == [p for p, _ in d1] and all(
                all(np.array_equal(np.asarray(x), np.asarray(y)) for A, B in zip(a0, a1) for x, y in zip(A, B))
                for (_, a0), (_, a1) in zip(d0, d1)) and all(len(a0) == len(a1) for (_, a0), (_, a1) in zip(d0, d1))
            tot += 1
            if not ok:
                bad += 1; print('MISMATCH', f, sig, delta, len(r0[0]), len(r1[0]), len(d0), len(d1))
print(f'{tot - bad}/{tot} identical')
sys.exit(1 if bad else 0)
