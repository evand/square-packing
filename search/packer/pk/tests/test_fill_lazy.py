"""bestofK placement: lazily tiled score maps give the same placements as the whole-grid maps (10-27)."""
import random
import numpy as np
from pk import fill as FL


def _case(seed):
    rng = random.Random(seed)
    S = 7.6
    sq = np.array([[i % 7 + 0.55, i // 7 + 0.55, 0.0] for i in range(49)], float)
    sq[rng.sample(range(49), 6), 2] = 12.0                   # a few tilted squares: two palette angles
    rm = np.zeros(49, bool); rm[rng.sample(range(49), 8)] = True
    keep = sq[~rm]; region = [(x - .5, y - .5, x + .5, y + .5) for x, y in sq[rm, :2]]
    return S, keep, int(rm.sum()), region


def test_lazy_equals_full():
    for seed in range(4):
        S, keep, m, region = _case(seed)
        kn = dict(FL.BESTOFK_DEFAULTS, K=8, m_max=10 ** 6)
        angs = FL.palette(S, keep, kn['TILT'])
        out = []
        for lazy in (0, 1):
            kn['LAZY'] = lazy
            X, _ = FL.bestofk_candidate(S, keep, m, random.Random(seed), kn, angs, region=region, p_region=0.75, ov_full=0.02)
            out.append(X)
        assert np.array_equal(out[0], out[1])
