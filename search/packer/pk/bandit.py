"""Thompson-sampling bandits over arms with a reward rate per CPU-second (Gamma posterior, decayed counts).

Moved verbatim from explore.py (10-13) so every search loop uses one implementation: explore.py (move kinds, parent
cells), compose2.py (piece-pair groups).  An arm is any hashable key; reward may be graded (fractional), cost is seconds.
"""
from __future__ import annotations
import collections, time


class KindBandit:
    """Thompson sampling over move kinds: reward = new below-k archive entry; rate per CPU-second with a Gamma posterior;
    counts decay with half-life `half` seconds (rates drift as the archive fills); each kind keeps a floor probability."""
    def __init__(self, kinds, half=600.0, floor=0.03, a0=1.0, b0=60.0):
        self.kinds = list(kinds); self.half, self.floor, self.a0, self.b0 = half, floor, a0, b0
        self.r = {k: 0.0 for k in kinds}; self.c = {k: 0.0 for k in kinds}; self.last = time.time()

    def decay(self):
        now = time.time(); f = 0.5 ** ((now - self.last) / self.half); self.last = now
        for k in self.kinds:
            self.r[k] *= f; self.c[k] *= f

    def update(self, kind, sec, reward):
        self.decay(); self.c[kind] += sec; self.r[kind] += reward

    def choose(self, rng):
        if rng.random() < self.floor * len(self.kinds):
            return rng.choice(self.kinds)
        draw = {k: rng.gammavariate(self.a0 + self.r[k], 1.0) / (self.b0 + self.c[k]) for k in self.kinds}
        return max(draw, key=draw.get)

    def summary(self):
        return {k: f'{self.r[k]:.1f}/{self.c[k]:.0f}s' for k in self.kinds}


class CellBandit:
    """Thompson sampling over descriptor cells for parent choice (reward: a child that is a new below-k basin; cost: CPU s),
    decayed like KindBandit; with probability `uniform` the cell is drawn uniformly (exploration of new cells)."""
    def __init__(self, half=900.0, uniform=0.3, a0=1.0, b0=60.0):
        self.half, self.uniform, self.a0, self.b0 = half, uniform, a0, b0
        self.r = collections.Counter(); self.c = collections.Counter(); self.last = time.time()

    def update(self, cell, sec, reward):
        now = time.time(); f = 0.5 ** ((now - self.last) / self.half); self.last = now
        for d in (self.r, self.c):
            for key in d: d[key] *= f
        self.c[cell] += sec; self.r[cell] += reward

    def choose(self, cells, rng):
        if rng.random() < self.uniform:
            return rng.choice(cells)
        draw = {c: rng.gammavariate(self.a0 + self.r[c], 1.0) / (self.b0 + self.c[c]) for c in cells}
        return max(draw, key=draw.get)
