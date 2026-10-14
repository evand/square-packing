"""Region clusters of the packings at one n (10-14; definition and evidence: ../analysis/similarity110.md, A2).

same packing   side key: sides within SAME (2e-9) are one minimum (the only key that works: contact fingerprints split
               91-99.9 % of same pairs; exact twins < 1e-9 apart are geometrically indistinguishable and harmless)
cluster        greedy leader balls in matching distance (rediscover.distance: assignment RMS of centre + angle moves over the
               8 box symmetries), radius r(n) = 0.42 / sqrt(n) (0.040 at 110), side keys visited deepest first; image
               fingerprint (pk/fingerprint) < 2 r as the prefilter.  At 110: 40 clusters, basins A / B apart at every
               radius 0.01-0.05, 72 % of explorer steps stay in a cluster (11.8x a shuffled null).

  labels(n, rows) -> {store id: cluster label}     rows = [(id, s, coords n x 3)], any subset (an as-of pool)
Cached per (n, cache tag) in runs/clusters/ (leaders with their coordinates + every labelled id): new rows are assigned
incrementally; delete the cache file to recompute deepest-first (orders agree at ARI 0.7-0.8, A2).
"""
from __future__ import annotations
import bisect, json, math, os
import numpy as np
from . import fingerprint as FP

SAME = 2e-9
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, '..', 'runs', 'clusters')


def radius(n):
    return 0.42 / math.sqrt(n)


def _sym(sq, s, k):
    x, y, a = sq[:, 0].copy(), sq[:, 1].copy(), sq[:, 2].copy()
    if k >= 4:
        x, a = s - x, -a
    for _ in range(k % 4):
        x, y, a = s - y, x, a + 90
    return np.stack([x, y, np.mod(a, 90)], 1)


def md(A, sA, B, sB):
    """matching distance = rediscover.distance (vectorised as in analysis/similarity110.py): min over the 8 box symmetries
    of the RMS optimal-assignment cost (centre moves + angle moves at the corner radius); A / B: n x 3 (x, y, deg)"""
    from scipy.optimize import linear_sum_assignment
    A = np.asarray(A, float); B = np.asarray(B, float)
    Af = A.copy(); Af[:, :2] *= sB / sA
    best = math.inf
    for k in range(8):
        Bk = _sym(B, sB, k)
        dx = Af[:, None, 0] - Bk[None, :, 0]; dy = Af[:, None, 1] - Bk[None, :, 1]
        da = np.mod(Af[:, None, 2] - Bk[None, :, 2] + 45, 90) - 45
        C = dx * dx + dy * dy + (0.5 * math.sqrt(2) * np.radians(da)) ** 2
        r, c = linear_sum_assignment(C)
        best = min(best, math.sqrt(C[r, c].mean()))
    return best


class Leaders:
    def __init__(self, n):
        self.n, self.r = n, radius(n)
        self.ids, self.s, self.sq, self.img = [], [], [], []
        self.bank = np.zeros((0, 2, FP.G, FP.G))

    def assign(self, pid, s, sq):
        """index of the leader within r (nearest by md among img-prefiltered candidates), else a new leader"""
        im = FP.image(s, sq)
        if len(self.ids):
            d = FP.nearest(im, self.bank, self.n)
            cand = np.nonzero(d < 2 * self.r)[0]
            if len(cand):
                cand = cand[np.argsort(d[cand])][:8]
                m = [(md(self.sq[j], self.s[j], sq, s), j) for j in cand]
                dm, j = min(m)
                if dm <= self.r:
                    return int(j)
        self.ids.append(pid); self.s.append(s); self.sq.append(np.asarray(sq)); self.img.append(im)
        self.bank = np.concatenate([self.bank, im[None]])
        return len(self.ids) - 1


def labels(n, rows, tag='live'):
    """rows: [(id, s, coords)] at n.  Returns {id: 'c<n>:<leader store id>'} (same side key -> same label)."""
    if not rows:
        return {}
    rows = sorted(rows, key=lambda r: (r[1], r[0]))
    os.makedirs(CACHE, exist_ok=True)
    cf = os.path.join(CACHE, f'n{n}_{tag}.json')
    C = json.load(open(cf)) if os.path.exists(cf) else None
    L = Leaders(n)
    keys, klab, lab = [], [], {}                      # sorted side keys and their leader index
    if C and 'lsq' in C:                              # leaders carry their own coordinates: rows may be any subset
        for pid, ls, lsq in zip(C['leaders'], C['ls'], C['lsq']):
            q = np.asarray(lsq); L.ids.append(pid); L.s.append(ls); L.sq.append(q); L.img.append(FP.image(ls, q))
        L.bank = np.stack(L.img) if L.img else L.bank
    else:
        C = None
    known = C['labels'] if C else {}
    for pid, s, sq in rows:
        i = bisect.bisect_left(keys, s - SAME)
        if i < len(keys) and abs(keys[i] - s) < SAME:
            lab[pid] = klab[i]; continue
        j = known.get(str(pid))
        if j is None:
            j = L.assign(pid, s, sq)
        keys.insert(i, s); klab.insert(i, j); lab[pid] = j
    known.update({str(k): v for k, v in lab.items()})
    json.dump(dict(n=n, best=min(rows[0][1], C['best'] if C else math.inf), leaders=L.ids, ls=L.s,
                   lsq=[q.tolist() for q in L.sq], labels=known), open(cf + '.tmp', 'w'))
    os.replace(cf + '.tmp', cf)
    return {pid: f'c{n}:{L.ids[j]}' for pid, j in lab.items()}
