"""Tilted bands of a packing (10-13): connected groups of tilted squares.

Ostashev's README (#575): records near k(k+1) are aligned corner blocks plus bands of tilted squares; relaxing a packing
whose tilted squares form two separate bands merged them into one V and lowered the side (241: bands 47 + 32 -> one V
of 85; de Winter's 211: 48 + 25 -> 76).  `bands(p)` finds the bands so such packings can be selected.

Two tilted squares are in the same band when their centres are within LINK (a touching tilted pair is at most sqrt 2
apart; LINK leaves room for a rattler gap).  Bands smaller than MIN are ignored (isolated tilted rattlers).

  python3 -m pk.bands ID|best:N ...            band sizes per packing
  python3 -m pk.bands --scan 100-324           best packing per n with >= 2 bands or >= 2 angle families
"""
from __future__ import annotations
import sys
import numpy as np

LINK = 1.6
MIN = 5


def tilted(sq, eps=1e-6):
    a = np.asarray(sq)[:, 2] % 90.0
    return (a > eps) & (a < 90.0 - eps)


def bands(sq, link=LINK, mn=MIN):
    """Sizes of the tilted bands (largest first) and the index lists."""
    sq = np.asarray(sq, float)
    idx = np.nonzero(tilted(sq))[0]
    if len(idx) == 0:
        return [], []
    P = sq[idx, :2]
    D = np.linalg.norm(P[:, None] - P[None], axis=2) < link
    seen, comps = np.zeros(len(idx), bool), []
    for i in range(len(idx)):
        if seen[i]:
            continue
        stack, comp = [i], []
        seen[i] = True
        while stack:
            j = stack.pop(); comp.append(j)
            for k in np.nonzero(D[j] & ~seen)[0]:
                seen[k] = True; stack.append(k)
        if len(comp) >= mn:
            comps.append(sorted(idx[comp].tolist()))
    comps.sort(key=len, reverse=True)
    return [len(c) for c in comps], comps


def families(sq, gap=2.5, mn=MIN, eps=1.0):
    """Angle families of the tilted squares: folded angles min(a, 90 - a) (so the two arms of a V fold together), split
    at gaps >= `gap` degrees; families of >= mn squares, sizes largest first.  De Winter's 211 is one connected group
    of 74 tilted squares but two families (18-20 deg: 45, 24-29: 23); Ostashev's merged 211 is one (22-25: 76)."""
    a = np.asarray(sq, float)[:, 2] % 90.0
    a = np.sort(np.minimum(a, 90.0 - a))
    a = a[a > eps]
    if len(a) == 0:
        return []
    cuts = np.nonzero(np.diff(a) >= gap)[0] + 1
    return sorted((len(g) for g in np.split(a, cuts) if len(g) >= mn), reverse=True)


def split(sq):
    """(connected bands, angle families): a packing is a band-merge candidate if either has >= 2 entries."""
    return bands(sq)[0], families(sq)


def _get(st, ref):
    if ref.startswith('best:'):
        r = st.query(n=int(ref[5:]), limit=1)
        return st.get(r[0]['id']) if r else None
    return st.get(int(ref))


if __name__ == '__main__':
    from .store import Store
    st = Store()
    if sys.argv[1] == '--scan':
        lo, hi = map(int, sys.argv[2].split('-'))
        for n in range(lo, hi + 1):
            p = _get(st, f'best:{n}')
            if p is None:
                continue
            b, f = split(p.sq)
            if len(b) >= 2 or len(f) >= 2:
                print(n, f'{p.s:.10f}', 'bands', b, 'families', f)
    else:
        for ref in sys.argv[1:]:
            p = _get(st, ref)
            b, f = split(p.sq)
            print(ref, p.n, f'{p.s:.10f}', 'bands', b, 'families', f)
