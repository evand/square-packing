"""Cover C (c0 = 0.9): per centre cell (folded to the lower-left quadrant), the inclusion-minimal caught point
sets over all probes (arrangement vertices, 0.25 deg).  Point names: r = (c0, 1) [2/3], s = (1, 1) [1/3],
m = (1, b/2) [1/3]; ' = mirror x -> a - x; ^ = mirror y -> b - y (m is on the midline: m^ = m)."""
import math, collections, numpy as np, ptcover as pc
from c10 import cover, line_b
from tightC import probes, cell
NAMES = ["r", "r^", "r'", "r'^", "s", "s^", "s'", "s'^", "m", "m'"]
def fold(n, fx, fy):
    base = n[0]; pr = "'" in n; up = "^" in n
    if fx: pr = not pr
    if fy and base != "m": up = not up
    return base + ("'" if pr else "") + ("^" if up else "")
if __name__ == '__main__':
    c0 = 0.9
    sets = collections.defaultdict(collections.Counter)
    for a in [1 + math.sqrt(2), 2.45, 2.5, 2.55, 4 - math.sqrt(2)]:
        b = line_b(a); P, W = cover(a, b, c0)
        for th in np.radians(np.arange(0, 90, 0.25)):
            M, k = probes(P, a, b, th)
            c, s = np.cos(th), np.sin(th)
            dx = P[None, :, 0] - M[:, None, 0]; dy = P[None, :, 1] - M[:, None, 1]
            ins = (np.abs(dx*c + dy*s) <= .5 + 1e-12) & (np.abs(-dx*s + dy*c) <= .5 + 1e-12)
            for r in range(len(M)):
                x, y = M[r]; fx, fy = x > a/2, y > b/2
                st = frozenset(fold(NAMES[t], fx, fy) for t in np.where(ins[r])[0])
                sets[cell(x, y, a, b, c0)][st] += 1
    W = dict(r=2/3, s=1/3, m=1/3)
    for cl, cnt in sorted(sets.items()):
        mins = [s for s in cnt if not any(o < s for o in cnt)]
        print(cl)
        for s in sorted(mins, key=lambda s: (len(s), sorted(s))):
            print(f"   {' + '.join(sorted(s)):26s} weight {sum(W[n[0]] for n in s):.3f}  probes {cnt[s]}")
