#!/usr/bin/env python3
"""Band surgery: delete squares with centre in a horizontal band [cy, cy+w) and a vertical band [cx, cx+w),
shift the rest to close the gaps; side s -> s - w.  Reports surviving count and overlap (pen) stats."""
import sys, itertools
from gen import load_site, pen, wall_viol, write


def cut(sq, cx, cy, w=1.0):
    out = []
    for x, y, a in sq:
        if cy <= y < cy + w or cx <= x < cx + w:
            continue
        out.append((x - w if x >= cx + w else x, y - w if y >= cy + w else y, a))
    return out


def badness(sq, s):
    worst = 0.0; tot = 0.0
    for i in range(len(sq)):
        worst = max(worst, wall_viol(sq[i], s)); 
        for j in range(i + 1, len(sq)):
            p = pen(sq[i], sq[j])
            if p > 0:
                tot += p * p; worst = max(worst, p)
    return worst, tot


if __name__ == '__main__':
    src = int(sys.argv[1]); target = int(sys.argv[2]); w = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    s, sq = load_site(src)
    res = []
    grid = [i * 0.1 for i in range(int((s - w) / 0.1) + 1)]
    for cx in grid:
        for cy in grid:
            c = cut(sq, cx, cy, w)
            if len(c) >= target:
                wv, tot = badness(c, s - w)
                res.append((tot, wv, len(c), round(cx, 2), round(cy, 2)))
    res.sort()
    for r in res[:20]:
        print('E=%.3e worst=%.3f n=%d cx=%.1f cy=%.1f' % r)
