#!/usr/bin/env python3
"""Cantrell-type gadgets: a left-wall channel and a bottom-wall channel of tilted rows, + MIS axis fill.

Left channel (w, L, theta): L rows, each w unit squares along u = (cos t, sin t), row's leftmost point on x = 0,
row i's top point at y = top - i / cos(t) (rows stacked downward at the minimal vertical period).
Bottom channel = transpose (x <-> y, angle -> 90 - angle) of a left channel, so it runs along the bottom wall,
stacked leftward from the right wall when top = s.
"""
import math, sys, itertools
from gen import pen, mis_fill, wall_viol, write, corners


def row(w, t_deg):
    t = math.radians(t_deg)
    return [(j * math.cos(t), j * math.sin(t), t_deg) for j in range(w)]


def left_channel(s, w, L, t_deg, top=None, x0=0.0):
    top = s if top is None else top
    out = []
    per = 1.0 / math.cos(math.radians(t_deg))
    for i in range(L):
        r = row(w, t_deg)
        pts = [p for q in r for p in corners(*q)]
        mx = min(p[0] for p in pts); my = max(p[1] for p in pts)
        dx, dy = x0 - mx, top - i * per - my
        out += [(x + dx, y + dy, a) for x, y, a in r]
    return out


def transpose(sq):
    return [(y, x, (90 - a) % 90) for x, y, a in sq]


def gadget(s, A, B):
    """A, B = (w, L, theta) for the left and bottom channels.  Squares of B that overlap A are dropped."""
    a = left_channel(s, *A)
    b = transpose(left_channel(s, *B))
    a = [q for q in a if wall_viol(q, s) <= 1e-9]
    b = [q for q in b if wall_viol(q, s) <= 1e-9 and all(pen(q, p) <= 1e-9 for p in a)]
    return a + b


def count(s, A, B, tol=1e-6):
    g = gadget(s, A, B)
    ax = mis_fill(s, g, tol)
    return len(g) + len(ax), g, ax


def scan(s, n, tol, ws=((3, 3), (2, 3), (3, 2), (2, 2)), angles=range(16, 41, 2)):
    out = []
    m = int(s)
    for w1, w2 in ws:
        for L1 in range(2, m + 1):
            for L2 in range(2, m + 1):
                for t1 in angles:
                    for t2 in angles:
                        tot, g, ax = count(s, (w1, L1, t1), (w2, L2, t2), tol)
                        out.append((tot, w1, L1, t1, w2, L2, t2))
    out.sort(reverse=True)
    return out


if __name__ == '__main__':
    import argparse, random
    ap = argparse.ArgumentParser()
    ap.add_argument('--s', type=float); ap.add_argument('--n', type=int); ap.add_argument('--tol', type=float, default=1e-6)
    ap.add_argument('--top', type=int, default=20); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    res = scan(a.s, a.n, a.tol)
    for k, r in enumerate(res[:a.top]):
        print(r)
        if a.out:
            tot, w1, L1, t1, w2, L2, t2 = r
            _, g, ax = count(a.s, (w1, L1, t1), (w2, L2, t2), a.tol)
            sq = (g + ax)[:a.n]
            rng = random.Random(k)
            while len(sq) < a.n:
                sq.append((rng.uniform(.5, a.s - .5), rng.uniform(.5, a.s - .5), rng.uniform(0, 90)))
            write(f'{a.out}_{k}.txt', a.s, sq)
