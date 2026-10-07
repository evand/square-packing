#!/usr/bin/env python3
"""Make n-square seeds: tilted gadget (transplanted from a record, optionally perturbed) + MIS axis fill.

  mkseed.py --n N --s S --src M [--place ll|scale] [--sym k] [--jitter a] [--seeds K] --out prefix
place ll: keep the source's coordinates (gadget anchored at the lower-left corner) after applying symmetry k;
place scale: scale the source's coordinates by S/s_src (for diagonal bands).
If the fill gives more than N squares, axis squares are dropped (those adjacent to most others last);
if fewer, the missing squares are put at random positions for the optimizer to place.
"""
import argparse, math, random
from gen import load_site, is_axis, mis_fill, write

D4 = [lambda x, y, s: (x, y), lambda x, y, s: (s - x, y), lambda x, y, s: (x, s - y), lambda x, y, s: (s - x, s - y),
      lambda x, y, s: (y, x), lambda x, y, s: (s - y, x), lambda x, y, s: (y, s - x), lambda x, y, s: (s - y, s - x)]
D4a = [1, -1, -1, 1, -1, 1, 1, -1]  # angle sign under each map (rotation part ignored: squares mod 90)


def make(n, s, src, place, sym, jitter, rng, sel=None):
    ss, sq = load_site(src)
    t = [q for q in sq if not is_axis(q)]
    out = []
    for x, y, a in t:
        x, y = D4[sym](x, y, ss)
        a = a * D4a[sym]
        if place == 'scale':
            x, y = x * s / ss, y * s / ss
        out.append((x + rng.gauss(0, jitter), y + rng.gauss(0, jitter), a + rng.gauss(0, jitter * 20)))
    # keep tilted squares whose centre is inside the box (others dropped)
    out = [q for q in out if 0.5 <= q[0] <= s - 0.5 and 0.5 <= q[1] <= s - 0.5]
    if sel is not None:
        out = sel(out)
    ax = mis_fill(s, out, tol=0.02)
    tot = len(out) + len(ax)
    allsq = out + ax
    if tot > n:
        allsq = out + ax[:n - len(out)] if len(out) <= n else out[:n]
    while len(allsq) < n:
        allsq.append((rng.uniform(0.5, s - 0.5), rng.uniform(0.5, s - 0.5), rng.uniform(0, 90)))
    return allsq, tot, len(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int); ap.add_argument('--s', type=float); ap.add_argument('--src', type=int)
    ap.add_argument('--place', default='ll'); ap.add_argument('--sym', type=int, default=0)
    ap.add_argument('--jitter', type=float, default=0.0); ap.add_argument('--seeds', type=int, default=1)
    ap.add_argument('--out')
    a = ap.parse_args()
    for k in range(a.seeds):
        rng = random.Random(k)
        sq, tot, nt = make(a.n, a.s, a.src, a.place, a.sym, a.jitter, rng)
        write(f'{a.out}_{k}.txt', a.s, sq)
        print(f'{a.out}_{k}.txt tilted={nt} fill_total={tot} need={a.n}')
