#!/usr/bin/env python3
"""Convert wand125's point_n61_L8/cover.txt (commit f8846cec) to our plain certificate format
(certificates/FORMAT.md) and check it exactly.  Written for this replay; shares no code with
their check_cover.py (which was read, not executed).

    python3 search/s61_wand125/convert_check.py THEIR_cover.txt OUT.txt

Their file is whitespace-separated integers: s_num s_den D W m, then m lines X Y w -- the same
layout as FORMAT.md.  We parse every token as an integer (refusing anything else), re-emit it in
canonical form (one header value per line, single spaces, trailing newline), and check, all in
integers / Fractions:
  * header positive, s = 8 exactly, s_den | s_num*D, exactly m point lines, no trailing data;
  * every point in [0,8]^2, every weight >= 0, no duplicate coordinates;
  * D4 invariance of the weighted multiset (generators x->8-x, y->8-y, x<->y);
  * total = sum w / W < 61, printed as a reduced fraction.
Also reports how much weight sits on the lines x or y in (1/2)Z and Z (where tight poses live)."""
import sys, hashlib, re
from fractions import Fraction as F
from collections import Counter

src, dst = sys.argv[1], sys.argv[2]
raw = open(src, 'rb').read()
toks = raw.split()
assert all(re.fullmatch(rb'-?[0-9]+', t) for t in toks), 'non-integer token'
v = [int(t) for t in toks]
sn, sd, D, W, m = v[:5]
assert min(sn, sd, D, W, m) > 0
assert F(sn, sd) == 8, (sn, sd)
assert (sn * D) % sd == 0
assert len(v) == 5 + 3 * m, 'point count / trailing data'
S = sn * D // sd                       # side in units of 1/D
pts = [tuple(v[5 + 3 * i: 8 + 3 * i]) for i in range(m)]
mass = Counter()
for x, y, w in pts:
    assert 0 <= x <= S and 0 <= y <= S, (x, y)
    assert w >= 0, w
    mass[(x, y)] += w
dup = m - len(mass)
d4 = all(mass.get(img, 0) == w for (x, y), w in mass.items()
         for img in ((S - x, y), (x, S - y), (y, x)))
tot = F(sum(w for _, _, w in pts), W)
with open(dst, 'w') as f:
    f.write(f'{sn} {sd}\n{D}\n{W}\n{m}\n')
    for x, y, w in pts:
        f.write(f'{x} {y} {w}\n')
out = open(dst, 'rb').read()
# round trip: our file parses to the same integers
assert [int(t) for t in out.split()] == v
zero = sum(1 for (_, _, w) in pts if w == 0)
half = D // 2
on_half = sum(w for x, y, w in pts if x % half == 0 or y % half == 0)
on_int = sum(w for x, y, w in pts if x % D == 0 or y % D == 0)
on_bdry = sum(1 for x, y, w in pts if x in (0, S) or y in (0, S))
xs = sorted({x for x, _, _ in pts})
print(f'source  {src} sha256 {hashlib.sha256(raw).hexdigest()}')
print(f'ours    {dst} sha256 {hashlib.sha256(out).hexdigest()} (byte-identical to source: {raw == out})')
print(f'header  s = {sn}/{sd}, D = {D}, W = {W}, m = {m}; distinct coordinates {len(mass)} (duplicates {dup})')
print(f'points  all in [0,8]^2; weights >= 0 (zero weights: {zero}); min w {min(w for *_, w in pts)}, max w {max(w for *_, w in pts)}')
print(f'D4      invariant: {d4}')
print(f'total   {tot} = {tot.numerator}/{tot.denominator}  ~ {float(tot):.14f};  61 - total = {61 - tot} ~ {float(61 - tot):.6e};  < 61: {tot < 61}')
print(f'lines   weight on x or y in (1/2)Z: {float(F(on_half, W)):.6f};  on x or y in Z: {float(F(on_int, W)):.6f};  points on the container boundary: {on_bdry}')
print(f'coords  {len(xs)} distinct x values; min x {F(xs[0], D)}, max x {F(xs[-1], D)}')
assert d4 and tot < 61 and dup == 0
print('CONVERT+CHECK OK')
