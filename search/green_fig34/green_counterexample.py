#!/usr/bin/env python3
"""Exact refutation of the reconstructed Green set (H1) at its claimed side t = (40 sqrt2 + 19)/17.

Exhibits one closed unit square Q, with rational centre and rational tan(theta/2), that lies in
[0, t]^2 and contains none of the 16 points.  All arithmetic is exact in Q(sqrt2) (sympy); a
sign test of a + b sqrt2 is done exactly via sympy's comparison of algebraic numbers.

The pose comes from the float scan (scan_sets.py): theta ~ 30.36 deg, centre ~ (2.7063, 2.2130),
in the middle band, inside the triangle (a+1, b+V), (a+2, b+V), (a+1+e, t-b-V) whose side from
(a+1+e, t-b-V) to (a+2, b+V) is sqrt((1-e)^2 + V^2) = 1.0110 > 1 (Lemma 3 of DS7 does not apply).
"""
import sympy as sp

r2 = sp.sqrt(2)
t = (40 * r2 + 19) / 17
s = sp.Integer(1)
V = (2 * r2 + 12) / 17
a = (16 * r2 - 6) / 17
b = r2 - sp.Rational(1, 2)
e = t - 2 * a - 2 * s
assert sp.simplify(e - (8 * r2 - 3) / 17) == 0
assert sp.simplify(e**2 + V**2 - 1) == 0           # the diagonal Green's bound makes tight
A = [a, a + e, a + e + s, t - a]
B = [a, a + s, a + 2 * s, t - a]
P = []
for y, row in ((b, A), (b + V, B), (t - b - V, A), (t - b, B)):
    P += [(x, y) for x in row]

u = sp.Rational(27133, 100000)                     # tan(theta/2), theta ~ 30.36 deg
cx, cy = sp.Rational(27063, 10000), sp.Rational(22130, 10000)
den = 1 + u * u
co, si = (1 - u * u) / den, 2 * u / den            # exact rational cos, sin
w = co + si                                         # theta in (0, 90): |cos| + |sin|
half = sp.Rational(1, 2)

ok = True
# containment of the closed square in [0,t]^2: cx, cy in [w/2, t - w/2]
for val in (cx - w / 2, t - w / 2 - cx, cy - w / 2, t - w / 2 - cy):
    v = sp.nsimplify(val)
    ok &= bool(v > 0)
    print(f"  wall slack {float(v):+.6f}")
worst = None
for (px, py) in P:
    dx, dy = px - cx, py - cy
    p1 = sp.expand(dx * co + dy * si)               # coordinate along e1
    p2 = sp.expand(-dx * si + dy * co)              # coordinate along e2
    m = max(abs(float(p1)), abs(float(p2)))
    outside = bool(sp.Abs(p1) > half) or bool(sp.Abs(p2) > half)
    ok &= outside
    worst = m if worst is None else min(worst, m)
    print(f"  point ({float(px):.4f}, {float(py):.4f}): |frame coords| max = {m:.6f}  outside: {outside}")
print(f"min over points of max(|x'|, |y'|) = {worst:.6f}  (> 1/2 means the closed unit square misses it)")
print("REFUTED: a closed unit square in [0,t]^2 at side t = (40 sqrt2+19)/17 contains none of the 16 points"
      if ok else "not refuted by this pose")
