#!/usr/bin/env python3
"""Exact unit tests of zm_mixed.cond_poly / bern / cond_region against first principles.

(1) G_chi(u; p) = N(u) * (ineq_k(p, c_chi(u), u)) with N(u) = 2(1+u^2)^2 (or 2(1+u^2) when both bounds are 'R'),
    where c_chi(u) is the centre bound of the choice chi and ineq_k is  X - 1/2, -X - 1/2, Y - 1/2, -Y - 1/2.
(2) cond_poly is affine in p (value at p equals the 3-point affine reconstruction).
(3) bern: at random t, p(u0 + h t) equals sum beta_i B_i(t).
(4) cond_region(...) vertices: each satisfies inequality k at the slot's worst corner at sampled u -- checked
    in the vertex test of test_poly.py; here: each per-chi region's vertices have all Bernstein coeffs <= 0.
"""
import sys, os, random
from fractions import Fraction as F
S = os.path.expanduser('~/math/square-packing/public/s12/search')
sys.path.insert(0, S)
import zm_mixed as ZM, zeromargin as zm

rng = random.Random(5)
HALF = F(1, 2)


def rf(lo, hi, den=10 ** 5): return lo + (hi - lo) * F(rng.randrange(den + 1), den)


def bound_val(k, m, u):
    kind, val = k
    c, s = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
    if kind == 'R': return val
    if kind == 'W': return (c + s) / 2
    return m - (c + s) / 2


def ineq(cond, px, py, cx, cy, u):
    c, s = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
    X = (px - cx) * c + (py - cy) * s; Y = -(px - cx) * s + (py - cy) * c
    return (X - HALF, -X - HALF, Y - HALF, -Y - HALF)[cond]


def pev(g, u): return sum(g[i] * u ** i for i in range(len(g)))


n = 0
m = F(7)
for trial in range(3000):
    px, py = rf(-2, 9), rf(-2, 9)
    ks = [('R', rf(0, 7)), ('W', None), ('M', None)]
    xk, yk = rng.choice(ks), rng.choice(ks)
    cond = rng.randrange(4)
    g = ZM.cond_poly(px, py, xk, yk, cond, m)
    both = xk[0] == 'R' and yk[0] == 'R'
    for _ in range(5):
        u = rf(0, 1)
        N = 2 * (1 + u * u) if both else 2 * (1 + u * u) ** 2
        # which centre coordinate each slot bounds: xk bounds cx, yk bounds cy
        cx, cy = bound_val(xk, m, u), bound_val(yk, m, u)
        assert pev(g, u) == N * ineq(cond, px, py, cx, cy, u), (px, py, xk, yk, cond, u)
        n += 1
    # affine
    g0 = ZM.cond_poly(F(0), F(0), xk, yk, cond, m)
    gx = ZM.cond_poly(F(1), F(0), xk, yk, cond, m)
    gy = ZM.cond_poly(F(0), F(1), xk, yk, cond, m)
    assert all(g[i] == g0[i] + px * (gx[i] - g0[i]) + py * (gy[i] - g0[i]) for i in range(5))
    # bern
    u0 = rf(0, F(1, 2)); h = rf(0, F(1, 2))
    be = ZM.bern(list(g), u0, h)
    t = rf(0, 1)
    from math import comb
    val = sum(be[i] * comb(4, i) * t ** i * (1 - t) ** (4 - i) for i in range(5))
    assert val == pev(g, u0 + h * t)
print('cond_poly = N * inequality at the bound centre: %d exact checks OK; affine & bern OK' % n)

# slot semantics: worst corner of each inequality over the centre rectangle for theta in [0, 90deg)
for trial in range(2000):
    px, py = rf(-2, 9), rf(-2, 9)
    x0, y0 = rf(0, 7), rf(0, 7); x1, y1 = x0 + rf(0, 1), y0 + rf(0, 1)
    u = rf(0, F(99, 100))
    slot = {0: (x0, y0), 1: (x1, y1), 2: (x1, y0), 3: (x0, y1)}
    for cond in range(4):
        worst = max(ineq(cond, px, py, cx, cy, u) for cx in (x0, x1) for cy in (y0, y1))
        assert ineq(cond, px, py, *slot[cond], u) == worst
print('slot table (Ax,Ay),(Bx,By),(Bx,Ay),(Ax,By) = worst corners for u in [0,1): OK')
