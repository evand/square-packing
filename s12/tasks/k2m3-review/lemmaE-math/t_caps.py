"""E' and E'' identities, derived independently and checked against exact polygon clipping."""
import random, sys
from fractions import Fraction as F
from geom import *

a, b = F(9, 5), F(26, 5)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)


def g_formula(d, s, c):
    if d <= 0: return F(0)
    if d <= s: return d * d / (2 * s * c)
    if d <= c: return (d - s / 2) / c
    if d <= s + c: return 1 - (s + c - d) ** 2 / (2 * s * c)
    return F(1)


def ghat(d, s, c):
    if d <= 0: return F(0)
    if d <= s: return d * d / (2 * s * c)
    return (d - s / 2) / c


def rand_u():
    k = rng.random()
    if k < 0.2: return F(rng.randint(1, 1000), 10 ** 6)                       # tiny tilt
    if k < 0.4: return F(41421356, 10 ** 8) - F(rng.randint(0, 1000), 10 ** 8)  # near 45 deg
    return F(rng.randint(1, 4142), 10 ** 4)


bad = 0
# ---- E': g, ghat, parabola
for it in range(4000):
    u = rand_u(); c, s = trig(u)
    cx = F(rng.randint(0, 4000), 1000); cy = F(5)
    P = square(cx, cy, u)
    d = a - min(p[0] for p in P)
    left = area(clip(P, F(1), F(0), -a))                    # x <= a
    g = g_formula(d, s, c)
    if left != g: bad += 1; print("g mismatch", u, cx, left, g)
    gh = ghat(d, s, c)
    par = d * d / (2 * s * c)
    if gh < g: bad += 1; print("ghat < g", u, d)
    if par < gh: bad += 1; print("parabola < ghat", u, d)
    # tan mode: area(Q ∩ {x >= a}) >= tangent at d* <= 1 when d >= c
    if d >= c:
        right = 1 - left
        for ds in (F(1), c, (c + 1) / 2, F(rng.randint(0, 1000), 1000)):
            if ds > s + c: continue
            k = s + c - ds
            T = k * k / (2 * s * c) - k / (s * c) * (d - ds)
            if T > right: bad += 1; print("tan above", u, d, ds)
print("E' done, bad =", bad)

# ---- E'': corner3 / xcut identities
nc3 = [0, 0]; nxc = [0, 0]
for it in range(20000):
    u = rand_u(); c, s = trig(u); C, S, N = 1 - u * u, 2 * u, 1 + u * u
    cx = a + F(rng.randint(-900, 900), 1000); cy = a + F(rng.randint(-900, 900), 1000)
    P = square(cx, cy, u)
    BL, BR, TR, TL = P
    # check vertex labelling for theta in (0,45]
    assert BL[1] == min(p[1] for p in P) and BR[0] == max(p[0] for p in P) and TL[0] == min(p[0] for p in P)
    ALL = area(clip(clip(P, F(1), F(0), -a), F(0), F(1), -a))
    gx = area(clip(P, F(1), F(0), -a)); gy = area(clip(P, F(0), F(1), -a))
    # corner3
    X_aa = ((a - cx) * C + (a - cy) * S) / N
    if BL[0] <= a <= BR[0] and TL[1] >= a and X_aa <= HALF:
        al = (a - BL[0]) * N; be = (a - BL[1]) * N; z = be * C - al * S
        quad = (2 * C * al * be + S * (be * be - al * al)) / (2 * C * N * N)
        if z <= 0:
            nc3[0] += 1
            if ALL != gy: bad += 1; print("corner3 z<=0", u, cx, cy, ALL, gy)
        if z >= 0:
            nc3[1] += 1
            if ALL != quad: bad += 1; print("corner3 z>=0", u, cx, cy, ALL, quad)
    # xcut
    if BL[0] >= a and TL[1] >= a:
        alp = (a - TL[0]) * N; bpp = (a - TL[1]) * N; zp = alp * C + bpp * S
        A = max(F(0), zp) ** 2 / (2 * S * C * N * N)
        if ALL != A: bad += 1; print("xcut A_LL", u, cx, cy, ALL, A)
        nxc[0] += 1
        if zp >= 0:
            nxc[1] += 1
            quad = (2 * C * alp * bpp + S * (bpp * bpp - alp * alp)) / (2 * C * N * N)
            if -gx + ALL != quad: bad += 1; print("xcut quad", u, cx, cy)
print("corner3 z<=0/z>=0 tested", nc3, "xcut", nxc, "bad =", bad)
