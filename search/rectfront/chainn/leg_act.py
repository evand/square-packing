"""Which candidate slope attains phi/psi (leg_sep.py), as a function of the angle."""
import math, numpy as np
from leg_sep import pots, k, r
def cands(H, c, s, lo):
    lam, a, b, u0, u1 = pots(H)
    x2 = (2*b - a)/(2*b)
    T0 = a + b*6/7 - (H/16 - 0.1)*8/7; T1 = H/16 - 0.1; y1 = T0/(b - T1)
    pts = {'lo': lo, 'kink': s/c if c > 1e-12 else 1e9, '8/7': 8/7, 'x2': x2, 'y1': y1, 'inf': 1e6}
    return {n: p for n, p in pts.items() if p >= lo}
def best(f, pts, mode):
    vals = {n: float(f(np.array([p]))[0]) for n, p in pts.items()}
    n = (min if mode == 'min' else max)(vals, key=vals.get)
    return n, vals[n]
H = float(__import__('sys').argv[1]) if len(__import__('sys').argv) > 1 else 1.86
lam, a, b, u0, u1 = pots(H)
prev = None
for deg in np.arange(0, 90.01, 0.5):
    th = math.radians(deg); C, S = math.cos(th), math.sin(th)
    row = []
    row.append(best(lambda x: u0(x) + k(C, S, x, H), cands(H, C, S, 0), 'min')[0])
    row.append(best(lambda x: u1(x) + k(C, S, x, H), cands(H, C, S, 0), 'min')[0])
    row.append(best(lambda x: u0(x) + k(C, S, x, H), cands(H, C, S, .5), 'min')[0])
    row.append(best(lambda x: u1(x) + k(C, S, x, H), cands(H, C, S, .5), 'min')[0])
    row.append(best(lambda x: u0(x) - k(S, C, x, H), cands(H, S, C, 0), 'max')[0])
    row.append(best(lambda x: u1(x) - k(S, C, x, H), cands(H, S, C, 0), 'max')[0])
    if row != prev: print(f"{deg:5.1f}  phi0={row[0]:5s} phi1={row[1]:5s} phi0h={row[2]:5s} phi1h={row[3]:5s} psi0={row[4]:5s} psi1={row[5]:5s}")
    prev = row
