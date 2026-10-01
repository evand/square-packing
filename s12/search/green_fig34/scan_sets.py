#!/usr/bin/env python3
"""Float scans (exploration, not proof) of the candidate Figure-34 sets.

Prints, for each candidate, smax = the largest square (any angle) avoiding the points in [0,t]^2
at its nominal side t (unavoidable for closed unit squares iff smax <= 1), the worst angle, and
the centre of the worst square with its nearest points.
"""
import sys, os, math, functools
print = functools.partial(print, flush=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from smax import smax, smax_theta
from family import green17, T17
from family19 import friedman19

R2, R65 = math.sqrt(2), math.sqrt(65)

def green_H1():
    """Green's set as reconstructed from the bound: s = 1, e = (8r2-3)/17, V = (2r2+12)/17."""
    return T17, green17(T17, (16 * R2 - 6) / 17, R2 - 0.5, 1.0, (2 * R2 + 12) / 17)

def best_Hstar():
    """The family optimum (all lemma constraints tight): s = 8/r65, e = 4/r65, V = 7/r65."""
    t = 2 * R2 + 13 / R65
    return t, green17(t, R2 - 3.5 / R65, R2 - 4 / R65, 8 / R65, 7 / R65)

def friedman():
    t = 6 * R2 - 4; d = 2 * R2 - 2
    u = (d + math.sqrt(1 - d * d) / 2) / 2
    return t, friedman19(t, u)

def report(name, t, P, nth=360):
    v, th, vals = smax(P, t, nth=nth, nloc=4)
    _, c = smax_theta(P, t, th, return_c=True)
    Pa = np.array(P)
    print(f"{name}: {len(P)} points, t = {t:.9f}, smax = {v:.9f} at theta = {math.degrees(th):.4f} deg, "
          f"centre ({c[0]:.4f}, {c[1]:.4f})")
    for lam in (0.999, 0.9999):
        Q = [(x * lam, y * lam) for x, y in P]
        v2, th2, _ = smax(Q, t * lam, nth=nth // 2, nloc=4)
        print(f"    scaled by {lam}: side {t*lam:.7f}, smax = {v2:.9f} (theta {math.degrees(th2):.3f})")

if __name__ == '__main__':
    which = sys.argv[1:] or ['H1', 'Hstar', 'F19']
    for w in which:
        t, P = {'H1': green_H1, 'Hstar': best_Hstar, 'F19': friedman}[w]()
        report(w, t, P)
