#!/usr/bin/env python3
"""Parity fixed point with routes (2026-10-24).  d[P][m] = least excess of a prefix (any length) from the wall to a
line of slope m with index parity P; predecessors and tilts kept, so each slope's optimal route can be read off."""
import math, sys
import numpy as np
from kwalk import Kmat

R2 = math.sqrt(2)

def fixpoint(H, phis):
    ms = np.tan(np.radians(phis)); ms[np.abs(ms) < 1e-12] = 0
    K, ts, arg = Kmat(H, ms)
    lam = 2 * R2 - H; a = 1 - lam; beta = 1 - H / 2; g = a + beta
    c = K - lam; z = int(np.argmin(np.abs(ms))); M = len(ms)
    d = [np.full(M, np.inf), c[z].copy()]
    pred = [np.full(M, -1), np.full(M, z)]           # predecessor slope index (parity 1 - P); -1 = none, wall = z at P=1 start
    start = [np.zeros(M, bool), np.ones(M, bool)]    # True: the route is the single wall step
    for _ in range(400):
        ch = False
        for P in (0, 1):
            v = d[1 - P][:, None] + c
            i = v.argmin(0); nv = v[i, np.arange(M)]
            better = nv < d[P] - 1e-15
            if better.any():
                ch = True
                d[P] = np.where(better, nv, d[P]); pred[P] = np.where(better, i, pred[P])
                start[P] = np.where(better, False, start[P])
        if not ch: break
    return dict(ms=ms, phis=phis, c=c, ts=ts, arg=arg, d=d, pred=pred, start=start, z=z, lam=lam, a=a, beta=beta, g=g)

def route(F, P, i):
    """List of (slope deg, tilt deg) from the wall to line i of parity P."""
    out = []; Q = P
    while True:
        if F['start'][Q][i]:
            out.append((F['phis'][i], math.degrees(F['ts'][F['arg'][F['z'], i]]))); break
        j = F['pred'][Q][i]
        out.append((F['phis'][i], math.degrees(F['ts'][F['arg'][j, i]])))
        i = j; Q = 1 - Q
        if len(out) > 30: out.append(('...', 0)); break
    return out[::-1]

if __name__ == '__main__':
    H = float(sys.argv[1]) if len(sys.argv) > 1 else 1.846
    F = fixpoint(H, np.arange(-85, 85 + 1e-9, 0.25))
    a, beta = F['a'], F['beta']; am = np.abs(F['ms'])
    uO = a + beta * am; uE = np.minimum(2 * a + beta * am, a + beta * (2 - am))
    print(f"H={H} lam={F['lam']:.5f} a={a:.5f} beta={beta:.5f}")
    for P, u in ((1, uO), (0, uE)):
        print(f"parity {P}:")
        for ph in range(0, 86, 5):
            i = int(np.argmin(np.abs(F['phis'] - ph)))
            r = route(F, P, i)
            print(f"  phi {ph:+3d}: d {F['d'][P][i]:.5f}  (Lean u {u[i]:.5f}, d-u {F['d'][P][i]-u[i]:+.4f})  len {len(r)}: "
                  + " -> ".join(f"[{s:+.2f} t{t:.1f}]" if s != '...' else '...' for s, t in r))
