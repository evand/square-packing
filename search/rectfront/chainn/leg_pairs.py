"""Reduced inequality list for the band conditions (legibility pass, 2026-10-25).
Each condition is f(x) + g(y) + base >= 0 with f, g one-variable; branches of the potentials make f a min of
convex pieces (lower bound: all candidate points) and -g ... see leg_sep.py.  Here: every (x-candidate, y-candidate)
pair, its min slack over angle and H, and where."""
import math, numpy as np, itertools
r = math.sqrt(2)
H9 = (14*r - 5)/8
def consts(H):
    lam = 2*r - H; a = 1 - lam; b = 1 - H/2; T1 = H/16 - 0.1; T0 = a + b*6/7 - T1*8/7
    x2 = (2*b - a)/(2*b); y1 = T0/(b - T1)
    return lam, a, b, T0, T1, x2, y1
def k(c, s, x, H): return (s - H/2)*x + max(c*x - s, 0)
# branches: name -> (function of x, interval)
def branches(H, which):
    lam, a, b, T0, T1, x2, y1 = consts(H)
    A0 = lambda x: 2*a + b*x; Lz = lambda x: a + b*(2 - x); T = lambda x: T0 + T1*x
    if which == 'u0+':   # lower bound of u0 = min(A0, max(Lz, T)): pieces A0 on [0,inf), Lz on [0, 8/7], T on [8/7, inf) (max >= each)
        return {'A0': (A0, 0, math.inf), 'Lz': (Lz, 0, 8/7), 'T': (T, 8/7, math.inf)}
    if which == 'u1+':
        return {'A1': (lambda x: a + b*x, 0, math.inf), 'aLz': (lambda x: a + Lz(x), 0, 8/7), 'aT': (lambda x: a + T(x), 8/7, math.inf)}
    if which == 'u0-':   # upper bound: u0 = A0 on [0,x2], Lz on [x2,8/7], T on [8/7,inf)
        return {'A0': (A0, 0, x2), 'Lz': (Lz, x2, 8/7), 'T': (T, 8/7, math.inf)}
    if which == 'u1-':
        return {'A1': (lambda x: a + b*x, 0, y1), 'aT': (lambda x: a + T(x), y1, math.inf)}
def cand_pts(c, s, lo, hi):
    pts = [('lo', lo)]
    if hi < math.inf: pts.append(('hi', hi))
    if c > 1e-14 and lo < s/c < hi: pts.append(('k', s/c))
    if hi == math.inf: pts.append(('inf', 1e7))
    return pts
def fcands(H, c, s, which, X0=0.0):
    out = []
    for nm, (f, lo, hi) in branches(H, which).items():
        lo = max(lo, X0)
        if lo >= hi: continue
        for pn, p in cand_pts(c, s, lo, hi):
            out.append((nm + ':' + pn, f(p) + k(c, s, p, H)))
    return out
def gcands(H, c, s, which):  # sup of u - k  ->  values of u - k at candidates
    out = []
    for nm, (f, lo, hi) in branches(H, which).items():
        for pn, p in cand_pts(c, s, lo, hi):
            out.append((nm + ':' + pn, f(p) - k(c, s, p, H)))
    return out
CONDS = {'Sa': ('u0+', 'u1-'), 'Sb': ('u1+', 'u0-'), 'C0': ('u0+', 'u0+'), 'C1': ('u1+', 'u1+')}
if __name__ == '__main__':
    Hs = np.linspace(H9, 1.87, 7); ths = np.linspace(0, math.pi/2, 721)
    for cond, (wf, wg) in CONDS.items():
        worst = {}
        for H in Hs:
            lam, a, b, *_ = consts(H)
            for th in ths:
                C, S = math.cos(th), math.sin(th); base = C + S - lam
                if cond in ('Sa', 'Sb'):
                    F = fcands(H, C, S, wf); G = gcands(H, S, C, wg)
                    for (nf, vf), (ng, vg) in itertools.product(F, G):
                        if abs(vf) > 1e5 or abs(vg) > 1e5:   # slope at infinity
                            v = 1.0 if (vf > 0 or vg < 0) else -1.0
                        else: v = vf + base - vg
                        key = (nf, ng); w = worst.get(key)
                        if w is None or v < w[0]: worst[key] = (v, H, math.degrees(th))
                else:
                    F = fcands(H, C, S, wf, 0.5); G = fcands(H, S, C, wg, 0.5)
                    for (nf, vf), (ng, vg) in itertools.product(F, G):
                        v = vf + vg + base - 2*(a + b) if abs(vf) < 1e5 and abs(vg) < 1e5 else 1.0
                        key = (nf, ng); w = worst.get(key)
                        if w is None or v < w[0]: worst[key] = (v, H, math.degrees(th))
        print(f"== {cond}: {len(worst)} pairs")
        for key, (v, H, d) in sorted(worst.items(), key=lambda kv: kv[1][0]):
            if v < 0.02: print(f"   {key[0]:10s} {key[1]:10s} min {v:+.6f} at H={H:.4f} th={d:.2f}")
