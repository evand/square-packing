"""Closed-form one-line bounds (legibility pass 2026-10-25).
For a branch l(x) = p + q x and a square (c, s) (k(c,s,x) = (s - H/2) x + (c x - s)^+, kink at t = s/c):
  I(l)  = inf_{x>=X} l(x) + k(x)   (convex; min at x = X or the aligned x = t)
  Sup(l) = sup_{y>=0} l(y) - k(y)  (concave; max at y = 0 or aligned y = t)
u0 = min(A0, max(Lz, T)) >= min(I A0, max(I Lz, I T));  u0 <= min(Sup A0, max(Sup Lz, Sup T)).
u1 = min(A1, a + max(Lz, T)) likewise."""
import math, numpy as np
r = math.sqrt(2); H9 = (14*r - 5)/8
def br(H):
    lam = 2*r - H; a = 1 - lam; b = 1 - H/2; T1 = H/16 - 0.1; T0 = a + b*6/7 - T1*8/7
    return lam, a, b, {'A0': (2*a, b), 'Lz': (a + 2*b, -b), 'T': (T0, T1), 'A1': (a, b)}
def f_at(p, q, c, s, H, x): return p + q*x + (s - H/2)*x + max(c*x - s, 0)
def g_at(p, q, c, s, H, x): return p + q*x - (s - H/2)*x - max(c*x - s, 0)
def I(pq, c, s, H, X=0.0):
    p, q = pq; vals = [f_at(p, q, c, s, H, X)]
    if c > 1e-15 and s/c > X: vals.append(f_at(p, q, c, s, H, s/c))
    assert q + s - H/2 + c >= -1e-12
    return min(vals)
def Sup(pq, c, s, H):
    p, q = pq; vals = [g_at(p, q, c, s, H, 0.0)]
    if c > 1e-15: vals.append(g_at(p, q, c, s, H, s/c))
    assert q - s + H/2 - c <= 1e-12, (pq, c, s)
    return max(vals)
def bounds(H, c, s, X=0.0):
    lam, a, b, B = br(H)
    phi0 = min(I(B['A0'], c, s, H, X), max(I(B['Lz'], c, s, H, X), I(B['T'], c, s, H, X)))
    phi1 = min(I(B['A1'], c, s, H, X), a + max(I(B['Lz'], c, s, H, X), I(B['T'], c, s, H, X)))
    psi0 = min(Sup(B['A0'], c, s, H), max(Sup(B['Lz'], c, s, H), Sup(B['T'], c, s, H)))
    psi1 = min(Sup(B['A1'], c, s, H), a + max(Sup(B['Lz'], c, s, H), Sup(B['T'], c, s, H)))
    return phi0, phi1, psi0, psi1
if __name__ == '__main__':
    for H in np.linspace(H9, 1.87, 5):
        lam, a, b, B = br(H); worst = {}
        for th in np.linspace(0, math.pi/2, 1801):
            C, S = math.cos(th), math.sin(th); base = C + S - lam
            p0, p1, _, _ = bounds(H, C, S); _, _, q0, q1 = bounds(H, S, C)
            h0, h1, _, _ = bounds(H, C, S, 0.5); h0_, h1_, _, _ = bounds(H, S, C, 0.5)
            v = dict(Sa=p0 + base - q1, Sb=p1 + base - q0, C0=h0 + h0_ + base - 2*(a + b), C1=h1 + h1_ + base - 2*(a + b))
            for kk, vv in v.items():
                if kk not in worst or vv < worst[kk][0]: worst[kk] = (vv, math.degrees(th))
        print(f"H={H:.5f} " + '  '.join(f"{kk} {vv:+.2e}@{d:.1f}" for kk, (vv, d) in worst.items()))

def SupI(pq, c, s, H, lo, hi):
    """sup over [lo, hi] of l - k: concave, so at clip(t, lo, hi)."""
    p, q = pq; t = s/c if c > 1e-15 else 1e12
    y = min(max(t, lo), hi)
    if hi == math.inf and t > 1e11: y = lo if q - s + H/2 <= 0 else 1e9
    return g_at(p, q, c, s, H, y)
def psi_clip(H, c, s, Y0=None, Y1=None):
    lam, a, b, B = br(H); T0, T1 = B['T']
    x2 = (2*b - a)/(2*b) if Y0 is None else Y0
    y1 = T0/(b - T1) if Y1 is None else Y1
    psi0 = max(SupI(B['A0'], c, s, H, 0, x2), SupI(B['Lz'], c, s, H, x2, 8/7), SupI(B['T'], c, s, H, 8/7, math.inf))
    psi1 = max(SupI(B['A1'], c, s, H, 0, y1), a + SupI(B['T'], c, s, H, y1, math.inf))
    return psi0, psi1
def check_clip(Y0=None, Y1=None, label=''):
    for H in np.linspace(H9, 1.87, 5):
        lam, a, b, B = br(H); worst = {}
        for th in np.linspace(0, math.pi/2, 1801):
            C, S = math.cos(th), math.sin(th); base = C + S - lam
            p0, p1, _, _ = bounds(H, C, S); q0, q1 = psi_clip(H, S, C, Y0, Y1)
            for kk, vv in dict(Sa=p0 + base - q1, Sb=p1 + base - q0).items():
                if kk not in worst or vv < worst[kk][0]: worst[kk] = (vv, math.degrees(th))
        print(label, f"H={H:.5f} " + '  '.join(f"{kk} {vv:+.2e}@{d:.1f}" for kk, (vv, d) in worst.items()))
if __name__ == '__main__':
    import sys
    check_clip(label='crossing')
