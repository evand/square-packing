#!/usr/bin/env python3
"""Final5 leaves as polynomials in half angles u = tan(t1/2), v = tan(t2/2) (2026-10-21).
L_bound(u, v, H) - target, for a case (eps1, eps2, x, br1, br2), with mu's sign and the sign of the kink
term fixed.  Prints degree and factorisation of the numerator at H in {1, 2} (and h0 symbolic)."""
import sympy as sp
u, v, H, r = sp.symbols('u v H r', positive=True)   # r = sqrt 2
C1, S1 = (1 - u**2) / (1 + u**2), 2 * u / (1 + u**2)
C2, S2 = (1 - v**2) / (1 + v**2), 2 * v / (1 + v**2)

def bnd(C, S, br):
    return 2 * S - C * (H - 1) if br == 'U' else (C * (S - C) + 1 + C + S - H) / C

def leaf(e1, e2, xi, br1, br2, mu_pos, kink_pos):
    a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
    v1 = C1 if e1 else -C1; v2 = C2 if e2 else -C2
    x = [(1, 0), a1, a2, (0, 1)][xi]
    cr = lambda p, q: p[0] * q[1] - p[1] * q[0]
    dt = lambda p, q: p[0] * q[0] + p[1] * q[1]
    mu = v2 / S2 - v1 / S1
    y = (x[0] + x[1]) / 2 if mu_pos else H - (x[0] + x[1]) / 2
    def ab(z, pos): return z if pos else -z
    k1 = cr(a1, x); k2 = cr(a2, x)
    # sign choices: kink_pos applies to the one nonzero kink when x is aligned; upright kinks have known sign
    h1 = (ab(k1, kink_pos) if xi == 2 else (-k1 if xi == 0 else (k1 if xi == 3 else 0))) + dt(a1, x)
    h2 = (ab(k2, kink_pos) if xi == 1 else (-k2 if xi == 0 else (k2 if xi == 3 else 0))) + dt(a2, x)
    m = mu * y + h1 / (2 * S1) + h2 / (2 * S2)
    s1 = bnd(C1, S1, br1) + (C1 * H if e1 else 0)
    X2 = bnd(C2, S2, br2) + (0 if e2 else C2 * H)
    return sp.together(s1 / S1 + m + X2 / S2)

if __name__ == '__main__':
    for case in [(True, False, 1, 'P', 'P', True, True), (True, False, 1, 'P', 'P', True, False),
                 (False, False, 1, 'U', 'U', True, True), (False, True, 1, 'U', 'P', True, True)]:
        Lb = leaf(*case)
        for Hv, tgt in [(1, 5), (2, 6 * r - 4)]:
            num, den = sp.fraction(sp.together(Lb.subs(H, Hv) - tgt))
            num = sp.expand(num)
            print(case, 'H=', Hv, 'deg(u,v)=', sp.Poly(num, u, v).total_degree(), ' den =', sp.factor(den))
            print('   num =', sp.factor(num))
