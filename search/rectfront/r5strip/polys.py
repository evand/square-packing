"""Final5 cases as polynomials (2026-10-21).  For each case/sub-case: Lb (lower bound on L) times
D = S1 S2 (C1 if b's branch P) (C2 if d's branch P), minus target times D; linear in H.
Sign constraints: smu * mu >= 0, s1 * cr1 >= 0, s2 * cr2 >= 0 (as polynomials, times S1 S2 for mu)."""
import sympy as sp
from cases_sym import cases, C1, S1, C2, S2, H
r = sp.Symbol('r', positive=True)  # sqrt 2

def poly_case(key, mu, c1, c2, Lb):
    D = S1 * S2 * (C1 if key[3] == 'P' else 1) * (C2 if key[4] == 'P' else 1)
    N = sp.expand(sp.cancel(Lb * D))
    cons = [sp.expand(key[5] * mu * S1 * S2), sp.expand(key[6] * c1), sp.expand(key[7] * c2)]
    return D, N, [c for c in cons if c != 0]

if __name__ == '__main__':
    for key, mu, c1, c2, Lb in cases():
        if key == (0, 1, 'a1', 'P', 'P', 1, 1, 1):
            D, N, cons = poly_case(key, mu, c1, c2, Lb)
            print('D =', D); print('N =', N); print('cons >= 0:', cons)
            Z = sp.expand(N - (4 + 6 * r - 4 * H) * D)
            print('N - zig*D =', sp.collect(Z, H))
