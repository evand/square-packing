"""Split functions (legibility pass 2026-10-25): find alpha(C, S) in a small basis with
   D(C,S) := inf_x [u_f(x) + k(C,S,x)] >= alpha(C,S) >= E(S,C) - base,   E(c,s) = sup_y [u_g(y) - k(c,s,y)]
(Sa: f=u0, g=u1; Sb: f=u1, g=u0), or for the cheap conditions alpha(C,S) + alpha(S,C) = 2g - base, D' >= alpha.
LP maximising the worst margin."""
import math, numpy as np, sys
from scipy.optimize import linprog
from leg_sep import pots, k, r
X = np.concatenate([np.linspace(0, 4, 40001), np.linspace(4, 60, 2001)])
def basis(C, S, deg):
    out = [1.0, C, S]
    if deg >= 2: out += [C*S, C*C]
    if deg >= 3: out += [C**3, C*C*S]
    return np.array(out)
def run(H, cond, deg, nth=361):
    lam, a, b, u0, u1 = pots(H); g = a + b
    uf, ug = {'Sa': (u0, u1), 'Sb': (u1, u0), 'C0': (u0, u0), 'C1': (u1, u1)}[cond]
    A, B = [], []
    for th in np.linspace(0, math.pi/2, nth):
        C, S = math.cos(th), math.sin(th); base = C + S - lam; phi = basis(C, S, deg)
        if cond in ('Sa', 'Sb'):
            D = (uf(X) + k(C, S, X, H)).min(); E = (ug(X) - k(S, C, X, H)).max()
            # alpha + m <= D ;  E - base + m <= alpha
            A.append(np.r_[phi, 1]); B.append(D)
            A.append(np.r_[-phi, 1]); B.append(base - E)
        else:
            X2 = X[X >= .5]; D = (uf(X2) + k(C, S, X2, H)).min()
            # alpha(C,S) = (2g - base)/2 + beta(C,S),  beta antisymmetric: beta = sum c_i (phi(C,S) - phi(S,C))
            ph = phi - basis(S, C, deg)
            A.append(np.r_[ph, 1]); B.append(D - (2*g - base)/2)
    n = len(A[0]); res = linprog(np.r_[np.zeros(n - 1), -1], A_ub=np.array(A), b_ub=np.array(B),
                                  bounds=[(None, None)]*(n - 1) + [(None, 1)], method='highs')
    return res.x
if __name__ == '__main__':
    H9 = (14*r - 5)/8
    for cond in ('Sa', 'Sb', 'C0', 'C1'):
        for deg in (1, 2):
            out = []
            for H in (H9, 1.86, 1.87):
                x = run(H, cond, deg); out.append(f"H={H:.4f} margin={x[-1]:+.2e}")
            print(cond, 'deg', deg, ' | '.join(out))

def basis2(C, S, deg):
    f1 = S*(1 - S)/C if C > 1e-12 else 0.0; f2 = C*(1 - C)/S if S > 1e-12 else 0.0
    out = [1.0, C, S, f1, f2]
    if deg >= 2: out += [C*S, C*C]
    return np.array(out)
if __name__ == '__main__' and 'b2' in sys.argv:
    basis = basis2
    H9 = (14*r - 5)/8
    for cond in ('Sa', 'Sb', 'C0', 'C1'):
        for deg in (1, 2):
            out = []
            for H in (H9, 1.86, 1.87):
                x = run(H, cond, deg); out.append(f"H={H:.4f} margin={x[-1]:+.2e}")
            print('b2', cond, 'deg', deg, ' | '.join(out))
