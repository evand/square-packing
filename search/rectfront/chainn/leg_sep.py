"""Legibility pass (2026-10-25): the band conditions are separable in the two slopes.
kap(x, -y) = C + S + k(C, S, x) + k(S, C, y),  k(c, s, x) = (s - H/2) x + max(c x - s, 0).
phi_j(c, s; X) = inf_{x >= X} u_j(x) + k(c, s, x);  psi_j(c, s) = sup_{y >= 0} u_j(y) - k(c, s, y).
Slacks: Sa = phi0(C,S) + C + S - lam - psi1(S,C); Sb = phi1(C,S) + C+S-lam - psi0(S,C);
C0 = phi0(C,S;1/2) + phi0(S,C;1/2) + C+S-lam - 2(a+b); C1 likewise with u1."""
import numpy as np, math
r = math.sqrt(2)
def pots(H):
    lam = 2*r - H; a = 1 - lam; b = 1 - H/2
    Lz = lambda x: a + b*(2 - x)
    T = lambda x: a + b*6/7 + (H/16 - 0.1)*(x - 8/7)
    L = lambda x: np.maximum(Lz(x), T(x))
    u0 = lambda x: np.minimum(2*a + b*x, L(x))
    u1 = lambda x: np.minimum(a + b*x, a + L(x))
    return lam, a, b, u0, u1
X = np.concatenate([np.linspace(0, 4, 400001), np.linspace(4, 60, 20001)])
def k(c, s, x, H): return (s - H/2)*x + np.maximum(c*x - s, 0)
def slacks(H, th):
    lam, a, b, u0, u1 = pots(H)
    C, S = math.cos(th), math.sin(th)
    X2 = X[X >= 0.5]
    phi0 = (u0(X) + k(C, S, X, H)).min(); phi1 = (u1(X) + k(C, S, X, H)).min()
    phi0h = (u0(X2) + k(C, S, X2, H)).min(); phi0h_ = (u0(X2) + k(S, C, X2, H)).min()
    phi1h = (u1(X2) + k(C, S, X2, H)).min(); phi1h_ = (u1(X2) + k(S, C, X2, H)).min()
    psi0 = (u0(X) - k(S, C, X, H)).max(); psi1 = (u1(X) - k(S, C, X, H)).max()
    base = C + S - lam
    return dict(Sa=phi0 + base - psi1, Sb=phi1 + base - psi0,
                C0=phi0h + phi0h_ + base - 2*(a + b), C1=phi1h + phi1h_ + base - 2*(a + b))
if __name__ == '__main__':
    H9 = (14*r - 5)/8
    for H in [H9, 1.855, 1.86, 1.87]:
        ths = np.linspace(0, math.pi/2, 181)
        res = {key: [] for key in ('Sa', 'Sb', 'C0', 'C1')}
        for th in ths:
            for kk, v in slacks(H, th).items(): res[kk].append(v)
        print(f"H={H:.5f}", '  '.join(f"{kk}: min {min(v):+.5f} at {math.degrees(ths[int(np.argmin(v))]):.1f}" for kk, v in res.items()))
