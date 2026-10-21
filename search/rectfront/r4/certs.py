#!/usr/bin/env python3
"""Exact certificates for the one-square lemma (2026-10-21): each case's inequality as
  target - sum(c_i * fact_i) = m1 (P^2+Q^2-1) + m2 (C^2+S^2-1)   (identity)."""
import sympy as sp
P, Q, C, S, H, x, y, lo, sg = sp.symbols('P Q C S H x y lo sg')
k = (P + Q) / 2; A = x - 1 - k
gm = S * Q + C * P
I = [P**2 + Q**2 - 1, C**2 + S**2 - 1]

def show(name, F):
    F = sp.expand(F)
    q, r = sp.reduced(F, I, P, Q, C, S) if F != 0 else ([0, 0], 0)
    print(name, ' remainder:', sp.simplify(r))
    q = list(q) + [0] * (2 - len(q))
    print('   m1 =', sp.factor(q[0]), '\n   m2 =', sp.factor(q[1]))

Up = 2 * S - C * (H - 1)
Psic = C * (S - C) + 1 + C + S - H
# U case, sigma <= 0: |sigma| = -sigma, sigma = S P - C Q
T_neg = S * x - C * y + (-(S * P - C * Q) + gm) / 2
d = lo - (y + (P - Q) / 2)
rep = (P + Q) * (S + C) + (-(S * P - C * Q)) + gm - 2 * (S + C)      # upright_repl slack (|sigma| = -sigma)
l3 = Q * A + P * d
F = Q * (T_neg - Up) - (S * l3 + Q * C * (H - 1 - lo) + d * (-(S * P - C * Q)) + Q * rep / 2 + Q * C * (1 - Q))
show('U-', F)
# U case, sigma > 0
T_pos = S * x - C * y + ((S * P - C * Q) + gm) / 2
s = S * P - C * Q
F = Q * (T_pos - Up) - (S * l3 + s * (y - k) + S * P * (H - 1 - lo) + s * (2 - H) + (1 - Q) * (C * Q + S * (1 - P)))
show('U+', F)
# D case: P T >= P(S-C) + P sigma+ + gamma (1+P+Q-H)
l3D = P * A + Q * (y + (P - Q) / 2 - lo - 1)
for nm, T, sp_ in [('D+', T_pos, s), ('D-', T_neg, 0)]:
    tau = (P * (S - C) + P * sp_ + gm * (1 + P + Q - H))
    F = P * T - tau - (S * l3D + S * Q * lo + gm * (H - k - y))
    show(nm + ' bound', F)
# D1: P S tau = Q Psic + sigma Up + polyD1
polyD1 = -C*P*S - C*Q**2*S + C*Q*S - C*Q + C*S + P*Q*S**2 - P*S**2 - Q*S**2 - Q*S + Q + S**2
tau = P * (S - C) + P * s + gm * (1 + P + Q - H)
show('D1 identity', S * tau - (Q * Psic + s * Up + polyD1))
# D2: C tau - P Psic = polyD2 + (2 - H)(C gamma - P), and C gamma - P = -S sigma
polyD2 = C * gm * (P + Q - 1) - P * (C + S - 1)
tau0 = P * (S - C) + gm * (1 + P + Q - H)
show('D2 identity', C * tau0 - P * Psic - (polyD2 + (2 - H) * (-S * s)))
