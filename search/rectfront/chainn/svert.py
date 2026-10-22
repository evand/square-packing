"""(S) opposite case: exact forms of the per-line vertex inequalities (multiplied out), factored.
Per-line pieces (k = (C+S)/2):  for x in [0, tan t]: T(x) = k + x(S - 1), T'(x) = k + x(S + 1 - H);
for y in [0, cot t]: T'(-y) = k + y(C + 1 - H).  V_m = T'(tan), V_y = T'(-cot), W = T(tan)."""
import sympy as sp
C, S, H, r = sp.symbols('C S H r', positive=True)
k = (C + S) / 2; b = 1 - H / 2; Z = 2 + 2 * r - 2 * H; mu = (Z - 1) / (2 * b)
tn, ct = S / C, C / S
Vm = k + tn * (S + 1 - H)          # T'(tan) via the first piece (continuous at the kink)
Vy = k + ct * (C + 1 - H)
W = k + tn * (S - 1)
Tpmu = k + mu * (S + 1 - H)        # T'(mu) if mu <= tan
Tmu = k + mu * (S - 1)             # T(mu)  if mu <= tan
red = lambda e: sp.factor(sp.simplify(sp.expand(e).subs(S**2, 1 - C**2)))
for name, e in [("Z1 Vm+Vy-Z", Vm + Vy - Z), ("Z2 Vm+k-Z", Vm + k - Z), ("Z3 T'(mu)+Vy-Z", Tpmu + Vy - Z),
                ("Z4 T'(mu)+k-Z", Tpmu + k - Z), ("G1 W+Vy-1", W + Vy - 1), ("G2 W+k-1", W + k - 1),
                ("G3 T(mu)+Vy-1", Tmu + Vy - 1), ("G4 T(mu)+k-1", Tmu + k - 1)]:
    print(name, ":", red(sp.together(e)))
