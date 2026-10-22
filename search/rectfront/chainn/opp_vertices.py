"""(S) opposite case: the eight per-line vertex inequalities, min over their region and unconditionally."""
import math, numpy as np
r2 = math.sqrt(2)
res = {}
def upd(key, v, where):
    if v < res.get(key, (9,))[0]: res[key] = (v, where)
for H in np.linspace(2 * r2 - 1, 1.99999, 120):
    b = 1 - H / 2; Z = 2 + 2 * r2 - 2 * H; mu = (Z - 1) / (2 * b)
    for t in np.linspace(1e-5, math.pi / 2 - 1e-5, 3000):
        C, S = math.cos(t), math.sin(t); k = (C + S) / 2
        if C + S > H: continue
        g = lambda x: abs(C + x * S) + abs(x * C - S)
        T = lambda x: g(x) / 2 - abs(x) * (1 - k)
        Tp = lambda x: g(x) / 2 - abs(x) * (H - 1 - k)
        tn = S / C; Vy = Tp(-C / S); Vm = Tp(tn)
        w = (round(math.degrees(t), 2), round(H, 4))
        cond_z = tn >= mu          # Z1, Z2 need the kink inside m >= mu
        for key, val, cond in (("Z1 T'(tan)+T'(-cot) >= Z", Vm + Vy - Z, cond_z),
                               ("Z2 T'(tan)+k >= Z", Vm + k - Z, cond_z),
                               ("Z3 T'(mu)+T'(-cot) >= Z", Tp(max(mu, 0)) + Vy - Z, True),
                               ("Z4 T'(mu)+k >= Z", Tp(max(mu, 0)) + k - Z, True),
                               ("G1 T(tan)+T'(-cot) >= 1", T(tn) + Vy - 1, mu > 0 and tn <= mu),
                               ("G2 T(tan)+k >= 1", T(tn) + k - 1, mu > 0 and tn <= mu),
                               ("G3 T(mu)+T'(-cot) >= 1", T(mu) + Vy - 1, mu > 0 and mu <= tn),
                               ("G4 T(mu)+k >= 1", T(mu) + k - 1, mu > 0 and mu <= tn)):
            if cond: upd(key + " [region]", val, w)
            if not key.startswith("Z3") and not key.startswith("Z4"): upd(key + " [unconditional]", val, w)
for k_ in sorted(res): print(f"{k_:45s} {res[k_][0]:+.4e} at (t, H) = {res[k_][1]}")
