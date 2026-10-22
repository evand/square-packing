"""Min of P / (S1 S2) over the sign region, per distinct leaf (reps.pkl)."""
import pickle, math, numpy as np, sympy as sp
from cases_sym import C1, S1, C2, S2
from polys import r
reps = pickle.load(open('reps.pkl', 'rb'))
ts = np.radians(np.linspace(0.05, 89.95, 500)); T1, T2 = np.meshgrid(ts, ts, indexing='ij')
a = (np.cos(T1), np.sin(T1), np.cos(T2), np.sin(T2), math.sqrt(2))
for i, (key, tag, tgt, h, P, cons) in enumerate(reps):
    f = sp.lambdify((C1, S1, C2, S2, r), P)
    ok = np.ones_like(T1, bool)
    for c in cons: ok &= sp.lambdify((C1, S1, C2, S2, r), c)(*a) * np.ones_like(T1) >= 0
    V = np.where(ok, f(*a) / (np.sin(T1) * np.sin(T2)), np.inf)
    j = np.unravel_index(np.argmin(V), V.shape)
    print(f"{i:3d} {V[j]:+.5f} at ({math.degrees(T1[j]):.1f},{math.degrees(T2[j]):.1f}) {key} {tgt[:5]} H={sp.N(h,5)}")
