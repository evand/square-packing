"""Per case of Final5 (H in [2sqrt2-1, 2)): min of Lb - 5 and of Lb - zig separately (zig = 4 + 6sqrt2 - 4H)."""
import math, numpy as np, sympy as sp
from cases_sym import cases, C1, S1, C2, S2, H
R2 = math.sqrt(2)
ts = np.radians(np.linspace(0.02, 89.98, 240)); T1, T2 = np.meshgrid(ts, ts, indexing='ij')
rows = []
for key, mu, c1, c2, Lb in cases():
    f = sp.lambdify((C1, S1, C2, S2, H), Lb); fm = sp.lambdify((C1, S1, C2, S2, H), mu)
    g1 = sp.lambdify((C1, S1, C2, S2, H), c1); g2 = sp.lambdify((C1, S1, C2, S2, H), c2)
    m5 = mz = mm = math.inf
    for h in np.linspace(2 * R2 - 1, 1.9999, 21):
        a = (np.cos(T1), np.sin(T1), np.cos(T2), np.sin(T2), h)
        ok = (key[5] * fm(*a) >= 0) & (key[6] * g1(*a) * np.ones_like(T1) >= 0) & (key[7] * g2(*a) * np.ones_like(T1) >= 0)
        if not ok.any(): continue
        L = f(*a) * np.ones_like(T1)
        m5 = min(m5, np.where(ok, L - 5, np.inf).min()); mz = min(mz, np.where(ok, L + 4 * h - 4 - 6 * R2, np.inf).min())
        mm = min(mm, np.where(ok, np.maximum(L - 5, L + 4 * h - 4 - 6 * R2), np.inf).min())
    if m5 < math.inf: rows.append((m5, mz, mm, key))
for m5, mz, mm, k in sorted(rows, key=lambda t: max(t[0], t[1])):
    tag = '5' if m5 >= -1e-9 else ('zig' if mz >= -1e-9 else 'BOTH')
    print(f"{tag:4s} min(L-5)={m5:+.4f} min(L-zig)={mz:+.4f} min(max)={mm:+.4f} e=({k[0]},{k[1]}) x={k[2]} br={k[3]}{k[4]} s=({k[5]},{k[6]},{k[7]})")
