"""Final5 per case, symbolic (2026-10-21).  L >= Lb := [B1 + C1 H [e1]]/S1 + midW + [B2 + C2 H [!e2]]/S2.
Write u = C/S (cot), so C = u S; prints Lb - 5 and Lb - (4 + 6 r - 4H) per case and sub-case (signs of mu, cr),
then numerically the min over t1, t2 in (0, 90), H in [2r - 1, 2) of max(Lb - 5, Lb + 4H - 4 - 6r)."""
import itertools, math
import numpy as np
import sympy as sp

C1, S1, C2, S2, H = sp.symbols('C1 S1 C2 S2 H', positive=True)
r = sp.sqrt(2)

def B(C, S, br):
    return 2 * S - C * (H - 1) if br == 'U' else (C * (S - C) + 1 + C + S - H) / C

def cr(a, b): return a[0] * b[1] - a[1] * b[0]
def dt(a, b): return a[0] * b[0] + a[1] * b[1]

def cases():
    for e1, e2, xi, b1, b2 in itertools.product((0, 1), (0, 1), ('up', 'a1', 'a2'), 'UP', 'UP'):
        a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
        v1 = C1 if e1 else -C1; v2 = C2 if e2 else -C2
        x = {'up': (1, 0), 'a1': a1, 'a2': a2}[xi]
        mu = v2 / S2 - v1 / S1
        for smu in (1, -1):
            y = (x[0] + x[1]) / 2 if smu > 0 else H - (x[0] + x[1]) / 2
            c1 = cr(a1, x); c2 = cr(a2, x)
            for s1 in ((1, -1) if c1 != 0 else (1,)):
                for s2 in ((1, -1) if c2 != 0 else (1,)):
                    m = mu * y + (s1 * c1 + dt(a1, x)) / (2 * S1) + (s2 * c2 + dt(a2, x)) / (2 * S2)
                    Lb = (B(C1, S1, b1) + (C1 * H if e1 else 0)) / S1 + m + (B(C2, S2, b2) + (0 if e2 else C2 * H)) / S2
                    yield (e1, e2, xi, b1, b2, smu, s1, s2), mu, c1, c2, Lb

R2 = math.sqrt(2)
if __name__ == '__main__':
    ts = np.radians(np.linspace(0.02, 89.98, 300))
    T1, T2 = np.meshgrid(ts, ts, indexing='ij')
    out = []
    for key, mu, c1, c2, Lb in cases():
        f = sp.lambdify((C1, S1, C2, S2, H), Lb, 'numpy')
        fmu = sp.lambdify((C1, S1, C2, S2, H), mu, 'numpy')
        fc1 = sp.lambdify((C1, S1, C2, S2, H), c1, 'numpy'); fc2 = sp.lambdify((C1, S1, C2, S2, H), c2, 'numpy')
        worst = (math.inf, None)
        for h in np.linspace(2 * R2 - 1, 1.9999, 25):
            args = (np.cos(T1), np.sin(T1), np.cos(T2), np.sin(T2), h)
            ok = (key[5] * fmu(*args) >= 0) & (key[6] * fc1(*args) * np.ones_like(T1) >= 0) & (key[7] * fc2(*args) * np.ones_like(T1) >= 0)
            if not ok.any(): continue
            L = f(*args) * np.ones_like(T1)
            sl = np.maximum(L - 5, L + 4 * h - 4 - 6 * R2)
            sl = np.where(ok, sl, np.inf)
            i = np.unravel_index(np.argmin(sl), sl.shape)
            if sl[i] < worst[0]:
                worst = (sl[i], (round(h, 4), round(math.degrees(T1[i]), 2), round(math.degrees(T2[i]), 2)), float(L[i]))
        out.append((worst, key))
    for w, k in sorted(out, key=lambda t: t[0][0]):
        if w[1] is None: continue
        print(f"{w[0]:+.4f} L={w[2]:.4f} at {w[1]}  e=({k[0]},{k[1]}) x={k[2]} br={k[3]}{k[4]} sgn(mu,cr1,cr2)=({k[5]},{k[6]},{k[7]})")
