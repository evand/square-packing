"""Which clips are needed (legibility pass 2026-10-25)."""
import math, numpy as np
import leg_closed as L
def glob(pq, c, s, H): return L.Sup(pq, c, s, H)
def test(variant):
    out = []
    for H in np.linspace(L.H9, 1.87, 5):
        lam, a, b, B = L.br(H); T0, T1 = B['T']; x2 = (2*b - a)/(2*b); y1 = T0/(b - T1); w = {}
        for th in np.linspace(0, math.pi/2, 1801):
            C, S = math.cos(th), math.sin(th); base = C + S - lam
            p0, p1, _, _ = L.bounds(H, C, S); c, s = S, C
            if variant['psi1'] == 'clipA':
                q1 = max(L.SupI(B['A1'], c, s, H, 0, y1), a + glob(B['T'], c, s, H))
            elif variant['psi1'] == 'clipT':
                q1 = max(glob(B['A1'], c, s, H), a + L.SupI(B['T'], c, s, H, y1, math.inf))
            else:
                q1 = max(L.SupI(B['A1'], c, s, H, 0, y1), a + L.SupI(B['T'], c, s, H, y1, math.inf))
            v = variant['psi0']
            pA = L.SupI(B['A0'], c, s, H, 0, x2) if 'A' in v else glob(B['A0'], c, s, H)
            pL = L.SupI(B['Lz'], c, s, H, x2 if 'l' in v else 0, 8/7 if 'L' in v else math.inf) if True else 0
            pT = L.SupI(B['T'], c, s, H, 8/7, math.inf) if 'T' in v else glob(B['T'], c, s, H)
            q0 = max(pA, pL, pT)
            for kk, vv in dict(Sa=p0 + base - q1, Sb=p1 + base - q0).items():
                if kk not in w or vv < w[kk][0]: w[kk] = (vv, math.degrees(th))
        out.append(f"{w['Sa'][0]:+.1e}/{w['Sb'][0]:+.1e}")
    print(variant, ' '.join(out))
for p1 in ('clipA', 'clipT', 'both'): test({'psi1': p1, 'psi0': 'AlLT'})
for p0 in ('AlLT', 'lLT', 'ALT', 'AlL', 'AlT', 'LT', 'AL', 'lL', 'A', 'L'): test({'psi1': 'both', 'psi0': p0})
