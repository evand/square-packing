#!/usr/bin/env python3
"""One-square lemma by witness case (2026-10-21).  For b with frame (P, Q) (tilt th), unit direction (C, S)
(tilt t), strip H: minimise T = S x - C y + h/2 over (y, lo, A) by LP in each witness case V / U / D, and
compare with Psi = S Phi (frontier branch) and Up = 2S - C(H - 1) (upright at the top), Weak = 2S - CH/2."""
import math, sys
import numpy as np
from scipy.optimize import linprog

def bounds_case(case, P, Q, C, S, H):
    k = (P + Q) / 2
    h = abs(S * P - C * Q) + S * Q + C * P
    # variables z = (y, lo, A);  T = S (1 + k + A) - C y + h/2
    c = [-C, 0, S]; c0 = S * (1 + k) + h / 2
    Aub, bub = [], []
    def le(row, rhs): Aub.append(row); bub.append(rhs)
    le([-1, 0, 0], -k); le([1, 0, 0], H - k)            # k <= y <= H - k
    le([0, -1, 0], 0); le([0, 1, 0], H - 1)             # 0 <= lo <= H - 1
    le([0, 0, -1], 1)                                   # x >= k  (A >= -1)
    if case == 'V':
        le([0, 0, -1], 0)
    elif case == 'U':                                   # y + (P-Q)/2 <= lo <= y + k; Q A + P (lo - y - (P-Q)/2) >= 0
        le([1, -1, 0], -(P - Q) / 2); le([-1, 1, 0], k)
        le([P, -P, -Q], -P * (P - Q) / 2)
    elif case == 'D':                                   # lo + 1 <= y + (P-Q)/2; y - k <= lo + 1; P A + Q (y + (P-Q)/2 - lo - 1) >= 0
        le([-1, 1, 0], (P - Q) / 2 - 1); le([1, -1, 0], 1 + k)
        le([-Q, Q, -P], Q * ((P - Q) / 2 - 1))
    r = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * 3, method='highs')
    if r.status != 0: return None
    return r.fun + c0

def targets(C, S, H):
    kt = (C + S) / 2
    Psi = S + 0.5 - (H - kt) * C + S * (0.5 - (H - kt - 1) * S) / C
    return Psi, 2 * S - C * (H - 1), 2 * S - C * H / 2

if __name__ == '__main__':
    rows = []
    ths = np.radians(np.linspace(0, 90, 91)); tsd = np.linspace(1, 89, 45)
    for H in [1.05, 1.3, 1.6, 1.8, 1.9, 1.95, 1.99]:
        for td in tsd:
            t = math.radians(td); C, S = math.cos(t), math.sin(t)
            Psi, Up, Weak = targets(C, S, H)
            for case in 'VUD':
                vals = [(bounds_case(case, math.cos(th), math.sin(th), C, S, H), th) for th in ths]
                vals = [(v, th) for v, th in vals if v is not None]
                if not vals: continue
                v, th = min(vals)
                rows.append((case, H, td, v - Psi, v - Up, v - Weak, v - min(Psi, Weak), math.degrees(th)))
    for case in 'VUD':
        R = [r for r in rows if r[0] == case]
        print(case, 'n=', len(R),
              ' min(v-Psi)=%.3g' % min(r[3] for r in R), ' min(v-Up)=%.3g' % min(r[4] for r in R),
              ' min(v-Weak)=%.3g' % min(r[5] for r in R), ' min(v-min(Psi,Weak))=%.3g' % min(r[6] for r in R))
        bad = [r for r in R if r[6] < -1e-9]
        for r in bad[:5]: print('   FAIL', r)
        nopsi = [r for r in R if r[3] < -1e-9]
        print('   cases below Psi: %d; of those below Up: %d; H range below Psi: %s' % (
            len(nopsi), sum(r[4] < -1e-9 for r in nopsi), sorted(set(r[1] for r in nopsi))))
