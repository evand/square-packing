"""Final5 leaves for H in [2r-1, 2) (2026-10-21): per case/sub-case, the polynomial inequalities at H endpoints.
Class (which_target.py): '5' -> N - 5D >= 0 at H = 2r-1, 2;  'zig' -> N - zD >= 0 at 2r-1, 2;
'BOTH' -> N - 5D at 2r-1, h0 and N - zD at h0, 2 (h0 = (6r-1)/4).  Grouped up to the swap 1 <-> 2."""
import math, pickle
import sympy as sp
from cases_sym import cases, C1, S1, C2, S2, H
from polys import poly_case, r

def classes():
    cl = {}
    for line in open('../../../runs/rectfront/r5strip_targets.txt'):
        p = line.split()
        tag = p[0]
        e = p[4][3:-1].split(','); x = p[5][2:]; br = p[6][3:]; s = p[7][3:-1].split(',')
        cl[(int(e[0]), int(e[1]), x, br[0], br[1], int(s[0]), int(s[1]), int(s[2]))] = tag
    return cl

H0 = (6 * r - 1) / 4; HL = 2 * r - 1
def leaves():
    cl = classes(); out = []
    for key, mu, c1, c2, Lb in cases():
        if key not in cl: continue
        D, N, cons = poly_case(key, mu, c1, c2, Lb)
        z = 4 + 6 * r - 4 * H
        tag = cl[key]
        if tag == '5': L = [(5, HL), (5, 2)]
        elif tag == 'zig': L = [(z, HL), (z, 2)]
        else: L = [(5, HL), (5, H0), (z, H0), (z, 2)]
        for tgt, h in L:
            P = sp.expand((N - tgt * D).subs(H, h))
            out.append((key, tag, str(tgt), h, P, cons))
    return out

if __name__ == '__main__':
    L = leaves()
    sw = {C1: C2, S1: S2, C2: C1, S2: S1}
    groups = {}
    for key, tag, tgt, h, P, cons in L:
        k1 = (sp.srepr(P), tuple(sorted(map(sp.srepr, cons))))
        k2 = (sp.srepr(sp.expand(P.xreplace(sw))), tuple(sorted(sp.srepr(sp.expand(c.xreplace(sw))) for c in cons)))
        g = k1 if k1 in groups else (k2 if k2 in groups else k1)
        groups.setdefault(g, []).append((key, tgt, h))
    print(len(L), 'leaves,', len(groups), 'distinct up to 1<->2')
    reps = [(v[0], [x for x in L if (x[0], x[2], x[3]) == v[0]][0]) for v in groups.values()]
    pickle.dump([(r_[1][0], r_[1][1], r_[1][2], r_[1][3], r_[1][4], r_[1][5]) for r_ in reps], open('reps.pkl', 'wb'))
