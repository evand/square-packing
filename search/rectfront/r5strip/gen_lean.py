"""Emit Lean test files for the distinct leaves (reps.pkl): plain nlinarith with a standard hint set."""
import pickle, sys, sympy as sp
from cases_sym import C1, S1, C2, S2
from polys import r
reps = pickle.load(open('reps.pkl', 'rb'))
def lean(e):
    s = sp.sstr(sp.expand(e)).replace('**', '^')
    return s
nfiles = int(sys.argv[1]) if len(sys.argv) > 1 else 4
outdir = sys.argv[2]
for f in range(nfiles):
    lines = ['import Mathlib', '', 'set_option maxHeartbeats 400000', '']
    for i, (key, tag, tgt, h, P, cons) in enumerate(reps):
        if i % nfiles != f: continue
        hs = ' '.join(f'(k{j} : 0 ≤ {lean(c)})' for j, c in enumerate(cons))
        kh = ''.join(f', mul_nonneg k{j} hC1, mul_nonneg k{j} hS1, mul_nonneg k{j} hC2, mul_nonneg k{j} hS2' for j in range(len(cons)))
        lines.append(f'-- leaf {i}: {key} target {tgt} at H = {h}')
        lines.append(f'theorem leaf{i} (C1 S1 C2 S2 r : ℝ) (hC1 : 0 ≤ C1) (hS1 : 0 ≤ S1) (hC2 : 0 ≤ C2) (hS2 : 0 ≤ S2)')
        lines.append(f'    (e1 : C1 ^ 2 + S1 ^ 2 = 1) (e2 : C2 ^ 2 + S2 ^ 2 = 1) (hr : r ^ 2 = 2) (hr1 : 1.4142 < r) (hr2 : r < 1.4143) {hs} :')
        lines.append(f'    0 ≤ {lean(P)} := by')
        lines.append(f'  nlinarith [mul_nonneg hC1 hS1, mul_nonneg hC2 hS2, mul_nonneg hC1 hC2, mul_nonneg hS1 hS2, mul_nonneg hC1 hS2, mul_nonneg hS1 hC2, sq_nonneg (C1 - S1), sq_nonneg (C2 - S2), sq_nonneg (C1 - C2), sq_nonneg (S1 - S2), sq_nonneg (C1 - S2), sq_nonneg (S1 - C2){kh}]')
        lines.append('')
    open(f'{outdir}/Leaves{f}.lean', 'w').write('\n'.join(lines))
