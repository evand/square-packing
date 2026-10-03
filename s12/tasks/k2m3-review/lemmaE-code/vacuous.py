"""Shows the missing guard: if any alternative list were empty, combos = [] and _cand_bin returns 'ok' vacuously.
(Monkeypatch only; in the shipped code _pieces cannot return [] -- see REPORT.)"""
from fractions import Fraction as F
from mass import *
import qx2_zm as QZ
cov, a, b = load(); ex = QZ.Exact(cov, a, b)
box = (F(3), F(3) + F(1, 80), F(1, 2), F(1, 2) + F(1, 80), F(0), F(1, 500))
print('tau=2 (false claim), unpatched:', ex.certify(box, tau=F(2)))
ex._pieces = lambda *A: []
print('tau=2 with _pieces -> []:', ex.certify(box, tau=F(2)))
