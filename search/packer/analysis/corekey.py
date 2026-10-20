"""Core-only identity vs the long-polish verdict on analysis/ident.py's calibration pairs (10-20; Evan: "fingerprint
for high precision dedupe might need to be over the rigid core only, or ignoring rattlers").  Core = fixlp.mobility < 1e-7
at the row's own point (no first-order motion at fixed side; rattlers and flat-motion squares excluded).  Pair test:
best of 8 symmetries + Hungarian on core centres; same core = equal core size and max centre / angle deviation <= tol."""
import sys, os, json, sqlite3
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rational'))
import numpy as np
from scipy.optimize import linear_sum_assignment
import fixlp
from core.coremap import D4, angd

db = sqlite3.connect('file:runs/store.sqlite?mode=ro', uri=True)
def row(pid):
    s, c = db.execute('SELECT s, coords FROM packing WHERE id=?', (pid,)).fetchone()
    return s, np.frombuffer(c, '<f8').reshape(-1, 3).copy()
def corepos(pid):
    s, sq = row(pid)
    return s, sq, fixlp.mobility(sq, s) < 1e-7
def compare(a, b):
    sa, qa, ca = a; sb, qb, cb = b
    A = qa[ca]; best = None
    for Bf in D4(qb, sb):
        B = Bf[cb]
        if len(B) != len(A):
            C = np.sqrt(((A[:, None, :2] - B[None, :, :2]) ** 2).sum(-1)); r, c = linear_sum_assignment(C)
        else:
            C = np.sqrt(((A[:, None, :2] - B[None, :, :2]) ** 2).sum(-1)); r, c = linear_sum_assignment(C)
        dp = C[r, c].max(); da = angd(A[r, 2], B[c, 2]).max()
        if best is None or dp < best[0]:
            best = (dp, da)
    return int(ca.sum()), int(cb.sum()), best[0], best[1]
def compare2(a, b):
    """core of each side matched into ALL squares of the other (core membership is threshold-sensitive): max centre /
    angle deviation over core squares, best of 8 symmetries, both directions."""
    sa, qa, ca = a; sb, qb, cb = b
    best = None
    for k, Bf in enumerate(D4(qb, sb)):
        C = np.sqrt(((qa[:, None, :2] - Bf[None, :, :2]) ** 2).sum(-1)); r, c = linear_sum_assignment(C)
        dev = C[r, c]; ang = angd(qa[r, 2], Bf[c, 2])
        core = ca[r] & cb[c]          # locked in both (squares on flat motions drift between copies)
        if core.sum() < 0.5 * min(ca.sum(), cb.sum()):          # a wrong symmetry: the cores barely overlap
            continue
        v = (dev[core].max(), ang[core].max(), dev[~core].max() if (~core).any() else 0.0)
        if best is None or v[0] < best[0]:
            best = v
    return best or (np.inf, np.inf, np.inf)


if __name__ != '__main__':
    pass
else:
  D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'runs/ident_274_cal.json'))
  out = []
  for d in D:
      if d['pi'] is None or d['pj'] is None:
          continue
      A_, B_ = corepos(d['i']), corepos(d['j']); na, nb = int(A_[2].sum()), int(B_[2].sum()); dp, da, dn = compare2(A_, B_)
      same_pol = abs(d['pi'] - d['pj']) < 1e-7
      out.append((d['d'], same_pol, na, nb, dp, da, abs(d['si'] - d['sj'])))
  print('fp_dist  polish_same  core_a core_b  max_dpos  max_dang  side_gap')
  for o in sorted(out, key=lambda o: (not o[1], o[4])):
      if o[1] or o[4] < 1e-3:
          print(f'{o[0]:.4f}  {str(o[1]):5s}  {o[2]:4d} {o[3]:4d}  {o[4]:.2e}  {o[5]:.2e}  {o[6]:.2e}')
  for tol in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
      tp = sum(o[1] and o[4] <= tol for o in out); fp = sum((not o[1]) and o[4] <= tol for o in out)
      print(f'tol {tol:g}: core-same & polish-same {tp} / {sum(o[1] for o in out)}; core-same but polish-distinct {fp} / {sum(not o[1] for o in out)}')
