"""Adversarial covers for the zmx2 audit (ZMX2_AUDIT.md sec 3): writes into $ZT (default runs/zmx2_audit).
  x9943   candidate x 0.9943: exact mu = 0.9999956 at the known dip (must be refused)
  g1      asymmetric germ hole at (2.5,2.5): x=3 zeroed on [2.88,2.92], x=2 on [2.92,2.96]
  g2a     D4-symmetric germ hole: x=2 on [2.92,2.94] and its 7 images
  w1      wall-germ hole: x=1 zeroed on [2.94,3.00]
  wrap_s  s_num = 5 + 2^125 (i128 wrap of s_num*D back to 5000)
  wrap_w  three extra points at (0,0), weights 2^127-1, 2^127-1, 2 (aggregate wraps to 0)"""
import os, sys
from fractions import Fraction as Fr
from aud import load, write
from holes import zero
C = sys.argv[1] if len(sys.argv) > 1 else 'runs/line-cover_m5_candidate_x1003.txt'
ZT = os.environ.get('ZT', 'runs/zmx2_audit')
os.makedirs(ZT, exist_ok=True)
cv0 = load(C)
f = Fr(9943, 10000)
cv = dict(cv0)
cv['pts'] = [(a, b, int(w * f)) for a, b, w in cv0['pts']]
cv['segs'] = [(a, b, c, d, int(w * f)) for a, b, c, d, w in cv0['segs']]
write(ZT + '/x9943.txt', cv)
write(ZT + '/g1.txt', zero(cv0, [(3000, 2880, 3000, 2900), (3000, 2900, 3000, 2920), (2000, 2920, 2000, 2940),
                                 (2000, 2940, 2000, 2960)], sym=False)[0])
write(ZT + '/g2a.txt', zero(cv0, [(2000, 2920, 2000, 2940)])[0])
write(ZT + '/w1.txt', zero(cv0, [(1000, y, 1000, y + 20) for y in range(2940, 3000, 20)], sym=False)[0])
lines = open(C).read().split('\n')
i = next(k for k, l in enumerate(lines) if l.strip() == '5 1')
l2 = list(lines); l2[i] = '%d 1' % (5 + 2 ** 125)
open(ZT + '/wrap_s.txt', 'w').write('\n'.join(l2))
j = next(k for k, l in enumerate(lines) if l.strip() == str(len(cv0['pts'])))
l3 = list(lines); l3[j] = str(len(cv0['pts']) + 3)
l3.insert(j + 1, '0 0 %d\n0 0 %d\n0 0 2' % (2 ** 127 - 1, 2 ** 127 - 1))
open(ZT + '/wrap_w.txt', 'w').write('\n'.join(l3))
print('wrote', ZT)
