import sys, numpy as np
from model import *
K = int(sys.argv[1]); target = int(sys.argv[2])
# 1. Ellsworth's layout for K=12: corners (0,0),(3,2),(5,5),(8,7),(10,10) -> indices 2x
if K == 12:
    E = {(0,0),(6,4),(10,10),(16,14),(20,20)}
    h,_,_ = build(K, E); print('Ellsworth layout value', solve(h)[0])
h, sq, bl = build(K)
layouts = []
while True:
    val, sol = solve(h)
    if val < target - 0.5: break
    L = frozenset(k for k,j in bl.items() if sol[j] > .5)
    layouts.append((val, L)); print(round(val), sorted((a/2, b/2) for a,b in L), flush=True)
    idx = [bl[k] for k in L]; others = [j for k,j in bl.items() if k not in L]
    # no-good cut: sum_{in L} y - sum_{not in L} y <= |L|-1
    allidx = idx+others; coef = [1.0]*len(idx)+[-1.0]*len(others)
    h.addRow(-highspy.kHighsInf, len(idx)-1, len(allidx), np.array(allidx), np.array(coef))
print('layouts reaching', target, ':', len(layouts))
