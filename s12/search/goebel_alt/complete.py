# Is every layout reaching 149 in class H? (5 blocks, integer corners, sorted x gaps in {2,3} from 0 to 10; same for y.)
# By the x<->y symmetry it suffices to show no 149-layout fails the x-condition.
import itertools, numpy as np
from model import *
inf = highspy.kHighsInf
def seqs():
    out=[]
    for st in itertools.product((2,3),repeat=4):
        if sum(st)==10: out.append(tuple(np.cumsum((0,)+st)*2))   # as indices
    return out
X = seqs(); print('x-sets (indices):', X)
def case(name, add):
    h, sq, bl = build(12)
    cols = {}
    for (a,b),j in bl.items(): cols.setdefault(a,[]).append(j)
    add(h, bl, cols); val,_ = solve(h); print(f'{name}: max = {val:.0f}', flush=True)
def addrow(h, idx, coef, lo, hi): h.addRow(lo, hi, len(idx), np.array(idx), np.array(coef, float))
allb = lambda bl: list(bl.values())
case('<=4 blocks', lambda h,bl,c: addrow(h, allb(bl), [1]*len(bl), -inf, 4))
case('>=6 blocks', lambda h,bl,c: addrow(h, allb(bl), [1]*len(bl), 6, inf))
def odd(h,bl,c):
    addrow(h, allb(bl), [1]*len(bl), 5, 5)
    j=[v for (a,b),v in bl.items() if a%2]; addrow(h, j,[1]*len(j),1,inf)
case('5 blocks, some odd x', odd)
def stack(h,bl,c):
    # 5 blocks, all even x, some column used twice: forbid by requiring max over a; do via extra binary
    addrow(h, allb(bl), [1]*len(bl), 5, 5)
    j=[v for (a,b),v in bl.items() if a%2]; addrow(h, j,[1]*len(j),-inf,0)
    n0 = h.getNumCol(); A=sorted(c)
    h.addVars(len(A), np.zeros(len(A)), np.ones(len(A)))
    h.changeColsIntegrality(len(A), np.arange(n0,n0+len(A)), np.array([highspy.HighsVarType.kInteger]*len(A)))
    for i,a in enumerate(A):   # w_a=1 -> cx_a >= 2
        addrow(h, c[a]+[n0+i], [1]*len(c[a])+[-2], 0, inf)
    addrow(h, list(range(n0,n0+len(A))), [1]*len(A), 1, inf)
case('5 blocks, two share an x', stack)
def other(h,bl,c):
    addrow(h, allb(bl), [1]*len(bl), 5, 5)
    j=[v for (a,b),v in bl.items() if a%2]; addrow(h, j,[1]*len(j),-inf,0)
    for a in c: addrow(h, c[a], [1]*len(c[a]), -inf, 1)
    for s in X:   # x-set != s : sum_{a not in s} cx_a >= 1
        j=[v for a in c if a not in s for v in c[a]]; addrow(h, j, [1]*len(j), 1, inf)
case('5 blocks, distinct even x, x-set not in the 6', other)
