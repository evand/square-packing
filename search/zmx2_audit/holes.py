from aud import *
from germ import germscan
import sys
def imgs(seg, S):
    X0,Y0,X1,Y1 = seg
    out=set()
    for f in [lambda x,y:(x,y),lambda x,y:(S-x,y),lambda x,y:(x,S-y),lambda x,y:(S-x,S-y),
              lambda x,y:(y,x),lambda x,y:(S-y,x),lambda x,y:(y,S-x),lambda x,y:(S-y,S-x)]:
        a=f(X0,Y0); b=f(X1,Y1)
        out.add(tuple(sorted([a,b])))
    return out
def zero(cv, keys, sym=True):
    S=int(cv['s']*cv['D'])
    Z=set()
    for k in keys:
        Z |= imgs(k,S) if sym else {tuple(sorted([(k[0],k[1]),(k[2],k[3])]))}
    n=0; segs=[]
    for (a,b,c,d,w) in cv['segs']:
        if tuple(sorted([(a,b),(c,d)])) in Z: w=0; n+=1
        segs.append((a,b,c,d,w))
    cv=dict(cv); cv['segs']=segs; return cv, n
def piece(cv, X, y):
    """the vertical piece of line x=X containing height y (units 1/D)"""
    for (a,b,c,d,w) in cv['segs']:
        if a==c==X and min(b,d)<=y<max(b,d): return (a,min(b,d),c,max(b,d)), w
