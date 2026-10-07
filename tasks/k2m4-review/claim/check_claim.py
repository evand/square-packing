# claim-area independent checks (read-only on public/): box file <-> Lean data, total, D4, polygon, duplicates
import re, sys
from fractions import Fraction as F
from collections import Counter
P='/home/evand/math/square-packing/public/s12/'
toks=[]
for ln in open(P+'search/qx2_data/K4_k008_box9.txt'):
    toks+=ln.split('#',1)[0].split()
assert toks[0]=='mixed' and toks[1]=='1'
t=list(map(int,toks[2:]))
sn,sd,Dc,W,npnt=t[:5]; i=5
assert (sn,sd,Dc,W,npnt)==(9,1,5,2000000000000,0), (sn,sd,Dc,W,npnt)
ns=t[i]; i+=1
segs=[tuple(t[i+5*j:i+5*j+5]) for j in range(ns)]; i+=5*ns
npg=t[i]; i+=1
assert npg==1
k=t[i]; pw=t[i+1]; verts=[(t[i+2+2*j],t[i+3+2*j]) for j in range(k)]; i+=2+2*k
assert i==len(t)
print('segments',ns,'poly',k,pw,verts)
# Lean data
L=open(P+'lean/Sqpack/Bentz4Data.lean').read()
body=L[L.index('def boxSegs'):L.index('def boxPolyVerts')]
lsegs=[tuple(map(int,m)) for m in re.findall(r'\((\d+), (\d+), (\d+), (\d+), (\d+)\)',body)]
print('lean segs == file segs (order too):', lsegs==segs)
lv=L[L.index('def boxPolyVerts'):L.index('def boxPolyW')]
print('lean verts', re.findall(r'\((\d+), (\d+)\)',lv), 'leanW', re.search(r'def boxPolyW : Nat := (\d+)',L).group(1))
# geometry
for (a,b,c,d,w) in segs:
    assert w>0
    assert (b==d and abs(a-c)==1) or (a==c and abs(b-d)==1)
    assert all(0<=v<=45 for v in (a,b,c,d))
# polygon: square [14/5,31/5]^2, mass = area
assert verts==[(14,14),(31,14),(31,31),(14,31)]
assert F(pw,W)==F(17,5)**2
tot=F(sum(s[4] for s in segs)+pw,W)
D=F(214770225571,200000000000)
print('total',tot,float(tot),'== 81-4D:',tot==81-4*D,' <77:',tot<77, 'D>1', D>1)
# duplicates
key=lambda s:(min((s[0],s[1]),(s[2],s[3])),max((s[0],s[1]),(s[2],s[3])))
c=Counter(key(s) for s in segs)
dups=[k_ for k_,v in c.items() if v>1]
print('duplicated unit segments:',len(dups),'max mult',max(c.values()))
for k_ in sorted(dups)[:16]: print('  ',k_,[s[4] for s in segs if key(s)==k_])
# D4 invariance of the measure (summed density per unit segment)
m=Counter()
for s in segs: m[key(s)]+=s[4]
S=45
def img(g):
    out=Counter()
    for (p,q),w in m.items():
        a,b=g(*p),g(*q); out[(min(a,b),max(a,b))]+=w
    return out
print('D4 invariant:', img(lambda x,y:(S-x,y))==m and img(lambda x,y:(y,x))==m)
# segments on the container walls?
wall=[s for s in segs if (s[0]==s[2] and s[0] in (0,45)) or (s[1]==s[3] and s[1] in (0,45))]
print('wall segments',len(wall), 'mass on walls', float(F(sum(s[4] for s in wall),W)))
