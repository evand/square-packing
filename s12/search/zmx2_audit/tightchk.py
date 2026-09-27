import sys, random
from fractions import Fraction as Fr
from aud import load, mu, inside
cv=load(sys.argv[1]); S=int(cv['s']*cv['D'])
cvr=dict(cv); cvr['pts']=[(X,S-Y,w) for (X,Y,w) in cv['pts']]; cvr['segs']=[(a,S-b,c,S-d,w) for (a,b,c,d,w) in cv['segs']]
rows=[l.split() for l in open(sys.argv[2])]
rng=random.Random(5); rng.shuffle(rows); rows=rows[:int(sys.argv[3])]
unit=cv['W']*20*2**30
nf=nc=0; gaps=[]
for pas,bx,bd in rows:
    b=[Fr(t) for t in bx.split(',')]; B=Fr(int(bd),unit); c=[cv,cvr][int(pas)]
    mn=None
    for k in range(20):
        if k<8: f=[Fr((k>>j)&1) for j in range(3)]
        else: f=[Fr(rng.randint(0,997),997) for _ in range(3)]
        P=tuple(b[2*i]+(b[2*i+1]-b[2*i])*f[i] for i in range(3))
        if not inside(c,*P): continue
        m=mu(c,*P); nc+=1
        if m<B: nf+=1; print('FAIL',pas,bx,P,float(m),float(B))
        mn=m if mn is None or m<mn else mn
    if mn is not None: gaps.append(float(mn-B))
gaps.sort()
print('tight leaves %d, poses %d, FAIL %d, min gap %.3g, median %.3g'%(len(rows),nc,nf,gaps[0],gaps[len(gaps)//2]))
