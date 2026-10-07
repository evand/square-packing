import sys, multiprocessing as mp
from fractions import Fraction as F
import ev
fn=sys.argv[1]; lo=F(sys.argv[2]); hi=F(sys.argv[3]); k=int(sys.argv[4]); step=F(sys.argv[5]); R=F(sys.argv[6])
COV=ev.Cover(fn)
n=int(R/step); T=[(i*step,j*step) for i in range(-n,n+1) for j in range(-n,n+1)]
def work(a):
    cx,cy,sg=a; u0=F(1,10**k); u=sg*u0; best=(9,None); nb=0
    for tx,ty in T:
        x,y=cx+u0*tx,cy+u0*ty
        if not COV.admissible(x,y,u): continue
        m=COV.mass(x,y,u); nb+=m<1
        if m<best[0]: best=(m,(str(tx),str(ty)))
    return (float(best[0]-1),str(cx),str(cy),sg,best[1],nb)
odd=[F(2*i+1,10) for i in range(0,40) if lo<=F(2*i+1,10)<=hi]
args=[(a,b,s) for a in odd for b in odd for s in (1,-1) if a<=b]
with mp.Pool(2) as p: res=sorted(p.map(work,args))
print('corners',len(args),'with sub-1 poses',sum(r[5]>0 for r in res))
for r in res[:6]: print(r)
