import sys, json, multiprocessing as mp
from fractions import Fraction as F
import ev
P='/home/evand/math/square-packing/public/s12/search/qx2_data/'
fn=sys.argv[1]; k=int(sys.argv[2]); step=F(sys.argv[3]); R=F(sys.argv[4]); src=sys.argv[5]; out=sys.argv[6]
COV=ev.Cover(P+fn)
n=int(R/step); T=[(i*step,j*step) for i in range(-n,n+1) for j in range(-n,n+1)]
def work(a):
    cx,cy,sg=F(a[0]),F(a[1]),a[2]; u0=F(1,10**k); u=sg*u0
    vals=[]
    for tx,ty in T:
        x,y=cx+u0*tx,cy+u0*ty
        if not COV.admissible(x,y,u): continue
        vals.append((COV.mass(x,y,u),tx,ty))
    vals.sort()
    below=[v for v in vals if v[0]<1]
    return (a[0],a[1],sg,str(vals[0][0]),float(vals[0][0]-1),str(vals[0][1]),str(vals[0][2]),len(below),len(vals))
r=[json.loads(l) for l in open(src)]
C=[(x[0],x[1],x[2]) for x in r if 0<float(F(x[5])-1)<1e-10 and F(x[0])<=F(x[1])]
print(len(C),'corners',flush=True)
with mp.Pool(2) as p:
    res=list(p.imap_unordered(work,C))
res.sort(key=lambda z:z[4])
with open(out,'w') as f:
    for z in res: f.write(json.dumps(z)+'\n'); print(z)
