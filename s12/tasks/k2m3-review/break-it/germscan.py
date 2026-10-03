import sys, itertools, multiprocessing as mp, json
from fractions import Fraction as F
import ev
P='/home/evand/math/square-packing/public/s12/search/qx2_data/'
fn=sys.argv[1]; k=int(sys.argv[2]); step=F(sys.argv[3]); R=F(sys.argv[4])
COV=ev.Cover(P+fn)
T=[]
n=int(R/step)
for i in range(-n,n+1):
    for j in range(-n,n+1): T.append((i*step,j*step))
def work(arg):
    cx,cy,sg=arg
    u=sg*F(1,10**k); best=(9,None)
    for tx,ty in T:
        x=cx+u*tx*sg*sg if False else cx+abs(u)*tx; y=cy+abs(u)*ty
        if not COV.admissible(x,y,u): continue
        m=COV.mass(x,y,u)
        if m<best[0]: best=(m,(str(tx),str(ty)))
    return (str(cx),str(cy),sg,float(best[0]),best[1], str(best[0]))
odd=[F(2*i+1,10) for i in range(2,18)]  # 0.5..3.5
args=[(a,b,s) for a in odd for b in odd for s in (1,-1)]
with mp.Pool(2) as pool:
    res=list(pool.imap_unordered(work,args,chunksize=4))
res.sort(key=lambda r:r[3])
out=open(f'germ_{fn[:6]}_k{k}.jsonl','w')
for r in res: out.write(json.dumps(r)+'\n')
for r in res[:25]: print(r[:5])
