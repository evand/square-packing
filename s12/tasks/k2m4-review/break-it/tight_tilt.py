# For each tight germ (exact limit mass 1), sample nearby tiny-tilt poses consistent with/near the germ and look for mass<1
import numpy as np, json, sys, time
from fastf import mass, feasible
part,nparts=int(sys.argv[1]),int(sys.argv[2])
G=[tuple(float(eval(v.replace('/','/'))) if i>0 else None for i,v in enumerate(r)) for r in json.load(open('germs_near.json'))]
G=[g for g in G]
rows=[r for r in json.load(open('germs_near.json')) if r[0]=='1']
rows=rows[part::nparts]
rng=np.random.default_rng(100+part)
worst=[];nv=0;t0=time.time();tot=0
K=240
for r in rows:
    x0,y0,sg,X,Y=[float(eval(v)) for v in r[1:]]
    tau=np.repeat([1e-4,1e-3,1e-2,3e-2],K//4)
    th=sg*tau
    Xs=X+rng.uniform(-0.25,0.25,K)*(rng.random(K)<0.7)
    Ys=Y+rng.uniform(-0.25,0.25,K)*(rng.random(K)<0.7)
    xi1=rng.uniform(-0.2,0.2,K)*(rng.random(K)<0.25)
    xi2=rng.uniform(-0.2,0.2,K)*(rng.random(K)<0.25)
    cx0=x0+.5; cy0=y0+.5
    dy=cx0-Xs; dx=Ys-cy0
    cx=cx0+xi1+tau*dx*0+th*dx; cy=cy0+xi2+th*dy
    # sign of th: formulas derived for theta>0 with offsets theta*dx; for theta<0 cut directions swap automatically
    ok=feasible(cx,cy,th)
    if ok.sum()==0: continue
    m=mass(cx[ok],cy[ok],th[ok]); tot+=ok.sum()
    rel=(m-1)/np.abs(th[ok])
    i=np.argmin(rel)
    nv+=int((m<1-1e-12).sum())
    worst.append((rel[i],m[i],cx[ok][i],cy[ok][i],th[ok][i],r))
worst.sort(key=lambda w:w[0])
print('germs',len(rows),'poses',tot,'violations',nv,'time',time.time()-t0)
for w in worst[:15]: print('slope %.6f mass %.15f pose %r %r %r germ %s'%(w[0],w[1],w[2],w[3],w[4],w[5]))
