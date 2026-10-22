import numpy as np, sys
from phase_lp2 import solve
H=float(sys.argv[1]) if len(sys.argv)>1 else 1.9
phis=np.linspace(-84,84,337); M=len(phis)
res,ms=solve(H,phis)
d=-res.ineqlin.marginals
lab=[]
ii,kk=np.meshgrid(np.arange(M),np.arange(M),indexing='ij')
for j in (0,1):
    for i,k in zip(ii.ravel(),kk.ravel()):
        lab+= [f"step j={j} {phis[i]:+.1f}->{phis[k]:+.1f} A", f"step j={j} {phis[i]:+.1f}->{phis[k]:+.1f} B"]
for i in range(M):
    lab+=[f"start ->{phis[i]:+.1f} A",f"start ->{phis[i]:+.1f} B",f"end {phis[i]:+.1f}-> A",f"end {phis[i]:+.1f}-> B"]
for q in np.argsort(-d)[:30]:
    if d[q]<1e-6: break
    print(round(d[q],4),lab[q])
