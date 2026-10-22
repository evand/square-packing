import math, numpy as np
from phase_lp import solve
from kwalk import Kmat
H=1.9
phis=np.linspace(-84,84,337)
res,ms=solve(H,phis,True)
M=len(ms)
d=-res.ineqlin.marginals
idx=np.argsort(-d)[:25]
# reconstruct constraint meaning: order of generation
lab=[]
for j in (0,1):
    for i in range(M):
        for k in range(M):
            u,v=j*M+i,(1-j)*M+k
            lab.append(('step',j,phis[i],phis[k],'A'))
            if u!=v: lab.append(('step',j,phis[i],phis[k],'B'))
for i in range(M):
    for j in (1,0):
        lab.append(('wall',j,phis[i],None,'A')); lab.append(('wall',j,phis[i],None,'B'))
print(len(lab),len(d))
for q in idx: print(round(d[q],4),lab[q])
p0=res.x[:M]; p1=res.x[M:2*M]
for ph in (-60,-45,-30,-15,0,15,30,45,60):
    i=int(np.argmin(abs(phis-ph))); print(ph, round(p0[i],3), round(p1[i],3))
