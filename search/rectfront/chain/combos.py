import numpy as np
r2=np.sqrt(2)
t=np.radians(np.linspace(1e-6,45,20001)); c,s=np.cos(t),np.sin(t); k=(c+s)/2
def terms(H):
    X=1+k
    U=.5+(k+.5-s*(k-.5))/c;  V=.5+(k+.5-s*(H-.5-k))/c
    Up=.5+(k+.5-c*(k-.5))/s; Vp=.5+(k+.5-c*(H-.5-k))/s
    return {'X':X,'U':U,"V'":Vp},{'X':X,'V':V,"U'":Up}
res={}
for Hn,Hf in (('lo',lambda:np.maximum(1,c+s)),('kink',lambda:np.full_like(t,.5+r2)),('2',lambda:np.full_like(t,2.0))):
    H=Hf(); m=np.minimum(1.5,2+r2-H); L,R=terms(H)
    for a in L:
        for b in R:
            g=L[a]+R[b]-2*m; i=np.argmin(g)
            res.setdefault((a,b),[]).append((Hn,round(g[i],5),round(np.degrees(t[i]),2)))
for kk,v in res.items(): print(kk,v)
