import numpy as np, sys
sys.argv=['x']; exec(open('lemP.py').read().split('worst=(9,)')[0])   # brute-force P(beta,H)
def hs(t,phi): d=phi-t; return (np.abs(np.cos(d))+np.abs(np.sin(d)))/2
def Pc(beta,H):
    kb=hs(beta,0); out=1+kb
    for phi in (beta,beta-np.pi/2,beta+np.pi/2):
        c,s=np.cos(phi),np.sin(phi)
        if c<=1e-12: continue
        for ya in (0.5,H-0.5):
            out=min(out,0.5+(kb+0.5-s*(kb-ya))/c)
    return out
r2=np.sqrt(2); md=0; ms=(9,)
for H in np.linspace(1,1.999,120):
    m=min(1.5,2+r2-H)
    for bd in np.linspace(-45,45,361):
        b=np.radians(bd)
        if 2*hs(b,0)>H: continue
        md=max(md,abs(P(b,H,721)-Pc(b,H)))
        s=Pc(b,H)+Pc(-b,H)-2*m
        if s<ms[0]: ms=(s,H,bd)
print('max |P_brute - P_closed|',md); print('min slack closed form',ms)
