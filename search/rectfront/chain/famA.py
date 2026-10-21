# Family A: separating normal phi is one of a's (phi = alpha or alpha-90, cos phi > 0).
# term = k_a + (1/2 + h_b(phi) - sin(phi)(k_b - y_a))/cos(phi), y_a in {k_a, H-k_a}.  Compare with Pc(beta) components.
import numpy as np
r2=np.sqrt(2)
def hs(t,phi): d=phi-t; return (np.abs(np.cos(d))+np.abs(np.sin(d)))/2
N=721
A=np.radians(np.linspace(0,90,N,endpoint=False))
worst={}
for H in np.linspace(1,1.999,60):
    for bd in np.linspace(-45,45,181):
        b=np.radians(bd); kb=hs(b,0)
        if 2*kb>H: continue
        comps={'X':1+kb}
        for nm,phi in (('b0',b),('b-90',b-np.pi/2),('b+90',b+np.pi/2)):
            c,s=np.cos(phi),np.sin(phi)
            if c>1e-12: comps[nm]=min(.5+(kb+.5-s*(kb-ya))/c for ya in (.5,H-.5))
        Pc=min(comps.values())
        ka=hs(A,0); ok=2*ka<=H
        for phi in (A,A-np.pi/2):
            c,s=np.cos(phi),np.sin(phi); g=ok&(c>1e-9)
            for ya in (ka,H-ka):
                t=ka+(.5+hs(b,phi)-s*(kb-ya))/np.where(g,c,1)
                # which single component does t dominate everywhere?
                for nm,v in comps.items():
                    d=np.where(g,t-v,np.inf).min()
                    worst.setdefault(nm,[]).append(d)
                d=np.where(g,t-Pc,np.inf).min(); worst.setdefault('Pc',[]).append(d)
for nm,v in worst.items(): print(nm,'min(term - comp)',min(v))
