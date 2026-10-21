# P(beta,H): b on floor (y_b=k_b), a left of b (any tilt, y_a extreme), bound via each separating normal phi (min over phi is
# NOT what we want: x_b >= term(phi*) for the actual separating phi*, so P = min over alpha,y_a of max?? no:
# disjoint => exists phi with dx >= term(phi) => dx >= min_phi term(phi).  So P = min over alpha,y_a,phi term.)
import numpy as np
r2=np.sqrt(2)
def hs(t,phi): d=phi-t; return (np.abs(np.cos(d))+np.abs(np.sin(d)))/2
def P(beta,H,N=1441):
    al=np.radians(np.linspace(0,90,N,endpoint=False)); ka=hs(al,0); kb=hs(beta,0)
    best=np.inf
    for ya in (ka,H-ka):
        dy=kb-ya
        for phi in (al,al-np.pi/2,np.full_like(al,beta),np.full_like(al,beta-np.pi/2),np.full_like(al,beta+np.pi/2)):
            c=np.cos(phi); good=(c>1e-9)&(2*ka<=H)
            t=ka+(hs(al,phi)+hs(beta,phi)-np.sin(phi)*dy)/np.where(good,c,1)
            best=min(best,np.where(good,t,np.inf).min())
    return best
worst=(9,)
for H in np.concatenate([np.linspace(1,1.99,100),[0.5+r2]]):
    m=min(1.5,2+r2-H)
    for bd in np.linspace(0,45,181):
        b=np.radians(bd)
        if 2*hs(b,0)>H: continue
        s=P(b,H)+P(-b,H)-2*m
        if s<worst[0]: worst=(s,H,bd,P(b,H),P(-b,H))
print(worst)
